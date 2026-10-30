"""Weighted switch-region defects from explicit pair-probability matrices.

All fractions refer to nucleotide pairing states in thermodynamic ensembles.
They are neither probabilities of an entire structure nor predicted expression.
"""
from __future__ import annotations

import math
from dataclasses import asdict
from numbers import Integral

import numpy as np

from .architecture import architecture_layout, target_switch_structures
from .models import FunctionalDefectWeights


def _pair_targets(structure: str, expected_lengths: list[int]) -> np.ndarray:
    parts = structure.split("+")
    if list(map(len, parts)) != expected_lengths or set(structure) - set(".()+"):
        raise ValueError("Target structure does not match the ordered strand lengths")
    partner = np.arange(sum(expected_lengths), dtype=int)
    stack = []
    for index, symbol in enumerate("".join(parts)):
        if symbol == "(":
            stack.append(index)
        elif symbol == ")":
            if not stack:
                raise ValueError("Target structure contains unbalanced parentheses")
            other = stack.pop()
            partner[index], partner[other] = other, index
    if stack:
        raise ValueError("Target structure contains unbalanced parentheses")
    return partner


def _pair_matrix(value, size: int, label: str) -> np.ndarray:
    if value is None:
        raise ValueError(f"Missing {label} pair probability matrix")
    try:
        if hasattr(value, "to_array"):
            value = value.to_array()
        original = np.asarray(value)
        if np.iscomplexobj(original) or original.dtype.kind == "b":
            raise ValueError("Matrix must contain real probabilities")
        matrix = np.asarray(original, dtype=float)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"Invalid {label} pair probability matrix") from exc
    if matrix.shape != (size, size):
        raise ValueError(f"{label} matrix dimensions must be {size} by {size}")
    if not np.isfinite(matrix).all() or np.any(matrix < -1e-12) or np.any(matrix > 1 + 1e-12):
        raise ValueError(f"{label} matrix entries must be finite probabilities between zero and one")
    if not np.allclose(matrix, matrix.T, rtol=0, atol=1e-8):
        raise ValueError(f"{label} complex pair probability matrix must be symmetric")
    if not np.allclose(matrix.sum(axis=1), 1.0, rtol=0, atol=1e-7):
        raise ValueError(f"{label} matrix must be complete, with each row summing to one")
    return np.clip(matrix, 0.0, 1.0)


def _regional_context(architecture, reporter_length, off_pairs, on_pairs) -> dict:
    """Validate full matrices and resolve switch errors against the targets."""
    if isinstance(reporter_length, bool) or not isinstance(reporter_length, Integral) or reporter_length < 0:
        raise ValueError("Reporter length must be a nonnegative integer")
    for field in ("toehold_length", "stem_length", "upper_stem_length", "upper_stem2_length",
                  "bulge_length", "rbs_prefix_length", "aug_spacer"):
        value = getattr(architecture, field, 0)
        if isinstance(value, bool) or not isinstance(value, Integral) or value < 0:
            raise ValueError("Architecture lengths must be nonnegative integers")
    layout = architecture_layout(architecture)
    trigger_length = int(layout["len_trigger"])
    architecture.validate_for_trigger(trigger_length)
    off_target, on_target = target_switch_structures(architecture, int(reporter_length))
    switch_length = len(off_target)
    off_partner = _pair_targets(off_target, [switch_length])
    on_partner = _pair_targets(on_target, [trigger_length, switch_length])
    off = _pair_matrix(off_pairs, switch_length, "OFF")
    on = _pair_matrix(on_pairs, trigger_length + switch_length, "ON")
    off_indices = np.arange(switch_length)
    on_indices = off_indices + trigger_length
    errors = {
        "off": 1.0 - off[off_indices, off_partner],
        "on": 1.0 - on[on_indices, on_partner[on_indices]],
    }
    return {
        "layout": layout, "switch_length": switch_length, "off_partner": off_partner,
        "errors": errors, "totals": {state: float(values.sum()) for state, values in errors.items()},
        "targets": {"off": off_target, "on": on_target},
    }


def _region_diagnostic(context: dict, state: str, indices: list[int], name: str) -> dict:
    if not indices or any(index < 0 or index >= context["switch_length"] for index in indices):
        raise ValueError(f"The {name} region is empty or outside the switch")
    local_errors = context["errors"][state][indices]
    expected_incorrect = float(local_errors.sum())
    total = context["totals"][state]
    return {
        "indices_switch": indices, "nucleotide_count": len(indices),
        "expected_incorrect_nt": expected_incorrect,
        "normalized_defect": float(local_errors.mean()),
        "share_of_switch_defect": expected_incorrect / total if total else None,
    }


def _aggregate_regions(context: dict, index_sets: dict, weights) -> dict:
    normalized_weights = weights.fractions()
    raw_weights = {name: float(value) for name, value in asdict(weights).items()}
    regions = {name: _region_diagnostic(context, name.split("_", 1)[0], indices, name)
               for name, indices in index_sets.items()}
    result = {name: region["normalized_defect"] for name, region in regions.items()}
    for state in ("off", "on"):
        names = [name for name in raw_weights if name.startswith(state + "_")]
        scale = max(raw_weights[name] for name in names)
        # Normalize each state independently, so a small but nonzero state
        # remains measurable even if its overall relative weight underflows.
        scaled = {name: raw_weights[name] / scale for name in names} if scale else {}
        result[state] = (math.fsum(scaled[name] * result[name] for name in names)
                         / math.fsum(scaled.values())) if scaled else None
    result["balanced"] = (math.fsum(normalized_weights[name] * result[name] for name in regions)
                          / math.fsum(normalized_weights.values()))
    result.update(regions=regions,
                  switch_structural={state: value / context["switch_length"] for state, value in context["totals"].items()},
                  switch_expected_incorrect_nt=context["totals"],
                  per_base={state: values.tolist() for state, values in context["errors"].items()},
                  indexing="switch_0_based", target_structures=context["targets"],
                  weights=raw_weights, normalized_weights=normalized_weights)
    return result




def functional_defect_metrics(
    architecture, reporter_length: int, off_pairs, on_pairs,
    weights: FunctionalDefectWeights | None = None,
) -> dict:
    """Measure OFF toehold accessibility/lock and two ON opening regions.

    The OFF lock measures the expected partner of both arms of the stems
    flanking AUG (6 + 9 pairs in the classic 30-nt architecture). AUG and its
    opposite bulge remain unpaired in the target: their sequestration is
    assessed indirectly through the flanking stems, never by demanding AUG
    pairing. ON RBS–AUG includes the RBS loop prefix and all three AUG bases.
    The second ON region extends from that same prefix through the linker's
    last base. Their overlap is intentional: both measurements respond to
    errors before AUG, while only the extended one responds after AUG. An ON
    trigger-duplex error alone does not contribute to either measurement.
    """
    weights = FunctionalDefectWeights() if weights is None else weights
    if not isinstance(weights, FunctionalDefectWeights):
        raise ValueError("Functional defect weights must be FunctionalDefectWeights")
    weights.validate()
    context = _regional_context(architecture, reporter_length, off_pairs, on_pairs)
    layout, partner = context["layout"], context["off_partner"]
    before_end = layout["aug_start"] - layout["len_aug_spacer"]
    after_start = layout["aug_start"] + 3
    right_arms = {
        "before_aug": list(range(before_end - layout["len_upper"], before_end)),
        "after_aug": list(range(after_start, after_start + layout["len_upper2"] + layout["len_stem"])),
    }
    stem_parts = {}
    for name, right_indices in right_arms.items():
        if any(partner[index] == index for index in right_indices):
            raise ValueError("OFF stem region contains a nucleotide without a target partner")
        indices = sorted(set(right_indices) | {int(partner[index]) for index in right_indices})
        stem_parts[name] = _region_diagnostic(context, "off", indices, name)
        stem_parts[name]["pair_count"] = len(right_indices)
    index_sets = {
        "off_toehold": list(range(layout["len_t7"], layout["len_t7"] + layout["len_toehold"])),
        "off_stem": sorted({index for part in stem_parts.values() for index in part["indices_switch"]}),
        "on_rbs_aug": list(range(layout["rbs_prefix_start"], after_start)),
        "on_rbs_linker": list(range(layout["rbs_prefix_start"], layout["rbs_linker_end"])),
    }
    result = _aggregate_regions(context, index_sets, weights)
    result.update(scheme="functional_v2", off_stem_parts=stem_parts)
    return result
