"""Assess switch-designed extensions against the original aptamer/trigger.

This stage does not redesign sequences or change switch performance scores.
It reuses the direct extension workflow's conservation metrics at independently
specified system conditions, with one original reference calculation.
"""
from __future__ import annotations

import copy
import math
from numbers import Integral

from .extension_thermodynamics import extension_model_settings, open_extension_model
from .localization import tr
from .models import RunResult
from .scoring import candidate_score_100, sort_candidates
from .sequences import normalize_dna
from .trigger_extension import (
    MfeResult, _nucleotide_probabilities, _reference_pair_rows,
    ensemble_metrics, extension_split, map_reference_pairs_to_extended,
    padded_reference_structure, parse_pairs, score_mfe_candidate,
)


def _validate_analysis(analysis, aptamer_length: int, trigger_length: int) -> None:
    import numpy as np

    if [len(part) for part in analysis.mfe_structure.split("+")] != [aptamer_length, trigger_length]:
        raise ValueError("The thermodynamics engine returned a structure incompatible with the analyzed strands")
    parse_pairs(analysis.mfe_structure)
    if not all(math.isfinite(float(value)) for value in (analysis.mfe_energy_kcal_mol, analysis.ensemble_free_energy_kcal_mol)):
        raise ValueError("The thermodynamics engine returned non-finite energies")
    for value in (analysis.mfe_structure_probability, analysis.target_structure_probability):
        if not math.isfinite(float(value)) or not 0 <= value <= 1:
            raise ValueError("The thermodynamics engine returned an invalid structure probability")
    matrix = np.asarray(analysis.pair_matrix)
    size = aptamer_length + trigger_length
    if (matrix.shape != (size, size) or not np.isfinite(matrix).all()
            or np.any(matrix < -1e-8) or np.any(matrix > 1 + 1e-8)
            or not np.allclose(matrix, matrix.T, atol=1e-7)
            or not np.allclose(matrix.sum(axis=1), 1, atol=1e-5)):
        raise ValueError("The thermodynamics engine returned an invalid pair probability matrix")


def analyze_inverse_extensions(result: RunResult, import_path: str | None = None,
                               progress_callback=None, cancel_check=None, *, backend=None) -> None:
    """Populate conservation metrics for up to X already evaluated switches.

    Thermodynamic calls are not interruptible. Cancellation is checked before
    and after each call; only complete candidate analyses are retained. A
    failure is recorded in this stage's metadata and never discards the switch
    design results or fabricates thermodynamic values.
    """
    settings = result.request.inverse_extension
    if settings is None:
        return
    stage = {
        "status": "running", "reference": {}, "model": {},
        "selected_trials": [], "selected_count": 0, "analyzed_count": 0,
        "unique_trigger_count": 0, "unique_trigger_total": 0, "candidates": [],
        "score_scope": "conservation_of_original_aptamer_trigger_complex",
        "switch_scores_unchanged": True,
    }
    result.inverse_extension_analysis = stage
    for candidate in result.candidates:
        candidate.inverse_extension_metrics = {}

    def cancelled() -> bool:
        if cancel_check and cancel_check():
            stage["status"] = "cancelled"
            stage["candidates"].sort(key=lambda row: row["original_switch_rank"])
            return True
        return False

    def progress(done: int, total: int, message: str) -> None:
        if progress_callback:
            progress_callback(done, total, message)

    try:
        if cancelled():
            return
        if backend is None:
            from .thermodynamics.setup import get_backend
            backend = get_backend()
        if isinstance(settings.top_candidates, bool) or not isinstance(settings.top_candidates, Integral) or settings.top_candidates < 1:
            raise ValueError("top_candidates must be a positive integer")
        if isinstance(settings.target_length, bool) or not isinstance(settings.target_length, Integral):
            raise ValueError("target_length must be an integer")
        thermo = settings.thermo
        model_settings = extension_model_settings(thermo, backend=backend)
        aptamer = normalize_dna(result.request.sequences.aptamer_dna, field_name="aptamer")
        original = normalize_dna(result.request.sequences.effective_trigger(), field_name="original trigger")
        added_length = int(settings.target_length) - len(original)
        if added_length <= 0:
            raise ValueError("The inverse extension must be longer than the original trigger")
        prefix_length, suffix_length = extension_split(added_length, settings.extension_side)
        selected = [candidate for candidate in sort_candidates(
            result.candidates,
            defect_weight=getattr(getattr(result.request, "scoring", None), "defect", None),
        ) if not candidate.excluded][:int(settings.top_candidates)]
        groups = {}
        for rank, candidate in enumerate(selected, 1):
            extended = normalize_dna(candidate.trigger_seq, field_name="extended trigger")
            if len(extended) != int(settings.target_length):
                raise ValueError(f"Trial {candidate.trial}: extended trigger has an unexpected length")
            if extended[prefix_length:prefix_length + len(original)] != original:
                raise ValueError(f"Trial {candidate.trial}: the original trigger core was changed")
            groups.setdefault(extended, []).append((rank, candidate))
        stage.update(model=model_settings, selected_trials=[c.trial for c in selected], selected_count=len(selected),
                     unique_trigger_total=len(groups), extension_side=settings.extension_side,
                     target_length=int(settings.target_length), original_trigger_length=len(original))
        if not selected:
            stage["status"] = "completed"
            return
        total = len(groups) + 1
        progress(0, total, tr("Référence originale aptamère + trigger pour l’extension inverse…"))
        with open_extension_model(thermo, import_path, backend=backend) as engine:
            stage.update(copy.deepcopy(engine.metadata))
            if cancelled():
                return
            reference = engine.analyze(aptamer, original)
            if cancelled():
                return
            _validate_analysis(reference, len(aptamer), len(original))
            reference_mfe = MfeResult(reference.mfe_structure, reference.mfe_energy_kcal_mol)
            reference_pairs = parse_pairs(reference.mfe_structure)
            mapped_pairs = map_reference_pairs_to_extended(reference_pairs, len(aptamer), prefix_length)
            padded_structure = padded_reference_structure(reference.mfe_structure, added_length, settings.extension_side)
            stage["reference"] = {
                "aptamer_sequence": aptamer, "trigger_sequence": original,
                "mfe_structure": reference.mfe_structure,
                "mfe_energy_kcal_mol": reference.mfe_energy_kcal_mol,
                "ensemble_free_energy_kcal_mol": reference.ensemble_free_energy_kcal_mol,
                "mfe_structure_probability": reference.mfe_structure_probability,
                "mfe_pairs": [list(pair) for pair in sorted(reference_pairs)],
                "pair_probabilities": _reference_pair_rows(reference, len(aptamer)),
                "nucleotide_probabilities": _nucleotide_probabilities(reference),
            }
            if not any(i < len(aptamer) <= j for i, j in reference_pairs):
                stage["warnings"] = [tr("La structure MFE de référence ne présente aucune paire entre l'aptamère et le trigger ; la conservation de cette référence ne démontre pas leur liaison.")]
            for index, (extended, candidates) in enumerate(groups.items(), 1):
                if cancelled():
                    return
                progress(index, total, tr("Comparaison de l’extension inverse {index}/{total}…", index=index, total=len(groups)))
                if cancelled():
                    return
                analysis = engine.analyze(aptamer, extended, target=padded_structure)
                if cancelled():
                    return
                _validate_analysis(analysis, len(aptamer), len(extended))
                extension = extended[:prefix_length] + (extended[-suffix_length:] if suffix_length else "")
                mfe = score_mfe_candidate(extension, MfeResult(analysis.mfe_structure, analysis.mfe_energy_kcal_mol),
                                          reference_mfe, mapped_pairs, len(aptamer), original,
                                          extension_side=settings.extension_side)
                metrics = ensemble_metrics(mfe, reference, analysis, len(aptamer), len(original))
                metrics.update(
                    status="completed", reference_mfe_structure=reference.mfe_structure,
                    reference_mfe_energy_kcal_mol=reference.mfe_energy_kcal_mol,
                    reference_pairs=_reference_pair_rows(reference, len(aptamer), analysis, prefix_length),
                    nucleotide_probabilities=_nucleotide_probabilities(analysis),
                )
                metrics.update(copy.deepcopy(engine.metadata))
                for switch_rank, candidate in candidates:
                    row = copy.deepcopy(metrics)
                    row.update(trial=candidate.trial, original_switch_rank=switch_rank)
                    candidate.inverse_extension_metrics = row
                    stage["candidates"].append(copy.deepcopy(row))
                stage["analyzed_count"] += len(candidates)
                stage["unique_trigger_count"] += 1
        stage["candidates"].sort(key=lambda row: row["original_switch_rank"])
        stage["status"] = "completed"
        progress(total, total, tr("Comparaison des extensions inverses terminée."))
    except Exception as exc:
        stage["status"] = "failed"
        stage["error"] = str(exc)
        stage["candidates"].sort(key=lambda row: row["original_switch_rank"])


def _average_ranks(values: list[float], *, reverse: bool = False) -> list[float]:
    ordered = sorted(range(len(values)), key=lambda i: values[i], reverse=reverse)
    ranks = [0.0] * len(values)
    start = 0
    while start < len(ordered):
        end = start + 1
        while end < len(ordered) and values[ordered[end]] == values[ordered[start]]:
            end += 1
        rank = (start + 1 + end) / 2.0
        for index in ordered[start:end]:
            ranks[index] = rank
        start = end
    return ranks


def ranked_inverse_candidates(result: RunResult, performance_weight: float = 0.5) -> list[dict]:
    """Rank analyzed switches by a dimensionless compromise of rank positions.

    High switch score and low conservation impact are good. Each axis uses
    average ranks for ties, normalized as (rank−1)/(n−1) over this analyzed
    subset. Lower compromise is better. These are relative ranking positions,
    not probabilities, physical energies, or experimentally calibrated scores.
    """
    if isinstance(performance_weight, bool):
        raise ValueError("performance_weight must be finite and between 0 and 1")
    try:
        weight = float(performance_weight)
    except (TypeError, ValueError) as exc:
        raise ValueError("performance_weight must be finite and between 0 and 1") from exc
    if not math.isfinite(weight) or not 0 <= weight <= 1:
        raise ValueError("performance_weight must be finite and between 0 and 1")
    selected = []
    for original_rank, candidate in enumerate([c for c in sort_candidates(
        result.candidates,
        defect_weight=getattr(getattr(result.request, "scoring", None), "defect", None),
    ) if not c.excluded], 1):
        metrics = getattr(candidate, "inverse_extension_metrics", {})
        if metrics.get("status") != "completed":
            continue
        try:
            raw_score = getattr(candidate, "score_raw", None)
            performance = float(candidate.score if raw_score is None else raw_score)
            impact = float(metrics["final_score"])
        except (TypeError, ValueError, KeyError):
            continue
        if math.isfinite(performance) and math.isfinite(impact):
            selected.append((original_rank, candidate, metrics, performance, impact))
    if not selected:
        return []
    performance_ranks = _average_ranks([row[3] for row in selected], reverse=True)
    structure_ranks = _average_ranks([row[4] for row in selected])
    denominator = max(1, len(selected) - 1)
    rows = []
    for index, (original_rank, candidate, metrics, performance, impact) in enumerate(selected):
        p = (performance_ranks[index] - 1) / denominator
        s = (structure_ranks[index] - 1) / denominator
        pareto = not any(other[3] >= performance and other[4] <= impact
                         and (other[3] > performance or other[4] < impact) for other in selected)
        rows.append({
            "trial": candidate.trial, "original_switch_rank": original_rank,
            "performance_rank": performance_ranks[index], "structure_rank": structure_ranks[index],
            "performance_rank_normalized": p, "structure_rank_normalized": s,
            "compromise_score": weight * p + (1 - weight) * s,
            "performance_weight": weight, "pareto": pareto,
            "switch_score": candidate_score_100(
                candidate, score_scale=result.request.score_scale,
            ),
            "switch_score_raw": performance, "conservation_score": impact,
            "candidate": candidate, "metrics": copy.deepcopy(metrics),
        })
    return sorted(rows, key=lambda row: (row["compromise_score"], row["original_switch_rank"], row["trial"]))
