"""Marginal pairing motifs and sampled, functional secondary-structure families.

Family populations come from Boltzmann samples, never from multiplying marginal
pair probabilities. Representatives are actual sampled structures. ON is
conditional on the formed trigger + switch complex, not on tube occupancy.
"""
from __future__ import annotations

from dataclasses import asdict
from datetime import datetime, timezone
from functools import lru_cache
import hashlib
import json
import math
from numbers import Integral
import os
from pathlib import Path
import tempfile
from types import SimpleNamespace

import numpy as np

from .defect_metrics import _pair_matrix, _pair_targets, functional_defect_metrics
from .localization import tr
from .models import SwitchArchitecture


SCHEMA = 1


@lru_cache(maxsize=8)
def ensemble_sample_variant_keys(serialized_inputs: str) -> dict[str, int]:
    """Index supported sample counts without opening unrelated saved matrices.

    This is a filename index only. ``load_ensemble`` must still validate each
    selected artifact before any of its scientific contents are displayed.
    """
    inputs = json.loads(serialized_inputs)
    return {ensemble_key({**inputs, "num_sample": count}): count for count in range(100, 5001)}


@lru_cache(maxsize=32)
def _recorded_ensemble_inputs(path: str, modified_ns: int, size: int) -> dict:
    """Cache only small input metadata, not saved probability matrices."""
    recorded = json.loads(Path(path).read_text(encoding="utf-8"))
    return recorded.get("inputs", {}) if isinstance(recorded, dict) else {}


def recorded_native_ensemble_inputs(path: Path, request, candidate) -> dict | None:
    """Match recorded model, sample count and seeds without loading a runtime.

    The caller must still pass the returned inputs to ``load_ensemble`` before
    displaying or copying the artifact; this function validates identity only.
    """
    from .thermodynamics.types import ModelProvenance
    off_length, on_length = len(candidate.switch_seq), len(candidate.switch_seq) + len(candidate.trigger_seq)
    maximum_bytes = min(512 * 1024 * 1024, max(32 * 1024 * 1024,
        32 * (off_length ** 2 + on_length ** 2) + 2 * 5000 * (off_length + on_length + 128)))
    if path.is_symlink() or not path.is_file():
        return None
    details = path.stat()
    if details.st_size > maximum_bytes:
        return None
    inputs = _recorded_ensemble_inputs(str(path), details.st_mtime_ns, details.st_size)
    if not isinstance(inputs, dict) or "thermodynamics" not in inputs:
        return None
    provenance = ModelProvenance(**inputs["thermodynamics"])
    seed = inputs["sampling"]["seeds"]["off"]
    expected = ensemble_inputs(request, candidate, inputs["num_sample"],
                               backend=SimpleNamespace(provenance=provenance), seed=seed)
    if inputs != expected or path.stem != ensemble_key(expected):
        return None
    return expected


class EnsembleCancelled(Exception):
    pass


def ensemble_inputs(request, candidate, num_sample: int = 1000, *, backend=None, seed: int = 0) -> dict:
    """Identify a saved legacy ensemble or a new explicitly selected native one."""
    if backend is not None:
        from .thermodynamics.ensemble_analysis import native_ensemble_inputs
        return native_ensemble_inputs(request, candidate, backend, num_sample=num_sample, seed=seed)
    if isinstance(num_sample, bool) or not isinstance(num_sample, Integral) or not 100 <= num_sample <= 5000:
        raise ValueError("Expected 100–5000 structures per state")
    for sequence in (candidate.trigger_seq, candidate.switch_seq):
        if not sequence or set(sequence) - set("ACGTU"):
            raise ValueError("Structure sampling requires concrete nucleotide sequences")
    inputs = {
        "schema": SCHEMA, "num_sample": int(num_sample),
        "trigger": candidate.trigger_seq, "switch": candidate.switch_seq,
        "architecture": asdict(request.architecture), "thermo": asdict(request.thermo),
        "reporter_length": len(request.reporter.sequence_after_start_rna),
    }
    inputs["functional_scheme"] = "functional_v2"
    return inputs


def _metrics_for_inputs(inputs, off, on):
    return functional_defect_metrics(SwitchArchitecture(**inputs["architecture"]), inputs["reporter_length"], off, on)


def ensemble_key(inputs: dict) -> str:
    return hashlib.sha256(json.dumps(inputs, sort_keys=True, allow_nan=False).encode()).hexdigest()


def analyze_structure_ensemble(request, candidate, import_path: str | None = None, *, num_sample: int = 1000,
                               cancel_check=None, progress_callback=None, backend=None, seed: int = 0,
                               thermodynamics_library_path: str | None = None) -> dict:
    """Calculate new structure families with an identified native backend."""
    from .thermodynamics.ensemble_analysis import analyze_structure_ensemble as analyze_native
    from .thermodynamics.setup import backend_for_request

    if backend is None:
        backend = backend_for_request(request, library_path=thermodynamics_library_path)
    return analyze_native(backend, request, candidate, num_sample=num_sample, seed=seed,
                          cancel_check=cancel_check, progress_callback=progress_callback)


def assemble_ensemble(inputs: dict, matrices: dict, samples: dict, provenance: dict) -> dict:
    """Shared display data from complete matrices and counted joint structures."""
    matrices = {state: _pair_matrix(matrices[state], len(inputs["switch"]) +
                                   (len(inputs["trigger"]) if state == "on" else 0), state)
                for state in ("off", "on")}
    for counts in samples.values():
        if (any(type(count) is not int or count <= 0 for count in counts.values())
                or sum(counts.values()) != inputs["num_sample"]):
            raise ValueError("Incomplete structure sample")
    # Reuse the exact same region definitions as the selection metric, without
    # applying its weights: grouping describes structures, not their score.
    regional = _metrics_for_inputs(inputs, matrices["off"], matrices["on"])
    states = {}
    for state in ("off", "on"):
        lengths = ([len(inputs["trigger"])] if state == "on" else []) + [len(inputs["switch"])]
        for structure in samples[state]:
            _pair_targets(structure, lengths)
        states[state] = {
            "strand_lengths": lengths, "switch_offset": lengths[0] if state == "on" else 0,
            "target": regional["target_structures"][state],
            "regions": {name: region["indices_switch"] for name, region in regional["regions"].items()
                        if name.startswith(state + "_")},
            "pairs": matrices[state].tolist(),
            "structures": [{"structure": structure, "count": count}
                           for structure, count in sorted(samples[state].items(), key=lambda item: (-item[1], item[0]))],
        }
    return {"schema": SCHEMA, "key": ensemble_key(inputs), "inputs": inputs,
            "created_at": datetime.now(timezone.utc).isoformat(), "model": dict(provenance), "states": states}


def cluster_structures(state: dict, *, tolerance: float = 0.05) -> dict:
    """Greedy, deterministic clustering by pairing partners on functional bases.

    Seeds are the most frequent remaining functional signatures (lexical ties).
    Each member differs from its seed on at most floor(tolerance * N) selected
    nucleotides. This is not a kinetic basin or a pairwise-diameter guarantee.
    Differences elsewhere may be arbitrary. All sampled counts are retained.
    """
    if not math.isfinite(tolerance) or not 0 <= tolerance <= 0.2:
        raise ValueError("Family tolerance must be between 0 and 0.2")
    lengths = state["strand_lengths"]
    offset = state["switch_offset"]
    regions = {name: np.asarray(indices, dtype=int) + offset for name, indices in state["regions"].items()}
    selected = sorted({int(i) for indices in regions.values() for i in indices})
    size = sum(lengths)
    if not selected or selected[0] < offset or selected[-1] >= size:
        raise ValueError("Invalid functional regions")
    target = _pair_targets(state["target"], lengths)
    signatures = {}
    total = 0
    for row in state["structures"]:
        count = row["count"]
        if type(count) is not int or count <= 0:
            raise ValueError("Sample counts must be positive integers")
        partners = _pair_targets(row["structure"], lengths)
        signature = tuple(int(i) for i in partners[selected])
        group = signatures.setdefault(signature, {"count": 0, "representative": row["structure"],
            "representative_count": 0, "errors": {name: 0. for name in regions}, "distinct_structures": 0})
        group["count"] += count
        group["distinct_structures"] += 1
        for name, indices in regions.items():
            group["errors"][name] += count * float(np.mean(partners[indices] != target[indices]))
        if (-count, row["structure"]) < (-group["representative_count"], group["representative"]):
            group["representative"], group["representative_count"] = row["structure"], count
        total += count
    if not total:
        raise ValueError("Empty structure sample")
    ordered = sorted(signatures.items(), key=lambda item: (-item[1]["count"], item[0]))
    vectors = np.asarray([signature for signature, _ in ordered], dtype=int)
    assigned = np.zeros(len(ordered), dtype=bool)
    max_changes = math.floor(tolerance * len(selected) + 1e-10)
    families = []
    for seed, (_, seed_group) in enumerate(ordered):
        if assigned[seed]:
            continue
        members = np.flatnonzero(~assigned & (np.count_nonzero(vectors != vectors[seed], axis=1) <= max_changes))
        assigned[members] = True
        count = sum(ordered[i][1]["count"] for i in members)
        families.append({
            "representative": seed_group["representative"], "count": count, "fraction": count / total,
            "signature_count": len(members),
            "distinct_structures": sum(ordered[i][1]["distinct_structures"] for i in members),
            "regional_defects": {name: math.fsum(ordered[i][1]["errors"][name] for i in members) / count
                                 for name in regions},
        })
    families.sort(key=lambda family: (-family["count"], family["representative"]))
    return {"families": families, "sample_count": total, "selected_count": len(selected),
            "max_changes": max_changes, "tolerance": tolerance}


def probability_helices(state: dict, *, threshold: float = 0.2) -> list[dict]:
    """Maximal stacks above a marginal threshold, with strand breaks respected.

    A row is a possible local motif, not evidence that all its pairs coexist.
    Competing rows share a nucleotide with a different partner. No product of
    probabilities is returned as a helix population.
    """
    if not math.isfinite(threshold) or not 0 < threshold <= 1:
        raise ValueError("Pair threshold must be in (0, 1]")
    lengths = state["strand_lengths"]
    matrix = _pair_matrix(state["pairs"], sum(lengths), "ensemble")
    chains = np.repeat(np.arange(len(lengths)), lengths)
    pairs = {(int(i), int(j)) for i, j in zip(*np.where(np.triu(matrix, 1) >= threshold))}
    helices = []
    for i, j in sorted(pairs):
        if (i-1, j+1) in pairs and chains[i-1] == chains[i] and chains[j+1] == chains[j]:
            continue
        stack = [(i, j)]
        a, b = i+1, j-1
        while (a, b) in pairs and chains[a] == chains[i] and chains[b] == chains[j]:
            stack.append((a, b))
            a, b = a+1, b-1
        values = [float(matrix[a, b]) for a, b in stack]
        helices.append({"pairs": stack, "length": len(stack), "mean_pair_probability": math.fsum(values) / len(values),
                        "min_pair_probability": min(values)})
    helices.sort(key=lambda h: (-h["mean_pair_probability"], -h["length"], h["pairs"]))
    endpoints = [dict([*h["pairs"], *((j, i) for i, j in h["pairs"])]) for h in helices]
    for index, helix in enumerate(helices):
        helix["competing_motifs"] = [other+1 for other, mapping in enumerate(endpoints)
                                    if any(i in mapping and mapping[i] != j for i, j in endpoints[index].items())]
    return helices


def ensemble_path(output_dir: Path, key: str) -> Path:
    if len(key) != 64 or set(key) - set("0123456789abcdef"):
        raise ValueError("Invalid ensemble key")
    return Path(output_dir) / "structure_ensembles" / f"{key}.json"


def save_ensemble(output_dir: Path, analysis: dict) -> Path:
    path = ensemble_path(output_dir, analysis["key"])
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary = tempfile.mkstemp(prefix=".ensemble-", suffix=".tmp", dir=path.parent)
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as stream:
            json.dump(analysis, stream, ensure_ascii=False, allow_nan=False, separators=(",", ":"))
        os.replace(temporary, path)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)
    return path


def load_ensemble(output_dir: Path, inputs: dict) -> dict | None:
    key = ensemble_key(inputs)
    path = ensemble_path(output_dir, key)
    if not path.exists():
        return None
    # Full-reporter ensembles contain two quadratic matrices and can exceed
    # the former fixed 32 MiB allowance without being corrupt. Budget for JSON
    # floats and sampled dot-brackets using the already selected input lengths;
    # retain an absolute ceiling for untrusted/imported files.
    off_length = len(inputs["switch"])
    on_length = off_length + len(inputs["trigger"])
    expected_bytes = (32 * (off_length ** 2 + on_length ** 2)
                      + 2 * inputs["num_sample"] * (off_length + on_length + 128))
    maximum_bytes = min(512 * 1024 * 1024, max(32 * 1024 * 1024, expected_bytes))
    if path.stat().st_size > maximum_bytes:
        raise ValueError("Saved ensemble is too large")
    analysis = json.loads(path.read_text(encoding="utf-8"))
    if (not isinstance(analysis, dict) or analysis.get("schema") != SCHEMA
            or analysis.get("key") != key or analysis.get("inputs") != inputs
            or not isinstance(analysis.get("model"), dict)):
        raise ValueError("Saved ensemble does not match this candidate")
    if "thermodynamics" in inputs and analysis["model"] != inputs["thermodynamics"]:
        raise ValueError("Saved ensemble has a different thermodynamic model")
    connected_sampling = inputs.get("sampling", {}).get("method") == "joint-connected-boltzmann"
    conventions = analysis["model"].get("ensemble", "")
    check_chemistry = (connected_sampling and isinstance(conventions, str)
                       and {"RNA-WC-GU", "DNA-WC", "hybrid-WC"}.issubset(conventions.split(";")))
    if connected_sampling:
        # These graph helpers have no energy/runtime dependency. A portable
        # archive must preserve the connected ensemble promised by its inputs.
        from .thermodynamics.structures import is_connected, parse_structure
    for name, lengths in (("off", [len(inputs["switch"])]),
                          ("on", [len(inputs["trigger"]), len(inputs["switch"])])):
        state = analysis["states"][name]
        if state["strand_lengths"] != lengths or state["switch_offset"] != (lengths[0] if name == "on" else 0):
            raise ValueError("Invalid saved strand order")
        _pair_matrix(state["pairs"], sum(lengths), name)
        if sum(row["count"] for row in state["structures"]) != inputs["num_sample"]:
            raise ValueError("Incomplete saved sample")
        if check_chemistry:
            sequences = ([inputs["trigger"]] if name == "on" else []) + [inputs["switch"]]
            strand_materials = (([inputs["thermo"]["trigger_material"]] if name == "on" else [])
                                + [inputs["thermo"]["switch_material"]])
            sequence = "".join(sequences).replace("U", "T")
            materials = [material for strand, material in zip(sequences, strand_materials) for _ in strand]
        # Validate every representative before a saved/imported file reaches
        # the viewer. Clustering is done separately with the chosen tolerance.
        _pair_targets(state["target"], lengths)
        for row in state["structures"]:
            if type(row["count"]) is not int or row["count"] <= 0:
                raise ValueError("Invalid saved sample count")
            partners = _pair_targets(row["structure"], lengths)
            if connected_sampling and not is_connected(parse_structure(row["structure"], lengths), lengths):
                raise ValueError("Saved native ensemble contains a disconnected complex")
            if check_chemistry:
                # Match the native model's declared pair mask, without loading
                # its DLL. Conceptual targets and probability matrices are not
                # sampled structures and are not subject to this check.
                for i, j in enumerate(partners):
                    if j <= i:
                        continue
                    pair = sequence[i] + sequence[j]
                    if pair not in ("AT", "TA", "CG", "GC") and not (
                            materials[i] == materials[j] == "rna" and pair in ("GT", "TG")):
                        raise ValueError("Saved native ensemble contains a chemically forbidden base pair")
        if len({row["structure"] for row in state["structures"]}) != len(state["structures"]):
            raise ValueError("Duplicate saved structures")
    expected = _metrics_for_inputs(inputs, analysis["states"]["off"]["pairs"], analysis["states"]["on"]["pairs"])
    for name in ("off", "on"):
        state = analysis["states"][name]
        regions = {key: value["indices_switch"] for key, value in expected["regions"].items() if key.startswith(name + "_")}
        if state["target"] != expected["target_structures"][name] or state["regions"] != regions:
            raise ValueError("Saved functional regions or targets do not match this candidate")
    return analysis
