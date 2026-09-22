"""Real NUPACK 5′, 3′ or two-ended trigger extension.

Numerical reference outputs from the original algorithm are preserved in
``tests/fixtures/extension_original_reference.json``.
For each requested final length we enumerate the complete sequence space, keep
the lexicographically best MFE preservation group, then rank that group using
the original ensemble probability loss and score. A survivor limit only limits
the ensemble stage; it never silently limits the exhaustive MFE search.

Each length is an independent execution of the source's fixed-length search.
This preserves candidates for every requested length. There is no fallback
sequence generator and no mock thermodynamics in this module.
"""

from __future__ import annotations

from .localization import tr

import itertools
import math
import time
from dataclasses import asdict, dataclass
from typing import Any, Callable, Sequence

from .nupack_engine import nupack_import_context
from .sequences import normalize_dna, normalize_rna

Pair = tuple[int, int]
ProgressCallback = Callable[[int, int, str], None]
EXTENSION_SIDES = {"5prime": "5′", "3prime": "3′", "both": "Moitié 5′ / moitié 3′"}


def extension_split(extension_length: int, extension_side: str = "5prime") -> tuple[int, int]:
    """Return (5′ bases, 3′ bases); an odd extra base goes to 5′ in split mode."""
    if extension_side not in EXTENSION_SIDES:
        raise ValueError(tr('Choisissez une extension en 5′, en 3′ ou moitié de chaque côté.'))
    if isinstance(extension_length, bool) or not isinstance(extension_length, int) or extension_length < 0:
        raise ValueError(tr('La longueur d’extension doit être un entier positif ou nul.'))
    if extension_side == "5prime":
        return extension_length, 0
    if extension_side == "3prime":
        return 0, extension_length
    return (extension_length + 1) // 2, extension_length // 2


@dataclass(frozen=True)
class MfeResult:
    structure: str
    energy_kcal_mol: float


@dataclass
class ComplexAnalysisResult:
    mfe_structure: str
    mfe_energy_kcal_mol: float
    ensemble_free_energy_kcal_mol: float
    mfe_structure_probability: float
    target_structure_probability: float
    pair_matrix: Any


@dataclass(frozen=True)
class MfeCandidate:
    extension_5prime: str
    extension_length: int
    trigger_length: int
    extended_trigger: str
    gc_fraction: float
    longest_homopolymer: int
    mfe_energy_kcal_mol: float
    mfe_energy_delta: float
    mfe_core_hamming_distance: int
    mfe_padded_hamming_distance: int
    missing_reference_mfe_pairs: int
    new_core_mfe_pairs: int
    extension_mfe_pairs: int
    extension_unpaired_fraction_mfe: float
    mfe_structure: str
    extension_3prime: str = ""


def estimate_extension_search(trigger_length: int, target_lengths: Sequence[int]) -> dict:
    """Validate lengths and report the exact, uncapped number of sequences."""
    if isinstance(trigger_length, bool) or int(trigger_length) != trigger_length or trigger_length <= 0:
        raise ValueError(tr('La longueur du trigger original doit être un entier positif.'))
    lengths = list(target_lengths)
    if not lengths:
        raise ValueError(tr('Choisissez au moins une longueur finale de trigger.'))
    if any(isinstance(value, bool) or not isinstance(value, int) for value in lengths):
        raise ValueError(tr('Les longueurs finales doivent être des entiers.'))
    lengths = sorted(set(lengths))
    if lengths[0] <= trigger_length:
        raise ValueError(tr('Chaque longueur finale doit dépasser celle du trigger original.'))
    by_length = [
        {
            "target_length": length,
            "extension_length": length - trigger_length,
            "total_candidates_expected": 4 ** (length - trigger_length),
        }
        for length in lengths
    ]
    return {
        "target_lengths": lengths,
        "extension_lengths": [row["extension_length"] for row in by_length],
        "total_candidates_expected": sum(row["total_candidates_expected"] for row in by_length),
        "by_length": by_length,
    }


def parse_pairs(structure: str) -> set[Pair]:
    open_to_close = {"(": ")", "[": "]", "{": "}", "<": ">"}
    close_to_open = {close: open_ for open_, close in open_to_close.items()}
    stacks: dict[str, list[int]] = {open_: [] for open_ in open_to_close}
    pairs: set[Pair] = set()
    index = 0
    for char in structure:
        if char == "+":
            continue
        if char in stacks:
            stacks[char].append(index)
        elif char in close_to_open:
            stack = stacks[close_to_open[char]]
            if not stack:
                raise ValueError(tr('Structure dot-bracket déséquilibrée : {v0}', v0=structure))
            pairs.add((stack.pop(), index))
        elif char != ".":
            raise ValueError(tr('Caractère de structure non pris en charge : {v0!r}', v0=char))
        index += 1
    if any(stacks.values()):
        raise ValueError(tr('Structure dot-bracket déséquilibrée : {v0}', v0=structure))
    return pairs


def padded_reference_structure(reference_structure: str, extension_length: int, extension_side: str = "5prime") -> str:
    parts = reference_structure.split("+")
    if len(parts) != 2:
        raise ValueError(tr('Une structure de référence à deux brins est attendue'))
    prefix, suffix = extension_split(extension_length, extension_side)
    return parts[0] + "+" + "." * prefix + parts[1] + "." * suffix


def map_reference_position_to_candidate(position: int, aptamer_length: int, extension_length: int) -> int:
    """Shift trigger positions only by the number of bases added in 5′."""
    return position if position < aptamer_length else position + extension_length


def map_reference_pairs_to_extended(pairs: set[Pair], aptamer_length: int, extension_length: int) -> set[Pair]:
    return {
        (
            map_reference_position_to_candidate(i, aptamer_length, extension_length),
            map_reference_position_to_candidate(j, aptamer_length, extension_length),
        )
        for i, j in pairs
    }


def mfe_key(candidate: MfeCandidate) -> tuple[int, ...]:
    return (
        candidate.mfe_core_hamming_distance,
        candidate.missing_reference_mfe_pairs,
        candidate.extension_mfe_pairs,
        candidate.new_core_mfe_pairs,
        candidate.mfe_padded_hamming_distance,
    )


def score_mfe_candidate(
    extension: str,
    candidate_mfe: MfeResult,
    reference_mfe: MfeResult,
    reference_pairs_extended: set[Pair],
    aptamer_length: int,
    trigger_sequence: str,
    *, extension_side: str = "5prime",
) -> MfeCandidate:
    """The five source MFE criteria, without system-specific constants."""
    extension_length = len(extension)
    prefix_length, suffix_length = extension_split(extension_length, extension_side)
    prefix, suffix = extension[:prefix_length], extension[prefix_length:]
    candidate_pairs = parse_pairs(candidate_mfe.structure)
    reference_padded = padded_reference_structure(reference_mfe.structure, extension_length, extension_side)
    parts = candidate_mfe.structure.split("+")
    if len(parts) != 2 or len(parts[0]) != aptamer_length or len(parts[1]) != extension_length + len(trigger_sequence):
        raise ValueError(tr('NUPACK a renvoyé une structure dont les longueurs de brins sont inattendues.'))
    core_structure = parts[0] + "+" + parts[1][prefix_length:prefix_length + len(trigger_sequence)]

    def is_extension(position: int) -> bool:
        return (aptamer_length <= position < aptamer_length + prefix_length
                or aptamer_length + prefix_length + len(trigger_sequence) <= position < aptamer_length + len(trigger_sequence) + extension_length)

    extension_pairs = sum(is_extension(i) or is_extension(j) for i, j in candidate_pairs)
    new_core_pairs = sum(
        not is_extension(i) and not is_extension(j)
        for i, j in candidate_pairs - reference_pairs_extended
    )
    longest = max((len(list(group)) for segment in (prefix, suffix) for _, group in itertools.groupby(segment)), default=0)
    return MfeCandidate(
        extension_5prime=prefix,
        extension_3prime=suffix,
        extension_length=extension_length,
        trigger_length=extension_length + len(trigger_sequence),
        extended_trigger=prefix + trigger_sequence + suffix,
        gc_fraction=round(sum(base in "GC" for base in extension) / extension_length, 6),
        longest_homopolymer=longest,
        mfe_energy_kcal_mol=round(candidate_mfe.energy_kcal_mol, 6),
        mfe_energy_delta=round(candidate_mfe.energy_kcal_mol - reference_mfe.energy_kcal_mol, 6),
        mfe_core_hamming_distance=sum(a != b for a, b in zip(reference_mfe.structure, core_structure)),
        mfe_padded_hamming_distance=sum(a != b for a, b in zip(reference_padded, candidate_mfe.structure)),
        missing_reference_mfe_pairs=len(reference_pairs_extended - candidate_pairs),
        new_core_mfe_pairs=new_core_pairs,
        extension_mfe_pairs=extension_pairs,
        # Keep the source's pair-count definition for numerical compatibility.
        extension_unpaired_fraction_mfe=round(1.0 - extension_pairs / extension_length, 6),
        mfe_structure=candidate_mfe.structure,
    )


def final_score(
    reference_pair_probability_loss: float,
    core_pair_probability_rmse: float,
    core_pair_probability_mae: float,
    extension_pair_probability_sum: float,
    extension_to_aptamer_probability_sum: float,
    extension_to_trigger_core_probability_sum: float,
    target_structure_probability: float,
    candidate: MfeCandidate,
) -> float:
    return round(
        1000.0 * reference_pair_probability_loss
        + 1000.0 * core_pair_probability_rmse
        + 500.0 * core_pair_probability_mae
        + 100.0 * extension_pair_probability_sum
        + 150.0 * extension_to_aptamer_probability_sum
        + 100.0 * extension_to_trigger_core_probability_sum
        + 10.0 * (1.0 - target_structure_probability)
        + 100.0 * candidate.mfe_core_hamming_distance
        + 100.0 * candidate.missing_reference_mfe_pairs
        + 25.0 * candidate.new_core_mfe_pairs
        + 25.0 * candidate.extension_mfe_pairs,
        9,
    )


def final_rank_key(candidate: dict) -> tuple:
    return (
        candidate["final_score"],
        candidate["reference_pair_probability_loss"],
        candidate["core_pair_probability_rmse"],
        candidate["core_pair_probability_mae"],
        candidate["extension_pair_probability_sum"],
        candidate["extension_to_aptamer_probability_sum"],
        candidate["extension_to_trigger_core_probability_sum"],
        -candidate["target_structure_probability"],
        candidate["mfe_core_hamming_distance"],
        candidate["missing_reference_mfe_pairs"],
        candidate["new_core_mfe_pairs"],
        candidate["extension_mfe_pairs"],
        abs(candidate["ensemble_free_energy_delta"]),
        abs(candidate["mfe_energy_delta"]),
        candidate["trigger_length"],
        candidate["extension_5prime"] + candidate.get("extension_3prime", ""),
    )


def _prefixed(sequence: str, chemistry: str) -> str:
    return ("r" if chemistry == "rna" else "d") + sequence


def _analyze_complex(nupack, model, aptamer: str, trigger: str, aptamer_chemistry: str,
                     trigger_chemistry: str, name: str, target: str | None = None) -> ComplexAnalysisResult:
    strands = [_prefixed(aptamer, aptamer_chemistry), _prefixed(trigger, trigger_chemistry)]
    complex_obj = nupack.Complex(
        [nupack.Strand(strands[0], name="aptamer"), nupack.Strand(strands[1], name="trigger")],
        name=name,
    )
    result = nupack.complex_analysis(complexes=[complex_obj], model=model, compute=["mfe", "pairs", "pfunc"])[complex_obj]
    if not result.mfe:
        raise RuntimeError(tr('NUPACK n’a renvoyé aucune structure MFE.'))
    mfe = result.mfe[0]
    structure = str(mfe.structure)
    return ComplexAnalysisResult(
        mfe_structure=structure,
        mfe_energy_kcal_mol=float(mfe.energy),
        ensemble_free_energy_kcal_mol=float(result.free_energy),
        mfe_structure_probability=float(nupack.structure_probability(strands=strands, structure=structure, model=model)),
        target_structure_probability=float(nupack.structure_probability(strands=strands, structure=target or structure, model=model)),
        pair_matrix=result.pairs.to_array(),
    )


def _nucleotide_probabilities(analysis: ComplexAnalysisResult) -> list[float]:
    partners = {position: other for i, j in parse_pairs(analysis.mfe_structure) for position, other in [(i, j), (j, i)]}
    return [float(analysis.pair_matrix[i, partners.get(i, i)]) for i in range(len(analysis.pair_matrix))]


def _reference_pair_rows(reference: ComplexAnalysisResult, aptamer_length: int,
                         candidate: ComplexAnalysisResult | None = None, extension_length: int = 0) -> list[dict]:
    """Compare reference positions shifted by the 5′ prefix length only."""
    rows = []
    for i, j in sorted(parse_pairs(reference.mfe_structure)):
        cross = i < aptamer_length <= j
        row = {
            "reference_position_a": i + 1,
            "reference_position_b": j + 1,
            "aptamer_index": i + 1 if cross else None,
            "trigger_index": j - aptamer_length + 1 if cross else None,
            "reference_probability": round(float(reference.pair_matrix[i, j]), 12),
        }
        if candidate is not None:
            a = map_reference_position_to_candidate(i, aptamer_length, extension_length)
            b = map_reference_position_to_candidate(j, aptamer_length, extension_length)
            row["candidate_probability"] = round(float(candidate.pair_matrix[a, b]), 12)
        rows.append(row)
    return rows


def _extension_positions(candidate: MfeCandidate, aptamer_length: int, trigger_length: int) -> list[int]:
    core_start = aptamer_length + len(candidate.extension_5prime)
    core_end = core_start + trigger_length
    return list(range(aptamer_length, core_start)) + list(range(core_end, core_end + len(candidate.extension_3prime)))


def ensemble_metrics(candidate: MfeCandidate, reference: ComplexAnalysisResult,
                     analysis: ComplexAnalysisResult, aptamer_length: int, trigger_length: int) -> dict:
    """Source ensemble score; diagonal unpaired probabilities are included."""
    import numpy as np

    prefix_length = len(candidate.extension_5prime)
    core_start = aptamer_length + prefix_length
    core_indices = list(range(aptamer_length)) + list(range(core_start, core_start + trigger_length))
    difference = analysis.pair_matrix[np.ix_(core_indices, core_indices)] - reference.pair_matrix
    rmse = float(np.sqrt(np.mean(difference * difference)))
    mae = float(np.mean(np.abs(difference)))
    loss = sum(
        max(0.0, float(reference.pair_matrix[i, j]) - float(analysis.pair_matrix[
            map_reference_position_to_candidate(i, aptamer_length, prefix_length),
            map_reference_position_to_candidate(j, aptamer_length, prefix_length),
        ]))
        for i, j in parse_pairs(reference.mfe_structure)
    )
    positions = _extension_positions(candidate, aptamer_length, trigger_length)
    paired = sum(1.0 - float(analysis.pair_matrix[pos, pos]) for pos in positions)
    to_aptamer = float(analysis.pair_matrix[np.ix_(positions, range(aptamer_length))].sum())
    to_core = float(analysis.pair_matrix[np.ix_(positions, range(core_start, core_start + trigger_length))].sum())
    result = asdict(candidate)
    result.update(
        rank=0,
        final_score=final_score(loss, rmse, mae, paired, to_aptamer, to_core, analysis.target_structure_probability, candidate),
        mfe_structure=analysis.mfe_structure,
        mfe_energy_kcal_mol=round(analysis.mfe_energy_kcal_mol, 9),
        mfe_energy_delta=round(analysis.mfe_energy_kcal_mol - reference.mfe_energy_kcal_mol, 9),
        ensemble_free_energy_kcal_mol=round(analysis.ensemble_free_energy_kcal_mol, 9),
        ensemble_free_energy_delta=round(analysis.ensemble_free_energy_kcal_mol - reference.ensemble_free_energy_kcal_mol, 9),
        mfe_structure_probability=round(analysis.mfe_structure_probability, 12),
        target_structure_probability=round(analysis.target_structure_probability, 12),
        core_pair_probability_rmse=round(rmse, 12),
        core_pair_probability_mae=round(mae, 12),
        reference_pair_probability_loss=round(loss, 12),
        extension_pair_probability_sum=round(paired, 12),
        extension_to_aptamer_probability_sum=round(to_aptamer, 12),
        extension_to_trigger_core_probability_sum=round(to_core, 12),
    )
    return result


def design_trigger_extensions(
    *, aptamer_dna: str, trigger_dna: str, target_lengths: Sequence[int], thermo,
    import_path: str, top: int = 100, max_ensemble_candidates: int = 50000,
    extension_side: str = "5prime",
    progress_callback: ProgressCallback | None = None,
    cancel_check: Callable[[], bool] | None = None,
) -> dict:
    """Run exhaustive MFE searches and ensemble ranking with real NUPACK.

    ``top`` and ``max_ensemble_candidates`` apply separately to each final
    length. Cancellation is observed between NUPACK calls; completed ensemble
    analyses can be returned, always labelled as cancelled/partial.
    """
    if not import_path:
        raise RuntimeError(tr('NUPACK est obligatoire pour calculer une extension réelle. Configurez son chemin.'))
    extension_split(0, extension_side)
    if any(isinstance(value, bool) or not isinstance(value, int) or value <= 0 for value in (top, max_ensemble_candidates)):
        raise ValueError(tr("Le nombre exporté et la limite de candidats d'ensemble doivent être des entiers positifs."))
    material = str(thermo.material)
    homogeneous_rna = material.startswith("rna") and "-" not in material
    aptamer_chemistry = "rna" if homogeneous_rna else "dna"
    trigger_chemistry = "rna" if homogeneous_rna else (str(getattr(thermo, "trigger_material", "dna")) if "-" in material else "dna")
    if trigger_chemistry not in {"dna", "rna"}:
        raise ValueError(tr('La chimie du trigger doit être DNA ou RNA.'))
    aptamer = (normalize_rna if aptamer_chemistry == "rna" else normalize_dna)(aptamer_dna, field_name="aptamer")
    trigger = (normalize_rna if trigger_chemistry == "rna" else normalize_dna)(trigger_dna, field_name="trigger")
    estimate = estimate_extension_search(len(trigger), target_lengths)
    model_settings = {
        "material": material,
        "temperature_c": float(thermo.temperature_c),
        "sodium_m": float(thermo.sodium_m),
        "magnesium_m": float(thermo.magnesium_m),
        "ensemble": "stacking",
        "aptamer_material": aptamer_chemistry,
        "trigger_material": trigger_chemistry,
    }
    if not all(math.isfinite(model_settings[key]) for key in ("temperature_c", "sodium_m", "magnesium_m")):
        raise ValueError(tr('Les conditions thermodynamiques doivent être finies.'))
    if model_settings["sodium_m"] < 0 or model_settings["magnesium_m"] < 0:
        raise ValueError(tr('Les concentrations ioniques doivent être positives ou nulles.'))
    expected = estimate["total_candidates_expected"]
    started = time.monotonic()
    warnings: list[str] = []
    if max(estimate["target_lengths"]) > 30:
        warnings.append(tr('Au-delà de 30 nt, le calcul peut prendre beaucoup de temps ; chaque base ajoutée multiplie la recherche par quatre.'))
    result: dict = {
        "status": "running",
        "input": {
            "aptamer_dna": aptamer, "trigger_dna": trigger,
            "target_lengths": estimate["target_lengths"],
            "extension_lengths": estimate["extension_lengths"],
            "mode": "exhaustive_mfe_then_pair_probability_ranking",
            "selection_scope": "independent_search_per_target_length",
            "top_per_length": top,
            "extension_side": extension_side,
            "odd_split_extra_base": "5prime" if extension_side == "both" else None,
        },
        "model": model_settings,
        "reference": {},
        "mfe_screen": {
            "total_candidates_expected": expected, "total_candidates_screened": 0,
            "mfe_survivor_cap": max_ensemble_candidates, "mfe_survivor_cap_hit": False,
            "exhaustive": False, "workers": 1, "by_length": [],
        },
        "candidates": [], "curves": {"mfe_progress": [], "ensemble_progress": [], "ranking": []},
        "warnings": warnings,
    }
    cancelled = lambda: bool(cancel_check and cancel_check())

    def progress(done: int, total: int, message: str) -> None:
        if progress_callback:
            progress_callback(done, total, message)

    progress(0, expected, tr('NUPACK · extension {v0} : {v1:,} extensions à examiner, recherche exhaustive.', v0=tr(EXTENSION_SIDES[extension_side]), v1=expected))
    if cancelled():
        result["status"] = "cancelled"
        return result
    try:
        with nupack_import_context(import_path) as nupack:
            model = nupack.Model(material=material, celsius=model_settings["temperature_c"],
                                 sodium=model_settings["sodium_m"], magnesium=model_settings["magnesium_m"], ensemble="stacking")
            progress(0, expected, tr('NUPACK · référence aptamère + trigger original : MFE, ensemble et probabilités de paires…'))
            reference = _analyze_complex(nupack, model, aptamer, trigger, aptamer_chemistry, trigger_chemistry, "reference")
            reference_mfe = MfeResult(reference.mfe_structure, reference.mfe_energy_kcal_mol)
            reference_pairs = parse_pairs(reference.mfe_structure)
            if not any(i < len(aptamer) <= j for i, j in reference_pairs):
                warnings.append(tr("La structure MFE de référence ne présente aucune paire entre l'aptamère et le trigger ; la conservation de cette référence ne démontre pas leur liaison."))
            result["reference"] = {
                "aptamer_sequence": aptamer, "trigger_sequence": trigger,
                "mfe_structure": reference.mfe_structure,
                "mfe_energy_kcal_mol": reference.mfe_energy_kcal_mol,
                "ensemble_free_energy_kcal_mol": reference.ensemble_free_energy_kcal_mol,
                "mfe_structure_probability": reference.mfe_structure_probability,
                "mfe_pairs": [list(pair) for pair in sorted(reference_pairs)],
                "pair_probabilities": _reference_pair_rows(reference, len(aptamer)),
                "nucleotide_probabilities": _nucleotide_probabilities(reference),
            }
            progress(0, expected, tr('Référence calculée : MFE {v0:.6f} kcal/mol ; début du criblage MFE.', v0=reference.mfe_energy_kcal_mol))
            survivors: list[MfeCandidate] = []
            seen = 0
            last_report = time.monotonic()
            screen_started = time.monotonic()
            was_cancelled = cancelled()
            for plan in estimate["by_length"]:
                if was_cancelled:
                    break
                extension_length = plan["extension_length"]
                prefix_length, suffix_length = extension_split(extension_length, extension_side)
                plan.update(extension_5prime_length=prefix_length, extension_3prime_length=suffix_length)
                progress(seen, expected, tr('Cible {v0} nt · ajout de {v1} bases en 5′ et {v2} bases en 3′.', v0=plan['target_length'], v1=prefix_length, v2=suffix_length))
                mapped_pairs = map_reference_pairs_to_extended(reference_pairs, len(aptamer), prefix_length)
                best_key = None
                best_count = 0
                group: list[MfeCandidate] = []
                length_seen = 0
                for bases in itertools.product("ACGU" if trigger_chemistry == "rna" else "ACGT", repeat=extension_length):
                    if cancelled():
                        was_cancelled = True
                        break
                    extension = "".join(bases)
                    extended_trigger = extension[:prefix_length] + trigger + extension[prefix_length:]
                    records = nupack.mfe(strands=[_prefixed(aptamer, aptamer_chemistry), _prefixed(extended_trigger, trigger_chemistry)], model=model)
                    if not records:
                        raise RuntimeError(tr('NUPACK n’a renvoyé aucune structure MFE.'))
                    candidate = score_mfe_candidate(extension, MfeResult(str(records[0].structure), float(records[0].energy)),
                                                    reference_mfe, mapped_pairs, len(aptamer), trigger, extension_side=extension_side)
                    key = mfe_key(candidate)
                    if best_key is None or key < best_key:
                        best_key, best_count, group = key, 1, [candidate]
                    elif key == best_key:
                        best_count += 1
                        if len(group) < max_ensemble_candidates:
                            group.append(candidate)
                    seen += 1
                    length_seen += 1
                    now = time.monotonic()
                    if seen == 1 or now - last_report >= 1.0 or length_seen == plan["total_candidates_expected"]:
                        last_report = now
                        progress(seen, expected, tr('MFE · cible {v0} nt : {v1:,}/{v2:,} ; meilleur groupe {v3:,}, clé {v4}.', v0=plan['target_length'], v1=length_seen, v2=plan['total_candidates_expected'], v3=best_count, v4=best_key))
                        result["curves"]["mfe_progress"].append({"screened": seen, "target_length": plan["target_length"], "best_group_count": best_count, "best_key": list(best_key)})
                cap_hit = best_count > len(group)
                stats = {**plan, "total_candidates_screened": length_seen,
                         "mfe_best_key": list(best_key) if best_key is not None else None,
                         "mfe_best_group_count": best_count, "mfe_survivors_for_ensemble": len(group),
                         "mfe_survivor_cap": max_ensemble_candidates, "mfe_survivor_cap_hit": cap_hit,
                         "exhaustive": length_seen == plan["total_candidates_expected"]}
                result["mfe_screen"]["by_length"].append(stats)
                if best_key is not None and any(best_key):
                    warnings.append(tr('Cible {v0} nt : aucune extension examinée ne conserve exactement la structure MFE de référence avec une extension non appariée (meilleure clé {v1}). Les candidats présentés sont les meilleurs compromis calculés.', v0=plan['target_length'], v1=best_key))
                if cap_hit:
                    warnings.append(tr("Cible {v0} nt : {v1:,} ex æquo MFE ; seuls les {v2:,} premiers en ordre alphabétique passent à l'ensemble. Le classement d'ensemble est partiel.", v0=plan['target_length'], v1=best_count, v2=max_ensemble_candidates))
                survivors.extend(group)
            screen = result["mfe_screen"]
            screen.update(total_candidates_screened=seen, exhaustive=seen == expected,
                          mfe_best_key=screen["by_length"][0]["mfe_best_key"] if len(screen["by_length"]) == 1 else None,
                          mfe_best_group_count=sum(row["mfe_best_group_count"] for row in screen["by_length"]),
                          mfe_survivors_for_ensemble=len(survivors),
                          mfe_survivor_cap_hit=any(row["mfe_survivor_cap_hit"] for row in screen["by_length"]),
                          mfe_screen_seconds=round(time.monotonic() - screen_started, 3))
            ensemble_started = time.monotonic()
            final_candidates: list[dict] = []
            # Keep detailed plots only for the top per length, bounded by top.
            details: dict[str, dict] = {}
            best_by_length: dict[int, list[dict]] = {}
            last_ensemble_report = time.monotonic()
            for index, candidate in enumerate(survivors, start=1):
                if was_cancelled or cancelled():
                    was_cancelled = True
                    break
                if index == 1 or index == len(survivors) or time.monotonic() - last_ensemble_report >= 1.0:
                    last_ensemble_report = time.monotonic()
                    progress(index - 1, len(survivors), tr('Ensemble NUPACK · candidat {v0:,}/{v1:,} ; 5′ {v2} ; 3′ {v3}…', v0=index, v1=len(survivors), v2=candidate.extension_5prime or '—', v3=candidate.extension_3prime or '—'))
                analysis = _analyze_complex(nupack, model, aptamer, candidate.extended_trigger, aptamer_chemistry,
                                            trigger_chemistry, f"extension_{candidate.extended_trigger}",
                                            padded_reference_structure(reference.mfe_structure, candidate.extension_length, extension_side))
                row = ensemble_metrics(candidate, reference, analysis, len(aptamer), len(trigger))
                row.update(extension_side=extension_side, extension_5prime_length=len(candidate.extension_5prime),
                           extension_3prime_length=len(candidate.extension_3prime))
                final_candidates.append(row)
                best_for_length = sorted(best_by_length.get(candidate.trigger_length, []) + [row], key=final_rank_key)[:top]
                best_by_length[candidate.trigger_length] = best_for_length
                if row in best_for_length:
                    details[candidate.extended_trigger] = {
                        "reference_pairs": _reference_pair_rows(reference, len(aptamer), analysis, len(candidate.extension_5prime)),
                        "nucleotide_probabilities": _nucleotide_probabilities(analysis),
                        "extension_base_probabilities": [
                            {"extension_index": ext_index, "base": candidate.extended_trigger[pos - len(aptamer)],
                             "side": "5prime" if pos < len(aptamer) + len(candidate.extension_5prime) else "3prime",
                             "trigger_index": pos - len(aptamer) + 1,
                             "paired_probability": round(1.0 - float(analysis.pair_matrix[pos, pos]), 12)}
                            for ext_index, pos in enumerate(_extension_positions(candidate, len(aptamer), len(trigger)), start=1)
                        ],
                    }
                    keep = {item["extended_trigger"] for item in best_for_length}
                    for sequence in list(details):
                        if len(sequence) == candidate.trigger_length and sequence not in keep:
                            del details[sequence]
                result["curves"]["ensemble_progress"].append({"analyzed": index, "trigger_length": candidate.trigger_length,
                                                              "final_score": row["final_score"], "core_pair_probability_rmse": row["core_pair_probability_rmse"],
                                                              "reference_pair_probability_loss": row["reference_pair_probability_loss"],
                                                              "extension_pair_probability_sum": row["extension_pair_probability_sum"]})
            counts: dict[int, int] = {}
            ranked = []
            for row in sorted(final_candidates, key=final_rank_key):
                length = row["trigger_length"]
                counts[length] = counts.get(length, 0) + 1
                if counts[length] <= top:
                    row.update(rank=len(ranked) + 1, rank_within_length=counts[length], **details[row["extended_trigger"]])
                    ranked.append(row)
            result["candidates"] = ranked
            result["curves"]["ranking"] = [
                {key: row[key] for key in ("rank", "trigger_length", "extension_5prime", "extension_3prime", "final_score", "core_pair_probability_rmse", "reference_pair_probability_loss", "extension_pair_probability_sum")}
                for row in ranked
            ]
            screen.update(ensemble_candidates_analyzed=len(final_candidates),
                          ensemble_pass_seconds=round(time.monotonic() - ensemble_started, 3),
                          ensemble_complete=len(final_candidates) == len(survivors) and not was_cancelled,
                          ensemble_ranking_exhaustive=not screen["mfe_survivor_cap_hit"] and len(final_candidates) == len(survivors) and screen["exhaustive"] and not was_cancelled)
            result["status"] = "cancelled" if was_cancelled else "completed"
            if was_cancelled:
                warnings.append(tr('Calcul annulé : les résultats éventuellement présents sont partiels.'))
            result["elapsed_seconds"] = round(time.monotonic() - started, 3)
            progress(len(final_candidates), len(survivors),
                     tr('Extension {v0} : {v1:,}/{v2:,} MFE, {v3:,} ensembles analysés, {v4:,} candidats exportés.', v0=tr('annulée') if was_cancelled else tr('terminée'), v1=seen, v2=expected, v3=len(final_candidates), v4=len(ranked)))
            return result
    except ImportError as exc:
        raise RuntimeError(tr("NUPACK ou une de ses dépendances est indisponible. Aucune extension n'a été simulée.")) from exc
