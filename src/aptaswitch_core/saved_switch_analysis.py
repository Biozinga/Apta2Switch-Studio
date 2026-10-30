"""Fixed-sequence native validation of a retained switch with its full reporter.

Every sweep changes one physical parameter at a time. Tube yields count strand
stoichiometry at equilibrium; they do not predict fluorescence, ligand binding,
kinetics, or the concentration of extracellular nucleic acids.
"""

from __future__ import annotations

import math
from numbers import Real
from typing import Callable

from .architecture import architecture_layout
from .models import SwitchArchitecture, ThermoConditions
from .sequences import normalize_dna, normalize_rna
from .thermodynamics.tube_reporting import MASS_BALANCE_RELATIVE_TOLERANCE, NEGATIVE_ROUNDOFF_RELATIVE_TOLERANCE


MAX_SWEEP_POINTS = 25
MAX_TOTAL_POINTS = 51
MAX_TRANSCRIPT_LENGTH = 2000
MAX_COMPLEX_SIZE = 3
_AXES = {"temperature_c": "temperature_values_c", "sodium_m": "sodium_values_m", "magnesium_m": "magnesium_values_m"}
_COMPLEX_STRAND_COUNTS = (("Trigger", "trigger_count"), ("Switch", "switch_count"), ("Aptamer", "aptamer_count"))


def model_limits(material: str, *, backend=None) -> dict:
    """Supported mixed-model ranges; installed parameter checks stay authoritative."""
    if backend is not None:
        from .thermodynamics.types import ModelProvenance
        if not isinstance(backend.provenance, ModelProvenance) or material != backend.provenance.model:
            raise ValueError("The calculation model differs from the selected backend")
        return dict(backend.condition_limits)
    ThermoConditions(material=material).validate_model()
    return {
        "sodium_min_m": .12, "sodium_max_m": 1.,
        "magnesium_min_m": 0., "magnesium_max_m": 0.,
        "temperature_min_c": 25., "temperature_max_c": 55.,
    }


def _number(value, name: str, lower: float, upper: float, *, positive: bool = False) -> float:
    if isinstance(value, bool) or not isinstance(value, Real) or not math.isfinite(value):
        raise ValueError(f"{name} doit être un nombre fini.")
    value = float(value)
    if value < lower or value > upper or (positive and value <= 0):
        raise ValueError(f"{name} doit être compris entre {lower:g} et {upper:g}.")
    return value


def build_saved_switch_context(snapshot: dict, config: dict, *, backend=None) -> dict:
    """Extend only the reporter suffix, keeping the saved candidate intact.

    Input reporter is explicitly *after* its start codon: the switch already
    contains the initiation AUG upstream of its linker. The full sequence
    must begin with the reporter context used to design this exact candidate.
    """
    if not isinstance(snapshot, dict) or not isinstance(config, dict):
        raise ValueError("Le candidat et les conditions doivent être des objets enregistrés.")
    request = snapshot.get("request", {})
    candidate = snapshot.get("candidate", {})
    thermo = request.get("thermo", {})
    ThermoConditions(switch_material=thermo.get("switch_material"),
                     trigger_material=thermo.get("trigger_material")).validate_switch_chemistry()
    material = backend.provenance.model if backend is not None else thermo.get("material")
    model_limits(material, backend=backend)
    try:
        architecture = SwitchArchitecture(**request["architecture"])
        layout = architecture_layout(architecture)
    except (KeyError, TypeError, ValueError) as error:
        raise ValueError("L’architecture enregistrée est absente ou invalide.") from error
    # Candidate-specific trigger is essential for inverse-extension designs.
    switch = normalize_rna(candidate.get("switch_seq"), field_name="Switch enregistré")
    trigger = normalize_dna(candidate.get("trigger_seq"), field_name="Trigger du candidat")
    reporter_prefix = normalize_rna(request.get("reporter", {}).get("sequence_after_start_rna"), field_name="Contexte du rapporteur")
    saved_reporter = request.get("reporter", {}).get("full_sequence_after_start_rna")
    reporter = normalize_rna(saved_reporter, field_name="Rapporteur complet après AUG")
    start = layout["reporter_start"]
    if len(switch) != start + len(reporter_prefix) or switch[start:] != reporter_prefix:
        raise ValueError("La séquence du candidat ne correspond pas au rapporteur et à l’architecture enregistrés.")
    if switch[layout["aug_start"]:layout["aug_start"] + 3] != "AUG":
        raise ValueError("Le codon d’initiation du candidat enregistré est incompatible avec son architecture.")
    if not reporter.startswith(reporter_prefix):
        raise ValueError("Le rapporteur complet doit commencer par le contexte exact utilisé pour ce candidat, sans répéter son AUG initial.")
    transcript = switch[:start] + reporter
    if len(transcript) > MAX_TRANSCRIPT_LENGTH:
        raise ValueError(f"Le transcrit complet est limité à {MAX_TRANSCRIPT_LENGTH} nucléotides pour cette analyse.")
    if max(len(trigger), len(request.get("sequences", {}).get("aptamer_dna", ""))) > MAX_TRANSCRIPT_LENGTH:
        raise ValueError("Un brin dépasse la taille maximale prise en charge pour cette analyse.")
    aptamer = normalize_dna(request.get("sequences", {}).get("aptamer_dna"), field_name="Aptamère enregistré")
    return {
        "candidate_name": str(snapshot.get("name", "")), "material": material,
        "switch_sequence": transcript, "full_switch_sequence": transcript,
        "designed_switch_sequence": switch,
        "trigger_sequence": trigger, "aptamer_sequence": aptamer,
        "switch_material": "rna", "trigger_material": "dna",
        "aptamer_material": "dna", "reporter_sequence_after_start": reporter,
        "designed_reporter_prefix": reporter_prefix, "reporter_start": start,
        "aug_start": layout["aug_start"], "rbs_start": layout["rbs_start"],
        "rbs_linker_start": layout["rbs_linker_start"], "rbs_linker_end": layout["rbs_linker_end"],
        "reporter_length": len(reporter), "switch_length": len(transcript),
        "reporter_has_complete_codons": len(reporter) % 3 == 0,
        "reporter_source": "saved_full_reporter",
    }


def validate_analysis_config(snapshot: dict, config: dict, *, backend=None) -> dict:
    """Validate all points before starting native calculations."""
    context = build_saved_switch_context(snapshot, config, backend=backend)
    thermo = snapshot["request"].get("thermo", {})
    limits = model_limits(context["material"], backend=backend)
    result = {}
    for key, lower, upper, default in (
        ("temperature_c", limits["temperature_min_c"], limits["temperature_max_c"], 37.),
        ("sodium_m", limits["sodium_min_m"], limits["sodium_max_m"], 1.),
        ("magnesium_m", 0., limits["magnesium_max_m"], 0.),
    ):
        result[key] = _number(config.get(key, thermo.get(key, default)), key, lower, upper)
    shared_concentration = thermo.get("concentration_m", 5e-6)
    for key in ("switch_concentration_m", "trigger_concentration_m", "aptamer_concentration_m"):
        fallback = thermo.get(key)
        if fallback is None and key == "aptamer_concentration_m":
            fallback = thermo.get("trigger_concentration_m")
        if fallback is None:
            fallback = shared_concentration
        result[key] = _number(config.get(key, fallback), key, 0., 1., positive=True)
    include = config.get("include_aptamer", False)
    if type(include) is not bool:
        raise ValueError("L’inclusion de l’aptamère doit être un booléen.")
    result["include_aptamer"] = include
    max_size = config.get("max_complex_size", MAX_COMPLEX_SIZE)
    if type(max_size) is not int or not 2 <= max_size <= MAX_COMPLEX_SIZE:
        raise ValueError("La taille maximale des complexes du rapporteur complet doit être un entier entre 2 et 3.")
    result["max_complex_size"] = max_size
    for axis, key in _AXES.items():
        values = config.get(key, [])
        if not isinstance(values, (list, tuple)) or len(values) > MAX_SWEEP_POINTS:
            raise ValueError(f"Chaque balayage doit contenir au maximum {MAX_SWEEP_POINTS} points.")
        if axis == "magnesium_m" and values and not limits["magnesium_max_m"]:
            raise ValueError("Le balayage Mg²⁺ est exclu pour les modèles ARN et ADN/ARN : Mg²⁺ = 0.")
        bounds = {
            "temperature_c": (limits["temperature_min_c"], limits["temperature_max_c"]),
            "sodium_m": (limits["sodium_min_m"], limits["sodium_max_m"]),
            "magnesium_m": (0., limits["magnesium_max_m"]),
        }[axis]
        result[key] = sorted(set(_number(value, axis, *bounds) for value in values))
    if 1 + sum(len(values) for key, values in result.items() if key in _AXES.values()) > MAX_TOTAL_POINTS:
        raise ValueError(f"L’analyse est limitée à {MAX_TOTAL_POINTS} points, référence comprise.")
    return result


def _cyclic_strand_order(names: list[str]) -> tuple[str, ...]:
    """Identify rotations, preserving distinct connected complex ensembles.

    For example, Trigger/Switch/Aptamer and Trigger/Aptamer/Switch have the
    same composition but different ensembles. Only rotating the strand cycle
    preserves its identity; sorting or reversing it would merge complexes.
    """
    priority = {name: index for index, (name, _) in enumerate(_COMPLEX_STRAND_COUNTS)}
    cycle = tuple(names)
    return min((cycle[index:] + cycle[:index] for index in range(len(cycle))),
               key=lambda order: tuple(priority[name] for name in order))


def _population(concentrations, totals: dict, max_size: int) -> dict:
    tolerance = max(totals.values()) * NEGATIVE_ROUNDOFF_RELATIVE_TOLERANCE
    grouped: dict[tuple, list[float]] = {}
    for complex_obj, value in concentrations.items():
        concentration = float(value)
        names = [strand.name for strand in complex_obj.strands]
        if (not math.isfinite(concentration) or concentration < -tolerance or not names
                or len(names) > max_size or any(name not in totals for name in names)):
            raise RuntimeError("Le moteur a retourné un complexe ou une concentration invalide.")
        key = _cyclic_strand_order(names)
        grouped.setdefault(key, []).append(max(0., concentration))
    values = {key: math.fsum(items) for key, items in grouped.items()}
    for name in ("Switch", "Trigger", "Aptamer"):
        recovered = math.fsum(names.count(name) * concentration for names, concentration in values.items())
        expected = totals.get(name, 0.)
        if not math.isclose(recovered, expected, rel_tol=MASS_BALANCE_RELATIVE_TOLERANCE, abs_tol=0.):
            raise RuntimeError("Les concentrations ne conservent pas les quantités de brins du tube.")
    rows = []
    composition_orders = {}
    for names in values:
        composition = tuple(names.count(name) for name in ("Switch", "Trigger", "Aptamer"))
        composition_orders[composition] = composition_orders.get(composition, 0) + 1
    for names, concentration in values.items():
        s, t, a = (names.count(name) for name in ("Switch", "Trigger", "Aptamer"))
        label = " + ".join(name if count == 1 else f"{count} {name}" for name, count in (("Switch", s), ("Trigger", t), ("Aptamer", a)) if count)
        rows.append({
            "switch_count": s, "trigger_count": t, "aptamer_count": a,
            "strand_order": list(names), "complex_key": "|".join(names),
            "composition_label": label,
            "label": label + (f" [{' → '.join(names)}]" if composition_orders[(s, t, a)] > 1 else ""),
            "concentration_m": concentration,
            "switch_fraction_pct": 100 * s * concentration / totals["Switch"],
            "trigger_fraction_pct": 100 * t * concentration / totals["Trigger"] if totals.get("Trigger") else 0.,
        })
    rows.sort(key=lambda row: -row["concentration_m"])
    free = values.get(("Switch",), 0.) / totals["Switch"] * 100
    target = values.get(("Trigger", "Switch"), 0.) / totals["Switch"] * 100
    trigger_target = values.get(("Trigger", "Switch"), 0.) / totals["Trigger"] * 100 if totals.get("Trigger") else 0.
    return {
        "total_strand_concentrations_m": dict(totals), "complexes": rows,
        "switch_monomer_pct": min(100., max(0., free)),
        "switch_associated_pct": min(100., max(0., 100 - free)),
        "target_switch_fraction_pct": min(100., max(0., target)),
        "target_trigger_fraction_pct": min(100., max(0., trigger_target)),
        "mass_balance_checked": True,
    }



def analyze_saved_switch(
    snapshot: dict, import_path: str | None, config: dict, *,
    progress_callback: Callable[[int, int, str], None] | None = None,
    cancel_check: Callable[[], bool] | None = None,
    backend=None,
    thermodynamics_library_path: str | None = None,
) -> dict:
    """Analyze a fixed full-reporter candidate using the configured native model.

    Original snapshot provenance is retained. New results identify the actual
    native calculation; the historical ``import_path`` is ignored.
    """
    from .thermodynamics.context_analysis import analyze_saved_switch as analyze_native
    from .thermodynamics.setup import get_backend

    if backend is None:
        backend = get_backend(thermodynamics_library_path)
    return analyze_native(backend, snapshot, config, progress_callback=progress_callback,
                          cancel_check=cancel_check)
