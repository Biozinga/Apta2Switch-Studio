"""Injectable replacements for the application's fixed system/tube analyses.

Nothing in this module selects an engine or updates a saved
design. Existing display fields are preserved, while ``thermodynamics`` records
the actual model rather than attributing its predictions to the former engine.
The aptamer and original trigger are both DNA.
"""
from __future__ import annotations

from dataclasses import asdict
import math
from numbers import Integral
from typing import Callable, Protocol

from ..models import ThermoConditions
from ..sequences import normalize_dna, normalize_rna
from .structures import format_structure, is_connected, parse_structure
from .tube_reporting import TubeAnalysisConditions, _population
from .tube_labels import labelled_populations
from .tubes import StandardStateConvention, TubeResult, solve_tube
from .types import (
    ComplexSpec, Conditions, ModelProvenance, PartitionResult, Strand,
    enumerate_complexes,
)


class FixedAnalysisBackend(Protocol):
    """The native adapter supplies these methods without an application import.

    ``minimum`` returns a record with ``context``, ``structure`` and
    ``energy_kcal_mol``. Context must identify the exact request and provenance.
    An optional ``minimum_structures(..., maximum_count=...)`` supplies multiple
    minimum witnesses and must explicitly report whether enumeration truncated.
    """

    provenance: ModelProvenance

    def partition(self, spec: ComplexSpec, conditions: Conditions) -> PartitionResult: ...

    def minimum(self, spec: ComplexSpec, conditions: Conditions): ...


def _inputs(backend: FixedAnalysisBackend, aptamer_sequence: str,
            trigger_sequence: str, thermo: ThermoConditions):
    provenance = getattr(backend, "provenance", None)
    if not isinstance(provenance, ModelProvenance):
        raise ValueError("The analysis requires explicit backend model provenance")
    if not isinstance(thermo, ThermoConditions):
        raise ValueError("thermo must be ThermoConditions")
    thermo.validate_system_chemistry()
    conditions = Conditions(thermo.temperature_c, thermo.sodium_m, thermo.magnesium_m)
    # The injected backend enforces its own temperature/salt domain.
    a = Strand(normalize_dna(aptamer_sequence, field_name="aptamer"), "dna", "A")
    t = Strand(normalize_dna(trigger_sequence, field_name="trigger"), "dna", "T")
    settings = {
        "material": provenance.model,
        "temperature_c": conditions.celsius,
        "sodium_m": conditions.sodium,
        "magnesium_m": conditions.magnesium,
        "ensemble": provenance.ensemble,
        "aptamer_material": a.material,
        "trigger_material": t.material,
    }
    return a, t, conditions, settings, provenance


def _check_context(record, spec: ComplexSpec, conditions: Conditions,
                   provenance: ModelProvenance):
    if (record.complex_spec != spec or record.conditions != conditions
            or record.provenance != provenance or record.constraints is not None):
        raise ValueError("Backend changed the requested coordinates, chemistry, conditions or ensemble")


def analyze_system(
    backend: FixedAnalysisBackend,
    aptamer_dna: str,
    trigger_dna: str,
    thermo: ThermoConditions,
    *,
    mfe_policy: str = "available",
    max_mfe_structures: int = 100,
) -> dict:
    """Return the initial two-strand minimum and its intermolecular contacts.

    ``available`` requests bounded alternative minima when the backend supports
    them; otherwise one witness is explicitly marked non-exhaustive. ``complete``
    requires exhaustive enumeration and fails if the backend cannot provide it
    within ``max_mfe_structures``. This is a display/enumeration policy, not a
    change to the minimum-energy objective or the thermodynamic ensemble.

    An empty connected ensemble is an error; no fabricated all-unpaired dimer
    or replacement chemistry is supplied to the application's 2D/3D renderer.
    """
    if mfe_policy not in ("available", "complete"):
        raise ValueError("mfe_policy must be 'available' or 'complete'")
    if (isinstance(max_mfe_structures, bool) or not isinstance(max_mfe_structures, Integral)
            or max_mfe_structures < 1):
        raise ValueError("max_mfe_structures must be a positive integer")
    a, t, conditions, settings, provenance = _inputs(backend, aptamer_dna, trigger_dna, thermo)
    spec = ComplexSpec((a, t))
    enumerate_minima = getattr(backend, "minimum_structures", None)
    if callable(enumerate_minima):
        record = enumerate_minima(spec, conditions, maximum_count=int(max_mfe_structures),
                                  absolute_tolerance=1e-8)
        if not isinstance(record.truncated, bool):
            raise ValueError("Minimum enumeration must explicitly report truncation")
        structures = tuple(record.structures)
        energies = tuple(float(value) for value in record.energies_kcal_mol)
        energy = float(record.minimum_energy_kcal_mol)
        tolerance = float(record.tolerance_kcal_mol)
        if not structures or len(structures) > max_mfe_structures:
            raise ValueError("Minimum enumeration returned an invalid number of structures")
        if (len(energies) != len(structures) or tolerance != 1e-8
                or record.maximum_count != max_mfe_structures):
            raise ValueError("Minimum enumeration changed its requested bound or tolerance")
        exhaustive, truncated = not record.truncated, record.truncated
        method = "minimum_structures"
    else:
        if mfe_policy == "complete":
            raise ValueError("This backend cannot enumerate all alternative minima")
        record = backend.minimum(spec, conditions)
        structures = (record.structure,)
        energy = float(record.energy_kcal_mol)
        energies, tolerance = (energy,), None
        exhaustive, truncated, method = False, None, "minimum"
    _check_context(record.context, spec, conditions, provenance)
    if not math.isfinite(energy) or any(not math.isfinite(value) for value in energies):
        raise ValueError("The backend returned a nonfinite minimum energy")
    numerical = 1e-9 + 1e-12 * max(abs(energy), *(abs(value) for value in energies))
    if (abs(energies[0] - energy) > numerical
            or any(value < energy - numerical or value > energy + (tolerance or 0.) + numerical
                   for value in energies)):
        raise ValueError("Minimum witnesses are inconsistent with their energy tolerance")
    lengths = tuple(len(strand.sequence) for strand in spec.strands)
    alternatives, states = [], []
    for structure, witness_energy in zip(structures, energies):
        state = parse_structure(structure, lengths)
        if not is_connected(state, lengths):
            raise ValueError("The backend returned a disconnected minimum for a connected complex")
        formatted = format_structure(state, lengths)
        if any(item["structure"] == formatted for item in alternatives):
            raise ValueError("Minimum enumeration returned a duplicate structure")
        alternatives.append({"structure": formatted, "energy_kcal_mol": witness_energy})
        states.append(state)
    if mfe_policy == "complete" and not exhaustive:
        raise ValueError("Minimum enumeration reached its limit; completeness is not established")
    n_a = len(a.sequence)
    paired_trigger_indices = sorted(
        partner - n_a for index, partner in enumerate(states[0].partners)
        if index < n_a <= partner
    )
    return {
        "aptamer_sequence": a.sequence,
        "trigger_sequence": t.sequence,
        "mfe_structure": alternatives[0]["structure"],
        "mfe_energy_kcal_mol": energy,
        "paired_trigger_indices": paired_trigger_indices,
        "model": settings,
        "thermodynamics": asdict(provenance),
        "strand_order": ["aptamer", "trigger"],
        "mfe_structures": alternatives,
        "mfe_enumeration": {
            "method": method, "exhaustive": exhaustive, "truncated": truncated,
            "returned_count": len(alternatives), "maximum_count": int(max_mfe_structures),
            "energy_tolerance_kcal_mol": tolerance,
        },
    }


def _report_equilibrium(tube: TubeResult, conditions: TubeAnalysisConditions, *,
                        mixed: bool, aliases: tuple[Strand, Strand] | None = None) -> dict:
    """Report a solved catalog, explicitly merging identical chemical inputs.

    If A and T are chemically identical they form one physical species with the
    sum of their totals. Their *origin labels*, which cannot affect folding,
    are then distributed binomially within each n-mer. This conserves each
    supplied label total without inventing different chemical complexes.
    """
    populations = {}
    if aliases is None:
        populations = {record.complex_spec: float(concentration)
                       for record, concentration in zip(tube.partitions, tube.concentrations_molar)}
    else:
        a, t = aliases
        populations = labelled_populations(tube, (a, t),
            (conditions.aptamer_concentration_m, conditions.trigger_concentration_m))
    result = _population(populations, conditions, mixed=mixed)
    result["equilibrium"] = {
        "standard_state": asdict(tube.standard_state),
        "chemical_aliases_merged": aliases is not None,
        "ordered_complex_count": len(tube.partitions),
        "ordered_complexes": [{
            "strand_order": [strand.name for strand in record.complex_spec.strands],
            "concentration_m": float(concentration),
            "empty_ensemble": record.is_empty,
            "ensemble_free_energy_kcal_mol": None if record.is_empty else record.free_energy_kcal_mol,
        } for record, concentration in zip(tube.partitions, tube.concentrations_molar)],
    }
    return result


def analyze_parasite_tubes(
    backend: FixedAnalysisBackend,
    aptamer_dna: str,
    trigger_dna: str,
    thermo: ThermoConditions,
    conditions: TubeAnalysisConditions = TubeAnalysisConditions(),
    *,
    standard_state: StandardStateConvention,
    cancel_check: Callable[[], bool] | None = None,
    progress_callback: Callable[[int, int, str], None] | None = None,
) -> dict:
    """Compute both complete tube catalogs through 1–3 strands using PFs.

    A single ordered partition cache serves both tubes. Equilibrium is always
    solved by :func:`solve_tube`; no concentration fit, MFE approximation, or
    symmetry correction is introduced here. ``standard_state`` is mandatory:
    a concentration standard must be established, not guessed from a label.

    Cancellation is checked between partition calculations and tube solves.
    An in-flight native call is left to the injected backend's cancellation
    facility. No partial populations are returned as completed predictions.
    """
    if not isinstance(standard_state, StandardStateConvention):
        raise ValueError("An explicit StandardStateConvention is required")
    if not isinstance(conditions, TubeAnalysisConditions):
        raise ValueError("conditions must be TubeAnalysisConditions")
    conditions = TubeAnalysisConditions(**asdict(conditions))
    a, t, native_conditions, settings, provenance = _inputs(backend, aptamer_dna, trigger_dna, thermo)
    if standard_state.provenance_label != provenance.standard_state:
        raise ValueError("Backend provenance does not match the declared concentration standard")
    result = {
        "status": "cancelled", "model": settings, "conditions": asdict(conditions),
        "aptamer_sequence": a.sequence, "trigger_sequence": t.sequence,
        "thermodynamics": asdict(provenance), "tubes": {},
    }
    if cancel_check and cancel_check():
        return result
    merged = a.chemical_identity == t.chemical_identity
    species = (t,) if merged else (a, t)
    catalog = enumerate_complexes(species, conditions.max_complex_size)
    total_steps = len(catalog) + 2
    if progress_callback:
        progress_callback(0, total_steps, "Connected ensembles for the two tube catalogs…")
    partitions = []
    for index, spec in enumerate(catalog):
        if cancel_check and cancel_check():
            return result
        record = backend.partition(spec, native_conditions)
        if not isinstance(record, PartitionResult):
            raise ValueError("Backend must return a PartitionResult")
        _check_context(record, spec, native_conditions, provenance)
        partitions.append(record)
        if cancel_check and cancel_check():
            return result
        if progress_callback:
            progress_callback(index + 1, total_steps, "Connected complex partitions calculated…")
    trigger_partitions = tuple(record for record in partitions
                               if all(s.chemical_identity == t.chemical_identity
                                      for s in record.complex_spec.strands))
    single = solve_tube((t,), (conditions.trigger_concentration_m,), trigger_partitions,
                        standard_state=standard_state)
    if cancel_check and cancel_check():
        return result
    if progress_callback:
        progress_callback(len(catalog) + 1, total_steps, "Trigger-only equilibrium calculated…")
    totals = ((conditions.aptamer_concentration_m + conditions.trigger_concentration_m,) if merged else
              (conditions.aptamer_concentration_m, conditions.trigger_concentration_m))
    mixed = solve_tube(species, totals, partitions, standard_state=standard_state)
    if cancel_check and cancel_check():
        return result
    populations = {
        "trigger_only": _report_equilibrium(single, conditions, mixed=False),
        "aptamer_trigger": _report_equilibrium(mixed, conditions, mixed=True,
                                                aliases=(a, t) if merged else None),
    }
    result.update(status="completed", tubes=populations)
    if progress_callback:
        progress_callback(total_steps, total_steps, "Both tube populations calculated; strand conservation checked.")
    return result


def analyze_strand(backend, sequence: str, material: str, thermo: ThermoConditions) -> dict:
    """An isolated strand's actual minimum, for oxDNA guides."""
    if material not in ("dna", "rna"):
        raise ValueError("Guide chemistry must be explicitly DNA or RNA")
    normalize = normalize_dna if material == "dna" else normalize_rna
    strand = Strand(normalize(sequence, field_name="strand"), material, "strand")
    spec = ComplexSpec((strand,))
    conditions = Conditions(thermo.temperature_c, thermo.sodium_m, thermo.magnesium_m)
    record = backend.minimum(spec, conditions)
    _check_context(record.context, spec, conditions, backend.provenance)
    if not math.isfinite(record.energy_kcal_mol):
        raise ValueError("The isolated-strand minimum energy is nonfinite")
    state = parse_structure(record.structure, (len(strand.sequence),))
    return {"sequence": strand.sequence, "material": material,
            "structure": format_structure(state, (len(strand.sequence),)),
            "mfe_energy_kcal_mol": record.energy_kcal_mol,
            "thermodynamics": asdict(backend.provenance),
            "conditions": {"temperature_c": conditions.celsius, "sodium_m": conditions.sodium,
                           "magnesium_m": conditions.magnesium, "ensemble": backend.provenance.ensemble}}
