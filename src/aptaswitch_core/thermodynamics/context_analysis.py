"""Fixed full-reporter tubes and independent condition sweeps with an explicit backend."""
from __future__ import annotations

from dataclasses import asdict
import math

from ..saved_switch_analysis import (
    _AXES, _population, build_saved_switch_context, validate_analysis_config,
)
from .application_analysis import _check_context
from .structures import is_connected, parse_structure
from .tube_labels import labelled_populations, merge_species
from .tubes import solve_tube
from .types import ComplexSpec, Conditions, PartitionResult, Strand, enumerate_complexes


class _Cancelled(Exception):
    pass


def minimum_record(backend, spec, conditions, *, maximum_count=50):
    """Existing MFE display fields, with an honest bounded-alternatives count."""
    record = backend.minimum_structures(spec, conditions, maximum_count=maximum_count)
    _check_context(record.context, spec, conditions, backend.provenance)
    structures, energies = record.structures, record.energies_kcal_mol
    if (not structures or len(structures) != len(energies) or len(structures) > maximum_count
            or record.maximum_count != maximum_count or type(record.truncated) is not bool
            or len(set(structures)) != len(structures)):
        raise ValueError("Invalid bounded minimum enumeration")
    minimum, tolerance = record.minimum_energy_kcal_mol, record.tolerance_kcal_mol
    if (not math.isfinite(minimum) or not math.isfinite(tolerance) or tolerance < 0
            or (record.truncated and len(structures) != maximum_count)):
        raise ValueError("Invalid minimum window or truncation count")
    numerical = 1e-9 + 1e-12 * abs(minimum)
    if abs(energies[0] - minimum) > numerical or any(e < minimum - numerical or e > minimum + tolerance + numerical
                                                   for e in energies):
        raise ValueError("Minimum energies do not satisfy the requested global window")
    lengths = tuple(len(s.sequence) for s in spec.strands)
    for structure, energy in zip(structures, energies):
        if not math.isfinite(energy) or not is_connected(parse_structure(structure, lengths), lengths):
            raise ValueError("Invalid minimum for a connected complex")
    return {
        "structure": structures[0], "energy_kcal_mol": energies[0],
        "strand_order": [s.name for s in spec.strands], "sequences": [s.sequence for s in spec.strands],
        "alternatives": [{"structure": s, "energy_kcal_mol": e} for s, e in zip(structures, energies)],
        "alternative_count": None if record.truncated else len(structures),
        "alternatives_truncated": record.truncated, "alternatives_returned": len(structures),
        "energy_tolerance_kcal_mol": record.tolerance_kcal_mol,
    }


def analyze_saved_switch(backend, snapshot, config, *, progress_callback=None, cancel_check=None):
    """Retain every completed point; never publish a partial point as complete.

    OFF/ON share identical physical partitions at each condition. All complexes
    through the selected max_size are solved, including unwanted complexes and
    distinct cyclic orders. Baseline structures cover every nonempty ensemble;
    sweeps retain isolated Switch and formed Trigger+Switch structures only.
    """
    validated = validate_analysis_config(snapshot, config, backend=backend)
    context = build_saved_switch_context(snapshot, validated, backend=backend)
    provenance = asdict(backend.provenance)
    result = dict(status="cancelled", context=context, config=validated, baseline=None,
                  sweeps={axis: [] for axis in _AXES}, thermodynamics=provenance,
                  model_limits=dict(backend.condition_limits),
                  notes=["Séquences fixes, rapporteur complet et chimie de chaque brin conservés.",
                         "Les sous-structures décrivent les complexes formés ; les concentrations décrivent les tubes.",
                         "Le modèle thermodynamique expérimental ne prédit ni fluorescence ni cinétique."])
    strands = {name: Strand(context[f"{name.lower()}_sequence"], context[f"{name.lower()}_material"], name)
               for name in ("Switch", "Trigger", "Aptamer")}
    totals = {"off": {"Switch": validated["switch_concentration_m"]},
              "on": {"Switch": validated["switch_concentration_m"], "Trigger": validated["trigger_concentration_m"]}}
    if validated["include_aptamer"]:
        for amounts in totals.values():
            amounts["Aptamer"] = validated["aptamer_concentration_m"]
    baseline = {axis: validated[axis] for axis in _AXES}
    planned = [(None, baseline)]
    for axis, key in _AXES.items():
        planned.extend((axis, {**baseline, axis: value}) for value in validated[key])
    points = {}

    def check():
        if cancel_check and cancel_check():
            raise _Cancelled()

    try:
        for index, (axis, settings) in enumerate(planned):
            check()
            if progress_callback:
                progress_callback(index, len(planned), f"Point {index + 1}/{len(planned)} · {settings['temperature_c']:g} °C · Na⁺ {settings['sodium_m']:g} M…")
            check()
            point_key = tuple(settings[name] for name in _AXES)
            if point_key not in points:
                conditions = Conditions(settings["temperature_c"], settings["sodium_m"], settings["magnesium_m"])
                partitions, minima, populations = {}, {}, {}

                def partition(spec):
                    if spec.cyclic_key not in partitions:
                        check()
                        record = backend.partition(spec, conditions)
                        check()
                        if not isinstance(record, PartitionResult):
                            raise ValueError("Backend must return an identified partition")
                        _check_context(record, spec, conditions, backend.provenance)
                        partitions[spec.cyclic_key] = record
                    return partitions[spec.cyclic_key]

                def mfe(names):
                    if names not in minima:
                        check()
                        spec = ComplexSpec(tuple(strands[name] for name in names))
                        if partition(spec).is_empty:
                            minima[names] = None
                        else:
                            minima[names] = minimum_record(backend, spec, conditions)
                        check()
                    return minima[names]

                for state, amounts in totals.items():
                    labelled = tuple(strands[name] for name in amounts)
                    species, concentrations = merge_species(labelled, amounts.values())
                    catalog = enumerate_complexes(species, validated["max_complex_size"])
                    records = tuple(partition(spec) for spec in catalog)
                    tube = solve_tube(species, concentrations, records, standard_state=backend.standard_state)
                    check()
                    population = _population(labelled_populations(tube, labelled, amounts.values()),
                                             amounts, validated["max_complex_size"])
                    population["equilibrium"] = {
                        "ordered_physical_complex_count": len(records),
                        "chemical_aliases_merged": len(species) != len(labelled),
                        "max_relative_mass_balance_error": tube.equilibrium.max_relative_mass_balance_error,
                        "standard_state": asdict(tube.standard_state),
                    }
                    if axis is None:
                        for row in population["complexes"]:
                            row["mfe"] = mfe(tuple(row["strand_order"]))
                            row["structure_status"] = "available" if row["mfe"] is not None else "empty_ensemble"
                    populations[state] = population
                off_mfe, on_mfe = mfe(("Switch",)), mfe(("Trigger", "Switch"))
                points[point_key] = {
                    "conditions": dict(settings), "thermodynamics": provenance,
                    "off": {"population": populations["off"], "mfe": off_mfe},
                    "on": {"population": populations["on"], "mfe": on_mfe},
                    "metrics": {"on_yield_pct": populations["on"]["target_switch_fraction_pct"],
                                "on_trigger_bound_pct": populations["on"]["target_trigger_fraction_pct"],
                                "off_switch_monomer_pct": populations["off"]["switch_monomer_pct"],
                                "off_switch_associated_pct": populations["off"]["switch_associated_pct"]},
                }
            if axis is None:
                result["baseline"] = points[point_key]
            else:
                result["sweeps"][axis].append(points[point_key])
            if progress_callback:
                progress_callback(index + 1, len(planned), "Point terminé ; conservation des brins vérifiée.")
        check()
    except _Cancelled:
        return result
    except Exception as error:
        result.update(status="failed", error=str(error))
        return result
    result["status"] = "completed"
    return result
