"""Native OFF/ON sampling for the existing four-region family presentation."""
from __future__ import annotations

from collections import Counter
from dataclasses import asdict

from ..localization import tr
from ..structure_ensemble import EnsembleCancelled, assemble_ensemble, ensemble_inputs
from .application_analysis import _check_context
from .structures import is_connected, parse_structure
from .types import ComplexSpec, Conditions, ModelProvenance, Strand


def native_ensemble_inputs(request, candidate, backend, *, num_sample=1000, seed=0):
    """Cache identity includes actual model, chemistry, sample count and RNG seeds."""
    if type(seed) is not int or not 0 <= seed < 2**64 - 1:
        raise ValueError("Ensemble seed must allow two consecutive uint64 seeds")
    if not isinstance(backend.provenance, ModelProvenance):
        raise ValueError("An explicit thermodynamic model identity is required")
    request.thermo.validate_switch_chemistry()
    inputs = ensemble_inputs(request, candidate, num_sample)
    # Calculated candidates must already carry their actual strand chemistry.
    # Do not silently convert stored RNA triggers or DNA switches.
    Strand(inputs["switch"], "rna", "Switch")
    Strand(inputs["trigger"], "dna", "Trigger")
    inputs["thermo"]["material"] = backend.provenance.model
    inputs["thermodynamics"] = asdict(backend.provenance)
    inputs["sampling"] = {"method": "joint-connected-boltzmann", "seeds": {"off": seed, "on": seed + 1}}
    return inputs


def analyze_structure_ensemble(backend, request, candidate, *, num_sample=1000, seed=0,
                               cancel_check=None, progress_callback=None):
    """Sample formed complexes, not a concentration-weighted mixture of species.

    The same native inside table serves probabilities and sampling. Family
    grouping is shared with the existing application; no marginal-product
    estimates, artificial consensus structure, or sequence redesign is used.
    """
    inputs = native_ensemble_inputs(request, candidate, backend, num_sample=num_sample, seed=seed)
    conditions = Conditions(request.thermo.temperature_c, request.thermo.sodium_m, request.thermo.magnesium_m)
    switch = Strand(inputs["switch"], "rna", "Switch")
    trigger = Strand(inputs["trigger"], "dna", "Trigger")
    matrices, samples = {}, {}

    def check():
        if cancel_check and cancel_check():
            raise EnsembleCancelled()

    for position, (state, strands) in enumerate((("off", (switch,)), ("on", (trigger, switch)))):
        check()
        if progress_callback:
            progress_callback(position, 2, tr("Ensemble {state} · échantillonnage de {count} structures…",
                                              state=state.upper(), count=num_sample))
        check()
        spec = ComplexSpec(strands)
        result = backend.ensemble(spec, conditions, num_sample, seed=inputs["sampling"]["seeds"][state])
        check()
        partition, pairs = result.analysis.partition, result.analysis.pairs
        _check_context(partition, spec, conditions, backend.provenance)
        _check_context(result.samples.context, spec, conditions, backend.provenance)
        if (pairs is None or pairs.partition != partition
                or result.samples.context.thermal_energy_kcal_mol != partition.thermal_energy_kcal_mol
                or result.samples.seed != inputs["sampling"]["seeds"][state]
                or len(result.samples.structures) != num_sample):
            raise ValueError("Native ensemble does not match its probabilities or sampling request")
        lengths = tuple(len(s.sequence) for s in strands)
        counts = Counter(result.samples.structures)
        for structure in counts:
            if not is_connected(parse_structure(structure, lengths), lengths):
                raise ValueError("A sampled complex is disconnected")
        matrices[state], samples[state] = pairs.matrix, counts
    check()
    result = assemble_ensemble(inputs, matrices, samples, inputs["thermodynamics"])
    if progress_callback:
        progress_callback(2, 2, tr("Ensembles OFF/ON calculés."))
    check()
    return result
