"""Functional design priorities from the actual OFF/ON tube evaluation.

The four local defects remain conditional structural probabilities. Separate
guard measurements record missing OFF/ON complexes and the expected ON duplex;
neither these measurements nor the global tube defect is an expression ratio.
"""
from __future__ import annotations

from dataclasses import asdict
import math

from ..defect_metrics import functional_defect_metrics, _pair_targets
from .design_constraints import BaseSite
from .optimization import ObjectiveValue


REGION_NAMES = ("off_toehold", "off_stem", "on_rbs_aug", "on_rbs_linker")
GUARD_NAMES = ("off_complex_shortage", "on_complex_shortage", "on_recognition_defect")


def objective_from_evaluation(request, evaluated, *, rbs_minimum=None, rbs_prediction=None, trigger_pairs=None) -> ObjectiveValue:
    """Reuse already-calculated matrices, with exact trigger+switch coordinates."""
    problem, assignment = evaluated.problem, evaluated.assignment
    strands = {s.name: s for s in assignment.strands}
    trigger_length = len(strands["Trigger"].sequence)
    switch_length = len(strands["Switch"].sequence)
    records = {tuple(s.name for s in a.partition.complex_spec.strands): a for a in evaluated.analyses}
    off, on = records[("Switch",)].pairs, records[("Trigger", "Switch")].pairs
    # A feasible complementary target must have a connected ensemble. Missing
    # probabilities here indicate inconsistent engine/constraint results.
    if off is None or on is None:
        raise ValueError("Functional optimization requires the connected OFF and ON probability matrices")
    regional = functional_defect_metrics(request.architecture,
        len(request.reporter.sequence_after_start_rna), off.matrix, on.matrix,
        weights=request.functional_defect_weights)
    weights = asdict(request.functional_defect_weights)
    fractions = request.functional_defect_weights.fractions()
    local_weights, global_weights = {}, {}
    off_target = _pair_targets(regional["target_structures"]["off"], [switch_length])
    on_target = _pair_targets(regional["target_structures"]["on"], [trigger_length, switch_length])
    off_sites = tuple(BaseSite("Switch", i) for i in range(switch_length))
    on_sites = tuple(BaseSite("Trigger", i) for i in range(trigger_length)) + off_sites

    def add_errors(destination, matrix, targets, sites, indices, factor):
        if not indices or factor <= 0:
            return
        scale = factor / len(indices)
        for i in indices:
            expected = targets[i]
            error = 1.0 - matrix[i][expected]
            destination[sites[i]] = destination.get(sites[i], 0.) + scale * error
            # In particular, a fixed toehold may be trapped by a mutable linker.
            # Guide mutations to the unintended partner instead of trying to
            # mutate only the fixed erroneous base itself.
            for j, probability in enumerate(matrix[i]):
                if j != i and j != expected and probability > 0:
                    destination[sites[j]] = destination.get(sites[j], 0.) + scale * probability

    for name in REGION_NAMES:
        is_on = name.startswith("on_")
        offset = trigger_length if is_on else 0
        indices = [i + offset for i in regional["regions"][name]["indices_switch"]]
        add_errors(local_weights, on.matrix if is_on else off.matrix,
                   on_target if is_on else off_target, on_sites if is_on else off_sites,
                   indices, fractions[name])
    add_errors(global_weights, off.matrix, off_target, off_sites, list(range(switch_length)), .5)
    add_errors(global_weights, on.matrix, on_target, on_sites, list(range(len(on_sites))), .5)
    layout = problem.layout
    rbs_energy = 0.0
    if rbs_minimum is not None:
        from .structures import parse_structure
        start, end = layout["rbs_linker_start"], layout["rbs_linker_end"]
        partners = parse_structure(rbs_minimum.structure, [end - start]).partners
        rbs_energy = rbs_minimum.energy_kcal_mol
        paired = [i for i, partner in enumerate(partners) if partner != i]
        # Guide the refinement toward actual paired sites of the isolated RNA
        # fragment. The local phase retains its four functional priorities.
        penalty = request.scoring.rbs * (max(0., -rbs_energy) / 100.)
        for i in paired:
            site = BaseSite("Switch", start + i)
            global_weights[site] = global_weights.get(site, 0.) + penalty / len(paired)
    recognition = list(range(trigger_length + layout["len_t7"],
                             trigger_length + layout["len_t7"] + layout["len_trigger"]))
    recognition_defect = math.fsum(1.0-on.matrix[i][on_target[i]] for i in recognition) / len(recognition)
    add_errors(local_weights, on.matrix, on_target, on_sites, recognition, .25)
    # The local phase keeps its original structural priorities. These extra
    # costs steer refinement under the already-fixed regional/tube guards.
    rbs_merit = rbs_prediction["merit"] if rbs_prediction else 0.0
    rbs_calculator_weight = request.scoring.rbs_calculator if rbs_prediction else 0.0
    if rbs_calculator_weight:
        window_start = max(0, layout["aug_start"] - 35)
        window_end = min(switch_length, layout["aug_start"] + 35)
        scale = rbs_calculator_weight * (1.0 - rbs_merit) / (window_end - window_start)
        for index in range(window_start, window_end):
            site = BaseSite("Switch", index)
            global_weights[site] = global_weights.get(site, 0.) + scale
    linearity = linearity_weight = 0.0
    if request.inverse_extension is not None and request.scoring.trigger_linearity > 0:
        pairs = trigger_pairs if trigger_pairs is not None else records[("Trigger",)].pairs
        if pairs is None:
            raise ValueError("Trigger refinement requires its isolated probability matrix")
        linearity = math.fsum(1.0 - pairs.matrix[i][i] for i in range(trigger_length)) / trigger_length
        linearity_weight = request.scoring.trigger_linearity
        for i in range(trigger_length):
            site = BaseSite("Trigger", i)
            global_weights[site] = global_weights.get(site, 0.) + linearity_weight * (1.0 - pairs.matrix[i][i]) / trigger_length
    contributions = {(tube.name, c.target.name): c
                     for tube in evaluated.tubes for c in tube.defect.contributions}

    def shortage(tube, target):
        contribution = contributions[(tube, target)]
        return max(0.0, 1.0-min(1.0, contribution.actual_concentration_molar /
                               contribution.target.concentration_molar))

    return ObjectiveValue(
        local_defects=tuple(regional[name] for name in REGION_NAMES),
        local_weights=tuple(weights[name] for name in REGION_NAMES),
        global_defect=evaluated.defect.normalized,
        guard_defects=(shortage("tube_off", "switch_off"), shortage("tube_on", "switch_on"), recognition_defect),
        site_weights=tuple(sorted(local_weights.items())),
        global_site_weights=tuple(sorted(global_weights.items())),
        rbs_linker_energy=rbs_energy,
        rbs_weight=request.scoring.rbs if rbs_minimum is not None else 0.0,
        rbs_calculator_merit=rbs_merit, rbs_calculator_weight=rbs_calculator_weight,
        trigger_linearity_defect=linearity, trigger_linearity_weight=linearity_weight,
    )
