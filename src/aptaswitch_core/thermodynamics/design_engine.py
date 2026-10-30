"""Independent multi-state design and existing candidate screening/ranking.

This orchestration requires an explicit thermodynamic backend. It returns an
identified in-memory batch; the application adapter persists that identity and
the actual optimizer outcomes alongside the candidate measurements.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass, replace
from collections import OrderedDict
from datetime import datetime
import math
from typing import Callable

from ..defect_metrics import functional_defect_metrics
from ..models import DesignCandidate, RunRequest
from ..scoring import compute_selection_score, selection_score_100, sort_candidates
from ..sequences import normalize_rna
from .defects import DesignDefect
from .accessibility import unpaired_window_accessibility
from .design_constraints import DesignAssignment, DesignProblem, build_design_problem
from .evaluation import ComplexAnalysisCache, DesignEvaluation, evaluate_design
from .design_objective import objective_from_evaluation
from .optimization import OptimizationResult, ScoredAssignment, TrialOptions, assignment_key, optimize
from .tubes import solve_tube
from .tube_labels import labelled_populations, merge_species
from .types import ComplexSpec, Conditions, ModelProvenance, Strand, enumerate_complexes


def _rbs_measurement(request, problem, assignment):
    """Keep the DNA-bound projection separate from the physical hybrid backend."""
    if request.scoring.rbs_calculator == 0:
        return {}
    from ..rbs_calculator import predict_rbs_switch
    layout = problem.layout
    return predict_rbs_switch(assignment.sequence("Switch"), layout["aug_start"],
        tuple(range(layout["len_t7"], layout["len_t7"] + layout["len_trigger"])))


@dataclass(frozen=True, slots=True)
class NativeDesignBatch:
    """Explicit model and optimizer outcomes, without legacy model relabelling."""
    candidates: tuple[DesignCandidate, ...]
    optimization_trials: tuple[OptimizationResult, ...]
    provenance: ModelProvenance
    conditions: Conditions
    status: str
    started_at: str
    completed_at: str
    analysis_cache: dict


def defect_breakdown(defect: DesignDefect) -> dict:
    """Preserve the application's reporting shape from typed tube quantities."""
    tubes = []
    for tube in defect.tubes:
        structural, concentration = tube.structural_nucleotide_molar, tube.shortage_nucleotide_molar
        tubes.append({
            "name": tube.name, "global": tube.normalized,
            "structural": tube.normalized_structural,
            "concentration": tube.normalized_concentration,
            "dimensional_defect": structural + concentration,
            "dimensional_structural": structural, "dimensional_concentration": concentration,
            "target_nucleotide_concentration_m": tube.target_nucleotide_molar,
            "complexes": [{"name": c.target.name, "structural": c.structural_nucleotide_molar,
                          "concentration": c.shortage_nucleotide_molar,
                          "total": c.total_nucleotide_molar} for c in tube.contributions],
        })
    return {"global": defect.normalized, "structural": defect.structural,
            "concentration": defect.concentration, "tubes": tubes}


class _CandidateCalculations:
    """One candidate's complete partitions, never shared across mutated bases."""
    def __init__(self, backend, conditions, evaluated=None):
        self.backend, self.conditions = backend, conditions
        self.analyses = {a.partition.complex_spec.ordered_key: a for a in evaluated.analyses} if evaluated else {}

    def analysis(self, spec, *, need_pairs=False):
        key = spec.ordered_key
        result = self.analyses.get(key)
        if result is None or (need_pairs and not result.partition.is_empty and result.pairs is None):
            result = self.backend.analyze(spec, self.conditions, need_pairs=need_pairs)
            if result.partition.complex_spec != spec or result.partition.conditions != self.conditions:
                raise ValueError("Backend changed the requested candidate context")
            if result.partition.provenance != self.backend.provenance or result.partition.constraints is not None:
                raise ValueError("Candidate analysis needs the selected unrestricted model")
            self.analyses[key] = result
        return result

    def complex_percent(self, species, totals, expected, denominator, max_size):
        physical, amounts = merge_species(species, totals)
        catalog = enumerate_complexes(physical, max_size)
        records = [self.analysis(spec).partition for spec in catalog]
        tube = solve_tube(physical, amounts, records, standard_state=self.backend.standard_state)
        names = tuple(s.name for s in expected)
        key = min(names[i:] + names[:i] for i in range(len(names)))
        actual = math.fsum(concentration for spec, concentration in labelled_populations(tube, species, totals).items()
                          if tuple(s.name for s in spec.strands) == key)
        return round(100 * actual / denominator, 4)


class _EvaluationCache:
    """Recent complete evaluations, bounded independently by entries and cells.

    Counting shared matrices more than once is deliberately conservative. This
    prevents the recent-assignment cache from retaining many full reporters
    after the complex cache has already evicted their probability matrices.
    """
    def __init__(self, *, max_entries, max_pair_cells):
        self.max_entries, self.max_pair_cells = max_entries, max_pair_cells
        self.pair_cells = 0
        self._values = OrderedDict()

    def remember(self, evaluation):
        key = assignment_key(evaluation.assignment)
        previous = self._values.pop(key, None)
        if previous is not None:
            self.pair_cells -= previous[1]
        cells = sum(a.partition.complex_spec.n_nucleotides ** 2
                    for a in evaluation.analyses if a.pairs is not None)
        if self.max_entries and cells <= self.max_pair_cells:
            self._values[key] = (evaluation, cells)
            self.pair_cells += cells
            while len(self._values) > self.max_entries or self.pair_cells > self.max_pair_cells:
                _, (_, evicted_cells) = self._values.popitem(last=False)
                self.pair_cells -= evicted_cells

    def get(self, assignment):
        key = assignment_key(assignment)
        result = self._values.get(key)
        if result is not None:
            self._values.move_to_end(key)
            return result[0]
        return None


def evaluate_candidate(request: RunRequest, problem: DesignProblem, assignment: DesignAssignment,
                       backend, *, trial: int, evaluated: DesignEvaluation | None = None) -> DesignCandidate:
    """Apply the current STOP, all-minimum ON, thermodynamic and ranking rules.

    STOP and ON exclusions intentionally retain no fabricated uncomputed metrics.
    An optimizer budget ending before its loss target is not a sequence-screen
    exclusion: its actual outcome is retained separately by ``NativeDesignBatch``.
    """
    request = request.normalized()
    if problem != build_design_problem(request):
        raise ValueError("Candidate constraints and target tubes differ from the request")
    problem.validate(assignment)
    if type(trial) is not int or trial < 1:
        raise ValueError("Candidate trial must be a positive integer")
    conditions = problem.conditions
    if evaluated is not None and (evaluated.problem != problem or evaluated.assignment != assignment or evaluated.conditions != conditions
                                 or any(a.partition.provenance != backend.provenance for a in evaluated.analyses)):
        raise ValueError("Cached design evaluation belongs to a different candidate or model")
    strands = {s.name: s for s in assignment.strands}
    trigger, switch = strands["Trigger"], strands["Switch"]
    layout = problem.layout
    frame = problem.reading_frame_report(assignment)
    if request.exclude_stop_candidates and frame["has_premature_stop"]:
        return DesignCandidate.excluded_for_stop(
            trial=trial, switch_seq=switch.sequence, trigger_seq=trigger.sequence,
            first_stop=frame["first_stop"], reporter_label=request.reporter.label)
    on = ComplexSpec((trigger, switch))
    if request.exclude_non_linear_candidates:
        screen = backend.screen_unpaired_region(
            on, conditions, start=len(trigger.sequence)+layout["rbs_linker_start"],
            end=len(trigger.sequence)+layout["rbs_linker_end"])
        if not screen.all_minima_unpaired:
            if screen.witness is None:
                raise ValueError("A failing ON minimum needs its structural witness")
            return DesignCandidate.excluded_for_non_linear(
                trial=trial, switch_seq=switch.sequence, trigger_seq=trigger.sequence,
                structure_on=screen.witness.structure, delta_g_on=screen.witness.energy_kcal_mol,
                has_stop_codon=bool(frame["has_premature_stop"]), first_stop=frame["first_stop"],
                reporter_label=request.reporter.label)
    minimum_on = backend.minimum(on, conditions)
    evaluated = evaluated or evaluate_design(problem, assignment, backend.analyze,
                                              standard_state=backend.standard_state)
    calculations = _CandidateCalculations(backend, conditions, evaluated)
    off = ComplexSpec((switch,))
    free_trigger = ComplexSpec((trigger,))
    off_pairs = calculations.analysis(off, need_pairs=True).pairs
    on_pairs = calculations.analysis(on, need_pairs=True).pairs
    if off_pairs is None or on_pairs is None:
        raise ValueError("OFF and ON require nonempty connected ensembles")
    reporter_length = len(normalize_rna(request.reporter.sequence_after_start_rna, field_name="reporter"))
    regional = functional_defect_metrics(request.architecture, reporter_length,
                                          off_pairs.matrix, on_pairs.matrix,
                                          weights=request.functional_defect_weights)
    breakdown = defect_breakdown(evaluated.defect)
    minimum_off = backend.minimum(off, conditions)
    minimum_trigger = backend.minimum(free_trigger, conditions)
    trigger_linearity = {}
    if request.inverse_extension is not None:
        pairs = calculations.analysis(free_trigger, need_pairs=True).pairs
        if pairs is None:
            raise ValueError("An isolated trigger must have a nonempty ensemble")
        per_base = [pairs.matrix[i][i] for i in range(len(trigger.sequence))]
        unpaired = math.fsum(per_base) / len(per_base)
        trigger_linearity = {
            "normalized_defect": 1-unpaired, "unpaired_fraction": unpaired,
            "per_base_unpaired": per_base, "structure_mfe": minimum_trigger.structure,
            "nucleotide_count": len(per_base), "expected_paired_nt": (1-unpaired)*len(per_base),
            "model": {**asdict(backend.provenance), "conditions": asdict(conditions)},
        }
    ddg = minimum_on.energy_kcal_mol-minimum_off.energy_kcal_mol-minimum_trigger.energy_kcal_mol
    region = switch.sequence[layout["rbs_linker_start"]:layout["rbs_linker_end"]]
    rbs_energy = backend.minimum(ComplexSpec((Strand(region, switch.material, "RBS_linker"),)),
                                 conditions).energy_kcal_mol
    accessibility = unpaired_window_accessibility(
        backend, on, conditions,
        start=len(trigger.sequence) + layout["rbs_prefix_start"],
        end=len(trigger.sequence) + layout["aug_start"] + 3,
        partition=calculations.analysis(on).partition)
    trigger_total = request.thermo.resolved_trigger_concentration_m
    switch_total = request.thermo.resolved_switch_concentration_m
    max_size = request.thermo.validation_max_complex_size
    on_yield = calculations.complex_percent((trigger, switch), (trigger_total, switch_total),
                                            (trigger, switch), switch_total, max_size)
    aptamer = Strand(request.sequences.aptamer_dna, "dna", "Aptamer")
    leak = calculations.complex_percent((aptamer, trigger, switch),
                                        (trigger_total, trigger_total, switch_total),
                                        (trigger, switch), switch_total, max_size)
    rbs_prediction = _rbs_measurement(request, problem, assignment)
    raw = compute_selection_score(on_yield_pct=on_yield, leak_pct=leak, ddg_activation=ddg,
                                  delta_g_rbs_linker=rbs_energy, weights=request.scoring,
                                  functional_defect=regional["balanced"],
                                  trigger_linearity_defect=trigger_linearity.get("normalized_defect"),
                                  rbs_calculator_merit=rbs_prediction.get("merit"))
    return DesignCandidate(
        trial=trial, switch_seq=switch.sequence, trigger_seq=trigger.sequence,
        score=selection_score_100(raw), defect=round(evaluated.loss, 6),
        on_yield_pct=on_yield, leak_pct=leak,
        delta_g_off=round(minimum_off.energy_kcal_mol, 6),
        delta_g_on=round(minimum_on.energy_kcal_mol, 6), ddg_activation=round(ddg, 6),
        delta_g_rbs_linker=round(rbs_energy, 6), has_bad_rbs_linker=rbs_energy < -8,
        has_stop_codon=bool(frame["has_premature_stop"]), first_stop=frame["first_stop"],
        structure_off=minimum_off.structure, structure_on=minimum_on.structure,
        reporter_label=request.reporter.label, defect_structural=breakdown["structural"],
        defect_concentration=breakdown["concentration"], defect_breakdown=breakdown,
        regional_defects=regional, trigger_linearity=trigger_linearity,
        rbs_accessibility=accessibility.to_dict(), rbs_calculator=rbs_prediction, score_raw=raw)


def design_candidates(request: RunRequest, backend, *, options: TrialOptions = TrialOptions(),
                      progress_callback: Callable | None = None,
                      cancel_check: Callable[[], bool] | None = None,
                      start_trial: int = 1, trial_callback: Callable | None = None,
                      control_check: Callable | None = None) -> NativeDesignBatch:
    """Diverse local-first trials, then globally refined viable candidates.

    Trial k uses seed ``options.seed + k - 1``. Thermodynamic failures propagate;
    they are never converted into bad scores or empty successful runs. The
    optimizer's finite budget and stop reason accompany every evaluated trial.
    """
    request = request.normalized()
    if type(start_trial) is not int or not 1 <= start_trial <= request.trials+1:
        raise ValueError('The first design trial is outside the requested run.')
    problem = build_design_problem(request)
    started = datetime.now().astimezone().isoformat()
    candidates, reports = [], []
    analysis_cache = ComplexAnalysisCache(backend)
    cancelled = False
    controlled_status = None
    for trial in range(start_trial, request.trials+1):
        if cancel_check and cancel_check():
            cancelled = True
            break
        # Soft controls never interrupt a thermodynamic call or discard a trial.
        from ..design_checkpoint import control_action
        action = control_action(control_check)
        if action is not None:
            controlled_status = 'paused' if action == 'pause' else 'stopped'
            break
        evaluations = _EvaluationCache(max_entries=max(8, 2 * options.lineages),
                                       max_pair_cells=analysis_cache.max_pair_cells)

        def loss(assignment):
            # Fail with a clear OSTIR diagnostic before the expensive tube
            # calculation when a custom RBS has no evaluable initiation site.
            rbs_prediction = _rbs_measurement(request, problem, assignment)
            # Short RNA-only calculation first; the bounded run cache reuses an
            # unchanged fragment across mutations and final candidate reporting.
            rbs_minimum = None
            if request.scoring.rbs > 0:
                sequence = assignment.sequence("Switch")
                region = sequence[problem.layout["rbs_linker_start"]:problem.layout["rbs_linker_end"]]
                rbs_minimum = analysis_cache.minimum(
                    ComplexSpec((Strand(region, "rna", "RBS_linker"),)), problem.conditions)
            evaluation = evaluate_design(problem, assignment, analysis_cache,
                                         standard_state=backend.standard_state)
            evaluations.remember(evaluation)
            trigger_pairs = None
            if request.inverse_extension is not None and request.scoring.trigger_linearity > 0:
                trigger = next(strand for strand in assignment.strands if strand.name == "Trigger")
                trigger_pairs = analysis_cache(ComplexSpec((trigger,)), problem.conditions, need_pairs=True).pairs
            return objective_from_evaluation(request, evaluation, rbs_minimum=rbs_minimum,
                rbs_prediction=rbs_prediction, trigger_pairs=trigger_pairs)

        def progress(state):
            if progress_callback:
                progress_callback(trial, request.trials, state)

        result = optimize(problem, loss, options=replace(options, seed=options.seed+trial-1),
                          progress_callback=progress, cancel_check=cancel_check)
        if result.stop_reason == "cancelled":
            reports.append(result)
            cancelled = True
            break
        if result.best_assignment is None:
            if result.stop_reason != "no_admissible_sequence":
                raise ValueError("A finished optimization did not produce a sequence")
            reports.append(result)
            if trial_callback:
                trial_callback(trial, None, result, analysis_cache.statistics())
            continue
        finalists = (ScoredAssignment(result.best_assignment, result.best_objective),) + tuple(
            elite for elite in result.elites if elite.assignment != result.best_assignment)
        first_candidate = chosen = None
        checked = 0
        screening_complete = False
        for elite in finalists:
            if cancel_check and cancel_check():
                cancelled = True
                break
            # A cache miss remains lazy: existing STOP/ON screens run before
            # evaluate_candidate reconstructs the full objective, if needed.
            evaluation = evaluations.get(elite.assignment)
            candidate = evaluate_candidate(request, problem, elite.assignment, analysis_cache,
                                           trial=trial, evaluated=evaluation)
            checked += 1
            if first_candidate is None:
                first_candidate, chosen = candidate, finalists[0]
            if not candidate.excluded:
                first_candidate, chosen = candidate, elite
                screening_complete = True
                break
        else:
            screening_complete = True
        result = replace(result, candidate_screened=checked)
        if screening_complete and first_candidate is not None:
            candidates.append(first_candidate)
            result = replace(result, selected_assignment=chosen.assignment,
                             selected_objective=chosen.objective)
        reports.append(result)
        # An excluded finalist is provisional while other finalists remain.
        # A hard interruption must resume this entire trial with the same seed.
        if screening_complete and first_candidate is not None and trial_callback:
            trial_callback(trial, first_candidate, result, analysis_cache.statistics())
        if cancelled or (cancel_check and cancel_check()):
            cancelled = True
            break
    return NativeDesignBatch(tuple(sort_candidates(candidates, defect_weight=request.scoring.defect)),
                              tuple(reports), backend.provenance, problem.conditions,
                              "cancelled" if cancelled else controlled_status or "completed",
                              started, datetime.now().astimezone().isoformat(), analysis_cache.statistics())
