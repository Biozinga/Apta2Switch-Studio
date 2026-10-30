"""Bounded population search: functional regions first, global/RBS refinement second.

This is a sequence optimization heuristic, not an ensemble sampler. Its fixed
nonregression guards are relative to an evaluated starting sequence, not
experimentally calibrated thresholds. Annealing temperatures are objective
units and have no relation to the physical folding temperature.
"""
from __future__ import annotations

from collections import OrderedDict
from dataclasses import dataclass
import math
from numbers import Real
import random
from typing import Callable, Literal

from .design_constraints import BaseSite, DesignAssignment, DesignProblem, NoMutableSitesError


StopReason = Literal["target_reached", "budget_exhausted", "no_mutable_sites", "cancelled", "no_admissible_sequence", "space_exhausted"]
BudgetKind = Literal["evaluations", "proposals"]
Phase = Literal["initialization", "local", "global"]
AssignmentKey = tuple[tuple[str, str, str], ...]


def _nonnegative_real(value, name: str) -> float:
    if isinstance(value, bool) or not isinstance(value, Real):
        raise ValueError(f"{name} must be a finite nonnegative real number")
    try:
        result = float(value)
    except (ValueError, OverflowError) as error:
        raise ValueError(f"{name} must be a finite nonnegative real number") from error
    if not math.isfinite(result) or result < 0:
        raise ValueError(f"{name} must be a finite nonnegative real number")
    return result


def _fraction(value, name: str) -> float:
    result = _nonnegative_real(value, name)
    if result > 1:
        raise ValueError(f"{name} must not exceed one")
    return result


def rbs_stability_settings(weight: float) -> dict:
    """Saved method identity; resuming cannot silently change the objective."""
    return {"scheme": "rbs_refinement_v1", "weight": _nonnegative_real(weight, "rbs_weight"),
            "refinement_objective": "global_defect + weight * max(0, -dG_RBS_linker) / 100",
            "scope": "isolated RNA RBS-prefix through linker; fixed local and tube guards",
            "joint_accessibility": "informative connected ON RBS-prefix through AUG; not added to score"}


def translation_search_settings(scoring, inverse_extension: bool) -> dict:
    """Freeze the additional refinement method and weights in each new run."""
    return {"scheme": "ostir_occupancy_and_trigger_refinement_v1",
            "rbs_method": "ostir_ideal_occupancy_v1", "ostir_version": "1.1.3", "vienna_version": "2.7.2",
            "rbs_calculator_weight": _nonnegative_real(scoring.rbs_calculator, "rbs_calculator_weight"),
            "trigger_linearity_weight": _nonnegative_real(scoring.trigger_linearity, "trigger_linearity_weight") if inverse_extension else 0.,
            "refinement_addition": "w_RBS_calculator * (1 - merit) + w_trigger * paired_fraction",
            "local_guards": "unchanged; local structural defects have priority"}


@dataclass(frozen=True, slots=True)
class ObjectiveValue:
    """Normalized defects and optional nucleotide-level mutation guidance."""
    local_defects: tuple[float, float, float, float]
    local_weights: tuple[float, float, float, float]
    global_defect: float
    guard_defects: tuple[float, float, float]
    site_weights: tuple[tuple[BaseSite, float], ...] = ()
    global_site_weights: tuple[tuple[BaseSite, float], ...] = ()
    rbs_linker_energy: float = 0.0
    rbs_weight: float = 0.0
    rbs_calculator_merit: float = 0.0
    rbs_calculator_weight: float = 0.0
    trigger_linearity_defect: float = 0.0
    trigger_linearity_weight: float = 0.0

    def __post_init__(self):
        for name, size, validate in (("local_defects", 4, _fraction), ("local_weights", 4, _nonnegative_real),
                                     ("guard_defects", 3, _fraction)):
            values = tuple(getattr(self, name))
            if len(values) != size:
                raise ValueError(f"{name} must contain {size} values")
            object.__setattr__(self, name, tuple(validate(value, name) for value in values))
        if not any(self.local_weights):
            raise ValueError("At least one local weight must be positive")
        object.__setattr__(self, "global_defect", _fraction(self.global_defect, "global_defect"))
        energy = self.rbs_linker_energy
        if isinstance(energy, bool) or not isinstance(energy, Real):
            raise ValueError("RBS–linker energy must be finite, in kcal/mol")
        try:
            energy = float(energy)
        except (ValueError, OverflowError) as error:
            raise ValueError("RBS–linker energy must be finite, in kcal/mol") from error
        if not math.isfinite(energy):
            raise ValueError("RBS–linker energy must be finite, in kcal/mol")
        object.__setattr__(self, "rbs_linker_energy", energy)
        object.__setattr__(self, "rbs_weight", _nonnegative_real(self.rbs_weight, "rbs_weight"))
        for name in ("rbs_calculator_merit", "trigger_linearity_defect"):
            object.__setattr__(self, name, _fraction(getattr(self, name), name))
        for name in ("rbs_calculator_weight", "trigger_linearity_weight"):
            object.__setattr__(self, name, _nonnegative_real(getattr(self, name), name))
        if not math.isfinite(self.refinement_loss):
            raise ValueError("The weighted RBS stability penalty must be finite")
        for name in ("site_weights", "global_site_weights"):
            rows = tuple(getattr(self, name))
            if any(not isinstance(site, BaseSite) for site, _ in rows):
                raise ValueError(f"{name} requires explicit BaseSite objects")
            if len({site for site, _ in rows}) != len(rows):
                raise ValueError(f"{name} cannot repeat a site")
            object.__setattr__(self, name, tuple((site, _nonnegative_real(weight, name)) for site, weight in rows))

    @property
    def local_loss(self) -> float:
        # Scaling avoids overflow for valid but unusually large relative weights.
        scale = max(self.local_weights)
        weights = tuple(weight / scale for weight in self.local_weights)
        return math.fsum(value * weight for value, weight in zip(self.local_defects, weights)) / math.fsum(weights)

    @property
    def rbs_stability_penalty(self) -> float:
        """Continuous cost; division by 100 matches percentage-score units.

        This is a search objective, not a probability or a calibrated cutoff.
        A zero existing RBS score weight also disables this search preference.
        """
        return self.rbs_weight * (max(0.0, -self.rbs_linker_energy) / 100.0)

    @property
    def refinement_loss(self) -> float:
        return (self.global_defect + self.rbs_stability_penalty
                + self.rbs_calculator_weight * (1.0 - self.rbs_calculator_merit)
                + self.trigger_linearity_weight * self.trigger_linearity_defect)


@dataclass(frozen=True, slots=True)
class ScoredAssignment:
    assignment: DesignAssignment
    objective: ObjectiveValue


@dataclass(frozen=True, slots=True)
class TrialOptions:
    """Total evaluation/proposal budgets shared by initialization and lineages.

    Every attempted generation or mutation spends one proposal, even when it
    fails the sequence screen or revisits the cache. Initialization uses at most
    20% of either budget (at least one); the remaining work is split into local
    and global phases. Small budgets adapt rather than allocating hidden work.
    The global target applies to the refinement cost, including the optional
    RBS stability penalty; the separately reported global defect stays physical.
    """
    seed: int = 0
    max_evaluations: int = 1000
    max_proposals: int = 5000
    cache_size: int = 1024
    population_size: int = 64
    lineages: int = 8
    local_fraction: float = .7
    local_tolerance: float = .005
    guard_tolerance: float = .01
    exploration_fraction: float = .2
    initial_temperature: float = .02
    final_temperature: float = .0005
    target_local_loss: float = .005
    target_global_loss: float = .005

    def __post_init__(self):
        if type(self.seed) is not int:
            raise ValueError("seed must be an explicit integer")
        for name, minimum in (("max_evaluations", 1), ("max_proposals", 1), ("cache_size", 0),
                              ("population_size", 1), ("lineages", 1)):
            if type(getattr(self, name)) is not int or getattr(self, name) < minimum:
                raise ValueError(f"{name} must be an integer of at least {minimum}")
        for name in ("local_fraction", "local_tolerance", "guard_tolerance", "exploration_fraction",
                     "target_local_loss", "target_global_loss"):
            object.__setattr__(self, name, _fraction(getattr(self, name), name))
        if not 0 < self.local_fraction < 1:
            raise ValueError("local_fraction must lie strictly between zero and one")
        for name in ("initial_temperature", "final_temperature"):
            object.__setattr__(self, name, _nonnegative_real(getattr(self, name), name))
        if self.final_temperature > self.initial_temperature:
            raise ValueError("final_temperature cannot exceed initial_temperature")


@dataclass(frozen=True, slots=True)
class OptimizationProgress:
    evaluations: int
    proposals: int
    cache_hits: int
    temperature: float
    current_assignment: DesignAssignment | None
    current_loss: float | None
    best_assignment: DesignAssignment | None
    best_loss: float | None
    phase: Phase
    population_evaluated: int
    lineages_started: int
    current_objective: ObjectiveValue | None
    best_objective: ObjectiveValue | None
    best_global_loss: float | None
    elites: tuple[ScoredAssignment, ...]
    guard_limits: tuple[float, float, float] | None
    local_limits: tuple[float, float, float, float] | None


@dataclass(frozen=True, slots=True)
class OptimizationResult(OptimizationProgress):
    stop_reason: StopReason
    budget_kind: BudgetKind | None
    options: TrialOptions
    selected_assignment: DesignAssignment | None = None
    selected_objective: ObjectiveValue | None = None
    candidate_screened: int = 0


def assignment_key(assignment: DesignAssignment) -> AssignmentKey:
    """Preserve strand order, names and chemistry in the invocation-local cache."""
    if not isinstance(assignment, DesignAssignment):
        raise ValueError("An explicit DesignAssignment is required")
    return tuple((strand.name, strand.material, strand.sequence) for strand in assignment.strands)


@dataclass(slots=True)
class _Lineage:
    current: ScoredAssignment
    local_best: ScoredAssignment
    limits: tuple[float, ...] | None = None


def optimize(
    problem: DesignProblem,
    loss: Callable[[DesignAssignment], ObjectiveValue],
    *,
    options: TrialOptions = TrialOptions(),
    initial_assignment: DesignAssignment | None = None,
    progress_callback: Callable[[OptimizationProgress], None] | None = None,
    cancel_check: Callable[[], bool] | None = None,
) -> OptimizationResult:
    """Optimize a fixed objective without importing a thermodynamic engine.

    Errors propagate. Cancellation preserves completed evaluations. The bounded
    cache and random generator belong exclusively to this call; neither physical
    conditions nor objective weights may change during it. Final elites are
    candidates in one fixed local neighborhood satisfying the fixed guards.
    """
    if not isinstance(problem, DesignProblem) or not isinstance(options, TrialOptions):
        raise ValueError("A DesignProblem and TrialOptions are required")
    if not callable(loss) or any(fn is not None and not callable(fn) for fn in (progress_callback, cancel_check)):
        raise ValueError("Objective and any callbacks must be callable")
    if initial_assignment is not None:
        problem.validate(initial_assignment)
    rng = random.Random(options.seed)
    cache: OrderedDict[AssignmentKey, ObjectiveValue] = OrderedDict()
    archive: dict[AssignmentKey, ScoredAssignment] = {}
    archive_size = max(32, options.population_size, 8 * options.lineages)
    evaluations = proposals = cache_hits = population_evaluated = lineages_started = 0
    phase: Phase = "initialization"
    temperature = options.initial_temperature
    current: ScoredAssignment | None = None
    guard_limits = local_limits = None
    objective_weights = None
    lineages: list[_Lineage] = []

    def cancelled() -> bool:
        return cancel_check is not None and cancel_check()

    def violation(item: ScoredAssignment) -> float:
        return max((max(0.0, value - ceiling) for value, ceiling in zip(item.objective.guard_defects, guard_limits)),
                   default=0.0) if guard_limits is not None else 0.0

    def local_key(item: ScoredAssignment):
        return (violation(item), item.objective.local_loss, item.objective.refinement_loss, assignment_key(item.assignment))

    def ceilings(item: ScoredAssignment):
        return tuple(min(1.0, value + options.local_tolerance) if weight > 0 else 1.0
                     for value, weight in zip(item.objective.local_defects, item.objective.local_weights))

    def within(item: ScoredAssignment, limits) -> bool:
        return violation(item) == 0 and all(value <= ceiling + 1e-12
                                            for value, ceiling in zip(item.objective.local_defects, limits))

    def best_local() -> ScoredAssignment | None:
        return min(archive.values(), key=local_key) if archive else None

    def final_candidates() -> list[ScoredAssignment]:
        anchor = best_local()
        if anchor is None or guard_limits is None:
            return []
        limits = local_limits if local_limits is not None else ceilings(anchor)
        return sorted((item for item in archive.values() if within(item, limits)),
                      key=lambda item: (item.objective.refinement_loss, item.objective.local_loss,
                                        assignment_key(item.assignment)))

    def diverse(rows, count: int, *, first=None) -> list[ScoredAssignment]:
        """Quality shortlist, then farthest-first on mutable coupled groups."""
        ordered = list(rows)
        if not ordered or count <= 0:
            return []
        selected = [first or ordered[0]]
        selected_keys = {assignment_key(selected[0].assignment)}
        shortlist = ordered[:max(count * 4, count)]
        while len(selected) < min(count, len(ordered)):
            candidates = [item for item in shortlist if assignment_key(item.assignment) not in selected_keys]
            if not candidates:
                candidates = [item for item in ordered if assignment_key(item.assignment) not in selected_keys]
            winner = max(candidates, key=lambda item: min(problem.design_distance(item.assignment, old.assignment)
                                                          for old in selected))
            selected.append(winner)
            selected_keys.add(assignment_key(winner.assignment))
        return selected

    def prune_archive():
        if len(archive) <= 2 * archive_size:
            return
        ordered = sorted(archive.values(), key=local_key)
        finalists = final_candidates()
        protected = [ordered[0]] + finalists[:options.lineages]
        # Protect current/best states of every lineage even after cache eviction.
        protected += [state for lineage in lineages for state in (lineage.current, lineage.local_best)]
        kept = {assignment_key(item.assignment): item for item in protected}
        others = [item for item in ordered if assignment_key(item.assignment) not in kept]
        for item in diverse(others, max(0, archive_size - len(kept))):
            kept[assignment_key(item.assignment)] = item
        archive.clear()
        archive.update(kept)

    def choose_best() -> ScoredAssignment | None:
        finalists = final_candidates()
        return finalists[0] if phase == "global" and finalists else best_local()

    def values(*, final=False) -> dict:
        finalists = final_candidates()
        best = finalists[0] if final and finalists else choose_best()
        elite_rows = diverse(finalists, options.lineages) if finalists else []
        elite_rows.sort(key=lambda item: (item.objective.refinement_loss, item.objective.local_loss,
                                          assignment_key(item.assignment)))
        return dict(evaluations=evaluations, proposals=proposals, cache_hits=cache_hits, temperature=temperature,
                    current_assignment=current.assignment if current else None,
                    current_loss=current.objective.local_loss if current else None,
                    best_assignment=best.assignment if best else None,
                    best_loss=best.objective.local_loss if best else None,
                    phase=phase, population_evaluated=population_evaluated, lineages_started=lineages_started,
                    current_objective=current.objective if current else None,
                    best_objective=best.objective if best else None,
                    best_global_loss=best.objective.global_defect if best else None,
                    elites=tuple(elite_rows), guard_limits=guard_limits, local_limits=local_limits)

    def set_guards():
        nonlocal guard_limits
        if guard_limits is None and archive:
            anchor = min(archive.values(), key=lambda item: (max(item.objective.guard_defects),
                         item.objective.global_defect, item.objective.local_loss, assignment_key(item.assignment)))
            guard_limits = tuple(min(1.0, value + options.guard_tolerance) for value in anchor.objective.guard_defects)

    def finish(reason: StopReason, budget: BudgetKind | None = None) -> OptimizationResult:
        nonlocal local_limits
        set_guards()
        if archive and local_limits is None:
            local_limits = ceilings(best_local())
        return OptimizationResult(**values(final=True), stop_reason=reason, budget_kind=budget, options=options)

    def at_target(item: ScoredAssignment | None) -> bool:
        return item is not None and guard_limits is not None and violation(item) == 0 and (
            item.objective.local_loss <= options.target_local_loss and item.objective.refinement_loss <= options.target_global_loss)

    def stop(*, target=True) -> OptimizationResult | None:
        if cancelled():
            return finish("cancelled")
        if target and at_target(next(iter(final_candidates()), choose_best())):
            return finish("target_reached")
        if evaluations >= options.max_evaluations:
            return finish("budget_exhausted", "evaluations")
        if proposals >= options.max_proposals:
            return finish("budget_exhausted", "proposals")
        return None

    def progress():
        if progress_callback is not None:
            progress_callback(OptimizationProgress(**values()))

    def admissible(assignment):
        return problem.screening is None or problem.passes_sequence_screen(assignment)

    def evaluate(assignment: DesignAssignment) -> ScoredAssignment:
        nonlocal evaluations, cache_hits, objective_weights
        key = assignment_key(assignment)
        if key in cache:
            cache_hits += 1
            objective = cache.pop(key)
            cache[key] = objective
        else:
            objective = loss(assignment)
            if not isinstance(objective, ObjectiveValue):
                raise ValueError("The objective must return an ObjectiveValue")
            weights = (objective.local_weights, objective.rbs_weight,
                       objective.rbs_calculator_weight, objective.trigger_linearity_weight)
            if objective_weights is not None and weights != objective_weights:
                raise ValueError("Objective weights cannot change during a trial")
            objective_weights = weights
            evaluations += 1
            if options.cache_size:
                cache[key] = objective
                if len(cache) > options.cache_size:
                    cache.popitem(last=False)
        item = ScoredAssignment(assignment, objective)
        archive[key] = item
        return item

    if cancelled():
        return finish("cancelled")
    pool_goal = min(options.population_size, max(1, options.max_evaluations // 5))
    generation_budget = max(1, options.max_proposals // 5)
    # Exact enumeration removes duplicate-generation stalls in genuinely small
    # spaces, without materializing an exponential space just to initialize.
    enumerated = problem.enumerate_assignments(pool_goal)
    pending = list(enumerated) if enumerated is not None else None
    if pending is not None:
        rng.shuffle(pending)
        if initial_assignment is not None:
            pending = [item for item in pending if assignment_key(item) != assignment_key(initial_assignment)]
        pool_goal = min(pool_goal, len(enumerated))
    seen: set[AssignmentKey] = set()
    while population_evaluated < pool_goal and proposals < generation_budget:
        if cancelled():
            return finish("cancelled")
        if initial_assignment is not None and proposals == 0:
            candidate = initial_assignment
        elif pending is not None:
            if not pending:
                break
            candidate = pending.pop()
        else:
            candidate = problem.generate(rng)
        proposals += 1
        key = assignment_key(candidate)
        if key in seen:
            continue
        seen.add(key)
        if cancelled():
            return finish("cancelled")
        if not admissible(candidate):
            continue
        current = evaluate(candidate)
        population_evaluated += 1
        progress()
    set_guards()
    if not archive:
        return finish("cancelled" if cancelled() else "no_admissible_sequence")
    if cancelled():
        return finish("cancelled")
    if at_target(next(iter(final_candidates()), choose_best())):
        return finish("target_reached")
    if not problem.mutable_components():
        return finish("no_mutable_sites")
    if enumerated is not None and all(assignment_key(item) in seen for item in enumerated):
        return finish("space_exhausted")
    stopped = stop()
    if stopped is not None:
        return stopped

    seeds = diverse(sorted(archive.values(), key=local_key), min(options.lineages, population_evaluated))
    lineages = [_Lineage(item, item) for item in seeds]
    lineages_started = len(lineages)
    remaining_evaluations = options.max_evaluations - evaluations
    remaining_proposals = options.max_proposals - proposals
    local_evaluations = max(1, int(remaining_evaluations * options.local_fraction))
    local_proposals = max(1, int(remaining_proposals * options.local_fraction))

    for phase in ("local", "global"):
        phase_start_evaluations, phase_start_proposals = evaluations, proposals
        evaluation_quota = local_evaluations if phase == "local" else options.max_evaluations - evaluations
        proposal_quota = local_proposals if phase == "local" else options.max_proposals - proposals
        if phase == "global":
            local_limits = ceilings(best_local())
            # Refinement only spends work on states that can contribute to the
            # final neighborhood. Each lineage additionally keeps its own fixed
            # regional ceilings, intersected with that shared neighborhood.
            global_seeds = diverse(final_candidates(), lineages_started)
            lineages = [_Lineage(item, item, tuple(min(a, b) for a, b in zip(
                ceilings(item), local_limits))) for item in global_seeds]
            current = choose_best()
        progress()
        turn = 0
        while (evaluations - phase_start_evaluations < evaluation_quota
               and proposals - phase_start_proposals < proposal_quota):
            stopped = stop()
            if stopped is not None:
                return stopped
            lineage = lineages[turn % len(lineages)]
            turn += 1
            current = lineage.current
            fraction = max((evaluations - phase_start_evaluations) / max(1, evaluation_quota - 1),
                           (proposals - phase_start_proposals) / max(1, proposal_quota - 1))
            fraction = min(1.0, fraction)
            temperature = (1.0 - fraction) * options.initial_temperature + fraction * options.final_temperature
            weights = current.objective.site_weights if phase == "local" else current.objective.global_site_weights
            proposals += 1
            try:
                candidate = problem.mutate_guided(current.assignment, rng, weights,
                                                 exploration_fraction=options.exploration_fraction)
            except NoMutableSitesError:
                # A chosen generic constraint component can be immovable even
                # when another is not; count the failed proposal and try again.
                continue
            if cancelled():
                return finish("cancelled")
            if not admissible(candidate):
                continue
            item = evaluate(candidate)
            if phase == "local":
                if local_key(item) < local_key(lineage.local_best):
                    lineage.local_best = item
                previous_violation, next_violation = violation(current), violation(item)
                allowed = next_violation <= previous_violation
                delta = item.objective.local_loss - current.objective.local_loss
                if delta == 0:
                    delta = item.objective.refinement_loss - current.objective.refinement_loss
                improves_guard = next_violation < previous_violation
            else:
                allowed = within(item, lineage.limits)
                delta = item.objective.refinement_loss - current.objective.refinement_loss
                improves_guard = False
            if allowed and (improves_guard or delta <= 0 or (temperature > 0 and
                    math.log(max(rng.random(), 1e-300)) < -delta / temperature)):
                lineage.current = current = item
            prune_archive()
            progress()
        stopped = stop()
        if stopped is not None:
            return stopped
    raise RuntimeError("Optimizer finished without consuming its shared budget")
