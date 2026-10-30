"""Evaluate a constrained sequence assignment against its complete target tubes.

The global diagnostic is the mean normalized tube ensemble defect. The design
objective adds functional-region priorities using these same matrices. Chemistry,
connected partition functions and their standard
state come from an explicit backend; this layer never selects a substitute
energy model. No application state or persistent cache is modified.
"""
from __future__ import annotations

from dataclasses import dataclass
from collections import OrderedDict
import math
from typing import Protocol

from .defects import DesignDefect, TubeDefect, TubeTarget, tube_defect
from .design_constraints import DesignAssignment, DesignProblem
from .structures import PairingConstraints
from .tubes import StandardStateConvention, TubeResult, solve_tube
from .types import (
    ComplexSpec, Conditions, PairProbabilities, PartitionResult, enumerate_complexes,
)


@dataclass(frozen=True, slots=True)
class ComplexAnalysis:
    partition: PartitionResult
    pairs: PairProbabilities | None = None

    def __post_init__(self):
        if not isinstance(self.partition, PartitionResult):
            raise ValueError("Analysis requires an explicit PartitionResult")
        if self.pairs is not None and (
            not isinstance(self.pairs, PairProbabilities)
            or self.pairs.partition != self.partition
        ):
            raise ValueError("Pair probabilities must belong to this exact partition")


class AnalyzeComplex(Protocol):
    def __call__(self, spec: ComplexSpec, conditions: Conditions, *,
                 need_pairs: bool) -> ComplexAnalysis: ...


class ComplexAnalysisCache:
    """Bounded run-local reuse of exact partitions, analyses, minima and screens.

    Neither mutations of another strand nor another tube require refolding an
    unchanged complex. Ordered strands, labels, chemistry and conditions all
    belong to the key; the cache is tied to one explicit backend provenance.
    Scalar/structure observables share the entry budget with analyses. Native
    handles are never retained; only completed result objects are reused.
    """
    def __init__(self, backend, *, max_entries=256, max_pair_cells=2_000_000):
        if type(max_entries) is not int or max_entries < 0 or type(max_pair_cells) is not int or max_pair_cells < 0:
            raise ValueError("Cache capacities must be nonnegative integers")
        self.backend, self.provenance = backend, backend.provenance
        self.max_entries, self.max_pair_cells = max_entries, max_pair_cells
        self._cache = OrderedDict()
        self._pair_cells = self.calls = self.hits = 0
        self._calls_by_kind = {kind: 0 for kind in ("analysis", "partition", "minimum", "screen")}

    @property
    def standard_state(self):
        return self.backend.standard_state

    def _check_provenance(self):
        if self.backend.provenance != self.provenance:
            raise ValueError("The thermodynamic backend changed during optimization")

    def _remember(self, key, result, cells=0):
        previous = self._cache.pop(key, None)
        if previous is not None:
            self._pair_cells -= previous[1]
        if self.max_entries and cells <= self.max_pair_cells:
            self._cache[key] = (result, cells)
            self._pair_cells += cells
            while len(self._cache) > self.max_entries or self._pair_cells > self.max_pair_cells:
                _, (_, evicted_cells) = self._cache.popitem(last=False)
                self._pair_cells -= evicted_cells

    def _observable(self, kind, spec, conditions, parameters, calculate, *, constraints=None):
        self._check_provenance()
        key = (kind, spec, conditions, parameters)
        cached = self._cache.get(key)
        if cached is not None:
            self._cache.move_to_end(key)
            self.hits += 1
            return cached[0]
        result = calculate()
        self.calls += 1
        self._calls_by_kind[kind] += 1
        context = getattr(result, "context", None)
        if context is not None and (
                context.complex_spec != spec or context.conditions != conditions
                or context.provenance != self.provenance or context.constraints != constraints):
            raise ValueError("A cached observable must preserve the exact unrestricted physical context")
        if kind == "partition" and (
                not isinstance(result, PartitionResult) or result.complex_spec != spec
                or result.conditions != conditions or result.provenance != self.provenance
                or result.constraints != constraints):
            raise ValueError("A cached partition must preserve the exact requested physical context")
        self._remember(key, result)
        return result

    def partition(self, spec, conditions, constraints=None):
        """Reuse an exact constrained PF, or the existing unrestricted analysis."""
        if constraints is not None:
            if not isinstance(constraints, PairingConstraints) or len(constraints.partners) != spec.n_nucleotides:
                raise ValueError("Constraints must describe every nucleotide of the complex")
            if all(partner == -1 for partner in constraints.partners):
                constraints = None
        if constraints is None:
            return self(spec, conditions, need_pairs=False).partition
        return self._observable("partition", spec, conditions, (constraints,),
                                lambda: self.backend.partition(spec, conditions, constraints),
                                constraints=constraints)

    def minimum(self, spec, conditions):
        """Reuse an unrestricted minimum under the selected backend convention."""
        return self._observable("minimum", spec, conditions, (),
                                lambda: self.backend.minimum(spec, conditions))

    def screen_unpaired_region(self, spec, conditions, *, start, end):
        """Reuse the exact default-tolerance unrestricted region screen."""
        if type(start) is not int or type(end) is not int or not 0 <= start < end <= spec.n_nucleotides:
            raise ValueError("Region must be a nonempty half-open interval in complex coordinates")
        return self._observable("screen", spec, conditions, (start, end),
                                lambda: self.backend.screen_unpaired_region(
                                    spec, conditions, start=start, end=end))

    def __call__(self, spec, conditions, *, need_pairs):
        if type(need_pairs) is not bool:
            raise ValueError("need_pairs must be bool")
        self._check_provenance()
        key = ("analysis", spec, conditions)
        cached = self._cache.get(key)
        result = cached[0] if cached is not None else None
        if result is not None and (not need_pairs or result.partition.is_empty or result.pairs is not None):
            self._cache.move_to_end(key)
            self.hits += 1
            return result
        result = self.backend.analyze(spec, conditions, need_pairs=need_pairs)
        self.calls += 1
        self._calls_by_kind["analysis"] += 1
        if (not isinstance(result, ComplexAnalysis) or result.partition.complex_spec != spec
                or result.partition.conditions != conditions or result.partition.provenance != self.provenance
                or result.partition.constraints is not None
                or (need_pairs and not result.partition.is_empty and result.pairs is None)):
            raise ValueError("A cached analysis must preserve the exact unrestricted physical context")
        cells = spec.n_nucleotides ** 2 if result.pairs is not None else 0
        self._remember(key, result, cells)
        return result

    analyze = __call__

    def statistics(self):
        return {"backend_calls": self.calls, "cache_hits": self.hits,
                "analysis_calls": self._calls_by_kind["analysis"],
                "partition_calls": self._calls_by_kind["partition"],
                "minimum_calls": self._calls_by_kind["minimum"],
                "screen_calls": self._calls_by_kind["screen"],
                "cached_complexes": sum(key[0] == "analysis" for key in self._cache),
                "cached_observables": len(self._cache), "cached_pair_cells": self._pair_cells,
                "max_entries": self.max_entries, "max_pair_cells": self.max_pair_cells}


@dataclass(frozen=True, slots=True)
class EvaluatedTube:
    name: str
    equilibrium: TubeResult
    targets: tuple[TubeTarget, ...]
    defect: TubeDefect


@dataclass(frozen=True, slots=True)
class DesignEvaluation:
    problem: DesignProblem
    assignment: DesignAssignment
    conditions: Conditions
    analyses: tuple[ComplexAnalysis, ...]
    tubes: tuple[EvaluatedTube, ...]
    defect: DesignDefect

    @property
    def loss(self) -> float:
        return self.defect.normalized


def evaluate_design(problem: DesignProblem, assignment: DesignAssignment,
                    analyze: AnalyzeComplex, *,
                    standard_state: StandardStateConvention) -> DesignEvaluation:
    """Calculate OFF/ON losses including every requested off-target complex.

    Strand totals are derived from on-target stoichiometries and concentrations.
    Each tube contains all cyclic orders up to its maximum size, plus any larger
    explicitly requested target. Monomers are always included. Matrix
    coordinates follow each target's supplied strand order, even when it differs
    from the enumeration's canonical order.

    The per-call cache shares identical ordered requests across tubes. Different
    rotations are evaluated separately when targets need different coordinates;
    no unrotated probability matrix is silently reused. Complete unconstrained
    matrices are requested only for on-target complexes. Empty target ensembles
    have a concentration shortage and no fictitious pairing distribution.
    """
    if not isinstance(problem, DesignProblem):
        raise ValueError("problem must be a DesignProblem")
    problem.validate(assignment)
    if not isinstance(problem.conditions, Conditions) or not problem.tubes:
        raise ValueError("Thermodynamic evaluation requires conditions and target tubes")
    if not callable(analyze):
        raise ValueError("analyze must be an explicit complex-analysis callable")
    if not isinstance(standard_state, StandardStateConvention):
        raise ValueError("An explicit standard-state convention is required")
    strands = {s.name: s for s in assignment.strands}
    # Alias concentrations must be merged deliberately before this interface.
    if len({s.chemical_identity for s in assignment.strands}) != len(strands):
        raise ValueError("Distinct strand labels with identical chemistry require an explicit alias merge")
    targets = {t.name: t for t in problem.targets}
    prepared, requested_pairs = [], set()
    for tube in problem.tubes:
        on_targets = tuple(TubeTarget(
            name, ComplexSpec(tuple(strands[n] for n in targets[name].strands)),
            targets[name].structure, concentration,
        ) for name, concentration in tube.on_targets)
        preferred = {}
        quantities = {}
        for target in on_targets:
            spec = target.complex_spec
            if spec.cyclic_key in preferred:
                raise ValueError("A tube cannot request two target structures for the same complex")
            preferred[spec.cyclic_key] = spec
            requested_pairs.add(spec.ordered_key)
            for strand in spec.strands:
                quantities.setdefault(strand.name, []).append(target.concentration_molar)
        species = tuple(s for s in assignment.strands if s.name in quantities)
        totals = tuple(math.fsum(quantities[s.name]) for s in species)
        catalog = {spec.cyclic_key: spec
                   for spec in enumerate_complexes(species, tube.max_complex_size)}
        catalog.update(preferred)
        prepared.append((tube.name, on_targets, species, totals, tuple(catalog.values())))

    cache = {}
    reference = None
    evaluated = []
    for name, on_targets, species, totals, catalog in prepared:
        records = []
        for spec in catalog:
            key = spec.ordered_key
            if key not in cache:
                result = analyze(spec, problem.conditions, need_pairs=key in requested_pairs)
                if not isinstance(result, ComplexAnalysis):
                    raise ValueError("Backend must return a ComplexAnalysis")
                partition = result.partition
                if partition.complex_spec != spec or partition.conditions != problem.conditions:
                    raise ValueError("Backend changed strand coordinates, chemistry, labels or conditions")
                if partition.constraints is not None:
                    raise ValueError("Design evaluation requires unrestricted connected ensembles")
                if reference is not None and (
                    partition.provenance != reference.provenance
                    or partition.thermal_energy_kcal_mol != reference.thermal_energy_kcal_mol
                ):
                    raise ValueError("All design tubes require the same model provenance and RT convention")
                reference = partition
                if key in requested_pairs and not partition.is_empty and result.pairs is None:
                    raise ValueError("A nonempty target ensemble requires full pair probabilities")
                cache[key] = result
            records.append(cache[key])
        equilibrium = solve_tube(species, totals, (r.partition for r in records),
                                 standard_state=standard_state)
        defect = tube_defect(name, equilibrium, on_targets,
                             (r.pairs for r in records if r.pairs is not None))
        evaluated.append(EvaluatedTube(name, equilibrium, on_targets, defect))
    design_defect = DesignDefect(tuple(t.defect for t in evaluated))
    return DesignEvaluation(problem, assignment, problem.conditions, tuple(cache.values()),
                            tuple(evaluated), design_defect)
