"""Engine-independent complex and tube ensemble defects.

Definitions follow the complex and test-tube ensemble defects described by
Wolfe et al. (2015); their implementation and result records are independent.
Per-base errors concern a specified pairing partner, not just paired/unpaired
status. Off-target complexes affect the shortage of on-target complexes and
must not be counted a second time as an additional structural penalty.
"""
from __future__ import annotations

from dataclasses import dataclass
import math
from numbers import Real
from typing import Iterable

from .structures import format_structure, is_connected, parse_structure
from .types import ComplexSpec, PairProbabilities
from .tubes import TubeResult


@dataclass(frozen=True, slots=True)
class ComplexDefect:
    complex_spec: ComplexSpec
    target_structure: str
    per_base: tuple[float, ...]

    @property
    def expected_incorrect_nucleotides(self) -> float:
        return math.fsum(self.per_base)

    @property
    def normalized(self) -> float:
        return self.expected_incorrect_nucleotides / len(self.per_base)


def _target(spec: ComplexSpec, structure: str):
    lengths = tuple(len(s.sequence) for s in spec.strands)
    partners = parse_structure(structure, lengths)
    if not is_connected(partners, lengths):
        raise ValueError("A complex target must be a connected structure")
    return partners, format_structure(partners, lengths)


def complex_defect(pairs: PairProbabilities, target_structure: str) -> ComplexDefect:
    """Expected incorrect nucleotides in this connected complex ensemble."""
    if not isinstance(pairs, PairProbabilities):
        raise ValueError("pairs must be a complete PairProbabilities result")
    spec = pairs.partition.complex_spec
    target, formatted = _target(spec, target_structure)
    return ComplexDefect(spec, formatted, tuple(
        1.0 - pairs.matrix[i][j] for i, j in enumerate(target.partners)))


@dataclass(frozen=True, slots=True)
class TubeTarget:
    name: str
    complex_spec: ComplexSpec
    structure: str
    concentration_molar: float

    def __post_init__(self):
        if not isinstance(self.name, str) or not self.name.strip():
            raise ValueError("A target name is required")
        if not isinstance(self.complex_spec, ComplexSpec):
            raise ValueError("Target complex_spec must be a ComplexSpec")
        _, formatted = _target(self.complex_spec, self.structure)
        if (isinstance(self.concentration_molar, bool)
                or not isinstance(self.concentration_molar, Real)
                or not math.isfinite(self.concentration_molar)
                or self.concentration_molar <= 0):
            raise ValueError("On-target concentration must be finite and positive")
        object.__setattr__(self, "name", self.name.strip())
        object.__setattr__(self, "structure", formatted)
        object.__setattr__(self, "concentration_molar", float(self.concentration_molar))


@dataclass(frozen=True, slots=True)
class TargetDefectContribution:
    target: TubeTarget
    actual_concentration_molar: float
    complex_defect: ComplexDefect | None
    structural_nucleotide_molar: float
    shortage_nucleotide_molar: float

    @property
    def total_nucleotide_molar(self) -> float:
        return self.structural_nucleotide_molar + self.shortage_nucleotide_molar


@dataclass(frozen=True, slots=True)
class TubeDefect:
    name: str
    contributions: tuple[TargetDefectContribution, ...]
    target_nucleotide_molar: float
    normalized_structural: float
    normalized_concentration: float

    @property
    def normalized(self) -> float:
        return self.normalized_structural + self.normalized_concentration

    @property
    def structural_nucleotide_molar(self) -> float:
        return math.fsum(c.structural_nucleotide_molar for c in self.contributions)

    @property
    def shortage_nucleotide_molar(self) -> float:
        return math.fsum(c.shortage_nucleotide_molar for c in self.contributions)


def tube_defect(name: str, tube: TubeResult, targets: Iterable[TubeTarget],
                pairs: Iterable[PairProbabilities]) -> TubeDefect:
    """Separate pairing errors from missing target complexes.

    Targets must account for all introduced strands, including residual free
    strands when the desired bound complexes have unequal strand totals.
    Matrices must come from the exact partitions used for the tube. No matrix
    is needed for a zero-concentration complex, whose structural contribution
    is zero (its concentration shortage still contributes in full).

    Coordinate rotations are intentionally not implicit: reorder a target or
    its matrix explicitly before calling this function.
    """
    if not isinstance(name, str) or not name.strip():
        raise ValueError("A tube name is required")
    if not isinstance(tube, TubeResult):
        raise ValueError("tube must be a validated TubeResult")
    targets = tuple(targets)
    if not targets or any(not isinstance(t, TubeTarget) for t in targets):
        raise ValueError("At least one TubeTarget is required")
    if len({t.name for t in targets}) != len(targets):
        raise ValueError("Target names must be unique in a tube")
    indices = {p.complex_spec.cyclic_key: i for i, p in enumerate(tube.partitions)}
    target_indices = []
    for target in targets:
        key = target.complex_spec.cyclic_key
        if key not in indices:
            raise ValueError("Target complex is missing from the tube ensemble")
        i = indices[key]
        if target.complex_spec.ordered_key != tube.partitions[i].complex_spec.ordered_key:
            raise ValueError("Target coordinates differ from the partition strand order")
        if i in target_indices:
            raise ValueError("A complex cannot have two target structures in one tube")
        target_indices.append(i)
    for s, total in enumerate(tube.total_concentrations_molar):
        requested = math.fsum(float(tube.stoichiometry[s, i]) * target.concentration_molar
                              for target, i in zip(targets, target_indices))
        if not math.isclose(requested, float(total), rel_tol=1e-10, abs_tol=0.0):
            raise ValueError("Target concentrations do not conserve introduced strand quantities")
    matrices = {}
    for matrix in pairs:
        if not isinstance(matrix, PairProbabilities):
            raise ValueError("Matrices must be complete PairProbabilities results")
        i = indices.get(matrix.partition.complex_spec.cyclic_key)
        if i is None or matrix.partition != tube.partitions[i]:
            raise ValueError("Pair probabilities do not match a tube partition exactly")
        if i in matrices:
            raise ValueError("Duplicate pair-probability result")
        matrices[i] = matrix
    # Normalize before multiplying tiny molar values, so a dilute tube does
    # not lose its dimensionless defect to underflow.
    scale = max(t.concentration_molar for t in targets)
    denominator = math.fsum(t.complex_spec.n_nucleotides * (t.concentration_molar / scale)
                            for t in targets)
    target_nt = scale * denominator
    if not math.isfinite(target_nt):
        raise ValueError("Target nucleotide concentration overflow")
    contributions, structural, concentration = [], [], []
    for target, i in zip(targets, target_indices):
        size, desired = target.complex_spec.n_nucleotides, target.concentration_molar
        actual = float(tube.concentrations_molar[i])
        if actual > 0 and i not in matrices:
            raise ValueError("A populated on-target complex needs its pair-probability matrix")
        defect = complex_defect(matrices[i], target.structure) if i in matrices else None
        fraction = 1.0 if actual >= desired else actual / desired
        weight = size * (desired / scale) / denominator
        d = defect.normalized if defect is not None else 0.0
        structural.append(weight * d * fraction)
        concentration.append(weight * (1.0 - fraction))
        contributions.append(TargetDefectContribution(
            target, actual, defect,
            (defect.expected_incorrect_nucleotides if defect is not None else 0.0)
            * min(actual, desired), size * max(desired - actual, 0.0)))
    return TubeDefect(name.strip(), tuple(contributions), target_nt,
                      math.fsum(structural), math.fsum(concentration))


@dataclass(frozen=True, slots=True)
class DesignDefect:
    tubes: tuple[TubeDefect, ...]

    def __post_init__(self):
        values = tuple(self.tubes)
        if not values or any(not isinstance(t, TubeDefect) for t in values):
            raise ValueError("A design defect requires at least one TubeDefect")
        if len({t.name for t in values}) != len(values):
            raise ValueError("Tube names must be unique")
        object.__setattr__(self, "tubes", values)

    @property
    def structural(self) -> float:
        return math.fsum(t.normalized_structural for t in self.tubes) / len(self.tubes)

    @property
    def concentration(self) -> float:
        return math.fsum(t.normalized_concentration for t in self.tubes) / len(self.tubes)

    @property
    def normalized(self) -> float:
        return self.structural + self.concentration
