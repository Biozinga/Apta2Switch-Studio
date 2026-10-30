"""Exhaustive, engine-independent reference for small connected complexes.

This is a scientific reference, not a production folding algorithm. Complexity
and storage grow exponentially; the explicit nucleotide limit raises before
enumeration rather than silently truncating an ensemble. Pair admissibility
and minimum hairpin separation are required model inputs: no RNA alphabet or
energy table is assumed, and no DNA/RNA transcription is performed.

The energy callback returns one total energy in kcal/mol per structure. It
must already include every association, rotational-symmetry or other constant
required by its stated ensemble. This module adds no such correction. Minimum
energies consequently include those statistical corrections when supplied;
they are not implicitly a native engine's uncorrected structural MFE.
"""
from __future__ import annotations

from collections.abc import Callable, Iterator
from dataclasses import dataclass
import math
from numbers import Real
import random

from .structures import PairingConstraints, format_structure, is_connected
from .types import ComplexSpec, Conditions, ModelProvenance, PairProbabilities, PartitionResult


PairAllowed = Callable[[ComplexSpec, int, int], bool]
StructureEnergy = Callable[[ComplexSpec, Conditions, PairingConstraints], float]
ThermalEnergy = Callable[[Conditions], float]


class EnumerationLimitError(ValueError):
    """An exhaustive reference request exceeds its explicit length limit."""


def _integer(value: object, name: str, minimum: int) -> int:
    if type(value) is not int or value < minimum:
        raise ValueError(f"{name} must be an integer >= {minimum}")
    return value


def _real(value: object, name: str, *, allow_positive_infinity: bool = False) -> float:
    if isinstance(value, bool) or not isinstance(value, Real):
        raise ValueError(f"{name} must be a real number")
    try:
        result = float(value)
    except (OverflowError, ValueError) as error:
        raise ValueError(f"{name} must be representable as a float") from error
    if not math.isfinite(result) and not (allow_positive_infinity and result == math.inf):
        raise ValueError(f"{name} must be finite" + (" or +inf" if allow_positive_infinity else ""))
    return result


def enumerate_structures(
    complex_spec: ComplexSpec, *, pair_allowed: PairAllowed,
    min_hairpin_size: int, constraints: PairingConstraints | None = None,
    max_nucleotides: int = 18,
) -> Iterator[PairingConstraints]:
    """Yield each admissible connected, pseudoknot-free matching exactly once.

    Coordinates always refer to the original concatenated strand order. The
    deterministic ``pair_allowed(spec, i, j)`` callback is queried only for
    i < j. A pair within one strand requires j-i-1 >= min_hairpin_size;
    interstrand pairs have no distance minimum, including adjacent bases
    across a nick. Forbidden forced pairs give an empty iterator.
    """
    if not isinstance(complex_spec, ComplexSpec):
        raise ValueError("complex_spec must be a ComplexSpec")
    if not callable(pair_allowed):
        raise ValueError("pair_allowed must be callable")
    _integer(min_hairpin_size, "min_hairpin_size", 0)
    _integer(max_nucleotides, "max_nucleotides", 1)
    length = complex_spec.n_nucleotides
    if length > max_nucleotides:
        raise EnumerationLimitError(
            f"Exhaustive reference supports at most {max_nucleotides} nucleotides; received {length}")
    if constraints is None:
        constraints = PairingConstraints.unconstrained(length)
    if not isinstance(constraints, PairingConstraints) or len(constraints.partners) != length:
        raise ValueError("constraints must describe every nucleotide of this complex")
    lengths = tuple(len(strand.sequence) for strand in complex_spec.strands)
    strand_ids = tuple(k for k, n in enumerate(lengths) for _ in range(n))
    allowed: set[tuple[int, int]] = set()
    for i in range(length):
        for j in range(i + 1, length):
            if strand_ids[i] == strand_ids[j] and j - i - 1 < min_hairpin_size:
                continue
            accepts = pair_allowed(complex_spec, i, j)
            if type(accepts) is not bool:
                raise ValueError("pair_allowed must return bool")
            if accepts:
                allowed.add((i, j))
    required = constraints.partners

    def interval(first: int, stop: int) -> Iterator[tuple[int, ...]]:
        if first == stop:
            yield ()
            return
        forced = required[first]
        if forced in (-1, first):
            for suffix in interval(first + 1, stop):
                yield (first,) + suffix
        partners = range(first + 1, stop) if forced == -1 else (forced,)
        for j in partners:
            if (not first < j < stop or (first, j) not in allowed
                    or required[j] not in (-1, first)):
                continue
            for inside in interval(first + 1, j):
                for outside in interval(j + 1, stop):
                    yield (j,) + inside + (first,) + outside

    def connected_states() -> Iterator[PairingConstraints]:
        for partners in interval(0, length):
            state = PairingConstraints(partners)
            if is_connected(state, lengths):
                yield state

    return connected_states()


@dataclass(frozen=True, slots=True)
class ExhaustiveResult:
    """Immutable evaluated ensemble; obtain it from ``ExhaustiveBackend.evaluate``.

    Structures with energy +inf are absent. Each remaining structure is stored
    once, even when its relative Boltzmann weight underflows to zero. As with
    other floating-point thermodynamic calculations, extremely small weights
    need not be representable. Reusing this result never re-invokes the oracle.
    """

    partition: PartitionResult
    structures: tuple[PairingConstraints, ...]
    energies_kcal_mol: tuple[float, ...]
    _scaled_weights: tuple[float, ...]

    def __post_init__(self) -> None:
        # Copy containers so callers cannot mutate the evaluated ensemble.
        for name in ("structures", "energies_kcal_mol", "_scaled_weights"):
            object.__setattr__(self, name, tuple(getattr(self, name)))

    @property
    def n_structures(self) -> int:
        return len(self.structures)

    @property
    def probabilities(self) -> tuple[float, ...]:
        self.partition.require_nonempty()
        total = math.fsum(self._scaled_weights)
        return tuple(weight / total for weight in self._scaled_weights)

    def pairs(self) -> PairProbabilities:
        """Return exact enumerated marginals, with unpaired on the diagonal."""
        self.partition.require_nonempty()
        contributions: dict[tuple[int, int], list[float]] = {}
        for state, weight in zip(self.structures, self._scaled_weights):
            for i, j in enumerate(state.partners):
                contributions.setdefault((i, j), []).append(weight)
        n = self.partition.complex_spec.n_nucleotides
        total = math.fsum(self._scaled_weights)
        matrix = tuple(tuple(math.fsum(contributions.get((i, j), ())) / total
                             for j in range(n)) for i in range(n))
        return PairProbabilities(self.partition, matrix)

    def mfe(self) -> tuple[float, tuple[str, ...]]:
        """Return minimum supplied total energy and every exactly tied structure.

        Dot-bracket strings use '+' at the original strand boundaries. There
        is no energy rounding or tolerance-based merging of distinct levels.
        """
        self.partition.require_nonempty()
        minimum = min(self.energies_kcal_mol)
        lengths = tuple(len(s.sequence) for s in self.partition.complex_spec.strands)
        return minimum, tuple(format_structure(state, lengths)
                              for state, energy in zip(self.structures, self.energies_kcal_mol)
                              if energy == minimum)

    def sample(self, count: int, *, seed: int) -> list[str]:
        """Seeded independent Boltzmann draws over complete enumerated states."""
        _integer(count, "count", 0)
        if type(seed) is not int:
            raise ValueError("seed must be an integer")
        if count == 0:
            return []
        self.partition.require_nonempty()
        lengths = tuple(len(s.sequence) for s in self.partition.complex_spec.strands)
        draws = random.Random(seed).choices(self.structures, weights=self._scaled_weights, k=count)
        return [format_structure(state, lengths) for state in draws]


class ExhaustiveBackend:
    """Reference backend with explicit, deterministic model callbacks.

    ``thermal_energy(conditions)`` supplies the model's exact RT in kcal/mol.
    ``energy(spec, conditions, complete_matching)`` must validate its own
    temperature/salt/model domain and include all ensemble corrections. The
    original typed inputs, order and chemistry reach both model callbacks.
    Model provenance is preserved verbatim. No external engine is imported.

    ``evaluate`` invokes energy once for each connected admissible structure.
    Keep its result to obtain multiple observables without re-evaluation;
    separate convenience-method calls evaluate independent ensembles.
    """

    def __init__(self, *, pair_allowed: PairAllowed, energy: StructureEnergy,
                 thermal_energy: ThermalEnergy, provenance: ModelProvenance,
                 min_hairpin_size: int, max_nucleotides: int = 18):
        if not all(callable(callback) for callback in (pair_allowed, energy, thermal_energy)):
            raise ValueError("pair_allowed, energy and thermal_energy must be callable")
        if not isinstance(provenance, ModelProvenance):
            raise ValueError("provenance must be ModelProvenance")
        _integer(min_hairpin_size, "min_hairpin_size", 0)
        _integer(max_nucleotides, "max_nucleotides", 1)
        self.pair_allowed, self.energy, self.thermal_energy = pair_allowed, energy, thermal_energy
        self.provenance = provenance
        self.min_hairpin_size, self.max_nucleotides = min_hairpin_size, max_nucleotides

    def evaluate(self, complex_spec: ComplexSpec, conditions: Conditions,
                 constraints: PairingConstraints | None = None) -> ExhaustiveResult:
        if not isinstance(conditions, Conditions):
            raise ValueError("conditions must be Conditions")
        conditions.validate_for(complex_spec)
        states = enumerate_structures(complex_spec, pair_allowed=self.pair_allowed,
            min_hairpin_size=self.min_hairpin_size, constraints=constraints,
            max_nucleotides=self.max_nucleotides)
        thermal = _real(self.thermal_energy(conditions), "thermal energy")
        if thermal <= 0:
            raise ValueError("thermal energy must be strictly positive")
        structures, energies = [], []
        for state in states:
            energy = _real(self.energy(complex_spec, conditions, state), "structure energy",
                           allow_positive_infinity=True)
            if energy != math.inf:
                structures.append(state)
                energies.append(energy)
        if not energies:
            weights, log_z, free_energy = (), -math.inf, math.inf
        else:
            minimum = min(energies)
            weights = tuple(math.exp((minimum - energy) / thermal) for energy in energies)
            log_z = -minimum / thermal + math.log(math.fsum(weights))
            free_energy = -thermal * log_z
            if not math.isfinite(log_z) or not math.isfinite(free_energy):
                raise ArithmeticError("Enumerated partition exceeds floating-point log/energy range")
        partition = PartitionResult(complex_spec, conditions, self.provenance,
                                    log_z, free_energy, thermal, constraints=constraints)
        return ExhaustiveResult(partition, tuple(structures), tuple(energies), weights)

    def partition(self, complex_spec: ComplexSpec, conditions: Conditions,
                  constraints: PairingConstraints | None = None) -> PartitionResult:
        return self.evaluate(complex_spec, conditions, constraints).partition

    def pairs(self, complex_spec: ComplexSpec, conditions: Conditions,
              constraints: PairingConstraints | None = None) -> PairProbabilities:
        return self.evaluate(complex_spec, conditions, constraints).pairs()

    def mfe(self, complex_spec: ComplexSpec, conditions: Conditions,
            constraints: PairingConstraints | None = None) -> tuple[float, tuple[str, ...]]:
        return self.evaluate(complex_spec, conditions, constraints).mfe()

    def sample(self, complex_spec: ComplexSpec, conditions: Conditions, count: int, *,
               seed: int, constraints: PairingConstraints | None = None) -> list[str]:
        return self.evaluate(complex_spec, conditions, constraints).sample(count, seed=seed)
