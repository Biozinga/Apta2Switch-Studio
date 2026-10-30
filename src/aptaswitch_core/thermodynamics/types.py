"""Immutable contracts for connected nucleic-acid complex ensembles.

All positions refer to the supplied strand order, from zero at its first base.
Rotations identify ordered complexes for caching, but never reorder result
coordinates. Concentrations and chemical identities must not depend on labels.
This module does not select an energy model or import a folding engine.
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import accumulate, product
import math
from numbers import Integral, Real
from typing import Iterable, Literal, TypeAlias

from .structures import PairingConstraints


Material: TypeAlias = Literal["dna", "rna"]
StrandIdentity: TypeAlias = tuple[Material, str]
ComplexIdentity: TypeAlias = tuple[StrandIdentity, ...]

_PROBABILITY_TOLERANCE = 1e-9


def _text(value: object, field: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{field} must be a nonempty string")
    return value.strip()


def _number(value: object, field: str, *, finite: bool = True) -> float:
    if isinstance(value, bool) or not isinstance(value, Real):
        raise ValueError(f"{field} must be a real number")
    try:
        result = float(value)
    except (OverflowError, ValueError) as error:
        raise ValueError(f"{field} must be representable as a float") from error
    if math.isnan(result) or (finite and not math.isfinite(result)):
        raise ValueError(f"{field} must be {'finite' if finite else 'non-NaN'}")
    return result


@dataclass(frozen=True, slots=True)
class Strand:
    """One homogeneous strand; names are labels, not chemical identities.

    Case and surrounding whitespace are normalized. Internal whitespace,
    ambiguous bases, separators, and the other chemistry's T/U are rejected.
    No transcription or chemistry inference is performed.
    """

    sequence: str
    material: Material
    name: str

    def __post_init__(self) -> None:
        if self.material not in ("dna", "rna"):
            raise ValueError("material must be explicitly 'dna' or 'rna'")
        sequence = _text(self.sequence, "sequence").upper()
        alphabet = "ACGT" if self.material == "dna" else "ACGU"
        if any(base not in alphabet for base in sequence):
            raise ValueError(f"{self.material} sequence must contain only {alphabet}")
        object.__setattr__(self, "sequence", sequence)
        object.__setattr__(self, "name", _text(self.name, "name"))

    @property
    def chemical_identity(self) -> StrandIdentity:
        return self.material, self.sequence


def _strands(values: Iterable[Strand], field: str) -> tuple[Strand, ...]:
    try:
        strands = tuple(values)
    except TypeError as error:
        raise ValueError(f"{field} must be an iterable of Strand objects") from error
    if not strands or any(not isinstance(strand, Strand) for strand in strands):
        raise ValueError(f"{field} must contain at least one Strand")
    names: dict[str, StrandIdentity] = {}
    for strand in strands:
        previous = names.setdefault(strand.name, strand.chemical_identity)
        if previous != strand.chemical_identity:
            raise ValueError(f"strand name {strand.name!r} denotes different chemical species")
    return strands


@dataclass(frozen=True, slots=True)
class ComplexSpec:
    """An ordered connected-ensemble request, whether or not it can bind.

    ``strands`` preserves caller order. Dataclass equality also preserves that
    order and its labels; use ``cyclic_key`` only for physical identification.
    Cached matrices/structures require an explicit coordinate rotation before
    reuse for a differently ordered request. Mirror orders are not merged.
    """

    strands: tuple[Strand, ...]

    def __post_init__(self) -> None:
        object.__setattr__(self, "strands", _strands(self.strands, "strands"))

    @property
    def ordered_key(self) -> ComplexIdentity:
        return tuple(strand.chemical_identity for strand in self.strands)

    @property
    def cyclic_key(self) -> ComplexIdentity:
        key = self.ordered_key
        return min(key[offset:] + key[:offset] for offset in range(len(key)))

    @property
    def symmetry_factor(self) -> int:
        """Number of rotations preserving chemical identities, including zero."""
        key = self.ordered_key
        return sum(key[offset:] + key[:offset] == key for offset in range(len(key)))

    @property
    def sequence(self) -> str:
        """Concatenated bases, without separators or any T/U conversion."""
        return "".join(strand.sequence for strand in self.strands)

    @property
    def cutpoints(self) -> tuple[int, ...]:
        """Zero-based indices at which subsequent strands start."""
        return tuple(accumulate(len(strand.sequence) for strand in self.strands[:-1]))

    @property
    def n_nucleotides(self) -> int:
        return sum(len(strand.sequence) for strand in self.strands)

    @property
    def is_mixed(self) -> bool:
        return len({strand.material for strand in self.strands}) > 1


def enumerate_complexes(
    species: Iterable[Strand], max_size: int = 5,
) -> tuple[ComplexSpec, ...]:
    """Enumerate each chemical necklace once, including homomers and monomers.

    Different cyclic orders with the same stoichiometry remain distinct.
    Catalog aliases/duplicates are rejected: callers must explicitly merge
    concentrations of chemically identical species before enumeration.
    Enumeration is exhaustive (exponential in ``max_size``), intended for the
    application's small catalogs and complexes of at most five strands.
    """
    if isinstance(max_size, bool) or not isinstance(max_size, Integral) or max_size < 1:
        raise ValueError("max_size must be a positive integer")
    catalog = _strands(species, "species")
    if len({strand.chemical_identity for strand in catalog}) != len(catalog):
        raise ValueError("species must have distinct chemical identities; merge aliases explicitly")
    catalog = tuple(sorted(catalog, key=lambda strand: strand.chemical_identity))
    complexes: list[ComplexSpec] = []
    for size in range(1, int(max_size) + 1):
        for order in product(range(len(catalog)), repeat=size):
            if any(order[offset:] + order[:offset] < order for offset in range(1, size)):
                continue
            complexes.append(ComplexSpec(tuple(catalog[index] for index in order)))
    return tuple(complexes)


@dataclass(frozen=True, slots=True)
class Conditions:
    """Temperature in degrees Celsius and salt concentrations in mol/L.

    Only physical numeric constraints are enforced here. Each energy model
    must separately enforce its experimentally supported temperature/salt
    domain. The current mixed-material contract excludes magnesium.
    """

    celsius: float = 37.0
    sodium: float = 0.12
    magnesium: float = 0.0

    def __post_init__(self) -> None:
        for field in ("celsius", "sodium", "magnesium"):
            object.__setattr__(self, field, _number(getattr(self, field), field))
        if self.celsius <= -273.15:
            raise ValueError("temperature must be strictly above absolute zero")
        if self.sodium < 0 or self.magnesium < 0:
            raise ValueError("salt concentrations must be nonnegative")

    @property
    def kelvin(self) -> float:
        return self.celsius + 273.15

    def validate_for(self, complex_spec: ComplexSpec) -> None:
        if not isinstance(complex_spec, ComplexSpec):
            raise ValueError("complex_spec must be a ComplexSpec")
        if complex_spec.is_mixed and self.magnesium != 0:
            raise ValueError("mixed DNA/RNA ensembles currently require magnesium=0")


@dataclass(frozen=True, slots=True)
class ModelProvenance:
    """Exact backend/model/parameter identity and ensemble conventions.

    ``parameter_digest`` must identify the parameter content, including its
    hash algorithm (for example ``sha256:<hex>``), rather than a file path.
    ``ensemble`` records structural and stacking conventions; ``standard_state``
    records concentration units and association/symmetry conventions. Their
    vocabulary belongs to the energy backend, not to these shared types.
    """

    backend: str
    version: str
    model: str
    parameter_digest: str
    ensemble: str
    standard_state: str

    def __post_init__(self) -> None:
        for field in (
            "backend", "version", "model", "parameter_digest", "ensemble", "standard_state",
        ):
            object.__setattr__(self, field, _text(getattr(self, field), field))
        algorithm, separator, digest = self.parameter_digest.partition(":")
        if not separator or not algorithm or not digest or any(
            character.isspace() for character in self.parameter_digest
        ):
            raise ValueError("parameter_digest must include its hash algorithm and digest")


class EmptyEnsembleError(ValueError):
    """No connected structure exists for a normalized observable or sample."""


@dataclass(frozen=True, slots=True)
class PartitionResult:
    """Connected-ensemble partition and free energy for one supplied order.

    ``log_partition`` is natural log Z. ``thermal_energy_kcal_mol`` is the
    backend's exact RT convention, not an independently assumed gas constant.
    The empty ensemble has exactly log Z = -inf and G = +inf. Other nonfinite
    values are invalid. Finite G must equal -RT log Z to numerical tolerance
    (1e-10 relative, 1e-8 kcal/mol absolute).

    ``constraints`` identifies a conditioned subset of the model's ensemble.
    None means the full ensemble; an all-unresolved matching is normalized to
    None. Conditional partitions must not be used as unrestricted formation
    energies when solving tube concentrations.
    """

    complex_spec: ComplexSpec
    conditions: Conditions
    provenance: ModelProvenance
    log_partition: float
    free_energy_kcal_mol: float
    thermal_energy_kcal_mol: float
    constraints: PairingConstraints | None = None

    def __post_init__(self) -> None:
        if not isinstance(self.conditions, Conditions):
            raise ValueError("conditions must be Conditions")
        self.conditions.validate_for(self.complex_spec)
        if not isinstance(self.provenance, ModelProvenance):
            raise ValueError("provenance must be ModelProvenance")
        if self.constraints is not None:
            if (not isinstance(self.constraints, PairingConstraints)
                    or len(self.constraints.partners) != self.complex_spec.n_nucleotides):
                raise ValueError("constraints must describe every nucleotide of this complex")
            if all(partner == -1 for partner in self.constraints.partners):
                object.__setattr__(self, "constraints", None)
        log_z = _number(self.log_partition, "log_partition", finite=False)
        energy = _number(self.free_energy_kcal_mol, "free_energy_kcal_mol", finite=False)
        thermal = _number(self.thermal_energy_kcal_mol, "thermal_energy_kcal_mol")
        if thermal <= 0:
            raise ValueError("thermal_energy_kcal_mol must be strictly positive")
        if log_z == -math.inf:
            if energy != math.inf:
                raise ValueError("an empty ensemble requires free_energy_kcal_mol=+inf")
        elif not math.isfinite(log_z) or not math.isfinite(energy):
            raise ValueError("nonempty partition results must be finite")
        elif not math.isclose(energy, -thermal * log_z, rel_tol=1e-10, abs_tol=1e-8):
            raise ValueError("free energy must equal -thermal_energy_kcal_mol * log_partition")
        object.__setattr__(self, "log_partition", log_z)
        object.__setattr__(self, "free_energy_kcal_mol", energy)
        object.__setattr__(self, "thermal_energy_kcal_mol", thermal)

    @property
    def is_empty(self) -> bool:
        return self.log_partition == -math.inf

    @property
    def partition_function(self) -> float:
        """Return Z; may underflow or raise OverflowError. Prefer log_partition."""
        return math.exp(self.log_partition)

    def require_nonempty(self) -> None:
        """Guard MFE, normalized probabilities, and Boltzmann sampling calls."""
        if self.is_empty:
            raise EmptyEnsembleError("the connected ensemble is empty")


@dataclass(frozen=True, slots=True)
class PairProbabilities:
    """Marginals of precisely ``partition``; diagonal means unpaired.

    Matrix coordinates follow ``partition.complex_spec.strands`` unchanged.
    An immutable tuple copy is made. Values must be finite and within [0, 1],
    symmetric, with each row summing to one (absolute tolerance 1e-9). No
    clamping, renormalization, averaging, or symmetrization is performed.
    """

    partition: PartitionResult
    matrix: tuple[tuple[float, ...], ...]

    def __post_init__(self) -> None:
        if not isinstance(self.partition, PartitionResult):
            raise ValueError("partition must be a PartitionResult")
        self.partition.require_nonempty()
        try:
            matrix = tuple(tuple(_number(value, "pair probability") for value in row)
                           for row in self.matrix)
        except TypeError as error:
            raise ValueError("matrix must be a square iterable of real numbers") from error
        length = self.partition.complex_spec.n_nucleotides
        if len(matrix) != length or any(len(row) != length for row in matrix):
            raise ValueError("matrix shape must match the complex nucleotide count")
        for i, row in enumerate(matrix):
            if any(value < 0 or value > 1 for value in row):
                raise ValueError("pair probabilities must lie in [0, 1]")
            if not math.isclose(math.fsum(row), 1.0, rel_tol=0.0,
                                abs_tol=_PROBABILITY_TOLERANCE):
                raise ValueError("pair probability rows must sum to one, including unpaired diagonal")
            if any(not math.isclose(row[j], matrix[j][i], rel_tol=0.0,
                                    abs_tol=_PROBABILITY_TOLERANCE) for j in range(i)):
                raise ValueError("pair probability matrix must be symmetric")
        object.__setattr__(self, "matrix", matrix)
