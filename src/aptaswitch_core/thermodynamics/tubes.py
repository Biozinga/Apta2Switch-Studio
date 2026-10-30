"""Convert compatible connected-complex partitions into tube concentrations.

Partitions already include association and cyclic-symmetry contributions. The
adapter neither inserts those corrections again nor merges distinct cyclic
orders. Its result describes exactly the supplied complex catalog; it cannot
establish that omitting larger complexes is physically justified.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
from numbers import Real
from typing import Iterable

import numpy as np
from numpy.typing import ArrayLike, NDArray

from .equilibrium import EquilibriumResult, solve_equilibrium
from .types import Conditions, ModelProvenance, PartitionResult, Strand


@dataclass(frozen=True, slots=True)
class StandardStateConvention:
    """An explicit, caller-audited convention for the supplied partitions.

    ``concentration_molar`` is the concentration corresponding to unit activity:
    use 1 for a 1 M standard, or an explicitly established water concentration
    for a mole-fraction standard. No value is inferred from ``provenance_label``.
    The label must match the backend's provenance exactly. This matching guards
    against accidental mixing; it is not evidence that a model implements the
    declared convention correctly. The caller must establish that separately.

    All three flags must explicitly be True. A partition missing association or
    symmetry terms, or including disconnected assemblies, needs a different
    upstream calculation rather than an implicit correction in this adapter.
    """

    provenance_label: str
    concentration_molar: float
    association_included: bool
    cyclic_symmetry_included: bool
    connected_ensembles: bool

    def __post_init__(self) -> None:
        if not isinstance(self.provenance_label, str) or not self.provenance_label.strip():
            raise ValueError("provenance_label must be a nonempty string.")
        object.__setattr__(self, "provenance_label", self.provenance_label.strip())
        if isinstance(self.concentration_molar, bool) or not isinstance(self.concentration_molar, Real):
            raise ValueError("concentration_molar must be an explicit positive concentration in mol/L.")
        concentration = float(self.concentration_molar)
        if not math.isfinite(concentration) or concentration <= 0:
            raise ValueError("concentration_molar must be finite and positive.")
        object.__setattr__(self, "concentration_molar", concentration)
        for flag in ("association_included", "cyclic_symmetry_included", "connected_ensembles"):
            if getattr(self, flag) is not True:
                raise ValueError(f"{flag} must explicitly be True for this tube adapter.")


@dataclass(frozen=True, slots=True)
class TubeResult:
    """A tube's species, original partition records, and validated equilibrium.

    Rows of ``stoichiometry`` and ``total_concentrations_molar`` follow species
    order. Columns, log K, and equilibrium concentrations follow partition
    order, including entries with empty connected ensembles. Arrays are
    read-only copies; the original immutable partition objects are preserved.
    """

    species: tuple[Strand, ...]
    partitions: tuple[PartitionResult, ...]
    total_concentrations_molar: NDArray[np.float64]
    stoichiometry: NDArray[np.int64]
    log_equilibrium_constants: NDArray[np.float64]
    standard_state: StandardStateConvention
    equilibrium: EquilibriumResult

    @property
    def concentrations_molar(self) -> NDArray[np.float64]:
        return self.equilibrium.concentrations_molar

    @property
    def log_concentrations_molar(self) -> NDArray[np.float64]:
        return self.equilibrium.log_concentrations_molar

    @property
    def conditions(self) -> Conditions:
        return self.partitions[0].conditions

    @property
    def provenance(self) -> ModelProvenance:
        return self.partitions[0].provenance

    @property
    def thermal_energy_kcal_mol(self) -> float:
        return self.partitions[0].thermal_energy_kcal_mol


def solve_tube(
    species: Iterable[Strand],
    total_concentrations_molar: ArrayLike,
    partitions: Iterable[PartitionResult],
    *,
    standard_state: StandardStateConvention,
    relative_tolerance: float = 1e-10,
    max_iterations: int = 1000,
) -> TubeResult:
    """Solve a fixed tube catalog from its compatible connected partitions.

    Every chemical strand type must occur once in ``species`` and have exactly
    one nonempty monomer partition. A partition may use other display labels
    for those same molecules: identity is always (material, sequence).
    Duplicate cyclic rotations are rejected; distinct cyclic orders with the
    same stoichiometry remain separate species. Input order is never changed.

    All records must have exactly identical model provenance, conditions, and
    backend RT. Do not combine constrained and unconstrained partition models
    under a shared provenance. The caller must supply the complete intended
    complex catalog and establish its standard-state convention explicitly.

    For a complex with n strands, the natural log formation constant in M is
    log Z_complex - sum(log Z_monomer) + (1-n)*log(c_standard).
    Using log Z directly retains the backend's exact RT convention without
    recomputing it using a different gas constant.
    """
    catalog = _catalog(species)
    records = _partitions(partitions)
    if not isinstance(standard_state, StandardStateConvention):
        raise ValueError("standard_state must be an explicit StandardStateConvention.")
    reference = records[0]
    if reference.provenance.standard_state != standard_state.provenance_label:
        raise ValueError("Partition provenance does not match the declared standard-state convention.")

    identities = {strand.chemical_identity: index for index, strand in enumerate(catalog)}
    matrix = np.zeros((len(catalog), len(records)), dtype=np.int64)
    monomers: dict[int, PartitionResult] = {}
    seen = set()
    for column, record in enumerate(records):
        if record.constraints is not None:
            raise ValueError("Constrained partition ensembles cannot define an unconstrained tube equilibrium.")
        if record.provenance != reference.provenance:
            raise ValueError("All tube partitions must have exactly the same model provenance and ensemble conventions.")
        if record.conditions != reference.conditions:
            raise ValueError("All tube partitions must have exactly the same temperature and salt conditions.")
        if record.thermal_energy_kcal_mol != reference.thermal_energy_kcal_mol:
            raise ValueError("All tube partitions must have exactly the same backend RT convention.")
        spec = record.complex_spec
        if spec.cyclic_key in seen:
            raise ValueError("Duplicate complex partitions, including cyclic rotations, cannot be counted twice.")
        seen.add(spec.cyclic_key)
        for strand in spec.strands:
            identity = strand.chemical_identity
            if identity not in identities:
                raise ValueError("A complex contains a chemical strand species absent from the tube catalog.")
            matrix[identities[identity], column] += 1
        if len(spec.strands) == 1:
            if record.is_empty:
                raise ValueError("Every monomer must have a nonempty partition ensemble.")
            monomers[identities[spec.strands[0].chemical_identity]] = record
    if len(monomers) != len(catalog):
        raise ValueError("A nonempty monomer partition is required for every tube strand species.")

    monomer_logs = [monomers[index].log_partition for index in range(len(catalog))]
    log_standard = math.log(standard_state.concentration_molar)
    constants = np.empty(len(records), dtype=float)
    for column, record in enumerate(records):
        if record.is_empty:
            constants[column] = -np.inf
        else:
            terms = [record.log_partition]
            terms.extend(-monomer_logs[identities[strand.chemical_identity]] for strand in record.complex_spec.strands)
            terms.append((1 - len(record.complex_spec.strands)) * log_standard)
            try:
                constants[column] = math.fsum(terms)
            except OverflowError as error:
                raise ValueError("Formation constants exceed the supported numerical range.") from error
            if not np.isfinite(constants[column]):
                raise ValueError("A nonempty complex must have a finite log formation constant.")

    # Snapshot without numeric coercion so the solver can reject original
    # booleans, strings, and complex values rather than silently changing them.
    original_totals = np.array(total_concentrations_molar, dtype=object, copy=True)
    equilibrium = solve_equilibrium(
        matrix, original_totals, log_equilibrium_constants=constants,
        relative_tolerance=relative_tolerance, max_iterations=max_iterations,
    )
    totals = np.array(original_totals, dtype=float)
    for array in (matrix, totals, constants):
        array.flags.writeable = False
    return TubeResult(catalog, records, totals, matrix, constants, standard_state, equilibrium)


def _catalog(species: Iterable[Strand]) -> tuple[Strand, ...]:
    try:
        catalog = tuple(species)
    except TypeError as error:
        raise ValueError("species must be an iterable of Strand objects.") from error
    if not catalog or any(not isinstance(strand, Strand) for strand in catalog):
        raise ValueError("species must contain at least one Strand.")
    if len({strand.chemical_identity for strand in catalog}) != len(catalog):
        raise ValueError("Tube species must have unique chemical identities; merge alias totals explicitly.")
    if len({strand.name for strand in catalog}) != len(catalog):
        raise ValueError("Tube species names cannot designate different chemical identities.")
    return catalog


def _partitions(partitions: Iterable[PartitionResult]) -> tuple[PartitionResult, ...]:
    try:
        records = tuple(partitions)
    except TypeError as error:
        raise ValueError("partitions must be an iterable of PartitionResult objects.") from error
    if not records or any(not isinstance(record, PartitionResult) for record in records):
        raise ValueError("partitions must contain at least one PartitionResult.")
    return records
