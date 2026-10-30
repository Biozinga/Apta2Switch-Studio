"""Exact joint accessibility of a window in an ordered connected complex.

The event is that every nucleotide in [start, end) is unpaired simultaneously.
Its weight comes from a partially constrained partition of the same complex,
never from independent marginal probabilities or an isolated RNA projection.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass, fields
import json
import math
from numbers import Real
from typing import Protocol

from .structures import PairingConstraints
from .types import ComplexSpec, Conditions, ModelProvenance, PartitionResult, Strand


class PartitionBackend(Protocol):
    provenance: ModelProvenance

    def partition(self, spec: ComplexSpec, conditions: Conditions,
                  constraints: PairingConstraints | None = None) -> PartitionResult: ...


def _finite_number(value, name):
    if isinstance(value, bool) or not isinstance(value, Real):
        raise ValueError(f"{name} must be a finite real number")
    try:
        result = float(value)
    except (OverflowError, ValueError) as error:
        raise ValueError(f"{name} must be representable as a finite real number") from error
    if not math.isfinite(result):
        raise ValueError(f"{name} must be a finite real number")
    return result


def _log_ratio(constrained, unrestricted):
    ratio = constrained - unrestricted
    if not math.isfinite(ratio):
        raise ValueError("Window log probability exceeds floating-point range")
    tolerance = 1e-10 + 1e-12 * max(abs(constrained), abs(unrestricted))
    if ratio > tolerance:
        raise ValueError("Window partition exceeds the unrestricted partition")
    return min(0., ratio)


@dataclass(frozen=True, slots=True)
class UnpairedWindowAccessibility:
    """Joint opening event with JSON-safe representation of impossible events.

    Coordinates refer to the supplied concatenated strand order. Opening free
    energy is -RT ln(probability), in kcal/mol. An empty constrained ensemble
    has probability zero and ``constrained_ensemble_empty=True``; its log
    probability (-infinity), opening energy (+infinity), and constrained log
    partition (-infinity) are represented by None for strict JSON export.
    A possible event whose probability underflows to zero retains its finite
    log probability and opening energy, with ``constrained_ensemble_empty=False``.
    """

    complex_spec: ComplexSpec
    conditions: Conditions
    provenance: ModelProvenance
    start: int
    end: int
    probability: float
    log_probability: float | None
    opening_free_energy_kcal_mol: float | None
    thermal_energy_kcal_mol: float
    unrestricted_log_partition: float
    constrained_log_partition: float | None
    constrained_ensemble_empty: bool

    def __post_init__(self):
        if (not isinstance(self.complex_spec, ComplexSpec)
                or not isinstance(self.conditions, Conditions)
                or not isinstance(self.provenance, ModelProvenance)):
            raise ValueError("Accessibility requires typed complex, conditions and model provenance")
        self.conditions.validate_for(self.complex_spec)
        if (type(self.start) is not int or type(self.end) is not int
                or not 0 <= self.start < self.end <= self.complex_spec.n_nucleotides):
            raise ValueError("Window must satisfy 0 <= start < end <= nucleotide count")
        if type(self.constrained_ensemble_empty) is not bool:
            raise ValueError("constrained_ensemble_empty must be bool")
        for name in ("probability", "thermal_energy_kcal_mol", "unrestricted_log_partition"):
            object.__setattr__(self, name, _finite_number(getattr(self, name), name))
        if not 0 <= self.probability <= 1 or self.thermal_energy_kcal_mol <= 0:
            raise ValueError("Accessibility requires probability in [0,1] and positive thermal energy")
        if self.constrained_ensemble_empty:
            if self.probability != 0 or any(value is not None for value in (
                    self.log_probability, self.opening_free_energy_kcal_mol, self.constrained_log_partition)):
                raise ValueError("Empty constrained ensembles require zero probability and null log/energy fields")
            return
        for name in ("log_probability", "opening_free_energy_kcal_mol", "constrained_log_partition"):
            object.__setattr__(self, name, _finite_number(getattr(self, name), name))
        if self.log_probability > 0 or self.opening_free_energy_kcal_mol < 0:
            raise ValueError("Log probability must be nonpositive and opening energy nonnegative")
        expected_log = _log_ratio(self.constrained_log_partition, self.unrestricted_log_partition)
        if not math.isclose(self.log_probability, expected_log, rel_tol=1e-12, abs_tol=1e-10):
            raise ValueError("Log probability differs from the constrained partition ratio")
        # No absolute tolerance here: a recorded zero is valid only when the
        # finite log probability actually underflows, or the ensemble is empty.
        if not math.isclose(self.probability, math.exp(self.log_probability), rel_tol=1e-10, abs_tol=0.):
            raise ValueError("Probability differs from its recorded log probability")
        energy = -self.thermal_energy_kcal_mol * self.log_probability
        if not math.isfinite(energy) or not math.isclose(
                self.opening_free_energy_kcal_mol, energy, rel_tol=1e-10, abs_tol=1e-8):
            raise ValueError("Opening energy must equal minus RT times log probability")

    def to_dict(self) -> dict:
        """Return the complete recorded context and JSON-safe numeric fields."""
        return json.loads(json.dumps(asdict(self), allow_nan=False))

    @classmethod
    def from_dict(cls, value: dict) -> UnpairedWindowAccessibility:
        """Validate persisted diagnostics without loading or invoking a backend."""
        if not isinstance(value, dict) or set(value) != {field.name for field in fields(cls)}:
            raise ValueError("Recorded accessibility fields are incomplete or unrecognized")
        try:
            spec, conditions, provenance = value["complex_spec"], value["conditions"], value["provenance"]
            if (not isinstance(spec, dict) or set(spec) != {"strands"}
                    or not isinstance(spec["strands"], (list, tuple))
                    or not isinstance(conditions, dict) or set(conditions) != {"celsius", "sodium", "magnesium"}
                    or not isinstance(provenance, dict)):
                raise ValueError("Recorded accessibility context is incomplete or malformed")
            return cls(**{**value,
                "complex_spec": ComplexSpec(tuple(Strand(**strand) for strand in spec["strands"])),
                "conditions": Conditions(**conditions), "provenance": ModelProvenance(**provenance)})
        except (TypeError, KeyError) as error:
            raise ValueError("Recorded accessibility context is incomplete or malformed") from error


def _require_context(partition, spec, conditions, provenance, constraints):
    if (not isinstance(partition, PartitionResult)
            or partition.complex_spec != spec
            or partition.conditions != conditions
            or partition.provenance != provenance
            or partition.constraints != constraints):
        raise ValueError("Partition must preserve the exact complex, strand order, chemistry, "
                         "labels, conditions, model and requested constraints")


def unpaired_window_accessibility(
    backend: PartitionBackend, spec: ComplexSpec, conditions: Conditions, *,
    start: int, end: int, partition: PartitionResult | None = None,
) -> UnpairedWindowAccessibility:
    """Compute P(all sites in [start,end) unpaired | this connected complex).

    The nonempty window uses zero-based concatenated coordinates and may cross
    strand boundaries. Sites outside it remain unresolved. No strand, chemistry,
    ordering, salt condition, or ensemble is substituted. The opening event is
    conditioned on the specified connected complex already being present; this
    is not a tube concentration, binding probability, or expression estimate.

    A supplied unrestricted partition is reused after checking its exact
    physical context. The constrained result must use the same model and RT.
    An empty unrestricted ensemble raises EmptyEnsembleError. Only tiny positive
    log-ratio endpoint roundoff is absorbed; larger inconsistencies raise.
    """
    if not isinstance(spec, ComplexSpec) or not isinstance(conditions, Conditions):
        raise ValueError("Explicit ComplexSpec and Conditions are required")
    conditions.validate_for(spec)
    if (type(start) is not int or type(end) is not int
            or not 0 <= start < end <= spec.n_nucleotides):
        raise ValueError("Window must satisfy 0 <= start < end <= nucleotide count")
    provenance = getattr(backend, "provenance", None)
    if not isinstance(provenance, ModelProvenance) or not callable(getattr(backend, "partition", None)):
        raise ValueError("Accessibility requires an explicit partition backend and model provenance")
    if partition is None:
        partition = backend.partition(spec, conditions)
    _require_context(partition, spec, conditions, provenance, None)
    partition.require_nonempty()
    constraints = PairingConstraints(tuple(i if start <= i < end else -1
                                           for i in range(spec.n_nucleotides)))
    constrained = backend.partition(spec, conditions, constraints)
    _require_context(constrained, spec, conditions, provenance, constraints)
    if backend.provenance != provenance:
        raise ValueError("The thermodynamic backend changed during the accessibility calculation")
    thermal = partition.thermal_energy_kcal_mol
    if constrained.thermal_energy_kcal_mol != thermal:
        raise ValueError("Constrained and unrestricted partitions use different thermal energy conventions")
    if constrained.is_empty:
        return UnpairedWindowAccessibility(spec, conditions, provenance, start, end,
            0., None, None, thermal, partition.log_partition, None, True)
    try:
        log_probability = _log_ratio(constrained.log_partition, partition.log_partition)
    except ValueError as error:
        raise ArithmeticError(str(error)) from error
    opening = 0. if log_probability == 0. else -thermal * log_probability
    if not math.isfinite(opening):
        raise ArithmeticError("Window opening energy exceeds floating-point range")
    return UnpairedWindowAccessibility(spec, conditions, provenance, start, end,
        math.exp(log_probability), log_probability, opening, thermal,
        partition.log_partition, constrained.log_partition, False)
