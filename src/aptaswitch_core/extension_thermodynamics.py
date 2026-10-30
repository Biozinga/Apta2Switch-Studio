"""Thermodynamic boundary shared by exhaustive and inverse extensions.

The search/ranking algorithms do not depend on an engine API. Selecting an
explicit backend never falls back to another engine. The initial system
uses DNA for both aptamer and trigger, independently of the switch's RNA chemistry.
"""
from __future__ import annotations

from contextlib import contextmanager
from dataclasses import asdict, dataclass
import math
from typing import Any

from .models import ThermoConditions

@dataclass(frozen=True)
class MfeResult:
    structure: str
    energy_kcal_mol: float


@dataclass
class ComplexAnalysisResult:
    mfe_structure: str
    mfe_energy_kcal_mol: float
    ensemble_free_energy_kcal_mol: float
    mfe_structure_probability: float
    target_structure_probability: float
    pair_matrix: Any


def extension_model_settings(thermo: ThermoConditions, backend=None) -> dict:
    """Resolve conditions and engine identity without starting a calculation."""
    if not isinstance(thermo, ThermoConditions):
        raise ValueError("thermo must be ThermoConditions")
    thermo.validate_model()
    if backend is None:
        from .thermodynamics.setup import get_backend
        backend = get_backend()
    from .thermodynamics.types import Conditions, ModelProvenance
    if not isinstance(getattr(backend, "provenance", None), ModelProvenance):
        raise ValueError("The extension backend requires explicit model provenance")
    thermo.validate_system_chemistry()
    conditions = Conditions(thermo.temperature_c, thermo.sodium_m, thermo.magnesium_m)
    limits = backend.condition_limits
    for name, value in (("temperature", conditions.celsius), ("sodium", conditions.sodium),
                        ("magnesium", conditions.magnesium)):
        unit = "c" if name == "temperature" else "m"
        if not limits[f"{name}_min_{unit}"] <= value <= limits[f"{name}_max_{unit}"]:
            raise ValueError(f"Extension {name} is outside the backend's declared model domain")
    material, ensemble = backend.provenance.model, backend.provenance.ensemble
    return dict(material=material, ensemble=ensemble,
                temperature_c=float(thermo.temperature_c), sodium_m=float(thermo.sodium_m),
                magnesium_m=float(thermo.magnesium_m),
                aptamer_material="dna", trigger_material="dna")


class _NativeExtensionModel:
    def __init__(self, backend, settings):
        from .thermodynamics.types import Conditions
        self.backend, self.settings = backend, settings
        self.provenance = backend.provenance
        self.conditions = Conditions(settings["temperature_c"], settings["sodium_m"], settings["magnesium_m"])
        self.metadata = {"thermodynamics": asdict(self.provenance)}

    def _spec(self, aptamer, trigger):
        from .thermodynamics.types import ComplexSpec, Strand
        return ComplexSpec((Strand(aptamer, "dna", "aptamer"), Strand(trigger, "dna", "trigger")))

    def _minimum(self, spec):
        from .thermodynamics.application_analysis import _check_context
        from .thermodynamics.structures import parse_structure, is_connected, format_structure
        record = self.backend.minimum(spec, self.conditions)
        _check_context(record.context, spec, self.conditions, self.provenance)
        if not math.isfinite(record.energy_kcal_mol):
            raise ValueError("The extension minimum energy is nonfinite")
        lengths = tuple(len(s.sequence) for s in spec.strands)
        state = parse_structure(record.structure, lengths)
        if not is_connected(state, lengths):
            raise ValueError("The extension minimum must belong to the connected complex")
        return MfeResult(format_structure(state, lengths), record.energy_kcal_mol), record.context

    def minimum(self, aptamer, trigger):
        return self._minimum(self._spec(aptamer, trigger))[0]

    def analyze(self, aptamer, trigger, *, target=None):
        import numpy as np
        from .thermodynamics.application_analysis import _check_context
        from .thermodynamics.types import PairProbabilities, PartitionResult
        spec = self._spec(aptamer, trigger)
        result = self.backend.analyze(spec, self.conditions, need_pairs=True)
        partition = result.partition
        if not isinstance(partition, PartitionResult):
            raise ValueError("Extension partition must identify its ensemble and RT")
        _check_context(partition, spec, self.conditions, self.provenance)
        partition.require_nonempty()
        if not isinstance(result.pairs, PairProbabilities) or result.pairs.partition != partition:
            raise ValueError("Extension pair probabilities must describe the same ensemble")
        mfe, context = self._minimum(spec)
        if not math.isclose(context.thermal_energy_kcal_mol, partition.thermal_energy_kcal_mol,
                            rel_tol=1e-12, abs_tol=0):
            raise ValueError("Extension minimum and partition use different RT conventions")
        def probability(structure):
            value = float(self.backend.structure_probability(spec, self.conditions, structure, partition=partition))
            if not math.isfinite(value) or not 0 <= value <= 1:
                raise ValueError("The extension structure probability is invalid")
            return value
        mfe_probability = probability(mfe.structure)
        target_probability = mfe_probability if target is None or target == mfe.structure else probability(target)
        return ComplexAnalysisResult(mfe.structure, mfe.energy_kcal_mol, partition.free_energy_kcal_mol,
            mfe_probability, target_probability, np.asarray(result.pairs.matrix))


@contextmanager
def open_extension_model(thermo, import_path=None, *, backend=None):
    """Open the selected native model; the historical import path is ignored."""
    if backend is None:
        from .thermodynamics.setup import get_backend
        backend = get_backend()
    settings = extension_model_settings(thermo, backend)
    yield _NativeExtensionModel(backend, settings)
