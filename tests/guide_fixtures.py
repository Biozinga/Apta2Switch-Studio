"""Typed recording fixture for guide routing, not a scientific folding model."""
from types import SimpleNamespace
from aptaswitch_core.thermodynamics.setup import MODEL_ID
from aptaswitch_core.thermodynamics.types import ModelProvenance

class GuideBackend:
    provenance = ModelProvenance('fixture', '1', MODEL_ID, 'sha256:' + '1'*64,
                                 'fixture-connected', '1 M')
    condition_limits = {'temperature_min_c': 25., 'temperature_max_c': 55.,
                        'sodium_min_m': .12, 'sodium_max_m': 1.,
                        'magnesium_min_m': 0., 'magnesium_max_m': 0.}

    def __init__(self, records):
        self.records = iter(records)
        self.calls = []

    def minimum(self, spec, conditions):
        self.calls.append((spec, conditions))
        structure, energy = next(self.records)
        context = SimpleNamespace(complex_spec=spec, conditions=conditions,
                                  provenance=self.provenance, constraints=None)
        return SimpleNamespace(context=context, structure=structure, energy_kcal_mol=energy)
