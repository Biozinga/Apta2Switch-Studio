"""Portable identified backend for UI tests; no scientific work or user settings."""
from types import SimpleNamespace
from unittest.mock import patch
from aptaswitch_core.thermodynamics.types import ModelProvenance
from aptaswitch_core.thermodynamics.setup import MODEL_ID
from aptaswitch_core.thermodynamics.tubes import StandardStateConvention


def patch_native_runtime(test):
    backend = SimpleNamespace(
        provenance=ModelProvenance('AptaSwitch-native-hybrid', 'ui-fixture-v1', MODEL_ID,
                                   'sha256:' + '1' * 64, 'ui-fixture-ensemble', '1 M'),
        condition_limits={'temperature_min_c': 25., 'temperature_max_c': 55.,
                          'sodium_min_m': .12, 'sodium_max_m': 1.,
                          'magnesium_min_m': 0., 'magnesium_max_m': 0.},
        standard_state=StandardStateConvention('1 M', 1.0, True, True, True),
    )
    for module in ('aptaswitch_core.thermodynamics.setup', 'aptaswitch_studio.scientific_engines_ui',
                   'aptaswitch_studio.extension_ui', 'aptaswitch_studio.extension_runtime', 'aptaswitch_studio.system_tube_runtime', 'aptaswitch_studio.ensemble_ui',
                   'aptaswitch_studio.saved_switch_ui', 'aptaswitch_studio.saved_switch_runtime'):
        test.enterContext(patch(module + '.get_backend', return_value=backend))
    test.enterContext(patch('aptaswitch_studio.ensemble_runtime.backend_for_request', return_value=backend))
    return backend
