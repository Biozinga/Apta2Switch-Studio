"""Native UI routing, unavailable runtimes and historical attribution."""
from dataclasses import asdict, replace
from types import SimpleNamespace
import unittest
from unittest.mock import patch

from aptaswitch_core.thermodynamics.setup import RuntimeValidation
from aptaswitch_studio import scientific_engines_ui as settings
from tests.native_ui_fixtures import patch_native_runtime
from tests.ui_helpers import FrenchAppTest as AppTest


class NativeThermodynamicsUITests(unittest.TestCase):
    def setUp(self):
        self.backend = patch_native_runtime(self)

    def test_identity_changes_when_binary_identity_changes_at_same_path(self):
        first = self.backend
        second = SimpleNamespace(provenance=replace(first.provenance, version='ui-fixture-v2'))
        with patch.object(settings, 'get_backend', side_effect=[first, second]):
            self.assertNotEqual(settings.native_backend_identity('/same/library'),
                                settings.native_backend_identity('/same/library'))

    def test_invalid_native_runtime_disables_identity_without_fallback(self):
        with patch.object(settings, 'get_backend', side_effect=RuntimeError('invalid digest')) as get_backend:
            self.assertIsNone(settings.native_backend_identity('/invalid/library'))
            self.assertIsNone(settings.native_backend_identity(None))
            self.assertEqual(get_backend.call_count, 1)

    def test_initial_structure_cache_uses_exact_identity(self):
        from aptaswitch_core.models import ThermoConditions
        from aptaswitch_studio.extension_ui import _analyze_initial_system
        _analyze_initial_system.clear()
        self.addCleanup(_analyze_initial_system.clear)
        with patch('aptaswitch_studio.extension_ui.get_backend', return_value=self.backend), \
             patch('aptaswitch_core.system_analysis.analyze_system', return_value={'mfe_structure': '..+..'}) as analyze:
            for identity in ('binary-A', 'binary-A', 'binary-B'):
                _analyze_initial_system('/same/library', 'AA', 'TT', ThermoConditions(), identity)
        self.assertEqual(analyze.call_count, 2)
        self.assertIs(analyze.call_args.kwargs['backend'], self.backend)
        self.assertIsNone(analyze.call_args.args[0])

    def test_saved_native_ensemble_cache_is_readable_without_runtime(self):
        from tests.test_structure_ensemble import fixture
        from aptaswitch_studio.ensemble_ui import _ensemble_input_options
        result = fixture()
        result.request = replace(result.request, engine='native-hybrid',
                                 thermo=replace(result.request.thermo, material=self.backend.provenance.model))
        result.thermodynamics = {**asdict(self.backend.provenance), 'engine': 'native-hybrid', 'conditions': {'celsius': 37}}
        with patch('aptaswitch_studio.ensemble_ui.get_backend', side_effect=AssertionError('No runtime while reopening')):
            options = _ensemble_input_options(result, result.candidates[0], 1000, allow_calculation=False)
        self.assertEqual(options[0]['thermodynamics'], asdict(self.backend.provenance))
        self.assertEqual(options[0]['thermo']['material'], self.backend.provenance.model)

    def test_finished_optimization_reports_unmet_target_honestly(self):
        script = """
from types import SimpleNamespace
from aptaswitch_studio.web_app import _render_optimization_status
_render_optimization_status(SimpleNamespace(optimization={'trials': [
    {'trial': 1, 'phase': 'global', 'evaluations': 37, 'best_loss': .2, 'best_global_loss': .3, 'stop_reason': 'budget_exhausted'},
    {'trial': 2, 'phase': 'local', 'evaluations': 12, 'best_loss': .001, 'best_global_loss': .002, 'stop_reason': 'target_reached'},
]}))
"""
        app = AppTest.from_string(script).run()
        self.assertFalse(app.exception)
        self.assertIn('1/2 essais ont atteint les cibles locale et globale', app.caption[0].value)
        self.assertEqual(list(app.table[0].value['Motif d’arrêt']),
                         ['Limite de recherche atteinte', 'Cibles locale et globale atteintes'])
        self.assertEqual(list(app.table[0].value['Défaut local (%)']), [20., .1])
        self.assertEqual(list(app.table[0].value['Défaut global (%)']), [30., .2])

    def test_optimization_status_distinguishes_exhausted_budgets(self):
        app = AppTest.from_string("""
from types import SimpleNamespace
from aptaswitch_studio.web_app import _render_optimization_status
_render_optimization_status(SimpleNamespace(optimization={'trials': [
    {'trial': 1, 'evaluations': 37, 'best_loss': .2, 'stop_reason': 'budget_exhausted', 'budget_kind': 'evaluations'},
    {'trial': 2, 'evaluations': 12, 'best_loss': .1, 'stop_reason': 'budget_exhausted', 'budget_kind': 'proposals'},
]}))
""").run()
        self.assertFalse(app.exception)
        self.assertEqual(list(app.table[0].value['Motif d’arrêt']),
                         ['Limite d’évaluations atteinte', 'Limite de propositions atteinte'])

    def test_no_admissible_sequence_label_is_limited_to_search_budget(self):
        for language, column, label in (
            ('fr', 'Motif d’arrêt', 'Aucune séquence admissible trouvée dans le budget'),
            ('en', 'Stopping reason', 'No admissible sequence found within the search budget'),
        ):
            with self.subTest(language=language):
                app = AppTest.from_string(f"""
from types import SimpleNamespace
from aptaswitch_core.localization import set_language
from aptaswitch_studio.web_app import _render_optimization_status
set_language({language!r})
_render_optimization_status(SimpleNamespace(optimization={{'trials': [
    {{'trial': 1, 'phase': 'initialization', 'evaluations': 0,
      'best_loss': None, 'best_global_loss': None, 'stop_reason': 'no_admissible_sequence'}},
]}}))
""").run()
                self.assertFalse(app.exception)
                self.assertEqual(list(app.table[0].value[column]), [label])

    def test_offline_unidentified_structure_has_no_invented_native_attribution(self):
        script = '''
import streamlit as st
from aptaswitch_studio.extension_ui import EXTENSION_DEFAULTS, render_initial_structure
for key, value in EXTENSION_DEFAULTS.items():
    st.session_state.setdefault(key, value)
render_initial_structure('ACGT', 'ACGT')
'''
        legacy = {'aptamer_sequence': 'ACGT', 'trigger_sequence': 'ACGT',
                  'mfe_structure': '((((+))))', 'mfe_energy_kcal_mol': -3.,
                  'paired_trigger_indices': [0, 1, 2, 3],
                  'model': {'material': 'external-archive', 'ensemble': 'archived',
                            'temperature_c': 37., 'sodium_m': .12, 'magnesium_m': 0.}}
        app = AppTest.from_string(script)
        app.session_state.system_analysis_signature = ('ACGT', 'ACGT', None, 'experimental-Sugimoto95-Vienna04-loopmean', 'dna', 37., .12, 0.)
        app.session_state.system_analysis_result = legacy
        with patch('aptaswitch_studio.extension_ui.configured_library_path', return_value=None), \
             patch('aptaswitch_studio.extension_ui._analyze_initial_system') as analyze:
            app.run()
        self.assertFalse(app.exception)
        analyze.assert_not_called()
        self.assertEqual(app.session_state.system_analysis_result, legacy)
        self.assertTrue(any('Moteur non renseigné' in caption.value for caption in app.caption))
        self.assertFalse(any(self.backend.provenance.backend in caption.value for caption in app.caption))
        conditions = next(caption.value for caption in app.caption if 'bases du trigger appariées' in caption.value)
        self.assertIn('external-archive', conditions)
        self.assertNotIn(self.backend.provenance.model, conditions)


if __name__ == '__main__':
    unittest.main()
