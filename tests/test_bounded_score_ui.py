"""Score display, inverse-only weight and saved scientific values agree."""
from dataclasses import replace
import unittest

from aptaswitch_core.models import InverseExtensionSettings
from aptaswitch_core.scoring import selection_score_100
from aptaswitch_studio.design_session import open_design_into_session, validate_design_state
from aptaswitch_studio.design_library import create_design
from tests.test_scoring_ui import ScoringUITests
from tests.test_structure_ensemble import fixture
from tests.ui_helpers import FrenchAppTest as AppTest


def score_app():
    import streamlit as st
    from aptaswitch_studio.score_ui import render_candidate_score
    from aptaswitch_studio.results_table_ui import result_rows
    result = st.session_state.saved_result
    st.table(result_rows(result))
    render_candidate_score(result, result.candidates[0])


class TriggerWeightUITests(unittest.TestCase):
    setUp = ScoringUITests.setUp
    app = ScoringUITests.app
    launch = ScoringUITests.launch

    def select_inverse(self, app):
        app.checkbox(key='extension_enabled').check().run()
        self.assertEqual(app.radio(key='design_section').value, 'extension')
        app.radio(key='inverse_extension_enabled').set_value(True).run()
        app.radio(key='design_section').set_value('switch').run()
        self.assertFalse(app.exception)

    def test_weight_only_visible_in_inverse_mode_and_reaches_engine_request(self):
        app = self.app(trigger_input='GTGCTGGTATAA')
        self.assertNotIn('weight_trigger_linearity', [n.key for n in app.number_input])
        self.select_inverse(app)
        control = app.number_input(key='weight_trigger_linearity')
        self.assertEqual(control.value, 3.)
        self.assertIn('fraction moyenne de bases appariées', control.proto.help)
        control.set_value(2.5).run()
        request = self.launch(app)
        self.assertEqual(request.score_scale, 'bounded_v1')
        self.assertEqual(request.scoring.trigger_linearity, 2.5)
        self.assertIsNotNone(request.inverse_extension)

    def test_weight_zero_is_allowed_and_hidden_control_retains_user_preference(self):
        app = self.app(trigger_input='GTGCTGGTATAA', inverse_extension_enabled=True, weight_trigger_linearity=0.)
        self.assertTrue(app.checkbox(key='extension_enabled').value)
        app.checkbox(key='extension_enabled').uncheck().run()
        self.assertFalse(app.session_state['inverse_extension_enabled'])
        self.assertNotIn('weight_trigger_linearity', [n.key for n in app.number_input])
        app.radio(key='nav').set_value('Réglages').run()
        app.radio(key='nav').set_value('Conception').run()
        self.select_inverse(app)
        self.assertEqual(app.number_input(key='weight_trigger_linearity').value, 0.)
        self.assertEqual(self.launch(app).scoring.trigger_linearity, 0.)

    def test_disabling_extension_or_selecting_exhaustive_removes_inverse_from_request(self):
        for action in ('disable', 'exhaustive'):
            with self.subTest(action=action):
                app = self.app()
                self.select_inverse(app)
                app.number_input(key='weight_trigger_linearity').set_value(2.5).run()
                if action == 'disable':
                    app.checkbox(key='extension_enabled').uncheck().run()
                    self.assertEqual(app.radio(key='design_section').value, 'switch')
                    self.assertEqual(len(app.radio(key='design_section').options), 2)
                else:
                    app.radio(key='design_section').set_value('extension').run()
                    app.radio(key='inverse_extension_enabled').set_value(False).run()
                    self.assertNotIn('inverse_extension_target_length', [control.key for control in app.number_input])
                    app.radio(key='design_section').set_value('switch').run()
                self.assertFalse(app.exception)
                self.assertFalse(app.session_state['inverse_extension_enabled'])
                self.assertNotIn('weight_trigger_linearity', [control.key for control in app.number_input])
                self.assertEqual(app.session_state['weight_trigger_linearity'], 2.5)
                self.assertIsNone(self.launch(app).inverse_extension)


class BoundedResultUITests(unittest.TestCase):

    def test_bounded_score_is_not_transformed_twice_and_trigger_diagnostics_are_saved(self):
        result = fixture()
        result.request = replace(result.request, inverse_extension=InverseExtensionSettings())
        candidate = result.candidates[0]
        candidate.score_raw = 199.
        candidate.score = selection_score_100(199)
        candidate.trigger_linearity = {
            'normalized_defect': .2, 'unpaired_fraction': .8,
            'per_base_unpaired': [.8] * len(candidate.trigger_seq), 'structure_mfe': '.' * len(candidate.trigger_seq),
        }
        app = AppTest.from_function(score_app)
        app.session_state.saved_result = result
        app.run()
        self.assertFalse(app.exception)
        values = {metric.label: metric.value for metric in app.metric}
        self.assertAlmostEqual(float(values['Score /100']), candidate.score, places=3)
        self.assertEqual(values['Trigger seul · bases libres'], '80.00 %')
        self.assertEqual(values['Trigger seul · défaut de linéarité'], '0.2000')
        self.assertEqual(app.code[0].value, '.' * len(candidate.trigger_seq))


    def test_weight_restored_and_validated_in_portable_design_form(self):
        from pathlib import Path
        import tempfile
        with tempfile.TemporaryDirectory() as folder:
            design = create_design(Path(folder), 'Current design', state={'form': {'weight_trigger_linearity': 2.5}})
            form = {}
            open_design_into_session(design.path, form)
        self.assertEqual(form['weight_trigger_linearity'], 2.5)
        validate_design_state({'form': form})
        form['weight_trigger_linearity'] = -1.
        with self.assertRaises(ValueError):
            validate_design_state({'form': form})
