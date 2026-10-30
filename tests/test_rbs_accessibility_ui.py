"""Live UI labels and numeric meaning of the recorded opening diagnostic."""
import unittest
from streamlit.testing.v1 import AppTest


class RBSAccessibilityUITests(unittest.TestCase):
    def test_candidate_without_opening_measurement_remains_viewable(self):
        app = AppTest.from_string('''
from aptaswitch_core.models import DesignCandidate
from aptaswitch_studio.accessibility_ui import render_rbs_accessibility
candidate = DesignCandidate.excluded_for_stop(trial=1, switch_seq='AUG',
    trigger_seq='CAT', first_stop=2, reporter_label='sfGFP')
del candidate.rbs_accessibility
render_rbs_accessibility(candidate)
''').run()
        self.assertFalse(app.exception)
        self.assertEqual(len(app.metric), 0)
        self.assertEqual(len(app.expander), 0)

    def app(self, language, *, impossible=False):
        return AppTest.from_string('''
from types import SimpleNamespace
from aptaswitch_core.localization import set_language
from aptaswitch_studio.accessibility_ui import render_rbs_accessibility
set_language(LANGUAGE)
measurement = {'probability': 0. if IMPOSSIBLE else .125,
    'opening_free_energy_kcal_mol': None if IMPOSSIBLE else 1.281,
    'start':60, 'end':83, 'constrained_ensemble_empty':IMPOSSIBLE}
render_rbs_accessibility(SimpleNamespace(trigger_seq='A'*30, rbs_accessibility=measurement))
'''.replace('LANGUAGE', repr(language)).replace('IMPOSSIBLE', repr(impossible))).run()

    def test_joint_measurement_is_clear_in_both_languages(self):
        for language, label, caption in (
            ('fr', 'Région entièrement libre (%)', 'sans malus supplémentaire'),
            ('en', 'Entire region unpaired (%)', 'no additional score penalty'),
        ):
            app = self.app(language)
            self.assertFalse(app.exception)
            self.assertEqual(app.metric[0].label, label)
            self.assertEqual(app.metric[0].value, '12.5')
            self.assertTrue(any(caption in item.value for item in app.caption))
            self.assertTrue(any('31' in item.value and '53' in item.value for item in app.caption))

    def test_impossible_opening_is_not_shown_as_missing_or_zero_energy(self):
        app = self.app('fr', impossible=True)
        self.assertFalse(app.exception)
        self.assertEqual(app.metric[0].value, '0')
        self.assertEqual(app.metric[1].value, '∞')
        self.assertIn('coût est infini', app.caption[-1].value)
