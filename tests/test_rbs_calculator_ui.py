"""Both result surfaces display the recorded, explicitly exploratory measurement."""
import unittest
from streamlit.testing.v1 import AppTest


class RBSCalculatorUITests(unittest.TestCase):
    def test_old_or_uncomputed_results_do_not_manufacture_metrics(self):
        app=AppTest.from_string('''
from types import SimpleNamespace
from aptaswitch_studio.rbs_calculator_ui import render_rbs_calculator
render_rbs_calculator(SimpleNamespace())
''').run()
        self.assertFalse(app.exception)
        self.assertFalse(app.metric)

    def test_both_languages_show_rates_and_context_without_recalculation(self):
        for lang,label,ref in [('fr','ON · ARN contraint','1,021 M'),('en','ON · constrained RNA','1.021 M')]:
            app=AppTest.from_string('''
from types import SimpleNamespace
from unittest.mock import patch
from aptaswitch_core.localization import set_language
from aptaswitch_studio.rbs_calculator_ui import render_rbs_calculator
set_language(LANGUAGE)
data={'off_rate':10.,'on_rate':1000.,'log10_on_off':2.,'merit':.4,
      'window_start':15,'window_end':85,'aug_start':50}
with patch('aptaswitch_core.rbs_calculator.predict_rbs_switch', side_effect=AssertionError('no compute on view')):
    render_rbs_calculator(SimpleNamespace(rbs_calculator=data))
'''.replace('LANGUAGE',repr(lang))).run()
            self.assertFalse(app.exception)
            self.assertEqual(app.metric[1].label,label)
            self.assertEqual(app.metric[2].value,'100×')
            self.assertEqual(app.metric[3].value,'40.00')
            self.assertTrue(any(ref in item.value for item in app.caption))
