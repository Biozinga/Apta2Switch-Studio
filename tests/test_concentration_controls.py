"""Concentration inputs preserve physical values through unit and page changes."""
import unittest

from streamlit.testing.v1 import AppTest


SCRIPT = """
import streamlit as st
from aptaswitch_core.localization import set_language
from aptaswitch_studio.concentration_ui import render_concentration_control, render_concentration_range
set_language('fr')
st.session_state.setdefault('physical_m', .175)
st.session_state.setdefault('range_m', (.15, .25))
st.radio('Page', ['controls', 'other'], key='page')
if st.session_state.page == 'controls':
    st.session_state.physical_m = render_concentration_control(
        'Na⁺', value_key='value', unit_key='unit', min_m=.12, max_m=1.,
        value_m=st.session_state.physical_m, default_unit='mM', help='Concentration du tube.',
    )
    st.session_state.range_m = render_concentration_range(
        'Na⁺', key='range', start_m=st.session_state.range_m[0], stop_m=st.session_state.range_m[1],
        min_m=.12, max_m=1.,
    )
"""


class ConcentrationControlsTests(unittest.TestCase):
    def test_all_units_preserve_concentration_and_adapt_numeric_limits(self):
        app = AppTest.from_string(SCRIPT).run()
        for unit, scale in (("M", 1.), ("µM", 1e6), ("nM", 1e9), ("mM", 1e3)):
            with self.subTest(unit=unit):
                app.selectbox(key="unit").set_value(unit).run()
                self.assertFalse(app.exception)
                self.assertEqual(app.session_state.physical_m, .175)
                self.assertEqual(app.number_input(key="value").value, .175 * scale)
                self.assertEqual(app.number_input(key="value").min, .12 * scale)
                self.assertEqual(app.number_input(key="value").max, 1. * scale)
        self.assertFalse(app.caption)
        self.assertEqual(app.number_input(key="value").proto.help, "Concentration du tube.")

    def test_shared_range_unit_preserves_both_bounds(self):
        app = AppTest.from_string(SCRIPT).run()
        for unit, scale in (("M", 1.), ("nM", 1e9), ("µM", 1e6), ("mM", 1e3)):
            with self.subTest(unit=unit):
                app.selectbox(key="range_unit").set_value(unit).run()
                self.assertFalse(app.exception)
                self.assertEqual(tuple(app.session_state.range_m), (.15, .25))
                self.assertEqual(app.number_input(key="range_min").value, .15 * scale)
                self.assertEqual(app.number_input(key="range_max").value, .25 * scale)

    def test_navigation_reconstructs_controls_in_last_selected_units(self):
        app = AppTest.from_string(SCRIPT).run()
        app.selectbox(key="unit").set_value("M").run()
        app.number_input(key="value").set_value(.185).run()
        app.selectbox(key="range_unit").set_value("µM").run()
        app.number_input(key="range_min").set_value(160000.).run()
        app.number_input(key="range_max").set_value(240000.).run()
        app.radio(key="page").set_value("other").run()
        self.assertFalse(app.number_input)
        app.radio(key="page").set_value("controls").run()
        self.assertFalse(app.exception)
        self.assertEqual(app.selectbox(key="unit").value, "M")
        self.assertEqual(app.number_input(key="value").value, .185)
        self.assertEqual(app.selectbox(key="range_unit").value, "µM")
        self.assertEqual(app.number_input(key="range_min").value, 160000.)
        self.assertEqual(app.number_input(key="range_max").value, 240000.)
        self.assertEqual(tuple(app.session_state.range_m), (.16, .24))


if __name__ == "__main__":
    unittest.main()
