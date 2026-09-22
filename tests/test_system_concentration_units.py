"""System salt units preserve physical conditions through analysis and extension."""

from __future__ import annotations

import tempfile
import unittest
from contextlib import ExitStack
from pathlib import Path
from unittest.mock import patch

from tests.ui_helpers import FrenchAppTest as AppTest

from aptaswitch_studio.extension_ui import _analyze_initial_system
from tests import test_system_sections as system_sections
from tests.test_system_sections import APP_PATH, APTAMER, TRIGGER, button


class SystemConcentrationUnitsTests(unittest.TestCase):
    def setUp(self):
        self.context = ExitStack()
        self.addCleanup(self.context.close)
        self.output = Path(self.context.enter_context(tempfile.TemporaryDirectory()))
        self.context.enter_context(patch("aptaswitch_core.nupack_setup.load_settings", return_value={}))
        self.context.enter_context(patch("aptaswitch_core.nupack_setup.detect_nupack", return_value=None))
        for save_function in ("save_settings", "save_nupack_settings"):
            self.context.enter_context(patch(
                f"aptaswitch_core.nupack_setup.{save_function}",
                side_effect=AssertionError("Tests must not save user settings"),
            ))
        for module in ("aptaswitch_core.nupack_setup", "aptaswitch_studio.extension_ui"):
            self.context.enter_context(patch(
                f"{module}.configured_nupack_import_path", return_value="/fixture/nupack",
            ))
        self.analyze = self.context.enter_context(patch(
            "aptaswitch_core.system_analysis.analyze_system",
            side_effect=system_sections.SystemSectionsTests.analysis_result,
        ))
        _analyze_initial_system.clear()
        self.addCleanup(_analyze_initial_system.clear)

    def app(self, **session_values) -> AppTest:
        app = AppTest.from_file(str(APP_PATH), default_timeout=20)
        app.session_state["dependency_detection_done"] = True
        app.session_state["output_dir"] = str(self.output)
        for key, value in session_values.items():
            app.session_state[key] = value
        app.run()
        self.assertFalse(app.exception)
        return app

    def fill_sequences(self, app: AppTest) -> None:
        app.text_area(key="aptamer_input").set_value(APTAMER)
        app.text_area(key="trigger_input").set_value(TRIGGER).run()
        self.assertFalse(app.exception)

    def test_system_salts_offer_same_units_as_switch_with_millimolar_defaults(self):
        app = self.app()
        for salt in ("sodium", "magnesium"):
            system_unit = app.selectbox(key=f"extension_{salt}_unit")
            design_unit = app.selectbox(key=f"{salt}_concentration_unit")
            self.assertEqual(system_unit.options, design_unit.options)
            self.assertEqual(set(system_unit.options), {"nM", "µM", "mM", "M"})
            self.assertEqual(system_unit.value, "mM")
        self.assertEqual(app.number_input(key="extension_sodium").value, 120.0)
        self.assertEqual(app.number_input(key="extension_magnesium").value, 0.0)
        self.analyze.assert_not_called()

    def test_unit_changes_preserve_salts_and_do_not_repeat_structure_calculation(self):
        app = self.app(extension_magnesium=7.0)
        self.fill_sequences(app)
        initial = app.session_state["system_analysis_result"]
        for unit, sodium, magnesium in (
            ("M", 0.12, 0.007),
            ("µM", 120000.0, 7000.0),
            ("nM", 120000000.0, 7000000.0),
            ("mM", 120.0, 7.0),
        ):
            with self.subTest(unit=unit):
                app.selectbox(key="extension_sodium_unit").set_value(unit)
                app.selectbox(key="extension_magnesium_unit").set_value(unit).run()
                self.assertFalse(app.exception)
                self.assertAlmostEqual(app.number_input(key="extension_sodium").value, sodium)
                self.assertAlmostEqual(app.number_input(key="extension_magnesium").value, magnesium)
                self.assertEqual(app.session_state["system_analysis_result"], initial)
                self.analyze.assert_called_once()
        thermo = self.analyze.call_args.args[3]
        self.assertEqual((thermo.sodium_m, thermo.magnesium_m), (0.12, 0.007))

    def test_molar_and_micromolar_entries_feed_both_structure_and_extension(self):
        app = self.app()
        app.selectbox(key="extension_sodium_unit").set_value("M")
        app.selectbox(key="extension_magnesium_unit").set_value("µM").run()
        app.number_input(key="extension_sodium").set_value(0.08)
        app.number_input(key="extension_magnesium").set_value(2500.0).run()
        self.fill_sequences(app)
        self.analyze.assert_called_once()
        thermo = self.analyze.call_args.args[3]
        self.assertEqual((thermo.sodium_m, thermo.magnesium_m), (0.08, 0.0025))
        with patch("aptaswitch_studio.extension_runtime.ExtensionRunJob.start") as start:
            button(app, "Étendre le trigger").click().run()
        self.assertFalse(app.exception)
        start.assert_called_once()
        request = app.session_state["extension_job"].request
        self.assertEqual(request.thermo, thermo)
        self.assertEqual((request.thermo.sodium_m, request.thermo.magnesium_m), (0.08, 0.0025))
        self.analyze.assert_called_once()

    def test_custom_units_and_values_survive_results_navigation(self):
        app = self.app()
        app.selectbox(key="extension_sodium_unit").set_value("M")
        app.selectbox(key="extension_magnesium_unit").set_value("µM").run()
        app.number_input(key="extension_sodium").set_value(0.075)
        app.number_input(key="extension_magnesium").set_value(3250.0).run()
        self.fill_sequences(app)
        app.radio(key="nav").set_value("Résultats").run()
        app.radio(key="nav").set_value("Conception").run()
        self.assertFalse(app.exception)
        self.assertEqual(app.selectbox(key="extension_sodium_unit").value, "M")
        self.assertEqual(app.selectbox(key="extension_magnesium_unit").value, "µM")
        self.assertEqual(app.number_input(key="extension_sodium").value, 0.075)
        self.assertEqual(app.number_input(key="extension_magnesium").value, 3250.0)
        self.analyze.assert_called_once()
        app.selectbox(key="extension_sodium_unit").set_value("mM")
        app.selectbox(key="extension_magnesium_unit").set_value("mM").run()
        self.assertFalse(app.exception)
        self.assertEqual(app.number_input(key="extension_sodium").value, 75.0)
        self.assertEqual(app.number_input(key="extension_magnesium").value, 3.25)
        self.analyze.assert_called_once()

    def test_previous_sessions_without_units_keep_millimolar_values(self):
        app = self.app(extension_sodium=42.0, extension_magnesium=7.0)
        self.fill_sequences(app)
        self.assertEqual(app.selectbox(key="extension_sodium_unit").value, "mM")
        self.assertEqual(app.selectbox(key="extension_magnesium_unit").value, "mM")
        self.assertEqual(app.number_input(key="extension_sodium").value, 42.0)
        self.assertEqual(app.number_input(key="extension_magnesium").value, 7.0)
        thermo = self.analyze.call_args.args[3]
        self.assertEqual((thermo.sodium_m, thermo.magnesium_m), (0.042, 0.007))


if __name__ == "__main__":
    unittest.main()
