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
        from tests.native_ui_fixtures import patch_native_runtime
        patch_native_runtime(self)
        self.context = ExitStack()
        self.addCleanup(self.context.close)
        self.output = Path(self.context.enter_context(tempfile.TemporaryDirectory()))
        self.context.enter_context(patch("aptaswitch_core.settings.load_settings", return_value={}))
        for save_function in ("save_settings",):
            self.context.enter_context(patch(
                f"aptaswitch_core.settings.{save_function}",
                side_effect=AssertionError("Tests must not save user settings"),
            ))
        for module in ("aptaswitch_core.thermodynamics.setup", "aptaswitch_studio.extension_ui"):
            self.context.enter_context(patch(
                f"{module}.configured_library_path", return_value="/fixture/native",
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
        system_options = {}
        for salt in ("sodium", "magnesium"):
            system_unit = app.selectbox(key=f"extension_{salt}_unit")
            system_options[salt] = system_unit.options
            self.assertEqual(set(system_unit.options), {"nM", "µM", "mM", "M"})
            self.assertEqual(system_unit.value, "mM")
        self.assertEqual(app.number_input(key="extension_sodium").value, 120.0)
        self.assertEqual(app.number_input(key="extension_magnesium").value, 0.0)
        app.radio(key="design_section").set_value("switch").run()
        self.assertFalse(app.exception)
        for salt in ("sodium", "magnesium"):
            design_unit = app.selectbox(key=f"{salt}_concentration_unit")
            self.assertEqual(system_options[salt], design_unit.options)
            self.assertEqual(design_unit.value, "mM")
        self.analyze.assert_not_called()

    def test_unit_changes_preserve_sodium_and_do_not_repeat_structure_calculation(self):
        app = self.app(extension_magnesium=7.0)
        self.fill_sequences(app)
        initial = app.session_state["system_analysis_result"]
        for unit, sodium in (("M", .12), ("µM", 120000.), ("nM", 120000000.), ("mM", 120.)):
            with self.subTest(unit=unit):
                app.selectbox(key="extension_sodium_unit").set_value(unit).run()
                self.assertFalse(app.exception)
                self.assertAlmostEqual(app.number_input(key="extension_sodium").value, sodium)
                self.assertEqual(app.session_state["system_analysis_result"], initial)
                self.analyze.assert_called_once()
        self.assertTrue(app.number_input(key="extension_magnesium").disabled)
        self.assertTrue(app.selectbox(key="extension_magnesium_unit").disabled)
        self.assertEqual(app.number_input(key="extension_magnesium").value, 0.)
        thermo = self.analyze.call_args.args[3]
        self.assertEqual((thermo.sodium_m, thermo.magnesium_m), (.12, 0.))

    def test_molar_sodium_entry_feeds_both_structure_and_extension(self):
        app = self.app()
        app.selectbox(key="extension_sodium_unit").set_value("M").run()
        app.number_input(key="extension_sodium").set_value(.18).run()
        self.fill_sequences(app)
        self.analyze.assert_called_once()
        thermo = self.analyze.call_args.args[3]
        self.assertEqual((thermo.sodium_m, thermo.magnesium_m), (.18, 0.))
        app.checkbox(key="extension_enabled").check().run()
        self.assertFalse(app.exception)
        self.assertEqual(app.radio(key="design_section").value, "extension")
        with patch("aptaswitch_studio.extension_runtime.ExtensionRunJob.start") as start:
            button(app, "Étendre le trigger").click().run()
        self.assertFalse(app.exception)
        start.assert_called_once()
        request = app.session_state["extension_job"].request
        self.assertEqual(request.thermo, thermo)
        self.analyze.assert_called_once()

    def test_custom_sodium_units_and_value_survive_results_navigation(self):
        app = self.app()
        app.selectbox(key="extension_sodium_unit").set_value("M").run()
        app.number_input(key="extension_sodium").set_value(.175).run()
        self.fill_sequences(app)
        app.radio(key="nav").set_value("Résultats").run()
        app.radio(key="nav").set_value("Conception").run()
        self.assertFalse(app.exception)
        self.assertEqual(app.selectbox(key="extension_sodium_unit").value, "M")
        self.assertEqual(app.number_input(key="extension_sodium").value, .175)
        self.assertEqual(app.number_input(key="extension_magnesium").value, 0.)
        self.analyze.assert_called_once()
        app.selectbox(key="extension_sodium_unit").set_value("mM").run()
        self.assertFalse(app.exception)
        self.assertEqual(app.number_input(key="extension_sodium").value, 175.)
        self.analyze.assert_called_once()


if __name__ == "__main__":
    unittest.main()
