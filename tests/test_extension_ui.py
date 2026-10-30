"""Portable UI regressions using a freshly recorded native technical fixture."""

from __future__ import annotations

import json
import base64
import re
import sys
import tempfile
import unittest
import xml.etree.ElementTree as ET
from contextlib import ExitStack
from pathlib import Path
from unittest.mock import patch

from tests.ui_helpers import FrenchAppTest as AppTest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from aptaswitch_core.architecture import default_architecture_for_trigger
from tests.design_fixtures import design_result_fixture
from aptaswitch_core.models import RunRequest, SequenceInput
from tests.export_fixtures import export_extension_result


APP_PATH = ROOT / "src" / "aptaswitch_studio" / "web_app.py"
ORIGINAL_TRIGGER = "GTCCAGGCTGGTATAATTAGATC"


def button_with_label(app: AppTest, label: str):
    return next(button for button in app.button if button.label == label)


class ExtensionUITests(unittest.TestCase):
    def setUp(self):
        from tests.native_ui_fixtures import patch_native_runtime
        patch_native_runtime(self)
        self.context = ExitStack()
        self.addCleanup(self.context.close)
        self.output = Path(self.context.enter_context(tempfile.TemporaryDirectory()))
        self.context.enter_context(patch("aptaswitch_core.settings.load_settings", return_value={}))
        self.context.enter_context(patch("aptaswitch_core.settings.save_settings", side_effect=AssertionError("Tests must not save user settings")))
        self.context.enter_context(patch("aptaswitch_core.thermodynamics.setup.configured_library_path", return_value="/fixture/native-hybrid"))
        self.native_library = self.context.enter_context(patch("aptaswitch_studio.extension_ui.configured_library_path", return_value="/fixture/native-hybrid"))

    def app(self) -> AppTest:
        app = AppTest.from_file(str(APP_PATH), default_timeout=20)
        app.session_state["dependency_detection_done"] = True
        self.output = Path(app.session_state["output_dir"])
        app.run()
        self.assertFalse(app.exception)
        return app

    def fill_sequences(self, app: AppTest) -> None:
        app.text_area(key="aptamer_input").set_value("ACGTACGT")
        app.text_area(key="trigger_input").set_value(ORIGINAL_TRIGGER)
        app.text_input(key="molecule_name").set_value("Mon système")
        app.run()
        self.assertFalse(app.exception)

    def open_extension(self, app: AppTest) -> None:
        if not app.checkbox(key="extension_enabled").value:
            app.checkbox(key="extension_enabled").check().run()
        else:
            app.radio(key="design_section").set_value("extension").run()
        self.assertFalse(app.exception)

    def open_system(self, app: AppTest) -> None:
        app.radio(key="design_section").set_value("system").run()
        self.assertFalse(app.exception)

    def completed_result(self) -> tuple[dict, dict[str, str]]:
        result = json.loads((ROOT / "tests" / "fixtures" / "extension_native_fixture.json").read_text())
        result["output_dir"] = str(self.output / result["run_id"])
        return result, export_extension_result(result)

    def show_completed_result(self, app: AppTest) -> dict:
        result, exports = self.completed_result()
        app.session_state["extension_result"] = result
        app.session_state["extension_exports"] = exports
        app.session_state["results_kind"] = "Extension de triggers"
        app.radio(key="nav").set_value("Résultats").run()
        self.assertFalse(app.exception)
        return result

    def test_extension_is_opt_in_and_has_a_separate_page(self):
        app = self.app()
        self.assertFalse(app.checkbox(key="extension_enabled").value)
        self.assertEqual(app.radio(key="design_section").value, "system")
        self.assertEqual(app.radio(key="design_section").options,
                         ["1 · Aptamère et trigger", "2 · Conception du toehold switch"])
        self.assertNotIn("extension_mode", [widget.key for widget in app.radio])
        self.assertNotIn("inverse_extension_target_length", [widget.key for widget in app.number_input])
        self.assertFalse(any(button.label == "Lancer la conception" for button in app.button))
        self.fill_sequences(app)
        with patch("aptaswitch_studio.extension_runtime.ExtensionRunJob.start") as extension_start, \
                patch("aptaswitch_studio.web_runtime.DesignRunJob.start") as design_start:
            self.open_extension(app)
            self.assertEqual(app.radio(key="design_section").value, "extension")
            self.assertFalse(app.radio(key="inverse_extension_enabled").value)
            self.assertNotIn("trigger_input", [widget.key for widget in app.text_area])
            self.assertNotIn("extension_temperature", [widget.key for widget in app.number_input])
            self.assertFalse(any(button.label == "Lancer la conception" for button in app.button))
            self.assertIn(ORIGINAL_TRIGGER, [item.value for item in app.code])
            button_with_label(app, "Modifier les séquences et les conditions").click().run()
            self.assertEqual(app.radio(key="design_section").value, "system")
            self.assertEqual(app.text_area(key="trigger_input").value, ORIGINAL_TRIGGER)
            self.assertNotIn("extension_mode", [widget.key for widget in app.radio])
            app.checkbox(key="extension_enabled").uncheck().run()
            self.assertEqual(app.radio(key="design_section").value, "system")
            self.assertEqual(len(app.radio(key="design_section").options), 2)
            extension_start.assert_not_called()
            design_start.assert_not_called()
        self.assertFalse(app.exception)

    def test_new_initial_tube_analysis_uses_three_strands(self):
        app = self.app()
        self.fill_sequences(app)
        self.open_system(app)
        heading = next(item for item in app.markdown if 'Complexes parasites · système initial' in item.value)
        self.assertIn('complexes de 1 à 3 brins', heading.proto.help)
        with patch('aptaswitch_studio.system_tube_runtime.SystemTubeJob.start') as start:
            app.button(key='analyze_system_tubes').click().run()
            self.assertFalse(app.exception)
            start.assert_called_once()
            self.assertEqual(app.session_state.system_tube_job.signature[4].max_complex_size, 3)

    def test_archived_five_strand_tubes_remain_visible_with_three_strand_new_inputs(self):
        from dataclasses import asdict, replace
        from aptaswitch_core.models import ThermoConditions
        from aptaswitch_core.tube_analysis import TubeAnalysisConditions
        from aptaswitch_core.thermodynamics.setup import MODEL_ID
        from aptaswitch_studio.design_session import _StoredTubeConditions, _StoredTubeJob
        from aptaswitch_studio.web_runtime import JobSnapshot
        thermo = ThermoConditions(material=MODEL_ID, switch_material='dna', trigger_material='dna')
        conditions = TubeAnalysisConditions()
        saved_conditions = _StoredTubeConditions(**{**asdict(conditions), 'max_complex_size': 5})
        result = {'status': 'completed', 'conditions': asdict(saved_conditions), 'tubes': {}}
        old_signature = ('AAAA', 'TTTT', '/fixture/native-hybrid',
                         replace(thermo, screening_max_complex_size=5, validation_max_complex_size=5),
                         saved_conditions, None)
        job = _StoredTubeJob(old_signature, JobSnapshot('done', 2, 2, '', result=result))
        app = AppTest.from_string('from aptaswitch_studio.extension_ui import _system_tube_monitor\n_system_tube_monitor()')
        app.session_state.system_tube_job = job
        app.session_state.system_tube_signature = ('AAAA', 'TTTT', '/fixture/native-hybrid', thermo, conditions, None)
        with patch('aptaswitch_studio.extension_ui.render_tube_analysis') as render:
            app.run()
            self.assertFalse(app.exception)
            render.assert_called_once_with(result, key='system_tubes')
            self.assertTrue(any('Résultat archivé' in item.value for item in app.caption))
            self.assertEqual(result['conditions']['max_complex_size'], 5)

    def test_extension_methods_are_exclusive_and_preserve_values_without_starting_work(self):
        app = self.app()
        self.fill_sequences(app)
        with patch("aptaswitch_studio.extension_runtime.ExtensionRunJob.start") as extension_start, \
                patch("aptaswitch_studio.web_runtime.DesignRunJob.start") as design_start:
            self.open_extension(app)
            app.radio(key="extension_mode").set_value("Longueur personnalisée").run()
            app.number_input(key="extension_target_length").set_value(29)
            app.radio(key="extension_side").set_value("both").run()
            app.radio(key="inverse_extension_enabled").set_value(True).run()
            self.assertNotIn("extension_mode", [widget.key for widget in app.radio])
            self.assertFalse(any(button.label == "Étendre le trigger" for button in app.button))
            app.number_input(key="inverse_extension_target_length").set_value(28)
            app.number_input(key="inverse_extension_top").set_value(12)
            app.radio(key="inverse_extension_side").set_value("3prime").run()
            button_with_label(app, "Continuer vers la conception du switch").click().run()
            self.assertEqual(app.radio(key="design_section").value, "switch")
            self.assertTrue(app.session_state["inverse_extension_enabled"])
            self.assertNotIn("inverse_extension_target_length", [widget.key for widget in app.number_input])
            self.assertNotIn("trigger_input", [widget.key for widget in app.text_area])
            self.assertEqual(app.session_state["trigger_input"], ORIGINAL_TRIGGER)
            self.open_extension(app)
            self.assertEqual(app.number_input(key="inverse_extension_target_length").value, 28)
            self.assertEqual(app.number_input(key="inverse_extension_top").value, 12)
            self.assertEqual(app.radio(key="inverse_extension_side").value, "3prime")
            app.radio(key="inverse_extension_enabled").set_value(False).run()
            self.assertNotIn("inverse_extension_target_length", [widget.key for widget in app.number_input])
            self.assertEqual(app.number_input(key="extension_target_length").value, 29)
            self.assertEqual(app.radio(key="extension_side").value, "both")
            app.radio(key="inverse_extension_enabled").set_value(True).run()
            app.checkbox(key="extension_enabled").uncheck().run()
            self.assertFalse(app.session_state["inverse_extension_enabled"])
            self.open_extension(app)
            self.assertFalse(app.radio(key="inverse_extension_enabled").value)
            self.assertEqual(app.number_input(key="extension_target_length").value, 29)
            self.assertIsNone(app.session_state["extension_job"])
            self.assertIsNone(app.session_state["design_job"])
            extension_start.assert_not_called()
            design_start.assert_not_called()
        self.assertFalse(app.exception)

    def test_initial_controls_use_reference_conditions_and_require_sequences(self):
        app = self.app()
        self.assertNotIn("extension_material", app.session_state)
        self.assertFalse(any(widget.label == "Chimie des brins" for widget in app.selectbox))
        self.assertEqual(app.number_input(key="extension_temperature").value, 37)
        self.assertEqual(app.number_input(key="extension_sodium").value, 120)
        self.assertEqual(app.number_input(key="extension_magnesium").value, 0)
        self.open_extension(app)
        self.assertTrue(button_with_label(app, "Étendre le trigger").disabled)
        self.assertEqual(app.radio(key="extension_mode").value, "Deux versions : 24 et 30 nt")
        self.open_system(app)
        self.fill_sequences(app)
        self.open_extension(app)
        self.assertFalse(button_with_label(app, "Étendre le trigger").disabled)

    def test_missing_bundled_runtime_disables_extension(self):
        self.native_library.return_value = None
        app = self.app()
        self.fill_sequences(app)
        self.open_extension(app)
        self.assertTrue(button_with_label(app, "Étendre le trigger").disabled)
        self.assertTrue(any("moteur thermodynamique intégré est indisponible" in caption.value for caption in app.caption))
        self.assertIsNone(app.session_state["extension_job"])

    def test_short_trigger_warning_requires_valid_input_and_survives_section_navigation(self):
        app = self.app()

        def short_warnings():
            return [warning.value for warning in app.warning
                    if "Attention" in warning.value and "trigger" in warning.value
                    and ("court" in warning.value or "faible" in warning.value)]

        self.assertFalse(short_warnings())
        self.fill_sequences(app)
        self.assertEqual(len(short_warnings()), 1)
        self.assertIn("activer l’extension", short_warnings()[0])
        app.radio(key="design_section").set_value("switch").run()
        self.assertNotIn("customize_tswitch", {item.key for item in app.checkbox})
        self.open_system(app)
        self.assertEqual(len(short_warnings()), 1)
        for value in ("ACGTX", "", "ACGU " * 6):
            with self.subTest(trigger=value):
                app.text_area(key="trigger_input").set_value(value).run()
                self.assertFalse(app.exception)
                self.assertFalse(short_warnings())

    def test_extension_controls_and_preview_stay_available_for_standard_length_triggers(self):
        app = self.app()
        self.fill_sequences(app)
        self.open_extension(app)
        app.radio(key="extension_mode").set_value("Longueur personnalisée").run()
        app.checkbox(key="extension_allow_long").check().run()
        app.number_input(key="extension_target_length").set_value(32)
        app.radio(key="extension_side").set_value("both").run()
        for trigger in (" acgu\n" * 6, "A" * 25, "A" * 29, "ACGT" * 7 + "AC",
                        "ACGT" * 7 + "ACG", ORIGINAL_TRIGGER):
            with self.subTest(trigger=trigger):
                self.open_system(app)
                app.text_area(key="trigger_input").set_value(trigger).run()
                self.assertFalse(app.text_area(key="trigger_input").disabled)
                self.open_extension(app)
                self.assertFalse(app.exception)
                self.assertTrue(any("Extension du trigger" in heading.value for heading in app.subheader))
                for widget in (
                    app.radio(key="extension_side"), app.radio(key="extension_mode"),
                    app.checkbox(key="extension_allow_long"),
                    app.number_input(key="extension_target_length"),
                    app.number_input(key="extension_top"),
                    app.number_input(key="extension_ensemble_cap"),
                    button_with_label(app, "Étendre le trigger"),
                ):
                    self.assertFalse(widget.disabled, widget.label)
                self.assertEqual(app.radio(key="extension_side").value, "both")
                self.assertEqual(app.number_input(key="extension_target_length").value, 32)
                preview_visible = any("Aperçu en direct de l’extension" in item.value for item in app.markdown)
                self.assertTrue(preview_visible)
                self.assertTrue(any("extensions à tester" in item.value for item in app.caption))
                self.assertTrue(any("au-delà de 30 nt" in item.value for item in app.warning))

    def test_start_callback_accepts_extensions_for_all_24_through_30_base_triggers(self):
        for length in range(24, 31):
            with self.subTest(length=length):
                app = AppTest.from_string("""
import streamlit as st
from aptaswitch_studio.extension_ui import EXTENSION_DEFAULTS, _start_extension
for key, value in EXTENSION_DEFAULTS.items():
    st.session_state.setdefault(key, value)
st.button("Force extension callback", on_click=_start_extension)
""")
                for key, value in {
                    "aptamer_input": "ACGTACGT", "trigger_input": "A" * length,
                    "molecule_name": "Test", "output_dir": str(self.output),
                    "extension_mode": "Longueur personnalisée", "extension_target_length": 31,
                    "extension_allow_long": True,
                }.items():
                    app.session_state[key] = value
                app.run()
                with patch("aptaswitch_studio.extension_runtime.ExtensionRunJob.start") as start:
                    button_with_label(app, "Force extension callback").click().run()
                self.assertFalse(app.exception)
                start.assert_called_once()
                request = app.session_state["extension_job"].request
                self.assertEqual(request.trigger_dna, "A" * length)
                self.assertEqual(request.target_lengths, (31,))
                self.assertEqual(app.session_state["extension_error"], "")

    def test_standard_targets_extend_24_through_29_bases_only_to_30(self):
        app = self.app()
        self.fill_sequences(app)
        for length in range(24, 30):
            with self.subTest(length=length):
                if app.session_state["extension_enabled"]:
                    self.open_system(app)
                app.text_area(key="trigger_input").set_value("A" * length).run()
                self.open_extension(app)
                self.assertFalse(app.exception)
                self.assertFalse(button_with_label(app, "Étendre le trigger").disabled)
                targets = [item.value for item in app.caption if "extensions à tester" in item.value]
                self.assertEqual(len(targets), 1)
                self.assertTrue(targets[0].startswith("30 nt :"))
                self.assertTrue(any("Trigger final · 30 nt" in item.value for item in app.markdown))
        with patch("aptaswitch_studio.extension_runtime.ExtensionRunJob.start") as start:
            button_with_label(app, "Étendre le trigger").click().run()
        start.assert_called_once()
        self.assertEqual(app.session_state["extension_job"].request.target_lengths, (30,))

    def test_30_base_trigger_can_extend_with_explicit_long_length_option(self):
        app = self.app()
        self.fill_sequences(app)
        app.text_area(key="trigger_input").set_value("A" * 30).run()
        self.open_extension(app)
        self.assertTrue(button_with_label(app, "Étendre le trigger").disabled)
        self.assertFalse(app.radio(key="extension_mode").disabled)
        self.assertTrue(any("Autoriser une longueur supérieure à 30 nt" in item.value for item in app.info))
        app.radio(key="extension_mode").set_value("Longueur personnalisée").run()
        self.assertFalse(app.checkbox(key="extension_allow_long").disabled)
        self.assertEqual(app.number_input(key="extension_target_length").max, 30)
        app.number_input(key="extension_target_length").set_value(30).run()
        self.assertTrue(button_with_label(app, "Étendre le trigger").disabled)
        app.checkbox(key="extension_allow_long").check().run()
        app.number_input(key="extension_target_length").set_value(31).run()
        self.assertFalse(button_with_label(app, "Étendre le trigger").disabled)
        with patch("aptaswitch_studio.extension_runtime.ExtensionRunJob.start") as start:
            button_with_label(app, "Étendre le trigger").click().run()
        self.assertFalse(app.exception)
        start.assert_called_once()
        self.assertEqual(app.session_state["extension_job"].request.target_lengths, (31,))

    def test_callback_requires_growth_and_explicit_permission_above_30(self):
        for original, target, allow_long in ((24, 24, False), (30, 29, False), (30, 31, False)):
            with self.subTest(original=original, target=target):
                app = AppTest.from_string("""
import streamlit as st
from aptaswitch_studio.extension_ui import EXTENSION_DEFAULTS, _start_extension
for key, value in EXTENSION_DEFAULTS.items():
    st.session_state.setdefault(key, value)
st.button("Force extension callback", on_click=_start_extension)
""")
                for key, value in {
                    "aptamer_input": "ACGTACGT", "trigger_input": "A" * original,
                    "molecule_name": "Test", "output_dir": str(self.output),
                    "extension_mode": "Longueur personnalisée", "extension_target_length": target,
                    "extension_allow_long": allow_long,
                }.items():
                    app.session_state[key] = value
                app.run()
                with patch("aptaswitch_studio.extension_runtime.ExtensionRunJob.start") as start:
                    button_with_label(app, "Force extension callback").click().run()
                self.assertFalse(app.exception)
                start.assert_not_called()
                self.assertIsNone(app.session_state["extension_job"])
                self.assertIn("supérieure" if target <= original else "dépassement", app.session_state["extension_error"])

    def test_ligand_selection_survives_navigation_and_is_saved_with_run(self):
        app = self.app()
        self.fill_sequences(app)
        app.checkbox(key="extension_ligand_base_1").check()
        app.checkbox(key="extension_ligand_base_4").check().run()
        self.assertEqual(app.session_state["extension_ligand_positions"], [1, 4])
        self.open_extension(app)
        app.radio(key="extension_mode").set_value("Longueur personnalisée").run()
        app.number_input(key="extension_target_length").set_value(len(ORIGINAL_TRIGGER) + 6)
        app.radio(key="extension_side").set_value("both").run()
        previews = [base64.b64decode(match).decode() for element in app.get("html")
                    for match in re.findall(r'data:image/svg\+xml;base64,([^" ]+)', element.proto.body)]
        preview = next(svg for svg in previews if 'data-domain="ligand"' in svg and 'data-domain="extension-5prime"' in svg)
        self.assertEqual(preview.count('data-domain="extension-5prime"'), 3)
        self.assertEqual(preview.count('data-domain="extension-3prime"'), 3)
        self.assertIn('data-domain="ligand" data-position="4"', preview)
        self.open_system(app)
        self.assertIn("T", app.checkbox(key="extension_ligand_base_4").label)
        app.radio(key="nav").set_value("Résultats").run()
        app.radio(key="nav").set_value("Conception").run()
        self.assertTrue(app.checkbox(key="extension_ligand_base_1").value)
        self.assertTrue(app.checkbox(key="extension_ligand_base_4").value)
        self.assertFalse(any("Cadre bleu" in item.value for item in app.caption))
        self.open_extension(app)
        with patch("aptaswitch_studio.extension_runtime.ExtensionRunJob.start"):
            button_with_label(app, "Étendre le trigger").click().run()
        self.assertFalse(app.exception)
        self.assertEqual(app.session_state["extension_job"].request.ligand_aptamer_positions, (1, 4))
        thermo = app.session_state["extension_job"].request.thermo
        self.assertEqual((thermo.switch_material, thermo.trigger_material), ("dna", "dna"))

    def test_ligand_selection_resets_only_when_aptamer_sequence_changes(self):
        app = self.app()
        self.fill_sequences(app)
        app.checkbox(key="extension_ligand_base_3").check().run()
        app.text_area(key="aptamer_input").set_value("acgu acgu").run()
        self.assertTrue(app.checkbox(key="extension_ligand_base_3").value)
        app.text_area(key="trigger_input").set_value("ACGT").run()
        self.assertTrue(app.checkbox(key="extension_ligand_base_3").value)
        app.text_area(key="aptamer_input").set_value("GGGGGGGG").run()
        self.assertEqual(app.session_state["extension_ligand_positions"], [])
        self.assertFalse(app.checkbox(key="extension_ligand_base_3").value)
        app.checkbox(key="extension_ligand_base_8").check().run()
        button_with_label(app, "Effacer la sélection du ligand").click().run()
        self.assertEqual(app.session_state["extension_ligand_positions"], [])
        self.assertFalse(app.checkbox(key="extension_ligand_base_8").value)
        self.assertFalse(app.exception)

    def test_extension_side_choice_persists_and_reaches_background_request(self):
        for side, expected in (("3prime", "+0 bases en 5′ et +5 bases en 3′"), ("both", "+3 bases en 5′ et +2 bases en 3′")):
            with self.subTest(side=side):
                app = self.app()
                self.fill_sequences(app)
                self.open_extension(app)
                self.assertEqual(app.radio(key="extension_side").value, "5prime")
                app.radio(key="extension_mode").set_value("Longueur personnalisée").run()
                app.number_input(key="extension_target_length").set_value(len(ORIGINAL_TRIGGER) + 5)
                app.radio(key="extension_side").set_value(side).run()
                self.assertTrue(any(expected in caption.value for caption in app.caption))
                if side == "both":
                    self.assertIn("impair", app.radio(key="extension_side").help)
                app.radio(key="nav").set_value("Résultats").run()
                app.radio(key="nav").set_value("Conception").run()
                self.assertEqual(app.radio(key="extension_side").value, side)
                with patch("aptaswitch_studio.extension_runtime.ExtensionRunJob.start"):
                    button_with_label(app, "Étendre le trigger").click().run()
                self.assertFalse(app.exception)
                self.assertEqual(app.session_state["extension_job"].request.extension_side, side)

    def test_both_ends_result_displays_suffix_and_reuses_complete_trigger(self):
        # Pure display fixture: preserve the stored structures and test suffix placement.
        result = json.loads((ROOT / "tests" / "fixtures" / "extension_native_fixture.json").read_text())
        result["input"].update(extension_side="both", target_lengths=[6], extension_lengths=[2])
        result["candidates"] = [row for row in result["candidates"] if row["trigger_length"] == 6]
        for index, row in enumerate(result["candidates"], 1):
            extension = row["extension_5prime"]
            row.update(rank=index, extension_side="both", extension_5prime=extension[:1],
                       extension_3prime=extension[1:], extension_5prime_length=1, extension_3prime_length=1,
                       extended_trigger=extension[:1] + "ACGT" + extension[1:])
        app = self.app()
        app.session_state["extension_result"] = result
        app.session_state["results_kind"] = "Extension de triggers"
        app.radio(key="nav").set_value("Résultats").run()
        self.assertFalse(app.exception)
        self.assertTrue(any("Moitié 5′ / moitié 3′" in caption.value for caption in app.caption))
        table = app.get("table")[0].value
        self.assertEqual(list(table["Extension 3′"]), [c["extension_3prime"] for c in result["candidates"]])
        self.assertEqual(list(table["Extension 5′"]), [c["extension_5prime"] for c in result["candidates"]])
        app.button(key="extension_use_2").click().run()
        self.assertFalse(app.exception)
        candidate = result["candidates"][0]
        self.assertEqual(app.radio(key="design_section").value, "switch")
        self.assertEqual(app.session_state["trigger_input"], candidate["extension_5prime"] + "ACGT" + candidate["extension_3prime"])
        self.open_system(app)
        self.assertEqual(app.text_area(key="trigger_input").value, candidate["extension_5prime"] + "ACGT" + candidate["extension_3prime"])

    def test_custom_length_is_limited_to_30_with_explicit_long_run_warning(self):
        app = self.app()
        self.fill_sequences(app)
        self.open_extension(app)
        app.radio(key="extension_mode").set_value("Longueur personnalisée").run()
        self.assertEqual(app.number_input(key="extension_target_length").max, 30)
        app.checkbox(key="extension_allow_long").check().run()
        app.number_input(key="extension_target_length").set_value(31).run()
        self.assertTrue(any("au-delà de 30 nt" in warning.value for warning in app.warning))
        self.assertFalse(button_with_label(app, "Étendre le trigger").disabled)
        app.checkbox(key="extension_allow_long").uncheck().run()
        self.assertEqual(app.number_input(key="extension_target_length").value, 30)
        self.assertEqual(app.number_input(key="extension_target_length").max, 30)
        self.assertFalse(app.exception)

    def test_start_opens_results_and_keeps_original_trigger_until_selection(self):
        app = self.app()
        self.fill_sequences(app)
        self.open_extension(app)
        with patch("aptaswitch_studio.extension_runtime.ExtensionRunJob.start") as start:
            button_with_label(app, "Étendre le trigger").click().run()
        self.assertFalse(app.exception)
        start.assert_called_once()
        self.assertEqual(app.radio(key="nav").value, "Résultats")
        self.assertEqual(app.radio(key="results_kind").value, "Extension de triggers")
        self.assertEqual(app.session_state["trigger_input"], ORIGINAL_TRIGGER)
        job = app.session_state["extension_job"]
        self.assertEqual(job.request.trigger_dna, ORIGINAL_TRIGGER)
        self.assertEqual(job.request.target_lengths, (24, 30))
        self.assertEqual(job.request.output_dir, self.output)
        self.assertTrue(any("Préparation de l’extension" in code.value for code in app.code))
        metrics = {metric.label for metric in app.metric}
        self.assertIn("Progression de l’étape", metrics)
        self.assertIn("Temps restant de l’étape", metrics)
        button_with_label(app, "Annuler l’extension").click().run()
        self.assertTrue(job._cancel.is_set())

    def test_real_result_renders_three_structures_curves_candidates_and_downloads(self):
        app = self.app()
        result = self.show_completed_result(app)
        downloads = [button.label for button in app.get("download_button")]
        self.assertEqual(downloads.count("Structure SVG"), 3)
        self.assertEqual(downloads.count("Structure PNG"), 3)
        self.assertEqual(downloads.count("Séquence linéaire SVG"), 3)
        self.assertEqual(downloads.count("Séquence linéaire PNG"), 3)
        self.assertEqual(downloads.count("Paysage de sélection SVG"), 2)
        self.assertEqual(downloads.count("Paysage de sélection PNG"), 2)
        self.assertEqual(downloads.count("Probabilités des paires SVG"), 2)
        self.assertEqual(downloads.count("Probabilités des paires PNG"), 2)
        self.assertIn("Tous les exports d’extension (.zip)", downloads)
        self.assertEqual(sum(len(popover.get("download_button")) for popover in app.get("popover")), len(downloads))
        self.assertFalse(any("dot-bracket" in expander.label for expander in app.expander))
        selected = [result["reference"]] + [
            next(candidate for candidate in result["candidates"] if candidate["trigger_length"] == length)
            for length in sorted({candidate["trigger_length"] for candidate in result["candidates"]})
        ]
        self.assertEqual(len(app.tabs), len(selected) * 2)
        for tab, data in zip(app.tabs, (candidate for candidate in selected for _ in range(2))):
            self.assertEqual([popover.proto.popover.label for popover in tab.get("popover")], ["Structure", "Exporter"])
            self.assertEqual([code.value for code in tab.code], [data["mfe_structure"]])
            drawing = next(base64.b64decode(match).decode() for element in tab.get("html")
                           for match in re.findall(r'data:image/svg\+xml;base64,([^" ]+)', element.proto.body))
            texts = " ".join(node.text or "" for node in ET.fromstring(drawing).iter("{http://www.w3.org/2000/svg}text"))
            self.assertNotIn(data["mfe_structure"].replace("+", " "), texts)
            self.assertNotIn("ΔG", texts)
            self.assertNotIn("kcal/mol", texts)
        self.assertEqual(len(app.get("table")), 1)
        table = app.get("table")[0].value
        self.assertEqual(len(table), len(result["candidates"]))
        self.assertIn("Trigger étendu 5′ → 3′", table.columns)
        self.assertIn("Structure 2D (dot-bracket)", table.columns)
        self.assertEqual(len(app.selectbox(key="extension_selected_5").options), 2)
        self.assertEqual(len(app.selectbox(key="extension_selected_6").options), 2)

    def test_advance_button_continues_without_cancelling(self):
        app = self.app()
        self.fill_sequences(app)
        self.open_extension(app)
        with patch("aptaswitch_studio.extension_runtime.ExtensionRunJob.start"):
            button_with_label(app, "Étendre le trigger").click().run()
        job = app.session_state["extension_job"]
        self.assertTrue(app.button(key="advance_extension").disabled)
        job._set_phase("mfe")
        job._progress(100, 1000, "MFE")
        app.run()
        self.assertFalse(app.button(key="advance_extension").disabled)
        app.button(key="advance_extension").click().run()
        self.assertTrue(job._advance.is_set())
        self.assertFalse(job._cancel.is_set())
        self.assertTrue(app.button(key="advance_extension").disabled)
        self.assertFalse(app.exception)

    def test_resume_cancelled_run_uses_original_sequences_not_edited_form(self):
        app = self.app()
        self.fill_sequences(app)
        self.open_extension(app)
        with patch("aptaswitch_studio.extension_runtime.ExtensionRunJob.start"):
            button_with_label(app, "Étendre le trigger").click().run()
        original_request = app.session_state["extension_job"].request
        result, exports = self.completed_result()
        result["request"] = original_request.to_dict()
        result["status"] = "cancelled"
        result["mfe_screen"].update(mfe_survivors_for_ensemble=2, ensemble_complete=False)
        app.session_state["extension_job"] = None
        app.session_state["extension_result"] = result
        app.session_state["extension_exports"] = exports
        app.session_state["trigger_input"] = "AAAACCCCGGGG"
        app.run()
        with patch("aptaswitch_studio.extension_runtime.ExtensionRunJob.start") as start:
            app.button(key="resume_extension").click().run()
        start.assert_called_once()
        job = app.session_state["extension_job"]
        self.assertEqual(job.request.trigger_dna, ORIGINAL_TRIGGER)
        self.assertEqual(app.session_state["trigger_input"], "AAAACCCCGGGG")
        self.assertEqual(job._resume_result["status"], "cancelled")
        self.assertFalse(app.exception)

    def test_selected_candidate_is_reused_for_switch_design(self):
        app = self.app()
        result = self.show_completed_result(app)
        result["input"]["ligand_aptamer_positions"] = [1, 8]
        candidates = [candidate for candidate in result["candidates"] if candidate["trigger_length"] == 5]
        app.selectbox(key="extension_selected_5").select_index(1).run()
        self.assertFalse(app.exception)
        for tab in app.tabs[2:4]:
            self.assertEqual([code.value for code in tab.code], [candidates[1]["mfe_structure"]])
        app.button(key="extension_use_1").click().run()
        self.assertFalse(app.exception)
        self.assertEqual(app.radio(key="nav").value, "Conception")
        self.assertEqual(app.session_state["design_section"], "switch")
        self.assertFalse(app.checkbox(key="extension_enabled").value)
        self.assertFalse(app.session_state["inverse_extension_enabled"])
        self.assertNotIn("trigger_input", [widget.key for widget in app.text_area])
        self.open_system(app)
        self.assertEqual(app.text_area(key="trigger_input").value, candidates[1]["extended_trigger"])
        self.assertEqual(app.text_area(key="aptamer_input").value, result["input"]["aptamer_dna"])
        self.assertEqual(app.text_input(key="molecule_name").value, result["request"]["molecule_name"])
        self.assertTrue(app.checkbox(key="extension_ligand_base_1").value)
        self.assertTrue(app.checkbox(key="extension_ligand_base_8").value)

    def test_navigation_preserves_extension_parameters_and_existing_switch_results(self):
        app = self.app()
        self.fill_sequences(app)
        app.number_input(key="extension_temperature").set_value(31.5)
        app.number_input(key="extension_sodium").set_value(190.0)
        app.run()
        self.open_extension(app)
        app.radio(key="extension_mode").set_value("Longueur personnalisée").run()
        app.number_input(key="extension_target_length").set_value(29)
        app.number_input(key="extension_top").set_value(25)
        app.run()
        design = design_result_fixture(RunRequest(
            sequences=SequenceInput(aptamer_dna="ACGTACGT", trigger_dna=ORIGINAL_TRIGGER),
            architecture=default_architecture_for_trigger(len(ORIGINAL_TRIGGER)),
            trials=1, output_dir=self.output,
        ))
        app.session_state["last_result"] = design
        self.show_completed_result(app)
        app.radio(key="results_kind").set_value("Design de switches").run()
        self.assertFalse(app.dataframe)
        self.assertTrue(any("Rang" in table.value.columns for table in app.table))
        self.assertEqual(app.session_state["last_result"].run_id, design.run_id)
        app.radio(key="nav").set_value("Conception").run()
        self.assertEqual(app.number_input(key="extension_target_length").value, 29)
        self.assertEqual(app.number_input(key="extension_top").value, 25)
        self.open_system(app)
        self.assertEqual(app.number_input(key="extension_temperature").value, 31.5)
        self.assertEqual(app.number_input(key="extension_sodium").value, 190)
        self.assertNotIn("extension_material", app.session_state)
        self.assertEqual(app.text_area(key="trigger_input").value, ORIGINAL_TRIGGER)
        self.assertFalse(app.exception)


if __name__ == "__main__":
    unittest.main()
