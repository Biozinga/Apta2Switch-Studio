"""Portable UI regressions using a saved real-NUPACK result and no user settings."""

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
from aptaswitch_studio.extension_runtime import export_extension_result


APP_PATH = ROOT / "src" / "aptaswitch_studio" / "web_app.py"
ORIGINAL_TRIGGER = "GTCCAGGCTGGTATAATTAGATC"


def button_with_label(app: AppTest, label: str):
    return next(button for button in app.button if button.label == label)


class ExtensionUITests(unittest.TestCase):
    def setUp(self):
        self.context = ExitStack()
        self.addCleanup(self.context.close)
        self.output = Path(self.context.enter_context(tempfile.TemporaryDirectory()))
        self.context.enter_context(patch("aptaswitch_core.nupack_setup.load_settings", return_value={}))
        self.context.enter_context(patch("aptaswitch_core.nupack_setup.detect_nupack", return_value=None))
        self.context.enter_context(patch("aptaswitch_core.nupack_setup.save_settings", side_effect=AssertionError("Tests must not save user settings")))
        self.context.enter_context(patch("aptaswitch_core.nupack_setup.save_nupack_settings", side_effect=AssertionError("Tests must not save NUPACK settings")))
        self.context.enter_context(patch("aptaswitch_core.nupack_setup.configured_nupack_import_path", return_value="/fixture/nupack"))
        self.nupack = self.context.enter_context(patch("aptaswitch_studio.extension_ui.configured_nupack_import_path", return_value="/fixture/nupack"))

    def app(self) -> AppTest:
        app = AppTest.from_file(str(APP_PATH), default_timeout=20)
        app.session_state["dependency_detection_done"] = True
        app.session_state["output_dir"] = str(self.output)
        app.run()
        self.assertFalse(app.exception)
        return app

    def fill_sequences(self, app: AppTest) -> None:
        app.text_area(key="aptamer_input").set_value("ACGTACGT")
        app.text_area(key="trigger_input").set_value(ORIGINAL_TRIGGER)
        app.text_input(key="molecule_name").set_value("Mon système")
        app.run()
        self.assertFalse(app.exception)

    def completed_result(self) -> tuple[dict, dict[str, str]]:
        result = json.loads((ROOT / "tests" / "fixtures" / "extension_nupack_smoke.json").read_text())
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

    def test_initial_controls_use_reference_conditions_and_require_sequences(self):
        app = self.app()
        self.assertTrue(button_with_label(app, "Étendre le trigger").disabled)
        self.assertEqual(app.radio(key="extension_mode").value, "Deux versions : 24 et 30 nt")
        self.assertEqual(app.selectbox(key="extension_material").value, "ADN · dna04")
        self.assertEqual(app.number_input(key="extension_temperature").value, 37)
        self.assertEqual(app.number_input(key="extension_sodium").value, 120)
        self.assertEqual(app.number_input(key="extension_magnesium").value, 0)
        self.fill_sequences(app)
        self.assertFalse(button_with_label(app, "Étendre le trigger").disabled)

    def test_missing_nupack_disables_extension(self):
        self.nupack.return_value = None
        app = self.app()
        self.fill_sequences(app)
        self.assertTrue(button_with_label(app, "Étendre le trigger").disabled)
        self.assertTrue(any("NUPACK est requis" in caption.value for caption in app.caption))
        self.assertIsNone(app.session_state["extension_job"])

    def test_short_trigger_warning_requires_valid_input_and_survives_customization(self):
        app = self.app()

        def short_warnings():
            return [warning.value for warning in app.warning
                    if "Attention" in warning.value and "trigger" in warning.value
                    and ("court" in warning.value or "faible" in warning.value)]

        self.assertFalse(short_warnings())
        self.fill_sequences(app)
        self.assertEqual(len(short_warnings()), 1)
        self.assertIn("étendre", short_warnings()[0])
        app.checkbox(key="customize_tswitch").check().run()
        self.assertEqual(len(short_warnings()), 1)
        for value in ("ACGTX", "", "ACGU " * 6):
            with self.subTest(trigger=value):
                app.text_area(key="trigger_input").set_value(value).run()
                self.assertFalse(app.exception)
                self.assertFalse(short_warnings())

    def test_extension_stays_visible_and_disabled_from_24_through_30_bases(self):
        app = self.app()
        self.fill_sequences(app)
        app.radio(key="extension_mode").set_value("Longueur personnalisée").run()
        app.checkbox(key="extension_allow_long").check().run()
        app.number_input(key="extension_target_length").set_value(32)
        app.radio(key="extension_side").set_value("both").run()
        for trigger, disabled in ((" acgu\n" * 6, True), ("ACGT" * 7 + "AC", True),
                                  ("ACGT" * 7 + "ACG", False), (ORIGINAL_TRIGGER, False)):
            with self.subTest(trigger=trigger, disabled=disabled):
                app.text_area(key="trigger_input").set_value(trigger).run()
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
                    self.assertEqual(widget.disabled, disabled, widget.label)
                self.assertFalse(app.text_area(key="trigger_input").disabled)
                self.assertFalse(app.selectbox(key="extension_material").disabled)
                self.assertEqual(app.radio(key="extension_side").value, "both")
                self.assertEqual(app.number_input(key="extension_target_length").value, 32)
                preview_visible = any("Aperçu en direct de l’extension" in item.value for item in app.markdown)
                self.assertEqual(preview_visible, not disabled)
                if disabled:
                    self.assertFalse(any("extensions à tester" in item.value for item in app.caption))
                    self.assertFalse(any("au-delà de 30 nt" in item.value or "recherche exhaustive" in item.value
                                         for item in app.warning))
                    self.assertTrue(any("24" in item.value and "30" in item.value for item in app.caption))

    def test_start_callback_rejects_standard_length_even_with_long_extension_enabled(self):
        for length in (24, 30):
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
                start.assert_not_called()
                self.assertIsNone(app.session_state["extension_job"])
                self.assertIn("24", app.session_state["extension_error"])
                self.assertIn("30", app.session_state["extension_error"])

    def test_ligand_selection_survives_navigation_and_is_saved_with_run(self):
        app = self.app()
        self.fill_sequences(app)
        app.checkbox(key="extension_ligand_base_1").check()
        app.checkbox(key="extension_ligand_base_4").check().run()
        self.assertEqual(app.session_state["extension_ligand_positions"], [1, 4])
        app.radio(key="extension_mode").set_value("Longueur personnalisée").run()
        app.number_input(key="extension_target_length").set_value(len(ORIGINAL_TRIGGER) + 6)
        app.radio(key="extension_side").set_value("both").run()
        previews = [base64.b64decode(match).decode() for element in app.get("html")
                    for match in re.findall(r'data:image/svg\+xml;base64,([^" ]+)', element.proto.body)]
        preview = next(svg for svg in previews if 'data-domain="ligand"' in svg and 'data-domain="extension-5prime"' in svg)
        self.assertEqual(preview.count('data-domain="extension-5prime"'), 3)
        self.assertEqual(preview.count('data-domain="extension-3prime"'), 3)
        self.assertIn('data-domain="ligand" data-position="4"', preview)
        app.selectbox(key="extension_material").set_value("ARN · rna06").run()
        self.assertIn("U", app.checkbox(key="extension_ligand_base_4").label)
        app.radio(key="nav").set_value("Résultats").run()
        app.radio(key="nav").set_value("Conception").run()
        self.assertTrue(app.checkbox(key="extension_ligand_base_1").value)
        self.assertTrue(app.checkbox(key="extension_ligand_base_4").value)
        self.assertFalse(any("Cadre bleu" in item.value for item in app.caption))
        with patch("aptaswitch_studio.extension_runtime.ExtensionRunJob.start"):
            button_with_label(app, "Étendre le trigger").click().run()
        self.assertFalse(app.exception)
        self.assertEqual(app.session_state["extension_job"].request.ligand_aptamer_positions, (1, 4))

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
                self.assertEqual(app.radio(key="extension_side").value, "5prime")
                self.fill_sequences(app)
                app.radio(key="extension_mode").set_value("Longueur personnalisée").run()
                app.number_input(key="extension_target_length").set_value(len(ORIGINAL_TRIGGER) + 5)
                app.radio(key="extension_side").set_value(side).run()
                self.assertTrue(any(expected in caption.value for caption in app.caption))
                if side == "both":
                    self.assertTrue(any("impair" in caption.value for caption in app.caption))
                app.radio(key="nav").set_value("Résultats").run()
                app.radio(key="nav").set_value("Conception").run()
                self.assertEqual(app.radio(key="extension_side").value, side)
                with patch("aptaswitch_studio.extension_runtime.ExtensionRunJob.start"):
                    button_with_label(app, "Étendre le trigger").click().run()
                self.assertFalse(app.exception)
                self.assertEqual(app.session_state["extension_job"].request.extension_side, side)

    def test_both_ends_result_displays_suffix_and_reuses_complete_trigger(self):
        from tests.test_trigger_extension import FakeNupack
        from aptaswitch_core import trigger_extension
        from aptaswitch_core.models import ThermoConditions
        with patch.object(trigger_extension, "nupack_import_context") as context:
            context.return_value.__enter__.return_value = FakeNupack()
            result = trigger_extension.design_trigger_extensions(
                aptamer_dna="ACGTACGT", trigger_dna="ACGT", target_lengths=[6],
                thermo=ThermoConditions(material="dna04"), import_path="/fixture", top=2, extension_side="both",
            )
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
        self.assertEqual(app.text_area(key="trigger_input").value, candidate["extension_5prime"] + "ACGT" + candidate["extension_3prime"])

    def test_custom_length_is_limited_to_30_with_explicit_long_run_warning(self):
        app = self.app()
        self.fill_sequences(app)
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
        self.assertTrue(any("Préparation de l’extension NUPACK" in code.value for code in app.code))
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
        self.assertEqual(app.text_area(key="trigger_input").value, candidates[1]["extended_trigger"])
        self.assertEqual(app.text_area(key="aptamer_input").value, result["input"]["aptamer_dna"])
        self.assertEqual(app.text_input(key="molecule_name").value, result["request"]["molecule_name"])
        self.assertEqual(app.session_state["binding_start"], 1)
        self.assertEqual(app.session_state["binding_end"], 5)
        self.assertTrue(app.checkbox(key="extension_ligand_base_1").value)
        self.assertTrue(app.checkbox(key="extension_ligand_base_8").value)

    def test_navigation_preserves_extension_parameters_and_existing_switch_results(self):
        app = self.app()
        self.fill_sequences(app)
        app.radio(key="extension_mode").set_value("Longueur personnalisée").run()
        app.number_input(key="extension_target_length").set_value(29)
        app.number_input(key="extension_temperature").set_value(31.5)
        app.number_input(key="extension_sodium").set_value(90.0)
        app.selectbox(key="extension_material").set_value("ARN · rna06")
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
        self.assertEqual(len(app.get("table")), 1)
        self.assertEqual(app.session_state["last_result"].run_id, design.run_id)
        app.radio(key="nav").set_value("Conception").run()
        self.assertEqual(app.number_input(key="extension_target_length").value, 29)
        self.assertEqual(app.number_input(key="extension_temperature").value, 31.5)
        self.assertEqual(app.number_input(key="extension_sodium").value, 90)
        self.assertEqual(app.number_input(key="extension_top").value, 25)
        self.assertEqual(app.selectbox(key="extension_material").value, "ARN · rna06")
        self.assertEqual(app.text_area(key="trigger_input").value, ORIGINAL_TRIGGER)
        self.assertFalse(app.exception)

    def test_rna_aptamer_result_cannot_be_reused_as_a_dna_switch_design(self):
        app = self.app()
        result, exports = self.completed_result()
        result["model"].update(material="rna06", aptamer_material="rna", trigger_material="rna")
        for name in ("aptamer_dna", "trigger_dna"):
            result["input"][name] = result["input"][name].replace("T", "U")
        for candidate in result["candidates"]:
            candidate["extended_trigger"] = candidate["extended_trigger"].replace("T", "U")
            candidate["extension_5prime"] = candidate["extension_5prime"].replace("T", "U")
        app.session_state["extension_result"] = result
        app.session_state["extension_exports"] = exports
        app.session_state["results_kind"] = "Extension de triggers"
        app.radio(key="nav").set_value("Résultats").run()
        self.assertFalse(app.exception)
        uses = [button for button in app.button if button.label == "Utiliser pour le design du switch"]
        self.assertEqual(len(uses), 2)
        self.assertTrue(all(button.disabled for button in uses))
        self.assertTrue(any("ADN" in caption.value and "ARN" in caption.value for caption in app.caption))


if __name__ == "__main__":
    unittest.main()
