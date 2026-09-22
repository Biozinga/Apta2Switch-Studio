"""UI regressions for the initial-system analysis and optional extension."""

from __future__ import annotations

import base64
import re
import tempfile
import unittest
import xml.etree.ElementTree as ET
from contextlib import ExitStack
from pathlib import Path
from unittest.mock import patch

from tests.ui_helpers import FrenchAppTest as AppTest
from streamlit.testing.v1.element_tree import Block

from aptaswitch_core.extension_preview import LIGAND_COLOR
from aptaswitch_studio.extension_ui import _analyze_initial_system


ROOT = Path(__file__).resolve().parents[1]
APP_PATH = ROOT / "src" / "aptaswitch_studio" / "web_app.py"
APTAMER = "ACGTACGT"
TRIGGER = "GTCCAGGCTGGTATAATTAGATC"


def initial_section(app: AppTest) -> Block:
    return next(node for node in app.main if isinstance(node, Block)
                and any(getattr(child, "type", None) == "subheader"
                        and child.value == "1 · Système initial : données et structure 2D"
                        for child in node.children.values()))


def initial_svgs(app: AppTest, *, linear: bool = False) -> list[str]:
    tabs = initial_section(app).tabs
    if not tabs:
        return []
    return [base64.b64decode(match).decode() for element in tabs[int(linear)].get("html")
            for match in re.findall(r'data:image/svg\+xml;base64,([^" ]+)', element.proto.body)]


def visible_svg_text(svg: str) -> str:
    return " ".join(element.text or "" for element in ET.fromstring(svg).iter("{http://www.w3.org/2000/svg}text"))


def button(app: AppTest, label: str):
    return next(item for item in app.button if item.label == label)


class SystemSectionsTests(unittest.TestCase):
    def setUp(self):
        self.context = ExitStack()
        self.addCleanup(self.context.close)
        self.output = Path(self.context.enter_context(tempfile.TemporaryDirectory()))
        self.context.enter_context(patch("aptaswitch_core.nupack_setup.load_settings", return_value={}))
        self.context.enter_context(patch("aptaswitch_core.nupack_setup.detect_nupack", return_value=None))
        self.context.enter_context(patch("aptaswitch_core.nupack_setup.save_settings", side_effect=AssertionError("Do not save user settings")))
        self.context.enter_context(patch("aptaswitch_core.nupack_setup.save_nupack_settings", side_effect=AssertionError("Do not save NUPACK settings")))
        self.context.enter_context(patch("aptaswitch_core.nupack_setup.configured_nupack_import_path", return_value="/fixture/nupack"))
        self.nupack = self.context.enter_context(patch("aptaswitch_studio.extension_ui.configured_nupack_import_path", return_value="/fixture/nupack"))
        self.analyze = self.context.enter_context(patch("aptaswitch_core.system_analysis.analyze_system", side_effect=self.analysis_result))
        _analyze_initial_system.clear()
        self.addCleanup(_analyze_initial_system.clear)

    @staticmethod
    def analysis_result(import_path, aptamer, trigger, thermo):
        chemistry = "rna" if thermo.material == "rna06" else "dna"
        if chemistry == "rna":
            aptamer, trigger = aptamer.replace("T", "U"), trigger.replace("T", "U")
        return {
            "aptamer_sequence": aptamer,
            "trigger_sequence": trigger,
            "mfe_structure": "((" + "." * (len(aptamer) - 2) + "+" + "." * (len(trigger) - 2) + "))",
            "mfe_energy_kcal_mol": -4.25,
            "paired_trigger_indices": [len(trigger) - 2, len(trigger) - 1],
            "model": {
                "material": thermo.material, "temperature_c": thermo.temperature_c,
                "sodium_m": thermo.sodium_m, "magnesium_m": thermo.magnesium_m,
                "ensemble": "stacking", "aptamer_material": chemistry, "trigger_material": chemistry,
            },
        }

    def app(self) -> AppTest:
        app = AppTest.from_file(str(APP_PATH), default_timeout=20)
        app.session_state["dependency_detection_done"] = True
        app.session_state["output_dir"] = str(self.output)
        app.run()
        self.assertFalse(app.exception)
        return app

    def fill_sequences(self, app: AppTest) -> None:
        app.text_area(key="aptamer_input").set_value(APTAMER)
        app.text_area(key="trigger_input").set_value(TRIGGER).run()
        self.assertFalse(app.exception)

    def test_input_analysis_and_extension_have_separate_ordered_sections(self):
        app = self.app()
        headings = [heading.value for heading in app.subheader]
        titles = ["1 · Système initial : données et structure 2D", "2 · Extension du trigger (facultatif)", "3 · Architecture du switch"]
        indices = [headings.index(title) for title in titles]
        self.assertEqual(indices, sorted(indices))
        sections = []
        for title in titles:
            matching = [node for node in app.main if isinstance(node, Block)
                        and any(getattr(child, "type", None) == "subheader" and child.value == title
                                for child in node.children.values())]
            self.assertEqual(len(matching), 1)
            sections.append(matching[0])
        initial, extension, switch = sections
        self.assertEqual({widget.key for widget in initial.text_area}, {"aptamer_input", "trigger_input"})
        self.assertEqual({widget.key for widget in initial.number_input}, {"extension_temperature", "extension_sodium", "extension_magnesium"})
        self.assertIn("extension_material", {widget.key for widget in initial.selectbox})
        self.assertFalse(initial.radio)
        self.assertIn("extension_side", {widget.key for widget in extension.radio})
        self.assertFalse(extension.text_area)
        self.assertNotIn("extension_temperature", {widget.key for widget in extension.number_input})
        self.assertNotIn("binding_start", {widget.key for widget in switch.number_input})
        self.assertIn("customize_tswitch", {widget.key for widget in switch.checkbox})

    def test_sequences_automatically_build_initial_structure_without_starting_extension(self):
        app = self.app()
        self.analyze.assert_not_called()
        self.fill_sequences(app)
        self.analyze.assert_called_once()
        result = app.session_state["system_analysis_result"]
        self.assertEqual(result["aptamer_sequence"], APTAMER)
        self.assertEqual(result["trigger_sequence"], TRIGGER)
        self.assertEqual(app.session_state["aptamer_indices"], [len(TRIGGER) - 2, len(TRIGGER) - 1])
        self.assertEqual(len(initial_svgs(app)), 1)
        self.assertNotIn(result["mfe_structure"], initial_svgs(app)[0])
        self.assertNotIn("ΔG", visible_svg_text(initial_svgs(app)[0]))
        self.assertNotIn("kcal/mol", visible_svg_text(initial_svgs(app)[0]))
        for tab in initial_section(app).tabs:
            self.assertEqual([code.value for code in tab.code], [result["mfe_structure"]])
            self.assertEqual([popover.proto.popover.label for popover in tab.get("popover")], ["Structure", "Exporter"])
            self.assertEqual(sum(len(popover.get("download_button")) for popover in tab.get("popover")), 2)
        self.assertEqual({item.label for item in app.get("download_button") if item.label.startswith("Structure initiale")},
                         {"Structure initiale SVG", "Structure initiale PNG"})
        self.assertIsNone(app.session_state["extension_job"])

    def test_ligand_changes_redraw_colors_without_repeating_folding(self):
        app = self.app()
        self.fill_sequences(app)
        previous_result = app.session_state["system_analysis_result"]
        self.assertNotIn('class="base-annotation"', initial_svgs(app)[0])
        app.checkbox(key="extension_ligand_base_3").check()
        app.checkbox(key="extension_ligand_base_7").check().run()
        self.assertFalse(app.exception)
        self.analyze.assert_called_once()
        svg = initial_svgs(app)[0]
        self.assertIn('class="base-annotation" data-position="3"', svg)
        self.assertIn('class="base-annotation" data-position="7"', svg)
        self.assertIn(f'fill="{LIGAND_COLOR}"', svg)
        self.assertEqual(app.session_state["system_analysis_result"], previous_result)
        self.assertEqual(app.session_state["extension_ligand_positions"], [3, 7])

    def test_initial_linear_tab_uses_real_sequences_and_preserves_annotations(self):
        app = self.app()
        self.fill_sequences(app)
        app.checkbox(key="extension_ligand_base_3").check().run()
        labels = [tab.label for tab in app.tabs]
        self.assertEqual(labels[:2], ["Structure 2D", "Séquence linéaire"])
        svg = initial_svgs(app, linear=True)[0]
        ns = {"s": "http://www.w3.org/2000/svg"}
        root = ET.fromstring(svg)
        groups = root.findall('.//s:g[@data-domain]', ns)
        self.assertEqual(''.join(group.find('s:text', ns).text for group in groups), APTAMER + TRIGGER)
        self.assertEqual([group.attrib['data-position'] for group in groups if group.attrib['data-domain'] == 'ligand'], ['3'])
        self.assertIn('width="20.0" height="26.0"', svg)
        self.assertNotIn("Séquences du système initial", visible_svg_text(svg))
        self.assertTrue(any(item.label == "Séquence initiale PNG" for item in app.get("download_button")))
        self.analyze.assert_called_once()

    def test_each_condition_change_recomputes_with_exact_model_and_salt_concentrations(self):
        app = self.app()
        self.fill_sequences(app)
        changes = [("extension_temperature", 26.75), ("extension_sodium", 42.0), ("extension_magnesium", 7.0)]
        for index, (key, value) in enumerate(changes, start=2):
            app.number_input(key=key).set_value(value).run()
            self.assertFalse(app.exception)
            self.assertEqual(self.analyze.call_count, index)
        thermo = self.analyze.call_args.args[3]
        self.assertEqual((thermo.material, thermo.temperature_c, thermo.sodium_m, thermo.magnesium_m), ("dna04", 26.75, .042, .007))
        app.selectbox(key="extension_material").set_value("ARN · rna06").run()
        self.assertFalse(app.exception)
        self.assertEqual(self.analyze.call_count, 5)
        thermo = self.analyze.call_args.args[3]
        self.assertEqual((thermo.material, thermo.temperature_c, thermo.sodium_m, thermo.magnesium_m), ("rna06", 26.75, .042, .007))
        self.assertEqual((thermo.switch_material, thermo.trigger_material), ("rna", "rna"))
        self.assertNotIn("T", app.session_state["system_analysis_result"]["aptamer_sequence"])
        self.assertNotIn("rna06", visible_svg_text(initial_svgs(app)[0]))
        self.assertTrue(any("rna06" in caption.value and "26.75" in caption.value for caption in app.caption))

    def test_extension_side_and_length_changes_preserve_initial_analysis(self):
        app = self.app()
        self.fill_sequences(app)
        initial = initial_svgs(app)[0]
        for side in ("3prime", "both", "5prime"):
            app.radio(key="extension_side").set_value(side).run()
        app.radio(key="extension_mode").set_value("Longueur personnalisée").run()
        app.number_input(key="extension_target_length").set_value(29).run()
        self.assertFalse(app.exception)
        self.analyze.assert_called_once()
        self.assertEqual(initial_svgs(app)[0], initial)

    def test_invalid_input_and_failed_new_sequence_clear_old_structure_and_binding(self):
        app = self.app()
        self.fill_sequences(app)
        app.text_area(key="trigger_input").set_value("ACGN").run()
        self.assertFalse(app.exception)
        self.assertIsNone(app.session_state["system_analysis_result"])
        self.assertEqual(app.session_state["aptamer_indices"], [])
        self.assertFalse(initial_svgs(app))
        self.analyze.assert_called_once()
        app.text_area(key="trigger_input").set_value(TRIGGER).run()
        self.assertTrue(initial_svgs(app))
        self.analyze.side_effect = RuntimeError("Le modèle ne peut pas calculer cette séquence")
        app.text_area(key="aptamer_input").set_value("GCGCGCGC").run()
        self.assertFalse(app.exception)
        self.assertIsNone(app.session_state["system_analysis_result"])
        self.assertEqual(app.session_state["aptamer_indices"], [])
        self.assertFalse(initial_svgs(app))
        self.assertTrue(any("Le modèle ne peut pas" in error.value for error in app.error))
        app.text_area(key="aptamer_input").set_value("").run()
        self.assertEqual(app.session_state["system_analysis_error"], "")
        self.assertFalse(initial_svgs(app))

    def test_missing_nupack_does_not_fabricate_initial_structure(self):
        self.nupack.return_value = None
        app = self.app()
        self.fill_sequences(app)
        self.analyze.assert_not_called()
        self.assertIsNone(app.session_state["system_analysis_result"])
        self.assertFalse(initial_svgs(app))
        self.assertTrue(any("NUPACK est requis pour la structure 2D" in info.value for info in app.info))
        self.assertTrue(button(app, "Étendre le trigger").disabled)
        self.nupack.return_value = "/fixture/nupack"
        app.run()
        self.assertTrue(initial_svgs(app))
        self.nupack.return_value = None
        app.run()
        self.assertIsNone(app.session_state["system_analysis_result"])
        self.assertEqual(app.session_state["aptamer_indices"], [])
        self.assertFalse(initial_svgs(app))

    def test_failed_analysis_can_be_retried_without_showing_a_fake_structure(self):
        self.analyze.side_effect = RuntimeError("NUPACK indisponible")
        app = self.app()
        self.fill_sequences(app)
        self.assertFalse(initial_svgs(app))
        self.assertFalse(any(item.label.startswith("Structure initiale") for item in app.get("download_button")))
        app.checkbox(key="extension_ligand_base_3").check().run()
        self.analyze.assert_called_once()
        self.analyze.side_effect = self.analysis_result
        button(app, "Réessayer l’analyse du système").click().run()
        self.assertFalse(app.exception)
        self.assertEqual(self.analyze.call_count, 2)
        self.assertEqual(len(initial_svgs(app)), 1)
        self.assertEqual(app.session_state["system_analysis_error"], "")

    def test_extension_request_reuses_initial_system_conditions_and_annotations(self):
        app = self.app()
        self.fill_sequences(app)
        app.selectbox(key="extension_material").set_value("ARN · rna06")
        app.number_input(key="extension_temperature").set_value(26.75)
        app.number_input(key="extension_sodium").set_value(42.0)
        app.number_input(key="extension_magnesium").set_value(7.0)
        app.checkbox(key="extension_ligand_base_3").check().run()
        initial_thermo = self.analyze.call_args.args[3]
        with patch("aptaswitch_studio.extension_runtime.ExtensionRunJob.start") as start:
            button(app, "Étendre le trigger").click().run()
        self.assertFalse(app.exception)
        start.assert_called_once()
        request = app.session_state["extension_job"].request
        self.assertEqual(request.thermo, initial_thermo)
        self.assertEqual(request.ligand_aptamer_positions, (3,))
        self.assertEqual(request.aptamer_dna, APTAMER)
        self.assertEqual(request.trigger_dna, TRIGGER)
        self.assertEqual(app.radio(key="nav").value, "Résultats")


if __name__ == "__main__":
    unittest.main()
