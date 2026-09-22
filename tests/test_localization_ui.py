"""Language changes translate the UI without changing scientific inputs or IDs."""

from __future__ import annotations

import copy
import json
import tempfile
import unittest
from contextlib import ExitStack
from pathlib import Path
from unittest.mock import patch

from streamlit.testing.v1 import AppTest

from aptaswitch_core.architecture import default_architecture_for_trigger
from aptaswitch_core.models import RunRequest, SequenceInput
from tests.design_fixtures import design_result_fixture

ROOT = Path(__file__).resolve().parents[1]
APP_PATH = ROOT / "src" / "aptaswitch_studio" / "web_app.py"
TRIGGER = "GTCCAGGCTGGTATAATTAGATCC"


class LocalizationUITests(unittest.TestCase):
    def setUp(self):
        self.context = ExitStack()
        self.addCleanup(self.context.close)
        self.workspace = Path(self.context.enter_context(tempfile.TemporaryDirectory()))
        self.settings = {"default_output_dir": str(self.workspace)}
        self.context.enter_context(patch(
            "aptaswitch_core.nupack_setup.load_settings",
            side_effect=lambda: copy.deepcopy(self.settings),
        ))
        self.save = self.context.enter_context(patch(
            "aptaswitch_core.nupack_setup.save_settings",
            side_effect=lambda settings: self.settings.update(copy.deepcopy(settings)),
        ))
        self.context.enter_context(patch(
            "aptaswitch_core.reporter_library.load_settings",
            side_effect=lambda: copy.deepcopy(self.settings),
        ))
        self.reporter_save = self.context.enter_context(patch(
            "aptaswitch_core.reporter_library.save_settings",
            side_effect=lambda settings: self.settings.update(copy.deepcopy(settings)),
        ))
        self.context.enter_context(patch("aptaswitch_core.nupack_setup.detect_nupack", return_value=None))
        for module in ("aptaswitch_core.nupack_setup", "aptaswitch_studio.extension_ui"):
            self.context.enter_context(patch(
                f"{module}.configured_nupack_import_path", return_value=None,
            ))
        self.design_start = self.context.enter_context(patch(
            "aptaswitch_studio.web_runtime.DesignRunJob.start",
            side_effect=AssertionError("Language changes must not start a design"),
        ))
        self.extension_start = self.context.enter_context(patch(
            "aptaswitch_studio.extension_runtime.ExtensionRunJob.start",
            side_effect=AssertionError("Language changes must not start an extension"),
        ))
        self.context.enter_context(patch("aptaswitch_studio.run_library_ui.default_export_dir", return_value=self.workspace))
        self.context.enter_context(patch("aptaswitch_studio.run_library_ui.legacy_export_dir", return_value=self.workspace))

    def app(self, **state):
        app = AppTest.from_file(str(APP_PATH), default_timeout=30)
        app.session_state["dependency_detection_done"] = True
        app.session_state["output_dir"] = str(self.workspace)
        for key, value in state.items():
            app.session_state[key] = value
        app.run()
        self.assertFalse(app.exception)
        return app

    def open_language_settings(self, app):
        app.radio(key="nav").set_value("Réglages").run()
        self.assertFalse(app.exception)
        # Settings may render inactive tabs eagerly. The key also makes the
        # selected tab observable when opened through the reporter shortcut.
        app.session_state["settings_tab"] = "Général"
        app.run()
        self.assertFalse(app.exception)
        return app.selectbox(key="ui_language")

    def switch_language(self, app, language):
        self.open_language_settings(app).set_value(language).run()
        self.assertFalse(app.exception)
        self.assertEqual(app.session_state["ui_language"], language)

    @staticmethod
    def visible_text(app):
        return "\n".join(str(getattr(item, "value", "")) for kind in (
            "markdown", "caption", "info", "success", "warning", "error", "subheader",
        ) for item in app.get(kind))

    def test_default_is_english_and_removed_badges_do_not_return(self):
        with patch("aptaswitch_core.nupack_setup.configured_nupack_import_path", return_value="/fixture/nupack"):
            app = self.app()
        self.assertEqual(app.session_state["ui_language"], "en")
        self.assertEqual(app.radio(key="nav").options, ["Design", "Results", "Settings"])
        self.assertEqual(app.radio(key="nav").value, "Conception")
        self.assertEqual(app.selectbox(key="extension_material").label, "Thermodynamic model")
        self.assertEqual(app.text_area(key="aptamer_input").label, "Aptamer sequence")
        self.assertEqual(app.text_area(key="trigger_input").label, "Trigger sequence")
        text = self.visible_text(app)
        for removed in (
            "NUPACK prêt", "NUPACK activé", "NUPACK ready", "NUPACK enabled",
            "Studio de conception moléculaire", "Molecular design studio",
            "Application locale · un utilisateur",
        ):
            self.assertNotIn(removed, text)
        self.save.assert_not_called()

    def test_switch_to_french_persists_and_preserves_form_and_navigation_values(self):
        app = self.app()
        app.text_area(key="aptamer_input").set_value("ACGTACGT")
        app.text_area(key="trigger_input").set_value(TRIGGER)
        app.text_input(key="molecule_name").set_value("My aptamer")
        app.number_input(key="trigger_concentration_value").set_value(7.25)
        app.number_input(key="weight_on").set_value(4.25)
        app.run()
        self.assertFalse(app.exception)
        material = app.session_state["material_label"]
        system_material = app.session_state["extension_material"]
        self.switch_language(app, "fr")
        self.assertEqual(self.settings["ui_language"], "fr")
        self.assertEqual(self.settings["default_output_dir"], str(self.workspace))
        self.assertEqual(app.radio(key="nav").options, ["Conception", "Résultats", "Réglages"])
        self.assertEqual(app.radio(key="nav").value, "Réglages")
        self.assertEqual(app.selectbox(key="ui_language").options, ["English", "Français"])
        self.assertEqual(app.selectbox(key="ui_language").label, "Langue")
        app.radio(key="nav").set_value("Conception").run()
        self.assertFalse(app.exception)
        self.assertEqual(app.text_area(key="aptamer_input").value, "ACGTACGT")
        self.assertEqual(app.text_area(key="trigger_input").value, TRIGGER)
        self.assertEqual(app.text_input(key="molecule_name").value, "My aptamer")
        self.assertEqual(app.number_input(key="trigger_concentration_value").value, 7.25)
        self.assertEqual(app.number_input(key="weight_on").value, 4.25)
        self.assertEqual(app.session_state["material_label"], material)
        self.assertEqual(app.session_state["extension_material"], system_material)
        self.assertEqual(app.selectbox(key="extension_material").label, "Modèle thermodynamique")
        self.assertEqual(self.app().session_state["ui_language"], "fr")
        self.design_start.assert_not_called()
        self.extension_start.assert_not_called()

    def test_conditions_and_scoring_help_translate_without_changing_numbers(self):
        app = self.app()
        english_help = app.number_input(key="weight_on").proto.help
        self.assertIn("switch", english_help.lower())
        self.assertNotIn("rendement", english_help.lower())
        self.assertNotIn("concentration initiale", english_help.lower())
        english_stop = app.checkbox(key="exclude_stop_candidates").label
        self.assertIn("Exclude", english_stop)
        self.assertTrue(app.checkbox(key="exclude_stop_candidates").value)
        self.assertTrue(app.checkbox(key="exclude_non_linear_candidates").value)
        self.assertTrue(any("experimental" in item.value.lower() for item in app.caption))
        self.switch_language(app, "fr")
        app.radio(key="nav").set_value("Conception").run()
        self.assertFalse(app.exception)
        french_help = app.number_input(key="weight_on").proto.help
        self.assertNotEqual(english_help, french_help)
        self.assertIn("pourcentage", french_help.lower())
        self.assertEqual(app.checkbox(key="exclude_stop_candidates").label, "Exclure les candidats avec un STOP prématuré")
        self.assertEqual(app.number_input(key="extension_sodium").value, 120.0)
        self.assertEqual(app.selectbox(key="extension_sodium_unit").value, "mM")
        self.assertEqual(app.number_input(key="weight_on").value, 3.0)
        self.switch_language(app, "en")
        app.radio(key="nav").set_value("Conception").run()
        self.assertEqual(app.number_input(key="weight_on").proto.help, english_help)
        self.assertEqual(app.checkbox(key="exclude_stop_candidates").label, english_stop)

    def test_reporter_shortcut_opens_library_without_changing_reporter_or_form(self):
        app = self.app()
        app.text_area(key="aptamer_input").set_value("ACGTACGT")
        app.run()
        app.selectbox(key="reporter_key").set_value("mng").run()
        selector = app.selectbox(key="reporter_key")
        self.assertTrue(any("Add" in option and "reporter" in option.lower() for option in selector.options))
        original = selector.value
        selector.set_value("__add_reporter__").run()
        self.assertFalse(app.exception)
        self.assertEqual(app.session_state["nav"], "Réglages")
        self.assertEqual(app.session_state["settings_tab"], "Reporters")
        self.assertEqual(app.session_state["reporter_key"], original)
        self.assertEqual(app.session_state["aptamer_input"], "ACGTACGT")
        self.assertTrue(any("Add a reporter" in item.value for item in app.markdown))
        self.design_start.assert_not_called()
        self.extension_start.assert_not_called()
        self.save.assert_not_called()

    def test_reporter_shortcut_preserves_selection_from_an_existing_session(self):
        app = self.app(reporter_key="mng")
        self.assertEqual(app.selectbox(key="reporter_key").value, "mng")
        app.selectbox(key="reporter_key").set_value("__add_reporter__").run()
        self.assertFalse(app.exception)
        self.assertEqual(app.session_state["nav"], "Réglages")
        self.assertEqual(app.session_state["settings_tab"], "Reporters")
        self.assertEqual(app.session_state["reporter_key"], "mng")
        self.save.assert_not_called()
        self.reporter_save.assert_not_called()

    def test_reporter_added_through_shortcut_is_selected_on_return_to_design(self):
        app = self.app()
        app.text_area(key="aptamer_input").set_value("ACGTACGT")
        app.text_area(key="trigger_input").set_value(TRIGGER)
        app.number_input(key="weight_on").set_value(4.75)
        app.run()
        app.selectbox(key="reporter_key").set_value("__add_reporter__").run()
        self.assertFalse(app.exception)
        app.text_input(key="new_reporter_label").set_value("My red reporter")
        app.text_area(key="new_reporter_sequence").set_value("acgtacgt")
        next(button for button in app.button if button.label == "Add").click().run()
        self.assertFalse(app.exception)
        self.reporter_save.assert_called_once()
        saved = self.settings["reporters"][0]
        self.assertEqual(saved["label"], "My red reporter")
        self.assertEqual(saved["sequence"], "ACGUACGU")
        self.assertEqual(app.session_state["reporter_key"], saved["key"])
        app.button(key="return_from_reporters").click().run()
        self.assertFalse(app.exception)
        self.assertEqual(app.radio(key="nav").value, "Conception")
        self.assertEqual(app.selectbox(key="reporter_key").value, saved["key"])
        self.assertIn("My red reporter", app.selectbox(key="reporter_key").options)
        self.assertEqual(app.text_area(key="aptamer_input").value, "ACGTACGT")
        self.assertEqual(app.text_area(key="trigger_input").value, TRIGGER)
        self.assertEqual(app.number_input(key="weight_on").value, 4.75)
        self.design_start.assert_not_called()
        self.extension_start.assert_not_called()

    def test_switch_results_table_and_navigation_translate_without_mutating_result(self):
        request = RunRequest(
            sequences=SequenceInput("ACGTACGT", TRIGGER),
            architecture=default_architecture_for_trigger(len(TRIGGER)),
            output_dir=self.workspace, trials=2,
        )
        result = design_result_fixture(request)
        app = self.app(last_result=result, nav="Résultats")
        table = app.get("table")[0].value
        for label in ("Rank", "Status", "First STOP", "Trial"):
            self.assertIn(label, table.columns)
        self.assertIn("no", table["STOP"].tolist())
        self.assertIn("Previous", app.button(key="previous_candidate").label)
        self.assertIn("Next", app.button(key="next_candidate").label)
        app.button(key="next_candidate").click().run()
        self.assertEqual(app.session_state["selected_trial"], 2)
        self.switch_language(app, "fr")
        app.radio(key="nav").set_value("Résultats").run()
        self.assertFalse(app.exception)
        table_fr = app.get("table")[0].value
        for label in ("Rang", "Statut", "Premier STOP", "Essai"):
            self.assertIn(label, table_fr.columns)
        self.assertIn("non", table_fr["STOP"].tolist())
        self.assertEqual(app.session_state["selected_trial"], 2)
        self.assertEqual(app.session_state["last_result"], result)
        self.assertEqual(len(app.get("dataframe")), 0)

    def test_extension_results_translate_table_headers_and_keep_sequences(self):
        result = json.loads((ROOT / "tests" / "fixtures" / "extension_nupack_smoke.json").read_text())
        app = self.app(extension_result=result, nav="Résultats", results_kind="Extension de triggers")
        self.assertEqual(len(app.get("table")), 1)
        table_en = app.get("table")[0].value
        self.assertNotIn("Longueur (nt)", table_en.columns)
        self.assertTrue(any("Length" in label or "length" in label for label in table_en.columns))
        self.switch_language(app, "fr")
        app.radio(key="nav").set_value("Résultats").run()
        self.assertFalse(app.exception)
        table_fr = app.get("table")[0].value
        self.assertIn("Longueur (nt)", table_fr.columns)
        self.assertEqual(table_en.shape, table_fr.shape)
        self.assertEqual(table_en.values.tolist(), table_fr.values.tolist())
        self.assertEqual(app.session_state["extension_result"], result)
        self.assertEqual(len(app.get("dataframe")), 0)


class LocalizationCatalogTests(unittest.TestCase):
    def test_translation_templates_keep_their_named_parameters(self):
        from importlib import import_module
        from string import Formatter

        formatter = Formatter()
        for module in ("en_core", "en_components", "en_main"):
            catalog = import_module(f"aptaswitch_core.locales.{module}").TRANSLATIONS
            self.assertTrue(catalog)
            for source, translated in catalog.items():
                with self.subTest(module=module, source=source):
                    try:
                        source_fields = {field for _, field, _, _ in formatter.parse(source) if field is not None}
                        translated_fields = {field for _, field, _, _ in formatter.parse(translated) if field is not None}
                    except ValueError:
                        # Literal scientific notation is not a format template.
                        continue
                    self.assertEqual(source_fields, translated_fields)

    def test_timestamped_console_lines_translate_while_preserving_the_timestamp(self):
        from aptaswitch_core.localization import language_context, tr

        messages = (
            ("[12:34:56] Préparation de l’extension NUPACK…", "[12:34:56] Preparing NUPACK extension…"),
            ("[12:34:56] Calcul en cours · 12/345", "[12:34:56] Calculation running · 12/345"),
            ("[12:34:56] Les exports seront enregistrés dans /tmp/Conception/échantillon {1}.",
             "[12:34:56] Exports will be saved in /tmp/Conception/échantillon {1}."),
        )
        for french, english in messages:
            with self.subTest(message=french):
                with language_context("en"):
                    self.assertEqual(tr(french), english)
                with language_context("fr"):
                    self.assertEqual(tr(english), french)

    def test_language_context_restores_previous_locale(self):
        from aptaswitch_core.localization import get_language, language_context, tr

        before = get_language()
        with language_context("en"):
            self.assertEqual(tr("Conception"), "Design")
            self.assertEqual(tr("ACGU+UCGC"), "ACGU+UCGC")
            with language_context("fr"):
                self.assertEqual(tr("Conception"), "Conception")
            self.assertEqual(get_language(), "en")
        self.assertEqual(get_language(), before)


if __name__ == "__main__":
    unittest.main()
