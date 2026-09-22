"""The app must not present missing inputs or demonstration data as a run."""

from __future__ import annotations

import tempfile
import unittest
from contextlib import ExitStack
from dataclasses import replace
from pathlib import Path
from unittest.mock import patch

from tests.ui_helpers import FrenchAppTest as AppTest

from aptaswitch_core.architecture import default_architecture_for_trigger
from aptaswitch_core.models import RunRequest, SequenceInput
from tests.design_fixtures import design_result_fixture


APP_PATH = Path(__file__).resolve().parents[1] / "src" / "aptaswitch_studio" / "web_app.py"
APTAMER = "ACGTACGT"
TRIGGER = "GTCCAGGCTGGTATAATTAGATCCACGTAC"
PREVIEW_TITLE = "4 · Aperçus du switch en direct"


class CallbackState(dict):
    """Support the attribute and mapping APIs used by Streamlit callbacks."""

    __getattr__ = dict.__getitem__
    __setattr__ = dict.__setitem__


class RealInputsOnlyTests(unittest.TestCase):
    def setUp(self):
        self.context = ExitStack()
        self.addCleanup(self.context.close)
        self.output = Path(self.context.enter_context(tempfile.TemporaryDirectory()))
        for name, value in (("load_settings", {}), ("detect_nupack", None)):
            self.context.enter_context(patch(f"aptaswitch_core.nupack_setup.{name}", return_value=value))
        for name in ("save_settings", "save_nupack_settings"):
            self.context.enter_context(patch(
                f"aptaswitch_core.nupack_setup.{name}",
                side_effect=AssertionError("Tests must not modify user settings"),
            ))
        self.nupack = self.context.enter_context(patch(
            "aptaswitch_core.nupack_setup.configured_nupack_import_path", return_value="/test/nupack",
        ))
        # These are UI tests. Do not perform automatic system folding or run a
        # scientific engine merely because both sequence inputs are present.
        self.context.enter_context(patch(
            "aptaswitch_studio.extension_ui.configured_nupack_import_path", return_value=None,
        ))
        self.start = self.context.enter_context(patch(
            "aptaswitch_studio.web_runtime.DesignRunJob.start", autospec=True,
        ))

    def app(self, aptamer="", trigger="", **session):
        app = AppTest.from_file(str(APP_PATH), default_timeout=20)
        for key, value in {
            "dependency_detection_done": True,
            "output_dir": str(self.output),
            "aptamer_input": aptamer,
            "trigger_input": trigger,
            **session,
        }.items():
            app.session_state[key] = value
        app.run()
        self.assertFalse(app.exception)
        return app

    @staticmethod
    def run_button(app):
        return next(button for button in app.button if button.label == "Lancer la conception")

    def assert_no_switch_preview(self, app):
        self.assertNotIn(PREVIEW_TITLE, [heading.value for heading in app.subheader])
        self.assertNotIn("preview_state", {radio.key for radio in app.radio})
        self.assertFalse(any(button.label == "Télécharger cet aperçu PNG"
                             for button in app.get("download_button")))

    @staticmethod
    def visible_words(app):
        values = [item.value for items in (app.markdown, app.caption, app.subheader, app.info, app.warning)
                  for item in items]
        values.extend(item.label for items in (app.button, app.checkbox, app.selectbox, app.number_input, app.tabs)
                      for item in items)
        return " ".join(values).lower()

    def test_absent_invalid_or_short_sequences_have_no_switch_preview_or_run(self):
        for aptamer, trigger in (
            ("", ""),
            (APTAMER, ""),
            ("", TRIGGER),
            ("ACGN", TRIGGER),
            (APTAMER, "ACGN"),
            (APTAMER, TRIGGER[:12]),
        ):
            with self.subTest(aptamer=aptamer, trigger=trigger):
                app = self.app(aptamer, trigger)
                self.assert_no_switch_preview(app)
                self.assertTrue(self.run_button(app).disabled)
                self.assertIsNone(app.session_state["design_job"])
        self.start.assert_not_called()

    def test_valid_inputs_use_their_real_trigger_and_only_start_nupack(self):
        # An existing browser session may still carry the removed mock choice.
        app = self.app(APTAMER, TRIGGER, engine="mock")
        self.assertIn(PREVIEW_TITLE, [heading.value for heading in app.subheader])
        self.assertFalse(self.run_button(app).disabled)
        self.assertEqual(app.session_state["engine"], "nupack")
        self.assertFalse(any("mock" in option.lower() for selector in app.selectbox for option in selector.options))
        captions = " ".join(item.value.lower() for item in app.caption)
        self.assertIn("cible", captions)
        self.assertIn("nupack", captions)

        self.run_button(app).click().run()
        self.assertFalse(app.exception)
        self.start.assert_called_once()
        request = self.start.call_args.args[0].request
        self.assertEqual(request.engine, "nupack")
        self.assertEqual(request.nupack_path, "/test/nupack")
        self.assertEqual(request.sequences.aptamer_dna, APTAMER)
        self.assertEqual(request.sequences.effective_trigger(), TRIGGER)
        self.assertEqual((request.architecture.toehold_length, request.architecture.stem_length), (12, 18))

    def test_clearing_an_input_removes_the_previous_preview(self):
        app = self.app(APTAMER, TRIGGER)
        self.assertIn(PREVIEW_TITLE, [heading.value for heading in app.subheader])
        app.text_area(key="aptamer_input").set_value("").run()
        self.assertFalse(app.exception)
        self.assert_no_switch_preview(app)
        self.assertTrue(self.run_button(app).disabled)

    def test_short_trigger_has_preview_only_after_explicit_customization(self):
        app = self.app(APTAMER, TRIGGER[:12])
        self.assert_no_switch_preview(app)
        self.assertTrue(self.run_button(app).disabled)
        app.checkbox(key="customize_tswitch").check().run()
        self.assertFalse(app.exception)
        self.assertIn(PREVIEW_TITLE, [heading.value for heading in app.subheader])
        self.assertFalse(self.run_button(app).disabled)

        app.text_area(key="trigger_input").set_value("A").run()
        self.assertFalse(app.exception)
        self.assert_no_switch_preview(app)
        self.assertTrue(self.run_button(app).disabled)

    def test_missing_nupack_cannot_fall_back_to_a_demonstration_run(self):
        self.nupack.return_value = None
        app = self.app(APTAMER, TRIGGER, engine="mock")
        self.assertTrue(self.run_button(app).disabled)
        self.assertEqual(app.session_state["engine"], "nupack")
        self.assertNotIn("mock", self.visible_words(app))
        self.assertNotIn("essai à blanc", self.visible_words(app))
        # Exercise the callback too: disabling the button must not be the only
        # protection against an old or manually submitted widget event.
        from aptaswitch_studio import web_app
        state = CallbackState(app.session_state.filtered_state)
        with patch.object(web_app.st, "session_state", state), patch.object(
            web_app, "configured_nupack_import_path", return_value=None,
        ):
            web_app._start_run()
        self.start.assert_not_called()
        self.assertIsNone(app.session_state["design_job"])
        self.assertIsNone(state.design_job)
        self.assertIn("NUPACK", state.run_error)

    def result(self, engine):
        request = RunRequest(
            sequences=SequenceInput(APTAMER, TRIGGER),
            architecture=default_architecture_for_trigger(len(TRIGGER)),
            trials=1,
            output_dir=self.output,
            engine="nupack",
        )
        result = design_result_fixture(request)
        # Reconstruct legacy session metadata without executing a removed
        # engine or passing it through today's request validation.
        result.request = replace(result.request, engine=engine)
        return result

    def test_stale_mock_result_is_not_presented_as_scientific_output(self):
        app = self.app(nav="Résultats", last_result=self.result("mock"))
        self.assertFalse(app.get("table"))
        self.assertNotIn("selected_trial", {selector.key for selector in app.selectbox})
        self.assertFalse(app.get("download_button"))
        self.assertNotIn("MOCK", [metric.value for metric in app.metric])

    def test_results_and_settings_have_no_kinetics_controls(self):
        app = self.app(nav="Résultats", last_result=self.result("nupack"))
        self.assertEqual(len(app.get("table")), 1)
        self.assertIn("selected_trial", {selector.key for selector in app.selectbox})
        for page in ("Résultats", "Réglages"):
            if page != "Résultats":
                app.radio(key="nav").set_value(page).run()
            self.assertFalse(app.exception)
            words = self.visible_words(app)
            self.assertNotIn("cinétique", words)
            self.assertNotIn("multistrand", words)
            self.assertFalse(any((item.key or "").startswith(("kinetic_", "ms_"))
                                 for items in (app.number_input, app.checkbox, app.button, app.text_input)
                                 for item in items))


if __name__ == "__main__":
    unittest.main()
