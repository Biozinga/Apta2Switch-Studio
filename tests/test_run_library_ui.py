"""Saved-run navigation must be read-only and independent of NUPACK."""

from __future__ import annotations

import json
import shutil
import sys
import tempfile
import unittest
from contextlib import ExitStack
from dataclasses import replace
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import Mock, patch

from tests.ui_helpers import FrenchAppTest as AppTest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from aptaswitch_core.architecture import default_architecture_for_trigger
from aptaswitch_core.exports import export_run
from aptaswitch_core.models import RunRequest, SequenceInput
from aptaswitch_studio.web_runtime import JobSnapshot
from tests.design_fixtures import design_result_fixture


APP_PATH = ROOT / "src" / "aptaswitch_studio" / "web_app.py"
TRIGGER = "GTCCAGGCTGGTATAATTAGATCC"


class RunLibraryUITests(unittest.TestCase):
    def setUp(self):
        self.context = ExitStack()
        self.addCleanup(self.context.close)
        self.workspace = Path(self.context.enter_context(tempfile.TemporaryDirectory())).resolve()
        self.runs = self.workspace / "Runs"
        self.legacy = self.workspace / "Legacy runs"
        self.custom = self.workspace / "Custom output"
        for path in (self.runs, self.legacy, self.custom):
            path.mkdir()
        for name, value in (
            ("load_settings", {}),
            ("detect_nupack", None),
            ("configured_nupack_import_path", None),
        ):
            self.context.enter_context(patch(f"aptaswitch_core.nupack_setup.{name}", return_value=value))
        for name in ("save_settings", "save_nupack_settings"):
            self.context.enter_context(patch(
                f"aptaswitch_core.nupack_setup.{name}",
                side_effect=AssertionError("Opening a run must not write user settings"),
            ))
        self.context.enter_context(patch("aptaswitch_studio.run_library_ui.default_export_dir", return_value=self.runs))
        self.context.enter_context(patch("aptaswitch_studio.run_library_ui.legacy_export_dir", return_value=self.legacy))
        self.design_start = self.context.enter_context(patch(
            "aptaswitch_studio.web_runtime.DesignRunJob.start",
            side_effect=AssertionError("Opening a saved run must not start a design"),
        ))
        self.extension_start = self.context.enter_context(patch(
            "aptaswitch_studio.extension_runtime.ExtensionRunJob.start",
            side_effect=AssertionError("Opening a saved run must not start an extension"),
        ))

    def saved_design(self, name="saved_design", *, root=None):
        output_root = root or self.runs / "designs"
        request = RunRequest(
            sequences=SequenceInput("ACGTACGT", TRIGGER, molecule_name="Saved system"),
            architecture=default_architecture_for_trigger(len(TRIGGER)),
            trials=2,
            output_dir=output_root,
        )
        result = replace(design_result_fixture(request), run_id=name, output_dir=output_root / name)
        exports = export_run(result)
        return result, exports

    def saved_extension(self):
        result = json.loads((ROOT / "tests" / "fixtures" / "extension_nupack_smoke.json").read_text())
        folder = self.runs / "extensions" / result["run_id"]
        folder.mkdir(parents=True)
        result["output_dir"] = str(folder)
        path = folder / f"{result['run_id']}_results.json"
        path.write_text(json.dumps(result), encoding="utf-8")
        (folder / f"{result['run_id']}.log").write_text("Saved NUPACK extension run\n", encoding="utf-8")
        return result, path

    def app(self):
        app = AppTest.from_file(str(APP_PATH), default_timeout=30)
        app.session_state["dependency_detection_done"] = True
        app.session_state["output_dir"] = str(self.custom)
        app.run()
        self.assertFalse(app.exception)
        # Keep a different form ready while inspecting historical results.
        app.text_area(key="aptamer_input").set_value("TGCAACGT")
        app.text_area(key="trigger_input").set_value(TRIGGER)
        app.text_input(key="molecule_name").set_value("Current form, not the saved system")
        app.number_input(key="trigger_concentration_value").set_value(8.5)
        app.run()
        app.radio(key="nav").set_value("Résultats").run()
        self.assertFalse(app.exception)
        return app

    def assert_form_preserved(self, app):
        self.assertEqual(app.session_state["aptamer_input"], "TGCAACGT")
        self.assertEqual(app.session_state["trigger_input"], TRIGGER)
        self.assertEqual(app.session_state["molecule_name"], "Current form, not the saved system")
        self.assertEqual(app.session_state["trigger_concentration_value"], 8.5)
        self.assertEqual(app.session_state["output_dir"], str(self.custom))

    @staticmethod
    def job(state="running", result=None):
        return SimpleNamespace(
            request=SimpleNamespace(engine="nupack"),
            snapshot=Mock(return_value=JobSnapshot(
                state=state, current=0, total=2, message="NUPACK en cours",
                result=result,
            )),
            consume_once=Mock(return_value=False),
            cancel=Mock(),
        )

    def test_library_finds_all_roots_and_opens_design_without_recomputing(self):
        expected, exports = self.saved_design()
        _, custom_exports = self.saved_design("custom_design", root=self.custom)
        _, legacy_exports = self.saved_design("legacy_design", root=self.legacy)
        _, extension_path = self.saved_extension()
        app = self.app()

        selector = app.selectbox(key="saved_run_path")
        for path in (exports["json"], custom_exports["json"], legacy_exports["json"], str(extension_path)):
            selector.set_value(path).run()
            self.assertFalse(app.exception)
            self.assertEqual(app.selectbox(key="saved_run_path").value, path)
            selector = app.selectbox(key="saved_run_path")

        app.selectbox(key="saved_run_path").set_value(exports["json"]).run()
        completed = self.job("done", expected)
        app.session_state["design_job"] = completed
        app.session_state["candidate_state"] = "ON + trigger"
        app.session_state["selected_trial"] = 999
        app.button(key="open_saved_run").click().run()

        self.assertFalse(app.exception)
        loaded = app.session_state["last_result"]
        self.assertEqual(loaded.run_id, expected.run_id)
        self.assertEqual(loaded.candidates, expected.candidates)
        self.assertEqual(app.radio(key="results_kind").value, "Design de switches")
        self.assertEqual(len(app.get("table")), 1)
        self.assertEqual(app.radio(key="candidate_state").value, "OFF")
        self.assertIsNone(app.session_state["design_job"])
        completed.consume_once.assert_not_called()
        self.assertIn("Télécharger la structure PNG", [button.label for button in app.get("download_button")])
        self.assert_form_preserved(app)
        self.design_start.assert_not_called()
        self.extension_start.assert_not_called()

    def test_open_moved_design_folder_rebases_available_exports(self):
        expected, _ = self.saved_design()
        moved = self.workspace / "Imported folder" / expected.run_id
        moved.parent.mkdir()
        shutil.move(expected.output_dir, moved)
        app = self.app()

        app.text_input(key="run_open_path").set_value(str(moved))
        app.button(key="open_run_path").click().run()

        self.assertFalse(app.exception)
        loaded = app.session_state["last_result"]
        self.assertEqual(loaded.run_id, expected.run_id)
        self.assertEqual(loaded.output_dir.resolve(), moved.resolve())
        exports = app.session_state["last_exports"]
        self.assertIn("json", exports)
        self.assertIn("csv", exports)
        self.assertTrue(all(Path(path).is_file() and Path(path).parent == moved for path in exports.values()))
        self.assert_form_preserved(app)

    def test_open_extension_retains_design_and_survives_page_navigation(self):
        previous, _ = self.saved_design()
        expected, path = self.saved_extension()
        app = self.app()
        app.session_state["last_result"] = previous
        completed = self.job("done", expected)
        app.session_state["extension_job"] = completed
        app.session_state["extension_selected_5"] = 999
        app.session_state["extension_selected_6"] = 999
        app.text_input(key="run_open_path").set_value(str(path))
        app.button(key="open_run_path").click().run()

        self.assertFalse(app.exception)
        self.assertEqual(app.radio(key="results_kind").value, "Extension de triggers")
        self.assertEqual(app.session_state["extension_result"]["candidates"], expected["candidates"])
        self.assertEqual(app.session_state["last_result"].run_id, previous.run_id)
        self.assertIsNone(app.session_state["extension_job"])
        completed.consume_once.assert_not_called()
        self.assertEqual(app.selectbox(key="extension_selected_5").value, 0)
        self.assertEqual(app.selectbox(key="extension_selected_6").value, 0)
        self.assertIn("Extensions criblées", [metric.label for metric in app.metric])
        self.assertTrue(any(button.label == "Structure PNG" for button in app.get("download_button")))
        self.assert_form_preserved(app)

        app.radio(key="nav").set_value("Conception").run()
        self.assert_form_preserved(app)
        app.radio(key="nav").set_value("Résultats").run()
        self.assertFalse(app.exception)
        self.assertEqual(app.session_state["extension_result"]["run_id"], expected["run_id"])
        self.assertEqual(app.radio(key="results_kind").value, "Extension de triggers")
        app.radio(key="results_kind").set_value("Design de switches").run()
        self.assertFalse(app.exception)
        self.assertEqual(app.session_state["last_result"].run_id, previous.run_id)
        self.assertEqual(len(app.get("table")), 1)
        self.design_start.assert_not_called()
        self.extension_start.assert_not_called()

    def test_invalid_path_or_json_keeps_previous_results(self):
        previous, exports = self.saved_design()
        app = self.app()
        app.session_state["last_result"] = previous
        app.session_state["last_exports"] = exports
        invalid_json = self.workspace / "not_a_run.json"
        invalid_json.write_text('{"message": "not a saved run"}', encoding="utf-8")

        for path in (self.workspace / "missing run", invalid_json):
            with self.subTest(path=path):
                app.text_input(key="run_open_path").set_value(str(path))
                app.button(key="open_run_path").click().run()
                self.assertFalse(app.exception)
                self.assertTrue(app.error)
                self.assertEqual(app.session_state["last_result"].run_id, previous.run_id)
                self.assertEqual(app.session_state["last_exports"], exports)
                self.assertEqual(app.radio(key="results_kind").value, "Design de switches")

    def test_job_started_after_render_prevents_opening_without_cancelling(self):
        _, exports = self.saved_design()
        previous, _ = self.saved_design("currently_displayed")
        for job_key in ("design_job", "extension_job"):
            with self.subTest(job=job_key):
                app = self.app()
                app.session_state["last_result"] = previous
                app.selectbox(key="saved_run_path").set_value(exports["json"]).run()
                # The page was rendered before another calculation started;
                # the callback must still inspect the current job state.
                running = self.job()
                app.session_state[job_key] = running
                app.button(key="open_saved_run").click().run()

                self.assertFalse(app.exception)
                self.assertEqual(app.session_state["last_result"].run_id, previous.run_id)
                self.assertIs(app.session_state[job_key], running)
                running.cancel.assert_not_called()
                self.assertTrue(app.button(key="open_saved_run").disabled)
                self.assertTrue(app.button(key="open_run_path").disabled)
                self.assertTrue(app.error or app.warning)

    def test_refresh_includes_run_created_after_first_render(self):
        self.saved_design()
        app = self.app()
        _, exports = self.saved_design("newly_completed")
        app.button(key="refresh_run_library").click().run()
        self.assertFalse(app.exception)
        app.selectbox(key="saved_run_path").set_value(exports["json"]).run()
        self.assertEqual(app.selectbox(key="saved_run_path").value, exports["json"])
        app.button(key="open_saved_run").click().run()
        self.assertFalse(app.exception)
        self.assertEqual(app.session_state["last_result"].run_id, "newly_completed")

    def test_old_default_migrates_once_and_custom_output_is_preserved(self):
        # Use the full app so the migration covers a live browser session as
        # well as the values later exposed by the settings page.
        with patch("aptaswitch_studio.web_runtime.default_export_dir", return_value=self.runs), \
             patch("aptaswitch_studio.web_runtime.legacy_export_dir", return_value=self.legacy):
            for old_path, expected_path in ((self.legacy, self.runs), (self.custom, self.custom)):
                with self.subTest(output_dir=old_path):
                    app = AppTest.from_file(str(APP_PATH), default_timeout=30)
                    app.session_state["dependency_detection_done"] = True
                    app.session_state["output_dir"] = str(old_path)
                    app.session_state["settings_default_output"] = str(old_path)
                    app.run()
                    self.assertFalse(app.exception)
                    self.assertEqual(app.session_state["output_dir"], str(expected_path))
                    self.assertEqual(app.session_state["settings_default_output"], str(expected_path))

            # A user may deliberately choose the old location after migration.
            # Subsequent reruns must not force it back to the new default.
            app.session_state["output_dir"] = str(self.legacy)
            app.run()
            self.assertFalse(app.exception)
            self.assertEqual(app.session_state["output_dir"], str(self.legacy))


if __name__ == "__main__":
    unittest.main()
