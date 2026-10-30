"""Control routing with explicit UI doubles; no scientific jobs are launched."""

from dataclasses import replace
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

from streamlit.testing.v1 import AppTest as LocalizedAppTest

from aptaswitch_core.models import RunRequest, SequenceInput
from aptaswitch_core.architecture import default_architecture_for_trigger
from aptaswitch_core.design_checkpoint import save_design_checkpoint
from aptaswitch_studio.run_library import load_run
from aptaswitch_studio.web_runtime import JobSnapshot
from tests.design_fixtures import design_result_fixture
from tests.ui_helpers import FrenchAppTest as AppTest, candidate_tables


ROOT = Path(__file__).resolve().parents[1]


class RecordingDesignJob:
    def __init__(self, *, continuation=True, state="running", result=None, phase="design"):
        self.supports_continuation = continuation
        self.calls = []
        self.consumed = False
        self.request = result.request if result else SimpleNamespace(trials=3, engine="native-hybrid")
        self.value = JobSnapshot(
            state=state, current=1, total=3, message="Essai terminé", result=result,
            phase=phase, can_request_control=state == "running" and phase == "design",
            exports={"checkpoint": str(result.output_dir / "design_checkpoint.json")} if result else {},
        )

    def snapshot(self):
        return self.value

    def pause(self):
        self.calls.append("pause")
        self.value = replace(self.value, control_requested="pause")

    def finish_early(self):
        self.calls.append("finish")
        self.value = replace(self.value, control_requested="finish")

    def cancel(self):
        self.calls.append("cancel")

    def consume_once(self):
        if self.consumed or self.value.state == "running":
            return False
        self.consumed = True
        return True


CONTROL_APP = """
import streamlit as st
from aptaswitch_studio.design_run_ui import render_design_run_controls
job = st.session_state.design_job
render_design_run_controls(job, job.snapshot())
"""

CONTINUATION_APP = """
import streamlit as st
from aptaswitch_studio.design_run_ui import render_design_continuation
render_design_continuation(st.session_state.last_result)
if st.session_state.get('run_error'):
    st.error(st.session_state.run_error)
"""


class DesignRunControlsUITests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix="aptaswitch-controls-ui-")
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)

    def result(self, *, status="paused", phase="design"):
        request = RunRequest(sequences=SequenceInput(
            aptamer_dna="ACGTACGT", trigger_dna="GTCCAGGCTGGTATAATTAGATCC"),
            architecture=default_architecture_for_trigger(24), trials=3, output_dir=self.root)
        result = design_result_fixture(request)
        result.status = status
        result.candidates = result.candidates[:1]
        result.optimization["trials"] = result.optimization["trials"][:1]
        result.optimization["checkpoint"] = {"phase": phase, "completed_trials": 1}
        result.output_dir.mkdir(parents=True)
        # Existence fixture only: resume API is mocked to test UI routing.
        (result.output_dir / "design_checkpoint.json").write_text("{}")
        return result

    def controls(self, job):
        app = AppTest.from_string(CONTROL_APP)
        app.session_state.design_job = job
        return app.run()

    def test_pause_and_finish_are_distinct_and_pending_request_remains_visible(self):
        for key, action, message in (("pause_design_run", "pause", "Pause demandée"),
                                     ("finish_design_run", "finish", "Fin demandée")):
            with self.subTest(action=action):
                job = RecordingDesignJob()
                app = self.controls(job)
                self.assertEqual({button.label for button in app.button},
                                 {"Terminer maintenant", "Mettre en pause"})
                app.button(key=key).click().run()
                self.assertEqual(job.calls, [action])
                self.assertTrue(all(button.disabled for button in app.button))
                self.assertTrue(any(message in item.value for item in app.info))
                app.run()
                self.assertTrue(all(button.disabled for button in app.button))
                self.assertTrue(any(message in item.value for item in app.info))
                self.assertFalse(app.exception)

    def test_independent_inverse_reanalysis_keeps_cancel_only(self):
        job = RecordingDesignJob(continuation=False)
        app = self.controls(job)
        self.assertEqual([button.label for button in app.button], ["Annuler le calcul"])
        app.button[0].click().run()
        self.assertEqual(job.calls, ["cancel"])
        self.assertFalse(app.exception)

    def test_finalizing_does_not_offer_cancel_or_claim_pause_is_available(self):
        app = self.controls(RecordingDesignJob(phase="inverse_extension"))
        self.assertTrue(all(button.disabled for button in app.button))
        self.assertNotIn("Annuler le calcul", [button.label for button in app.button])
        self.assertTrue(any("Finalisation en cours" in info.value for info in app.info))
        self.assertFalse(app.exception)

    def test_resume_and_finish_use_checkpoint_and_not_changed_form_settings(self):
        result = self.result()
        for key, finish in (("resume_design_run", False), ("finish_saved_design_run", True)):
            with self.subTest(finish=finish):
                app = AppTest.from_string(CONTINUATION_APP)
                app.session_state.last_result = result
                app.session_state.trigger_input = "AAAA"
                app.session_state.trials = 99
                app.session_state.temperature_c = 55.0
                app.session_state.last_exports = {"csv": "old.csv"}
                app.run()
                resumed = RecordingDesignJob(result=result)
                with patch("aptaswitch_studio.design_run_ui.active_job_running", return_value=False), \
                     patch("aptaswitch_studio.design_run_ui.configured_library_path", return_value="/new/library"), \
                     patch("aptaswitch_studio.web_runtime.resume_design_job", return_value=resumed) as resume:
                    app.button(key=key).click().run()
                resume.assert_called_once_with(result.output_dir / "design_checkpoint.json",
                                               finish_early=finish, library_path="/new/library")
                self.assertIs(app.session_state.design_job, resumed)
                self.assertEqual(app.session_state.last_exports, {})
                self.assertEqual(resumed.request.trials, 3)
                self.assertFalse(any(button.key == "resume_design_run" for button in app.button))
                self.assertFalse(app.exception)

    def test_completed_trials_with_exports_pending_still_offer_continuation(self):
        app = AppTest.from_string(CONTINUATION_APP)
        app.session_state.last_result = self.result(status="completed", phase="exports")
        app.run()
        self.assertIn("Reprendre le calcul", [button.label for button in app.button])
        self.assertTrue(any("interrompu" in item.value for item in app.markdown))
        self.assertFalse(app.exception)

    def test_missing_checkpoint_disables_both_actions(self):
        result = self.result()
        (result.output_dir / "design_checkpoint.json").unlink()
        app = AppTest.from_string(CONTINUATION_APP)
        app.session_state.last_result = result
        app.run()
        self.assertTrue(all(button.disabled for button in app.button))
        self.assertTrue(any("introuvable" in item.value for item in app.error))
        self.assertFalse(app.exception)

    def test_done_checkpoint_has_no_resume_and_stopped_count_is_explicit(self):
        app = AppTest.from_string(CONTINUATION_APP)
        app.session_state.last_result = self.result(status="stopped", phase="done")
        app.run()
        self.assertFalse(app.button)
        self.assertTrue(any("1/3 essais terminés" in item.value for item in app.info))
        self.assertFalse(app.exception)

    def test_pause_monitor_keeps_partial_result_without_failure_or_final_progress(self):
        result = self.result()
        job = RecordingDesignJob(state="paused", result=result)
        app = AppTest.from_string("""
import streamlit as st
from aptaswitch_studio.web_app import _design_job_monitor
from aptaswitch_studio.design_run_ui import render_design_continuation
_design_job_monitor()
if st.session_state.get('last_result'):
    render_design_continuation(st.session_state.last_result)
""")
        app.session_state.design_job = job
        app.session_state.run_error = ""
        with patch("aptaswitch_studio.web_app.persist_active_design", return_value=True) as persist:
            app.run()
        persist.assert_called_once()
        self.assertIs(app.session_state.last_result, result)
        self.assertEqual(app.session_state.run_error, "")
        self.assertFalse(app.error)
        self.assertFalse(app.get("progress"))
        self.assertIn("Reprendre le calcul", [button.label for button in app.button])
        self.assertFalse(app.exception)

    def test_results_page_displays_paused_candidates_but_not_final_checkpoint_export(self):
        result = self.result()
        app = AppTest.from_file(str(ROOT / "src/aptaswitch_studio/web_app.py"), default_timeout=20).run()
        app.session_state.design_job = RecordingDesignJob(state="paused", result=result)
        app.radio(key="nav").set_value("Résultats").run()
        self.assertFalse(app.exception)
        self.assertFalse(app.error)
        self.assertEqual(len(candidate_tables(app)), 1)
        self.assertEqual(len(candidate_tables(app)[0].value), 1)
        self.assertIn("Reprendre le calcul", [button.label for button in app.button])
        self.assertFalse(app.get("progress"))
        self.assertFalse(any(button.label in {"CHECKPOINT", "Tous les exports (.zip)"}
                             for button in app.get("download_button")))

    def test_disk_checkpoint_hides_json_and_stale_csv_final_exports(self):
        app = AppTest.from_file(str(ROOT / "src/aptaswitch_studio/web_app.py"), default_timeout=20).run()
        result = self.result()
        output = Path(app.session_state.output_dir)
        result.request = replace(result.request, output_dir=output)
        result.output_dir = output / "designs" / result.run_id
        checkpoint = save_design_checkpoint(result, phase="design")
        csv = result.output_dir / f"{result.run_id}_ranked.csv"
        csv.write_text("outdated export from a previous calculation boundary\n")
        loaded = load_run(checkpoint)
        self.assertEqual(loaded.exports["json"], str(checkpoint))
        self.assertEqual(loaded.exports["csv"], str(csv))
        before = {path: path.read_bytes() for path in (checkpoint, csv)}

        app.radio(key="nav").set_value("Résultats").run()
        app.selectbox(key="saved_run_path").set_value(str(checkpoint)).run()
        app.button(key="open_saved_run").click().run()

        self.assertFalse(app.exception)
        self.assertFalse(app.error)
        self.assertEqual(app.session_state.last_exports["json"], str(checkpoint))
        self.assertEqual(app.session_state.last_exports["csv"], str(csv))
        self.assertIn("Reprendre le calcul", [button.label for button in app.button])
        self.assertFalse(any((button.key or "").startswith("switch_run_")
                             for button in app.get("download_button")))
        self.assertNotIn("Exporter les résultats", [item.proto.popover.label for item in app.get("popover")])
        self.assertEqual({path: path.read_bytes() for path in before}, before)

    def test_control_and_resume_labels_follow_real_french_and_english_catalogs(self):
        result = self.result()
        language_setup = """
import streamlit as st
from aptaswitch_core.localization import set_language
set_language(st.session_state.ui_language)
"""
        for language, running_labels, continuation_labels, pending in (
            ("fr", {"Terminer maintenant", "Mettre en pause"},
             {"Reprendre le calcul", "Terminer avec les résultats actuels"}, "Pause demandée"),
            ("en", {"Finish now", "Pause"},
             {"Resume calculation", "Finish with current results"}, "Pause requested"),
        ):
            with self.subTest(language=language):
                app = LocalizedAppTest.from_string(language_setup + CONTROL_APP)
                app.session_state.ui_language = language
                app.session_state.design_job = RecordingDesignJob()
                app.run()
                self.assertEqual({button.label for button in app.button}, running_labels)
                app.button(key="pause_design_run").click().run()
                self.assertTrue(any(pending in item.value for item in app.info))
                self.assertFalse(app.exception)

                app = LocalizedAppTest.from_string(language_setup + CONTINUATION_APP)
                app.session_state.ui_language = language
                app.session_state.last_result = result
                app.run()
                self.assertEqual({button.label for button in app.button}, continuation_labels)
                self.assertFalse(app.exception)

    def test_monitor_failure_preserves_last_checkpoint_for_recovery(self):
        result = self.result(status="running")
        job = RecordingDesignJob(state="failed", result=result)
        job.value = replace(job.value, error="Échec simulé après sauvegarde")
        app = AppTest.from_string("""
import streamlit as st
from aptaswitch_studio.web_app import _design_job_monitor
from aptaswitch_studio.design_run_ui import render_design_continuation
_design_job_monitor()
if st.session_state.get('last_result'):
    render_design_continuation(st.session_state.last_result)
""")
        app.session_state.design_job = job
        with patch("aptaswitch_studio.web_app.persist_active_design", return_value=True):
            app.run()
        self.assertIs(app.session_state.last_result, result)
        self.assertEqual(app.session_state.run_error, "Échec simulé après sauvegarde")
        self.assertIn("Reprendre le calcul", [button.label for button in app.button])
        self.assertFalse(app.exception)

    def test_running_final_trial_progress_does_not_claim_complete(self):
        job = RecordingDesignJob()
        job.value = replace(job.value, current=3)
        app = AppTest.from_string("""
from aptaswitch_studio.web_app import _design_job_monitor
_design_job_monitor()
""")
        app.session_state.design_job = job
        app.run()
        self.assertEqual(app.get("progress")[0].proto.value, 99)
        self.assertFalse(app.exception)


if __name__ == "__main__":
    unittest.main()
