"""Saved analysis and FFS jobs must not race other native calculation entry points."""

from pathlib import Path
from types import MappingProxyType, SimpleNamespace
import tempfile
import unittest
from unittest.mock import Mock, patch

from aptaswitch_studio import design_session, extension_ui
from aptaswitch_studio import saved_switch_runtime


class SessionState(dict):
    def __getattr__(self, key):
        try:
            return self[key]
        except KeyError as error:
            raise AttributeError(key) from error

    def __setattr__(self, key, value):
        self[key] = value


class SavedSwitchJobGuardTests(unittest.TestCase):
    def setUp(self):
        from tests.native_ui_fixtures import patch_native_runtime
        patch_native_runtime(self)

    def state(self, key):
        job = Mock()
        job.snapshot.return_value = SimpleNamespace(state="running")
        return SessionState(**{**extension_ui.EXTENSION_DEFAULTS, key: job, "ui_language": "fr"})

    def test_saved_workers_participate_in_application_job_guard(self):
        for name in ("saved_switch_analysis_job", "saved_switch_ffs_job"):
            with self.subTest(name=name):
                state = self.state(name)
                self.assertTrue(design_session.active_job_running(state))
                state[name].snapshot.return_value = SimpleNamespace(state="done")
                self.assertFalse(design_session.active_job_running(state))

    def test_initial_tubes_and_extension_start_or_resume_are_blocked(self):
        for name in ("saved_switch_analysis_job", "saved_switch_ffs_job"):
            with self.subTest(name=name):
                state = self.state(name)
                with patch.object(extension_ui.st, "session_state", state), \
                     patch.object(extension_ui, "ExtensionRunJob") as extension, \
                     patch.object(extension_ui, "SystemTubeJob") as tubes:
                    extension_ui._start_system_tubes()
                    extension_ui._start_extension()
                    self.assertTrue(state["extension_error"])
                    extension_ui._resume_extension({})
                extension.assert_not_called()
                tubes.assert_not_called()

    def test_automatic_initial_mfe_is_deferred(self):
        for name in ("saved_switch_analysis_job", "saved_switch_ffs_job"):
            with self.subTest(name=name):
                state = self.state(name)
                with patch.object(extension_ui.st, "session_state", state), \
                     patch.object(extension_ui.st, "markdown"), patch.object(extension_ui.st, "info"), \
                     patch.object(extension_ui, "configured_library_path", return_value="/fake"), \
                     patch.object(extension_ui, "_analyze_initial_system") as compute:
                    extension_ui.render_initial_structure("AACCGGTT", "AACCGG")
                compute.assert_not_called()
                self.assertIsNone(state["system_analysis_signature"])


    def test_browser_reload_finds_running_saved_analysis_in_current_design(self):
        with tempfile.TemporaryDirectory() as directory:
            design = Path(directory).resolve() / "Design"
            candidate = design / "tswitch_enregistres" / "Retained"
            candidate.mkdir(parents=True)
            job = Mock()
            job.saved = SimpleNamespace(path=candidate, name="Retained")
            job.snapshot.return_value = SimpleNamespace(state="running")
            with patch.dict(saved_switch_runtime._JOBS, {str(candidate): job}, clear=True):
                fresh_session = {"active_design_path": str(design)}
                self.assertTrue(design_session.active_job_running(fresh_session))
                self.assertIs(fresh_session["saved_switch_analysis_job"], job)
                self.assertFalse(design_session.active_job_running({"active_design_path": str(design.parent / "Other")}))
                self.assertTrue(design_session.active_job_running(MappingProxyType({"active_design_path": str(design)})))
                with self.assertRaises(ValueError):
                    design_session.reset_design_session({"active_design_path": str(design)})
                job.snapshot.return_value = SimpleNamespace(state="done")
                self.assertFalse(design_session.active_job_running({"active_design_path": str(design)}))

    def test_similarly_named_design_is_not_the_current_project(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory).resolve()
            path = root / "Design copy" / "tswitch_enregistres" / "Retained"
            job = Mock()
            job.saved = SimpleNamespace(path=path)
            job.snapshot.return_value = SimpleNamespace(state="running")
            with patch.dict(saved_switch_runtime._JOBS, {str(path): job}, clear=True):
                self.assertFalse(design_session.active_job_running({"active_design_path": str(root / "Design")}))


if __name__ == "__main__":
    unittest.main()
