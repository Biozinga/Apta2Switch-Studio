"""Persisted context work survives navigation and respects candidate revisions."""

from __future__ import annotations

import copy
import json
from pathlib import Path
import tempfile
import threading
import unittest
from unittest.mock import patch

from aptaswitch_core.models import ReporterContext, RunRequest, SequenceInput, SwitchArchitecture
from aptaswitch_studio import saved_switch_library as library
from aptaswitch_studio import saved_switch_runtime as runtime
from design_fixtures import design_result_fixture
from test_thermodynamics_context_analysis import FixedWeightBackend


class SavedSwitchRuntimeTests(unittest.TestCase):
    def setUp(self):
        self.backend = FixedWeightBackend()
        self.enterContext(patch.object(runtime, "get_backend", return_value=self.backend))
        self.folder = Path(self.enterContext(tempfile.TemporaryDirectory())).resolve()
        self.root = self.folder / "design" / library.LIBRARY_FOLDER
        request = RunRequest(sequences=SequenceInput("AACCGGTTAACC", "AACC"),
                             architecture=SwitchArchitecture(toehold_length=2, stem_length=2),
                             reporter=ReporterContext("test", "Test", "ACG" * 10, "ACG" * 10 + "GCUUAA"),
                             trials=1, output_dir=self.folder / "original")
        self.result = design_result_fixture(request)
        self.candidate = self.result.candidates[0]
        self.saved = self.save()
        self.config = {"temperature_c": 29., "sodium_m": .15, "magnesium_m": 0.,
                       "temperature_values_c": [25., 30.]}
        self.addCleanup(runtime._JOBS.clear)

    def save(self, *, overwrite=False):
        snapshot = library.candidate_snapshot(self.result, self.candidate, "Retained")
        return library.save_switch(self.root, snapshot["name"], snapshot, overwrite=overwrite)

    def test_context_is_detached_and_restored_with_the_assembled_full_sequence(self):
        original = copy.deepcopy(self.saved.snapshot)
        data = runtime.save_context(self.saved, self.config)
        self.config["temperature_c"] = 40
        loaded = runtime.load_context(library.load_switch(self.saved.path))
        self.assertEqual(loaded["config"]["temperature_c"], 29)
        self.assertEqual(loaded["context"]["full_switch_sequence"], self.candidate.switch_seq + "GCUUAA")
        self.assertEqual(self.saved.snapshot, original)
        self.assertEqual(loaded["snapshot_id"], data["snapshot_id"])

    def test_results_in_same_second_are_sorted_by_completion_not_random_folder_suffix(self):
        for suffix, microseconds in (("zzzz", "100000"), ("aaaa", "900000")):
            run_id = "20260930_200000_" + suffix
            folder = self.saved.analyses_dir / run_id
            folder.mkdir(parents=True)
            (folder / "result.json").write_text(json.dumps({
                "snapshot_id": self.saved.snapshot["snapshot_id"], "run_id": run_id,
                "status": "completed", "completed_at": f"2026-09-30T20:00:00.{microseconds}+00:00"}))
        self.assertEqual([row["run_id"] for row in runtime.saved_analyses(self.saved)],
                         ["20260930_200000_aaaa", "20260930_200000_zzzz"])

    def test_atomic_context_failure_retains_previous_valid_file_without_temporary_files(self):
        runtime.save_context(self.saved, self.config)
        before = self.saved.context_path.read_bytes()
        with patch.object(runtime.os, "replace", side_effect=OSError("replace failed")), self.assertRaises(OSError):
            runtime.save_context(self.saved, {**self.config, "temperature_c": 45.})
        self.assertEqual(self.saved.context_path.read_bytes(), before)
        self.assertEqual(list(self.saved.path.glob(".analysis-*")), [])

    def test_stale_revision_cannot_overwrite_context_or_start_native_calculation(self):
        replacement = self.save(overwrite=True)
        with patch.object(runtime, "analyze_saved_switch") as compute, self.assertRaises(ValueError):
            runtime.start_analysis(self.saved, self.config, "/native")
        compute.assert_not_called()
        self.assertFalse(replacement.context_path.exists())

    def test_results_remain_loadable_after_moving_the_complete_candidate_directory(self):
        job = runtime.SavedSwitchAnalysisJob(self.saved, self.config, "/native")
        with patch.object(runtime, "analyze_saved_switch", return_value={"status": "completed", "baseline": {"on_yield_pct": 60}}):
            job._run()
        moved = self.folder / "moved"
        self.saved.path.rename(moved)
        reopened = library.load_switch(moved)
        records = runtime.saved_analyses(reopened)
        self.assertEqual(len(records), 1)
        self.assertTrue(Path(records[0]["result_path"]).is_relative_to(moved))
        self.assertEqual(records[0]["baseline"]["on_yield_pct"], 60)

    def test_cancelled_sweep_keeps_computed_partial_results_and_progress(self):
        job = runtime.SavedSwitchAnalysisJob(self.saved, self.config, "/native")

        def calculate(snapshot, import_path, config, *, progress_callback, cancel_check, backend):
            progress_callback(1, 3, "Reference complete")
            job.cancel()
            self.assertTrue(cancel_check())
            return {"status": "cancelled", "baseline": {"on_yield_pct": 50}, "sweeps": {"temperature_c": []}}

        with patch.object(runtime, "analyze_saved_switch", side_effect=calculate):
            job._run()
        state = job.snapshot()
        self.assertEqual(state.state, "cancelled")
        self.assertEqual(state.current, 1)
        self.assertEqual(state.total, 3)
        self.assertEqual(state.result["baseline"]["on_yield_pct"], 50)
        self.assertEqual(runtime.saved_analyses(self.saved)[0]["status"], "cancelled")
        self.assertTrue((job.folder / "progress.json").is_file())
        self.assertTrue(job.consume_once())
        self.assertFalse(job.consume_once())

    def test_replace_during_work_never_publishes_old_result_into_new_revision(self):
        job = runtime.SavedSwitchAnalysisJob(self.saved, self.config, "/native")

        def calculate(*args, **kwargs):
            self.save(overwrite=True)
            return {"status": "completed", "baseline": {"on_yield_pct": 60}}

        with patch.object(runtime, "analyze_saved_switch", side_effect=calculate):
            job._run()
        self.assertEqual(job.snapshot().state, "failed")
        self.assertIn("replaced", job.snapshot().error)
        self.assertFalse((job.folder / "result.json").exists())
        self.assertEqual(runtime.saved_analyses(library.load_switch(self.saved.path)), [])

    def test_same_candidate_reconnects_active_job_and_cannot_replace_its_config(self):
        started, release = threading.Event(), threading.Event()

        def calculate(*args, **kwargs):
            started.set()
            release.wait(3)
            return {"status": "completed"}

        with patch.object(runtime, "analyze_saved_switch", side_effect=calculate):
            first, launched = runtime.start_analysis(self.saved, self.config, "/native")
            try:
                self.assertTrue(launched)
                self.assertTrue(started.wait(2))
                second, launched_again = runtime.start_analysis(self.saved, {**self.config, "temperature_c": 45}, "/native")
                self.assertIs(first, second)
                self.assertFalse(launched_again)
                self.assertIs(runtime.active_analysis(self.saved), first)
                self.assertEqual(runtime.load_context(self.saved)["config"]["temperature_c"], 29)
            finally:
                release.set()
                first._thread.join(3)
            self.assertFalse(first._thread.is_alive())
            self.assertIsNone(runtime.active_analysis(self.saved))

    def test_finished_job_can_be_restarted_without_erasing_previous_analysis(self):
        with patch.object(runtime, "analyze_saved_switch", return_value={"status": "completed"}):
            first, first_started = runtime.start_analysis(self.saved, self.config, "/native")
            first._thread.join(3)
            second, second_started = runtime.start_analysis(self.saved, self.config, "/native")
            second._thread.join(3)
        self.assertTrue(first_started and second_started)
        self.assertNotEqual(first.run_id, second.run_id)
        self.assertEqual(len(runtime.saved_analyses(self.saved)), 2)

    def test_failed_native_calculation_is_persisted_with_original_configuration(self):
        job = runtime.SavedSwitchAnalysisJob(self.saved, self.config, "/native")
        self.config["sodium_m"] = .9
        with patch.object(runtime, "analyze_saved_switch", side_effect=RuntimeError("native failed")):
            job._run()
        records = runtime.saved_analyses(self.saved)
        self.assertEqual(records[0]["status"], "failed")
        self.assertEqual(records[0]["error"], "native failed")
        self.assertEqual(records[0]["config"]["sodium_m"], .15)

    def test_invalid_context_and_external_symlinks_are_rejected(self):
        runtime.save_context(self.saved, self.config)
        data = json.loads(self.saved.context_path.read_text())
        data["snapshot_id"] = "another-candidate"
        self.saved.context_path.write_text(json.dumps(data))
        with self.assertRaises(ValueError):
            runtime.load_context(self.saved)
        external = self.folder / "external-context.json"
        external.write_text(json.dumps(data))
        self.saved.context_path.unlink()
        self.saved.context_path.symlink_to(external)
        with self.assertRaises(ValueError):
            runtime.load_context(self.saved)

    def test_invalid_analysis_records_do_not_hide_valid_results(self):
        job = runtime.SavedSwitchAnalysisJob(self.saved, self.config, "/native")
        with patch.object(runtime, "analyze_saved_switch", return_value={"status": "completed"}):
            job._run()
        corrupt = self.saved.analyses_dir / "corrupt"
        corrupt.mkdir()
        (corrupt / "result.json").write_text("{")
        stale = self.saved.analyses_dir / "stale"
        stale.mkdir()
        (stale / "result.json").write_text(json.dumps({"snapshot_id": "different", "status": "completed"}))
        self.assertEqual(len(runtime.saved_analyses(self.saved)), 1)


if __name__ == "__main__":
    unittest.main()
