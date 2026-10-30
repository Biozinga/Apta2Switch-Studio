"""Full-context sampling survives navigation without crossing candidate revisions."""
from __future__ import annotations

from copy import deepcopy
from dataclasses import replace
import json
from pathlib import Path
import tempfile
import threading
import unittest
from unittest.mock import patch
import uuid

from aptaswitch_core.structure_ensemble import EnsembleCancelled, ensemble_inputs, ensemble_key, ensemble_path
from aptaswitch_studio import design_session, ensemble_runtime
from aptaswitch_studio import saved_switch_ensemble_runtime as runtime
from aptaswitch_studio import saved_switch_library as library
from test_saved_switch_ensemble import context_fixture
from test_thermodynamics_context_analysis import FixedWeightBackend
from aptaswitch_core.thermodynamics.context_analysis import analyze_saved_switch
from aptaswitch_core.thermodynamics.setup import MODEL_ID


class ContextEnsembleRuntimeTests(unittest.TestCase):
    def setUp(self):
        self.folder = Path(self.enterContext(tempfile.TemporaryDirectory())).resolve()
        self.root = self.folder / "design" / library.LIBRARY_FOLDER
        _, snapshot, _ = context_fixture()
        self.backend = FixedWeightBackend()
        # These arbitrary weights exercise UI worker orchestration. The mocked
        # runtime uses the supported interface identity, never physical results.
        self.backend.provenance = replace(
            self.backend.provenance, backend="AptaSwitch-native-hybrid", model=MODEL_ID,
            version="synthetic-ui-fixture-v1", parameter_digest="fixture:arbitrary-weights",
            ensemble="synthetic fixed-weight orchestration fixture; no physical prediction",
        )
        self.enterContext(patch.object(runtime, "get_backend", return_value=self.backend))
        self.analysis = analyze_saved_switch(self.backend, snapshot, {})
        self.saved = library.save_switch(self.root, snapshot["name"], snapshot)
        self.analysis.update(snapshot_id=snapshot["snapshot_id"], run_id="reference_1")
        self.analysis_folder = self.saved.analyses_dir / self.analysis["run_id"]
        self.analysis_folder.mkdir()
        self.analysis_path = self.analysis_folder / "result.json"
        self.analysis_path.write_text(json.dumps(self.analysis))
        self.source, self.candidate = runtime.load_context_ensemble_source(self.saved, self.analysis)
        inputs = ensemble_inputs(self.source.request, self.candidate, backend=self.backend)
        self.artifact = {"key": ensemble_key(inputs), "inputs": inputs, "marker": "computed families"}
        self.target = ensemble_path(self.analysis_folder, self.artifact["key"])
        self.addCleanup(runtime._JOBS.clear)

    def job(self):
        return runtime.create_context_ensemble_job(
            self.saved, self.analysis, self.source, self.candidate, "/fake",
        )

    def test_unknown_engine_is_rejected_before_loading_a_runtime(self):
        source = replace(self.source, request=replace(self.source.request, engine="external-archive"))
        with patch.object(runtime, "get_backend", side_effect=AssertionError("Runtime must not load")) as load:
            for factory in (runtime.create_context_ensemble_job, runtime.SavedSwitchEnsembleJob):
                with self.subTest(factory=factory.__name__), self.assertRaisesRegex(ValueError, "native-hybrid"):
                    factory(self.saved, self.analysis, source, self.candidate, "/fake")
            load.assert_not_called()

    def run_job(self, side_effect=None):
        job = self.job()
        with patch.object(ensemble_runtime, "analyze_structure_ensemble",
                          side_effect=side_effect, return_value=self.artifact):
            job.start()
            job._thread.join(5)
        self.assertFalse(job._thread.is_alive())
        return job

    def test_publishes_in_context_folder_without_mutating_original_analysis(self):
        before = self.analysis_path.read_bytes()
        job = self.run_job()
        self.assertEqual(job.snapshot().state, "done")
        stored = json.loads(self.target.read_text())
        self.assertEqual(stored["inputs"], self.artifact["inputs"])
        self.assertEqual(stored["context_source"], job.analysis_identity)
        self.assertEqual(self.analysis_path.read_bytes(), before)
        self.assertIsNone(runtime.active_context_ensemble(self.saved))
        self.assertFalse((self.saved.path / "structure_ensembles").exists())

    def test_saved_absolute_path_is_ignored_and_scientific_inputs_are_checked(self):
        supplied = {**self.analysis, "result_path": str(self.folder / "other" / "result.json")}
        source, candidate = runtime.load_context_ensemble_source(self.saved, supplied)
        self.assertEqual(source.output_dir, self.analysis_folder)
        with self.assertRaises(ValueError):
            runtime.create_context_ensemble_job(self.saved, self.analysis, source,
                replace(candidate, switch_seq="ACGT"), "/fake")
        with self.assertRaises(ValueError):
            runtime.create_context_ensemble_job(self.saved, self.analysis,
                replace(source, output_dir=self.folder / "outside"), candidate, "/fake")

    def test_navigation_recovers_worker_and_prevents_replacement_and_duplicate_start(self):
        started, release = threading.Event(), threading.Event()
        duplicate = self.job()

        def compute(*args, **kwargs):
            started.set()
            release.wait(5)
            return self.artifact

        with patch.object(ensemble_runtime, "analyze_structure_ensemble", side_effect=compute):
            job = self.job()
            job.start()
            try:
                self.assertTrue(started.wait(3))
                session = {"active_design_path": str(self.root.parent)}
                self.assertTrue(design_session.active_job_running(session))
                self.assertIs(session["ensemble_job"], job)
                self.assertEqual(session["ensemble_job_key"], job.key)
                self.assertIs(runtime.active_context_ensemble(self.saved.path), job)
                self.assertFalse(design_session.active_job_running({"active_design_path": str(self.folder / "another")}))
                with self.assertRaises(ValueError):
                    duplicate.start()
                with self.assertRaises(ValueError):
                    library.save_switch(self.root, self.saved.name, self.saved.snapshot, overwrite=True)
            finally:
                release.set()
                job._thread.join(5)
        self.assertEqual(job.snapshot().state, "done")

    def test_candidate_replaced_before_start_cannot_start_native_work(self):
        job = self.job()
        replacement = deepcopy(self.saved.snapshot)
        replacement["snapshot_id"] = uuid.uuid4().hex
        library.save_switch(self.root, self.saved.name, replacement, overwrite=True)
        with patch.object(ensemble_runtime, "analyze_structure_ensemble") as compute:
            with self.assertRaises(ValueError):
                job.start()
            with self.assertRaises(ValueError):
                self.job()
        compute.assert_not_called()

    def test_external_revision_during_work_never_receives_old_families(self):
        def compute(*args, **kwargs):
            data = deepcopy(self.saved.snapshot)
            data["snapshot_id"] = uuid.uuid4().hex
            self.saved.manifest_path.write_text(json.dumps(data))
            return self.artifact

        job = self.run_job(compute)
        self.assertEqual(job.snapshot().state, "failed")
        self.assertIn("replaced", job.snapshot().error)
        self.assertFalse(self.target.exists())

    def test_cancelled_worker_publishes_no_partial_ensemble(self):
        job = self.run_job(EnsembleCancelled())
        self.assertEqual(job.snapshot().state, "cancelled")
        self.assertFalse(self.target.exists())
        self.assertIsNone(runtime.active_context_ensemble(self.saved))

    def test_moving_candidate_keeps_full_context_and_family_location(self):
        self.run_job()
        moved = self.folder / "moved"
        self.saved.path.rename(moved)
        restored = library.load_switch(moved)
        source, candidate = runtime.load_context_ensemble_source(restored, self.analysis)
        path = ensemble_path(source.output_dir, ensemble_key(ensemble_inputs(source.request, candidate, backend=self.backend)))
        self.assertTrue(path.is_relative_to(moved))
        self.assertEqual(json.loads(path.read_text())["inputs"], self.artifact["inputs"])

    def test_foreign_identity_traversal_and_linked_artifacts_are_rejected(self):
        for fields in ({"snapshot_id": uuid.uuid4().hex}, {"run_id": "../reference_1"},
                       {"run_id": str(self.analysis_folder)}, {"run_id": "."}):
            with self.subTest(fields=fields), self.assertRaises(ValueError):
                runtime.load_context_ensemble_source(self.saved, {**self.analysis, **fields})
        outside = self.folder / "foreign"
        outside.mkdir()
        cache = self.target.parent
        cache.symlink_to(outside, target_is_directory=True)
        with self.assertRaises(ValueError):
            runtime.load_context_ensemble_source(self.saved, self.analysis)
        cache.unlink()
        cache.mkdir()
        self.target.symlink_to(outside / "foreign.json")
        with self.assertRaises(ValueError):
            runtime.load_context_ensemble_source(self.saved, self.analysis)

    def test_symlink_created_during_calculation_is_not_followed_on_publication(self):
        foreign = self.folder / "foreign.json"
        foreign.write_text("untouched")

        def compute(*args, **kwargs):
            self.target.parent.mkdir()
            self.target.symlink_to(foreign)
            return self.artifact

        job = self.run_job(compute)
        self.assertEqual(job.snapshot().state, "failed")
        self.assertEqual(foreign.read_text(), "untouched")

    def test_failed_publication_and_thread_start_release_registry(self):
        with patch.object(runtime, "_atomic_json", side_effect=OSError("disk full")):
            job = self.run_job()
        self.assertEqual(job.snapshot().state, "failed")
        self.assertIn("disk full", job.snapshot().error)
        self.assertIsNone(runtime.active_context_ensemble(self.saved))
        job = self.job()
        with patch.object(job._thread, "start", side_effect=RuntimeError("thread failed")):
            with self.assertRaises(RuntimeError):
                job.start()
        self.assertIsNone(runtime.active_context_ensemble(self.saved))


if __name__ == "__main__":
    unittest.main()
