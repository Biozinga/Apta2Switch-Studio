"""Recoverable design workers, trial-boundary controls, and process ownership."""
from __future__ import annotations

from copy import deepcopy
from dataclasses import replace
import importlib
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import threading
import unittest
from unittest.mock import patch

from aptaswitch_core.architecture import default_architecture_for_trigger
from aptaswitch_core.models import RunRequest, SequenceInput
from aptaswitch_studio import design_job_registry as registry
from aptaswitch_studio.run_library import LoadedRun
from aptaswitch_studio.web_runtime import DesignRunJob
from tests.design_fixtures import design_result_fixture


class DesignContinuationRuntimeTests(unittest.TestCase):
    def setUp(self):
        self.folder = tempfile.TemporaryDirectory()
        self.addCleanup(self.folder.cleanup)
        self.root = Path(self.folder.name).resolve()
        self.request = RunRequest(SequenceInput('ACGTACGT', 'ACGT'*6),
            default_architecture_for_trigger(24), trials=3, output_dir=self.root)
        self.jobs = []
        self.addCleanup(self.clean_jobs)

    def clean_jobs(self):
        for job in self.jobs:
            registry.release_design_job(job, remove=True)

    def result(self, *, status='running', finished=1, phase='design'):
        result = design_result_fixture(self.request)
        result.output_dir = self.root/'designs'/'test_fixture'
        result.status = status
        result.candidates = []  # Trials can legitimately produce no candidate.
        result.optimization['trials'] = result.optimization['trials'][:finished]
        result.optimization['checkpoint'] = {'version': 1, 'phase': phase,
                                              'completed_trials': finished}
        return result

    def job(self, *args, **kwargs):
        job = DesignRunJob(*(args or (self.request,)), **kwargs)
        self.jobs.append(job)
        return job

    def finish(self, job):
        job._thread.join(timeout=3)
        self.assertFalse(job._thread.is_alive())

    def persist(self, result):
        path = result.output_dir/'design_checkpoint.json'
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(result.to_dict()))
        return path

    def test_pause_completes_current_trial_without_cancellation_and_is_reattachable(self):
        started, finish_trial = threading.Event(), threading.Event()
        def workflow(req, *, progress_callback, cancel_check, control_check, resume_result, checkpoint_callback):
            checkpoint_callback(self.result(finished=0))
            started.set()
            self.assertTrue(finish_trial.wait(2))
            self.assertFalse(cancel_check())
            self.assertEqual(control_check(), 'pause')
            result = self.result(status='paused', finished=1)
            self.persist(result)
            checkpoint_callback(result)
            return result, {'checkpoint': str(result.output_dir/'design_checkpoint.json')}
        with patch('aptaswitch_studio.web_runtime.run_and_export', side_effect=workflow):
            job = self.job()
            job.start()
            self.assertTrue(started.wait(2))
            self.assertIs(registry.find_design_job(self.root), job)
            self.assertIs(registry.find_design_job(self.root/'designs'), job)
            job.pause()
            self.assertEqual(job.snapshot().control_requested, 'pause')
            self.assertEqual(job.snapshot().state, 'running')
            finish_trial.set()
            self.finish(job)
        snapshot = job.snapshot()
        self.assertEqual(snapshot.state, 'paused')
        self.assertEqual(snapshot.current, 1)
        self.assertEqual(snapshot.total, 3)
        self.assertEqual(len(snapshot.result.candidates), 0)
        self.assertTrue(job.consume_once())
        self.assertFalse(job.consume_once())
        self.assertIs(registry.find_design_job(job.checkpoint_path), job)

    def test_finish_requested_before_start_is_visible_before_another_trial(self):
        def workflow(req, **kwargs):
            self.assertEqual(kwargs['control_check'](), 'finish')
            self.assertFalse(kwargs['cancel_check']())
            result = self.result(status='stopped', finished=0, phase='done')
            return result, {}
        with patch('aptaswitch_studio.web_runtime.run_and_export', side_effect=workflow):
            job = self.job()
            job.finish_early()
            job.pause()  # A terminal finish request wins over pause.
            job.start()
            self.finish(job)
        self.assertEqual(job.snapshot().state, 'done')
        self.assertEqual(job.snapshot().result.status, 'stopped')

    def test_checkpoint_change_does_not_lose_pending_control_during_finalization(self):
        job = self.job()
        job.pause()
        job._checkpoint(self.result(phase='inverse_extension'))
        job._thread_started = True
        self.assertTrue(job.supports_continuation)
        self.assertEqual(job.snapshot().phase, 'inverse_extension')
        self.assertFalse(job.snapshot().can_request_control)
        self.assertEqual(job._control_check(), 'pause')
        job.finish_early()  # Cannot change controls while postprocessing runs.
        self.assertEqual(job._control_check(), 'pause')

    def test_resume_finalization_can_queue_finish_before_start_without_enabling_live_pause(self):
        result = self.result(status='paused', phase='inverse_extension')
        job = DesignRunJob(resume_result=result)
        self.jobs.append(job)
        self.assertEqual(job.snapshot().phase, 'inverse_extension')
        self.assertFalse(job.snapshot().can_request_control)
        job.finish_early()
        self.assertEqual(job.control_requested, 'finish')

    def test_registry_survives_reload_and_rejects_duplicate_launch(self):
        started, release = threading.Event(), threading.Event()
        def workflow(req, **kwargs):
            started.set()
            self.assertTrue(release.wait(2))
            return self.result(status='completed', phase='done'), {}
        with patch('aptaswitch_studio.web_runtime.run_and_export', side_effect=workflow):
            job = self.job()
            job.start()
            self.assertTrue(started.wait(2))
            importlib.reload(registry)
            self.assertIs(registry.find_design_job(self.root), job)
            duplicate = self.job()
            with self.assertRaises(ValueError):
                duplicate.start()
            release.set()
            self.finish(job)

    def test_failed_thread_start_releases_registry_and_process_locks(self):
        job = self.job()
        with patch.object(job._thread, 'start', side_effect=RuntimeError('cannot start')):
            with self.assertRaisesRegex(RuntimeError, 'cannot start'):
                job.start()
        self.assertIsNone(registry.find_design_job(self.root))
        lock = registry.DesignFileLock(self.root/'.design.lock')
        self.assertTrue(lock.acquire())
        lock.close()

    def test_second_process_cannot_acquire_active_checkpoint_lock(self):
        path = self.root/'designs'/'run'/'design.worker.lock'
        lock = registry.DesignFileLock(path)
        self.assertTrue(lock.acquire())
        try:
            script = ('from pathlib import Path; from aptaswitch_studio.design_job_registry import DesignFileLock; '
                      'import sys; lock=DesignFileLock(Path(sys.argv[1])); print(lock.acquire()); lock.close()')
            child = subprocess.run([sys.executable, '-c', script, str(path)],
                capture_output=True, text=True, check=True, timeout=10)
            self.assertEqual(child.stdout.strip(), 'False')
        finally:
            lock.close()
        other = registry.DesignFileLock(path)
        self.assertTrue(other.acquire())
        other.close()

    def test_inverse_reanalysis_retains_cancel_and_has_no_fake_pause(self):
        source = self.result()
        job = self.job(analysis_source=source)
        self.assertFalse(job.supports_continuation)
        job.pause()
        job.finish_early()
        self.assertIsNone(job.control_requested)
        job.cancel()
        self.assertTrue(job._cancel.is_set())
        with patch('aptaswitch_core.workflow.reanalyze_inverse_and_export',
                   return_value=(self.result(status='cancelled'), {})) as reanalyze:
            job.start()
            self.finish(job)
        self.assertNotIn('control_check', reanalyze.call_args.kwargs)

    def test_resume_reloads_after_lock_and_uses_frozen_request_with_updated_library(self):
        saved = self.result(status='paused')
        path = self.persist(saved)
        loaded = LoadedRun('design', saved, {}, path)
        def loader(source):
            self.assertEqual(Path(source), path)
            return deepcopy(loaded)
        def workflow(req, **kwargs):
            self.assertEqual(req.trials, 3)
            self.assertEqual(req.thermodynamics_library_path, '/new/unchanged-model.dylib')
            self.assertEqual(kwargs['resume_result'].optimization['checkpoint']['completed_trials'], 1)
            self.assertEqual(kwargs['control_check'](), 'finish')
            done = deepcopy(kwargs['resume_result'])
            done.status = 'stopped'
            return done, {}
        with patch('aptaswitch_studio.run_library.load_run', side_effect=loader) as reader, \
             patch('aptaswitch_studio.web_runtime.run_and_export', side_effect=workflow):
            job = registry.resume_design_job(path, finish_early=True, library_path='/new/unchanged-model.dylib')
            self.jobs.append(job)
            self.finish(job)
        self.assertEqual(reader.call_count, 2)
        self.assertEqual(job.snapshot().state, 'done')

    def test_resume_refuses_checkpoint_completed_between_click_and_lock(self):
        saved = self.result(status='paused')
        path = self.persist(saved)
        final = self.result(status='completed', phase='done')
        with patch('aptaswitch_studio.run_library.load_run', side_effect=[
            LoadedRun('design', saved, {}, path), LoadedRun('design', final, {}, path)]):
            with self.assertRaises(ValueError):
                registry.resume_design_job(path)
        self.assertIsNone(registry.find_design_job(self.root))

    def test_listing_uses_strict_reader_and_excludes_terminal_or_active_checkpoints(self):
        saved = self.result(status='running')
        path = self.persist(saved)
        with patch('aptaswitch_studio.run_library.load_run',
                   return_value=LoadedRun('design', saved, {}, path)) as reader:
            found = registry.list_design_checkpoints(self.root)
            self.assertEqual(len(found), 1)
            self.assertEqual(found[0].state, 'interrupted')
            self.assertEqual(found[0].completed_trials, 1)
            self.assertEqual(found[0].total_trials, 3)
            reader.assert_called_once_with(path)
            saved.optimization['checkpoint']['phase'] = 'done'
            saved.status = 'completed'
            self.assertEqual(registry.list_design_checkpoints(self.root), [])
        lock = registry.DesignFileLock(path.parent/'design.worker.lock')
        self.assertTrue(lock.acquire())
        try:
            with patch('aptaswitch_studio.run_library.load_run') as reader:
                self.assertEqual(registry.list_design_checkpoints(self.root), [])
                reader.assert_not_called()
        finally:
            lock.close()

    def test_failed_inverse_checkpoint_can_be_found_and_resumed(self):
        saved = self.result(status='failed', finished=3, phase='inverse_extension')
        saved.inverse_extension_analysis = {'status': 'failed', 'error': 'temporary failure'}
        path = self.persist(saved)
        def loader(source):
            return LoadedRun('design', deepcopy(saved), {}, path)
        def workflow(req, **kwargs):
            self.assertEqual(kwargs['resume_result'].optimization['checkpoint']['phase'], 'inverse_extension')
            done = deepcopy(kwargs['resume_result'])
            done.status = 'completed'
            done.inverse_extension_analysis = {'status': 'completed'}
            return done, {}
        with patch('aptaswitch_studio.run_library.load_run', side_effect=loader), \
             patch('aptaswitch_studio.web_runtime.run_and_export', side_effect=workflow):
            found = registry.list_design_checkpoints(self.root)
            self.assertEqual(len(found), 1)
            job = registry.resume_design_job(found[0].path)
            self.jobs.append(job)
            self.finish(job)
            self.assertEqual(job.snapshot().state, 'done')

    def test_request_is_detached_from_mutable_ui_options(self):
        request = replace(self.request, optimization_options={'seed': 12})
        job = self.job(request)
        request.optimization_options['seed'] = 999
        self.assertEqual(job.request.optimization_options['seed'], 12)
        saved = self.result(status='paused')
        with self.assertRaisesRegex(ValueError, 'saved design request'):
            DesignRunJob(replace(saved.request, trials=9), resume_result=saved)


@unittest.skipUnless(os.environ.get('APTASWITCH_HYBRID_NATIVE_LIBRARY'), 'Explicit native library required')
class NativeDesignRuntimeContinuationTests(unittest.TestCase):
    def test_real_worker_pause_disk_reload_and_complete_remaining_trial(self):
        from aptaswitch_core.models import ReporterContext, SwitchArchitecture, ThermoConditions
        from aptaswitch_core.thermodynamics.optimization import TrialOptions
        from dataclasses import asdict
        with tempfile.TemporaryDirectory() as folder:
            request = RunRequest(SequenceInput('GGTT', 'AACC'),
                SwitchArchitecture(2, 2, upper_stem_length=1, upper_stem2_length=0,
                                   rbs_loop='AGGA', t7_leader='', rbs_prefix_length=0),
                reporter=ReporterContext(key='short', label='Short test reporter', sequence_after_start_rna='GCU'),
                thermo=ThermoConditions(temperature_c=37, validation_max_complex_size=2),
                trials=2, output_dir=Path(folder),
                thermodynamics_library_path=os.environ['APTASWITCH_HYBRID_NATIVE_LIBRARY'],
                exclude_stop_candidates=False, exclude_non_linear_candidates=False,
                optimization_options=asdict(TrialOptions(seed=17, max_evaluations=1, max_proposals=1)))
            job = DesignRunJob(request)
            original = job._checkpoint
            def pause_after_one(result):
                original(result)
                if registry.completed_trials(result) == 1 and result.status == 'running':
                    job.pause()
            job._checkpoint = pause_after_one
            continuation = None
            try:
                job.start()
                job._thread.join(30)
                self.assertFalse(job._thread.is_alive())
                paused = job.snapshot()
                self.assertEqual(paused.state, 'paused', paused.error)
                self.assertEqual(paused.current, 1)
                first_trial = deepcopy(paused.result.optimization['trials'][0])
                found = registry.list_design_checkpoints(request.output_dir)
                self.assertEqual(len(found), 1)
                continuation = registry.resume_design_job(found[0].path)
                continuation._thread.join(30)
                self.assertFalse(continuation._thread.is_alive())
                final = continuation.snapshot()
                self.assertEqual(final.state, 'done', final.error)
                self.assertEqual(final.result.status, 'completed')
                self.assertEqual(final.result.run_id, paused.result.run_id)
                self.assertEqual(final.current, 2)
                self.assertEqual(final.result.optimization['trials'][0], first_trial)
                self.assertEqual([trial['trial'] for trial in final.result.optimization['trials']], [1, 2])
                self.assertEqual(registry.list_design_checkpoints(request.output_dir), [])
            finally:
                registry.release_design_job(job, remove=True)
                if continuation is not None:
                    registry.release_design_job(continuation, remove=True)


if __name__ == '__main__':
    unittest.main()
