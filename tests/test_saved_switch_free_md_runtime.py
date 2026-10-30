"""Free-MD ownership, pause/relaunch/resume and immutable saved conditions."""
from __future__ import annotations
import importlib
import json
from pathlib import Path
import tempfile
import threading
from types import SimpleNamespace
import unittest
from unittest.mock import patch

from aptaswitch_core.models import RunRequest, SequenceInput, SwitchArchitecture
from aptaswitch_studio import saved_switch_free_md_runtime as runtime
from aptaswitch_studio import saved_switch_library as library
from aptaswitch_studio.design_session import active_job_running, reset_design_session
from aptaswitch_studio.web_runtime import JobSnapshot
from design_fixtures import design_result_fixture


class FreeMDRuntimeTests(unittest.TestCase):
    def setUp(self):
        from tests.native_ui_fixtures import patch_native_runtime
        self.native_backend_fixture = patch_native_runtime(self)
        self.root = Path(self.enterContext(tempfile.TemporaryDirectory())).resolve()
        self.design = self.root / 'Design'
        request = RunRequest(sequences=SequenceInput('AACCGGTTAACC', 'AACC'),
            architecture=SwitchArchitecture(toehold_length=2, stem_length=2), trials=1, output_dir=self.root)
        self.result = design_result_fixture(request)
        self.saved = self.save()
        self.frozen = {'options': {'replicas': 4, 'monitored_bases': [2, 3], 'duration_us': .1},
                       'runtime': {'home': '/oxdna'}, 'aptamer_sequence': 'AACCGGTTAACC'}
        self.prepare = self.enterContext(patch.object(runtime, 'prepare_free_md_request', return_value=self.frozen))
        self.native = self.enterContext(patch.object(runtime, '_native_running', return_value=False))
        self.addCleanup(runtime._JOBS.clear)

    def save(self, overwrite=False):
        snapshot = library.candidate_snapshot(self.result, self.result.candidates[0], 'Retained',
            original_sequences={'aptamer_dna': 'AACCGGTTAACC', 'trigger_dna': 'CC'})
        return library.save_switch(library.saved_switch_root(self.design), 'Retained', snapshot, overwrite=overwrite)

    def paused(self):
        job = runtime.SavedSwitchFreeMDJob(self.saved, self.frozen)
        def calculate(request, *, output_dir, progress_callback, cancel_check):
            (output_dir / 'resume_state.json').write_text('{}')
            (output_dir / 'worker_request.json').write_text(json.dumps(request))
            progress_callback(1, 4, 'sampled')
            return {'status': 'cancelled', 'sampling_complete': False,
                    'replicas': [], 'comparison': {'status': 'partial'}}
        with patch.object(runtime, 'run_free_md', side_effect=calculate):
            job._run()
        return job

    def test_results_and_frozen_inputs_survive_move_without_native_compute(self):
        before = self.saved.manifest_path.read_bytes()
        job = self.paused()
        self.assertEqual(job.snapshot().state, 'done')
        self.assertEqual(self.saved.manifest_path.read_bytes(), before)
        moved = self.root / 'Moved'
        self.design.rename(moved)
        saved = library.load_switch(moved / library.LIBRARY_FOLDER / self.saved.name)
        results = runtime.saved_free_md_results(saved)
        self.assertEqual(len(results), 1)
        self.assertEqual(results[0]['requested_parameters'], self.frozen)
        self.assertEqual(results[0]['status'], 'cancelled')
        self.assertTrue(results[0]['resumable'])
        self.assertTrue(Path(results[0]['run_dir']).is_relative_to(moved))

    def test_paused_run_allows_new_run_then_resumes_its_own_frozen_request(self):
        paused = self.paused()
        with patch.object(runtime.SavedSwitchFreeMDJob, 'start'):
            second, started = runtime.start_free_md(self.saved, {}, {'replicas': 10})
        self.assertTrue(started)
        self.assertNotEqual(second.folder, paused.folder)
        self.assertEqual(len(runtime.saved_free_md_results(self.saved)), 1)
        second.cancel()
        with patch.object(runtime, 'run_free_md', return_value={'status': 'cancelled'}):
            second._run()
        self.prepare.reset_mock()
        with patch('aptaswitch_core.oxdna_free_md.load_free_md_resume_request', return_value=self.frozen), \
                patch.object(runtime.SavedSwitchFreeMDJob, 'start'):
            resumed, started = runtime.resume_free_md(self.saved, paused.folder)
        self.assertTrue(started)
        self.assertEqual(resumed.folder, paused.folder)
        self.assertEqual(resumed.request['options']['replicas'], 4)
        self.prepare.assert_not_called()

    def test_active_run_blocks_candidate_overwrite_navigation_and_ffs(self):
        from aptaswitch_studio.saved_switch_ffs_runtime import start_ffs
        with patch.object(runtime.SavedSwitchFreeMDJob, 'start'):
            job, _ = runtime.start_free_md(self.saved, {}, {})
        with self.assertRaises(ValueError):
            self.save(overwrite=True)
        session = {'active_design_path': str(self.design)}
        self.assertTrue(active_job_running(session))
        self.assertIs(session['saved_switch_free_md_job'], job)
        with self.assertRaises(ValueError):
            reset_design_session(session)
        with self.assertRaises(ValueError):
            start_ffs(self.saved, {}, {})

    def test_context_analysis_and_free_md_cannot_run_on_same_candidate(self):
        from aptaswitch_studio import saved_switch_runtime as context
        with patch.object(runtime.SavedSwitchFreeMDJob, 'start'):
            job, _ = runtime.start_free_md(self.saved, {}, {})
        with self.assertRaises(ValueError):
            context.start_analysis(self.saved, {})
        runtime._JOBS.clear()
        running = SimpleNamespace(snapshot=lambda: JobSnapshot(state='running', current=0, total=1,
            message='', elapsed_seconds=0., log=()))
        with patch.object(context, 'active_analysis', return_value=running):
            with self.assertRaises(ValueError):
                runtime.start_free_md(self.saved, {}, {})

    def test_competing_context_starts_preserve_winning_frozen_configuration(self):
        from aptaswitch_studio import saved_switch_runtime as context
        barrier, release = threading.Barrier(2), threading.Event()
        outcomes, errors = [], []
        prepare = context._context_record
        def prepare_together(*args, **kwargs):
            record = prepare(*args, **kwargs)
            barrier.wait(timeout=3)
            return record
        def compute(*args, **kwargs):
            release.wait(3)
            return {'status': 'completed'}
        def start(temperature):
            try:
                outcomes.append(context.start_analysis(self.saved, {'temperature_c': temperature}))
            except Exception as exc:
                errors.append(exc)
        try:
            with patch.object(context, '_context_record', side_effect=prepare_together), \
                    patch.object(context, 'analyze_saved_switch', side_effect=compute):
                threads = [threading.Thread(target=start, args=(temperature,)) for temperature in (30., 40.)]
                for thread in threads:
                    thread.start()
                for thread in threads:
                    thread.join(4)
                self.assertFalse(any(thread.is_alive() for thread in threads))
                self.assertEqual(errors, [])
                self.assertEqual(len(outcomes), 2)
                self.assertEqual(sum(started for _, started in outcomes), 1)
                self.assertIs(outcomes[0][0], outcomes[1][0])
                self.assertEqual(context.load_context(self.saved)['config'], outcomes[0][0].config)
                release.set()
                outcomes[0][0]._thread.join(4)
        finally:
            release.set()
            context._JOBS.clear()

    def test_native_worker_recovered_and_pause_targets_same_folder(self):
        job = self.paused()
        self.native.return_value = True
        recovered = runtime.active_free_md(self.saved)
        self.assertIsInstance(recovered, runtime.DetachedFreeMDJob)
        self.assertEqual(recovered.folder, job.folder)
        recovered.cancel()
        self.assertTrue((job.folder / 'CANCEL').is_file())

    def test_foreign_identity_and_outside_resume_are_rejected(self):
        job = self.paused()
        metadata = json.loads((job.folder / 'request.json').read_text())
        metadata['candidate_identity'] = 'another candidate'
        (job.folder / 'request.json').write_text(json.dumps(metadata))
        self.assertEqual(runtime.saved_free_md_results(self.saved), [])
        with self.assertRaises(ValueError):
            runtime.resume_free_md(self.saved, self.root)
