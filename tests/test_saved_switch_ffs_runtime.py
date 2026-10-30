"""FFS jobs preserve candidate identity, cancellation and portable artifacts."""
from __future__ import annotations

import builtins
import json
import importlib
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from aptaswitch_core.models import RunRequest, SequenceInput, SwitchArchitecture
from aptaswitch_studio import saved_switch_ffs_runtime as runtime
from aptaswitch_studio import saved_switch_library as library
from aptaswitch_studio.design_session import active_job_running, reset_design_session
from design_fixtures import design_result_fixture


class SavedSwitchFFSRuntimeTests(unittest.TestCase):
    def setUp(self):
        self.folder = Path(self.enterContext(tempfile.TemporaryDirectory())).resolve()
        self.root = self.folder / 'design' / library.LIBRARY_FOLDER
        request = RunRequest(sequences=SequenceInput('AACCGGTTAACC', 'AACC'),
                             architecture=SwitchArchitecture(toehold_length=2, stem_length=2),
                             trials=1, output_dir=self.folder / 'original')
        self.result = design_result_fixture(request)
        self.saved = self.save()
        self.options = {'temperature_c': 37., 'salt_m': .12, 'kinetics_trials': 1000,
                        'kinetics_source_samples': 200, 'kinetics_replicas': 4, 'kinetics_workers': 4}
        self.request = {'sampling': {'budget_seconds': None}, 'input_marker': 'exact candidate',
                        'kinetics': {'version': 3, 'mode': 'opening', 'trials': 1000, 'source_samples': 200,
                                     'replicas': 4, 'workers': 4}}
        self.intent = {'runtime': {'home': '/oxdna'}, 'options': self.options}
        self.prepare = self.enterContext(patch.object(runtime, 'prepare_ffs_request', return_value=self.request))
        self.addCleanup(runtime._JOBS.clear)

    def save(self, overwrite=False):
        snapshot = library.candidate_snapshot(self.result, self.result.candidates[0], 'Retained')
        return library.save_switch(self.root, snapshot['name'], snapshot, overwrite=overwrite)

    def start(self):
        return runtime.start_ffs(self.saved, {'home': '/oxdna'}, self.options)

    def run_job(self, status='completed'):
        job = runtime.SavedSwitchFFSJob(self.saved, self.intent)

        def calculate(request, *, output_dir, progress_callback, cancel_check):
            (output_dir / 'progress.json').write_text(json.dumps({'phase': 'native sampling', 'stages': [{'crossings': 3}]}))
            progress_callback(2, 4, 'sampling')
            (output_dir / 'movie.dat').write_text('trajectory')
            (output_dir / 'preparation').mkdir()
            (output_dir / 'preparation' / 'system.top').write_text('topology')
            return {'status': status, 'movie_file': 'movie.dat', 'topology_file': 'preparation/system.top', 'stages': [{'complete': True}]}

        with patch.object(runtime, 'run_oxdna_ffs', side_effect=calculate):
            job._run()
        return job

    def test_request_result_and_movie_survive_folder_move_without_any_engine(self):
        before = self.saved.manifest_path.read_bytes()
        job = self.run_job()
        self.assertEqual(job.snapshot().state, 'done')
        self.assertEqual(json.loads((job.folder / 'request.json').read_text())['request'], self.request)
        self.assertEqual(json.loads((job.folder / 'job_progress.json').read_text())['message'], 'sampling')
        self.assertEqual(json.loads((job.folder / 'progress.json').read_text()),
                         {'phase': 'native sampling', 'stages': [{'crossings': 3}]})
        self.assertEqual(self.saved.manifest_path.read_bytes(), before)
        moved = self.folder / 'moved'
        self.root.parent.rename(moved)
        saved = library.load_switch(moved / library.LIBRARY_FOLDER / self.saved.name)
        with patch.object(runtime, 'run_oxdna_ffs', side_effect=AssertionError('must not compute')):
            records = runtime.saved_ffs_results(saved)
        self.assertEqual(len(records), 1)
        result = records[0]
        self.assertEqual(result['status'], 'completed')
        self.assertTrue(Path(result['run_dir']).is_relative_to(moved))
        self.assertEqual((Path(result['run_dir']) / result['movie_file']).read_text(), 'trajectory')
        self.assertEqual((Path(result['run_dir']) / result['topology_file']).read_text(), 'topology')

    def test_partial_and_cancelled_runs_remain_terminal_without_claiming_success(self):
        for status in ('cancelled', 'failed'):
            job = self.run_job(status)
            with self.subTest(status=status):
                self.assertEqual(job.snapshot().state, 'failed' if status == 'failed' else 'done')
                self.assertEqual(job.snapshot().result['status'], status)
                self.assertTrue(job.consume_once())
                self.assertFalse(job.consume_once())
                self.assertEqual(runtime.saved_ffs_results(self.saved)[0]['status'], status)

    def test_proximity_four_archive_is_detected_from_its_frozen_protocol(self):
        self.request['kinetics'].update(version=4, mode='proximity', source_samples=0,
            equilibration_steps=200000, initial_distance_min=2., initial_distance_max=2.5,
            equilibration_method='isolated_monomers_before_close_placement_v1')
        job = self.run_job()
        path = job.folder / 'request.json'
        metadata = json.loads(path.read_text())
        metadata.pop('method')
        path.write_text(json.dumps(metadata))
        records = runtime.saved_ffs_results(self.saved)
        self.assertEqual(len(records), 1)
        self.assertEqual(records[0]['simulation_mode'], 'proximity')
        self.assertIsNone(records[0]['rate_estimate'])

    def test_fresh_full_simulations_do_not_reuse_previous_sources(self):
        self.assertEqual(runtime.saved_ffs_results(self.saved), [])
        jobs = [self.run_job(), self.run_job()]
        self.assertNotEqual(jobs[0].folder, jobs[1].folder)
        for job in jobs:
            metadata = json.loads((job.folder / 'request.json').read_text())
            self.assertEqual(metadata['method'], 'NATIVE_FFS')
            self.assertEqual(metadata['request']['kinetics']['mode'], 'opening')
            self.assertNotIn('source_run_dir', metadata['request'])
        records = runtime.saved_ffs_results(self.saved)
        self.assertEqual(len(records), 2)
        self.assertEqual({r['method'] for r in records}, {'NATIVE_FFS'})

    def test_pause_resumes_frozen_conditions_when_form_has_changed(self):
        original = self.resumable_job()
        self.options['temperature_c'] = 42.
        self.options['kinetics_replicas'] = 10
        self.prepare.reset_mock()
        with patch.object(runtime.SavedSwitchFFSJob, 'start'):
            resumed, started = runtime.resume_ffs(self.saved, original.folder)
        self.assertTrue(started)
        with patch.object(runtime, 'run_oxdna_ffs', return_value={'status': 'completed'}) as compute:
            resumed._run()
        self.prepare.assert_not_called()
        self.assertEqual(compute.call_args.args[0]['kinetics']['replicas'], 4)
        self.assertEqual(resumed.snapshot().result['method'], 'NATIVE_FFS')

    def test_first_contact_only_archive_is_not_presented_as_complete_pathway(self):
        original = self.resumable_job()
        path = original.folder / 'worker_request.json'
        request = json.loads(path.read_text())
        request['kinetics']['mode'] = 'encounter'
        path.write_text(json.dumps(request))
        self.assertEqual(runtime.saved_ffs_results(self.saved), [])
        self.assertTrue(path.is_file())

    def test_obsolete_protocol_cannot_resume_in_native_ffs(self):
        self.request['kinetics'].update(mode='opening', version=1)
        original = self.resumable_job()
        records = runtime.saved_ffs_results(self.saved)
        self.assertEqual(len(records), 1)
        self.assertEqual(records[0]['simulation_mode'], 'opening')
        self.assertFalse(records[0]['resumable'])
        self.assertTrue((original.folder / 'resume_state.json').is_file())
        with self.assertRaises(ValueError):
            runtime.resume_ffs(self.saved, original.folder)

    def test_cancel_before_execution_never_starts_native_engine(self):
        job = runtime.SavedSwitchFFSJob(self.saved, self.intent)
        job.cancel()
        with patch.object(runtime, 'run_oxdna_ffs') as compute:
            job._run()
        compute.assert_not_called()
        self.assertEqual(job.snapshot().result['status'], 'cancelled')

    def test_cancellation_is_forwarded_during_execution(self):
        job = runtime.SavedSwitchFFSJob(self.saved, self.intent)

        def calculate(request, *, output_dir, progress_callback, cancel_check):
            self.assertFalse(cancel_check())
            job.cancel()
            self.assertTrue(cancel_check())
            return {'status': 'cancelled'}

        with patch.object(runtime, 'run_oxdna_ffs', side_effect=calculate):
            job._run()
        self.assertEqual(job.snapshot().result['status'], 'cancelled')

    def test_cancel_during_request_preparation_prevents_native_launch(self):
        job = runtime.SavedSwitchFFSJob(self.saved, self.intent)
        def prepare(*args):
            job.cancel()
            return self.request
        with patch.object(runtime, 'prepare_ffs_request', side_effect=prepare), \
             patch.object(runtime, 'run_oxdna_ffs') as compute:
            job._run()
        compute.assert_not_called()
        self.assertEqual(job.snapshot().result['status'], 'cancelled')

    def test_active_preparation_is_not_also_listed_as_interrupted(self):
        with patch.object(runtime.SavedSwitchFFSJob, 'start'):
            job, _ = self.start()

        def prepare(*args):
            self.assertTrue((job.folder / 'request.json').is_file())
            self.assertFalse((job.folder / 'worker.lock').exists())
            self.assertEqual(runtime.saved_ffs_results(self.saved), [])
            return self.request

        with patch.object(runtime, 'prepare_ffs_request', side_effect=prepare), \
                patch.object(runtime, 'run_oxdna_ffs', return_value={'status': 'completed'}):
            job._run()
        self.assertEqual(runtime.saved_ffs_results(self.saved)[0]['status'], 'completed')


    def test_active_job_is_recovered_and_prevents_navigation_and_overwrite(self):
        with patch.object(runtime.SavedSwitchFFSJob, 'start'):
            job, started = self.start()
            again, second = self.start()
        self.assertTrue(started)
        self.assertFalse(second)
        self.assertIs(again, job)
        self.prepare.assert_not_called()
        session = {'active_design_path': str(self.root.parent)}
        self.assertTrue(active_job_running(session))
        self.assertIs(session['saved_switch_ffs_job'], job)
        with self.assertRaises(ValueError):
            reset_design_session(session)
        original = self.saved.manifest_path.read_bytes()
        with self.assertRaises(ValueError):
            self.save(overwrite=True)
        self.assertEqual(self.saved.manifest_path.read_bytes(), original)
        with patch.object(runtime, 'run_oxdna_ffs', return_value={'status': 'completed'}):
            job._run()
        self.assertIsNone(runtime.active_ffs(self.saved))
        self.assertFalse(active_job_running(session))

    def test_replaced_candidate_refuses_stale_start_and_stale_result_publication(self):
        job = runtime.SavedSwitchFFSJob(self.saved, self.intent)
        self.save(overwrite=True)
        with self.assertRaises(ValueError):
            self.start()
        with patch.object(runtime, 'run_oxdna_ffs') as compute:
            job._run()
        compute.assert_not_called()
        self.assertEqual(job.snapshot().state, 'failed')
        self.assertFalse((job.folder / 'result.json').exists())

    def test_start_needs_no_pdb_and_prepares_the_immutable_candidate_in_worker(self):
        self.assertFalse((self.saved.path / 'calculations' / 'structures3d').exists())
        with patch.object(runtime.SavedSwitchFFSJob, 'start'):
            job, _ = self.start()
        self.prepare.assert_not_called()
        with patch.object(runtime, 'run_oxdna_ffs', return_value={'status': 'completed'}):
            job._run()
        self.prepare.assert_called_once_with(self.saved.snapshot, {'home': '/oxdna'}, self.options)
        self.assertEqual(job.snapshot().result['status'], 'completed')

    def test_active_context_analysis_still_blocks_start_and_resume(self):
        original = self.resumable_job()
        checkpoint = original.folder / 'resume_state.json'
        before = checkpoint.read_bytes()
        self.prepare.reset_mock()
        with patch('aptaswitch_studio.saved_switch_runtime.active_analysis', return_value=object()), \
                patch.object(runtime.SavedSwitchFFSJob, 'start') as start:
            for launch in (self.start, lambda: runtime.resume_ffs(self.saved, original.folder)):
                with self.subTest(launch=launch), self.assertRaises(ValueError):
                    launch()
        start.assert_not_called()
        self.prepare.assert_not_called()
        self.assertEqual(checkpoint.read_bytes(), before)
        self.assertEqual(runtime._JOBS, {})

    def test_result_loading_rejects_foreign_identity_and_unsafe_movie_paths(self):
        job = self.run_job()
        manifest = job.folder / 'result.json'
        original = json.loads(manifest.read_text())
        external = self.folder / 'private.dat'
        external.write_text('PRIVATE')
        (job.folder / 'alias.dat').symlink_to(external)
        for movie in (str(external), '../private.dat', 'alias.dat'):
            manifest.write_text(json.dumps({**original, 'movie_file': movie}))
            with self.subTest(movie=movie):
                self.assertIsNone(runtime.saved_ffs_results(self.saved)[0]['movie_file'])
        manifest.write_text(json.dumps({**original, 'candidate_identity': 'another candidate'}))
        self.assertEqual(runtime.saved_ffs_results(self.saved), [])

    def test_symlinked_output_folder_is_rejected_before_any_write(self):
        outside = self.folder / 'outside'
        outside.mkdir()
        (self.saved.path / 'calculations' / 'oxdna_ffs').symlink_to(outside, target_is_directory=True)
        with self.assertRaises(ValueError):
            self.start()
        self.assertEqual(list(outside.iterdir()), [])

    def test_native_exception_is_persisted_and_unregisters_the_worker(self):
        with patch.object(runtime.SavedSwitchFFSJob, 'start'):
            job, _ = self.start()
        with patch.object(runtime, 'run_oxdna_ffs', side_effect=RuntimeError('native failure')):
            job._run()
        self.assertEqual(job.snapshot().state, 'failed')
        self.assertEqual(runtime.saved_ffs_results(self.saved)[0]['error'], 'native failure')
        self.assertIsNone(runtime.active_ffs(self.saved))

    def resumable_job(self):
        job = self.run_job('cancelled')
        (job.folder / 'preparation' / 'validated.dat').write_text('positions and velocities')
        frozen = {**self.request, 'runtime': {'home': '/oxdna'},
                  'prepared': {'configuration_file': str(job.folder / 'preparation' / 'validated.dat'),
                               'topology_file': str(job.folder / 'preparation' / 'system.top')}}
        (job.folder / 'worker_request.json').write_text(json.dumps(frozen))
        (job.folder / 'resume_state.json').write_text(json.dumps({'format': 'oxdna-native-ffs', 'version': 3}))
        return job

    def test_module_reload_and_lost_browser_session_preserve_calculation(self):
        with patch.object(runtime.SavedSwitchFFSJob, 'start'):
            job, _ = self.start()
        importlib.reload(runtime)
        self.assertIs(runtime.active_ffs(self.saved), job)
        session = {'active_design_path': str(self.root.parent)}
        self.assertTrue(active_job_running(session))
        self.assertIs(session['saved_switch_ffs_job'], job)
        self.assertFalse(job._cancel.is_set())

    def test_progress_needs_no_application_imports_after_job_creation(self):
        job = runtime.SavedSwitchFFSJob(self.saved, self.intent)
        original_import = builtins.__import__

        def reject_application_import(name, globals=None, locals=None, fromlist=(), level=0):
            package = (globals or {}).get('__package__', '') or ''
            if name.startswith('aptaswitch_') or (level and package.startswith('aptaswitch_')):
                raise AssertionError(f'Background progress tried to import {name!r}')
            return original_import(name, globals, locals, fromlist, level)

        # A Streamlit source reload may evict local modules while this thread
        # lives on. Its heartbeat must use bound functions, not import the UI
        # or rerun scientific snapshot validation.
        with patch.object(builtins, '__import__', side_effect=reject_application_import):
            job._progress(2, 4, 'sampling after a source reload')

        progress = json.loads((job.folder / 'job_progress.json').read_text())
        self.assertEqual(progress['message'], 'sampling after a source reload')
        self.assertEqual((progress['current'], progress['total']), (2, 4))
        self.assertEqual((job.snapshot().current, job.snapshot().total), (2, 4))
        self.assertFalse(job._cancel.is_set())

    def test_progress_and_publication_reject_changed_manifest_with_same_snapshot_id(self):
        job = runtime.SavedSwitchFFSJob(self.saved, self.intent)
        job._progress(1, 4, 'original candidate')
        destination = job.folder / 'job_progress.json'
        original_progress = destination.read_bytes()
        manifest = json.loads(self.saved.manifest_path.read_text())
        manifest['request']['trials'] += 1
        self.assertEqual(manifest['snapshot_id'], self.saved.snapshot['snapshot_id'])
        self.saved.manifest_path.write_text(json.dumps(manifest))

        with self.assertRaises(ValueError):
            job._progress(2, 4, 'must not publish for a changed candidate')
        self.assertEqual(destination.read_bytes(), original_progress)
        # The shared writer must also enforce full immutable identity. A
        # snapshot UUID alone cannot protect final result publication.
        with self.assertRaises(ValueError):
            runtime._atomic_json(self.saved, destination, {'message': 'foreign result'})
        self.assertEqual(destination.read_bytes(), original_progress)

    def test_resume_reuses_run_and_frozen_request_without_a_current_runtime(self):
        original = self.resumable_job()
        self.prepare.reset_mock()
        (original.folder / 'bank_00000.dat').write_text('keep this crossing')
        with patch.object(runtime.SavedSwitchFFSJob, 'start'):
            resumed, started = runtime.resume_ffs(self.saved, original.folder)
            duplicate, second = runtime.resume_ffs(self.saved, original.folder)
        self.assertTrue(started)
        self.assertFalse(second)
        self.assertIs(duplicate, resumed)
        self.assertEqual(resumed.run_id, original.run_id)
        self.assertEqual(resumed.folder, original.folder)
        with patch.object(runtime, 'run_oxdna_ffs', return_value={'status': 'completed'}) as compute:
            resumed._run()
        self.prepare.assert_not_called()
        self.assertTrue(compute.call_args.args[0]['resume'])
        self.assertEqual(compute.call_args.args[0]['kinetics'], self.request['kinetics'])
        self.assertNotIn('upgrade_diversity', compute.call_args.args[0])
        self.assertEqual((original.folder / 'bank_00000.dat').read_text(), 'keep this crossing')

    def test_abrupt_stop_without_final_result_is_visible_and_resumable(self):
        job = self.resumable_job()
        (job.folder / 'result.json').unlink()
        results = runtime.saved_ffs_results(self.saved)
        self.assertEqual(len(results), 1)
        self.assertEqual(results[0]['status'], 'interrupted')
        self.assertTrue(results[0]['resumable'])
        self.assertIsNone(results[0]['result_path'])

    def test_only_kinetics_runs_are_offered_in_the_current_workflow(self):
        job = self.run_job('cancelled')
        metadata_path = job.folder / 'request.json'
        metadata = json.loads(metadata_path.read_text())
        metadata.pop('method')
        metadata['request'].pop('kinetics')
        metadata_path.write_text(json.dumps(metadata))
        self.assertEqual(runtime.saved_ffs_results(self.saved), [])
        self.assertIsNone(runtime.active_ffs(self.saved))

    def test_preparation_without_a_checkpoint_is_not_resumable(self):
        job = self.resumable_job()
        (job.folder / 'resume_state.json').unlink()
        self.assertFalse(runtime.saved_ffs_results(self.saved)[0]['resumable'])
        with self.assertRaises(ValueError):
            runtime.resume_ffs(self.saved, job.folder)

    def test_surviving_native_worker_is_adopted_without_starting_another(self):
        job = self.resumable_job()
        (job.folder / 'result.json').unlink()
        with patch.object(runtime, '_native_running', return_value=True), \
                patch.object(runtime.SavedSwitchFFSJob, 'start') as start:
            active = runtime.active_ffs(self.saved)
            resumed, created = runtime.resume_ffs(self.saved, job.folder)
            self.assertIs(resumed, active)
            self.assertFalse(created)
            self.assertEqual(active.folder, job.folder)
            self.assertEqual(active.snapshot().state, 'running')
            self.assertEqual(runtime.saved_ffs_results(self.saved), [])
            start.assert_not_called()
            active.cancel()
        self.assertTrue((job.folder / 'CANCEL').is_file())

    def test_reattached_progress_counts_completed_repetitions(self):
        job = self.resumable_job()
        (job.folder / 'progress.json').write_text(json.dumps({
            'status': 'running', 'current': 2, 'kinetics': {'replica': 2, 'replicas': 4}}))
        with patch.object(runtime, '_native_running', return_value=True):
            snapshot = runtime.active_ffs(self.saved).snapshot()
        self.assertEqual((snapshot.current, snapshot.total), (1, 4))

    def test_resume_refuses_foreign_folder_and_completed_run(self):
        job = self.resumable_job()
        result = json.loads((job.folder / 'result.json').read_text())
        (job.folder / 'result.json').write_text(json.dumps({**result, 'status': 'completed'}))
        with self.assertRaises(ValueError):
            runtime.resume_ffs(self.saved, job.folder)
        with self.assertRaises(ValueError):
            runtime.resume_ffs(self.saved, self.folder)

    def test_proximity_result_keeps_its_mode_and_suppresses_residual_rates(self):
        self.request['kinetics']['mode'] = 'proximity'
        job = runtime.SavedSwitchFFSJob(self.saved, self.intent)
        with patch.object(runtime, 'run_oxdna_ffs', return_value={'status': 'completed',
            'kinetics': {'rate': {'status': 'model_estimate', 'value_model': 12.}},
            'rate_estimate': {'value_model': 12.},
            'conditional_probability': {'status': 'conditional_estimate', 'value': .25, 'unit': 'probability'}}):
            job._run()
        result = job.snapshot().result
        self.assertEqual(result['simulation_mode'], 'proximity')
        self.assertEqual(result['kinetics']['mode'], 'proximity')
        self.assertNotIn('rate', result['kinetics'])
        self.assertIsNone(result['rate_estimate'])
        record = runtime.saved_ffs_results(self.saved)[0]
        self.assertEqual(record['simulation_mode'], 'proximity')
        self.assertEqual(record['conditional_probability']['value'], .25)

    def test_proximity_pause_resumes_frozen_mode_after_form_change(self):
        self.request['kinetics'].update(mode='proximity', source_samples=0, equilibration_steps=0,
                                        initial_distance_min=2., initial_distance_max=2.5)
        original = self.resumable_job()
        self.options['kinetics_mode'] = 'opening'
        self.prepare.reset_mock()
        with patch.object(runtime.SavedSwitchFFSJob, 'start'):
            resumed, started = runtime.resume_ffs(self.saved, original.folder)
        self.assertTrue(started)
        with patch.object(runtime, 'run_oxdna_ffs', return_value={'status': 'completed'}) as compute:
            resumed._run()
        self.prepare.assert_not_called()
        self.assertEqual(compute.call_args.args[0]['kinetics']['mode'], 'proximity')
        self.assertEqual(resumed.snapshot().result['simulation_mode'], 'proximity')
        self.assertEqual(resumed.run_id, original.run_id)

    def test_proximity_failure_before_preparation_preserves_chosen_mode(self):
        self.options['kinetics_mode'] = 'proximity'
        self.prepare.side_effect = ValueError('Preparation failed')
        job = runtime.SavedSwitchFFSJob(self.saved, self.intent)
        job._run()
        self.assertEqual(job.snapshot().result['simulation_mode'], 'proximity')
        record = runtime.saved_ffs_results(self.saved)[0]
        self.assertEqual(record['simulation_mode'], 'proximity')
        self.assertEqual(record['status'], 'failed')
        self.assertFalse(record['resumable'])

    def test_proximity_resume_rejects_invalid_frozen_initial_distance(self):
        self.request['kinetics'].update(mode='proximity', source_samples=0, equilibration_steps=0,
                                        initial_distance_min=2., initial_distance_max=2.5)
        original = self.resumable_job()
        path = original.folder / 'worker_request.json'
        request = json.loads(path.read_text())
        request['kinetics']['initial_distance_max'] = 8.
        path.write_text(json.dumps(request))
        record = runtime.saved_ffs_results(self.saved)[0]
        self.assertEqual(record['simulation_mode'], 'proximity')
        self.assertFalse(record['resumable'])
        with self.assertRaises(ValueError):
            runtime.resume_ffs(self.saved, original.folder)


if __name__ == '__main__':
    unittest.main()
