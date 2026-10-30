"""Boundary controls and durable resume with the analytical backend contract."""
from copy import deepcopy
from dataclasses import asdict, replace
import json
import os
from pathlib import Path
import random
import tempfile
import unittest
from unittest.mock import patch

from aptaswitch_core.design_checkpoint import CHECKPOINT_NAME, save_design_checkpoint
from aptaswitch_core.engines import run_design
from aptaswitch_core.models import (InverseExtensionSettings, SequenceInput, RunRequest,
                                   ReporterContext, SwitchArchitecture, ThermoConditions)
from aptaswitch_core.thermodynamics.optimization import ScoredAssignment, TrialOptions, optimize
from aptaswitch_core.thermodynamics.design_constraints import UnsatisfiableDesignError
from aptaswitch_core.workflow import run_and_export
from aptaswitch_studio.run_library import load_run
from tests.test_thermodynamics_design_engine import AnalyticalBackend, request as design_request


class DesignCheckpointTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.directory = Path(temporary.name).resolve()
        self.request = replace(design_request(), trials=3, output_dir=self.directory/'calculations',
            optimization_options=asdict(TrialOptions(seed=31, max_evaluations=1, max_proposals=1)))
        self.backend = self.make_backend(self.request)
        self.enterContext(patch('aptaswitch_core.thermodynamics.setup.backend_for_request',
                                side_effect=lambda request: self.backend))

    def make_backend(self, request):
        backend = AnalyticalBackend(request)
        backend.provenance = replace(backend.provenance, model=request.thermo.material)
        return backend

    def pause_after_first_trial(self, *, request=None, snapshots=None):
        action = None

        def checkpoint(result):
            nonlocal action
            if snapshots is not None:
                snapshots.append(deepcopy(result))
            if len(result.optimization['trials']) == 1:
                action = 'pause'

        result, exports = run_and_export(request or self.request,
            control_check=lambda: action, checkpoint_callback=checkpoint)
        return result, exports

    def test_initial_checkpoint_exists_before_first_thermodynamic_call(self):
        seen = []
        analyze = self.backend.analyze

        def check_then_analyze(*args, **kwargs):
            files = list(self.directory.glob('calculations/designs/*/'+CHECKPOINT_NAME))
            self.assertEqual(len(files), 1)
            checkpoint = json.loads(files[0].read_text())
            if not seen:
                self.assertEqual(checkpoint['results'], [])
                self.assertEqual(checkpoint['optimization']['checkpoint']['completed_trials'], 0)
                self.assertEqual(checkpoint['status'], 'running')
                self.assertEqual(load_run(files[0]).result.run_id, checkpoint['run_id'])
            seen.append(True)
            return analyze(*args, **kwargs)

        self.backend.analyze = check_then_analyze
        self.pause_after_first_trial()
        self.assertTrue(seen)

    def test_pause_finishes_current_trial_even_when_requested_during_evaluation(self):
        action = None
        analyze = self.backend.analyze

        def request_pause(*args, **kwargs):
            nonlocal action
            action = 'pause'
            return analyze(*args, **kwargs)

        self.backend.analyze = request_pause
        result, exports = run_and_export(self.request, control_check=lambda: action)
        self.assertEqual(result.status, 'paused')
        self.assertEqual([candidate.trial for candidate in result.candidates], [1])
        self.assertEqual([report['trial'] for report in result.optimization['trials']], [1])
        self.assertEqual(set(exports), {'checkpoint'})
        self.assertFalse(list(result.output_dir.glob('*_results.json')))

    def optimizer_with_two_finalists(self, seeds):
        """Keep the real optimizer result, with a second valid screening finalist."""
        def two_finalists(problem, loss, **kwargs):
            report = optimize(problem, loss, **kwargs)
            seeds.append(report.options.seed)
            other = problem.generate(random.Random(report.options.seed+100))
            self.assertNotEqual(other, report.best_assignment)
            return replace(report, elites=(ScoredAssignment(other, loss(other)),))
        return two_finalists

    def test_cancel_mid_screening_repeats_unfinished_trial_from_disk_with_same_seed(self):
        request = replace(self.request, trials=1)
        seeds = []
        cancelled = False
        screen = self.backend.screen_unpaired_region
        self.backend.bad_on = True

        def reject_then_cancel(*args, **kwargs):
            nonlocal cancelled
            candidate = screen(*args, **kwargs)
            cancelled = True
            return candidate

        self.backend.screen_unpaired_region = reject_then_cancel
        with patch('aptaswitch_core.thermodynamics.design_engine.optimize',
                   side_effect=self.optimizer_with_two_finalists(seeds)):
            result = run_design(request, cancel_check=lambda: cancelled)
            self.assertEqual(result.status, 'cancelled')
            self.assertEqual(result.candidates, [])
            self.assertEqual(result.optimization['trials'], [])
            self.assertEqual(result.optimization['checkpoint']['completed_trials'], 0)
            interrupted = result.optimization['interrupted_trial']
            self.assertEqual(interrupted['trial'], 1)
            self.assertEqual(interrupted['stage'], 'screening')
            self.assertEqual(interrupted['candidate_screened'], 1)
            self.assertEqual(interrupted['stop_reason'], 'budget_exhausted')
            self.assertIsNone(interrupted['selected_assignment'])
            loaded = load_run(result.output_dir/CHECKPOINT_NAME).result
            self.backend.bad_on = False
            self.backend.screen_unpaired_region = screen
            resumed = run_design(loaded.request, resume_result=loaded)
        self.assertEqual(seeds, [31, 31])
        self.assertEqual(resumed.run_id, result.run_id)
        self.assertEqual(resumed.status, 'completed')
        self.assertEqual([report['trial'] for report in resumed.optimization['trials']], [1])
        self.assertEqual([candidate.trial for candidate in resumed.candidates], [1])
        self.assertFalse(resumed.candidates[0].excluded)
        self.assertNotIn('interrupted_trial', resumed.optimization)

    def test_cancel_after_admission_or_exhausted_finalists_keeps_completed_trial(self):
        for all_excluded, cancel_at in ((False, 1), (True, 2)):
            with self.subTest(all_excluded=all_excluded):
                request = replace(self.request, trials=1)
                self.backend = self.make_backend(request)
                self.backend.bad_on = all_excluded
                checked = 0
                screen = self.backend.screen_unpaired_region

                def screen_then_cancel(*args, **kwargs):
                    nonlocal checked
                    checked += 1
                    return screen(*args, **kwargs)

                self.backend.screen_unpaired_region = screen_then_cancel
                with patch('aptaswitch_core.thermodynamics.design_engine.optimize',
                           side_effect=self.optimizer_with_two_finalists([])):
                    result = run_design(request, cancel_check=lambda: checked >= cancel_at)
                self.assertEqual(result.status, 'cancelled')
                self.assertEqual(result.optimization['checkpoint']['completed_trials'], 1)
                self.assertEqual(result.optimization['trials'][0]['candidate_screened'], cancel_at)
                self.assertEqual(len(result.candidates), 1)
                self.assertEqual(result.candidates[0].excluded, all_excluded)
                self.assertNotIn('interrupted_trial', result.optimization)

    def test_disk_resume_preserves_identity_seeds_and_matches_uninterrupted_run(self):
        snapshots = []
        paused, paths = self.pause_after_first_trial(snapshots=snapshots)
        initial_report = deepcopy(paused.optimization['trials'][0])
        loaded = load_run(Path(paths['checkpoint'])).result
        self.assertEqual(loaded.request.output_dir.name, 'designs')
        progress = []
        resumed, exports = run_and_export(loaded.request, resume_result=loaded,
                                         progress_callback=lambda *event: progress.append(event))
        self.assertEqual(progress[0][0], 1)
        self.assertEqual(resumed.run_id, paused.run_id)
        self.assertEqual(resumed.output_dir, paused.output_dir)
        self.assertEqual(resumed.started_at, paused.started_at)
        self.assertEqual(resumed.status, 'completed')
        self.assertEqual(resumed.optimization['trials'][0], initial_report)
        self.assertEqual([report['options']['seed'] for report in resumed.optimization['trials']], [31, 32, 33])
        self.assertEqual(len({candidate.trial for candidate in resumed.candidates}), 3)
        self.assertEqual(resumed.optimization['checkpoint']['phase'], 'done')
        baseline, _ = run_and_export(replace(self.request, output_dir=self.directory/'baseline'))
        self.assertEqual([candidate.to_dict() for candidate in resumed.candidates],
                         [candidate.to_dict() for candidate in baseline.candidates])
        self.assertEqual(resumed.optimization['trials'], baseline.optimization['trials'])
        self.assertEqual(load_run(Path(exports['json'])).result.to_dict(), resumed.to_dict())
        self.assertEqual([len(result.optimization['trials']) for result in snapshots], [0, 1, 1, 1])

    def test_category_output_directory_and_rebuilt_library_path_resume(self):
        request = replace(self.request, output_dir=self.directory/'other'/'designs',
                          thermodynamics_library_path='/old/app/library')
        paused, paths = self.pause_after_first_trial(request=request)
        loaded = load_run(Path(paths['checkpoint'])).result
        rebuilt = replace(loaded.request, thermodynamics_library_path='/rebuilt/app/library')
        result, _ = run_and_export(rebuilt, resume_result=loaded)
        self.assertEqual(result.status, 'completed')
        self.assertEqual(result.output_dir, paused.output_dir)
        self.assertEqual(result.request.thermodynamics_library_path, '/rebuilt/app/library')

    def test_finish_saves_partial_ranked_exports_and_does_not_run_next_trial(self):
        action = None

        def checkpoint(result):
            nonlocal action
            if len(result.optimization['trials']) == 1:
                action = 'finish'

        result, exports = run_and_export(self.request, control_check=lambda: action,
                                         checkpoint_callback=checkpoint)
        self.assertEqual(result.status, 'stopped')
        self.assertEqual(len(result.candidates), 1)
        self.assertEqual(len(result.optimization['trials']), 1)
        self.assertEqual(result.optimization['checkpoint']['phase'], 'done')
        self.assertTrue(Path(exports['json']).is_file())
        self.assertTrue(Path(exports['csv']).is_file())
        self.assertEqual(load_run(Path(exports['json'])).result.status, 'stopped')

    def test_stop_before_any_trial_still_has_valid_empty_exports(self):
        result, exports = run_and_export(self.request, control_check=lambda: 'finish')
        self.assertEqual(result.status, 'stopped')
        self.assertEqual(result.candidates, [])
        self.assertEqual(self.backend.calls, [])
        self.assertEqual(load_run(Path(exports['json'])).result.status, 'stopped')

    def test_changed_request_or_native_provenance_is_rejected_without_new_folds(self):
        paused, paths = self.pause_after_first_trial()
        source = load_run(Path(paths['checkpoint'])).result
        before = len(self.backend.calls)
        with self.assertRaisesRegex(ValueError, 'settings differ'):
            run_and_export(replace(source.request, trials=4), resume_result=source)
        self.backend.provenance = replace(self.backend.provenance, version='different-native-build')
        with self.assertRaisesRegex(ValueError, 'engine or parameters differ'):
            run_and_export(source.request, resume_result=source)
        self.assertEqual(len(self.backend.calls), before)
        self.assertEqual(json.loads(Path(paths['checkpoint']).read_text())['status'], 'paused')

    def test_stale_or_modified_checkpoint_is_not_allowed_to_overwrite_latest(self):
        snapshots = []
        paused, _ = self.pause_after_first_trial(snapshots=snapshots)
        with self.assertRaisesRegex(ValueError, 'newer or different'):
            run_design(snapshots[0].request, resume_result=snapshots[0])
        changed = deepcopy(paused)
        changed.candidates[0].switch_seq = 'A'+changed.candidates[0].switch_seq[1:]
        with self.assertRaisesRegex(ValueError, 'integrity'):
            run_design(changed.request, resume_result=changed)

    def test_duplicate_trial_records_fail_even_with_resealed_file(self):
        paused, _ = self.pause_after_first_trial()
        paused.optimization['trials'].append(deepcopy(paused.optimization['trials'][0]))
        save_design_checkpoint(paused)
        with self.assertRaisesRegex(ValueError, 'trial sequence'):
            run_design(paused.request, resume_result=paused)

    def inverse_request(self, *, trials=3):
        request = replace(self.request, sequences=SequenceInput('AACCGGTTAACC', 'ACGT'*6),
                          inverse_extension=InverseExtensionSettings(30, 'both'), trials=trials)
        self.backend = self.make_backend(request)
        return request

    def test_inverse_pause_after_design_resumes_comparison_without_redesign(self):
        request = self.inverse_request(trials=1)
        action = None

        def checkpoint(result):
            nonlocal action
            if result.optimization['checkpoint']['phase'] == 'inverse_extension':
                action = 'pause'

        def analyze(result, **kwargs):
            result.inverse_extension_analysis = {'status': 'completed', 'test_comparison': True}

        with patch('aptaswitch_core.inverse_extension.analyze_inverse_extensions', side_effect=analyze) as comparison:
            paused, paths = run_and_export(request, control_check=lambda: action,
                                           checkpoint_callback=checkpoint)
            self.assertEqual(paused.status, 'paused')
            self.assertEqual(paused.optimization['checkpoint']['phase'], 'inverse_extension')
            comparison.assert_not_called()
            loaded = load_run(Path(paths['checkpoint'])).result
            with patch('aptaswitch_core.engines.design_candidates', side_effect=AssertionError('must not redesign')):
                result, exports = run_and_export(loaded.request, resume_result=loaded)
            comparison.assert_called_once()
        self.assertEqual(result.status, 'completed')
        self.assertTrue(result.inverse_extension_analysis['test_comparison'])
        self.assertTrue(Path(exports['inverse_csv']).is_file())

    def test_finish_partial_inverse_run_still_performs_comparison(self):
        request = self.inverse_request()
        action = None

        def checkpoint(result):
            nonlocal action
            if len(result.optimization['trials']) == 1:
                action = 'finish'

        def analyze(result, **kwargs):
            self.assertEqual(len(result.candidates), 1)
            result.inverse_extension_analysis = {'status': 'completed'}

        with patch('aptaswitch_core.inverse_extension.analyze_inverse_extensions', side_effect=analyze) as comparison:
            result, exports = run_and_export(request, control_check=lambda: action,
                                             checkpoint_callback=checkpoint)
        comparison.assert_called_once()
        self.assertEqual(result.status, 'stopped')
        self.assertEqual(result.inverse_extension_analysis['status'], 'completed')
        self.assertTrue(Path(exports['inverse_csv']).is_file())

    def test_export_failure_can_resume_exports_without_repeating_design(self):
        with patch('aptaswitch_core.workflow.export_run', side_effect=OSError('disk temporarily unavailable')):
            with self.assertRaises(OSError):
                run_and_export(self.request)
        path, = self.directory.glob('calculations/designs/*/'+CHECKPOINT_NAME)
        loaded = load_run(path).result
        self.assertEqual(loaded.optimization['checkpoint']['phase'], 'exports')
        with patch('aptaswitch_core.engines.design_candidates', side_effect=AssertionError('must not redesign')):
            result, exports = run_and_export(loaded.request, resume_result=loaded)
        self.assertEqual(result.status, 'completed')
        self.assertEqual(result.optimization['checkpoint']['phase'], 'done')
        self.assertTrue(Path(exports['json']).is_file())

    def test_inverse_analysis_recorded_failure_is_resumable_without_redesign(self):
        request = self.inverse_request(trials=1)
        # Exercise the real analyzer's exception capture, rather than an
        # analyzer mock raising directly past the workflow's status check.
        with patch('aptaswitch_core.inverse_extension.open_extension_model',
                   side_effect=RuntimeError('temporary comparison failure')), \
             patch('aptaswitch_core.inverse_extension.extension_model_settings', return_value={'model': 'fixture'}), \
             patch('aptaswitch_core.workflow.export_run') as export:
            with self.assertRaisesRegex(RuntimeError, 'temporary comparison failure'):
                run_and_export(request)
            export.assert_not_called()
        path, = self.directory.glob('calculations/designs/*/'+CHECKPOINT_NAME)
        loaded = load_run(path).result
        self.assertEqual(loaded.status, 'failed')
        self.assertEqual(loaded.inverse_extension_analysis['status'], 'failed')
        self.assertEqual(loaded.optimization['checkpoint']['phase'], 'inverse_extension')
        original_reports = deepcopy(loaded.optimization['trials'])

        def analyze(result, **kwargs):
            result.inverse_extension_analysis = {'status': 'completed'}

        with patch('aptaswitch_core.engines.design_candidates', side_effect=AssertionError('must not redesign')), \
             patch('aptaswitch_core.inverse_extension.analyze_inverse_extensions', side_effect=analyze) as comparison:
            result, exports = run_and_export(loaded.request, resume_result=loaded)
        comparison.assert_called_once()
        self.assertEqual(result.run_id, loaded.run_id)
        self.assertEqual(result.status, 'completed')
        self.assertEqual(result.optimization['trials'], original_reports)
        self.assertEqual(result.optimization['checkpoint']['phase'], 'done')
        self.assertTrue(Path(exports['json']).is_file())

    def test_impossible_stop_does_not_record_a_completed_checkpoint_trial(self):
        request = replace(self.request, exclude_stop_candidates=True,
                          architecture=replace(self.request.architecture, frame_linker='UAA'))
        self.backend = self.make_backend(request)
        snapshots = []
        with patch('aptaswitch_core.thermodynamics.design_engine.optimize',
                   side_effect=AssertionError('An impossible fixed STOP must fail before trials')):
            with self.assertRaisesRegex(UnsatisfiableDesignError, 'STOP'):
                self.pause_after_first_trial(request=request, snapshots=snapshots)
        self.assertEqual(self.backend.calls, [])
        for result in snapshots:
            self.assertEqual(result.candidates, [])
            self.assertEqual(result.optimization['checkpoint']['completed_trials'], 0)
            self.assertEqual(result.optimization['trials'], [])
        self.assertFalse(list(self.directory.rglob('*_results.json')))


@unittest.skipUnless(os.environ.get('APTASWITCH_HYBRID_NATIVE_LIBRARY'), 'Native library not configured')
class NativeCheckpointRoundTripTests(unittest.TestCase):
    def test_native_pause_disk_reload_and_resume_equals_uninterrupted_design(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory).resolve()
            request = RunRequest(SequenceInput('GGTT', 'AACC'),
                SwitchArchitecture(2, 2, upper_stem_length=1, upper_stem2_length=0,
                                   rbs_loop='AGGA', t7_leader='', rbs_prefix_length=0),
                reporter=ReporterContext(key='short', label='Test reporter', sequence_after_start_rna='GCU'),
                thermo=ThermoConditions(temperature_c=37., validation_max_complex_size=2),
                trials=2, output_dir=root/'calculations', exclude_stop_candidates=False,
                exclude_non_linear_candidates=False,
                thermodynamics_library_path=os.environ['APTASWITCH_HYBRID_NATIVE_LIBRARY'],
                optimization_options=asdict(TrialOptions(seed=19, max_evaluations=1, max_proposals=1)))
            action = None

            def checkpoint(result):
                nonlocal action
                if len(result.optimization['trials']) == 1:
                    action = 'pause'

            paused, paths = run_and_export(request, control_check=lambda: action,
                                           checkpoint_callback=checkpoint)
            self.assertEqual(paused.status, 'paused')
            loaded = load_run(Path(paths['checkpoint'])).result
            resumed, exports = run_and_export(loaded.request, resume_result=loaded)
            baseline, _ = run_and_export(replace(request, output_dir=root/'baseline'))
            self.assertEqual(resumed.run_id, paused.run_id)
            self.assertEqual(resumed.optimization['trials'], baseline.optimization['trials'])
            self.assertEqual([candidate.to_dict() for candidate in resumed.candidates],
                             [candidate.to_dict() for candidate in baseline.candidates])
            self.assertEqual(load_run(Path(exports['json'])).result.to_dict(), resumed.to_dict())


if __name__ == '__main__':
    unittest.main()
