"""Automatic preparation retries preserve physical inputs and native checks."""
import copy
import json
import shutil
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import numpy as np

from aptaswitch_core.oxdna_geometry import OxDNAPreparationCancelled
from aptaswitch_core.oxdna_guided import prepare_guided_system
from aptaswitch_core.oxdna_guided_worker import validation_from_bonds
from aptaswitch_core.scientific_process import CalculationCancelled
from aptaswitch_core.localization import language_context


class GuidedRetryTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        runtime = {'python': 'mock-native-python'}
        for chemistry in ('DNA', 'RNA', 'NA'):
            parameter = self.root/(chemistry+'.txt')
            parameter.write_text('immutable '+chemistry+' parameters\n')
            runtime['seq_dep_file_'+chemistry] = str(parameter)
        self.request = {
            'runtime': runtime, 'seed': 21, 'temperature_c': 37., 'salt_m': .12,
            'box_size': 100., 'initial_separation': 6., 'use_average_seq': False,
            'switch_sequence': 'AGGGGAAAACCCCA', 'switch_guide': '.((((....)))).',
            'switch_material': 'rna', 'trigger_sequence': 'ACGT'*7+'AC',
            'trigger_guide': '.'*30, 'trigger_material': 'dna',
            'plan': {'stem_pairs': [[1, 12], [2, 11], [3, 10], [4, 9]],
                     'basin_min_stem': 2, 'basin_distance': 4.},
            'guide_model': 'fixture-guide-model', 'guide_source': 'fixture-mfe',
            'guide_thermodynamics': {'parameter_digest': 'immutable-fixture-digest'},
        }
        self.output = self.root/'preparation'
        self.calls = []

    def worker_input(self, command):
        request = json.loads(Path(command[-1]).read_text())
        self.calls.append(request)
        return request

    def completed(self, request, *, missing=0, extra=(), strict=True):
        folder = Path(request['configuration_file']).parent
        expected = sorted(tuple(pair) for guide in request['guides'].values() for pair in guide['pairs'])
        bonds = set(expected[missing:]) | set(extra)
        validation = validation_from_bonds(request, bonds, 6.)
        validation['geometry'] = {'passed': True}
        validation['loop_threading'] = {'passed': True}
        configuration = folder/'validated.dat'
        shutil.copyfile(request['configuration_file'], configuration)
        result = {'status': 'completed', 'configuration_file': str(configuration),
                  'box_size': request['box_size'], 'validation': validation,
                  'strict_force_field_initialization_passed': strict,
                  'production_external_forces': False,
                  'guide_accepted': validation['guide_accepted'],
                  'off_basin_eligible': validation['off_basin_eligible'],
                  'ffs_eligible': validation['ffs_eligible']}
        (folder/'guided_result.json').write_text(json.dumps(result))
        return result

    def saved_report(self):
        return json.loads((self.output/'preparation.json').read_text())

    def test_rejected_then_accepted_uses_new_seed_and_internal_geometry(self):
        progress = []
        def worker(command, **kwargs):
            request = self.worker_input(command)
            self.completed(request, missing=1 if len(self.calls) == 1 else 0)

        with patch('aptaswitch_core.oxdna_guided.run_stage', side_effect=worker), language_context('fr'):
            result = prepare_guided_system(self.request, self.output, progress_callback=progress.append)
        self.assertTrue(result['ffs_eligible'])
        self.assertFalse(result['retry_exhausted'])
        self.assertEqual(result['selected_attempt'], 2)
        self.assertEqual([row['status'] for row in result['attempts']], ['rejected', 'accepted'])
        self.assertNotEqual(self.calls[0]['seed'], self.calls[1]['seed'])
        arrays = [np.loadtxt(row['configuration_file'], skiprows=3) for row in self.calls]
        internal = lambda rows: np.linalg.norm(rows[:14, None, :3]-rows[None, :14, :3], axis=-1)
        self.assertGreater(float(np.max(np.abs(internal(arrays[0])-internal(arrays[1])))), .01)
        self.assertEqual(self.saved_report()['selected_attempt'], 2)
        self.assertTrue(progress)
        self.assertTrue(all(message.startswith('Préparation 3D guidée · essai ') for message in progress))

    def test_exhausts_ten_attempts_and_retains_best_strict_diagnostic(self):
        def worker(command, **kwargs):
            request = self.worker_input(command)
            self.completed(request, missing=1 if len(self.calls) == 4 else 2)

        with patch('aptaswitch_core.oxdna_guided.run_stage', side_effect=worker):
            result = prepare_guided_system(self.request, self.output)
        self.assertEqual(len(self.calls), 10)
        self.assertEqual(len({row['seed'] for row in self.calls}), 10)
        self.assertTrue(result['retry_exhausted'])
        self.assertFalse(result['ffs_eligible'])
        self.assertEqual(result['selected_attempt'], 4)
        self.assertEqual(Path(result['configuration_file']).parent.name, 'attempt_04')
        self.assertEqual(result['validation']['guides']['switch']['missing_pairs'], [[1, 12]])
        self.assertEqual(len(list(self.output.glob('attempt_*/guided_result.json'))), 10)

    def test_exact_guide_cannot_override_a_failed_geometry_check(self):
        def worker(command, **kwargs):
            request = self.worker_input(command)
            result = self.completed(request)
            if len(self.calls) == 1:
                result['validation']['geometry']['passed'] = False
                path = Path(request['configuration_file']).parent/'guided_result.json'
                path.write_text(json.dumps(result))

        with patch('aptaswitch_core.oxdna_guided.run_stage', side_effect=worker):
            result = prepare_guided_system(self.request, self.output)
        self.assertEqual([row['status'] for row in result['attempts']], ['rejected', 'accepted'])
        self.assertEqual(result['selected_attempt'], 2)

    def test_threaded_seed_is_reconstructed_before_native_relaxation(self):
        def worker(command, **kwargs):
            self.completed(self.worker_input(command))

        from aptaswitch_core.oxdna_guided import configuration_loop_threading
        calls = 0
        def check(*args):
            nonlocal calls
            calls += 1
            if calls == 1:
                return {'passed': False, 'piercings': [{'segment': [2, 3]}], 'unresolved': []}
            return configuration_loop_threading(*args)

        with patch('aptaswitch_core.oxdna_guided.run_stage', side_effect=worker), \
             patch('aptaswitch_core.oxdna_guided.configuration_loop_threading', side_effect=check):
            result = prepare_guided_system(self.request, self.output)
        self.assertEqual(len(self.calls), 1)
        self.assertEqual(result['selected_attempt'], 2)
        self.assertEqual(result['attempts'][0]['failure_kind'], 'loop_threading')
        self.assertFalse((self.output/'attempt_01/guided_result.json').exists())

    def test_exact_pairs_and_valid_sterics_do_not_override_loop_threading(self):
        def worker(command, **kwargs):
            req = self.worker_input(command)
            result = self.completed(req)
            if len(self.calls) == 1:
                result['validation']['loop_threading']['passed'] = False
                (Path(req['configuration_file']).parent/'guided_result.json').write_text(json.dumps(result))
        with patch('aptaswitch_core.oxdna_guided.run_stage', side_effect=worker):
            result = prepare_guided_system(self.request, self.output)
        self.assertEqual(result['selected_attempt'], 2)
        self.assertFalse(result['attempts'][0]['validation']['loop_threading']['passed'])

    def test_attempt_timeout_then_success_is_recorded_without_losing_inputs(self):
        def worker(command, **kwargs):
            request = self.worker_input(command)
            if len(self.calls) == 1:
                raise TimeoutError('native preparation timeout')
            self.completed(request)

        with patch('aptaswitch_core.oxdna_guided.run_stage', side_effect=worker):
            result = prepare_guided_system(self.request, self.output)
        self.assertTrue(result['ffs_eligible'])
        self.assertEqual(result['attempts'][0]['failure_kind'], 'timeout')
        self.assertTrue((self.output/'attempt_01'/'system.dat').is_file())
        self.assertEqual(result['attempts'][1]['status'], 'accepted')

    def test_cancellation_between_retries_does_not_start_next_attempt(self):
        def worker(command, **kwargs):
            self.completed(self.worker_input(command), missing=1)

        def cancelled():
            path = self.output/'preparation.json'
            if not path.exists():
                return False
            rows = json.loads(path.read_text())['attempts']
            return bool(rows and rows[-1]['status'] == 'rejected')

        with patch('aptaswitch_core.oxdna_guided.run_stage', side_effect=worker):
            with self.assertRaises(OxDNAPreparationCancelled):
                prepare_guided_system(self.request, self.output, cancel_check=cancelled)
        report = self.saved_report()
        self.assertEqual(len(self.calls), 1)
        self.assertEqual(report['status'], 'cancelled')
        self.assertFalse(report['ffs_eligible'])
        self.assertNotIn('configuration_file', report)
        self.assertFalse((self.output/'attempt_02').exists())

    def test_cancelled_after_worker_completion_never_exposes_accepted_configuration(self):
        cancelled = False

        def worker(command, **kwargs):
            nonlocal cancelled
            self.completed(self.worker_input(command))
            cancelled = True

        with patch('aptaswitch_core.oxdna_guided.run_stage', side_effect=worker):
            with self.assertRaises(OxDNAPreparationCancelled):
                prepare_guided_system(self.request, self.output, cancel_check=lambda: cancelled)
        report = self.saved_report()
        self.assertEqual(report['status'], 'cancelled')
        self.assertEqual(report['attempts'][0]['status'], 'cancelled')
        self.assertFalse(report['ffs_eligible'])
        self.assertNotIn('configuration_file', report)

    def test_process_cancellation_has_no_retry_even_if_result_was_written(self):
        def worker(command, **kwargs):
            self.completed(self.worker_input(command))
            raise CalculationCancelled('cancelled at process completion')

        with patch('aptaswitch_core.oxdna_guided.run_stage', side_effect=worker):
            with self.assertRaises(OxDNAPreparationCancelled):
                prepare_guided_system(self.request, self.output)
        self.assertEqual(len(self.calls), 1)
        self.assertFalse(self.saved_report()['ffs_eligible'])

    def test_cancellation_wins_over_simultaneous_native_failure(self):
        cancelled = False

        def worker(command, **kwargs):
            nonlocal cancelled
            self.worker_input(command)
            cancelled = True
            raise RuntimeError('native process terminated at cancellation')

        with patch('aptaswitch_core.oxdna_guided.run_stage', side_effect=worker):
            with self.assertRaises(OxDNAPreparationCancelled):
                prepare_guided_system(self.request, self.output, cancel_check=lambda: cancelled)
        report = self.saved_report()
        self.assertEqual(len(self.calls), 1)
        self.assertEqual(report['status'], 'cancelled')
        self.assertEqual(report['attempts'][0]['status'], 'cancelled')
        self.assertFalse(report['ffs_eligible'])

    def test_missing_runtime_fails_once(self):
        def worker(command, **kwargs):
            self.worker_input(command)
            raise FileNotFoundError('missing native executable')

        with patch('aptaswitch_core.oxdna_guided.run_stage', side_effect=worker):
            with self.assertRaises(FileNotFoundError):
                prepare_guided_system(self.request, self.output)
        report = self.saved_report()
        self.assertEqual(len(self.calls), 1)
        self.assertEqual(report['attempts'][0]['failure_kind'], 'environment')
        self.assertEqual(report['status'], 'failed')

    def test_worker_environment_failure_is_not_retried(self):
        def worker(command, **kwargs):
            request = self.worker_input(command)
            folder = Path(request['configuration_file']).parent
            (folder/'guided_result.json').write_text(json.dumps(
                {'status': 'failed', 'failure_kind': 'environment', 'error': 'oxpy is unavailable'}))
            raise RuntimeError('worker failed')

        with patch('aptaswitch_core.oxdna_guided.run_stage', side_effect=worker):
            with self.assertRaisesRegex(RuntimeError, 'oxpy is unavailable'):
                prepare_guided_system(self.request, self.output)
        self.assertEqual(len(self.calls), 1)
        self.assertEqual(self.saved_report()['attempts'][0]['failure_kind'], 'environment')

    def test_contradictory_off_and_guide_targets_fail_before_native_attempts(self):
        self.request['plan']['stem_pairs'] = [[0, 13], [2, 10]]
        with patch('aptaswitch_core.oxdna_guided.run_stage') as worker:
            with self.assertRaises(ValueError):
                prepare_guided_system(self.request, self.output)
        worker.assert_not_called()
        self.assertEqual(self.saved_report()['attempts'], [])

    def test_retries_preserve_sequences_guides_conditions_and_parameter_contents(self):
        before = copy.deepcopy(self.request)

        def worker(command, **kwargs):
            request = self.worker_input(command)
            self.completed(request, missing=1 if len(self.calls) < 3 else 0)

        with patch('aptaswitch_core.oxdna_guided.run_stage', side_effect=worker):
            prepare_guided_system(self.request, self.output)
        self.assertEqual(self.request, before)
        keys = ('switch_sequence', 'trigger_sequence', 'switch_guide', 'trigger_guide',
                'switch_material', 'trigger_material', 'temperature_c', 'salt_m',
                'use_average_seq', 'plan', 'guide_model', 'guide_source', 'guide_thermodynamics')
        for attempt in self.calls:
            for key in keys:
                self.assertEqual(attempt[key], before[key])
            for chemistry in ('DNA', 'RNA', 'NA'):
                key = 'seq_dep_file_'+chemistry
                self.assertEqual(Path(attempt['runtime'][key]).read_bytes(), Path(before['runtime'][key]).read_bytes())

    def test_extra_native_pair_is_rejected_then_retried(self):
        def worker(command, **kwargs):
            self.completed(self.worker_input(command), extra=[(15, 20)] if len(self.calls) == 1 else [])

        with patch('aptaswitch_core.oxdna_guided.run_stage', side_effect=worker):
            report = prepare_guided_system(self.request, self.output)
        self.assertEqual(report['selected_attempt'], 2)
        self.assertFalse(report['attempts'][0]['validation']['guide_exact_match'])
        self.assertTrue(report['ffs_eligible'])

    def test_relaxed_only_result_cannot_become_accepted_or_a_strict_diagnostic(self):
        def worker(command, **kwargs):
            self.completed(self.worker_input(command), strict=False)

        with patch('aptaswitch_core.oxdna_guided.run_stage', side_effect=worker):
            with self.assertRaises(RuntimeError):
                prepare_guided_system(self.request, self.output)
        report = self.saved_report()
        self.assertEqual(len(self.calls), 10)
        self.assertFalse(report['ffs_eligible'])
        self.assertNotIn('configuration_file', report)


if __name__ == '__main__':
    unittest.main()
