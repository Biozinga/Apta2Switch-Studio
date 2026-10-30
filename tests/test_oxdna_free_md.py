"""Frozen comparison inputs, geometry, periodic metrics and export boundaries."""
from copy import deepcopy
import io
import json
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch
import zipfile

import numpy as np

from aptaswitch_core.oxdna_free_md import (available_free_md, free_md_sequences,
    parse_monitored_bases, prepare_free_md_request, prepare_free_md_system,
    load_free_md_resume_request, free_md_worker_running, make_free_md_archive, run_free_md)
from aptaswitch_core.oxdna_free_md_metrics import frame_metrics, summarize_run, LENGTH_NM
from aptaswitch_core.oxdna_free_md_worker import _identity, worker_lock
from aptaswitch_core.thermodynamics.types import ModelProvenance
from aptaswitch_core.thermodynamics.setup import MODEL_ID


def snapshot():
    return {'snapshot_id': 'saved', 'saved_at': 'date', 'candidate': {},
            'request': {'thermo': {'trigger_material': 'dna', 'temperature_c': 25., 'sodium_m': .5}},
            'sequences': {'aptamer': 'GCGAAACGC', 'original_trigger': 'GCGTTTCGC',
                          'trigger': 'AAGCGTTTCGC'}}


class FixtureBackend:
    provenance = ModelProvenance('AptaSwitch-native-hybrid', 'fixture', MODEL_ID,
                                 'sha256:fixture', 'connected-fixture', '1M')

    def minimum(self, spec, conditions):
        aptamer, trigger = spec.strands
        guide = '(' * len(aptamer.sequence) + '+' + '.' * (len(trigger.sequence) - len(aptamer.sequence)) + ')' * len(aptamer.sequence)
        return SimpleNamespace(structure=guide, energy_kcal_mol=-20.,
            context=SimpleNamespace(complex_spec=spec, conditions=conditions,
                                    provenance=self.provenance, solver='AptaSwitch-connected-DP'))


class FreeMDCoreTests(unittest.TestCase):
    def setUp(self):
        self.root = Path(self.enterContext(tempfile.TemporaryDirectory())).resolve()
        parameter = self.root / 'params.txt'
        parameter.write_text('test\n')
        self.runtime = {'python': '/test/python', 'seq_dep_file_DNA': str(parameter)}
        self.backend = FixtureBackend()
        self.enterContext(patch('aptaswitch_core.thermodynamics.setup.get_backend', return_value=self.backend))

    def request(self):
        return prepare_free_md_request(snapshot(), self.runtime,
            {'monitored_bases': [1, 3], 'duration_us': 1., 'thermodynamics_library_path': '/test/native.dylib'})

    def prepared(self):
        request = self.request()
        return prepare_free_md_system(request, self.root / 'run')

    def test_only_actual_dna_extension_available_and_sequences_frozen(self):
        source = snapshot()
        self.assertTrue(available_free_md(source))
        self.assertEqual(free_md_sequences(source)['extended'], 'AAGCGTTTCGC')
        for missing in ('aptamer', 'original_trigger', 'trigger'):
            invalid = deepcopy(source)
            invalid['sequences'].pop(missing)
            self.assertFalse(available_free_md(invalid))
        for trigger in ('GCGTTTCGC', 'GCGUUUCGC', 'AAAA', 'NNGCGTTTCGC'):
            invalid = deepcopy(source)
            invalid['sequences']['trigger'] = trigger
            self.assertFalse(available_free_md(invalid))
        source['request']['thermo']['trigger_material'] = 'rna'
        self.assertFalse(available_free_md(source))

    def test_indices_are_aptamer_positions_and_ranges_validated(self):
        self.assertEqual(parse_monitored_bases(' 1; 3 - 5, 4 9 ', 9), [1, 3, 4, 5, 9])
        for value in ('', '0', '10', '2-1', '1--3', '-1', '1.5', '1:a', '1-999999999'):
            with self.subTest(value=value), self.assertRaises(ValueError):
                parse_monitored_bases(value, 9)

    def test_balanced_allocation_model_time_and_independent_form_copy(self):
        request = self.request()
        self.assertEqual(request['sampling']['production_steps'], 109974706)
        self.assertEqual(request['replicas'], 10)
        self.assertEqual(request['model'], 'DNA2')
        self.assertFalse(request['production_external_forces'])
        self.runtime['python'] = 'changed'
        self.assertEqual(request['runtime']['python'], '/test/python')
        for options in ({'replicas': 3}, {'duration_us': float('nan')}, {'workers': 101},
                        {'monitored_bases': [True]}, {'monitored_bases': [10]}, {'seed': -1},
                        {'trigger_sequence': 'AAAA'}, {'temperature_c': 24.9}, {'temperature_c': 55.1},
                        {'salt_m': .119}, {'salt_m': 1.001}, {'import_path': '/legacy'}):
            with self.subTest(options=options), self.assertRaises(ValueError):
                prepare_free_md_request(snapshot(), self.runtime, {'monitored_bases': [1], **options})

    def test_complex_seed_contains_exactly_two_strands_and_no_dummy_linker(self):
        request = self.prepared()
        for variant in request['variants']:
            expected = len(variant['aptamer_sequence'] + variant['trigger_sequence'])
            top = Path(variant['topology_file']).read_text().splitlines()
            self.assertEqual(top[0], f'{expected} 2 5->3')
            self.assertEqual(len(Path(variant['configuration_file']).read_text().splitlines()) - 3, expected)
            self.assertEqual(len(variant['guide_pairs']), 9)
            self.assertTrue(all(i < 9 <= j for i, j in variant['guide_pairs']))
        self.assertEqual(request['variants'][0]['box_size'], request['variants'][1]['box_size'])
        self.assertEqual(len(request['prepared_file_sha256']), 5)

    def test_resume_reuses_inputs_and_rejects_parameter_or_condition_mutation(self):
        request = self.prepared()
        output = self.root / 'run'
        (output / 'worker_request.json').write_text(json.dumps(request))
        (output / 'resume_state.json').write_text(json.dumps({'identity': _identity(request)}))
        resumed = load_free_md_resume_request(output)
        self.assertTrue(resumed.pop('resume'))
        self.assertEqual(resumed, request)
        with worker_lock(output):
            self.assertTrue(free_md_worker_running(output))
            with self.assertRaises(ValueError):
                load_free_md_resume_request(output)
        self.assertFalse(free_md_worker_running(output))
        changed = deepcopy(request)
        changed['sampling']['temperature_c'] = 30
        (output / 'worker_request.json').write_text(json.dumps(changed))
        with self.assertRaises(ValueError):
            load_free_md_resume_request(output)
        (output / 'worker_request.json').write_text(json.dumps(request))
        Path(request['runtime']['seq_dep_file_DNA']).write_text('changed')
        with self.assertRaises(ValueError):
            load_free_md_resume_request(output)

    def test_preparation_refuses_changed_model_or_conditions_before_folding(self):
        request = self.request()
        mutated = deepcopy(request)
        mutated['guide_thermodynamics']['parameter_digest'] = 'sha256:changed'
        with patch.object(self.backend, 'minimum') as minimum:
            with self.assertRaisesRegex(ValueError, 'changé'):
                prepare_free_md_system(mutated, self.root / 'changed_model')
            minimum.assert_not_called()
        mutated = deepcopy(request)
        mutated['sampling']['temperature_c'] = 30.
        with self.assertRaisesRegex(ValueError, 'incompatible'):
            prepare_free_md_system(mutated, self.root / 'changed_conditions')

    def test_resume_rejects_legacy_folding_model_without_relabelling(self):
        request = self.request()
        request.pop('guide_thermodynamics')
        request['guide_model'] = 'dna'
        output = self.root / 'legacy_pause'
        report = run_free_md(request, output, cancel_check=lambda: True)
        self.assertEqual(report['status'], 'cancelled')
        with self.assertRaisesRegex(ValueError, 'incompatible'):
            load_free_md_resume_request(output)

    def test_archive_does_not_follow_symlinks_or_include_control_files(self):
        run = self.root / 'archive'
        run.mkdir()
        (run / 'samples.jsonl').write_text('{}\n')
        (run / 'outside').symlink_to(self.root / 'params.txt')
        (run / 'worker.lock').touch()
        (run / 'CANCEL').touch()
        with zipfile.ZipFile(io.BytesIO(make_free_md_archive(run))) as archive:
            self.assertEqual(archive.namelist(), ['samples.jsonl'])

    def test_pause_before_native_preparation_freezes_options_and_remains_resumable(self):
        request = self.request()
        output = self.root / 'early_pause'
        result = run_free_md(request, output, cancel_check=lambda: True)
        self.assertEqual(result['status'], 'cancelled')
        self.assertNotIn('error', result)
        self.assertFalse((output / 'worker_request.json').exists())
        resumed = load_free_md_resume_request(output)
        self.assertTrue(resumed.pop('resume'))
        self.assertTrue(resumed.pop('_resume_preparation'))
        self.assertEqual(resumed, request)
        old_seed = resumed['sampling']['seed']
        request['sampling']['seed'] = 42
        self.assertEqual(load_free_md_resume_request(output)['sampling']['seed'], old_seed)
        # A second pause preserves the original request and first-attempt artifacts.
        again = run_free_md(load_free_md_resume_request(output), output, cancel_check=lambda: True)
        self.assertEqual(again['status'], 'cancelled')
        self.assertTrue((output / 'preparation').is_dir())
        self.assertTrue((output / 'preparation_002').is_dir())
        (output / 'preparation_request.json').write_text('{}')
        with self.assertRaises(ValueError):
            load_free_md_resume_request(output)


class FreeMDMetricsTests(unittest.TestCase):
    def setUp(self):
        self.rows = np.zeros((6, 15))
        self.rows[:, :3] = [[9.5, 1, 1], [.2, 1, 1], [.8, 1.5, 1], [1.2, 1.7, 1.6],
                            [9.5, 1, 2], [.2, 1, 2]]
        self.box = np.array([10., 10., 10.])

    def test_minimum_periodic_distance_pair_occupancy_and_unwrapped_rg(self):
        metrics = frame_metrics(self.rows, self.box, [(0, 4), (4, 0), (1, 2)], 4,
                                [1, 2], self.rows[:4, :3], 1.)
        self.assertAlmostEqual(metrics['minimum_trigger_site_distance_nm'], LENGTH_NM)
        self.assertEqual(metrics['site_contact_fraction'], 1)
        self.assertEqual(metrics['site_paired_fraction'], 1)
        self.assertEqual(metrics['site_trigger_paired_fraction'], .5)
        self.assertEqual(metrics['interstrand_pairs'], 1)
        self.assertLess(metrics['rmsd_aptamer_nm'], 1e-12)
        unwrapped = self.rows[:4, :3].copy()
        unwrapped[1:, 0] += 10
        expected = np.sqrt(np.mean(np.sum((unwrapped - unwrapped.mean(axis=0)) ** 2, axis=1))) * LENGTH_NM
        self.assertAlmostEqual(metrics['rg_aptamer_nm'], expected)

    def test_rigid_rotation_translation_and_wrapping_do_not_produce_rmsd(self):
        rows = self.rows.copy()
        rows[:, :3] = (rows[:, :3] @ np.array([[0, -1, 0], [1, 0, 0], [0, 0, 1]]) + [4, 7, 2]) % 10
        metrics = frame_metrics(rows, self.box, [], 4, [1, 4], self.rows[:4, :3])
        self.assertLess(metrics['rmsd_aptamer_nm'], 1e-12)
        self.assertLess(metrics['rmsd_site_nm'], 1e-12)

    def test_partial_comparison_counts_replicas_and_does_not_extend_short_trajectories(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            topology = root / 'system.top'
            topology.write_text('6 2 5->3\nACGT\nAC\n')
            variants = [{'name': name, 'aptamer_sequence': 'ACGT', 'trigger_sequence': 'AC',
                         'topology_file': str(topology)} for name in ('original', 'extended')]
            request = {'replicas': 4, 'variants': variants, 'monitored_bases': [1, 2],
                       'sampling': {'dt': .003, 'production_steps': 10}}
            for index, steps in enumerate(([0, 10], [0, 10], [0], [])):
                directory = root / f'replica_{index:03d}'
                directory.mkdir()
                state = {'variant': variants[index % 2]['name'], 'samples_count': len(steps),
                         'production_steps': steps[-1] if steps else 0,
                         'status': 'completed' if len(steps) == 2 else 'cancelled'}
                (directory / 'state.json').write_text(json.dumps(state))
                samples = []
                for step in steps:
                    sample = frame_metrics(self.rows, self.box, [], 4, [1, 2], self.rows[:4, :3])
                    sample.update(step=step, time_us=step * .003 * 3.031e-6,
                                  variant=state['variant'], replica_index=index)
                    samples.append(sample)
                # A partial append must not be treated as a committed observation.
                (directory / 'samples.jsonl').write_text(''.join(json.dumps(s) + '\n' for s in samples) + '{broken')
            result = summarize_run(root, request)
            self.assertFalse(result['sampling_complete'])
            comparison = result['comparison']
            self.assertEqual(comparison['status'], 'partial')
            self.assertEqual(comparison['groups']['original']['available'], 2)
            curve = [row for row in comparison['curves'] if row['variant'] == 'original'
                     and row['metric'] == 'site_contact_fraction']
            self.assertEqual([row['n'] for row in curve], [2, 1])
            self.assertIsNone(curve[-1]['sd'])
            self.assertEqual(result['replicas'][0]['site'][0]['rmsf_nm'], 0)
            self.assertIsNone(result['replicas'][2]['site'][0]['rmsf_nm'])
            self.assertTrue((root / 'comparison_curves.csv').is_file())
            self.assertTrue((root / 'replica_000' / 'metrics.csv').is_file())


if __name__ == '__main__':
    unittest.main()
