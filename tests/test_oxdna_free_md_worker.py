"""Free MD checkpoint integrity and native DNA2 integration smoke checks."""
from copy import deepcopy
import json
import os
from pathlib import Path
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'src'))
from aptaswitch_core.oxdna_free_md_worker import (
    DNAAdapter, Replica, _identity, _validate, guide_evidence, parse_frame,
    read_json, run_worker, worker_lock, write_json,
)


def frame_text(rows, step=0):
    return f't = {step}\nb = 40 40 40\nE = 0 0 0\n' + '\n'.join(
        ' '.join(f'{value:.17g}' for value in row) for row in rows) + '\n'


def request_fixture(root, production_steps=2500, sample_every=1500):
    params = root / 'dna_parameters.txt'
    params.write_text('fixture parameters\n')
    variants = []
    for name, trigger in [('original', 'CC'), ('extended', 'ACC')]:
        rows = np.zeros((2 + len(trigger), 15))
        rows[:, 3], rows[:, 8] = 1, 1
        rows[:, :3] = np.arange(len(rows))[:, None] + 10
        rows[:, 9:] = .1
        conf, top = root / f'{name}.dat', root / f'{name}.top'
        conf.write_text(frame_text(rows))
        top.write_text(f'{len(rows)} 2 5->3\nGG\n{trigger}\n')
        variants.append({'name': name, 'aptamer_sequence': 'GG', 'trigger_sequence': trigger,
                         'guide': '((+' + '.' * (len(trigger) - 2) + '))',
                         'guide_pairs': [[0, len(rows) - 1], [1, len(rows) - 2]],
                         'configuration_file': str(conf), 'topology_file': str(top), 'box_size': 40.})
    return {'format': 'aptaswitch-oxdna-free-md', 'version': 1, 'snapshot_id': 'fixture',
            'saved_switch_identity': 'test-saved-switch', 'runtime': {'python': sys.executable, 'seq_dep_file_DNA': str(params)},
            'variants': variants, 'monitored_bases': [1, 2], 'replicas': 2, 'workers': 1,
            'sampling': {'dt': .003, 'production_steps': production_steps, 'sample_every': sample_every,
                         'equilibration_steps': 20, 'temperature_c': 25., 'salt_m': .5, 'seed': 19,
                         'duration_us': production_steps * .003 * 3.031e-6}}


class FakeAdapter:
    records = []
    pause_at = None

    def __init__(self, request, variant, configuration, **kwargs):
        self.rows, self.box = parse_frame(Path(configuration).read_text())
        self.manager = SimpleNamespace(current_step=int(float(Path(configuration).read_text().splitlines()[0].split('=')[1])))
        self.variant, self.phase = variant, kwargs['phase']
        self.root = Path(kwargs['directory']).parent
        self.production = 0
        self.records.append(dict(kwargs, start_frame=self.frame()))

    def frame(self):
        return frame_text(self.rows, self.manager.current_step)

    def run(self, steps):
        self.manager.current_step += steps
        self.rows[:, 0] += steps / 1e6
        if self.phase == 'production':
            self.production += steps
            if self.pause_at is not None and self.production >= self.pause_at:
                (self.root / 'CANCEL').touch()

    def bonds(self):
        return {tuple(pair) for pair in self.variant['guide_pairs']}

    def checkpoint(self, target):
        Path(target).write_text(self.frame())

    def energy(self):
        return -.7

    def close(self):
        pass


class WorkerTests(unittest.TestCase):
    def setUp(self):
        tmp = tempfile.TemporaryDirectory(prefix='aptaswitch_free_md_')
        self.addCleanup(tmp.cleanup)
        self.root = Path(tmp.name)
        self.request = request_fixture(self.root)
        FakeAdapter.records, FakeAdapter.pause_at = [], None
        module = SimpleNamespace(frame_metrics=lambda *args: {'interstrand_pairs': len(args[2])},
                                 summarize_run=lambda *args: {})
        self.mock_module = patch.dict(sys.modules, {'aptaswitch_core.oxdna_free_md_metrics': module})
        self.mock_module.start()
        self.addCleanup(self.mock_module.stop)

    def test_balanced_replica_assignment_and_independent_seeds(self):
        request = deepcopy(self.request)
        request['replicas'] = 6
        replicas = [Replica(request, self.root / 'balance', i, adapter_factory=FakeAdapter) for i in range(6)]
        self.assertEqual([r.state['variant'] for r in replicas], ['original', 'extended'] * 3)
        self.assertEqual(len({r.state['seed'] for r in replicas}), 6)
        request['replicas'] = 3
        with self.assertRaisesRegex(ValueError, 'even'):
            _validate(request)

    def test_bound_complex_is_retained_and_restraints_are_only_preparatory(self):
        replica = Replica(self.request, self.root / 'bound', 0, adapter_factory=FakeAdapter)
        state = replica.run()
        self.assertEqual(state['status'], 'completed')
        self.assertEqual(state['samples_count'], 3)  # initial, 1500, final 2500
        self.assertEqual(state['production_steps'], 2500)
        self.assertTrue(all(stage['validation']['interstrand_pairs'] == 2 for stage in state['preparation']))
        phases = {item['phase']: item for item in FakeAdapter.records}
        for phase in ('release', 'equilibration', 'production'):
            self.assertIsNone(phases[phase]['forces'])
            self.assertFalse(phases[phase]['initial'])
            self.assertTrue(phases[phase]['binary'])
        data = [json.loads(line) for line in (replica.output / 'samples.jsonl').read_text().splitlines()]
        self.assertEqual([row['step'] for row in data], [0, 1500, 2500])
        self.assertAlmostEqual(data[-1]['time_us'], 2500 * .003 * 3.031e-6)
        self.assertEqual(data[-1]['hydrogen_bonds'], [[0, 3], [1, 2]])
        self.assertEqual((replica.output / 'trajectory.dat').read_text().count('t = '), 3)

    def test_pause_mid_sample_resumes_without_duplicate_samples_or_rethermalizing(self):
        FakeAdapter.pause_at = 1000
        replica = Replica(self.request, self.root / 'pause', 0, adapter_factory=FakeAdapter)
        first = deepcopy(replica.run())
        self.assertEqual(first['status'], 'cancelled')
        self.assertEqual(first['production_steps'], 1000)
        self.assertEqual(first['samples_count'], 1)
        saved = (replica.output / first['checkpoint_file']).read_text()
        # Emulate a crash after two noncommitted append operations.
        with (replica.output / 'samples.jsonl').open('a') as handle:
            handle.write('{"partial":')
        with (replica.output / 'trajectory.dat').open('a') as handle:
            handle.write('t = 900000\nb = ')
        FakeAdapter.pause_at = None
        (replica.root / 'CANCEL').unlink()
        resumed = Replica(self.request, replica.root, 0, adapter_factory=FakeAdapter)
        self.assertEqual((replica.output / 'samples.jsonl').stat().st_size, first['samples_bytes'])
        state = resumed.run()
        self.assertEqual(FakeAdapter.records[-1]['start_frame'], saved)
        self.assertFalse(FakeAdapter.records[-1]['initial'])
        self.assertTrue(FakeAdapter.records[-1]['binary'])
        data = [json.loads(line) for line in (replica.output / 'samples.jsonl').read_text().splitlines()]
        self.assertEqual([row['step'] for row in data], [0, 1500, 2500])
        self.assertEqual(state['samples_count'], 3)
        self.assertEqual(state['production_steps'], 2500)

    def test_checkpoint_identity_protects_coordinates_parameters_and_conditions(self):
        output = self.root / 'identity'
        Replica(self.request, output, 0, adapter_factory=FakeAdapter)
        workers = deepcopy(self.request)
        workers['workers'] = 2
        self.assertEqual(_identity(self.request), _identity(workers))
        changed = deepcopy(self.request)
        changed['monitored_bases'] = [2]
        with self.assertRaisesRegex(ValueError, 'different conditions'):
            Replica(changed, output, 0, adapter_factory=FakeAdapter)
        Path(self.request['variants'][0]['configuration_file']).write_text('changed seed bytes')
        with self.assertRaisesRegex(ValueError, 'different conditions'):
            Replica(self.request, output, 0, adapter_factory=FakeAdapter)

    def test_missing_committed_data_fails_instead_of_restarting(self):
        replica = Replica(self.request, self.root / 'truncated', 0, adapter_factory=FakeAdapter)
        replica.run()
        (replica.output / 'samples.jsonl').write_text('')
        with self.assertRaisesRegex(ValueError, 'Truncated'):
            Replica(self.request, replica.root, 0, adapter_factory=FakeAdapter)

    def test_lock_prevents_duplicate_supervisors_and_releases_after_exit(self):
        output = self.root / 'lock'
        output.mkdir()
        with worker_lock(output):
            with self.assertRaisesRegex(RuntimeError, 'already owns'):
                with worker_lock(output):
                    self.fail('duplicate supervisor acquired lock')
        with worker_lock(output):
            pass

    def test_guide_evidence_reports_missing_and_additional_pairs_without_rejecting_bound_complex(self):
        evidence = guide_evidence(self.request['variants'][0], {(0, 3), (0, 2)})
        self.assertEqual(evidence['interstrand_pairs'], 2)
        self.assertEqual(evidence['missing_pairs'], [[1, 2]])
        self.assertEqual(evidence['additional_pairs'], [[0, 2]])
        self.assertFalse(evidence['equilibrium_validated'])

    def test_periodic_self_guard_rejects_unsafe_geometry_before_retaining_a_sample(self):
        request = deepcopy(self.request)
        request['box_policy'] = {'geometry_monitor_required': True,
                                 'interaction_cutoff': 39., 'periodic_margin': 2.}
        replica = Replica(request, self.root / 'unsafe', 0, adapter_factory=FakeAdapter)
        with self.assertRaisesRegex(ValueError, 'Free MD stopped'):
            replica.run()
        self.assertFalse((replica.output / 'samples.jsonl').exists())
        self.assertEqual(read_json(replica.state_path)['samples_count'], 0)

    def test_master_resume_manifest_is_written_and_rejects_changed_conditions(self):
        (self.root / 'CANCEL').touch()
        report = run_worker(self.request, self.root)
        self.assertEqual(report['status'], 'cancelled')
        manifest = read_json(self.root / 'resume_state.json')
        self.assertEqual(manifest['identity'], _identity(self.request))
        changed = deepcopy(self.request)
        changed['sampling']['temperature_c'] += 1
        with self.assertRaisesRegex(ValueError, 'different conditions'):
            run_worker(changed, self.root)
        self.assertEqual(read_json(self.root / 'resume_state.json'), manifest)


@unittest.skipUnless(os.environ.get('APTA_FREE_MD_NATIVE_HOME'), 'requires external oxpy installation')
class NativeDNA2Tests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory(prefix='aptaswitch_free_md_native_')
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        self.request = request_fixture(self.root, production_steps=20, sample_every=10)
        native = Path(os.environ['APTA_FREE_MD_NATIVE_HOME'])
        self.request['runtime']['seq_dep_file_DNA'] = str(native / 'third_party/oxDNA/oxDNA2_sequence_dependent_parameters.txt')
        # Use an ideal antiparallel duplex. This fixture has no artificial loop
        # or nonphysical extension bonds, and all outputs remain in a tmpdir.
        from aptaswitch_core.oxdna_guided import _helix_pair
        pairs = [_helix_pair(i, 'dna') for i in range(2)]
        self.rows = np.asarray([pair[0] for pair in pairs] + [pair[1] for pair in reversed(pairs)])
        self.rows[:, :3] += 20
        self.variant = self.request['variants'][0]
        self.configuration = Path(self.variant['configuration_file'])
        self.configuration.write_text(frame_text(self.rows))

    def adapter(self, configuration, **kwargs):
        return DNAAdapter(self.request, self.variant, configuration, directory=self.root,
                          seed=2831, phase='production', **kwargs)

    def test_native_dna2_free_md_and_binary_restart_preserve_phase_point_and_rng(self):
        adapter = self.adapter(self.configuration, initial=True)
        checkpoint = self.root / 'restart.bin'
        try:
            self.assertEqual(len(adapter.bonds()), 2)
            adapter.run(120)
            adapter.checkpoint(checkpoint)
            saved, _ = parse_frame(adapter.frame())
            adapter.run(120)
            continuous, _ = parse_frame(adapter.frame())
            self.assertTrue(np.isfinite(adapter.energy()))
        finally:
            adapter.close()
        resumed = self.adapter(checkpoint, binary=True)
        try:
            reloaded, _ = parse_frame(resumed.frame())
            np.testing.assert_allclose(saved, reloaded, rtol=0, atol=1e-12)
            resumed.run(120)
            restarted, _ = parse_frame(resumed.frame())
            np.testing.assert_allclose(continuous, restarted, rtol=0, atol=2e-11)
            self.assertEqual(resumed.manager.current_step, 240)
            settings = (self.root / 'native_production_input.txt').read_text()
            self.assertIn('interaction_type = DNA2\n', settings)
            self.assertIn('external_forces = false\n', settings)
            self.assertIn('refresh_vel = false\n', settings)
            self.assertIn('reload_from = ', settings)
            self.assertNotIn('\nseed = ', settings)
        finally:
            resumed.close()

    def test_native_preparation_releases_traps_and_production_pause_resumes_mid_sample(self):
        request = deepcopy(self.request)
        request['sampling'].update(production_steps=2500, sample_every=1500)

        class PausingAdapter(DNAAdapter):
            def __init__(inner, *args, **kwargs):
                inner.phase = kwargs['phase']
                super().__init__(*args, **kwargs)

            def run(inner, steps):
                super().run(steps)
                if inner.phase == 'production':
                    (self.root / 'CANCEL').touch()

        replica = Replica(request, self.root, 0, adapter_factory=PausingAdapter)
        paused = deepcopy(replica.run())
        self.assertEqual(paused['status'], 'cancelled')
        self.assertEqual(paused['phase'], 'production')
        self.assertEqual(paused['production_steps'], 1000)
        self.assertEqual(paused['samples_count'], 1)
        self.assertEqual([row['phase'] for row in paused['preparation']],
                         ['minimize_guided', 'relax_guided', 'guide_physical', 'release', 'equilibration'])
        self.assertFalse(paused['preparation'][-1]['external_forces'])
        self.assertFalse(paused['equilibrium_validated'])
        reference = (replica.output / 'reference.dat').read_bytes()
        (self.root / 'CANCEL').unlink()
        resumed = Replica(request, self.root, 0)
        result = resumed.run()
        self.assertEqual(result['status'], 'completed')
        self.assertEqual(result['production_steps'], 2500)
        self.assertEqual((replica.output / 'reference.dat').read_bytes(), reference)
        samples = [json.loads(line) for line in (replica.output / 'samples.jsonl').read_text().splitlines()]
        self.assertEqual([sample['step'] for sample in samples], [0, 1500, 2500])
        self.assertTrue(all(np.isfinite(sample['minimum_trigger_site_distance_nm']) for sample in samples))
        settings = (replica.output / 'native_production_input.txt').read_text()
        self.assertIn('interaction_type = DNA2\n', settings)
        self.assertIn('external_forces = false\n', settings)
        self.assertIn('refresh_vel = false\n', settings)
        self.assertEqual((replica.output / 'trajectory.dat').read_text().count('t = '), 3)

    def test_two_native_subprocesses_produce_balanced_report_and_completed_resume_is_idle(self):
        from aptaswitch_core.oxdna_guided import _helix_pair
        request = deepcopy(self.request)
        request['workers'] = 2
        pairs = [_helix_pair(i, 'dna') for i in range(3)]
        rows = np.asarray([pairs[0][0], pairs[1][0], pairs[2][1], pairs[1][1], pairs[0][1]])
        rows[:, :3] += 20
        Path(request['variants'][1]['configuration_file']).write_text(frame_text(rows))
        report = run_worker(request, self.root)
        self.assertEqual(report['status'], 'completed')
        self.assertTrue(report['sampling_complete'])
        self.assertEqual(report['comparison']['groups']['original']['completed'], 1)
        self.assertEqual(report['comparison']['groups']['extended']['completed'], 1)
        self.assertEqual(read_json(self.root / 'progress.json')['fraction'], 1.)
        for index in range(2):
            directory = self.root / f'replica_{index:03d}'
            self.assertEqual(read_json(directory / 'state.json')['samples_count'], 3)
            self.assertTrue((directory / 'metrics.csv').is_file())
            self.assertTrue((directory / 'checkpoint.dat').is_file())
        with patch('aptaswitch_core.oxdna_free_md_worker.subprocess.Popen',
                   side_effect=AssertionError('Completed replicas must not be launched again')):
            self.assertEqual(run_worker(request, self.root)['status'], 'completed')

    @unittest.skipUnless(os.environ.get('APTA_FREE_MD_THERMODYNAMICS_LIBRARY'), 'requires integrated Vienna thermodynamics library')
    def test_real_native_thermodynamics_request_geometry_native_worker_and_high_level_resume(self):
        from aptaswitch_core.oxdna_free_md import (
            prepare_free_md_request, run_free_md, load_free_md_resume_request,
        )
        snapshot = {'snapshot_id': 'native-free-md-e2e', 'saved_at': '2026-10-14T00:00:00Z',
                    'request': {'thermo': {'trigger_material': 'dna', 'temperature_c': 25., 'sodium_m': .5}},
                    'candidate': {'id': 'native-test'},
                    'sequences': {'aptamer': 'GGGCGGGC', 'original_trigger': 'GCCCGCCC',
                                  'trigger': 'TTGCCCGCCCAA'}}
        request = prepare_free_md_request(snapshot, self.request['runtime'], {
            'monitored_bases': '3-6', 'duration_us': .001, 'replicas': 2, 'workers': 2,
            'seed': 371, 'thermodynamics_library_path': os.environ['APTA_FREE_MD_THERMODYNAMICS_LIBRARY']})
        # Only the test's internal request is shortened; production UI bounds
        # and user experiments remain unchanged. Run 100 actual MD steps/arm.
        request['sampling'].update(production_steps=100, sample_every=50, equilibration_steps=20,
                                   actual_duration_us=100 * .003 * 3.031e-6)
        output = self.root / 'end_to_end'
        paused = run_free_md(request, output, cancel_check=lambda: True)
        self.assertEqual(paused['status'], 'cancelled')
        self.assertFalse(paused['sampling_complete'])
        self.assertEqual(read_json(output / 'resume_state.json')['stage'], 'preparation')
        self.assertFalse(list(output.glob('replica_*/state.json')))
        pending = load_free_md_resume_request(output)
        self.assertTrue(pending['_resume_preparation'])
        self.assertEqual(pending['sampling']['seed'], request['sampling']['seed'])
        self.assertEqual(pending['variants'], request['variants'])
        report = run_free_md(pending, output)
        failures = str(report.get('error', '')) + str(report.get('failures', ''))
        if report.get('status') != 'completed':
            failures += '\n' + '\n'.join(path.read_text()[-4000:] for path in output.glob('*.log'))
        self.assertEqual(report['status'], 'completed', failures)
        self.assertTrue(report['sampling_complete'])
        resolved = report['request']
        self.assertNotIn('stage', read_json(output / 'resume_state.json'))
        self.assertEqual(read_json(output / 'resume_state.json')['identity'], _identity(resolved))
        self.assertNotEqual(resolved['guide_provenance']['version'], 'unknown')
        self.assertEqual(resolved['guide_source'], 'native_hybrid_mfe_two_strand_dna_complex')
        self.assertEqual(resolved['guide_provenance']['backend'], 'AptaSwitch-native-hybrid')
        self.assertTrue(all(v['guide_solver'] == 'AptaSwitch-connected-DP' for v in resolved['variants']))
        self.assertEqual(resolved['variants'][0]['box_size'], resolved['variants'][1]['box_size'])
        for index, variant in enumerate(resolved['variants']):
            aptamer, trigger = variant['aptamer_sequence'], variant['trigger_sequence']
            rows, box = parse_frame(Path(variant['configuration_file']).read_text())
            self.assertEqual(len(rows), len(aptamer) + len(trigger))
            self.assertEqual(Path(variant['topology_file']).read_text().splitlines(),
                             [f'{len(rows)} 2 5->3', f'{aptamer} type=DNA', f'{trigger} type=DNA'])
            self.assertTrue(any(i < len(aptamer) <= j for i, j in variant['guide_pairs']))
            self.assertEqual(variant['guide'].count('+'), 1)
            state = read_json(output / f'replica_{index:03d}' / 'state.json')
            self.assertEqual(state['production_steps'], 100)
            self.assertEqual(state['samples_count'], 3)
            self.assertIn('periodic_self_check', state)
            self.assertGreaterEqual(state['periodic_self_check']['minimum_clearance'],
                                    state['periodic_self_check']['required_clearance'])
            samples = [json.loads(line) for line in (output / f'replica_{index:03d}' / 'samples.jsonl').read_text().splitlines()]
            self.assertEqual([row['step'] for row in samples], [0, 50, 100])
            self.assertGreater(samples[0]['interstrand_pairs'], 0)
        trajectories = {path: path.read_bytes() for path in output.glob('replica_*/trajectory.dat')}
        resumed = run_free_md(load_free_md_resume_request(output), output)
        self.assertEqual(resumed['status'], 'completed', resumed.get('error'))
        self.assertEqual({path: path.read_bytes() for path in trajectories}, trajectories)


if __name__ == '__main__':
    unittest.main()
