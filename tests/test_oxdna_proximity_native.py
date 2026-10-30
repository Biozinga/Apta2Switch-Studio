"""Opt-in real oxNA checks of the artificial, unpaired proximity ensemble.

Run with the external oxpy Python and APTA_NATIVE_FFS_REQUEST pointing to a
prepared worker request. Original run inputs are read only; all native logs,
configurations, banks and checkpoints are written into temporary directories.
These short checks validate initialization and restart, not opening kinetics.
"""
from copy import deepcopy
import json
import math
import os
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'src'))

from aptaswitch_core.oxdna_ffs_worker import _same_checkpoint, validate_request
from aptaswitch_core.oxdna_kinetics_worker import FFSReplica
from aptaswitch_core.oxdna_native_ffs import NativeFFSAdapter


def particle_rows(frame):
    return [list(map(float, row.split())) for row in frame.splitlines()[3:]]


@unittest.skipUnless(os.environ.get('APTA_NATIVE_FFS_REQUEST'),
                     'requires external oxpy and a prepared request')
class ActualNativeProximityTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix='aptaswitch_proximity_native_')
        self.addCleanup(self.temporary.cleanup)
        self.directory = Path(self.temporary.name)
        request = deepcopy(json.loads(Path(os.environ['APTA_NATIVE_FFS_REQUEST']).read_text()))
        request['kinetics'] = {
            'version': 3, 'mode': 'proximity', 'source_samples': 0,
            'trials': 3, 'replicas': 10, 'workers': 1,
            'equilibration_steps': 0, 'chunk_steps': 20,
            'initial_distance_min': 2.0, 'initial_distance_max': 2.5,
        }
        request['plan']['rate_estimation'] = False
        # This is a NEW artificial proximity experiment using only the frozen
        # monomer inputs. Older full-kinetic requests keep their terminal-pair
        # coordinate on resume; never rewrite their on-disk plan.
        plan = request['plan']
        if plan.get('distance_scope') == 'trigger_3prime_to_toehold_start':
            n, total = plan['switch_length'], plan['nucleotide_count']
            exposed = sorted({pair[0] for pair in plan['toehold_pairs']})
            plan.update(distance_scope='trigger_to_toehold',
                        distance_pairs=[[i, j] for i in exposed for j in range(n, total)])
        self.request = validate_request(request)

    def test_ten_distinct_unpaired_native_seeds_without_motion_or_restraints(self):
        placements = set()
        for index in range(10):
            with self.subTest(replica=index):
                replica = FFSReplica(self.request, self.directory / f'replica_{index:03d}', index)
                drawn = replica.state['frame']
                initial_step = int(float(drawn.splitlines()[0].split('=', 1)[1]))
                replica.prepare_proximity()
                state = replica.state
                observation = state['observation']
                placements.add(tuple(tuple(row[:3]) for row in particle_rows(state['frame'])))
                self.assertEqual(observation['associations'], 0)
                self.assertEqual(observation['binding'], 0)
                self.assertGreaterEqual(observation['distance'], 2.0 - 1e-9)
                self.assertLessEqual(observation['distance'], 2.5 + 1e-9)
                self.assertGreaterEqual(observation['separation'], 2.0 - 1e-9)
                self.assertEqual(observation['step'], initial_step)
                self.assertEqual(state['steps'], 0)
                self.assertEqual(state['burn_steps'], 0)
                self.assertEqual(state['flux']['exposure_steps'], 0)
                self.assertEqual(state['flux']['bank'], [])
                self.assertFalse(state['initial_velocities'])
                self.assertEqual(state['phase'], 'production')
                self.assertEqual(len(state['production_bank']), 1)

                # Native initialization may refresh momenta and normalize axes;
                # it must not move any nucleotide or deform the prepared strands.
                before, after = particle_rows(drawn), particle_rows(state['frame'])
                self.assertEqual(len(before), len(after))
                for original, native in zip(before, after):
                    self.assertEqual(len(native), 15)
                    self.assertTrue(all(math.isfinite(value) for value in native))
                    for first, second in zip(original[:9], native[:9]):
                        self.assertAlmostEqual(first, second, places=9)
                self.assertTrue(any(abs(value) > 0 for row in after for value in row[9:]))
                inputs = list(replica.output.rglob('native_input.txt'))
                self.assertTrue(inputs)
                for source in inputs:
                    self.assertIn('external_forces = false\n', source.read_text())

        self.assertEqual(len(placements), 10)

    def test_initialized_seed_resumes_without_redraw_or_velocity_refresh(self):
        output = self.directory / 'resume'
        first = FFSReplica(self.request, output, 0)
        first.prepare_proximity()
        first.save(force=True)
        frame = first.state['frame']
        rng = first.rng.getstate()
        bank = deepcopy(first.state['banks'])
        bank_files = {path.name: path.read_bytes() for path in (output / 'banks').iterdir()}

        with patch('aptaswitch_core.oxdna_kinetics_worker.random_proximity_configuration',
                   side_effect=AssertionError('A saved initial state must not be redrawn')):
            resumed = FFSReplica(self.request, output, 0,
                adapter_factory=lambda *args, **kwargs: self.fail('A ready seed must not be reinitialized'))
            resumed.prepare_proximity()
        self.assertEqual(resumed.state['frame'], frame)
        self.assertEqual(resumed.rng.getstate(), rng)
        self.assertEqual(resumed.state['banks'], bank)
        self.assertEqual({path.name: path.read_bytes() for path in (output / 'banks').iterdir()}, bank_files)
        self.assertFalse(resumed.state['initial_velocities'])

        configuration = self.directory / 'restart.dat'
        configuration.write_text(frame)
        adapter = NativeFFSAdapter(resumed.request, configuration, 371, False,
            directory=self.directory / 'native_restart', conditions=[])
        try:
            self.assertTrue(_same_checkpoint(frame, adapter.frame()))
            self.assertEqual(adapter.observe().associations, 0)
        finally:
            adapter.close()

    def test_close_seed_advances_twenty_native_steps_with_finite_phase_point(self):
        replica = FFSReplica(self.request, self.directory / 'short_md', 0)
        replica.prepare_proximity()
        configuration = self.directory / 'short_md.dat'
        configuration.write_text(replica.state['frame'])
        adapter = NativeFFSAdapter(replica.request, configuration, 733, False,
            directory=self.directory / 'native_short_md', conditions=[])
        try:
            before = adapter.observe()
            result = adapter.advance(20)
            self.assertEqual(result.steps, 20)
            self.assertIsNone(result.condition)
            self.assertEqual(result.observation.step, before.step + 20)
            self.assertTrue(math.isfinite(adapter.manager.system_energy()))
            self.assertTrue(all(math.isfinite(value) for row in particle_rows(adapter.frame()) for value in row))
        finally:
            adapter.close()


if __name__ == '__main__':
    unittest.main()
