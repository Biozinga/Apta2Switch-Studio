"""Opt-in short real-NA check; original prepared calculations remain read only."""
from copy import deepcopy
import json
import os
from pathlib import Path
import tempfile
import unittest
import numpy as np

from aptaswitch_core.oxdna_kinetics import PROXIMITY_EQUILIBRATION_METHOD
from aptaswitch_core.oxdna_kinetics_worker import FFSReplica


@unittest.skipUnless(os.environ.get('APTA_NATIVE_FFS_REQUEST'), 'requires oxpy and a read-only prepared request')
class NativeMonomerEquilibrationTests(unittest.TestCase):
    def test_isolated_native_monomers_move_before_unbound_close_source_with_no_production(self):
        original = Path(os.environ['APTA_NATIVE_FFS_REQUEST'])
        before = original.read_bytes()
        req = deepcopy(json.loads(before))
        req['kinetics'].update(version=4, mode='proximity', source_samples=0, equilibration_steps=5,
            equilibration_method=PROXIMITY_EQUILIBRATION_METHOD,
            initial_distance_min=2., initial_distance_max=2.5, chunk_steps=5)
        n, total = req['plan']['switch_length'], req['plan']['nucleotide_count']
        exposed = sorted({pair[0] for pair in req['plan']['toehold_pairs']})
        req['plan'].update(rate_estimation=False, distance_scope='trigger_to_toehold',
                           distance_pairs=[[i, j] for i in exposed for j in range(n, total)])
        with tempfile.TemporaryDirectory(prefix='aptaswitch_proximity4_native_') as temporary:
            output = Path(temporary) / 'replica_001'
            replica = FFSReplica(req, output, 0)
            initial = [part['frame'] for part in replica.state['monomer_equilibration']['monomers']]
            replica.prepare_proximity()
            eq = replica.summary()['monomer_equilibration']
            self.assertEqual(eq['completed_steps'], 10)
            self.assertEqual(eq['switch_steps'], 5)
            self.assertEqual(eq['trigger_steps'], 5)
            self.assertTrue(eq['completed'])
            self.assertEqual(replica.state['flux']['exposure_steps'], 0)
            self.assertEqual(replica.state['stages'], [])
            self.assertEqual(replica.state['observation']['associations'], 0)
            self.assertFalse(replica.state['initial_velocities'])
            self.assertEqual(replica.state['phase'], 'production')
            self.assertEqual(len(replica.state['banks']), 1)
            placed_rows = np.asarray([list(map(float, row.split())) for row in replica.state['frame'].splitlines()[3:]])
            offset = 0
            for index, part in enumerate(replica.state['monomer_equilibration']['monomers']):
                old = [row.split()[:3] for row in initial[index].splitlines()[3:]]
                new = [row.split()[:3] for row in part['frame'].splitlines()[3:]]
                self.assertNotEqual(old, new)
                native_rows = np.asarray([list(map(float, row.split())) for row in part['frame'].splitlines()[3:]])
                np.testing.assert_allclose(placed_rows[offset:offset+len(native_rows), 12:15], native_rows[:, 12:15],
                                           rtol=0, atol=1e-10)
                offset += len(native_rows)
                text = (output / 'monomer_equilibration' / part['role'] / 'native_input.txt').read_text()
                self.assertIn('sim_type = MD\n', text)
                self.assertIn('interaction_type = NA\n', text)
                self.assertIn('external_forces = false\n', text)
                self.assertTrue((output / 'monomer_equilibration' / (part['role'] + '.top')).read_text().splitlines()[0].endswith(' 1 5->3'))
        self.assertEqual(original.read_bytes(), before)
