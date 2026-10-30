"""Native monomer burn-in precedes artificial proximity and resumes exactly."""
from copy import deepcopy
import json
import random
from pathlib import Path
from types import SimpleNamespace
import tempfile
import unittest
from unittest.mock import patch
import numpy as np

from aptaswitch_core.oxdna_ffs_worker import Observation, SamplingStopped
from aptaswitch_core.oxdna_kinetics import (PROXIMITY_EQUILIBRATION_METHOD, source_fingerprint,
                                         validate_kinetic_protocol, proximity_ffs_report)
from aptaswitch_core.oxdna_kinetics_worker import FFSReplica, IsolatedMonomerAdapter, KineticSampler, read_json
from aptaswitch_core.oxdna_initialization import random_proximity_configuration
from test_oxdna_kinetics_worker import request
from test_oxdna_kinetics import stage


class MonomerScript:
    def __init__(self):
        self.starts, self.chunks = [], []
        self.after_chunk = None
        self.frame_guard = None

    def __call__(self, request, configuration, seed, initial, *, topology_file, directory):
        frame = Path(configuration).read_text()
        role = Path(topology_file).stem
        self.starts.append({'role': role, 'initial': initial, 'frame': frame,
                            'topology': Path(topology_file).read_text()})
        parent = self
        initial_frame = frame

        class Manager:
            current = initial_frame

            def frame(self):
                if parent.frame_guard:
                    parent.frame_guard(self.current)
                return self.current

            def advance(self, steps):
                lines = self.current.splitlines()
                lines[0] = 't = ' + str(int(lines[0].split('=')[1]) + steps)
                rows = []
                for index, row in enumerate(lines[3:]):
                    values = list(map(float, row.split()))
                    values[0] += .001 * steps * (index + 1)
                    values[9:12] = [.1, .2, .3]
                    values[12:15] = [.4, .5, .6]
                    rows.append(' '.join(map(str, values)))
                self.current = '\n'.join([*lines[:3], *rows]) + '\n'
                parent.chunks.append((role, steps))
                if parent.after_chunk:
                    parent.after_chunk()
                return steps

            def close(self):
                pass

        return Manager()


class SeedValidator:
    def __init__(self):
        self.starts = []

    def __call__(self, request, configuration, seed, initial, **kwargs):
        frame = Path(configuration).read_text()
        self.starts.append({'frame': frame, 'initial': initial})

        class Adapter:
            def observe(self):
                return Observation(0, 1, 1, 2.25, 0, 0, 2.25)

            def frame(self):
                return frame

            def close(self):
                pass

        return Adapter()


class EquilibratedProximityTests(unittest.TestCase):
    def setUp(self):
        self.root = Path(self.enterContext(tempfile.TemporaryDirectory())).resolve()
        self.req = request(self.root)
        self.req.update(switch_sequence='AU', trigger_sequence='AT', switch_material='rna', trigger_material='dna')
        self.req['kinetics'].update(version=4, mode='proximity', source_samples=0, equilibration_steps=5,
            equilibration_method=PROXIMITY_EQUILIBRATION_METHOD, initial_distance_min=2., initial_distance_max=2.5,
            chunk_steps=2)
        self.req['plan'].update(distance_scope='trigger_to_toehold', distance_pairs=[[0, 2], [0, 3]],
            separation_pairs=[[0, 2], [0, 3], [1, 2], [1, 3]], rate_estimation=False)
        self.monomers, self.validator = MonomerScript(), SeedValidator()

    def replica(self):
        return FFSReplica(self.req, self.root / 'replica_001', 0,
                          adapter_factory=self.validator, monomer_adapter_factory=self.monomers)

    def test_new_protocol_freezes_positive_isolated_burn_in_and_preserves_legacy(self):
        validate_kinetic_protocol(self.req['kinetics'])
        for changes in ({'equilibration_steps': 0}, {'equilibration_steps': True},
                        {'equilibration_method': 'joint_pair_md'}, {'source_samples': 75}):
            with self.subTest(changes=changes), self.assertRaises(ValueError):
                validate_kinetic_protocol(self.req['kinetics'] | changes)
        legacy = self.req['kinetics'] | {'version': 3, 'equilibration_steps': 0}
        validate_kinetic_protocol(legacy)
        self.assertNotEqual(source_fingerprint(self.req), source_fingerprint(self.req | {'kinetics': legacy}))
        changed = deepcopy(self.req)
        changed['kinetics']['equilibration_steps'] += 1
        self.assertNotEqual(source_fingerprint(self.req), source_fingerprint(changed))

    def test_both_monomers_equilibrate_without_partner_then_only_close_source_is_retained(self):
        replica = self.replica()
        self.assertEqual(self.validator.starts, [])
        replica.prepare_proximity()
        self.assertEqual(self.monomers.chunks, [('switch', 2), ('switch', 2), ('switch', 1),
                                                ('trigger', 2), ('trigger', 2), ('trigger', 1)])
        self.assertEqual([row['topology'] for row in self.monomers.starts],
                         ['2 1 5->3\nAU type=RNA\n', '2 1 5->3\nAT type=DNA\n'])
        self.assertEqual(replica.state['burn_steps'], 10)
        self.assertEqual(replica.state['steps'], 10)
        self.assertEqual(replica.state['flux']['exposure_steps'], 0)
        self.assertEqual(replica.state['flux']['bank'], [])
        self.assertEqual(len(replica.state['banks']), 1)
        self.assertEqual(replica.state['phase'], 'production')
        self.assertFalse(self.validator.starts[0]['initial'])
        self.assertTrue(replica.state['initialization']['native_unbound_verified'])
        self.assertEqual(replica.summary()['monomer_equilibration']['completed_steps'], 10)
        self.assertFalse(replica.summary()['monomer_equilibration']['equilibrium_validated'])
        seed = replica.state['banks'][replica.state['initial_source']]
        initial = (replica.output / seed['trajectory_file']).read_text()
        self.assertTrue(initial.startswith('t = 0\n'))
        self.assertEqual(initial, self.validator.starts[0]['frame'])
        self.assertEqual(initial.count('t ='), 1)
        # L is body-frame in oxDNA, unlike r, orientation axes and velocity.
        rows = np.asarray([list(map(float, row.split())) for row in initial.splitlines()[3:]])
        for index, part in enumerate(replica.state['monomer_equilibration']['monomers']):
            original = np.asarray([list(map(float, row.split())) for row in part['frame'].splitlines()[3:]])
            placed = rows[index*2:(index+1)*2]
            rotation = np.asarray(replica.state['initialization']['strand_rotations'][index])
            np.testing.assert_allclose(placed[:, :3] - placed[0, :3],
                                       (original[:, :3] - original[0, :3]) @ rotation.T,
                                       rtol=0, atol=1e-12)
            for begin in (3, 6, 9):
                np.testing.assert_allclose(placed[:, begin:begin+3], original[:, begin:begin+3] @ rotation.T,
                                           rtol=0, atol=1e-12)
            np.testing.assert_allclose(placed[:, 12:15], original[:, 12:15], rtol=0, atol=0)
            for old, new in zip(original, placed):
                old_axes = np.column_stack((old[3:6], np.cross(old[6:9], old[3:6]), old[6:9]))
                new_axes = np.column_stack((new[3:6], np.cross(new[6:9], new[3:6]), new[6:9]))
                np.testing.assert_allclose(new_axes @ new[12:15], rotation @ (old_axes @ old[12:15]),
                                           rtol=0, atol=1e-12)
        legacy = deepcopy(self.req)
        legacy['kinetics'].update(version=3, equilibration_steps=0)
        legacy['prepared']['configuration_file'] = str(replica.output / 'equilibrated_monomers.dat')
        legacy_frame, legacy_metadata = random_proximity_configuration(legacy, random.Random(53))
        legacy_rows = np.asarray([list(map(float, row.split())) for row in legacy_frame.splitlines()[3:]])
        for index, part in enumerate(replica.state['monomer_equilibration']['monomers']):
            original = np.asarray([list(map(float, row.split())) for row in part['frame'].splitlines()[3:]])
            rotation = np.asarray(legacy_metadata['strand_rotations'][index])
            np.testing.assert_allclose(legacy_rows[index*2:(index+1)*2, 12:15], original[:, 12:15] @ rotation.T,
                                       rtol=0, atol=1e-12)

    def test_pause_in_first_monomer_retains_current_full_frame_and_resumes_remaining_steps(self):
        replica = self.replica()
        self.monomers.after_chunk = lambda: (self.root / 'CANCEL').touch()
        with self.assertRaises(SamplingStopped):
            replica.prepare_proximity()
        before = read_json(replica.output / 'state.json')
        monomer = before['monomer_equilibration']['monomers'][0]
        self.assertEqual(monomer['steps'], 2)
        self.assertEqual(before['phase'], 'equilibration_switch')
        self.assertEqual(len(before['frame'].splitlines()[3:]), 4)
        self.assertEqual((replica.output / 'live.dat').read_text(), before['frame'])
        self.assertEqual(self.validator.starts, [])
        self.assertEqual(before['banks'], {})
        (self.root / 'CANCEL').unlink()
        self.monomers = MonomerScript()
        with patch('aptaswitch_core.oxdna_kinetics_worker.random_separated_configuration',
                   side_effect=AssertionError('Never replace a saved monomer phase point')):
            resumed = self.replica()
        resumed.prepare_proximity()
        self.assertEqual(self.monomers.starts[0]['frame'], monomer['frame'])
        self.assertFalse(self.monomers.starts[0]['initial'])
        self.assertEqual(self.monomers.chunks, [('switch', 2), ('switch', 1),
                                                ('trigger', 2), ('trigger', 2), ('trigger', 1)])
        self.assertEqual(resumed.state['burn_steps'], 10)
        calls = (len(self.validator.starts), len(self.monomers.starts))
        resumed.prepare_proximity()
        self.assertEqual(calls, (len(self.validator.starts), len(self.monomers.starts)))
        changed = deepcopy(self.req)
        changed['kinetics']['equilibration_steps'] = 6
        with self.assertRaisesRegex(ValueError, 'different simulation'):
            FFSReplica(changed, resumed.output, 0, adapter_factory=self.validator,
                       monomer_adapter_factory=self.monomers)

    def test_progress_and_conditional_report_record_preparation_without_source_flux(self):
        replica = self.replica()
        replica.prepare_proximity()
        replica.state.update(completed=True, stages=[stage(50)])
        controller = KineticSampler(self.req, self.root)
        progress = controller.aggregate({0: replica.summary()}, 0)
        self.assertEqual(progress['monomer_equilibration']['completed_steps'], 10)
        self.assertEqual(progress['monomer_equilibration']['total_steps'], 40)
        self.assertNotIn('flux_crossings', progress['kinetics'])
        report = proximity_ffs_report([replica.state])
        self.assertEqual(report['unit'], 'probability')
        self.assertEqual(report['initial_condition']['equilibration_steps_per_monomer'], 5)
        self.assertEqual(report['initial_condition']['ensemble'], 'isolated_monomer_md_then_unbound_proximity_v2')
        self.assertNotIn('value_model', report)
        legacy = {'replica_id': 'legacy', 'completed': True, 'initial_source': 'seed', 'stages': [stage(50)]}
        with self.assertRaisesRegex(ValueError, 'cannot be averaged'):
            proximity_ffs_report([replica.state, legacy])

    def test_periodic_self_image_guard_rejects_an_extended_monomer_before_any_source_or_movie(self):
        frame = 't = 2\nb = 20 20 20\nE = 0 0 0\n' + ''.join(
            f'{i} 0 0 1 0 0 0 0 1 0 0 0 .4 .5 .6\n' for i in range(16))
        adapter = object.__new__(IsolatedMonomerAdapter)
        adapter.manager = SimpleNamespace(current_step=2, update_observable_data=lambda: None)
        adapter.observable = SimpleNamespace(get_output_string=lambda step: frame)
        adapter.nucleotide_count = 16
        adapter.box_policy = {'interaction_cutoff': 5.}
        with self.assertRaisesRegex(ValueError, 'periodic image'):
            adapter.frame()
        replica = self.replica()
        self.monomers.after_chunk = lambda: setattr(self.monomers, 'frame_guard', lambda frame: adapter.frame())
        result = replica.run()
        self.assertEqual(result['status'], 'failed')
        self.assertIn('periodic image', result['error'])
        self.assertIsNone(result['initial_source'])
        self.assertIsNone(result['movie'])
        self.assertEqual(result['banks'], {})
        report = proximity_ffs_report([result])
        self.assertEqual(report['status'], 'incomplete')
        self.assertIsNone(report['value'])
