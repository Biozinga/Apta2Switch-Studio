"""Separated source placement preserves monomers and never randomizes a resume."""
import random
from types import SimpleNamespace
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import numpy as np
from aptaswitch_core.oxdna_guided import seed_strand, _configuration

from aptaswitch_core.oxdna_initialization import random_separated_configuration, random_proximity_configuration
from aptaswitch_core.oxdna_kinetics_worker import FFSReplica
from test_oxdna_kinetics_worker import NativeScript, request


def coordinates(frame):
    return np.asarray([[float(x) for x in row.split()] for row in frame.splitlines()[3:]])


class IndependentPlacementTests(unittest.TestCase):
    def setUp(self):
        self.root = Path(self.enterContext(tempfile.TemporaryDirectory()))
        self.request = request(self.root)
        self.request['plan'].update(switch_length=3, nucleotide_count=5)
        self.request['switch_guide'] = '...'
        # A covalent chain crosses the x box boundary in the prepared input.
        self.rows = np.asarray([[x, y, 3., 0, 0, 1, 1, 0, 0, 1, 2, 3, 3, 2, 1]
                                for x, y in ((99.4, 2.), (.1, 2.2), (.8, 2.6), (4., 6.), (4.7, 6.2))])
        self.source = self.root/'initial.dat'
        self.source.write_text('t = 0\nb = 100 100 100\nE = 0 0 0\n'+
                               '\n'.join(' '.join(map(str, row)) for row in self.rows)+'\n')

    def test_rotation_preserves_chirality_internal_geometry_and_momenta(self):
        frozen = self.source.read_bytes()
        frame, provenance = random_separated_configuration(self.request, random.Random(37))
        result = coordinates(frame)
        for original, transformed, rotation in zip((self.rows[:3], self.rows[3:]),
                                                   (result[:3], result[3:]),
                                                   provenance['strand_rotations']):
            delta = original[:, None, :3] - original[None, :, :3]
            delta -= 100 * np.floor(delta / 100 + .5)
            after = transformed[:, None, :3] - transformed[None, :, :3]
            np.testing.assert_allclose(np.linalg.norm(delta, axis=2), np.linalg.norm(after, axis=2), atol=1e-12)
            matrix = np.asarray(rotation)
            self.assertAlmostEqual(np.linalg.det(matrix), 1.)
            for offset in range(3, 15, 3):
                np.testing.assert_allclose(transformed[:, offset:offset+3], original[:, offset:offset+3] @ matrix.T)
        self.assertEqual(self.source.read_bytes(), frozen)
        self.assertFalse(provenance['internal_conformations_changed'])
        self.assertFalse(provenance['equilibrium_sample'])

    def test_ten_reproducible_independent_far_placements_check_all_periodic_pairs(self):
        frames = []
        for index in range(10):
            frame, metadata = random_separated_configuration(self.request, random.Random(43+104729*(index+1)))
            again = random_separated_configuration(self.request, random.Random(43+104729*(index+1)))
            self.assertEqual((frame, metadata), again)
            result = coordinates(frame)
            delta = result[:3, None, :3] - result[None, 3:, :3]
            delta -= 100 * np.floor(delta / 100 + .5)
            self.assertGreaterEqual(float(np.linalg.norm(delta, axis=2).min()), 10.)
            self.assertGreaterEqual(metadata['periodic_self_check']['minimum_clearance'],
                                    metadata['periodic_self_check']['required_clearance'])
            self.assertTrue(metadata['geometry_check']['passed'])
            self.assertTrue(metadata['loop_threading_check']['passed'])
            self.assertEqual(metadata['geometry_check']['repulsive_core_contacts'], 0)
            frames.append(frame)
        self.assertEqual(len(set(frames)), 10)

    def test_impossible_box_rejected_before_random_draw_and_attempts_bounded(self):
        self.request['initial_separation'] = 90.
        with patch('aptaswitch_core.oxdna_initialization._rotation', side_effect=AssertionError('should not draw')):
            with self.assertRaisesRegex(ValueError, 'box is too small'):
                random_separated_configuration(self.request, random.Random(3))
        self.request['initial_separation'] = 10.
        with patch('aptaswitch_core.oxdna_initialization.check_periodic_self_clearance',
                   side_effect=ValueError('periodic collision')) as guard:
            with self.assertRaisesRegex(ValueError, 'after 3 attempts'):
                random_separated_configuration(self.request, random.Random(3), max_attempts=3)
            self.assertEqual(guard.call_count, 3)

    def test_invalid_prepared_monomer_rejected_before_random_placement(self):
        lines = self.source.read_text().splitlines()
        lines[4] = lines[3]  # coincident covalent bases, not repairable by rigid rotation
        self.source.write_text('\n'.join(lines)+'\n')
        with patch('aptaswitch_core.oxdna_initialization._rotation', side_effect=AssertionError('no draw')):
            with self.assertRaisesRegex(ValueError, 'prepared monomer'):
                random_separated_configuration(self.request, random.Random(3))

    def test_threaded_source_is_not_randomly_rotated_into_ten_replicas(self):
        with patch('aptaswitch_core.oxdna_initialization.configuration_loop_threading',
                   return_value={'passed': False, 'piercings': [{'segment': [1, 2]}]}), \
             patch('aptaswitch_core.oxdna_initialization._rotation', side_effect=AssertionError('no draw')):
            with self.assertRaisesRegex(ValueError, 'threads a guide loop'):
                random_separated_configuration(self.request, random.Random(3))

    def test_resume_keeps_placement_and_rng_without_any_new_draw(self):
        req = request(self.root)
        folder = self.root/'replica_001'
        original = FFSReplica(req, folder, 0, adapter_factory=NativeScript([]))
        original.state['steps'] = 172
        original.save(force=True)
        with patch('aptaswitch_core.oxdna_kinetics_worker.random_separated_configuration',
                   side_effect=AssertionError('resume must not reposition strands')):
            resumed = FFSReplica(req, folder, 0, adapter_factory=NativeScript([]))
        self.assertEqual(resumed.state['frame'], original.state['frame'])
        self.assertEqual(resumed.state['initialization'], original.state['initialization'])
        self.assertEqual(resumed.state['steps'], 172)
        self.assertEqual(resumed.rng.getstate(), original.rng.getstate())


if __name__ == '__main__':
    unittest.main()


class ProximityPlacementTests(unittest.TestCase):
    def setUp(self):
        self.root = Path(self.enterContext(tempfile.TemporaryDirectory()))
        self.request = request(self.root)
        switch = seed_strand('AAAGGGGAAAACCCC', '...((((....))))', 'rna')
        trigger = seed_strand('TTTAAAT', '.......', 'dna')
        self.n, self.total = len(switch), len(switch)+len(trigger)
        self.rows = np.vstack([switch, trigger])
        self.rows[:, :3] = (self.rows[:, :3]+99) % 100
        self.rows[:, 9:] = [1, 2, 3, 3, 2, 1]
        self.source = self.root/'initial.dat'
        _configuration(self.source, self.rows, 100)
        self.request['plan'].update(switch_length=self.n, nucleotide_count=self.total,
                                    toehold_pairs=[[i, self.n+2-i] for i in range(3)])
        self.request.update(switch_guide='...((((....))))', trigger_guide='.......')
        self.request['kinetics'].update(mode='proximity', initial_distance_min=2., initial_distance_max=2.5)

    def test_independent_reproducible_encounters_respect_all_periodic_pair_distances(self):
        frozen = self.source.read_bytes()
        frames = []
        for seed in range(12):
            frame, metadata = random_proximity_configuration(self.request, random.Random(seed))
            self.assertEqual((frame, metadata), random_proximity_configuration(self.request, random.Random(seed)))
            rows = coordinates(frame)
            delta = rows[:self.n, None, :3] - rows[None, self.n:, :3]
            delta -= 100*np.floor(delta/100+.5)
            distances = np.linalg.norm(delta, axis=2)
            self.assertGreaterEqual(float(distances.min()), 2.)
            self.assertLessEqual(float(distances[:3].min()), 2.5)
            self.assertAlmostEqual(metadata['minimum_interstrand_distance'], float(distances.min()))
            self.assertAlmostEqual(metadata['minimum_toehold_distance'], float(distances[:3].min()))
            i, j = metadata['anchor_pair']
            self.assertIn(i, (0, 1, 2))
            self.assertIn(j, range(self.n, self.total))
            self.assertGreaterEqual(metadata['periodic_self_check']['minimum_clearance'],
                                    metadata['periodic_self_check']['required_clearance'])
            self.assertNotEqual(metadata['strand_rotations'][0], metadata['strand_rotations'][1])
            self.assertTrue(metadata['artificial'])
            self.assertFalse(metadata['equilibrium_sample'])
            self.assertFalse(metadata['first_crossing_sample'])
            self.assertFalse(metadata['native_interstrand_bonds_checked'])
            frames.append(frame)
        self.assertEqual(len(set(frames)), len(frames))
        self.assertEqual(self.source.read_bytes(), frozen)

    def test_proximity_preserves_hairpin_geometry_chirality_and_momenta(self):
        frame, metadata = random_proximity_configuration(self.request, random.Random(37))
        rows = coordinates(frame)
        for start, end, rotation in ((0, self.n, metadata['strand_rotations'][0]),
                                     (self.n, self.total, metadata['strand_rotations'][1])):
            original, transformed = self.rows[start:end], rows[start:end]
            delta = original[:, None, :3]-original[None, :, :3]
            delta -= 100*np.floor(delta/100+.5)
            after = transformed[:, None, :3]-transformed[None, :, :3]
            np.testing.assert_allclose(np.linalg.norm(delta, axis=2), np.linalg.norm(after, axis=2), atol=1e-12)
            matrix = np.asarray(rotation)
            self.assertAlmostEqual(np.linalg.det(matrix), 1.)
            for offset in range(3, 15, 3):
                np.testing.assert_allclose(transformed[:, offset:offset+3], original[:, offset:offset+3] @ matrix.T)
        self.assertFalse(metadata['internal_conformations_changed'])

    def test_close_nontoehold_nucleotide_rejects_otherwise_valid_anchor(self):
        rows = np.zeros((3, 15))
        rows[:, 3] = rows[:, 8] = 1.
        rows[1, 1] = .761070781051
        _configuration(self.source, rows, 100)
        self.request['plan'].update(switch_length=2, nucleotide_count=3, toehold_pairs=[[0, 2]])
        self.request.update(switch_guide='..', trigger_guide='.')
        # First proposal approaches through the adjacent non-toehold base;
        # second approaches the same toehold from the opposite direction.
        draws = iter([.5, .25, .5, .5, .5, .5, .5, .75, .5, .5, .5, .5])
        rng = SimpleNamespace(random=lambda: next(draws), choice=lambda values: values[0], randrange=lambda stop: 0)
        with patch('aptaswitch_core.oxdna_initialization._rotation', return_value=np.eye(3)):
            frame, metadata = random_proximity_configuration(self.request, rng, max_attempts=2)
        self.assertEqual(metadata['attempts'], 2)
        result = coordinates(frame)
        self.assertLess(result[2, 1], result[0, 1])
        self.assertGreaterEqual(metadata['minimum_interstrand_distance'], 2.)

    def test_invalid_proximity_conditions_fail_before_a_random_draw(self):
        for key, value in (('mode', 'opening'), ('initial_distance_min', True),
                           ('initial_distance_min', 0.), ('initial_distance_min', 3.),
                           ('initial_distance_max', float('nan')), ('initial_distance_max', float('inf'))):
            previous = self.request['kinetics'][key]
            self.request['kinetics'][key] = value
            with self.subTest(key=key, value=value), \
                    patch('aptaswitch_core.oxdna_initialization._rotation', side_effect=AssertionError('unexpected draw')):
                with self.assertRaises(ValueError):
                    random_proximity_configuration(self.request, random.Random(3))
            self.request['kinetics'][key] = previous
        self.request['plan']['toehold_pairs'] = []
        with self.assertRaisesRegex(ValueError, 'toehold pairs'):
            random_proximity_configuration(self.request, random.Random(3))

    def test_proximity_rejects_unsafe_images_with_a_bounded_attempt_count(self):
        with patch('aptaswitch_core.oxdna_initialization.check_periodic_self_clearance',
                   side_effect=ValueError('periodic collision')) as guard:
            with self.assertRaisesRegex(ValueError, 'after 3 attempts'):
                random_proximity_configuration(self.request, random.Random(3), max_attempts=3)
            self.assertEqual(guard.call_count, 3)
        for attempts in (0, -1, True, 1.5):
            with self.subTest(attempts=attempts), self.assertRaises(ValueError):
                random_proximity_configuration(self.request, random.Random(3), max_attempts=attempts)
