"""Preserve explicit box geometry and strict native checkpoint validation."""
from pathlib import Path
import tempfile
import unittest

import numpy as np

from aptaswitch_core.oxdna_geometry import (
    _number, _place_strands, _read_configuration, configuration_geometry,
)
from aptaswitch_core.oxdna_guided import seed_strand


def particles(count):
    rows = np.zeros((count, 15))
    rows[:, 0] = np.arange(count) * .76
    rows[:, 3] = rows[:, 8] = 1
    rows[:, 9:] = .25
    return rows


def configuration(rows, box=30):
    return (f't = 17\nb = {box} {box} {box}\nE = -1 -2 1\n'
            + '\n'.join(' '.join(map(str, row)) for row in rows) + '\n')


class OxDNAGeometryTests(unittest.TestCase):
    def setUp(self):
        self.root = Path(self.enterContext(tempfile.TemporaryDirectory()))
        self.path = self.root / 'configuration.dat'

    def test_checkpoint_preserves_positions_axes_and_momenta(self):
        rows = particles(3)
        self.path.write_text(configuration(rows))
        parsed, box = _read_configuration(self.path, 3, box_size=30)
        np.testing.assert_array_equal(parsed, rows)
        np.testing.assert_array_equal(box, [30, 30, 30])

    def test_nonfinite_nonorthogonal_or_truncated_configuration_is_rejected(self):
        for defect in ('position', 'momentum', 'nonorthogonal', 'nonunit', 'truncated'):
            rows = particles(3)
            if defect == 'position':
                rows[1, 1] = float('nan')
            elif defect == 'momentum':
                rows[1, 14] = float('inf')
            elif defect == 'nonorthogonal':
                rows[1, 6:9] = rows[1, 3:6]
            elif defect == 'nonunit':
                rows[1, 3:6] *= 2
            else:
                rows = rows[:-1]
            self.path.write_text(configuration(rows))
            with self.subTest(defect=defect), self.assertRaises(ValueError):
                _read_configuration(self.path, 3, box_size=30)

    def test_multiple_frames_and_changed_or_malformed_headers_are_rejected(self):
        text = configuration(particles(3))
        for malformed in (text + text, text.replace('t = 17', 't = nan'),
                          text.replace('E = -1 -2 1', 'E = -1 -2'),
                          text.replace('b = 30 30 30', 'b = 30 30 31'),
                          text.replace('b = 30 30 30', 'b = 0 0 0')):
            self.path.write_text(malformed)
            with self.subTest(header=malformed[:60]), self.assertRaises(ValueError):
                _read_configuration(self.path, 3, box_size=30)

    def test_separation_preserves_both_shapes_axes_and_periodic_gap(self):
        first, second = particles(8), particles(5)
        originals = first.copy(), second.copy()
        rows, box = _place_strands(first, second, 6, 30)
        np.testing.assert_array_equal(first, originals[0])
        np.testing.assert_array_equal(second, originals[1])
        for placed, original in ((rows[:8], first), (rows[8:], second)):
            np.testing.assert_allclose(np.diff(placed[:, :3], axis=0), np.diff(original[:, :3], axis=0))
            np.testing.assert_array_equal(placed[:, 3:9], original[:, 3:9])
        self.assertAlmostEqual(rows[8:, 0].min() - rows[:8, 0].max(), 6)
        self.assertGreaterEqual(box - (rows[8:, 0].max() - rows[:8, 0].min()), 6)
        np.testing.assert_array_equal(rows[:, 9:], 0)
        self.assertEqual(box, 30)

    def test_small_box_is_rejected_without_resizing(self):
        with self.assertRaisesRegex(ValueError, 'too small'):
            _place_strands(particles(8), particles(5), 6, 10)

    def test_automatic_placement_rotates_long_strands_without_deforming_them(self):
        first, second = particles(45), particles(45)
        placed, side = _place_strands(first, second, 10, 42,
                                     periodic_margin=8, auto_size=True, reorient=True)
        self.assertEqual(side, 42)  # Not both contour lengths end-to-end along X.
        for before, after in ((first, placed[:45]), (second, placed[45:])):
            original_d = np.linalg.norm(before[:, None, :3]-before[None, :, :3], axis=2)
            rotated_d = np.linalg.norm(after[:, None, :3]-after[None, :, :3], axis=2)
            np.testing.assert_allclose(original_d, rotated_d, atol=1e-12)
            # Signed scalar triple product is unchanged by a proper rotation.
            old = np.dot(np.cross(before[0,3:6], before[0,6:9]), before[1,:3]-before[0,:3])
            new = np.dot(np.cross(after[0,3:6], after[0,6:9]), after[1,:3]-after[0,:3])
            self.assertAlmostEqual(old, new)
            np.testing.assert_allclose(np.linalg.norm(after[:,3:6], axis=1), 1)
        self.assertGreaterEqual(side-np.ptp(placed[:,0]), 10)
        self.assertTrue(all(side-np.ptp(placed[:,axis]) >= 8 for axis in (1,2)))

    def test_only_preparation_expands_box_for_actual_placement(self):
        placed, side = _place_strands(particles(20), particles(20), 10, 22,
                                     periodic_margin=8, auto_size=True, reorient=True)
        self.assertGreater(side, 22)
        self.assertTrue(all(side-np.ptp(placed[:,axis]) >= 8 for axis in range(3)))

    def test_numeric_parameters_reject_nonfinite_or_nonpositive_values(self):
        for value in (True, None, float('nan'), float('inf'), 'invalid'):
            with self.subTest(value=value), self.assertRaises(ValueError):
                _number(value, 'parameter')
        for value in (0, -1):
            with self.subTest(value=value), self.assertRaises(ValueError):
                _number(value, 'parameter', positive=True)
        self.assertEqual(_number(-5, 'temperature'), -5)


class ConfigurationGeometryTests(unittest.TestCase):
    def test_nonbonded_collisions_use_each_interaction_site_not_centers(self):
        # Each pair has separated centers and exactly one coincident site pair.
        fixtures = {
            'rna': {
                'backbone-backbone': ((-.8, 0, 0), True),
                'base-base': ((.8, 0, 0), True),
                'backbone-base': ((-.8, 0, .2), False),
                'base-backbone': ((.8, 0, -.2), False),
            },
            'dna': {
                'backbone-backbone': ((-.68, .6816, 0), True),
                'base-base': ((.8, 0, 0), True),
                'backbone-base': ((-.74, .3408, 0), False),
                'base-backbone': ((.74, -.3408, 0), False),
            },
        }
        for acid, contacts in fixtures.items():
            for site, (position, reverse_base) in contacts.items():
                with self.subTest(acid=acid, site=site):
                    rows = particles(2)
                    rows[1, :3] = position
                    if reverse_base:
                        rows[1, 3:6] *= -1
                    self.assertGreater(np.linalg.norm(rows[1, :3]), .675)
                    result = configuration_geometry(configuration(rows), (1, 1), (acid, acid))
                    self.assertFalse(result['passed'])
                    self.assertEqual(result['invalid_backbone_bonds'], [])
                    self.assertEqual(result['repulsive_core_contacts'], 1)
                    contact, = result['contact_examples']
                    self.assertEqual(contact['nucleotides'], [1, 2])
                    self.assertEqual(contact['sites'], site)
                    self.assertAlmostEqual(contact['distance'], 0)

    def test_adjacent_backbones_are_excluded_but_strand_boundaries_are_not(self):
        for acid in ('rna', 'dna'):
            with self.subTest(acid=acid):
                rows = particles(2)
                rows[1, :3] = [0, 0, .6]
                bonded = configuration_geometry(configuration(rows), (2,), (acid,))
                self.assertTrue(bonded['passed'])
                self.assertIsNone(bonded['minimum_site_distances']['backbone-backbone'])
                nonbonded = configuration_geometry(configuration(rows), (1, 1), (acid, acid))
                self.assertFalse(nonbonded['passed'])
                self.assertEqual(nonbonded['repulsive_core_contacts'], 1)
                self.assertEqual(nonbonded['contact_examples'][0]['sites'], 'backbone-backbone')

    def test_adjacent_base_contacts_are_checked_even_with_valid_backbone_bonds(self):
        for acid in ('rna', 'dna'):
            for site in ('base-base', 'base-backbone', 'backbone-base'):
                with self.subTest(acid=acid, site=site):
                    rows = particles(2)
                    if site == 'base-base':
                        rows[1, 3:6] = [.5, np.sqrt(3)/2, 0]
                        rows[1, :3] = .4*(rows[0, 3:6]-rows[1, 3:6])
                    else:
                        rows[1, :3] = ((.8, 0, -.2) if acid == 'rna'
                                       else (.74, -.3408, 0))
                        if site == 'backbone-base':
                            rows[1, :3] *= -1
                    result = configuration_geometry(configuration(rows), (2,), (acid,))
                    self.assertFalse(result['passed'])
                    self.assertEqual(result['invalid_backbone_bonds'], [])
                    self.assertEqual(result['repulsive_core_contacts'], 1)
                    self.assertEqual(result['contact_examples'][0]['sites'], site)

    def test_fene_limits_reject_compressed_and_stretched_bonds_for_each_chemistry(self):
        for acid in ('rna', 'dna'):
            for distance in (.49, 1.1):
                with self.subTest(acid=acid, distance=distance):
                    rows = particles(2)
                    rows[1, :3] = [0, 0, distance]
                    result = configuration_geometry(configuration(rows), (2,), (acid,))
                    self.assertFalse(result['passed'])
                    self.assertEqual(result['invalid_backbone_bonds'], [[1, 2]])
                    self.assertEqual(result['repulsive_core_contacts'], 0)
                    self.assertGreater(result['maximum_backbone_extension'], .25)

    def test_fene_limits_use_the_requested_rna_or_dna_rest_length(self):
        for distance, valid_acid in ((.509, 'dna'), (1.009, 'rna')):
            for acid in ('rna', 'dna'):
                with self.subTest(acid=acid, distance=distance):
                    rows = particles(2)
                    rows[1, :3] = [0, 0, distance]
                    result = configuration_geometry(configuration(rows), (2,), (acid,))
                    self.assertEqual(result['passed'], acid == valid_acid)
                    self.assertEqual(result['invalid_backbone_bonds'],
                                     [] if acid == valid_acid else [[1, 2]])
                    self.assertEqual(result['repulsive_core_contacts'], 0)

    def test_nonbonded_site_collisions_are_detected_across_periodic_boundary(self):
        rows = particles(2)
        rows[1, :3] = [30.8, 0, -.2]
        result = configuration_geometry(configuration(rows), (1, 1), ('rna', 'rna'))
        self.assertFalse(result['passed'])
        self.assertEqual(result['repulsive_core_contacts'], 1)
        contact, = result['contact_examples']
        self.assertEqual(contact['sites'], 'base-backbone')
        self.assertAlmostEqual(contact['distance'], 0)

    def test_periodic_images_cannot_repair_an_unwrapped_covalent_bond(self):
        rows = particles(2)
        rows[1, :3] = [0, 0, 29.24]
        result = configuration_geometry(configuration(rows), (2,), ('rna',))
        self.assertFalse(result['passed'])
        self.assertEqual(result['invalid_backbone_bonds'], [[1, 2]])
        self.assertEqual(result['repulsive_core_contacts'], 0)

    def test_duplex_stems_from_seed_strand_pass_without_claiming_equilibration(self):
        for acid in ('rna', 'dna'):
            with self.subTest(acid=acid):
                hairpin = seed_strand('GGGGAAAACCCC', '((((....))))', acid)
                duplex = np.vstack((hairpin[:4], hairpin[-4:]))
                result = configuration_geometry(configuration(duplex), (4, 4), (acid, acid))
                self.assertTrue(result['passed'])
                self.assertEqual(result['invalid_backbone_bonds'], [])
                self.assertEqual(result['repulsive_core_contacts'], 0)
                self.assertFalse(result['equilibrium_validated'])
                self.assertFalse(result['topology_validated'])

    def test_geometry_is_invariant_under_a_proper_rigid_rotation(self):
        rotation, _ = np.linalg.qr(np.array([[2., 1, 3], [-1, 2, 2], [3, -2, 1]]))
        self.assertAlmostEqual(np.linalg.det(rotation), 1)
        for acid in ('rna', 'dna'):
            hairpin = seed_strand('GGGGAAAACCCC', '((((....))))', acid)
            duplex = np.vstack((hairpin[:4], hairpin[-4:]))
            for defective in (False, True):
                with self.subTest(acid=acid, defective=defective):
                    rows = duplex.copy()
                    if defective:
                        rows[4] = rows[0]
                    rotated = rows.copy()
                    for start in range(0, 15, 3):
                        rotated[:, start:start+3] = rows[:, start:start+3] @ rotation
                    rotated[:, :3] += [7, 8, 9]
                    before = configuration_geometry(configuration(rows), (4, 4), (acid, acid))
                    after = configuration_geometry(configuration(rotated), (4, 4), (acid, acid))
                    self.assertEqual(after['passed'], before['passed'])
                    self.assertEqual(after['invalid_backbone_bonds'], before['invalid_backbone_bonds'])
                    self.assertEqual(after['repulsive_core_contacts'], before['repulsive_core_contacts'])
                    self.assertAlmostEqual(after['maximum_backbone_extension'],
                                           before['maximum_backbone_extension'])
                    for site, distance in before['minimum_site_distances'].items():
                        self.assertAlmostEqual(after['minimum_site_distances'][site], distance)

    def test_soft_contacts_outside_repulsive_cores_are_allowed(self):
        fixtures = (
            ('backbone-backbone', (0, 0, .7), False, .7),
            ('base-base', (1.14, 0, 0), True, .34),
            ('base-backbone', (.8, 0, .32), False, .52),
            ('backbone-base', (-.8, 0, -.32), False, .52),
        )
        for site, position, reverse_base, distance in fixtures:
            with self.subTest(site=site):
                rows = particles(2)
                rows[1, :3] = position
                if reverse_base:
                    rows[1, 3:6] *= -1
                result = configuration_geometry(configuration(rows), (1, 1), ('rna', 'rna'))
                self.assertTrue(result['passed'])
                self.assertEqual(result['repulsive_core_contacts'], 0)
                self.assertAlmostEqual(result['minimum_site_distances'][site], distance)


if __name__ == '__main__':
    unittest.main()
