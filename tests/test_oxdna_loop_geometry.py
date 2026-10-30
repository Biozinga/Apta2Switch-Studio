"""A guide-loop piercing is a local geometric diagnostic, not a knot proof."""
import unittest

import numpy as np

from aptaswitch_core.oxdna_geometry import configuration_geometry
from aptaswitch_core.oxdna_loop_geometry import configuration_loop_threading


def rows_from_backbones(points, chemistry='rna', a1=(1., 0., 0.), a3=(0., 0., 1.)):
    """Specify actual interaction sites independently of their COM offsets."""
    rows = np.zeros((len(points), 15))
    rows[:, 3:6], rows[:, 6:9] = a1, a3
    offsets = (-.4 * rows[:, 3:6] + .2 * rows[:, 6:9] if chemistry == 'rna'
               else -.34 * rows[:, 3:6] + .3408 * np.cross(rows[:, 6:9], rows[:, 3:6]))
    rows[:, :3] = np.asarray(points, dtype=float) - offsets
    return rows


def configuration(rows, box=(30., 30., 30.)):
    return ('t = 0\nb = ' + ' '.join(map(str, box)) + '\nE = 0 0 0\n'
            + '\n'.join(' '.join(map(str, row)) for row in rows) + '\n')


SQUARE = np.array([[-2., -2., 0.], [2., -2., 0.], [2., 2., 0.], [-2., 2., 0.]])


def square_with_segment(segment, chemistry='rna'):
    return rows_from_backbones(np.vstack((SQUARE, segment)), chemistry)


def check_square(segment, chemistry='rna'):
    return configuration_loop_threading(configuration(square_with_segment(segment, chemistry)),
                                        [4, 2], [chemistry, chemistry], ['(..)', '..'])


class ConfigurationLoopThreadingTests(unittest.TestCase):
    def assert_clear(self, result, faces=1):
        self.assertTrue(result['passed'], result)
        self.assertEqual(result['faces_checked'], faces)
        self.assertEqual(result['piercings'], [])
        self.assertEqual(result['unresolved'], [])
        self.assertFalse(result['knot_invariant'])
        self.assertIsInstance(result['method'], str)
        self.assertTrue(result['method'])

    def assert_pierced(self, result, pair, segment, faces=1):
        self.assertFalse(result['passed'], result)
        self.assertEqual(result['faces_checked'], faces)
        self.assertEqual(result['unresolved'], [])
        self.assertEqual({tuple(hit['loop_closing_pair']) for hit in result['piercings']},
                         {tuple(pair)})
        self.assertEqual({tuple(hit['segment']) for hit in result['piercings']},
                         {tuple(segment)})
        self.assertFalse(result['knot_invariant'])

    def test_threading_can_fail_even_when_all_local_geometry_checks_pass(self):
        angles = np.arange(12) * 2 * np.pi / 12
        circle = np.column_stack((np.zeros(12), 1.5 * np.cos(angles), 1.5 * np.sin(angles)))
        rows = np.vstack((rows_from_backbones(circle),
                          rows_from_backbones([[-.38, .15, .2], [.38, .15, .2]],
                                              a1=(0., 1., 0.))))
        frame = configuration(rows)
        local = configuration_geometry(frame, [12, 2], ['rna', 'rna'])
        self.assertTrue(local['passed'], local)
        self.assertEqual(local['repulsive_core_contacts'], 0)
        self.assertEqual(local['invalid_backbone_bonds'], [])
        result = configuration_loop_threading(frame, [12, 2], ['rna', 'rna'],
                                              ['(' + '.' * 10 + ')', '..'])
        self.assert_pierced(result, [1, 12], [13, 14])

    def test_threaded_bond_is_found_for_both_backbone_site_conventions(self):
        for chemistry in ('rna', 'dna'):
            with self.subTest(chemistry=chemistry):
                result = check_square([[.25, .15, -1.], [.25, .15, 1.]], chemistry)
                self.assert_pierced(result, [1, 4], [5, 6])

    def test_bond_through_square_center_is_not_lost_on_a_triangulation_diagonal(self):
        self.assert_pierced(check_square([[0., 0., -1.], [0., 0., 1.]]), [1, 4], [5, 6])

    def test_mixed_chemistries_use_backbone_sites_rather_than_centers(self):
        # The .06 difference in x offsets puts the COM bond on the wrong side
        # of a COM polygon while the real backbone bond crosses its aperture.
        for loop_acid, thread_acid, x in (('rna', 'dna', -1.97), ('dna', 'rna', 1.97)):
            with self.subTest(loop_acid=loop_acid, thread_acid=thread_acid):
                rows = np.vstack((rows_from_backbones(SQUARE, loop_acid),
                                  rows_from_backbones([[x, .15, -1.], [x, .15, 1.]], thread_acid)))
                result = configuration_loop_threading(configuration(rows), [4, 2],
                                                      [loop_acid, thread_acid], ['(..)', '..'])
                self.assert_pierced(result, [1, 4], [5, 6])

    def test_transverse_bond_outside_the_loop_is_clear(self):
        self.assert_clear(check_square([[3., .15, -1.], [3., .15, 1.]]))

    def test_endpoint_touch_coplanar_bond_and_polygon_boundary_graze_are_excluded(self):
        segments = {
            'endpoint on interior': [[.25, .15, 0.], [.25, .15, 1.]],
            'endpoint from below': [[.25, .15, -1.], [.25, .15, 0.]],
            'coplanar across aperture': [[-3., .15, 0.], [3., .15, 0.]],
            'edge graze': [[2., .15, -1.], [2., .15, 1.]],
            'vertex graze': [[2., 2., -1.], [2., 2., 1.]],
            'parallel above surface': [[-1., .15, .2], [1., .15, .2]],
        }
        for label, segment in segments.items():
            with self.subTest(case=label):
                self.assert_clear(check_square(segment))

    def test_v_shaped_tangent_touch_is_not_a_transverse_piercing(self):
        rows = rows_from_backbones(np.vstack((SQUARE, [[-.5, .15, 1.], [.25, .15, 0.],
                                                      [.8, .15, 1.]])))
        result = configuration_loop_threading(configuration(rows), [4, 3], ['rna', 'rna'],
                                              ['(..)', '...'])
        self.assert_clear(result)

    def test_foreign_chain_crossing_at_an_interior_vertex_reports_both_bonds(self):
        # A shared endpoint is a true crossing when its two neighbors occupy
        # opposite sides of the surface, unlike a free endpoint or V tangent.
        points = np.vstack((SQUARE, [[.25, .15, -1.], [.25, .15, 0.], [.25, .15, 1.]]))
        result = configuration_loop_threading(configuration(rows_from_backbones(points)),
                                              [4, 3], ['rna', 'rna'], ['(..)', '...'])
        self.assertFalse(result['passed'], result)
        self.assertEqual(result['faces_checked'], 1)
        self.assertEqual(result['unresolved'], [])
        self.assertEqual({tuple(hit['loop_closing_pair']) for hit in result['piercings']},
                         {(1, 4)})
        self.assertEqual({tuple(hit['segment']) for hit in result['piercings']},
                         {(5, 6), (6, 7)})
        self.assertFalse(result['knot_invariant'])

    def test_rigid_rotations_and_translations_preserve_the_piercing(self):
        rotation, _ = np.linalg.qr(np.array([[2., 1., 3.], [-1., 2., 2.], [3., -2., 1.]]))
        self.assertAlmostEqual(np.linalg.det(rotation), 1.)
        rows = square_with_segment([[.25, .15, -1.], [.25, .15, 1.]])
        for transform, translation in ((np.eye(3), [9., -7., 5.]),
                                       (rotation, [7., 6., 8.]),
                                       (rotation.T, [-9., 4., -8.])):
            with self.subTest(translation=translation):
                changed = rows.copy()
                changed[:, :3] = changed[:, :3] @ transform + translation
                changed[:, 3:6] = changed[:, 3:6] @ transform
                changed[:, 6:9] = changed[:, 6:9] @ transform
                result = configuration_loop_threading(configuration(changed), [4, 2],
                                                      ['rna', 'rna'], ['(..)', '..'])
                self.assert_pierced(result, [1, 4], [5, 6])

    def test_periodic_images_of_a_separate_strand_are_checked(self):
        for translation in ([30., 0., 0.], [-30., 30., 0.], [0., 0., 30.]):
            with self.subTest(translation=translation):
                segment = np.array([[.25, .15, -1.], [.25, .15, 1.]]) + translation
                self.assert_pierced(check_square(segment), [1, 4], [5, 6])

    def test_no_covalent_bond_is_invented_between_two_strands(self):
        rows = square_with_segment([[.25, .15, -1.], [.25, .15, 1.]])
        result = configuration_loop_threading(configuration(rows), [4, 1, 1],
                                              ['rna', 'rna', 'rna'], ['(..)', '.', '.'])
        self.assert_clear(result)

    def test_unpaired_strands_have_no_guide_loop_faces(self):
        rows = square_with_segment([[.25, .15, -1.], [.25, .15, 1.]])
        result = configuration_loop_threading(configuration(rows), [4, 2], ['rna', 'rna'],
                                              ['....', '..'])
        self.assert_clear(result, faces=0)

    def test_nested_bulge_and_hairpin_faces_are_checked_without_a_stacked_face(self):
        # Pair 1-11 closes a bulge whose boundary jumps over its child 3-9.
        # Pair 3-9 is stacked onto 4-8, which closes the right-hand hairpin.
        points = np.array([[-3., -2., 0.], [-3., 2., 0.], [-1., 2., 0.],
                           [0., 2., 0.], [2., 2., 0.], [3., 0., 0.],
                           [2., -2., 0.], [0., -2., 0.], [-1., -2., 0.],
                           [-1.7, -3., 0.], [-2.3, -3., 0.]])
        for xy, pair in (([-2.25, .15], [1, 11]), ([1.4, .2], [4, 8])):
            with self.subTest(pair=pair):
                segment = [[*xy, -1.], [*xy, 1.]]
                rows = rows_from_backbones(np.vstack((points, segment)))
                result = configuration_loop_threading(configuration(rows), [11, 2],
                                                      ['rna', 'rna'], ['(.((...)).)', '..'])
                self.assert_pierced(result, pair, [12, 13], faces=2)

    def test_same_strand_can_pierce_its_own_nonadjacent_loop(self):
        # The first two nucleotides lie before the guide-paired loop.
        points = np.vstack(([[.25, .15, -1.], [.25, .15, 1.]], SQUARE))
        result = configuration_loop_threading(configuration(rows_from_backbones(points)),
                                              [6], ['rna'], ['..(..)'])
        self.assert_pierced(result, [3, 6], [1, 2])

    def test_degenerate_loop_is_unresolved_and_never_silently_passes(self):
        for points in ([[-2., 0., 0.], [-1., 0., 0.], [1., 0., 0.], [2., 0., 0.]],
                        [[0., 0., 0.]] * 4):
            with self.subTest(points=points):
                result = configuration_loop_threading(configuration(rows_from_backbones(points)),
                                                      [4], ['rna'], ['(..)'])
                self.assertFalse(result['passed'])
                self.assertEqual(result['piercings'], [])
                self.assertTrue(result['unresolved'])
                self.assertTrue(all(entry.get('reason') for entry in result['unresolved']))
                self.assertFalse(result['knot_invariant'])

    def test_self_crossing_planar_loop_is_unresolved(self):
        bow_tie = SQUARE[[0, 2, 1, 3]]
        result = configuration_loop_threading(configuration(rows_from_backbones(bow_tie)),
                                              [4], ['rna'], ['(..)'])
        self.assertFalse(result['passed'])
        self.assertEqual(result['piercings'], [])
        self.assertTrue(result['unresolved'])

    def test_guides_are_strictly_validated(self):
        rows = rows_from_backbones(SQUARE)
        for guide in ('(...', '...)', ')(..', '([)]', '[..]', '(..x', '(.. )', '...', ''):
            with self.subTest(guide=guide), self.assertRaises(ValueError):
                configuration_loop_threading(configuration(rows), [4], ['rna'], [guide])

    def test_strand_metadata_and_guide_count_are_strictly_validated(self):
        frame = configuration(rows_from_backbones(SQUARE))
        cases = (
            ([], [], []),
            ([4], ['rna'], []),
            ([4], ['rna'], ['(..)', '']),
            ([4], [], ['(..)']),
            ([4], ['protein'], ['(..)']),
            ([0, 4], ['rna', 'rna'], ['', '(..)']),
            ([True, 3], ['rna', 'rna'], ['.', '...']),
            ([4.0], ['rna'], ['(..)']),
            ([3], ['rna'], ['...']),
        )
        for lengths, acids, guides in cases:
            with self.subTest(lengths=lengths, acids=acids, guides=guides), self.assertRaises(ValueError):
                configuration_loop_threading(frame, lengths, acids, guides)


if __name__ == '__main__':
    unittest.main()
