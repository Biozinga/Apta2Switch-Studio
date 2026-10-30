"""Generic construction checks; the native NA model still validates each seed."""
import unittest
import random

import numpy as np

from aptaswitch_core.oxdna_guided import _backbone, guide_pairs, seed_strand


def sequence_for(guide, *, material='rna', wobble=False):
    sequence = ['A']*len(guide)
    stack = []
    for index, mark in enumerate(guide):
        if mark == '(':
            stack.append(index)
        elif mark == ')':
            first = stack.pop()
            sequence[first] = 'G'
            sequence[index] = 'U' if material == 'rna' and wobble else 'C'
    return ''.join(sequence)


class SeedVariantsTests(unittest.TestCase):
    GUIDES = (
        '....((((....))))....',
        '((.((....)).....................))',
        '((.....................((....)).))',
        '((((....))))((((....))))',
        '(((...(((....)))...)))..(((....)))',
        '((..((....))..((....))..))',
        '((..((....))..((....))..((....))..))',
        '((((...)))(((...))))',
        '((..((....))...................................((....))..))',
    )

    def test_diverse_architectures_retain_indices_pairs_and_proper_frames(self):
        for material in ('rna', 'dna'):
            for guide in self.GUIDES:
                sequence = sequence_for(guide, material=material)
                for seed in (None, 0, 7):
                    with self.subTest(material=material, guide=guide, seed=seed):
                        rows = seed_strand(sequence, guide, material, seed=seed)
                        self.assertEqual(rows.shape, (len(sequence), 15))
                        self.assertTrue(np.isfinite(rows).all())
                        a1, a3 = rows[:, 3:6], rows[:, 6:9]
                        frames = np.stack([a1, np.cross(a3, a1), a3], axis=-1)
                        np.testing.assert_allclose(np.linalg.det(frames), 1, atol=1e-12)
                        np.testing.assert_allclose(np.linalg.norm(a1, axis=1), 1, atol=1e-12)
                        np.testing.assert_allclose(np.linalg.norm(a3, axis=1), 1, atol=1e-12)
                        np.testing.assert_allclose(np.sum(a1*a3, axis=1), 0, atol=1e-12)
                        for first, second in guide_pairs(sequence, guide, material):
                            self.assertAlmostEqual(np.linalg.norm(rows[first, :3]-rows[second, :3]), 1.2)
                            np.testing.assert_allclose(rows[first, 3:9], -rows[second, 3:9], atol=1e-12)

    def test_short_connectors_are_not_stretched_by_bulges_or_other_hairpins(self):
        for guide in self.GUIDES:
            for seed in (None, 1, 11):
                with self.subTest(guide=guide, seed=seed):
                    rows = seed_strand(sequence_for(guide), guide, 'rna', seed=seed)
                    distances = np.linalg.norm(np.diff(_backbone(rows, 'rna'), axis=0), axis=1)
                    # RNA FENE bounds in the native model are r0 +/- .25.
                    # This only checks backbone geometry, not native energies.
                    self.assertGreater(float(distances.min()), .511070781051)
                    self.assertLess(float(distances.max()), 1.011070781051)

    def test_variants_are_reproducible_and_change_internal_geometry(self):
        for guide in ('....((((....))))....', '((..((....))..((....))..))'):
            sequence = sequence_for(guide)
            first = seed_strand(sequence, guide, 'rna', seed=19)
            repeat = seed_strand(sequence, guide, 'rna', seed=19)
            second = seed_strand(sequence, guide, 'rna', seed=20)
            np.testing.assert_array_equal(first, repeat)
            distances = lambda rows: np.linalg.norm(rows[:, None, :3]-rows[None, :, :3], axis=-1)
            self.assertGreater(float(np.max(np.abs(distances(first)-distances(second)))), .01)

    def test_unpaired_variants_keep_covalent_contour_and_chemistry(self):
        for material, sequence in (('dna', 'ACGT'*7+'AC'), ('rna', 'ACGU'*6)):
            guide = '.'*len(sequence)
            plain = seed_strand(sequence, guide, material)
            rows = seed_strand(sequence, guide, material, seed=6)
            np.testing.assert_allclose(np.linalg.norm(np.diff(_backbone(rows, material), axis=0), axis=1), .76)
            self.assertFalse(np.allclose(plain, rows))
            self.assertEqual(guide_pairs(sequence, guide, material), [])

    def test_wobble_pairs_remain_the_requested_pairs(self):
        guide = '...(((....)))...'
        sequence = sequence_for(guide, wobble=True)
        pairs = guide_pairs(sequence, guide, 'rna')
        rows = seed_strand(sequence, guide, 'rna', seed=2)
        self.assertEqual(len(pairs), 3)
        for first, second in pairs:
            self.assertEqual(sequence[first]+sequence[second], 'GU')
            self.assertAlmostEqual(np.linalg.norm(rows[first, :3]-rows[second, :3]), 1.2)

    def test_invalid_seeds_fail_before_construction(self):
        for seed in (-1, True, 1.2, '7', 2**31):
            with self.subTest(seed=seed), self.assertRaises(ValueError):
                seed_strand('GAAAC', '(...)', 'rna', seed=seed)

    def test_nested_junctions_and_long_loops_keep_discrete_backbone_lengths(self):
        rng = random.Random(127)

        def nested(depth):
            stem = rng.randint(1, 3)
            if not depth:
                inside = '.'*rng.randint(3, 12)
            elif rng.random() < .5:
                inside = '.'*rng.randint(0, 8)+nested(depth-1)+'.'*rng.randint(0, 8)
            else:
                inside = ('.'*rng.randint(0, 3)+nested(depth-1)+'.'*rng.randint(0, 3)
                          +nested(depth-1)+'.'*rng.randint(0, 3))
            return '('*stem+inside+')'*stem

        for seed in range(50):
            guide = '.'*rng.randint(0, 8)+nested(3)+'.'*rng.randint(0, 8)
            rows = seed_strand(sequence_for(guide), guide, 'rna', seed=seed)
            distances = np.linalg.norm(np.diff(_backbone(rows, 'rna'), axis=0), axis=1)
            with self.subTest(seed=seed, guide=guide):
                self.assertTrue(np.isfinite(rows).all())
                self.assertGreater(float(distances.min()), .75999)
                self.assertLess(float(distances.max()), .76108)


if __name__ == '__main__':
    unittest.main()
