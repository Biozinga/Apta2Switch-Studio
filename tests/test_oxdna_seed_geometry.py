"""Exterior seeds must not superpose distant parts of the covalent chain."""
import unittest

import numpy as np

from aptaswitch_core.oxdna_geometry import configuration_geometry
from aptaswitch_core.oxdna_guided import _backbone, guide_pairs, seed_strand


# Recorded native-guide architecture: long 5-prime tail, internal bulges,
# and a second exterior stem after a 30-base connector. No runtime is used.
SWITCH_SEQUENCE = (
    'GGGUUAUACCAGCACCAUGUAUCUGUGCGUCGCGAGAGAGGAGAGUGACGAUGAGAUAUAUGAAUAAU'
    'AAAAAGAAGAACAAAGUGAGCAAGGGCGAGGAGGAUAACAUGGCC'
)
SWITCH_GUIDE = (
    '...............(((((((((...((((((...........))))))...)))))))))'
    '..............................(((.((........).).)))'
)


def sequence_for(guide):
    """Canonical GC pairs, with adenine at unpaired positions."""
    return ''.join('G' if mark == '(' else 'C' if mark == ')' else 'A' for mark in guide)


def geometry(rows, material):
    # Keep all images far from this isolated strand; this checks construction,
    # while preparation chooses and validates the actual simulation box.
    box = float(np.ptp(rows[:, :3], axis=0).max()+100.)
    frame = (f't = 0\nb = {box} {box} {box}\nE = 0 0 0\n'
             + '\n'.join(' '.join(map(str, row)) for row in rows))
    return configuration_geometry(frame, [len(rows)], [material])


class ExteriorSeedGeometryTests(unittest.TestCase):
    def assert_guide_and_backbone(self, rows, sequence, guide, material):
        self.assertTrue(np.isfinite(rows).all())
        a1, a3 = rows[:, 3:6], rows[:, 6:9]
        frames = np.stack([a1, np.cross(a3, a1), a3], axis=-1)
        np.testing.assert_allclose(np.linalg.det(frames), 1., atol=1e-12)
        np.testing.assert_allclose(np.linalg.norm(a1, axis=1), 1., atol=1e-12)
        np.testing.assert_allclose(np.linalg.norm(a3, axis=1), 1., atol=1e-12)
        np.testing.assert_allclose(np.sum(a1*a3, axis=1), 0., atol=1e-12)
        for first, second in guide_pairs(sequence, guide, material):
            self.assertAlmostEqual(np.linalg.norm(rows[first, :3]-rows[second, :3]), 1.2)
            np.testing.assert_allclose(rows[second, :3]-rows[first, :3], 1.2*a1[first], atol=1e-12)
            np.testing.assert_allclose(rows[first, 3:9], -rows[second, 3:9], atol=1e-12)
        distances = np.linalg.norm(np.diff(_backbone(rows, material), axis=0), axis=1)
        # Native FENE finite domain, using each chemistry's equilibrium length.
        rest = .761070781051 if material == 'rna' else .7564
        self.assertTrue(np.all(np.abs(distances-rest) < .25))

    def test_multiple_exterior_stems_do_not_return_to_previous_coordinates(self):
        hairpin = '(((....)))'
        for material in ('rna', 'dna'):
            for connector in ('', '.', '.'*10, '.'*30):
                guide = connector.join([hairpin]*4)
                sequence = sequence_for(guide)
                for seed in (None, 0, 7):
                    with self.subTest(material=material, connector=len(connector), seed=seed):
                        rows = seed_strand(sequence, guide, material, seed=seed)
                        self.assert_guide_and_backbone(rows, sequence, guide, material)
                        screen = geometry(rows, material)
                        self.assertTrue(screen['passed'], screen)
                        # The formerly superposed first/third stems remain
                        # separated even with zero unpaired connector bases.
                        third = 2*(len(hairpin)+len(connector))
                        self.assertGreater(np.linalg.norm(rows[0, :3]-rows[third, :3]), 2.)

    def test_long_tail_and_exterior_connector_do_not_overlap_each_other(self):
        guide = '.'*15+'(((....)))'+'.'*30+'(((....)))'+'.'*15
        sequence = sequence_for(guide)
        for material in ('rna', 'dna'):
            for seed in (None, 0, 7):
                with self.subTest(material=material, seed=seed):
                    rows = seed_strand(sequence, guide, material, seed=seed)
                    self.assert_guide_and_backbone(rows, sequence, guide, material)
                    screen = geometry(rows, material)
                    self.assertTrue(screen['passed'], screen)
                    backbone = _backbone(rows, material)
                    # Thirty unpaired bases give 31 nominal covalent links;
                    # spacing does not stretch a single bond across the gap.
                    connector = backbone[24:56]
                    np.testing.assert_allclose(np.linalg.norm(np.diff(connector, axis=0), axis=1), .76)

    def test_recorded_native_guide_has_no_core_contacts_including_internal_bulge(self):
        for seed in (None, *range(20)):
            with self.subTest(seed=seed):
                rows = seed_strand(SWITCH_SEQUENCE, SWITCH_GUIDE, 'rna', seed=seed)
                self.assert_guide_and_backbone(rows, SWITCH_SEQUENCE, SWITCH_GUIDE, 'rna')
                # This pair used to coincide to 0.04 units in the first seed.
                self.assertGreater(np.linalg.norm(rows[11, :3]-rows[112, :3]), 10.)
                screen = geometry(rows, 'rna')
                self.assertTrue(screen['passed'], screen)

    def test_internal_bulges_extrude_away_from_their_paired_partners(self):
        guides = (
            '(((...(((....)))...)))',
            '((.((....)).....................))',
            '((.....................((....)).))',
            SWITCH_GUIDE,
        )
        for material in ('rna', 'dna'):
            for guide in guides:
                sequence = sequence_for(guide)
                pairs = guide_pairs(sequence, guide, material)
                partners = {i: j for a, b in pairs for i, j in ((a, b), (b, a))}
                anchors = sorted(partners)
                reference = seed_strand(sequence, guide, material)
                for seed in (None, 0, 7, 19):
                    with self.subTest(material=material, guide=guide, seed=seed):
                        rows = seed_strand(sequence, guide, material, seed=seed)
                        self.assert_guide_and_backbone(rows, sequence, guide, material)
                        backbone = _backbone(rows, material)
                        for left, right in zip(anchors, anchors[1:]):
                            # Same-side partners identify the two arms of an
                            # internal loop; omit exterior chains and hairpins.
                            if right-left <= 1 or not (
                                    partners[left] > right and partners[right] > right
                                    or partners[left] < left and partners[right] < left):
                                continue
                            chord = backbone[right]-backbone[left]
                            axis = chord/np.linalg.norm(chord)
                            midpoint = (backbone[left]+backbone[right])/2
                            outward = midpoint-backbone[[partners[left], partners[right]]].mean(axis=0)
                            outward -= axis*np.dot(outward, axis)
                            outward /= np.linalg.norm(outward)
                            offsets = backbone[left+1:right]-midpoint
                            offsets -= (offsets@axis)[:, None]*axis
                            cosines = offsets@outward/np.linalg.norm(offsets, axis=1)
                            # Variants may rotate the arc by at most .65 rad,
                            # while every interior point stays on the outer
                            # side of the chord rather than following the stem.
                            self.assertTrue(np.all(cosines >= np.cos(.65)-1e-12), cosines)
                            np.testing.assert_allclose(
                                np.linalg.norm(np.diff(backbone[left:right+1], axis=0), axis=1), .76)
                        # In a single-root architecture, random loop planes
                        # never move or rotate the paired stem anchors.
                        if guide != SWITCH_GUIDE:
                            np.testing.assert_array_equal(rows[anchors], reference[anchors])


if __name__ == '__main__':
    unittest.main()
