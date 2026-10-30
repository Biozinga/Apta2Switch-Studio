"""Independent combinatorial checks for the thermodynamic validation helpers."""
from __future__ import annotations

from itertools import combinations
import math
import unittest

from tests.native_reference.reference_checks import (
    boltzmann_summary, connected, enum_structures, pair_matrix,
)


class EnumerationTests(unittest.TestCase):
    def test_unpaired_sequence_and_minimal_hairpin(self):
        self.assertEqual(list(enum_structures("AAAA")), ["...."])
        self.assertEqual(set(enum_structures("GAAAC")), {".....", "(...)"})
        self.assertEqual(list(enum_structures("GAAC")), ["...."])
        self.assertIn("(..)", set(enum_structures("GAAC", min_loop_size=2)))

    def test_interstrand_pair_has_no_minimum_separation(self):
        self.assertEqual(list(enum_structures("GC")), [".."])
        self.assertEqual(set(enum_structures("GC", (1,))), {"..", "()"})
        self.assertEqual(set(enum_structures("GCGC", (2,))),
                         {"....", "(..)", ".().", "(())"})
        # A break elsewhere does not relax the hairpin rule within strand 1.
        self.assertNotIn("()..", set(enum_structures("GCAA", (2,))))

    def test_wobble_and_dna_alias_are_explicit(self):
        self.assertIn("()", set(enum_structures("GU", (1,))))
        self.assertEqual(set(enum_structures("gt", (1,))), set(enum_structures("GU", (1,))))
        self.assertEqual(set(enum_structures("GATTC", (2,))), set(enum_structures("GAUUC", (2,))))
        self.assertEqual(list(enum_structures("AC", (1,))), [".."])

    def test_enumeration_matches_independent_pair_subset_oracle(self):
        # The oracle enumerates subsets of all candidate edges, not recursive
        # interval decompositions. It independently checks occupancy/crossings.
        for sequence, cuts, minimum in (("GCGCUA", (), 0), ("GCGCUA", (2, 4), 3),
                                         ("GAUACU", (), 3), ("GCGC", (1, 2, 3), 3)):
            strands = [sum(point <= index for point in cuts) for index in range(len(sequence))]
            edges = [(i, j) for i, j in combinations(range(len(sequence)), 2)
                     if sequence[i] + sequence[j] in {"AU", "UA", "CG", "GC", "GU", "UG"}
                     and (strands[i] != strands[j] or j - i - 1 >= minimum)]
            expected = set()
            for count in range(len(sequence) // 2 + 1):
                for selected in combinations(edges, count):
                    if len({index for pair in selected for index in pair}) != 2 * count:
                        continue
                    if any(i < k < j < l or k < i < l < j
                           for (i, j), (k, l) in combinations(selected, 2)):
                        continue
                    symbols = ["."] * len(sequence)
                    for i, j in selected:
                        symbols[i], symbols[j] = "(", ")"
                    expected.add("".join(symbols))
            actual = list(enum_structures(sequence, cuts, minimum))
            with self.subTest(sequence=sequence, cuts=cuts, minimum=minimum):
                self.assertEqual(len(actual), len(set(actual)))
                self.assertEqual(set(actual), expected)

    def test_invalid_sequences_loop_sizes_and_cutpoints(self):
        for sequence in ("", "N", "AC&G", "A C", 123, None):
            with self.subTest(sequence=sequence), self.assertRaises(ValueError):
                list(enum_structures(sequence))
        for minimum in (-1, True, 3.0, "3"):
            with self.subTest(minimum=minimum), self.assertRaises(ValueError):
                list(enum_structures("GAAAC", min_loop_size=minimum))
        for cuts in ((0,), (4,), (2, 2), (3, 1), (1.5,), (True,), [2]):
            with self.subTest(cuts=cuts), self.assertRaises(ValueError):
                list(enum_structures("GCGC", cuts))


class StructureTests(unittest.TestCase):
    def test_connectivity_uses_the_complete_strand_graph(self):
        self.assertTrue(connected("....", ()))
        self.assertFalse(connected("....", (2,)))
        self.assertTrue(connected("(())", (2,)))
        self.assertTrue(connected("()()", (1, 3)))  # strand 1 -- 2 -- 3
        self.assertFalse(connected("(())", (1, 3)))  # strand 2 is disconnected
        self.assertFalse(connected("()()", (1, 2, 3)))  # two separate dimers

    def test_pair_matrix_has_symmetric_pairs_and_unpaired_diagonal(self):
        matrix = pair_matrix("(.())")
        self.assertEqual(matrix[0][4], 1.0)
        self.assertEqual(matrix[4][0], 1.0)
        self.assertEqual(matrix[2][3], 1.0)
        self.assertEqual(matrix[1][1], 1.0)
        for i, row in enumerate(matrix):
            self.assertEqual(sum(row), 1.0)
            for j in range(len(matrix)):
                self.assertEqual(matrix[i][j], matrix[j][i])

    def test_invalid_dot_brackets_and_boundaries_are_rejected(self):
        for structure in ("", "(()", "())", ")(", "([)]", "..&..", "..x..", None):
            with self.subTest(structure=structure), self.assertRaises(ValueError):
                pair_matrix(structure)
            with self.subTest(structure=structure), self.assertRaises(ValueError):
                connected(structure, ())
        with self.assertRaises(ValueError):
            connected("....", (4,))


class BoltzmannTests(unittest.TestCase):
    def test_single_state_and_equal_energy_degeneracy(self):
        single = boltzmann_summary(["(...)"], [-2.0], 0.6)
        self.assertEqual(single["free_energy"], -2.0)
        self.assertEqual(single["pair_probabilities"], pair_matrix("(...)"))
        summary = boltzmann_summary(["....", "()..", "..()"], [2.0] * 3, 0.5)
        self.assertEqual(summary["n_structures"], 3)
        self.assertAlmostEqual(summary["free_energy"], 2.0 - 0.5 * math.log(3))
        self.assertAlmostEqual(summary["pair_probabilities"][0][0], 2 / 3)
        self.assertAlmostEqual(summary["pair_probabilities"][0][1], 1 / 3)

    def test_known_two_state_partition_and_generator_inputs(self):
        summary = boltzmann_summary(iter(["..", "()"]), iter([0.0, math.log(3)]), 1.0)
        self.assertAlmostEqual(summary["free_energy"], -math.log(4 / 3))
        self.assertEqual(summary["pair_probabilities"], [[0.75, 0.25], [0.25, 0.75]])

    def test_shift_invariance_and_overflow_resistance(self):
        states, energies = ["....", "()..", "..()"], [-4.0, -2.0, 1.0]
        reference = boltzmann_summary(states, energies, 0.5)
        for shift in (-10000.0, 10000.0):
            shifted = boltzmann_summary(states, [value + shift for value in energies], 0.5)
            self.assertAlmostEqual(shifted["free_energy"], reference["free_energy"] + shift)
            self.assertEqual(shifted["pair_probabilities"], reference["pair_probabilities"])
        extreme = boltzmann_summary(["..", "()"], [1e308, -1e308], 0.5)
        self.assertEqual(extreme["free_energy"], -1e308)
        self.assertEqual(extreme["pair_probabilities"], pair_matrix("()"))

    def test_probability_rows_normalize_and_input_order_does_not_matter(self):
        states = list(enum_structures("GCGCUA", (2, 4)))
        energies = [0.1 * index for index in range(len(states))]
        summary = boltzmann_summary(states, energies, 0.616)
        reverse = boltzmann_summary(reversed(states), reversed(energies), 0.616)
        self.assertEqual(summary, reverse)
        matrix = summary["pair_probabilities"]
        for i, row in enumerate(matrix):
            self.assertAlmostEqual(math.fsum(row), 1.0)
            self.assertTrue(all(0 <= value <= 1 for value in row))
            for j in range(len(matrix)):
                self.assertEqual(matrix[i][j], matrix[j][i])

    def test_invalid_ensembles_and_thermodynamic_inputs(self):
        cases = (([], [], 1), ([".."], [], 1), (["..", "..."], [0, 0], 1),
                 (["..", ".."], [0, 0], 1), ([".."], [math.nan], 1),
                 ([".."], [math.inf], 1), ([".."], [-math.inf], 1),
                 ([".."], [True], 1), ([".."], ["0"], 1),
                 ([".."], [0], 0), ([".."], [0], -1),
                 ([".."], [0], math.nan), ([".."], [0], True))
        for states, energies, temperature in cases:
            with self.subTest(states=states, energies=energies, kT=temperature), self.assertRaises(ValueError):
                boltzmann_summary(states, energies, temperature)


if __name__ == "__main__":
    unittest.main()
