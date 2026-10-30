"""Official OFF ABI checked against independent integer-energy enumeration.

These optional integration tests require explicit APTASWITCH_VIENNA_OFF_LIBRARY
and the isolated ViennaRNA 2.7.2 Python bindings. The application backend itself
does not require Python ViennaRNA. No existing settings or runtime are selected.
"""
import math
import os
import unittest

from aptaswitch_core.thermodynamics.structures import PairingConstraints, parse_structure
from aptaswitch_core.thermodynamics.types import ComplexSpec, Conditions, Strand
from tests.native_reference.vienna_reference import direct_vienna
from tests.native_reference.hybrid_native import HybridNativeGrammar


LIBRARY = os.environ.get("APTASWITCH_VIENNA_OFF_LIBRARY")


@unittest.skipUnless(LIBRARY, "explicit official OFF test library required")
class ViennaOffNativeTests(unittest.TestCase):
    def test_integer_minimum_pair_minima_and_complete_alternatives_match_enumeration(self):
        spec = ComplexSpec((Strand("GGGGAAAACCCC", "rna", "s"),))
        for celsius, sodium in ((25., .12), (37., 1.), (55., .12)):
            conditions = Conditions(celsius=celsius, sodium=sodium)
            reference = direct_vienna(spec, conditions, 0, 91)
            energies = {s: e/100 for s, e in reference["integer_structure_energies_centikcal"].items()}
            minimum = min(energies.values())
            pair_minima = [[math.inf] * spec.n_nucleotides for _ in range(spec.n_nucleotides)]
            for structure, energy in energies.items():
                for i, j in enumerate(parse_structure(structure, (spec.n_nucleotides,)).partners):
                    if j != i:
                        pair_minima[i][j] = min(pair_minima[i][j], energy)
            with self.subTest(celsius=celsius, sodium=sodium), HybridNativeGrammar(
                    LIBRARY, spec, conditions, sodium_correction="published-pair", factory="rna") as native:
                observed, matrix = native.minimum_pair_energies()
                self.assertAlmostEqual(observed, minimum, places=5)
                for expected_row, actual_row in zip(pair_minima, matrix):
                    for expected, actual in zip(expected_row, actual_row):
                        if math.isfinite(expected):
                            self.assertAlmostEqual(actual, expected, places=5)
                        else:
                            self.assertEqual(actual, math.inf)
                observed, structures, truncated = native.minimum_structures(atol_kcal=.021, count_limit=100)
                self.assertAlmostEqual(observed, minimum, places=5)
                self.assertFalse(truncated)
                self.assertEqual({s for _, s in structures}, {s for s, e in energies.items() if e <= minimum+.021})
                for energy, structure in structures:
                    self.assertAlmostEqual(energy, energies[structure], places=5)
                _, limited, truncated = native.minimum_structures(atol_kcal=10., count_limit=2)
                self.assertEqual(len(limited), 2)
                self.assertTrue(truncated)

    def test_fully_specified_state_preserves_continuous_pf_and_seeded_samples(self):
        spec = ComplexSpec((Strand("GCGAAAGC", "rna", "s"),))
        conditions = Conditions(celsius=25., sodium=.12)
        state = parse_structure("((....))", (8,))
        with HybridNativeGrammar(LIBRARY, spec, conditions, sodium_correction="published-pair",
                                 constraints=state, factory="rna") as official, HybridNativeGrammar(
                    LIBRARY, spec, conditions, sodium_correction="published-pair", constraints=state) as legacy:
            self.assertAlmostEqual(official.log_partition(), legacy.log_partition(), places=11)
            self.assertEqual(official.sample(20, seed=417), ["((....))"]*20)
            matrix, _ = official.pairs()
            for i, j in enumerate(state.partners):
                self.assertAlmostEqual(matrix[i][j], 1., places=11)

    def test_partial_constraints_match_legacy_pf_and_preserve_partner_indices(self):
        spec = ComplexSpec((Strand("GGGGAAAACCCC", "rna", "s"),))
        conditions = Conditions(celsius=55., sodium=.12)
        state = PairingConstraints.unconstrained(12).assign(0, 11).assign(4, 4)
        with HybridNativeGrammar(LIBRARY, spec, conditions, sodium_correction="published-pair",
                                 constraints=state, factory="rna") as official, HybridNativeGrammar(
                    LIBRARY, spec, conditions, sodium_correction="published-pair", constraints=state) as legacy:
            self.assertAlmostEqual(official.log_partition(), legacy.log_partition(), places=11)
            matrix, _ = official.pairs()
            self.assertAlmostEqual(matrix[0][11], 1., places=11)
            self.assertAlmostEqual(matrix[4][4], 1., places=11)

    def test_factory_rejects_wrong_chemistry_and_multistrand_inputs(self):
        for strands in ((Strand("GCGAAAGC", "dna", "s"),),
                        (Strand("GCGAAAGC", "rna", "s"), Strand("GCTTT", "dna", "t"))):
            with self.subTest(strands=strands), self.assertRaises(ValueError):
                HybridNativeGrammar(LIBRARY, ComplexSpec(strands), Conditions(sodium=1), factory="rna")


if __name__ == "__main__":
    unittest.main()
