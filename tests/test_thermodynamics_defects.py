"""Independent finite-ensemble and analytical mass-action defect checks."""
import math
import unittest
from dataclasses import replace

from aptaswitch_core.thermodynamics.defects import (
    DesignDefect, TubeTarget, complex_defect, tube_defect,
)
from aptaswitch_core.thermodynamics.structures import PairingConstraints, parse_structure
from aptaswitch_core.thermodynamics.tubes import StandardStateConvention, solve_tube
from aptaswitch_core.thermodynamics.types import (
    ComplexSpec, Conditions, ModelProvenance, PairProbabilities, PartitionResult, Strand,
)


CONDITIONS = Conditions()
MODEL = ModelProvenance("finite", "1", "fixture", "fixture:defects",
                        "connected finite states", "fixture-1M")


def partition(spec, log_z=0.0):
    return PartitionResult(spec, CONDITIONS, MODEL, log_z, -log_z, 1.0)


def finite_pairs(p, states):
    n = p.complex_spec.n_nucleotides
    matrix = [[0.0] * n for _ in range(n)]
    lengths = tuple(len(s.sequence) for s in p.complex_spec.strands)
    for structure, probability in states:
        for i, j in enumerate(parse_structure(structure, lengths).partners):
            matrix[i][j] += probability
    return PairProbabilities(p, tuple(tuple(row) for row in matrix))


class ComplexDefectTests(unittest.TestCase):
    def test_average_disagreement_is_not_probability_of_whole_structure(self):
        spec = ComplexSpec((Strand("GAAAC", "rna", "A"),))
        p = partition(spec)
        matrix = finite_pairs(p, [("(...)", .75), (".....", .25)])
        result = complex_defect(matrix, "(...)")
        self.assertEqual(result.per_base, (.25, 0, 0, 0, .25))
        self.assertEqual(result.expected_incorrect_nucleotides, .5)
        self.assertEqual(result.normalized, .1)
        self.assertNotEqual(result.normalized, 1 - .75)

    def test_wrong_partner_counts_even_when_nucleotide_is_paired(self):
        spec = ComplexSpec((Strand("GGCC", "rna", "A"),))
        result = complex_defect(finite_pairs(partition(spec), [("()()", 1)]), "(())")
        self.assertEqual(result.normalized, 1)

    def test_target_cuts_and_connectivity_are_checked(self):
        a, b = Strand("GG", "rna", "A"), Strand("CC", "dna", "B")
        matrix = finite_pairs(partition(ComplexSpec((a, b))), [("((+))", 1)])
        for target in ("...", "(+.))", "..+.."):
            with self.subTest(target=target), self.assertRaises(ValueError):
                complex_defect(matrix, target)


class TubeDefectTests(unittest.TestCase):
    def setUp(self):
        self.a, self.b = Strand("GG", "rna", "A"), Strand("CC", "dna", "B")
        self.sa, self.sb = ComplexSpec((self.a,)), ComplexSpec((self.b,))
        self.ab = ComplexSpec((self.a, self.b))
        # An independently selected equilibrium: A=.6uM, B=1.6uM, AB=.4uM.
        k = .4e-6 / (.6e-6 * 1.6e-6)
        self.pa, self.pb, self.pab = partition(self.sa), partition(self.sb), partition(self.ab, math.log(k))
        self.standard = StandardStateConvention(
            provenance_label="fixture-1M", concentration_molar=1.0,
            association_included=True, cyclic_symmetry_included=True, connected_ensembles=True)
        self.tube = solve_tube((self.a, self.b), (1e-6, 2e-6),
                              (self.pa, self.pb, self.pab), standard_state=self.standard)
        self.targets = (TubeTarget("AB", self.ab, "((+))", 1e-6),
                        TubeTarget("free B", self.sb, "..", 1e-6))
        self.b_pairs = finite_pairs(self.pb, [("..", 1)])
        self.ab_pairs = finite_pairs(self.pab, [("((+))", .75), ("(.+.)", .25)])

    def test_analytical_tube_with_residual_free_strand(self):
        result = tube_defect("ON", self.tube, self.targets, (self.b_pairs, self.ab_pairs))
        self.assertAlmostEqual(result.target_nucleotide_molar / 1e-6, 6)
        self.assertAlmostEqual(result.structural_nucleotide_molar / 1e-6, .2)
        self.assertAlmostEqual(result.shortage_nucleotide_molar / 1e-6, 2.4)
        self.assertAlmostEqual(result.normalized_structural, .2 / 6)
        self.assertAlmostEqual(result.normalized_concentration, 2.4 / 6)
        self.assertAlmostEqual(result.normalized, 2.6 / 6)
        self.assertGreater(result.contributions[1].actual_concentration_molar, 1e-6)
        self.assertEqual(result.contributions[1].shortage_nucleotide_molar, 0)

    def test_empty_target_has_full_shortage_and_no_fabricated_matrix(self):
        empty = partition(self.ab, -math.inf)
        tube = solve_tube((self.a, self.b), (1e-6, 2e-6),
                          (self.pa, self.pb, empty), standard_state=self.standard)
        result = tube_defect("empty", tube, self.targets, (self.b_pairs,))
        self.assertIsNone(result.contributions[0].complex_defect)
        self.assertEqual(result.normalized_structural, 0)
        self.assertAlmostEqual(result.normalized_concentration, 4 / 6)

    def test_perfect_target_zero_defect_remains_well_defined(self):
        tube = solve_tube((self.a,), (1e-6,), (self.pa,), standard_state=self.standard)
        result = tube_defect("OFF", tube, (TubeTarget("A", self.sa, "..", 1e-6),),
                             (finite_pairs(self.pa, [("..", 1)]),))
        self.assertEqual(result.normalized, 0)
        self.assertEqual(result.target_nucleotide_molar, 2e-6)
        self.assertIsNotNone(result.target_nucleotide_molar)

    def test_design_is_mean_of_normalized_tubes_not_molar_sums(self):
        on = tube_defect("ON", self.tube, self.targets, (self.b_pairs, self.ab_pairs))
        off_tube = solve_tube((self.a,), (1.0,), (self.pa,), standard_state=self.standard)
        off = tube_defect("OFF", off_tube, (TubeTarget("A", self.sa, "..", 1.0),),
                          (finite_pairs(self.pa, [("..", 1)]),))
        result = DesignDefect((off, on))
        self.assertAlmostEqual(result.normalized, on.normalized / 2)
        self.assertAlmostEqual(result.structural, on.normalized_structural / 2)
        with self.assertRaises(ValueError):
            DesignDefect((on, on))

    def test_missing_population_matrix_and_unbalanced_targets_fail(self):
        with self.assertRaisesRegex(ValueError, "needs its pair"):
            tube_defect("ON", self.tube, self.targets, (self.b_pairs,))
        with self.assertRaisesRegex(ValueError, "conserve"):
            tube_defect("ON", self.tube, self.targets[:1], (self.ab_pairs,))

    def test_wrong_conditions_or_partition_cannot_silently_change_defect(self):
        wrong = replace(self.pb, conditions=Conditions(celsius=25))
        with self.assertRaisesRegex(ValueError, "match a tube partition"):
            tube_defect("ON", self.tube, self.targets,
                        (finite_pairs(wrong, [("..", 1)]), self.ab_pairs))
        with self.assertRaisesRegex(ValueError, "Duplicate"):
            tube_defect("ON", self.tube, self.targets,
                        (self.b_pairs, self.b_pairs, self.ab_pairs))

    def test_conditional_matrix_is_not_a_full_tube_observable(self):
        conditional = replace(self.pb, constraints=PairingConstraints((0, 1)))
        with self.assertRaisesRegex(ValueError, "match a tube partition"):
            tube_defect("ON", self.tube, self.targets,
                        (finite_pairs(conditional, [("..", 1)]), self.ab_pairs))

    def test_rotated_coordinates_are_not_implicitly_reused(self):
        rotated = TubeTarget("AB", ComplexSpec((self.b, self.a)), "((+))", 1e-6)
        with self.assertRaisesRegex(ValueError, "coordinates"):
            tube_defect("ON", self.tube, (rotated, self.targets[1]), (self.b_pairs, self.ab_pairs))

    def test_zero_negative_nonfinite_and_bool_targets_rejected(self):
        for value in (0, -1, math.nan, math.inf, True):
            with self.subTest(value=value), self.assertRaises(ValueError):
                TubeTarget("A", self.sa, "..", value)


if __name__ == "__main__":
    unittest.main()
