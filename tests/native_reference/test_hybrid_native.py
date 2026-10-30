"""Native physical DP against complete independent loop-oracle enumeration.

Set APTASWITCH_HYBRID_NATIVE_LIBRARY explicitly to the isolated build. These
checks prove model agreement, not empirical accuracy of its mixed-loop means.
"""
from collections import Counter
import math
import os
import unittest

from aptaswitch_core.thermodynamics.enumeration import ExhaustiveBackend
from aptaswitch_core.thermodynamics.hybrid_motifs import HybridMotifModel
from aptaswitch_core.thermodynamics.structures import format_structure, parse_structure
from aptaswitch_core.thermodynamics.types import ComplexSpec, Conditions, EmptyEnsembleError, Strand
from tests.native_reference.hybrid_native import HybridNativeGrammar


LIBRARY = os.environ.get("APTASWITCH_HYBRID_NATIVE_LIBRARY")


def spec(*items):
    return ComplexSpec(tuple(Strand(sequence, material, f"s{i}")
                             for i, (sequence, material) in enumerate(items)))


@unittest.skipUnless(LIBRARY, "explicit isolated native hybrid library required")
class HybridNativeTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.model = HybridMotifModel()
        cls.salt_model = HybridMotifModel(sodium_correction="published-pair")

    def reference(self, s, c, model=None, constraints=None):
        model = model or self.model
        return ExhaustiveBackend(pair_allowed=model.pair_allowed, energy=model.energy,
                                 thermal_energy=model.thermal_energy, provenance=model.provenance,
                                 min_hairpin_size=3).evaluate(s, c, constraints=constraints)

    def compare(self, s, c, model=None, constraints=None):
        model = model or self.model
        reference = self.reference(s, c, model, constraints)
        with HybridNativeGrammar(LIBRARY, s, c, sodium_correction=model.sodium_correction,
                                 constraints=constraints) as native:
            self.assertEqual(native.thermal_energy_kcal_mol, reference.partition.thermal_energy_kcal_mol)
            if reference.partition.is_empty:
                self.assertEqual(native.log_partition(), -math.inf)
                with self.assertRaises(EmptyEnsembleError):
                    native.pairs()
                return
            self.assertAlmostEqual(native.log_partition(), reference.partition.log_partition, places=10)
            minimum, structure = native.minimum()
            self.assertAlmostEqual(minimum, reference.mfe()[0], places=10)
            self.assertIn(structure, reference.mfe()[1])
            actual, adjustment = native.pairs()
            self.assertLessEqual(adjustment, 1e-12)
            for expected, observed in zip(reference.pairs().matrix, actual):
                for a, b in zip(expected, observed):
                    self.assertAlmostEqual(a, b, places=10)

    def test_complete_ensembles_rna_dna_hybrid_and_one_to_five_strands(self):
        cases = [spec(("GCGAAAGC", "rna")), spec(("GCGAAAGC", "dna")),
                 spec(("GCAAA", "rna"), ("TTTGC", "dna")),
                 spec(("TTTGC", "dna"), ("GCAAA", "rna")),
                 spec(("GAAA", "rna"), ("AAAA", "dna"))]
        cases += [spec(*[("GC", "rna" if i % 2 else "dna") for i in range(n)]) for n in (2,3,4,5)]
        for celsius in (25, 37, 55):
            for s in cases:
                with self.subTest(order=s.ordered_key, celsius=celsius):
                    self.compare(s, Conditions(celsius=celsius, sodium=1))

    def test_each_complete_state_matches_energy_including_small_interior_tables(self):
        s, c = spec(("GCAAA", "rna"), ("TTTGC", "dna")), Conditions(sodium=1)
        reference = self.reference(s, c)
        self.assertGreater(len(reference.structures), 50)
        for state, energy in zip(reference.structures, reference.energies_kcal_mol):
            with HybridNativeGrammar(LIBRARY, s, c, constraints=state) as native:
                self.assertAlmostEqual(native.minimum()[0], energy, places=10)
                self.assertAlmostEqual(native.log_partition(), -energy/native.thermal_energy_kcal_mol, places=10)

    def test_sodium_corrects_all_paired_states_once(self):
        cases = [spec(("GCGAAAGC", "rna")), spec(("GCGAAAGC", "dna")),
                 spec(("GCAAA", "rna"), ("TTTGC", "dna")),
                 spec(("GC", "rna"), ("GC", "dna"), ("GC", "rna"))]
        for celsius in (25, 37, 55):
            for sodium in (.12, .5, 1.):
                for s in cases:
                    with self.subTest(order=s.ordered_key, celsius=celsius, sodium=sodium):
                        self.compare(s, Conditions(celsius=celsius, sodium=sodium), self.salt_model)

    def test_partial_constraints_match_conditioned_oracle(self):
        s, c = spec(("GCAAA", "rna"), ("TTTGC", "dna")), Conditions(sodium=1)
        reference = self.reference(s, c)
        first = reference.structures[-1]
        from aptaswitch_core.thermodynamics.structures import PairingConstraints
        for i in (0, 2, 4, 9):
            partial = PairingConstraints.unconstrained(s.n_nucleotides).assign(i, first.partners[i])
            self.compare(s, c, constraints=partial)

    def test_mixed_multiloop_with_nicks_inside_children_matches_oracle(self):
        sequence = list("A"*18)
        for i in (0, 1, 9):
            sequence[i] = "G"
        for i in (8, 16, 17):
            sequence[i] = "C"
        sequence = "".join(sequence)
        s = spec((sequence[:6], "rna"), (sequence[6:12], "dna"), (sequence[12:], "rna"))
        state = parse_structure("((......)(......))", (6,6,6))
        for celsius in (25,37,55):
            self.compare(s, Conditions(celsius=celsius, sodium=.12), self.salt_model, state)

    def test_seeded_samples_follow_complete_joint_distribution(self):
        s, c = spec(("GCAAA", "rna"), ("TTTGC", "dna")), Conditions(sodium=1)
        reference = self.reference(s, c)
        with HybridNativeGrammar(LIBRARY, s, c) as native:
            samples = native.sample(8000, seed=7421)
            self.assertEqual(native.sample(20, seed=7421), samples[:20])
        observed = Counter(samples)
        expected = dict(zip((format_structure(state, (5,5)) for state in reference.structures), reference.probabilities))
        self.assertTrue(set(observed) <= set(expected))
        self.assertLess(sum(abs(observed[structure]/8000-p) for structure,p in expected.items()) / 2, .04)

    def test_invalid_domain_and_memory_budget_do_not_produce_results(self):
        s = spec(("GC", "rna"), ("GC", "dna"))
        for c in (Conditions(sodium=.12), Conditions(sodium=1, magnesium=.001),
                  Conditions(celsius=60, sodium=1)):
            with self.subTest(c=c), self.assertRaises(ValueError):
                HybridNativeGrammar(LIBRARY, s, c)
        with self.assertRaises(MemoryError):
            HybridNativeGrammar(LIBRARY, s, Conditions(sodium=1), memory_limit_bytes=8)


if __name__ == "__main__":
    unittest.main()
