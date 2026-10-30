"""Independent finite joint laws check sampling, not a folding implementation."""
from collections import Counter
from dataclasses import replace
import math
import unittest

from aptaswitch_core.thermodynamics.sampling import ConditionalSampler, ProbabilityConsistencyError
from aptaswitch_core.thermodynamics.structures import (
    PairingConstraints, parse_structure, format_structure, is_connected,
)
from aptaswitch_core.thermodynamics.types import (
    Strand, ComplexSpec, Conditions, ModelProvenance, PartitionResult, PairProbabilities,
)


class FiniteOracle:
    """Exact conditional marginals of an explicitly supplied joint law."""

    def __init__(self, complex_spec, conditions, structures):
        self.complex_spec, self.conditions = complex_spec, conditions
        self.lengths = tuple(len(s.sequence) for s in complex_spec.strands)
        self.entries = [(parse_structure(s, self.lengths), w) for s, w in structures.items()]
        self.provenance = ModelProvenance(backend="test-enumeration", version="1", model="finite-law",
            parameter_digest="fixture:synthetic", ensemble="connected-pseudoknot-free", standard_state="1 M")

    def __call__(self, constraints):
        entries = [(s, w) for s, w in self.entries if all(
            j == -1 or s.partners[i] == j for i, j in enumerate(constraints.partners))]
        z = math.fsum(w for _, w in entries)
        if z <= 0:
            raise AssertionError("Sampler reached a zero-probability state")
        n = sum(self.lengths)
        matrix = [[0.0] * n for _ in range(n)]
        for state, weight in entries:
            for i, j in enumerate(state.partners):
                matrix[i][j] += weight / z
        partition = PartitionResult(complex_spec=self.complex_spec, conditions=self.conditions,
            provenance=self.provenance, thermal_energy_kcal_mol=1,
            log_partition=math.log(z), free_energy_kcal_mol=-math.log(z), constraints=constraints)
        return PairProbabilities(partition=partition, matrix=tuple(tuple(row) for row in matrix))


class SamplingTests(unittest.TestCase):
    def setUp(self):
        self.spec = ComplexSpec((Strand("GC", "rna", "A"), Strand("GC", "rna", "B")))
        self.conditions = Conditions()
        self.expected = {"((+))": .5, "(.+.)": .3, ".(+).": .2}
        self.oracle = FiniteOracle(self.spec, self.conditions, self.expected)

    def sampler(self, **kwargs):
        return ConditionalSampler(self.spec, self.conditions, self.oracle, seed=27, **kwargs)

    def test_complete_joint_law_matches_exactly(self):
        sampler = self.sampler()
        observed = Counter()
        def visit(state, probability):
            if state.complete:
                observed[format_structure(state, (2, 2))] += probability
                return
            branch = sampler.conditional(state)
            for j, p in branch.choices:
                visit(state.assign(branch.index, j), probability * p)
        visit(PairingConstraints.unconstrained(4), 1)
        self.assertEqual(set(observed), set(self.expected))
        for s, p in self.expected.items():
            self.assertAlmostEqual(observed[s], p, places=14)

    def test_seed_and_cache_eviction_do_not_change_draws(self):
        a, b = self.sampler(cache_size=0), self.sampler(cache_size=2)
        samples = a.sample(200)
        self.assertEqual(samples, b.sample(200))
        self.assertEqual(set(samples), set(self.expected))
        self.assertLessEqual(len(b.cache), 2)

    def test_unconditional_rows_are_rejected(self):
        always = self.oracle(PairingConstraints.unconstrained(4))
        sampler = ConditionalSampler(self.spec, self.conditions, lambda state: always, seed=0)
        with self.assertRaises(ProbabilityConsistencyError):
            sampler.conditional(PairingConstraints.unconstrained(4).assign(0, 3))

    def test_correct_resolved_rows_do_not_excuse_a_wrong_constraint_identity(self):
        state = PairingConstraints.unconstrained(4).assign(0, 3)
        correct = self.oracle(state)
        mislabeled = replace(correct, partition=replace(correct.partition, constraints=None))
        self.assertEqual(mislabeled.matrix[0][3], 1)
        sampler = ConditionalSampler(self.spec, self.conditions, lambda state: mislabeled, seed=0)
        with self.assertRaisesRegex(ProbabilityConsistencyError, "constrained ensemble"):
            sampler.conditional(state)

    def test_correct_constraint_label_does_not_excuse_unconditional_rows(self):
        state = PairingConstraints.unconstrained(4).assign(0, 3)
        always = self.oracle(PairingConstraints.unconstrained(4))
        mislabeled = replace(always, partition=replace(always.partition, constraints=state))
        sampler = ConditionalSampler(self.spec, self.conditions, lambda state: mislabeled, seed=0)
        with self.assertRaisesRegex(ProbabilityConsistencyError, "not enforced"):
            sampler.conditional(state)

    def test_disconnected_oracle_is_rejected(self):
        oracle = FiniteOracle(self.spec, self.conditions, {"..+..": 1})
        sampler = ConditionalSampler(self.spec, self.conditions, oracle, seed=0)
        with self.assertRaisesRegex(ProbabilityConsistencyError, "disconnected"):
            sampler.sample_one()

    def test_model_change_is_rejected(self):
        sampler = self.sampler(cache_size=0)
        sampler.conditional(PairingConstraints.unconstrained(4))
        self.oracle.provenance = replace(self.oracle.provenance, parameter_digest="fixture:changed")
        with self.assertRaisesRegex(ProbabilityConsistencyError, "energy model"):
            sampler.conditional(PairingConstraints.unconstrained(4))

    def test_five_strand_chain_preserves_cuts(self):
        spec = ComplexSpec(tuple(Strand("GC", "rna", str(i)) for i in range(5)))
        state = parse_structure(".()()()().", (2,) * 5)
        structure = format_structure(state, (2,) * 5)
        self.assertTrue(is_connected(state, (2,) * 5))
        sampler = ConditionalSampler(spec, self.conditions,
            FiniteOracle(spec, self.conditions, {structure: 1}), seed=0)
        self.assertEqual(sampler.sample(3), [structure] * 3)

    def test_zero_samples_do_not_call_oracle(self):
        sampler = self.sampler()
        self.assertEqual(sampler.sample(0), [])
        self.assertEqual(sampler.oracle_calls, 0)
        for bad in (-1, True, 1.5):
            with self.assertRaises(ValueError):
                sampler.sample(bad)


class StructureTests(unittest.TestCase):
    def test_round_trip_and_boundaries(self):
        state = parse_structure("((+))", (2, 2))
        self.assertEqual(state, parse_structure("(())", (2, 2)))
        self.assertEqual(format_structure(state, (2, 2)), "((+))")
        for invalid in ("(+)))", "((+)", "((&&))", "([+])", "((+&))"):
            with self.assertRaises(ValueError):
                parse_structure(invalid, (2, 2))

    def test_partial_matching_validation(self):
        for invalid in ((2, 3, 0, 1), (1, -1), (True,), (-2,), (2,)):
            with self.assertRaises(ValueError):
                PairingConstraints(invalid)
        with self.assertRaises(ValueError):
            PairingConstraints.unconstrained(3).assign(0, 0).assign(0, 1)


if __name__ == "__main__":
    unittest.main()
