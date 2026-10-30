"""Independent finite-state proofs for the exhaustive reference backend.

All energies below are synthetic; these tests validate combinatorics,
normalization and contracts, not a physical DNA/RNA thermodynamic model.
"""
from collections import Counter
from dataclasses import FrozenInstanceError
from itertools import product
import math
import unittest
from unittest.mock import Mock

from aptaswitch_core.thermodynamics.enumeration import (
    EnumerationLimitError, ExhaustiveBackend, enumerate_structures,
)
from aptaswitch_core.thermodynamics.sampling import ConditionalSampler
from aptaswitch_core.thermodynamics.structures import PairingConstraints, format_structure, parse_structure
from aptaswitch_core.thermodynamics.types import (
    ComplexSpec, Conditions, EmptyEnsembleError, ModelProvenance, Strand,
)


PROVENANCE = ModelProvenance(
    backend="test-exhaustive", version="1", model="explicit-synthetic-energies",
    parameter_digest="fixture:synthetic", ensemble="connected;pseudoknot-free;one state per structure",
    standard_state="energy includes association and rotational symmetry",
)


def spec(*sequences, material="rna"):
    return ComplexSpec(tuple(Strand(s, material, f"s{i}") for i, s in enumerate(sequences)))


def gc_pairs(current, i, j):
    return current.sequence[i] + current.sequence[j] in ("GC", "CG")


def all_pairs(current, i, j):
    return True


def backend(**overrides):
    options = dict(pair_allowed=gc_pairs, min_hairpin_size=3, energy=lambda c, t, s: 0.0,
                   thermal_energy=lambda t: 1.0, provenance=PROVENANCE)
    options.update(overrides)
    return ExhaustiveBackend(**options)


def flat(state):
    return format_structure(state, (len(state.partners),))


def brute_strings(current, pair_allowed, minimum):
    """Independent 3^n word search, stack parser and graph connectivity check.

    It does not use the production interval recursion, parser or connectivity
    function. It is intentionally limited to the tiny fixtures in these tests.
    """
    n = current.n_nucleotides
    strand_ids = [k for k, strand in enumerate(current.strands) for _ in strand.sequence]
    states = {}
    for letters in product(".()", repeat=n):
        stack, pairs, partners = [], [], list(range(n))
        valid = True
        for j, letter in enumerate(letters):
            if letter == "(":
                stack.append(j)
            elif letter == ")":
                if not stack:
                    valid = False
                    break
                i = stack.pop()
                if (not pair_allowed(current, i, j)
                        or (strand_ids[i] == strand_ids[j] and j - i - 1 < minimum)):
                    valid = False
                    break
                pairs.append((i, j))
                partners[i], partners[j] = j, i
        if not valid or stack:
            continue
        components = [{i} for i in range(len(current.strands))]
        for i, j in pairs:
            first = next(c for c in components if strand_ids[i] in c)
            second = next(c for c in components if strand_ids[j] in c)
            if first is not second:
                first.update(second)
                components.remove(second)
        if len(components) == 1:
            states["".join(letters)] = tuple(partners)
    return states


class EnumerationStructureTests(unittest.TestCase):
    def test_unrestricted_monomer_counts_are_motzkin_numbers(self):
        for length, expected in enumerate((1, 2, 4, 9, 21, 51, 127), start=1):
            states = tuple(enumerate_structures(spec("A" * length), pair_allowed=all_pairs,
                                              min_hairpin_size=0))
            self.assertEqual(len(states), expected)
            self.assertEqual(len(set(states)), expected)
            self.assertTrue(all(s.complete for s in states))

    def test_hairpin_separation_does_not_apply_across_a_nick(self):
        self.assertEqual(tuple(map(flat, enumerate_structures(spec("GC"), pair_allowed=gc_pairs,
                                                             min_hairpin_size=3))), ("..",))
        self.assertEqual(tuple(map(flat, enumerate_structures(spec("G", "C"), pair_allowed=gc_pairs,
                                                             min_hairpin_size=100))), ("()",))
        for current, expected in ((spec("GAAAC"), {".....", "(...)"}),
                                  (spec("GAAC"), {"...."})):
            self.assertEqual({flat(s) for s in enumerate_structures(current, pair_allowed=gc_pairs,
                                                                  min_hairpin_size=3)}, expected)
        # A nick elsewhere must not enable a short pair within the first strand.
        current = spec("GC", "GC")
        forced = PairingConstraints.unconstrained(4).assign(0, 1)
        self.assertEqual(tuple(enumerate_structures(current, pair_allowed=gc_pairs,
                                                   min_hairpin_size=3, constraints=forced)), ())

    def test_one_to_five_strands_match_independent_word_enumeration(self):
        fixtures = ((spec("GCGAAAGC"), 5), (spec("GC", "GC"), 3),
                    (spec("GC", "GC", "GC"), 4), (spec(*("GC",) * 4), 5),
                    (spec(*("GC",) * 5), 6))
        for current, count in fixtures:
            with self.subTest(n_strands=len(current.strands)):
                expected = brute_strings(current, gc_pairs, 3)
                actual = tuple(enumerate_structures(current, pair_allowed=gc_pairs, min_hairpin_size=3))
                self.assertEqual(len(actual), count)
                self.assertEqual({flat(s): s.partners for s in actual}, expected)

    def test_all_partial_assignments_match_filtering_independent_full_ensemble(self):
        current = spec("GC", "GC", "GC")
        expected = brute_strings(current, gc_pairs, 3)
        partials = {PairingConstraints.unconstrained(6)}
        # All subsets of every complete matching, including resolved unpaired
        # positions and forced pairs, rather than only left-to-right prefixes.
        for partners in expected.values():
            decisions = [(i, j) for i, j in enumerate(partners) if j >= i]
            for mask in range(1 << len(decisions)):
                state = PairingConstraints.unconstrained(6)
                for k, decision in enumerate(decisions):
                    if mask & (1 << k):
                        state = state.assign(*decision)
                partials.add(state)
        # Also test a geometrically forbidden pair and disconnected complete state.
        partials.add(PairingConstraints.unconstrained(6).assign(0, 1))
        partials.add(PairingConstraints(tuple(range(6))))
        for constraint in partials:
            wanted = {s for s, p in expected.items()
                      if all(j == -1 or p[i] == j for i, j in enumerate(constraint.partners))}
            actual = {flat(s) for s in enumerate_structures(current, pair_allowed=gc_pairs,
                                                           min_hairpin_size=3, constraints=constraint)}
            with self.subTest(constraint=constraint):
                self.assertEqual(actual, wanted)

    def test_callback_defines_chemistry_and_caller_order_is_preserved(self):
        current = ComplexSpec((Strand("T", "dna", "d"), Strand("A", "rna", "r")))
        seen = []
        def mixed_allowed(value, i, j):
            seen.append((value, i, j))
            return value.strands[0].material == "dna" and value.sequence == "TA"
        result = backend(pair_allowed=mixed_allowed, min_hairpin_size=3).evaluate(current, Conditions())
        self.assertEqual(result.n_structures, 1)
        self.assertIs(result.partition.complex_spec, current)
        self.assertEqual(result.mfe()[1], ("(+)",))
        self.assertEqual(seen, [(current, 0, 1)])
        reverse = ComplexSpec(tuple(reversed(current.strands)))
        self.assertTrue(backend(pair_allowed=mixed_allowed).partition(reverse, Conditions()).is_empty)

    def test_length_limit_fails_before_callbacks_without_returning_a_partial_ensemble(self):
        pair, energy, thermal = Mock(return_value=True), Mock(return_value=0), Mock(return_value=1)
        model = backend(pair_allowed=pair, energy=energy, thermal_energy=thermal, max_nucleotides=3)
        with self.assertRaisesRegex(EnumerationLimitError, "at most 3.*received 4"):
            model.evaluate(spec("AAAA"), Conditions())
        pair.assert_not_called()
        energy.assert_not_called()
        thermal.assert_not_called()
        self.assertEqual(backend(pair_allowed=all_pairs, min_hairpin_size=0,
                                 max_nucleotides=4).evaluate(spec("AAAA"), Conditions()).n_structures, 9)

    def test_invalid_enumeration_configuration_is_explicit(self):
        for keyword, bad_values in (("min_hairpin_size", (-1, True, 2.0, None)),
                                    ("max_nucleotides", (0, True, 3.0, None))):
            for bad in bad_values:
                with self.subTest(keyword=keyword, value=bad), self.assertRaises(ValueError):
                    backend(**{keyword: bad})
        for keyword in ("pair_allowed", "energy", "thermal_energy", "provenance"):
            with self.subTest(keyword=keyword), self.assertRaises(ValueError):
                backend(**{keyword: None})
        for constraint in ("....", PairingConstraints.unconstrained(3)):
            with self.assertRaises(ValueError):
                backend().evaluate(spec("GC", "GC"), Conditions(), constraint)
        with self.assertRaises(ValueError):
            backend(pair_allowed=lambda c, i, j: 1).evaluate(spec("G", "C"), Conditions())


class EnumerationEnergyTests(unittest.TestCase):
    def setUp(self):
        self.current = spec("GC", "GC")
        self.conditions = Conditions()
        self.weights = {".().": 1.0, "(..)": 2.0, "(())": 3.0}
        self.model = backend(energy=lambda c, t, s: -math.log(self.weights[flat(s)]))

    def test_known_partition_marginals_and_minimum_total_energy(self):
        result = self.model.evaluate(self.current, self.conditions)
        self.assertAlmostEqual(result.partition.log_partition, math.log(6), places=14)
        self.assertAlmostEqual(result.partition.free_energy_kcal_mol, -math.log(6), places=14)
        self.assertIs(result.partition.provenance, PROVENANCE)
        self.assertIs(result.partition.conditions, self.conditions)
        self.assertIsNone(result.partition.constraints)
        self.assertEqual(result.partition.thermal_energy_kcal_mol, 1)
        expected = ((1 / 6, 0, 0, 5 / 6), (0, 1 / 3, 2 / 3, 0),
                    (0, 2 / 3, 1 / 3, 0), (5 / 6, 0, 0, 1 / 6))
        for row, wanted in zip(result.pairs().matrix, expected):
            for value, target in zip(row, wanted):
                self.assertAlmostEqual(value, target, places=14)
        self.assertEqual(result.mfe(), (-math.log(3), ("((+))",)))
        self.assertEqual(self.model.partition(self.current, self.conditions), result.partition)
        self.assertEqual(self.model.pairs(self.current, self.conditions), result.pairs())
        self.assertEqual(self.model.mfe(self.current, self.conditions), result.mfe())

    def test_constant_association_and_symmetry_are_not_applied_twice(self):
        correction = 1.7 + math.log(self.current.symmetry_factor)
        model = backend(energy=lambda c, t, s: correction - math.log(self.weights[flat(s)]))
        result = model.evaluate(self.current, self.conditions)
        self.assertEqual(self.current.symmetry_factor, 2)
        self.assertAlmostEqual(result.partition.log_partition, math.log(6) - correction, places=14)
        self.assertAlmostEqual(result.mfe()[0], correction - math.log(3), places=14)
        self.assertEqual(result.mfe()[1], ("((+))",))
        for row, wanted in zip(result.pairs().matrix, self.model.pairs(self.current, self.conditions).matrix):
            for p, expected in zip(row, wanted):
                self.assertAlmostEqual(p, expected, places=14)

    def test_each_connected_structure_is_evaluated_once_and_reused(self):
        seen = []
        def energy(current, conditions, state):
            self.assertIs(current, self.current)
            self.assertIs(conditions, self.conditions)
            seen.append(state)
            return -math.log(self.weights[flat(state)])
        result = backend(energy=energy).evaluate(self.current, self.conditions)
        self.assertEqual(len(seen), 3)
        self.assertEqual(len(set(seen)), 3)
        result.pairs(), result.mfe(), result.sample(20, seed=6), result.probabilities
        self.assertEqual(len(seen), 3)
        with self.assertRaises(FrozenInstanceError):
            result.structures = ()
        with self.assertRaises(FrozenInstanceError):
            result.structures[0].partners = ()

    def test_partial_and_complete_conditioning_obeys_partition_ratio(self):
        unconstrained = PairingConstraints.unconstrained(4)
        result = self.model.evaluate(self.current, self.conditions)
        for partner, probability in ((0, 1 / 6), (3, 5 / 6), (1, 0.0), (2, 0.0)):
            constraint = unconstrained.assign(0, partner)
            restricted = self.model.evaluate(self.current, self.conditions, constraint)
            self.assertIs(restricted.partition.constraints, constraint)
            ratio = math.exp(restricted.partition.log_partition - result.partition.log_partition)
            self.assertAlmostEqual(ratio, probability, places=14)
        for structure in self.weights:
            matching = parse_structure(structure, (2, 2))
            restricted = self.model.evaluate(self.current, self.conditions, matching)
            self.assertEqual(restricted.n_structures, 1)
            for i, j in enumerate(matching.partners):
                self.assertEqual(restricted.pairs().matrix[i][j], 1)
        self.assertIsNone(self.model.partition(self.current, self.conditions, unconstrained).constraints)

    def test_enumerated_oracle_drives_the_conditional_sampler_joint_law(self):
        sampler = ConditionalSampler(self.current, self.conditions,
            lambda state: self.model.pairs(self.current, self.conditions, state), seed=17)
        joint = {}
        def visit(state, probability):
            if state.complete:
                joint[flat(state)] = probability
                return
            branch = sampler.conditional(state)
            for partner, conditional_probability in branch.choices:
                visit(state.assign(branch.index, partner), probability * conditional_probability)
        visit(PairingConstraints.unconstrained(4), 1.0)
        self.assertEqual(set(joint), set(self.weights))
        for structure, weight in self.weights.items():
            self.assertAlmostEqual(joint[structure], weight / 6, places=14)
        self.assertEqual(set(sampler.sample(100)), {".(+).", "(.+.)", "((+))"})

    def test_positive_infinity_excludes_states_but_nan_and_negative_infinity_fail(self):
        def selective(c, t, s):
            return 0.0 if flat(s) == "(())" else math.inf
        result = backend(energy=selective).evaluate(self.current, self.conditions)
        self.assertEqual(result.n_structures, 1)
        self.assertEqual(result.probabilities, (1.0,))
        self.assertEqual(result.partition.log_partition, 0)
        for bad in (math.nan, -math.inf, True, "0", None):
            with self.subTest(energy=bad), self.assertRaises(ValueError):
                backend(energy=lambda c, t, s: bad).evaluate(self.current, self.conditions)

    def test_empty_ensembles_do_not_produce_normalized_observables(self):
        for model, current in ((backend(), spec("A", "A")),
                               (backend(energy=lambda c, t, s: math.inf), self.current)):
            result = model.evaluate(current, self.conditions)
            self.assertEqual(result.n_structures, 0)
            self.assertEqual((result.partition.log_partition, result.partition.free_energy_kcal_mol),
                             (-math.inf, math.inf))
            for observable in (result.pairs, result.mfe, lambda: result.sample(1, seed=3),
                               lambda: result.probabilities):
                with self.assertRaises(EmptyEnsembleError):
                    observable()
            self.assertEqual(result.sample(0, seed=3), [])
        oracle = Mock(return_value=0)
        backend(energy=oracle).evaluate(spec("A", "A"), self.conditions)
        oracle.assert_not_called()

    def test_equal_energy_degeneracy_is_retained_without_duplicate_states(self):
        result = backend().evaluate(self.current, self.conditions)
        self.assertEqual(result.n_structures, 3)
        self.assertEqual(result.probabilities, (1 / 3,) * 3)
        self.assertAlmostEqual(result.partition.log_partition, math.log(3), places=14)
        self.assertEqual(set(result.mfe()[1]), {".(+).", "(.+.)", "((+))"})

    def test_logsumexp_handles_large_energies_and_relative_underflow(self):
        for offset in (-10000.0, 10000.0):
            result = backend(energy=lambda c, t, s: offset).evaluate(self.current, self.conditions)
            self.assertAlmostEqual(result.partition.log_partition, -offset + math.log(3), places=11)
            self.assertAlmostEqual(result.partition.free_energy_kcal_mol, offset - math.log(3), places=11)
            self.assertEqual(result.probabilities, (1 / 3,) * 3)
        result = backend(energy=lambda c, t, s: -10000 if flat(s) == "(())" else 10000).evaluate(
            self.current, self.conditions)
        self.assertEqual(result.n_structures, 3)
        self.assertEqual(sum(result.probabilities), 1)
        self.assertEqual(result.sample(100, seed=11), ["((+))"] * 100)
        with self.assertRaises(ArithmeticError):
            backend(energy=lambda c, t, s: 1e308, thermal_energy=lambda t: 1e-300).evaluate(
                self.current, self.conditions)

    def test_thermal_convention_and_invalid_conditions_are_enforced(self):
        conditions = Conditions(celsius=25, sodium=1)
        thermal = Mock(return_value=2.0)
        result = backend(energy=lambda c, t, s: 4.0, thermal_energy=thermal).evaluate(self.current, conditions)
        thermal.assert_called_once_with(conditions)
        self.assertEqual(result.partition.thermal_energy_kcal_mol, 2)
        self.assertAlmostEqual(result.partition.log_partition, math.log(3) - 2, places=14)
        for bad in (0, -1, math.inf, -math.inf, math.nan, True, "1", None):
            with self.subTest(thermal=bad), self.assertRaises(ValueError):
                backend(thermal_energy=lambda t: bad).evaluate(self.current, self.conditions)
        mixed = ComplexSpec((Strand("G", "rna", "r"), Strand("C", "dna", "d")))
        with self.assertRaises(ValueError):
            backend().evaluate(mixed, Conditions(magnesium=.001))
        for conditions in (None, 37, "37"):
            with self.assertRaises(ValueError):
                backend().evaluate(self.current, conditions)

    def test_exact_seeded_sampling_matches_known_binomial_marginals(self):
        result = self.model.evaluate(self.current, self.conditions)
        draws = result.sample(6000, seed=7351)
        self.assertEqual(draws, result.sample(6000, seed=7351))
        self.assertEqual(draws, self.model.sample(self.current, self.conditions, 6000, seed=7351))
        counts = Counter(draws)
        self.assertEqual(set(counts), {".(+).", "(.+.)", "((+))"})
        for structure, weight in self.weights.items():
            label = structure[:2] + "+" + structure[2:]
            probability = weight / 6
            sigma = math.sqrt(6000 * probability * (1 - probability))
            self.assertLess(abs(counts[label] - 6000 * probability), 6 * sigma)
        for count, seed in ((-1, 3), (True, 3), (1.0, 3), (1, None), (1, True), (1, 3.0)):
            with self.subTest(count=count, seed=seed), self.assertRaises(ValueError):
                result.sample(count, seed=seed)


if __name__ == "__main__":
    unittest.main()
