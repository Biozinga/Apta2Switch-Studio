"""Numerical grammar proof with synthetic, deliberately asymmetric chemistry."""
from collections import Counter
from functools import lru_cache
from itertools import combinations, product
import math
import unittest

from aptaswitch_core.thermodynamics.enumeration import ExhaustiveBackend, enumerate_structures
from aptaswitch_core.thermodynamics.loops import decompose_structure
from aptaswitch_core.thermodynamics.structures import PairingConstraints, format_structure, parse_structure
from aptaswitch_core.thermodynamics.types import ComplexSpec, Conditions, EmptyEnsembleError, ModelProvenance, Strand
from tests.native_reference.connected_grammar import ConnectedGrammar, DNA, RNA, MIXED


PROVENANCE = ModelProvenance(backend="test-grammar", version="1", model="synthetic-asymmetric-faces",
    parameter_digest="fixture:explicit-coefficients-in-test", ensemble="connected;pseudoknot-free;dangles=0",
    standard_state="synthetic association;uniform cyclic factor in labelled state energies")


def make_spec(lengths, *, mixed=False, material="rna"):
    return ComplexSpec(tuple(Strand("A" * n, "dna" if mixed and k % 2 else material, str(k))
                             for k, n in enumerate(lengths)))


def material_mask(materials):
    value = 0
    for material in materials:
        value |= DNA if material == "dna" else RNA
    return value


class UnitWeights:
    def pair_allowed(self, spec, i, j):
        return True
    def hairpin(self, spec, i, j):
        return 1.0
    def interior(self, spec, i, j, k, l):
        return 1.0
    def unpaired(self, spec, mode, family, i):
        return 1.0
    def branch(self, spec, mode, family, i, j):
        return 1.0
    def closing(self, spec, mode, family, i, j):
        return 1.0


class AsymmetricWeights:
    """An explicit artificial model; a whole-complex material average is wrong."""
    bias = {DNA: .23, RNA: .79, MIXED: 1.67}

    def pair_allowed(self, spec, i, j):
        # A-A pairs are intentional for this topology/weight proof. This class
        # is not a physical DNA, RNA or hybrid nearest-neighbor parameter set.
        return True

    @staticmethod
    def _materials(spec):
        return tuple(s.material for s in spec.strands for _ in s.sequence)

    def _mask_at(self, spec, indices):
        materials = self._materials(spec)
        return material_mask(materials[i] for i in indices)

    def hairpin_energy(self, spec, i, j):
        family = self._mask_at(spec, range(i, j + 1))
        return .19 + .07 * (j - i - 1) + .011 * i + .017 * j + self.bias[family]

    def interior_energy(self, spec, i, j, k, l):
        family = self._mask_at(spec, tuple(range(i, k + 1)) + tuple(range(l, j + 1)))
        return (.23 + .03 * (k - i - 1) + .07 * (j - l - 1)
                + .009 * (i + 2 * k + 3 * l + 5 * j) + self.bias[family])

    def unpaired_energy(self, spec, mode, family, i):
        local = self._mask_at(spec, (i,))
        return .013 * (i + 1) + .021 * local + (.041 if mode == "multi" else .019) * self.bias[family]

    def branch_energy(self, spec, mode, family, i, j):
        local = self._mask_at(spec, (i, j))
        return .09 + .011 * (i + 1) + .017 * (j + 1) + .13 * local + (
            .29 if mode == "multi" else .11) * self.bias[family]

    def closing_energy(self, spec, mode, family, i, j):
        return self.branch_energy(spec, mode, family, j, i) + (
            .57 + .19 * self.bias[family] if mode == "multi" else 0)

    def hairpin(self, *args):
        return math.exp(-self.hairpin_energy(*args))
    def interior(self, *args):
        return math.exp(-self.interior_energy(*args))
    def unpaired(self, *args):
        return math.exp(-self.unpaired_energy(*args))
    def branch(self, *args):
        return math.exp(-self.branch_energy(*args))
    def closing(self, *args):
        return math.exp(-self.closing_energy(*args))

    @staticmethod
    def global_energy(spec):
        return .27 * (len(spec.strands) - 1) + math.log(spec.symmetry_factor)

    def energy(self, spec, conditions, state, *, forced_family=None):
        # This evaluator traverses complete faces from the independent loop
        # decomposition; it does not call the DP's productions or forests.
        lengths = tuple(len(s.sequence) for s in spec.strands)
        loops = decompose_structure(spec, format_structure(state, lengths))
        values = [self.global_energy(spec)]
        for face in loops:
            if face.kind == "hairpin":
                values.append(self.hairpin_energy(spec, *face.closing_pair))
            elif face.kind in ("stack", "bulge", "interior"):
                if face.unpaired_count > 30:
                    return math.inf
                values.append(self.interior_energy(spec, *face.closing_pair, *face.child_pairs[0]))
            else:
                mode = "multi" if face.kind == "multiloop" else "exterior"
                family = forced_family if forced_family is not None else material_mask(face.materials)
                if face.closing_pair is not None:
                    values.append(self.closing_energy(spec, mode, family, *face.closing_pair))
                values.extend(self.branch_energy(spec, mode, family, i, j) for i, j in face.child_pairs)
                values.extend(self.unpaired_energy(spec, mode, family, i)
                              for a, b in face.unpaired_runs for i in range(a, b))
        return math.fsum(values)


class ConnectedGrammarTests(unittest.TestCase):
    def compare(self, spec, *, minimum=0, constraints=None, model=None):
        model = model or AsymmetricWeights()
        reference = ExhaustiveBackend(pair_allowed=model.pair_allowed, energy=model.energy,
            thermal_energy=lambda conditions: 1.0, provenance=PROVENANCE,
            min_hairpin_size=minimum, max_nucleotides=80)
        result = reference.evaluate(spec, Conditions(), constraints)
        grammar = ConnectedGrammar(spec, model, min_hairpin_size=minimum, constraints=constraints,
                                   global_factor=math.exp(-model.global_energy(spec)))
        actual = grammar.partition_function()
        expected = result.partition.partition_function
        self.assertAlmostEqual(actual, expected, delta=max(1e-14, abs(expected) * 2e-12))
        if actual:
            outside = grammar.pairs()
            self.assertLessEqual(outside.maximum_endpoint_adjustment, 1e-12)
            self.assert_matrix_close(outside.matrix, result.pairs().matrix)
            minimum, structure = grammar.minimum()
            self.assertAlmostEqual(minimum, result.mfe()[0], delta=2e-12)
            selected = parse_structure(structure, tuple(len(s.sequence) for s in spec.strands))
            self.assertIn(selected, result.structures)
            self.assertAlmostEqual(model.energy(spec, Conditions(), selected), minimum, delta=2e-12)
        else:
            for observable in (grammar.pairs, grammar.minimum, lambda: grammar.sample(1, seed=3)):
                with self.assertRaises(EmptyEnsembleError):
                    observable()
        return actual, result

    def assert_matrix_close(self, actual, expected):
        self.assertEqual(len(actual), len(expected))
        for i, (row, wanted) in enumerate(zip(actual, expected)):
            self.assertAlmostEqual(math.fsum(row), 1, delta=1e-12)
            for j, (p, q) in enumerate(zip(row, wanted)):
                self.assertEqual(p, actual[j][i])
                self.assertAlmostEqual(p, q, delta=2e-12)

    def test_unit_weights_count_every_small_ordered_connected_ensemble(self):
        checked = 0
        for n in range(1, 8):
            for cut_count in range(min(4, n - 1) + 1):
                for cuts in combinations(range(1, n), cut_count):
                    ends = (0,) + cuts + (n,)
                    spec = make_spec(tuple(b - a for a, b in zip(ends, ends[1:])), mixed=True)
                    expected = len(tuple(enumerate_structures(spec, pair_allowed=lambda c, i, j: True,
                                                             min_hairpin_size=0)))
                    actual = ConnectedGrammar(spec, UnitWeights(), min_hairpin_size=0).partition_function()
                    self.assertEqual(actual, expected, (n, cuts))
                    checked += 1
        self.assertEqual(checked, 119)

    def test_asymmetric_weights_match_exhaustive_homogeneous_and_mixed_ensembles(self):
        for spec in (make_spec((8,)), make_spec((8,), material="dna"),
                     make_spec((4, 4), mixed=True), make_spec((2, 2, 2), mixed=True),
                     make_spec((2,) * 4, mixed=True), make_spec((2,) * 5, mixed=True)):
            with self.subTest(strands=spec.strands):
                self.compare(spec)

    def test_all_partial_constraints_match_their_exhaustive_subsets(self):
        spec = make_spec((2, 2, 2), mixed=True)
        states = tuple(enumerate_structures(spec, pair_allowed=lambda c, i, j: True, min_hairpin_size=0))
        constraints = {PairingConstraints.unconstrained(6), PairingConstraints(tuple(range(6)))}
        for state in states:
            decisions = [(i, j) for i, j in enumerate(state.partners) if j >= i]
            for mask in range(1 << len(decisions)):
                partial = PairingConstraints.unconstrained(6)
                for k, decision in enumerate(decisions):
                    if mask & (1 << k):
                        partial = partial.assign(*decision)
                constraints.add(partial)
        for partial in constraints:
            with self.subTest(partial=partial):
                self.compare(spec, constraints=partial)

    def test_hairpin_minimum_and_nick_exception_with_complete_constraints(self):
        for unpaired in range(5):
            word = "(" + "." * unpaired + ")"
            state = parse_structure(word, (len(word),))
            for spec in (make_spec((len(word),)), make_spec((1, len(word) - 1), mixed=True)):
                actual, result = self.compare(spec, minimum=3, constraints=state)
                self.assertEqual(actual > 0, len(spec.strands) > 1 or unpaired >= 3)

    def test_interior_limit_and_nicks_on_parent_or_child_faces(self):
        for total in (0, 1, 29, 30, 31):
            a, b = total // 2, total - total // 2
            word = "(" + "." * a + "()" + "." * b + ")"
            n, child_left = len(word), a + 1
            state = parse_structure(word, (n,))
            # A nick on the outer face makes it exterior, so the interior
            # size cap no longer applies. A nick inside the child does not.
            choices = ((), (1,), (child_left + 1,), tuple(sorted({1, child_left + 1})))
            for cuts in choices:
                ends = (0,) + cuts + (n,)
                spec = make_spec(tuple(y - x for x, y in zip(ends, ends[1:])), mixed=True)
                actual, result = self.compare(spec, constraints=state)
                self.assertEqual(actual > 0, total <= 30 or 1 in cuts)

    def test_child_interior_material_does_not_contaminate_parent_face(self):
        spec = make_spec((5, 4, 16), mixed=True)
        partners = list(range(25))
        for i, j in ((0, 24), (2, 16), (4, 8), (18, 22)):
            partners[i], partners[j] = j, i
        state = PairingConstraints(tuple(partners))
        actual, result = self.compare(spec, constraints=state)
        model = AsymmetricWeights()
        contaminated = math.exp(-model.energy(spec, Conditions(), state, forced_family=MIXED))
        self.assertGreater(abs(actual - contaminated), actual * .01)
        loops = decompose_structure(spec, format_structure(state, (5, 4, 16)))
        parent = next(face for face in loops if face.closing_pair == (0, 24))
        self.assertEqual(parent.materials, {"rna"})

    def test_global_correction_is_applied_once_and_root_physical_nicks_are_forbidden(self):
        spec = make_spec((2, 2))
        raw = ConnectedGrammar(spec, UnitWeights(), min_hairpin_size=0).partition_function()
        corrected = ConnectedGrammar(spec, UnitWeights(), min_hairpin_size=0,
                                     global_factor=.125).partition_function()
        self.assertEqual(corrected, raw * .125)
        disconnected = PairingConstraints(tuple(range(4)))
        self.assertEqual(ConnectedGrammar(spec, UnitWeights(), min_hairpin_size=0,
                                          constraints=disconnected).partition_function(), 0)

    def test_complete_traceback_joint_law_has_one_derivation_per_structure(self):
        for spec, constraint in ((make_spec((3, 3), mixed=True), None),
                                 (make_spec((2, 2, 2), mixed=True), None),
                                 (make_spec((2,) * 4, mixed=True), None),
                                 (make_spec((2,) * 5, mixed=True), None),
                                 (make_spec((2, 2, 2), mixed=True),
                                  PairingConstraints.unconstrained(6).assign(0, 5))):
            model = AsymmetricWeights()
            grammar = ConnectedGrammar(spec, model, min_hairpin_size=0, constraints=constraint,
                                       global_factor=math.exp(-model.global_energy(spec)))
            _, reference = self.compare(spec, constraints=constraint, model=model)
            @lru_cache(maxsize=None)
            def derivations(node):
                results = []
                own_pair = (node[1:],) if node[0] == "B" else ()
                denominator = grammar._inside(node)
                for coefficient, children in grammar._productions(node):
                    contribution = coefficient * math.prod(grammar._inside(child) for child in children)
                    if contribution == 0:
                        continue
                    for subtrees in product(*(derivations(child) for child in children)):
                        pairs = own_pair + tuple(pair for subtree, _ in subtrees for pair in subtree)
                        probability = contribution / denominator * math.prod(p for _, p in subtrees)
                        results.append((tuple(sorted(pairs)), probability))
                return tuple(results)
            law = derivations(("R",))
            counts = Counter(pairs for pairs, probability in law)
            self.assertEqual(set(counts.values()), {1})
            probabilities = dict(law)
            self.assertEqual(len(probabilities), reference.n_structures)
            self.assertAlmostEqual(math.fsum(probabilities.values()), 1, delta=2e-12)
            for state, expected in zip(reference.structures, reference.probabilities):
                pairs = tuple((i, j) for i, j in enumerate(state.partners) if i < j)
                self.assertAlmostEqual(probabilities[pairs], expected, delta=2e-12)

    def test_weighted_tracebacks_are_seeded_and_match_binomial_frequencies(self):
        spec, model = make_spec((2, 2), mixed=True), AsymmetricWeights()
        _, reference = self.compare(spec, model=model)
        grammar = ConnectedGrammar(spec, model, min_hairpin_size=0,
                                   global_factor=math.exp(-model.global_energy(spec)))
        draws = grammar.sample(6000, seed=98341)
        self.assertEqual(draws, grammar.sample(6000, seed=98341))
        counts = Counter(draws)
        expected_labels = {format_structure(state, (2, 2)) for state in reference.structures}
        self.assertEqual(set(counts), expected_labels)
        for state, p in zip(reference.structures, reference.probabilities):
            deviation = abs(counts[format_structure(state, (2, 2))] - 6000 * p)
            self.assertLess(deviation, 6 * math.sqrt(6000 * p * (1 - p)))
        for state in reference.structures:
            forced = ConnectedGrammar(spec, model, min_hairpin_size=0, constraints=state)
            self.assertEqual(forced.sample(3, seed=22), [format_structure(state, (2, 2))] * 3)
        for count, seed in ((-1, 0), (True, 0), (1, True), (1, None)):
            with self.assertRaises(ValueError):
                grammar.sample(count, seed=seed)

    def test_outside_is_invariant_to_global_factor_and_minimum_uses_supplied_rt(self):
        spec, model = make_spec((2, 2, 2), mixed=True), AsymmetricWeights()
        base = ConnectedGrammar(spec, model, min_hairpin_size=0)
        shifted = ConnectedGrammar(spec, model, min_hairpin_size=0, global_factor=.123)
        self.assert_matrix_close(base.pairs().matrix, shifted.pairs().matrix)
        energy, structure = base.minimum()
        corrected, other_structure = shifted.minimum()
        self.assertEqual(structure, other_structure)
        self.assertAlmostEqual(corrected, energy - math.log(.123), delta=2e-12)
        scaled, scaled_structure = base.minimum(thermal_energy_kcal_mol=.616)
        self.assertAlmostEqual(scaled, .616 * energy, delta=2e-12)
        self.assertEqual(scaled_structure, structure)
        for thermal in (0, -1, math.nan, math.inf, True):
            with self.assertRaises(ValueError):
                base.minimum(thermal_energy_kcal_mol=thermal)

    def test_zero_and_tied_weights_have_consistent_all_observables(self):
        spec = make_spec((2, 2))
        grammar = ConnectedGrammar(spec, UnitWeights(), min_hairpin_size=0)
        self.assertEqual(grammar.minimum()[0], 0)
        self.assertIn(grammar.minimum()[1], grammar.sample(200, seed=94))
        zero = ConnectedGrammar(spec, UnitWeights(), min_hairpin_size=0, global_factor=0)
        self.assertEqual(zero.partition_function(), 0)
        self.assertEqual(zero.sample(0, seed=12), [])
        for observable in (zero.pairs, zero.minimum, lambda: zero.sample(1, seed=4)):
            with self.assertRaises(EmptyEnsembleError):
                observable()


if __name__ == "__main__":
    unittest.main()
