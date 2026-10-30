"""Optional native proofs; set APTASWITCH_DP_FIXTURE_LIBRARY to a built library."""
from collections import Counter
import ctypes as C
from itertools import combinations
import math
import os
import unittest

from aptaswitch_core.thermodynamics.enumeration import ExhaustiveBackend, enumerate_structures
from aptaswitch_core.thermodynamics.structures import PairingConstraints, format_structure, parse_structure
from aptaswitch_core.thermodynamics.types import Conditions, EmptyEnsembleError
from tests.native_reference.connected_grammar import ConnectedGrammar
from tests.native_reference.connected_native import FixtureNativeGrammar, NativeInput, NativeModel
from tests.native_reference.test_connected_grammar import AsymmetricWeights, UnitWeights, make_spec, PROVENANCE

LIBRARY = os.environ.get("APTASWITCH_DP_FIXTURE_LIBRARY")


@unittest.skipUnless(LIBRARY, "optional native fixture library path is not configured")
class ConnectedNativeTests(unittest.TestCase):
    def native(self, spec, **options):
        return FixtureNativeGrammar(LIBRARY, spec, fixture=options.pop("fixture", "asymmetric"),
                                    min_hairpin_size=options.pop("minimum", 0), **options)

    def compare(self, spec, *, constraints=None, minimum=0, thermal=1.0, pair_allowed=None):
        model = AsymmetricWeights()
        allowed = pair_allowed or model.pair_allowed
        reference = ExhaustiveBackend(pair_allowed=allowed, energy=model.energy,
            thermal_energy=lambda conditions: thermal, provenance=PROVENANCE,
            min_hairpin_size=minimum, max_nucleotides=80).evaluate(spec, Conditions(), constraints)
        with self.native(spec, constraints=constraints, minimum=minimum, thermal_energy_kcal_mol=thermal,
                         pair_allowed=allowed, log_global_factor=-model.global_energy(spec) / thermal) as native:
            if reference.partition.is_empty:
                self.assertEqual(native.log_partition(), -math.inf)
                for operation in (native.pairs, native.minimum, native.minimum_pair_energies, native.minimum_structures,
                                  lambda: native.sample(1, seed=0)):
                    with self.assertRaises(EmptyEnsembleError):
                        operation()
            else:
                self.assertAlmostEqual(native.log_partition(), reference.partition.log_partition, delta=3e-12)
                matrix, adjustment = native.pairs()
                self.assertLessEqual(adjustment, 1e-12)
                self.assert_matrix_close(matrix, reference.pairs().matrix)
                energy, structure = native.minimum()
                self.assertAlmostEqual(energy, reference.mfe()[0], delta=3e-12)
                state = parse_structure(structure, tuple(len(s.sequence) for s in spec.strands))
                self.assertIn(state, reference.structures)
                self.assertAlmostEqual(model.energy(spec, Conditions(), state), energy, delta=3e-12)
                minimum_energy, pair_energies = native.minimum_pair_energies()
                self.assertAlmostEqual(minimum_energy, energy, delta=3e-12)
                for i, row in enumerate(pair_energies):
                    for j, pair_energy in enumerate(row):
                        expected = min((e for state, e in zip(reference.structures, reference.energies_kcal_mol)
                                        if i != j and state.partners[i] == j), default=math.inf)
                        if expected == math.inf:
                            self.assertEqual(pair_energy, math.inf)
                        else:
                            self.assertAlmostEqual(pair_energy, expected, delta=3e-12)
                            self.assertGreaterEqual(pair_energy, energy - 3e-12)
                        self.assertEqual(pair_energy, pair_energies[j][i])
                # The max/outside operation must preserve the log-sum cache.
                self.assertAlmostEqual(native.log_partition(), reference.partition.log_partition, delta=3e-12)
                for structure in native.sample(5, seed=9421):
                    self.assertIn(parse_structure(structure, tuple(len(s.sequence) for s in spec.strands)),
                                  reference.structures)
        return reference

    def assert_matrix_close(self, actual, expected):
        self.assertEqual(len(actual), len(expected))
        for i, (row, wanted) in enumerate(zip(actual, expected)):
            self.assertAlmostEqual(math.fsum(row), 1, delta=1e-12)
            for j, (p, q) in enumerate(zip(row, wanted)):
                self.assertEqual(p, actual[j][i])
                self.assertAlmostEqual(p, q, delta=3e-12)

    def test_all_119_unit_weight_ensembles_match_python_grammar(self):
        checked = 0
        for n in range(1, 8):
            for cut_count in range(min(4, n - 1) + 1):
                for cuts in combinations(range(1, n), cut_count):
                    ends = (0,) + cuts + (n,)
                    spec = make_spec(tuple(b - a for a, b in zip(ends, ends[1:])), mixed=True)
                    expected = ConnectedGrammar(spec, UnitWeights(), min_hairpin_size=0)
                    z = expected.partition_function()
                    with self.native(spec, fixture="unit") as native:
                        self.assertAlmostEqual(native.log_partition(), math.log(z) if z else -math.inf, delta=1e-12)
                        if z:
                            probabilities = expected.pairs().matrix
                            self.assert_matrix_close(native.pairs()[0], probabilities)
                            self.assertEqual(native.minimum()[0], 0)
                            minimum, energies = native.minimum_pair_energies()
                            self.assertEqual(minimum, 0)
                            for i, row in enumerate(energies):
                                for j, energy in enumerate(row):
                                    self.assertEqual(energy, 0 if i != j and probabilities[i][j] > 0 else math.inf)
                    checked += 1
        self.assertEqual(checked, 119)

    def test_nonuniform_rna_dna_and_mixed_weights_at_distinct_thermal_conventions(self):
        for spec in (make_spec((8,)), make_spec((8,), material="dna"),
                     make_spec((4, 4), mixed=True), make_spec((2,) * 3, mixed=True),
                     make_spec((2,) * 4, mixed=True), make_spec((2,) * 5, mixed=True)):
            for thermal in (.5925, .6163207755, .652, 1.0):
                with self.subTest(spec=spec, thermal=thermal):
                    self.compare(spec, thermal=thermal)

    def test_all_partial_trimolecular_constraints_match_exact_subsets(self):
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
            with self.subTest(constraints=partial):
                self.compare(spec, constraints=partial)

    def test_hairpin_minimum_interior_30_and_nick_boundary_cases(self):
        for u in range(5):
            state = parse_structure("(" + "." * u + ")", (u + 2,))
            self.compare(make_spec((u + 2,)), minimum=3, constraints=state)
            self.compare(make_spec((1, u + 1), mixed=True), minimum=3, constraints=state)
        for u in (0, 1, 29, 30, 31):
            a, b = u // 2, u - u // 2
            word = "(" + "." * a + "()" + "." * b + ")"
            n = len(word)
            state = parse_structure(word, (n,))
            for cuts in ((), (1,), (a + 2,), tuple(sorted({1, a + 2}))):
                ends = (0,) + cuts + (n,)
                self.compare(make_spec(tuple(y - x for x, y in zip(ends, ends[1:])), mixed=True), constraints=state)

    def test_child_interior_chemistry_does_not_leak_into_parent_loop(self):
        spec = make_spec((5, 4, 16), mixed=True)
        partners = list(range(25))
        for i, j in ((0, 24), (2, 16), (4, 8), (18, 22)):
            partners[i], partners[j] = j, i
        self.compare(spec, constraints=PairingConstraints(tuple(partners)))

    def test_arbitrary_pair_masks_and_zero_partition(self):
        spec = make_spec((2, 2), mixed=True)
        self.compare(spec, pair_allowed=lambda c, i, j: (i, j) == (0, 3))
        self.compare(spec, pair_allowed=lambda c, i, j: False)
        self.compare(make_spec((5,)), pair_allowed=lambda c, i, j: False)

    def test_log_global_constants_do_not_change_marginals_or_seeded_law(self):
        spec = make_spec((2, 2, 2), mixed=True)
        with self.native(spec) as native:
            baseline, matrix, minimum, samples = native.log_partition(), native.pairs()[0], native.minimum(), native.sample(100, seed=88)
            pair_minimum, pair_energies = native.minimum_pair_energies()
        for constant in (-10000.0, 10000.0, -1e6, 1e6):
            with self.native(spec, log_global_factor=constant) as native:
                self.assertAlmostEqual(native.log_partition(), baseline + constant, delta=1e-10)
                self.assertEqual(native.pairs()[0], matrix)
                self.assertEqual(native.sample(100, seed=88), samples)
                energy, structure = native.minimum()
                self.assertAlmostEqual(energy, minimum[0] - constant, delta=1e-10)
                self.assertEqual(structure, minimum[1])
                shifted_minimum, shifted_pairs = native.minimum_pair_energies()
                self.assertAlmostEqual(shifted_minimum, pair_minimum - constant, delta=2e-10)
                for row, shifted_row in zip(pair_energies, shifted_pairs):
                    for value, shifted in zip(row, shifted_row):
                        if math.isinf(value):
                            self.assertEqual(shifted, math.inf)
                        else:
                            self.assertAlmostEqual(shifted, value - constant, delta=2e-10)
        with self.native(spec, log_global_factor=-math.inf) as native:
            self.assertEqual(native.log_partition(), -math.inf)
            self.assertEqual(native.sample(0, seed=0), [])
            with self.assertRaises(EmptyEnsembleError):
                native.pairs()
            with self.assertRaises(EmptyEnsembleError):
                native.minimum_pair_energies()

    def test_pair_minima_equal_independent_forced_pair_calculations(self):
        for lengths in ((5,), (3, 3), (2, 2, 2), (2,) * 4, (2,) * 5):
            spec = make_spec(lengths, mixed=True)
            n = spec.n_nucleotides
            for base in (PairingConstraints.unconstrained(n), PairingConstraints.unconstrained(n).assign(0, n - 1)):
                with self.native(spec, constraints=base, thermal_energy_kcal_mol=.6163207755) as native:
                    # Exercise max/outside before any PF or standalone minimum.
                    minimum, matrix = native.minimum_pair_energies()
                    self.assertEqual(minimum, native.minimum()[0])
                    for i in range(n):
                        for j in range(i + 1, n):
                            if base.partners[i] == j:
                                wanted = minimum
                            else:
                                try:
                                    forced = base.assign(i, j)
                                except ValueError:
                                    wanted = math.inf
                                else:
                                    with self.native(spec, constraints=forced, thermal_energy_kcal_mol=.6163207755) as constrained:
                                        try:
                                            wanted = constrained.minimum()[0]
                                        except EmptyEnsembleError:
                                            wanted = math.inf
                            if wanted == math.inf:
                                self.assertEqual(matrix[i][j], math.inf)
                            else:
                                self.assertAlmostEqual(matrix[i][j], wanted, delta=3e-12)
                    self.assertTrue(math.isfinite(native.log_partition()))

    def test_minimum_first_requests_preserve_subsequent_partition_pairs_and_samples(self):
        operations = (
            ("minimum", lambda native: native.minimum()),
            ("minimum pair energies", lambda native: native.minimum_pair_energies()),
            ("minimum structures", lambda native: native.minimum_structures(atol_kcal=.3, count_limit=5)),
        )
        specs = (make_spec((7,)), make_spec((7,), material="dna"), make_spec((3, 3, 2), mixed=True))
        for fixture in ("unit", "asymmetric"):
            for spec in specs:
                with self.native(spec, fixture=fixture) as baseline:
                    log_partition = baseline.log_partition()
                    pairs = baseline.pairs()
                    samples = baseline.sample(32, seed=8291)
                    minima = {name: operation(baseline) for name, operation in operations}
                for name, operation in operations:
                    with self.subTest(fixture=fixture, spec=spec, first=name), \
                         self.native(spec, fixture=fixture) as native:
                        # Fresh minimum-only requests need no partition table.
                        # Its later lazy initialization must preserve every
                        # observable and the seeded traceback law on this handle.
                        self.assertEqual(operation(native), minima[name])
                        self.assertEqual(native.log_partition(), log_partition)
                        self.assertEqual(native.pairs(), pairs)
                        self.assertEqual(native.sample(32, seed=8291), samples)
                        self.assertEqual(operation(native), minima[name])
                        self.assertEqual(native.log_partition(), log_partition)
                        self.assertEqual(native.pairs(), pairs)

    def test_minimum_structures_match_exhaustive_global_windows_and_caps(self):
        model = AsymmetricWeights()
        for lengths in ((6,), (3, 3), (2, 2, 2), (2,) * 4, (2,) * 5):
            spec = make_spec(lengths, mixed=True)
            reference = ExhaustiveBackend(pair_allowed=model.pair_allowed, energy=model.energy,
                thermal_energy=lambda conditions: 1, provenance=PROVENANCE,
                min_hairpin_size=0).evaluate(spec, Conditions())
            expected = {format_structure(s, lengths): e for s, e in zip(reference.structures, reference.energies_kcal_mol)}
            global_factor = -model.global_energy(spec)
            with self.native(spec, log_global_factor=global_factor) as native:
                minimum, first = native.minimum()
                for tolerance in (0, .05, .3, 1.5, 100):
                    eligible = {s: e for s, e in expected.items() if e - minimum <= tolerance + 3e-12}
                    for cap in sorted({0, 1, 2, len(eligible), len(eligible) + 1}):
                        with self.subTest(lengths=lengths, tolerance=tolerance, cap=cap):
                            actual_minimum, values, truncated = native.minimum_structures(atol_kcal=tolerance, count_limit=cap)
                            self.assertEqual(actual_minimum, minimum)
                            self.assertEqual(len(values), min(cap, len(eligible)))
                            self.assertEqual(truncated, len(eligible) > cap)
                            self.assertEqual(len({s for _, s in values}), len(values))
                            if values:
                                self.assertEqual(values[0], (minimum, first))
                            for energy, structure in values:
                                self.assertIn(structure, eligible)
                                self.assertAlmostEqual(energy, expected[structure], delta=3e-12)
                                self.assertLessEqual(energy - minimum, tolerance + 3e-12)
                            if not truncated:
                                self.assertEqual({s for _, s in values}, set(eligible))

    def test_minimum_structures_degeneracy_constraints_and_global_offsets(self):
        spec = make_spec((2, 2, 2), mixed=True)
        constraints = PairingConstraints.unconstrained(6).assign(0, 5)
        for state in (None, constraints):
            all_states = tuple(enumerate_structures(spec, pair_allowed=lambda c, i, j: True,
                                                   min_hairpin_size=0, constraints=state))
            for factor in (0, -1e6, 1e6):
                with self.native(spec, fixture="unit", constraints=state, log_global_factor=factor) as native:
                    for cap in (0, 1, len(all_states), len(all_states) + 1):
                        minimum, values, truncated = native.minimum_structures(count_limit=cap)
                        self.assertEqual(minimum, -factor)
                        self.assertEqual(truncated, cap < len(all_states))
                        self.assertEqual(len(values), min(cap, len(all_states)))
                        self.assertTrue(all(e == -factor for e, _ in values))
                        if not truncated:
                            self.assertEqual({s for _, s in values}, {format_structure(s, (2, 2, 2)) for s in all_states})
                    first = native.minimum()
                    self.assertEqual(native.minimum_structures(count_limit=1)[1], (first,))
        with self.native(make_spec((1,)), fixture="unit") as native:
            self.assertEqual(native.minimum_structures(count_limit=1), (0.0, ((0.0, "."),), False))

    def test_minimum_structures_reject_invalid_tolerances_and_buffers(self):
        with self.native(make_spec((2, 2))) as native:
            for tolerance in (-1, math.nan, math.inf, -math.inf, True):
                with self.assertRaises(ValueError):
                    native.minimum_structures(atol_kcal=tolerance)
            for cap in (-1, 2**64 - 1, 2**64, True):
                with self.assertRaises(ValueError):
                    native.minimum_structures(count_limit=cap)
            output, energies = C.create_string_buffer(5), (C.c_double * 1)()
            count, truncated, minimum = C.c_uint64(), C.c_uint8(), C.c_double()
            for cells, size, cap in ((0, 5, 1), (1, 4, 1), (1, 5, 2**64 - 1)):
                with self.assertRaises(ValueError):
                    native._call(native.library.aptadp_minimum_structures, native.handle, 0, cap,
                        output, size, energies, cells, C.byref(count), C.byref(truncated), C.byref(minimum))
        with self.native(make_spec((2,)), thermal_energy_kcal_mol=.1) as native:
            with self.assertRaisesRegex(ValueError, "atol/RT"):
                native.minimum_structures(atol_kcal=1e308)

    def test_seeded_native_tracebacks_follow_the_complete_weighted_law(self):
        spec, model = make_spec((2, 2), mixed=True), AsymmetricWeights()
        reference = self.compare(spec)
        with self.native(spec) as native:
            samples = native.sample(12000, seed=1938433)
            self.assertEqual(samples, native.sample(12000, seed=1938433))
            counts = Counter(samples)
            self.assertEqual(set(counts), {format_structure(state, (2, 2)) for state in reference.structures})
            for state, p in zip(reference.structures, reference.probabilities):
                self.assertLess(abs(counts[format_structure(state, (2, 2))] - 12000 * p),
                                6 * math.sqrt(12000 * p * (1 - p)))

    def test_native_boundary_errors_and_memory_budget_are_explicit(self):
        spec = make_spec((2, 2), mixed=True)
        for options in ({"memory_limit_bytes": 1}, {"thermal_energy_kcal_mol": math.nan},
                        {"thermal_energy_kcal_mol": 0}, {"log_global_factor": math.inf},
                        {"minimum": 2**32}, {"max_interior_size": 31}, {"thermal_energy_kcal_mol": True}):
            with self.subTest(options=options), self.assertRaises(ValueError):
                self.native(spec, **options)
        with self.native(spec) as native:
            self.assertGreater(native.memory_estimate(), 0)
            out, adjustment = (C.c_double * 16)(), C.c_double()
            with self.assertRaisesRegex(ValueError, "buffer"):
                native._call(native.library.aptadp_pair_probabilities, native.handle, out, 15, C.byref(adjustment))
            with self.assertRaisesRegex(ValueError, "buffer"):
                native._call(native.library.aptadp_minimum_pair_energies, native.handle, out, 15, C.byref(adjustment))
            with self.assertRaisesRegex(ValueError, "NULL"):
                native._call(native.library.aptadp_minimum_pair_energies, native.handle, out, 16, None)
            native.close()
            with self.assertRaisesRegex(ValueError, "NULL"):
                native.log_partition()

    def test_c_abi_rejects_invalid_structs_materials_cuts_masks_and_constraints(self):
        with self.native(make_spec((2,))) as owner:
            lib = owner.library
            def fresh():
                request = NativeInput(C.sizeof(NativeInput), 2, b"GC", (C.c_uint8 * 2)(2, 2),
                    (C.c_uint8 * 3)(0, 0, 0), (C.c_uint8 * 4)(0, 1, 1, 0),
                    (C.c_int32 * 2)(-1, -1), 0, 30, 0)
                return request, NativeModel(C.sizeof(NativeModel), 1, 1.0, 0.0)
            mutations = (
                lambda r, m: setattr(r, "struct_size", 0),
                lambda r, m: setattr(m, "struct_size", 0),
                lambda r, m: setattr(m, "fixture_provider", 99),
                lambda r, m: setattr(r, "sequence", None),
                lambda r, m: setattr(r, "materials", (C.c_uint8 * 2)(3, 2)),
                lambda r, m: setattr(r, "materials", (C.c_uint8 * 2)(1, 2)),
                lambda r, m: setattr(r, "nicks_before", (C.c_uint8 * 3)(1, 0, 0)),
                lambda r, m: setattr(r, "nicks_before", (C.c_uint8 * 3)(0, 2, 0)),
                lambda r, m: setattr(r, "pair_allowed", (C.c_uint8 * 4)(0, 1, 0, 0)),
                lambda r, m: setattr(r, "pair_allowed", (C.c_uint8 * 4)(1, 1, 1, 0)),
                lambda r, m: setattr(r, "constraints", (C.c_int32 * 2)(1, -1)),
                lambda r, m: setattr(r, "constraints", (C.c_int32 * 2)(2, -1)),
            )
            for mutate in mutations:
                request, model = fresh()
                mutate(request, model)
                handle, error = C.c_void_p(), C.create_string_buffer(300)
                status = lib.aptadp_create_v1(C.byref(request), C.byref(model), C.byref(handle), error, len(error))
                self.assertEqual(status, 1, error.value)
                self.assertFalse(handle.value)
                self.assertTrue(error.value)


if __name__ == "__main__":
    unittest.main()
