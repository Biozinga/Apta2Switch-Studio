"""Synthetic thermodynamic identities, independent of any folding engine."""

import itertools
import math
import unittest
from unittest.mock import patch

import numpy as np
from numpy.testing import assert_allclose, assert_array_equal

from aptaswitch_core.thermodynamics.equilibrium import (
    GAS_CONSTANT_KCAL_MOL_K,
    EquilibriumConvergenceError,
    formation_energies_to_log_constants,
    solve_equilibrium,
)


def cyclic_compositions(strand_count, max_size=5):
    words = [
        word
        for size in range(1, max_size + 1)
        for word in itertools.product(range(strand_count), repeat=size)
        if word == min(word[offset:] + word[:offset] for offset in range(size))
    ]
    return np.array([[word.count(index) for word in words] for index in range(strand_count)])


class EquilibriumTests(unittest.TestCase):
    def assert_valid(self, result, matrix, totals):
        assert_allclose(matrix @ result.concentrations_molar, totals, rtol=1e-10, atol=0)
        self.assertLessEqual(result.max_relative_mass_balance_error, 1e-10)
        self.assertLessEqual(result.max_log_mass_action_error, 1e-10)
        self.assertFalse(result.concentrations_molar.flags.writeable)
        self.assertFalse(result.log_concentrations_molar.flags.writeable)
        self.assertFalse(result.mass_balance_residual_molar.flags.writeable)

    def test_analytical_homodimer_includes_two_strands(self):
        total, association = 3e-6, 2e7
        free = 2 * total / (1 + math.sqrt(1 + 8 * association * total))
        matrix = np.array([[1, 2]])
        result = solve_equilibrium(matrix, [total], log_equilibrium_constants=[0, math.log(association)])
        assert_allclose(result.concentrations_molar, [free, (total - free) / 2], rtol=1e-10)
        self.assert_valid(result, matrix, [total])

    def test_analytical_heterodimer_unequal_totals(self):
        ca, cb, association = 2e-6, 7e-6, 8e6
        scale = ca + cb + 1 / association
        bound = 2 * ca * cb / (scale + math.sqrt(scale * scale - 4 * ca * cb))
        matrix = np.array([[1, 0, 1], [0, 1, 1]])
        result = solve_equilibrium(matrix, [ca, cb], log_equilibrium_constants=[0, 0, math.log(association)])
        assert_allclose(result.concentrations_molar, [ca - bound, cb - bound, bound], rtol=1e-10)
        self.assert_valid(result, matrix, [ca, cb])

    def test_analytical_trimer_reconstructed_from_known_free_species(self):
        # Choose free monomers first; mass action independently supplies a tube
        # with known equilibrium, including strand multiplicity in AAB.
        free_a, free_b, association = 7e-8, 3e-7, 5e14
        trimer = association * free_a**2 * free_b
        expected = np.array([free_a, free_b, trimer])
        matrix = np.array([[1, 0, 2], [0, 1, 1]])
        totals = matrix @ expected
        result = solve_equilibrium(matrix, totals, log_equilibrium_constants=[0, 0, math.log(association)])
        assert_allclose(result.concentrations_molar, expected, rtol=1e-10)
        self.assert_valid(result, matrix, totals)

    def test_23_and_95_cyclic_complexes_with_distinct_orders(self):
        for strand_count, expected_count in ((2, 23), (3, 95)):
            with self.subTest(strand_count=strand_count):
                matrix = cyclic_compositions(strand_count)
                self.assertEqual(matrix.shape[1], expected_count)
                sizes = matrix.sum(axis=0)
                constants = (sizes - 1) * math.log(1e6) - 0.03 * np.arange(1, expected_count + 1)
                constants[:strand_count] = 0
                totals = np.array([2e-6, 3e-6, 1e-6])[:strand_count]
                result = solve_equilibrium(matrix, totals, log_equilibrium_constants=constants)
                self.assert_valid(result, matrix, totals)
                log_c = result.log_concentrations_molar
                assert_allclose(log_c - matrix.T @ log_c[:strand_count], constants, atol=1e-12)
                # Deliberately different energies for cyclic orders with the
                # same counts must survive; the solver must not merge them.
                if strand_count == 3:
                    columns = np.flatnonzero(np.all(matrix == np.array([1, 1, 1])[:, None], axis=0))
                    self.assertEqual(len(columns), 2)
                    first, second = columns
                    self.assertAlmostEqual(log_c[first] - log_c[second], constants[first] - constants[second], places=12)

    def test_extreme_affinities_do_not_exponentiate_constants(self):
        for constant in (-1000.0, 1000.0):
            with self.subTest(log_k=constant):
                matrix = np.array([[1, 2]])
                total = 2e-6
                result = solve_equilibrium(matrix, [total], log_equilibrium_constants=[0, constant])
                self.assert_valid(result, matrix, [total])
                if constant > 0:
                    self.assertAlmostEqual(result.log_concentrations_molar[0], (math.log(total / 2) - constant) / 2, places=10)
                    self.assertAlmostEqual(result.concentrations_molar[1] / (total / 2), 1, places=10)
                else:
                    self.assertAlmostEqual(result.concentrations_molar[0] / total, 1, places=10)
                    self.assertEqual(result.concentrations_molar[1], 0)  # float underflow only
                    self.assertTrue(np.isfinite(result.log_concentrations_molar[1]))

    def test_strong_heterodimer_handles_limiting_and_equivalent_strands(self):
        matrix = np.array([[1, 0, 1], [0, 1, 1]])
        for totals in ([2e-6, 2e-6], [2e-6, 7e-6]):
            with self.subTest(totals=totals):
                result = solve_equilibrium(matrix, totals, log_equilibrium_constants=[0, 0, 1000])
                self.assert_valid(result, matrix, totals)
                self.assertAlmostEqual(result.concentrations_molar[2] / min(totals), 1, places=10)
                if totals[0] == totals[1]:
                    assert_allclose(result.log_concentrations_molar[:2], (math.log(totals[0]) - 1000) / 2, atol=1e-10)
                else:
                    self.assertAlmostEqual(result.concentrations_molar[1] / (totals[1] - totals[0]), 1, places=10)
                    self.assertTrue(np.isfinite(result.log_concentrations_molar[0]))

    def test_strong_trimer_preserves_stoichiometric_free_strand_ratio(self):
        # Exact subtraction of total_A - 2*total_B must imply free_A=2*free_B,
        # even when neither free species contributes a representable amount
        # to the original mass sums.
        matrix = np.array([[1, 0, 2], [0, 1, 1]])
        totals = [4e-6, 2e-6]
        result = solve_equilibrium(matrix, totals, log_equilibrium_constants=[0, 0, 1000])
        expected_log_b = (math.log(totals[1]) - math.log(4) - 1000) / 3
        assert_allclose(result.log_concentrations_molar[:2], [expected_log_b + math.log(2), expected_log_b], atol=1e-10, rtol=0)
        self.assert_valid(result, matrix, totals)

    def test_almost_equal_totals_keep_the_small_exact_excess(self):
        matrix = np.array([[1, 0, 1], [0, 1, 1]])
        totals = [2e-6, 2e-6 * (1 + 1e-12)]
        result = solve_equilibrium(matrix, totals, log_equilibrium_constants=[0, 0, 1000])
        expected_log_b = math.log(totals[1] - totals[0])
        expected_log_a = math.log(totals[0]) - 1000 - expected_log_b
        assert_allclose(result.log_concentrations_molar[:2], [expected_log_a, expected_log_b], atol=1e-10, rtol=0)
        self.assert_valid(result, matrix, totals)

    def test_seeded_multicomplex_tubes_recover_known_equilibria(self):
        rng = np.random.default_rng(20261005)
        for strands in (2, 3):
            matrix = cyclic_compositions(strands)
            size = matrix.sum(axis=0)
            for trial in range(5):
                with self.subTest(strands=strands, trial=trial):
                    log_free = rng.uniform(math.log(1e-9), math.log(1e-7), strands)
                    constants = (size - 1) * math.log(1e6) + rng.uniform(-10, 10, matrix.shape[1])
                    constants[:strands] = 0
                    expected = np.exp(constants + matrix.T @ log_free)
                    totals = matrix @ expected
                    result = solve_equilibrium(matrix, totals, log_equilibrium_constants=constants)
                    assert_allclose(result.concentrations_molar, expected, rtol=1e-10, atol=0)
                    self.assert_valid(result, matrix, totals)

    def test_absent_strand_excludes_every_containing_complex(self):
        matrix = cyclic_compositions(3)
        constants = (matrix.sum(axis=0) - 1) * math.log(1e7)
        for totals in ([2e-6, 0, 3e-6], [0, 0, 0]):
            with self.subTest(totals=totals):
                result = solve_equilibrium(matrix, totals, log_equilibrium_constants=constants)
                excluded = np.any(matrix[np.array(totals) == 0] > 0, axis=0)
                assert_array_equal(result.concentrations_molar[excluded], 0)
                assert_array_equal(result.log_concentrations_molar[excluded], -np.inf)
                self.assert_valid(result, matrix, totals)

    def test_impossible_complex_has_exact_zero_without_losing_monomers(self):
        result = solve_equilibrium([[1, 0, 1], [0, 1, 1]], [1e-6, 2e-6], log_equilibrium_constants=[0, 0, -np.inf])
        assert_allclose(result.concentrations_molar, [1e-6, 2e-6, 0], rtol=1e-14)
        self.assertEqual(result.log_concentrations_molar[2], -np.inf)
        self.assertEqual(result.iterations, 0)

    def test_monomers_only_and_permuted_columns(self):
        matrix = np.array([[0, 1], [1, 0]])
        result = solve_equilibrium(matrix, [1e-200, 2e-9], log_equilibrium_constants=[0, 0])
        self.assert_valid(result, matrix, [1e-200, 2e-9])
        assert_allclose(result.concentrations_molar, [2e-9, 1e-200], rtol=1e-13)

    def test_molar_dimensional_covariance(self):
        matrix = cyclic_compositions(3)
        sizes = matrix.sum(axis=0)
        constants = (sizes - 1) * math.log(1e6)
        totals = np.array([2e-6, 3e-6, 1e-6])
        reference = solve_equilibrium(matrix, totals, log_equilibrium_constants=constants)
        # All calls still use M. These are scaled physical systems whose K
        # changes by s**(1-n); no hidden uM conversion is allowed in the API.
        for factor in (1e-12, 1e6):
            result = solve_equilibrium(matrix, totals * factor, log_equilibrium_constants=constants + (1 - sizes) * math.log(factor))
            assert_allclose(result.concentrations_molar / factor, reference.concentrations_molar, rtol=1e-10)
            self.assert_valid(result, matrix, totals * factor)

    def test_free_energy_conversion_temperature_and_water_standard(self):
        matrix = np.array([[1, 0, 1], [0, 1, 1]])
        formation_energy = -8.0
        for temperature in (283.15, 310.15, 343.15):
            for standard in (1.0, 55.1):
                with self.subTest(temperature=temperature, standard=standard):
                    expected = -formation_energy / (GAS_CONSTANT_KCAL_MOL_K * temperature) - math.log(standard)
                    constants = formation_energies_to_log_constants(matrix, [0, 0, formation_energy], temperature_kelvin=temperature, standard_concentration_molar=standard)
                    assert_allclose(constants, [0, 0, expected], atol=1e-14)
                    direct = solve_equilibrium(matrix, [1e-6, 1e-6], log_equilibrium_constants=constants)
                    energy = solve_equilibrium(matrix, [1e-6, 1e-6], formation_free_energies_kcal_mol=[0, 0, formation_energy], temperature_kelvin=temperature, standard_concentration_molar=standard)
                    assert_array_equal(direct.concentrations_molar, energy.concentrations_molar)
                    association = math.exp(expected)
                    free = 2e-6 / (1 + math.sqrt(1 + 4 * association * 1e-6))
                    assert_allclose(energy.concentrations_molar, [free, free, 1e-6 - free], rtol=1e-10)

    def test_equivalent_water_and_one_molar_energies_do_not_double_correct(self):
        matrix = cyclic_compositions(3)
        constants = (matrix.sum(axis=0) - 1) * math.log(2e6)
        totals = [2e-6, 3e-6, 1e-6]
        rt = GAS_CONSTANT_KCAL_MOL_K * 310.15
        reference = solve_equilibrium(matrix, totals, log_equilibrium_constants=constants)
        for standard in (1.0, 55.1):
            formation_energy = -rt * (constants + (matrix.sum(axis=0) - 1) * math.log(standard))
            result = solve_equilibrium(matrix, totals, formation_free_energies_kcal_mol=formation_energy, standard_concentration_molar=standard)
            assert_allclose(result.concentrations_molar, reference.concentrations_molar, rtol=1e-10)

    def test_per_species_validation_for_totals_spanning_many_orders(self):
        matrix = np.array([[1, 0, 1], [0, 1, 1]])
        totals = [1e-30, 1e-3]
        result = solve_equilibrium(matrix, totals, log_equilibrium_constants=[0, 0, 100])
        self.assert_valid(result, matrix, totals)

    def test_iteration_exhaustion_refuses_partial_result(self):
        with self.assertRaises(EquilibriumConvergenceError):
            solve_equilibrium([[1, 2]], [2e-6], log_equilibrium_constants=[0, 1000], max_iterations=1)

    def test_optimizer_success_alone_does_not_validate_conservation(self):
        with patch("aptaswitch_core.thermodynamics.equilibrium.least_squares") as optimizer:
            optimizer.return_value.success = True
            optimizer.return_value.x = np.array([math.log(2e-6)])
            optimizer.return_value.nfev = 1
            with self.assertRaisesRegex(EquilibriumConvergenceError, "conservation"):
                solve_equilibrium([[1, 2]], [2e-6], log_equilibrium_constants=[0, 20])

    def test_invalid_input_is_rejected(self):
        for matrix, totals, constants in (
            ([[1, 2]], [-1], [0, 1]), ([[1, 2]], [np.nan], [0, 1]),
            ([[1, 2]], [np.inf], [0, 1]), ([[1, 2]], [1], [0, np.nan]),
            ([[1, 2]], [1], [0, np.inf]), ([[1, 2]], [1], [1, 2]),
            ([[1, 1, 2]], [1], [0, 0, 1]), ([[2]], [1], [1]),
            ([[1, -2]], [1], [0, 1]), ([[1, 1.5]], [1], [0, 1]),
            ([[1, 0]], [1], [0, 1]), ([[1, np.inf]], [1], [0, 1]),
            ([[1, 2]], [1, 2], [0, 1]), ([[1, 2]], [1], [0]),
        ):
            with self.subTest(matrix=matrix, totals=totals, constants=constants), self.assertRaises(ValueError):
                solve_equilibrium(matrix, totals, log_equilibrium_constants=constants)
        for kwargs in ({}, {"log_equilibrium_constants": [0, 1], "formation_free_energies_kcal_mol": [0, -1]}):
            with self.assertRaises(ValueError):
                solve_equilibrium([[1, 2]], [1], **kwargs)
        for kwargs in ({"relative_tolerance": 0}, {"relative_tolerance": np.nan}, {"max_iterations": 0}, {"max_iterations": 1.5}, {"max_iterations": True}, {"standard_concentration_molar": 55.1}):
            with self.subTest(kwargs=kwargs), self.assertRaises(ValueError):
                solve_equilibrium([[1, 2]], [1], log_equilibrium_constants=[0, 1], **kwargs)
        for kwargs in ({"temperature_kelvin": 0}, {"temperature_kelvin": np.nan}, {"standard_concentration_molar": -1}, {"standard_concentration_molar": np.inf}):
            with self.subTest(kwargs=kwargs), self.assertRaises(ValueError):
                formation_energies_to_log_constants([[1, 2]], [0, -1], **kwargs)

    def test_scalar_types_are_validated_before_float_conversion(self):
        invalid = (True, np.bool_(False), "1", 1 + 0j, 1 + 2j)
        for value in invalid:
            with self.subTest(value=value):
                with self.assertRaises(ValueError):
                    solve_equilibrium([[1, 2]], [value], log_equilibrium_constants=[0, 1])
                with self.assertRaises(ValueError):
                    solve_equilibrium([[1, value]], [1], log_equilibrium_constants=[0, 1])
                with self.assertRaises(ValueError):
                    solve_equilibrium([[1, 2]], [1], log_equilibrium_constants=[0, value])
                with self.assertRaises(ValueError):
                    formation_energies_to_log_constants([[1, 2]], [0, value])
                for field in ("temperature_kelvin", "standard_concentration_molar"):
                    with self.assertRaises(ValueError):
                        formation_energies_to_log_constants([[1, 2]], [0, -1], **{field: value})
                with self.assertRaises(ValueError):
                    solve_equilibrium([[1, 2]], [1], log_equilibrium_constants=[0, 1], relative_tolerance=value)
        for array in (np.array([True]), np.array(["1"]), np.array([1 + 0j])):
            with self.subTest(dtype=array.dtype), self.assertRaises(ValueError):
                solve_equilibrium([[1, 2]], array, log_equilibrium_constants=[0, 1])
        for kwargs in ({"standard_concentration_molar": True},
                       {"standard_concentration_molar": 1 + 0j},
                       {"temperature_kelvin": 310.15 + 0j}):
            with self.assertRaises(ValueError):
                solve_equilibrium([[1, 2]], [1], log_equilibrium_constants=[0, 1], **kwargs)
        result = solve_equilibrium(np.array([[1]], dtype=np.int64),
                                   np.array([2e-6], dtype=np.float64),
                                   log_equilibrium_constants=np.array([0], dtype=np.int32))
        assert_allclose(result.concentrations_molar, [2e-6], rtol=1e-14)


if __name__ == "__main__":
    unittest.main()
