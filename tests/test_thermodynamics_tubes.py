"""Tube conversion identities; all partitions and energies here are synthetic."""

from dataclasses import FrozenInstanceError, replace
import math
import unittest

import numpy as np
from numpy.testing import assert_allclose, assert_array_equal

from aptaswitch_core.thermodynamics.equilibrium import EquilibriumConvergenceError
from aptaswitch_core.thermodynamics.structures import PairingConstraints
from aptaswitch_core.thermodynamics.tubes import StandardStateConvention, solve_tube
from aptaswitch_core.thermodynamics.types import (
    ComplexSpec,
    Conditions,
    ModelProvenance,
    PartitionResult,
    Strand,
    enumerate_complexes,
)


PROVENANCE = ModelProvenance(
    backend="synthetic", version="1", model="analytic-fixture",
    parameter_digest="sha256:" + "a" * 64,
    ensemble="synthetic-connected-ordered", standard_state="audited-fixture-v1",
)
CONDITIONS = Conditions(celsius=31, sodium=0.3)


def convention(concentration=1.0, **changes):
    values = dict(
        provenance_label=PROVENANCE.standard_state,
        concentration_molar=concentration, association_included=True,
        cyclic_symmetry_included=True, connected_ensembles=True,
    )
    values.update(changes)
    return StandardStateConvention(**values)


def partition(strands, log_z, *, provenance=PROVENANCE, conditions=CONDITIONS, rt=0.6):
    return PartitionResult(
        complex_spec=ComplexSpec(tuple(strands)), conditions=conditions,
        provenance=provenance, log_partition=log_z,
        free_energy_kcal_mol=math.inf if log_z == -math.inf else -rt * log_z,
        thermal_energy_kcal_mol=rt,
    )


class StandardStateTests(unittest.TestCase):
    def test_concentration_is_explicit_and_flags_are_required(self):
        self.assertEqual(convention().concentration_molar, 1)
        self.assertEqual(convention(55.1).concentration_molar, 55.1)
        for value in (0, -1, math.inf, math.nan, True, "1 M"):
            with self.subTest(value=value), self.assertRaises(ValueError):
                convention(value)
        for flag in ("association_included", "cyclic_symmetry_included", "connected_ensembles"):
            for value in (False, 1, "yes", None):
                with self.subTest(flag=flag, value=value), self.assertRaises(ValueError):
                    convention(**{flag: value})
        for label in ("", " ", None, 1):
            with self.assertRaises(ValueError):
                convention(provenance_label=label)
        with self.assertRaises(FrozenInstanceError):
            convention().concentration_molar = 55.1


class TubeTests(unittest.TestCase):
    def setUp(self):
        # Equal alphabet strings but different chemistry must stay two species.
        self.a = Strand("ACG", "rna", "switch")
        self.b = Strand("ACG", "dna", "trigger")
        self.c = Strand("UU", "rna", "aptamer")
        self.monomers = (partition((self.a,), 3.0), partition((self.b,), 5.0))

    def heterodimer(self, association=8e6, standard=1):
        return partition((self.a, self.b), 8 + math.log(association) + math.log(standard))

    def test_heterodimer_uses_formation_partition_ratio(self):
        association = 8e6
        ca, cb = 2e-6, 7e-6
        records = self.monomers + (self.heterodimer(association),)
        result = solve_tube((self.a, self.b), (ca, cb), records, standard_state=convention())
        scale = ca + cb + 1 / association
        bound = 2 * ca * cb / (scale + math.sqrt(scale**2 - 4 * ca * cb))
        assert_allclose(result.concentrations_molar, [ca - bound, cb - bound, bound], rtol=1e-10)
        assert_allclose(result.log_equilibrium_constants, [0, 0, math.log(association)], atol=1e-14)
        assert_array_equal(result.stoichiometry, [[1, 0, 1], [0, 1, 1]])
        self.assertEqual(result.conditions, CONDITIONS)
        self.assertEqual(result.provenance, PROVENANCE)
        self.assertEqual(result.thermal_energy_kcal_mol, 0.6)

    def test_homodimer_symmetry_and_association_are_not_counted_twice(self):
        association, total = 3e7, 2e-6
        monomer = self.monomers[0]
        dimer = partition((self.a, self.a), 2 * monomer.log_partition + math.log(association))
        self.assertEqual(dimer.complex_spec.symmetry_factor, 2)
        result = solve_tube((self.a,), [total], (dimer, monomer), standard_state=convention())
        free = 2 * total / (1 + math.sqrt(1 + 8 * association * total))
        assert_allclose(result.concentrations_molar, [(total - free) / 2, free], rtol=1e-10)
        assert_allclose(result.log_equilibrium_constants, [math.log(association), 0], atol=1e-14)

    def test_water_and_molar_partitions_convert_to_same_physical_tube(self):
        results = []
        for standard in (1, 55.1):
            records = self.monomers + (self.heterodimer(standard=standard),)
            results.append(solve_tube((self.a, self.b), [2e-6, 7e-6], records, standard_state=convention(standard)))
        assert_allclose(results[0].concentrations_molar, results[1].concentrations_molar, rtol=1e-12)

    def test_no_concentration_is_guessed_from_a_free_text_label(self):
        source = replace(PROVENANCE, standard_state="1 M; association and symmetry included")
        records = tuple(replace(record, provenance=source) for record in self.monomers + (self.heterodimer(),))
        result = solve_tube(
            (self.a, self.b), [2e-6, 7e-6], records,
            standard_state=convention(55.1, provenance_label=source.standard_state),
        )
        # The explicitly supplied numerical convention governs the conversion;
        # a misleading label cannot silently override the caller's declaration.
        self.assertAlmostEqual(result.log_equilibrium_constants[-1], math.log(8e6 / 55.1), places=12)

    def test_offsets_in_monomer_reference_energies_cancel(self):
        records = self.monomers + (self.heterodimer(),)
        reference = solve_tube((self.a, self.b), [2e-6, 7e-6], records, standard_state=convention())
        offsets = {self.a.chemical_identity: 7.5, self.b.chemical_identity: -3.5}
        shifted = tuple(partition(record.complex_spec.strands,
                                  record.log_partition + sum(offsets[s.chemical_identity] for s in record.complex_spec.strands))
                        for record in records)
        result = solve_tube((self.a, self.b), [2e-6, 7e-6], shifted, standard_state=convention())
        assert_allclose(result.concentrations_molar, reference.concentrations_molar, rtol=1e-12)

    def test_log_partition_and_native_rt_are_used_without_exp_or_new_gas_constant(self):
        records = (partition((self.a,), 5000, rt=0.7123),
                   partition((self.b,), 4200, rt=0.7123),
                   partition((self.a, self.b), 9200 + math.log(8e6), rt=0.7123))
        result = solve_tube((self.a, self.b), [2e-6, 7e-6], records, standard_state=convention())
        self.assertAlmostEqual(result.log_equilibrium_constants[2], math.log(8e6), places=11)
        self.assertEqual(result.thermal_energy_kcal_mol, 0.7123)
        self.assertTrue(np.all(np.isfinite(result.concentrations_molar)))

    def test_alias_names_map_by_chemistry_and_orders_are_unchanged(self):
        alias = replace(self.a, name="RNA alias")
        dimer = partition((self.b, alias), 8 + math.log(8e6))
        records = [dimer, self.monomers[1], partition((alias,), 3)]
        totals = np.array([7e-6, 2e-6])
        result = solve_tube((self.b, self.a), totals, records, standard_state=convention())
        assert_array_equal(result.stoichiometry, [[1, 1, 0], [1, 0, 1]])
        self.assertEqual(result.species, (self.b, self.a))
        self.assertIs(result.partitions[0], dimer)
        self.assertEqual(result.partitions[0].complex_spec.strands, (self.b, alias))
        self.assertEqual(result.partitions[0].complex_spec.cutpoints, (3,))
        records.clear()
        totals[:] = 0
        self.assertEqual(len(result.partitions), 3)
        assert_array_equal(result.total_concentrations_molar, [7e-6, 2e-6])
        for array in (result.total_concentrations_molar, result.stoichiometry,
                      result.log_equilibrium_constants, result.concentrations_molar,
                      result.log_concentrations_molar):
            self.assertFalse(array.flags.writeable)

    def test_distinct_cyclic_orders_remain_distinct_species(self):
        species = (self.a, self.b, self.c)
        monomers = self.monomers + (partition((self.c,), 4),)
        abc = partition(species, 12 + math.log(4e12))
        acb = partition((self.a, self.c, self.b), 12 + math.log(1e12))
        result = solve_tube(species, [2e-6, 3e-6, 4e-6], monomers + (abc, acb), standard_state=convention())
        assert_array_equal(result.stoichiometry[:, -2:], [[1, 1], [1, 1], [1, 1]])
        self.assertAlmostEqual(result.concentrations_molar[-2] / result.concentrations_molar[-1], 4, places=11)
        self.assertIs(result.partitions[-2], abc)
        self.assertIs(result.partitions[-1], acb)

    def test_rotations_and_chemical_alias_duplicates_are_rejected(self):
        abc = partition((self.a, self.b, self.c), 20)
        bca = partition((self.b, self.c, self.a), 20)
        alias = partition((replace(self.a, name="alias"),), 3)
        for species, records in (
            ((self.a, self.b), self.monomers + (self.heterodimer(), self.heterodimer())),
            ((self.a, self.b), self.monomers + (alias,)),
            ((self.a, self.b, self.c), self.monomers + (partition((self.c,), 4), abc, bca)),
        ):
            with self.subTest(records=records), self.assertRaisesRegex(ValueError, "Duplicate"):
                solve_tube(species, [1e-6] * len(species), records, standard_state=convention())

    def test_empty_multimer_and_absent_strand_keep_aligned_zero_entries(self):
        impossible = partition((self.a, self.b), -math.inf)
        result = solve_tube((self.a, self.b), [2e-6, 7e-6], self.monomers + (impossible,), standard_state=convention())
        assert_allclose(result.concentrations_molar, [2e-6, 7e-6, 0], rtol=1e-14)
        self.assertEqual(result.log_equilibrium_constants[-1], -math.inf)
        self.assertEqual(result.log_concentrations_molar[-1], -math.inf)
        records = self.monomers + (self.heterodimer(),)
        for totals in ([2e-6, 0], [0, 0]):
            with self.subTest(totals=totals):
                result = solve_tube((self.a, self.b), totals, records, standard_state=convention())
                assert_allclose(result.concentrations_molar, [totals[0], 0, 0], rtol=1e-14)
                self.assertIs(result.partitions[-1], records[-1])

    def test_all_monomers_are_required_and_cannot_be_empty(self):
        for records in (self.monomers[:1] + (self.heterodimer(),),
                        (partition((self.a,), -math.inf), self.monomers[1], self.heterodimer())):
            with self.subTest(records=records), self.assertRaisesRegex(ValueError, "monomer"):
                solve_tube((self.a, self.b), [0, 0], records, standard_state=convention())

    def test_every_provenance_field_must_match(self):
        changes = dict(backend="another", version="2", model="other-model",
                       parameter_digest="sha256:" + "b" * 64,
                       ensemble="different-ensemble", standard_state="other-standard")
        for field, value in changes.items():
            record = replace(self.heterodimer(), provenance=replace(PROVENANCE, **{field: value}))
            with self.subTest(field=field), self.assertRaisesRegex(ValueError, "provenance"):
                solve_tube((self.a, self.b), [2e-6, 7e-6], self.monomers + (record,), standard_state=convention())

    def test_conditions_and_rt_mismatch_are_rejected(self):
        for field, value in (("celsius", 32), ("sodium", 0.31), ("magnesium", 0.001)):
            changed = replace(self.monomers[1], conditions=replace(CONDITIONS, **{field: value}))
            with self.subTest(field=field), self.assertRaisesRegex(ValueError, "conditions"):
                solve_tube((self.a, self.b), [2e-6, 7e-6], (self.monomers[0], changed, self.heterodimer()), standard_state=convention())
        changed = partition((self.b,), 5, rt=0.600000000000001)
        with self.assertRaisesRegex(ValueError, "RT"):
            solve_tube((self.a, self.b), [2e-6, 7e-6], (self.monomers[0], changed, self.heterodimer()), standard_state=convention())

    def test_conditional_partitions_cannot_enter_an_unconstrained_tube(self):
        constrained = replace(self.monomers[0], constraints=PairingConstraints((0, -1, -1)))
        with self.assertRaisesRegex(ValueError, "Constrained"):
            solve_tube((self.a, self.b), [2e-6, 7e-6], (constrained, self.monomers[1], self.heterodimer()), standard_state=convention())
        constrained_dimer = replace(self.heterodimer(), constraints=PairingConstraints((5, -1, -1, -1, -1, 0)))
        with self.assertRaisesRegex(ValueError, "Constrained"):
            solve_tube((self.a, self.b), [2e-6, 7e-6], self.monomers + (constrained_dimer,), standard_state=convention())
        unconstrained = replace(self.monomers[0], constraints=PairingConstraints.unconstrained(3))
        self.assertIsNone(unconstrained.constraints)
        result = solve_tube((self.a, self.b), [2e-6, 7e-6], (unconstrained, self.monomers[1], self.heterodimer()), standard_state=convention())
        self.assertLess(result.equilibrium.max_relative_mass_balance_error, 1e-10)

    def test_convention_mismatch_and_missing_structured_convention_are_rejected(self):
        for state in ("1 M", None, convention(provenance_label="unreviewed label")):
            with self.subTest(state=state), self.assertRaises(ValueError):
                solve_tube((self.a, self.b), [2e-6, 7e-6], self.monomers, standard_state=state)
        with self.assertRaises(TypeError):
            solve_tube((self.a, self.b), [2e-6, 7e-6], self.monomers)

    def test_catalog_must_be_unique_and_contain_every_complex_strand(self):
        for species in ((self.a, self.a), (self.a, replace(self.a, name="alias")),
                        (self.a, replace(self.b, name=self.a.name)), (), None, ("ACG",)):
            with self.subTest(species=species), self.assertRaises(ValueError):
                solve_tube(species, [2e-6, 7e-6], self.monomers, standard_state=convention())
        with self.assertRaisesRegex(ValueError, "absent from"):
            solve_tube((self.a, self.b), [2e-6, 7e-6], self.monomers + (partition((self.c,), 4),), standard_state=convention())
        for records in ((), None, ("partition",)):
            with self.assertRaises(ValueError):
                solve_tube((self.a, self.b), [2e-6, 7e-6], records, standard_state=convention())

    def test_invalid_totals_and_solver_failure_are_not_hidden(self):
        records = self.monomers + (self.heterodimer(),)
        for totals in ([1e-6], [1e-6, -1e-6], [math.nan, 1e-6], [math.inf, 1e-6]):
            with self.subTest(totals=totals), self.assertRaises(ValueError):
                solve_tube((self.a, self.b), totals, records, standard_state=convention())
        with self.assertRaises(EquilibriumConvergenceError):
            solve_tube((self.a, self.b), [2e-6, 7e-6], records, standard_state=convention(), max_iterations=1)

    def test_totals_cannot_be_coerced_from_bool_string_or_complex(self):
        for totals in ([True, 7e-6], [2e-6, False], ["2e-6", 7e-6],
                       [2e-6, 7e-6 + 0j], np.array([2e-6, 7e-6 + 1j]),
                       np.array([True, False]), np.array(["2e-6", "7e-6"])):
            with self.subTest(totals=totals), self.assertRaises(ValueError):
                solve_tube((self.a, self.b), totals, self.monomers, standard_state=convention())

    def test_full_23_and_95_complex_catalogs_recover_known_tubes(self):
        for species in ((self.a, self.b), (self.a, self.b, self.c)):
            for standard in (1, 55.1):
                with self.subTest(strands=len(species), standard=standard):
                    free = {strand.chemical_identity: (index + 1) * 1e-7 for index, strand in enumerate(species)}
                    mono_log = {strand.chemical_identity: 3 + index for index, strand in enumerate(species)}
                    records, expected = [], []
                    totals = np.zeros(len(species))
                    indices = {strand.chemical_identity: index for index, strand in enumerate(species)}
                    for index, spec in enumerate(enumerate_complexes(species, max_size=5)):
                        size = len(spec.strands)
                        log_k = 0 if size == 1 else (size - 1) * math.log(1e6) - index * 0.03
                        log_z = log_k + sum(mono_log[strand.chemical_identity] for strand in spec.strands) + (size - 1) * math.log(standard)
                        records.append(partition(spec.strands, log_z))
                        concentration = math.exp(log_k) * math.prod(free[strand.chemical_identity] for strand in spec.strands)
                        expected.append(concentration)
                        for strand in spec.strands:
                            totals[indices[strand.chemical_identity]] += concentration
                    result = solve_tube(species, totals, records, standard_state=convention(standard))
                    self.assertEqual(len(records), 23 if len(species) == 2 else 95)
                    assert_allclose(result.concentrations_molar, expected, rtol=1e-10, atol=0)
                    assert_allclose(result.stoichiometry @ result.concentrations_molar, totals, rtol=1e-10, atol=0)


if __name__ == "__main__":
    unittest.main()
