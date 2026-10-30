"""Chemical identity, cyclic complex counting, and thermodynamic contracts."""

from collections import Counter
from dataclasses import FrozenInstanceError, replace
import math
import unittest

from aptaswitch_core.thermodynamics.structures import PairingConstraints

from aptaswitch_core.thermodynamics import (
    ComplexSpec,
    Conditions,
    EmptyEnsembleError,
    ModelProvenance,
    PairProbabilities,
    PartitionResult,
    Strand,
    enumerate_complexes,
)


def provenance(**overrides):
    fields = dict(
        backend="reference", version="1", model="synthetic-test",
        parameter_digest="sha256:" + "a" * 64,
        ensemble="connected;pseudoknot-free;dangles=0;ordered",
        standard_state="1 M;association included;rotational symmetry included",
    )
    fields.update(overrides)
    return ModelProvenance(**fields)


def partition(**overrides):
    fields = dict(
        complex_spec=ComplexSpec((Strand("GC", "rna", "r"),)),
        conditions=Conditions(), provenance=provenance(),
        log_partition=2.0, free_energy_kcal_mol=-1.2,
        thermal_energy_kcal_mol=0.6,
    )
    fields.update(overrides)
    return PartitionResult(**fields)


class StrandTests(unittest.TestCase):
    def test_normalizes_case_and_surrounding_whitespace_only(self):
        strand = Strand("  acgt\n", "dna", " probe ")
        self.assertEqual((strand.sequence, strand.material, strand.name), ("ACGT", "dna", "probe"))
        self.assertEqual(Strand("acgu", "rna", "r").sequence, "ACGU")
        with self.assertRaises(FrozenInstanceError):
            strand.sequence = "AAAA"

    def test_chemistry_is_explicit_and_t_u_are_never_converted(self):
        for sequence, material in (("ACGU", "dna"), ("ACGT", "rna")):
            with self.subTest(material=material), self.assertRaises(ValueError):
                Strand(sequence, material, "x")
        for material in ("DNA", "RNA", "hybrid", "", None, True):
            with self.subTest(material=material), self.assertRaises(ValueError):
                Strand("ACG", material, "x")
        self.assertNotEqual(Strand("ACG", "dna", "x").chemical_identity,
                            Strand("ACG", "rna", "x").chemical_identity)

    def test_rejects_noncanonical_or_empty_input(self):
        for sequence in ("", " \n", "A C", "A\nC", "N", "AC-G", "AC&G", "AC+G", "A1", None, 123):
            with self.subTest(sequence=sequence), self.assertRaises(ValueError):
                Strand(sequence, "dna", "x")
        for name in ("", "  ", None, 1):
            with self.subTest(name=name), self.assertRaises(ValueError):
                Strand("ACG", "rna", name)

    def test_identity_does_not_depend_on_display_name(self):
        first = Strand("GCG", "rna", "first")
        alias = Strand("GCG", "rna", "alias")
        self.assertEqual(first.chemical_identity, alias.chemical_identity)
        self.assertNotEqual(first, alias)


class ComplexTests(unittest.TestCase):
    def setUp(self):
        self.a = Strand("ACG", "rna", "A")
        self.b = Strand("ACG", "dna", "B")
        self.c = Strand("UU", "rna", "C")

    def test_retains_caller_order_and_coordinates(self):
        original = [self.c, self.b, self.a]
        spec = ComplexSpec(original)
        original.reverse()
        self.assertEqual(spec.strands, (self.c, self.b, self.a))
        self.assertEqual(spec.sequence, "UUACGACG")
        self.assertEqual(spec.cutpoints, (2, 5))
        self.assertEqual(spec.n_nucleotides, 8)
        self.assertTrue(spec.is_mixed)
        with self.assertRaises(FrozenInstanceError):
            spec.strands = ()

    def test_rotations_identify_complex_without_reordering_or_merging_mirrors(self):
        abc = ComplexSpec((self.a, self.b, self.c))
        bca = ComplexSpec((self.b, self.c, self.a))
        acb = ComplexSpec((self.a, self.c, self.b))
        self.assertEqual(abc.cyclic_key, bca.cyclic_key)
        self.assertNotEqual(abc, bca)
        self.assertNotEqual(abc.ordered_key, bca.ordered_key)
        self.assertNotEqual(abc.cyclic_key, acb.cyclic_key)
        self.assertEqual(abc.strands, (self.a, self.b, self.c))

    def test_rotational_symmetry_uses_chemistry_and_sequence_not_names(self):
        alias = Strand(self.a.sequence, self.a.material, "another-name")
        fixtures = [((self.a,), 1), ((self.a,) * 5, 5),
                    ((self.a, alias), 2), ((self.a, self.b), 1),
                    ((self.a, self.b, self.a, self.b), 2),
                    ((self.a, self.a, self.b, self.b), 1)]
        for strands, expected in fixtures:
            with self.subTest(strands=strands):
                self.assertEqual(ComplexSpec(strands).symmetry_factor, expected)

    def test_empty_invalid_and_conflicting_names_are_rejected(self):
        for strands in ((), None, "ACG", ("ACG",),
                        (self.a, Strand("AAA", "rna", "A")),
                        (self.a, Strand("ACG", "dna", "A"))):
            with self.subTest(strands=strands), self.assertRaises(ValueError):
                ComplexSpec(strands)
        self.assertEqual(ComplexSpec((self.a, self.a)).cutpoints, (3,))

    def test_exact_necklace_counts_and_unique_chemical_keys(self):
        for species, expected in (
            ((self.a, self.b), {1: 2, 2: 3, 3: 4, 4: 6, 5: 8}),
            ((self.a, self.b, self.c), {1: 3, 2: 6, 3: 11, 4: 24, 5: 51}),
        ):
            with self.subTest(n_species=len(species)):
                complexes = enumerate_complexes(species)
                self.assertEqual(Counter(len(c.strands) for c in complexes), expected)
                self.assertEqual(len(complexes), 23 if len(species) == 2 else 95)
                self.assertEqual(len({c.cyclic_key for c in complexes}), len(complexes))
                self.assertTrue(all(c.ordered_key == c.cyclic_key for c in complexes))
                self.assertEqual(enumerate_complexes(reversed(species)), complexes)

    def test_necklaces_keep_different_orders_with_same_stoichiometry(self):
        complexes = enumerate_complexes((self.a, self.b), max_size=4)
        aabb = [c for c in complexes if Counter(c.strands) == Counter((self.a, self.a, self.b, self.b))]
        self.assertEqual(len(aabb), 2)
        self.assertEqual({c.symmetry_factor for c in aabb}, {1, 2})
        triples = enumerate_complexes((self.a, self.b, self.c), max_size=3)
        keys = {c.cyclic_key for c in triples}
        self.assertIn(ComplexSpec((self.a, self.b, self.c)).cyclic_key, keys)
        self.assertIn(ComplexSpec((self.a, self.c, self.b)).cyclic_key, keys)

    def test_catalog_aliases_and_name_collisions_cannot_duplicate_concentrations(self):
        for species in ((self.a, self.a),
                        (self.a, Strand("ACG", "rna", "alias")),
                        (self.a, Strand("AAA", "rna", "A")),
                        (self.a, Strand("ACG", "dna", "A")), ()):
            with self.subTest(species=species), self.assertRaises(ValueError):
                enumerate_complexes(species)

    def test_enumeration_size_validation_and_one_species(self):
        for size in (True, 0, -1, 2.0, "2", None):
            with self.subTest(size=size), self.assertRaises(ValueError):
                enumerate_complexes((self.a,), size)
        complexes = enumerate_complexes((self.a,), 5)
        self.assertEqual(tuple(c.symmetry_factor for c in complexes), (1, 2, 3, 4, 5))


class ConditionsTests(unittest.TestCase):
    def test_defaults_and_physical_domain_do_not_invent_model_limits(self):
        conditions = Conditions()
        self.assertEqual((conditions.celsius, conditions.sodium, conditions.magnesium), (37, .12, 0))
        self.assertEqual(conditions.kelvin, 310.15)
        self.assertAlmostEqual(Conditions(celsius=-100, sodium=0, magnesium=.01).kelvin, 173.15)
        self.assertEqual(Conditions(celsius=150, sodium=2).celsius, 150)
        with self.assertRaises(FrozenInstanceError):
            conditions.sodium = 1

    def test_nonfinite_non_numeric_and_unphysical_values_are_rejected(self):
        for field in ("celsius", "sodium", "magnesium"):
            for value in (math.nan, math.inf, -math.inf, True, "37", None):
                with self.subTest(field=field, value=value), self.assertRaises(ValueError):
                    Conditions(**{field: value})
        for options in ({"celsius": -273.15}, {"celsius": -274},
                        {"sodium": -.001}, {"magnesium": -.001}):
            with self.subTest(options=options), self.assertRaises(ValueError):
                Conditions(**options)

    def test_magnesium_is_rejected_at_mixed_ensemble_boundary(self):
        dna = Strand("ACG", "dna", "d")
        rna = Strand("ACG", "rna", "r")
        magnesium = Conditions(magnesium=.001)
        magnesium.validate_for(ComplexSpec((dna, dna)))
        magnesium.validate_for(ComplexSpec((rna,)))
        mixed = ComplexSpec((dna, rna))
        Conditions().validate_for(mixed)
        with self.assertRaisesRegex(ValueError, "magnesium=0"):
            magnesium.validate_for(mixed)
        with self.assertRaises(ValueError):
            partition(complex_spec=mixed, conditions=magnesium)


class ResultTests(unittest.TestCase):
    def test_provenance_is_immutable_and_requires_explicit_conventions(self):
        source = provenance(backend=" reference ")
        self.assertEqual(source.backend, "reference")
        with self.assertRaises(FrozenInstanceError):
            source.version = "2"
        for field in ("backend", "version", "model", "parameter_digest", "ensemble", "standard_state"):
            with self.subTest(field=field), self.assertRaises(ValueError):
                provenance(**{field: " "})
        for digest in ("/path/parameters.json", "sha256:", ":abcd", "sha256:ab cd"):
            with self.subTest(digest=digest), self.assertRaises(ValueError):
                provenance(parameter_digest=digest)

    def test_partition_records_exact_thermal_convention(self):
        result = partition()
        self.assertAlmostEqual(result.partition_function, math.exp(2))
        self.assertFalse(result.is_empty)
        result.require_nonempty()
        for gas_constant in (1.98717e-3, 1.98720425864083e-3):
            rt = gas_constant * result.conditions.kelvin
            changed = replace(result, thermal_energy_kcal_mol=rt,
                              free_energy_kcal_mol=-rt * result.log_partition)
            self.assertEqual(changed.thermal_energy_kcal_mol, rt)
        with self.assertRaises(FrozenInstanceError):
            result.log_partition = 3

    def test_partition_rejects_divergent_free_energy(self):
        with self.assertRaisesRegex(ValueError, "free energy must equal"):
            partition(free_energy_kcal_mol=-1.19)
        partition(free_energy_kcal_mol=-1.2 + 1e-9)
        with self.assertRaises(ValueError):
            partition(free_energy_kcal_mol=-1.2 + 1e-6)

    def test_partition_preserves_constraint_identity_and_normalizes_unresolved(self):
        full = partition()
        self.assertIsNone(full.constraints)
        explicit = partition(constraints=PairingConstraints.unconstrained(2))
        self.assertIsNone(explicit.constraints)
        self.assertEqual(explicit, full)
        self.assertEqual(hash(explicit), hash(full))
        state = PairingConstraints.unconstrained(2).assign(0, 0)
        conditional = partition(constraints=state)
        self.assertIs(conditional.constraints, state)
        self.assertNotEqual(conditional, full)
        self.assertEqual(conditional.provenance, full.provenance)
        with self.assertRaises(FrozenInstanceError):
            conditional.constraints = None
        for invalid in ((-1, -1), "..", PairingConstraints.unconstrained(1), True):
            with self.subTest(constraints=invalid), self.assertRaises(ValueError):
                partition(constraints=invalid)

    def test_empty_ensemble_is_explicit_and_has_no_normalized_observables(self):
        result = partition(log_partition=-math.inf, free_energy_kcal_mol=math.inf)
        self.assertTrue(result.is_empty)
        self.assertEqual(result.partition_function, 0)
        with self.assertRaises(EmptyEnsembleError):
            result.require_nonempty()
        with self.assertRaises(EmptyEnsembleError):
            PairProbabilities(result, ((1, 0), (0, 1)))

    def test_only_the_empty_ensemble_may_have_nonfinite_energy_and_log_z(self):
        for field in ("log_partition", "free_energy_kcal_mol", "thermal_energy_kcal_mol"):
            for value in (math.nan, math.inf, -math.inf, True, None, "1"):
                with self.subTest(field=field, value=value), self.assertRaises(ValueError):
                    partition(**{field: value})
        for rt in (0, -1):
            with self.subTest(rt=rt), self.assertRaises(ValueError):
                partition(thermal_energy_kcal_mol=rt)
        with self.assertRaises(ValueError):
            partition(log_partition=-math.inf, free_energy_kcal_mol=-math.inf)

    def test_stable_log_partition_does_not_conflate_underflow_with_empty(self):
        tiny = partition(log_partition=-1000, free_energy_kcal_mol=600)
        self.assertEqual(tiny.partition_function, 0)
        self.assertFalse(tiny.is_empty)
        tiny.require_nonempty()
        large = partition(log_partition=1000, free_energy_kcal_mol=-600)
        with self.assertRaises(OverflowError):
            _ = large.partition_function

    def test_probabilities_copy_input_and_keep_partition_context(self):
        rows = [[.25, .75], [.75, .25]]
        result = partition()
        probabilities = PairProbabilities(result, rows)
        rows[0][0] = 0
        self.assertIs(probabilities.partition, result)
        self.assertEqual(probabilities.matrix, ((.25, .75), (.75, .25)))
        with self.assertRaises(FrozenInstanceError):
            probabilities.matrix = ()
        with self.assertRaises(TypeError):
            probabilities.matrix[0][0] = 0

    def test_probabilities_reject_wrong_shape_nonfinite_and_invalid_values(self):
        invalid = ((), ((1,),), ((1, 0), (0,)), None,
                   ((math.nan, 0), (0, 1)), ((math.inf, 0), (0, 1)),
                   ((1.1, -.1), (-.1, 1.1)), ((True, 0), (0, 1)))
        for matrix in invalid:
            with self.subTest(matrix=matrix), self.assertRaises(ValueError):
                PairProbabilities(partition(), matrix)

    def test_probabilities_require_unpaired_diagonal_row_mass_and_symmetry(self):
        for matrix in (((0, .75), (.75, 0)), ((.4, .6), (.5, .5))):
            with self.subTest(matrix=matrix), self.assertRaises(ValueError):
                PairProbabilities(partition(), matrix)
        self.assertEqual(PairProbabilities(partition(), ((1, 0), (0, 1))).matrix,
                         ((1.0, 0.0), (0.0, 1.0)))

    def test_results_require_typed_context(self):
        for field in ("complex_spec", "conditions", "provenance"):
            with self.subTest(field=field), self.assertRaises(ValueError):
                partition(**{field: None})
        with self.assertRaises(ValueError):
            PairProbabilities(None, ((1, 0), (0, 1)))


if __name__ == "__main__":
    unittest.main()
