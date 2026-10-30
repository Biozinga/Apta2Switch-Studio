"""Homogeneous adapter contracts and optional ViennaRNA 2.7.2 integration.

The exhaustive reference enumerates structures independently, but deliberately
uses Vienna energies: these tests check ensemble bookkeeping, not the physical
accuracy of the parameter sets. Fully constrained PF energies avoid comparing
integer MFE energies with floating-point PF parameters away from 37 C.
"""
from collections import Counter
import hashlib
import math
from types import SimpleNamespace
import unittest
from unittest.mock import Mock, patch

from aptaswitch_core.thermodynamics import ComplexSpec, Conditions, EmptyEnsembleError, Strand
from aptaswitch_core.thermodynamics.structures import (
    PairingConstraints, is_connected, parse_structure,
)
from aptaswitch_core.thermodynamics.vienna import UnsupportedChemistryError, ViennaBackend

try:
    import RNA
except ImportError:
    RNA = None
from tests.native_reference.reference_checks import (
    boltzmann_summary, connected, enum_structures,
)

NATIVE = RNA is not None and RNA.__version__ == "2.7.2"
CONDITIONS = tuple(Conditions(celsius=t, sodium=s) for t in (25, 37, 55) for s in (.12, 1))


def spec(sequences, material="rna"):
    # Explicit fixture generation for two different chemistries, never an
    # adapter-side transcription of an input sequence.
    return ComplexSpec(tuple(Strand(sequence, material, f"s{i}")
                             for i, sequence in enumerate(sequences)))


def fake_native(*, mfe=0.0, mfe_structure="..", energy=0.0, upper=None):
    if upper is None:
        upper = ((0, 0, 0), (0, 0, .4), (0, 0, 0))
    fc = SimpleNamespace(
        mfe=Mock(return_value=(mfe_structure, mfe)), pf=Mock(return_value=("..", energy)),
        bpp=Mock(return_value=upper), exp_params=SimpleNamespace(kT=616.3207755),
        exp_params_rescale=Mock(), hc_add_up=Mock(), hc_add_bp=Mock(),
    )
    native = SimpleNamespace(
        __version__="2.7.2", INF=10000000,
        params_load_RNA_Turner2004=Mock(), params_load_DNA_Mathews2004=Mock(),
        parameter_set_rna_turner2004="rna-test", parameter_set_dna_mathews2004="dna-test",
        md=SimpleNamespace, fold_compound=Mock(return_value=fc),
        CONSTRAINT_CONTEXT_ALL_LOOPS=63, CONSTRAINT_CONTEXT_ENFORCE=64,
    )
    return native, fc


class ViennaGuardsTests(unittest.TestCase):
    def backend(self, **kwargs):
        native, fc = fake_native(**kwargs)
        with patch("aptaswitch_core.thermodynamics.vienna.identify_vienna_runtime",
                   return_value=SimpleNamespace(version_label="fixture:simulated-native")), \
                patch("aptaswitch_core.thermodynamics.vienna.importlib.import_module", return_value=native):
            backend = ViennaBackend()
        return backend, native, fc

    def test_missing_and_unvalidated_native_versions_fail_explicitly(self):
        target = "aptaswitch_core.thermodynamics.vienna.importlib.import_module"
        with patch(target, side_effect=ImportError), self.assertRaises(RuntimeError):
            ViennaBackend()
        for version in ("2.7.1", "2.8.0"):
            with patch(target, return_value=SimpleNamespace(__version__=version)), self.assertRaises(RuntimeError):
                ViennaBackend()

    def test_mixed_chemistry_is_rejected_before_native_encoding(self):
        backend, native, _ = self.backend()
        mixed = ComplexSpec((Strand("ACG", "rna", "r"), Strand("ACG", "dna", "d")))
        for operation in (backend.partition, backend.pairs):
            with self.assertRaises(UnsupportedChemistryError):
                operation(mixed, Conditions())
        native.fold_compound.assert_not_called()
        native.params_load_RNA_Turner2004.assert_not_called()
        native.params_load_DNA_Mathews2004.assert_not_called()

    def test_unsupported_salts_and_wrong_constraint_length_fail_before_encoding(self):
        backend, native, _ = self.backend()
        for conditions in (Conditions(sodium=0), Conditions(magnesium=.001)):
            with self.assertRaises(ValueError):
                backend.partition(spec(("GC",)), conditions)
        with self.assertRaises(ValueError):
            backend.partition(spec(("GC",)), Conditions(), PairingConstraints.unconstrained(3))
        native.fold_compound.assert_not_called()

    def test_nonfinite_native_mfe_is_not_an_empty_ensemble(self):
        for bad in (math.nan, math.inf, -math.inf):
            with self.subTest(value=bad):
                backend, _, fc = self.backend(mfe=bad)
                with self.assertRaises(ArithmeticError):
                    backend.partition(spec(("G", "C")), Conditions())
                fc.pf.assert_not_called()

    def test_nan_or_negative_infinite_pf_is_a_numerical_failure(self):
        for bad in (math.nan, -math.inf):
            with self.subTest(value=bad):
                backend, _, _ = self.backend(energy=bad)
                with self.assertRaises(ArithmeticError):
                    backend.partition(spec(("G", "C")), Conditions())

    def test_positive_infinite_pf_is_native_empty_even_when_mfe_is_finite(self):
        backend, _, fc = self.backend(energy=math.inf)
        self.assertTrue(backend.partition(spec(("A", "A")), Conditions()).is_empty)
        with self.assertRaises(EmptyEnsembleError):
            backend.pairs(spec(("A", "A")), Conditions())
        fc.bpp.assert_not_called()

    def test_finite_connected_mfe_witness_makes_positive_infinite_pf_a_failure(self):
        for sequences, structure in ((("G", "C"), "()"), (("GC",), "..")):
            backend, _, _ = self.backend(mfe_structure=structure, energy=math.inf)
            with self.subTest(sequences=sequences), self.assertRaises(ArithmeticError):
                backend.partition(spec(sequences), Conditions())
        # An MFE traceback that violates the requested subset proves nothing
        # about its connected partition. Neither does an absent traceback.
        backend, _, _ = self.backend(mfe_structure="()", energy=math.inf)
        unpaired = PairingConstraints.unconstrained(2).assign(0, 0)
        self.assertTrue(backend.partition(spec(("G", "C")), Conditions(), unpaired).is_empty)
        backend, _, _ = self.backend(mfe_structure="", energy=math.inf)
        self.assertTrue(backend.partition(spec(("A", "A")), Conditions()).is_empty)

    def test_bad_native_pair_values_and_dimensions_are_rejected(self):
        for bad in (math.nan, math.inf, -math.inf, -.01, 1.01):
            with self.subTest(value=bad):
                backend, _, _ = self.backend(upper=((0, 0, 0), (0, 0, bad), (0, 0, 0)))
                with self.assertRaises(ArithmeticError):
                    backend.pairs(spec(("G", "C")), Conditions())
        backend, _, _ = self.backend(upper=((0, 0), (0, 0)))
        with self.assertRaises(ArithmeticError):
            backend.pairs(spec(("G", "C")), Conditions())

    def test_empty_mfe_sentinel_does_not_call_pf_or_normalize_pairs(self):
        backend, _, fc = self.backend(mfe=100000.0)
        result = backend.partition(spec(("A", "A")), Conditions())
        self.assertTrue(result.is_empty)
        self.assertEqual((result.log_partition, result.free_energy_kcal_mol), (-math.inf, math.inf))
        with self.assertRaises(EmptyEnsembleError):
            backend.pairs(spec(("A", "A")), Conditions())
        fc.pf.assert_not_called()
        fc.bpp.assert_not_called()

    def test_parameter_provenance_and_thermal_convention_are_explicit(self):
        backend, native, _ = self.backend()
        for material, text, model in (("rna", "rna-test", "RNA-Turner2004"),
                                      ("dna", "dna-test", "DNA-Mathews2004")):
            result = backend.partition(spec(("GC",), material), Conditions())
            self.assertEqual(result.provenance.model, model)
            self.assertEqual(result.provenance.version, "fixture:simulated-native")
            self.assertEqual(result.provenance.parameter_digest, "sha256:" + hashlib.sha256(text.encode()).hexdigest())
            self.assertEqual(result.thermal_energy_kcal_mol, .6163207755)
            self.assertIn("native cyclic symmetry", result.provenance.standard_state)
        native.params_load_RNA_Turner2004.assert_called_once()
        native.params_load_DNA_Mathews2004.assert_called_once()

    def test_model_flags_are_explicit_and_initialized_after_loading_tables(self):
        backend, native, _ = self.backend()
        events = []
        native.params_load_DNA_Mathews2004.side_effect = lambda: events.append("load-dna")

        def make_details(**kwargs):
            events.append("details")
            return SimpleNamespace(**kwargs)

        native.md = Mock(side_effect=make_details)
        for with_pairs in (False, True):
            operation = backend.pairs if with_pairs else backend.partition
            result = operation(spec(("G", "C"), "dna"), Conditions(celsius=25, sodium=.12))
            md = native.fold_compound.call_args.args[1]
            self.assertEqual(events[-2:], ["load-dna", "details"])
            self.assertEqual(md.temperature, 25)
            self.assertEqual(md.salt, .12)
            for field, expected in {
                "pf_smooth": 1, "betaScale": 1, "special_hp": 1,
                "noGUclosure": 0, "noGU": 0, "noLP": 0, "dangles": 0,
                "logML": 0, "circ": 0, "gquad": 0, "energy_set": 0,
                "max_bp_span": -1, "min_loop_size": 3,
                "compute_bpp": int(with_pairs),
            }.items():
                self.assertEqual(getattr(md, field), expected, field)
            partition = result.partition if with_pairs else result
            self.assertIn("pf_smooth=1", partition.provenance.ensemble)


@unittest.skipUnless(NATIVE, "optional ViennaRNA 2.7.2 bindings are unavailable")
class ViennaNativeTests(unittest.TestCase):
    def setUp(self):
        self.backend = ViennaBackend()

    def assert_matrix_close(self, actual, expected, tolerance=2e-6):
        self.assertEqual(len(actual), len(expected))
        for row, wanted in zip(actual, expected):
            self.assertEqual(len(row), len(wanted))
            for value, target in zip(row, wanted):
                self.assertAlmostEqual(value, target, delta=tolerance)

    def test_results_resist_changed_process_model_defaults(self):
        fixtures = (spec(("GCGAAAGC",)), spec(("GCGAAAGC",), "dna"),
                    spec(("GCGC", "GCGU")), spec(("GCGC", "GCGT"), "dna"))
        conditions = Conditions(celsius=25, sodium=.12)
        expected = [self.backend.pairs(current, conditions) for current in fixtures]
        changed = {
            "pf_smooth": 0, "betaScale": 2.0, "dangles": 2, "special_hp": 0,
            "noLP": 1, "noGU": 1, "noGUclosure": 1, "logML": 1,
            "circ": 1, "gquad": 1, "energy_set": 1, "max_bp_span": 4,
            "min_loop_size": 5, "backtrack_type": "C", "compute_bpp": 0,
            "saltDPXInit": 700,
        }
        original = {field: getattr(RNA.cvar, field) for field in changed}
        try:
            for field, value in changed.items():
                setattr(RNA.cvar, field, value)
            # This must change real native defaults, not just a Python mock.
            self.assertEqual(RNA.md().pf_smooth, 0)
            self.assertEqual(RNA.md().energy_set, 1)
            for current, wanted in zip(fixtures, expected):
                with self.subTest(material=current.strands[0].material, n=len(current.strands)):
                    actual = ViennaBackend().pairs(current, conditions)
                    self.assertEqual(actual.partition, wanted.partition)
                    self.assert_matrix_close(actual.matrix, wanted.matrix, tolerance=1e-14)
            # Constructing an adapter must not reset unrelated callers' defaults.
            self.assertEqual(RNA.md().pf_smooth, 0)
            self.assertEqual(RNA.md().energy_set, 1)
        finally:
            for field, value in original.items():
                setattr(RNA.cvar, field, value)

    def test_connected_ensembles_one_to_five_strands_reconstruct_pf_and_pairs(self):
        fixtures = (("GCGAAAGC",), ("GCGC", "GCGU"), ("GC",) * 3,
                    ("GC",) * 4, ("GC",) * 5)
        for material in ("rna", "dna"):
            for sequences, count in zip(fixtures, (5, 25, 4, 5, 6)):
                if material == "dna":
                    sequences = tuple(sequence.replace("U", "T") for sequence in sequences)
                current = spec(sequences, material)
                lengths = tuple(map(len, sequences))
                states = [s for s in enum_structures(current.sequence, current.cutpoints)
                          if connected(s, current.cutpoints)]
                self.assertEqual(len(states), count)
                for conditions in CONDITIONS:
                    with self.subTest(material=material, n=len(sequences), conditions=conditions):
                        full = self.backend.pairs(current, conditions)
                        energies = []
                        for structure in states:
                            constraint = parse_structure(structure, lengths)
                            isolated = self.backend.pairs(current, conditions, constraint)
                            energies.append(isolated.partition.free_energy_kcal_mol)
                            for i, j in enumerate(constraint.partners):
                                self.assertAlmostEqual(isolated.matrix[i][j], 1, delta=1e-8)
                        expected = boltzmann_summary(states, energies, full.partition.thermal_energy_kcal_mol)
                        # SWIG returns PF energies as a C float; the structural
                        # probabilities themselves retain native double precision.
                        self.assertAlmostEqual(full.partition.free_energy_kcal_mol,
                                               expected["free_energy"], delta=4e-6)
                        self.assert_matrix_close(full.matrix, expected["pair_probabilities"])

    def test_native_rotational_symmetry_is_included_exactly_once(self):
        conditions = Conditions(celsius=37, sodium=1.021)
        for material in ("rna", "dna"):
            for number in (2, 3, 4, 5):
                with self.subTest(material=material, n=number):
                    current = spec(("GC",) * number, material)
                    result = self.backend.pairs(current, conditions)
                    md = RNA.md()
                    md.temperature, md.salt, md.dangles = 37, 1.021, 0
                    md.noGU, md.noLP, md.uniq_ML = False, False, 1
                    fc = RNA.fold_compound("&".join(s.sequence for s in current.strands), md)
                    states = [s for s in enum_structures(current.sequence, current.cutpoints)
                              if connected(s, current.cutpoints)]
                    energies = [fc.eval_structure(s) for s in states]
                    labeled = boltzmann_summary(states, energies, result.partition.thermal_energy_kcal_mol)
                    correction = result.partition.thermal_energy_kcal_mol * math.log(number)
                    self.assertEqual(current.symmetry_factor, number)
                    self.assertAlmostEqual(result.partition.free_energy_kcal_mol,
                                           labeled["free_energy"] + correction, delta=3e-6)
                    self.assert_matrix_close(result.matrix, labeled["pair_probabilities"])

    def test_conditional_marginals_equal_ratios_of_partition_functions(self):
        for material, sequences in (("rna", ("GCGC", "GCGU")), ("dna", ("GCGC", "GCGT"))):
            current = spec(sequences, material)
            for conditions in CONDITIONS:
                for first_choice in (None, (0, 0), (0, 7)):
                    state = PairingConstraints.unconstrained(current.n_nucleotides)
                    if first_choice is not None:
                        state = state.assign(*first_choice)
                    with self.subTest(material=material, conditions=conditions, first=first_choice):
                        parent = self.backend.pairs(current, conditions, state)
                        self.assertEqual(parent.partition.constraints, None if first_choice is None else state)
                        i = state.partners.index(-1)
                        total = 0.0
                        for j, p in enumerate(parent.matrix[i]):
                            if state.partners[j] != -1:
                                continue
                            # A crossing partner is outside the typed matching.
                            try:
                                child = state.assign(i, j)
                            except ValueError:
                                self.assertEqual(p, 0)
                                continue
                            result = self.backend.partition(current, conditions, child)
                            self.assertEqual(result.constraints, child)
                            ratio = math.exp(result.log_partition - parent.partition.log_partition)
                            self.assertAlmostEqual(p, ratio, delta=2e-6)
                            total += ratio
                        self.assertAlmostEqual(total, 1, delta=3e-6)

    def test_empty_connected_ensembles_have_no_marginals_or_samples(self):
        for material in ("rna", "dna"):
            for number in range(2, 6):
                current = spec(("A",) * number, material)
                self.assertTrue(self.backend.partition(current, Conditions()).is_empty)
                with self.assertRaises(EmptyEnsembleError):
                    self.backend.pairs(current, Conditions())
                with self.assertRaises(EmptyEnsembleError):
                    self.backend.sampler(current, Conditions(), seed=5).sample_one()

    def test_impossible_constraints_do_not_enable_pairs_or_get_ignored(self):
        for material in ("rna", "dna"):
            for sequences, forced in ((("AAAAA",), (0, 4)), (("GC",), (0, 1)),
                                      (("A", "A"), (0, 1)), (("G", "C"), (0, 0))):
                current = spec(sequences, material)
                state = PairingConstraints.unconstrained(current.n_nucleotides).assign(*forced)
                with self.subTest(material=material, sequences=sequences, forced=forced):
                    result = self.backend.partition(current, Conditions(), state)
                    self.assertTrue(result.is_empty)
                    self.assertIs(result.constraints, state)
                    with self.assertRaises(EmptyEnsembleError):
                        self.backend.pairs(current, Conditions(), state)
            # Nicks exempt interstrand pairs from the hairpin length minimum.
            pair = PairingConstraints.unconstrained(2).assign(0, 1)
            self.assertEqual(self.backend.pairs(spec(("G", "C"), material), Conditions(), pair).matrix,
                             ((0.0, 1.0), (1.0, 0.0)))

    def test_chemistry_parameter_loading_cannot_leak_between_calls(self):
        rna, dna = spec(("GCGAAAGC",)), spec(("GCGAAAGC",), "dna")
        first = self.backend.pairs(rna, Conditions())
        other = self.backend.pairs(dna, Conditions())
        again = self.backend.pairs(rna, Conditions())
        self.assertEqual(first, again)
        self.assertNotEqual(first.partition.provenance.parameter_digest, other.partition.provenance.parameter_digest)
        self.assertNotEqual(first.partition.free_energy_kcal_mol, other.partition.free_energy_kcal_mol)

    def test_uninitialized_salt_slot_cannot_silently_remove_long_loop_states(self):
        RNA = self.backend.RNA
        for material in ("rna", "dna"):
            for celsius in (25, 37, 55):
                for sodium in (.12, 1.0):
                    conditions = Conditions(celsius=celsius, sodium=sodium)
                    # At 31 nt the affected motifs cannot occur.
                    self.assertFalse(self.backend.partition(
                        spec(("G" + "A" * 29 + "C",), material), conditions).is_empty)
                    md = RNA.md()
                    md.temperature, md.salt = celsius, sodium
                    has_bug = RNA.doubleP_getitem(RNA.exp_param(md).expSaltLoop, 31) == 0
                    current = spec(("G" + "A" * 30 + "C",), material)
                    if has_bug:
                        with self.assertRaisesRegex(ArithmeticError, "salt-loop-31"):
                            self.backend.partition(current, conditions)
                    else:
                        self.assertFalse(self.backend.partition(current, conditions).is_empty)
                # Native baseline bypasses the faulty salt-loop table.
                self.assertFalse(self.backend.partition(
                    spec(("G" + "A" * 30 + "C",), material),
                    Conditions(celsius=celsius, sodium=RNA.MODEL_DEFAULT_SALT)).is_empty)

    def test_native_conditional_sampler_is_reproducible_connected_and_matches_marginals(self):
        current = spec(("GC",) * 3)
        conditions = Conditions()
        expected = self.backend.pairs(current, conditions)
        sampler = self.backend.sampler(current, conditions, seed=712)
        samples = sampler.sample(1000)
        self.assertEqual(samples, self.backend.sampler(current, conditions, seed=712).sample(1000))
        counts = Counter()
        for structure in samples:
            state = parse_structure(structure, (2, 2, 2))
            self.assertTrue(is_connected(state, (2, 2, 2)))
            counts.update(enumerate(state.partners))
        for i, row in enumerate(expected.matrix):
            for j, p in enumerate(row):
                self.assertAlmostEqual(counts[i, j] / len(samples), p, delta=.07)
        self.assertGreater(sampler.cache_hits, 0)


if __name__ == "__main__":
    unittest.main()
