"""Independent motif assembly versus fully constrained native partitions.

These validate ensemble bookkeeping and parameter/API usage, not experimental
accuracy. No reference structure energy is obtained from integer eval_structure.
"""
import math
from types import SimpleNamespace
import unittest
from unittest.mock import patch

from aptaswitch_core.thermodynamics.enumeration import ExhaustiveBackend, enumerate_structures
from aptaswitch_core.thermodynamics.loops import decompose_structure
from aptaswitch_core.thermodynamics.structures import PairingConstraints, format_structure, parse_structure
from aptaswitch_core.thermodynamics.types import ComplexSpec, Conditions, Strand
from aptaswitch_core.thermodynamics.vienna import UnsupportedChemistryError, ViennaBackend
from aptaswitch_core.thermodynamics.vienna_motifs import ViennaMotifModel

try:
    import RNA
except ImportError:
    RNA = None
NATIVE = RNA is not None and RNA.__version__ == "2.7.2"
CONDITIONS = tuple(Conditions(celsius=t, sodium=s) for t in (25, 37, 55) for s in (.12, 1))


def spec(sequences, material="rna"):
    return ComplexSpec(tuple(Strand(s, material, f"s{i}") for i, s in enumerate(sequences)))


def compatible_structure(structure, material="rna"):
    lengths = tuple(map(len, structure.split("+")))
    state = parse_structure(structure, lengths)
    sequence = ["A"] * sum(lengths)
    for i, j in enumerate(state.partners):
        if j > i:
            sequence[i], sequence[j] = "G", "C"
    strands, offset = [], 0
    for length in lengths:
        strands.append("".join(sequence[offset:offset+length]))
        offset += length
    return spec(strands, material), state


class MotifGuardsTests(unittest.TestCase):
    def model(self, material="rna"):
        native = SimpleNamespace(__version__="2.7.2", parameter_set_rna_turner2004="rna",
                                 parameter_set_dna_mathews2004="dna")
        with patch("aptaswitch_core.thermodynamics.vienna_motifs.identify_vienna_runtime",
                   return_value=SimpleNamespace(version_label="fake-2.7.2")), \
                patch("aptaswitch_core.thermodynamics.vienna_motifs.importlib.import_module", return_value=native):
            return ViennaMotifModel(material)

    def test_optional_import_and_exact_version_guard(self):
        target = "aptaswitch_core.thermodynamics.vienna_motifs.importlib.import_module"
        with patch(target, side_effect=ImportError), self.assertRaises(RuntimeError):
            ViennaMotifModel()
        with patch(target, return_value=SimpleNamespace(__version__="2.7.3")), self.assertRaises(RuntimeError):
            ViennaMotifModel()
        with self.assertRaises(ValueError):
            ViennaMotifModel("mixed")

    def test_pair_rules_are_explicit_including_dna_gt(self):
        for chemistry, sequence, expected in (("rna", "ACGU", {(0, 3), (1, 2), (2, 3)}),
                                             ("dna", "ACGT", {(0, 3), (1, 2), (2, 3)})):
            model, current = self.model(chemistry), spec((sequence,), chemistry)
            observed = {(i, j) for i in range(4) for j in range(i+1, 4) if model.pair_allowed(current, i, j)}
            self.assertEqual(observed, expected)
            for i, j in observed:
                self.assertTrue(model.pair_allowed(current, j, i))
            self.assertFalse(model.pair_allowed(current, 0, 0))
            with self.assertRaises(ValueError):
                model.pair_allowed(current, True, 1)

    def test_chemistry_salt_and_partial_structure_guards(self):
        model = self.model()
        mixed = ComplexSpec((Strand("G", "rna", "r"), Strand("C", "dna", "d")))
        with self.assertRaises(UnsupportedChemistryError):
            model.energy(mixed, Conditions(), "(+)")
        for bad in (Conditions(sodium=0), Conditions(magnesium=.001), []):
            with self.assertRaises(ValueError):
                model.thermal_energy(bad)
        with self.assertRaises(ValueError):
            model.energy(spec(("GC",)), Conditions(), PairingConstraints.unconstrained(2))

    def test_zero_weight_topologies_never_reach_native_apis(self):
        model = self.model()
        cases = [(spec(("G", "C")), ".+."), (spec(("A", "A")), "(+)"),
                 (spec(("GAC",)), "(.)"), compatible_structure("(" + "."*31 + "(...))")]
        for current, state in cases:
            with self.subTest(state=state):
                self.assertEqual(model.energy(current, Conditions(), state), math.inf)

    def test_native_weight_validation(self):
        self.assertEqual(ViennaMotifModel._energy_of_weight(0, .6), math.inf)
        for value in (-1, math.nan, math.inf):
            with self.assertRaises(ArithmeticError):
                ViennaMotifModel._energy_of_weight(value, .6)


@unittest.skipUnless(NATIVE, "optional ViennaRNA 2.7.2 bindings are unavailable")
class NativeMotifTests(unittest.TestCase):
    def setUp(self):
        self.native = ViennaBackend()

    def test_complete_small_ensembles_reconstruct_pf_pairs_and_weight_minimum(self):
        fixtures = (("GCGAAAGC",), ("GCGC", "GCGU"), ("GC",)*3, ("GC",)*4, ("GC",)*5)
        for material in ("rna", "dna"):
            model = ViennaMotifModel(material)
            exhaustive = ExhaustiveBackend(pair_allowed=model.pair_allowed, energy=model.energy,
                thermal_energy=model.thermal_energy, provenance=model.provenance, min_hairpin_size=3)
            for sequences in fixtures:
                if material == "dna":
                    sequences = tuple(s.replace("U", "T") for s in sequences)
                current = spec(sequences, material)
                for conditions in CONDITIONS:
                    with self.subTest(material=material, sequences=sequences, conditions=conditions):
                        observed = exhaustive.evaluate(current, conditions)
                        expected = self.native.pairs(current, conditions)
                        self.assertEqual(observed.partition.provenance.parameter_digest,
                                         expected.partition.provenance.parameter_digest)
                        self.assertAlmostEqual(observed.partition.free_energy_kcal_mol,
                                               expected.partition.free_energy_kcal_mol, delta=4e-6)
                        self.assertAlmostEqual(observed.partition.thermal_energy_kcal_mol,
                                               expected.partition.thermal_energy_kcal_mol, delta=1e-12)
                        for row, wanted in zip(observed.pairs().matrix, expected.matrix):
                            for actual, target in zip(row, wanted):
                                self.assertAlmostEqual(actual, target, delta=2e-8)
                        native_energies = [self.native.partition(current, conditions, state).free_energy_kcal_mol
                                           for state in observed.structures]
                        for actual, target in zip(observed.energies_kcal_mol, native_energies):
                            self.assertAlmostEqual(actual, target, delta=4e-6)
                        self.assertAlmostEqual(observed.mfe()[0], min(native_energies), delta=4e-6)

    def test_each_loop_family_and_nicked_faces_match_constrained_pf(self):
        structures = ("(...)" , "(....)", "(......)", "((...))", "(.(...))", "(....(...))",
                      "(.(...).)", "(.(...)..)", "(..(...)..)", "(..(...)...)",
                      "(.(...)....)", "(...(...)...)", "((...)(...))",
                      "(..(...)...(....)....)", "(.(...)+(...).)")
        kinds = set()
        for material in ("rna", "dna"):
            model = ViennaMotifModel(material)
            for structure in structures:
                current, state = compatible_structure(structure, material)
                kinds.update(loop.kind for loop in decompose_structure(current, structure))
                for conditions in CONDITIONS:
                    with self.subTest(material=material, structure=structure, conditions=conditions):
                        energy = model.energy(current, conditions, state)
                        expected = self.native.partition(current, conditions, state).free_energy_kcal_mol
                        self.assertAlmostEqual(energy, expected, delta=4e-6)
        self.assertEqual(kinds, {"exterior", "hairpin", "stack", "bulge", "interior", "multiloop"})

    def test_native_structure_evaluator_is_not_used(self):
        model = ViennaMotifModel()
        current = spec(("GCGAAAGC",))
        state = parse_structure("((....))", (8,))
        with patch.object(RNA, "fold_compound", side_effect=AssertionError("No native structure evaluator")):
            self.assertTrue(math.isfinite(model.energy(current, Conditions(), state)))

    def test_parameter_snapshots_survive_global_material_reloads(self):
        models = [ViennaMotifModel(m) for m in ("rna", "dna")]
        c = Conditions(celsius=25, sodium=.12)
        values = [m.energy(spec(("GCGAAAGC",), m.material), c, "((....))") for m in models]
        self.assertNotEqual(values[0], values[1])
        for model, target in zip(models, values):
            self.assertEqual(model.energy(spec(("GCGAAAGC",), model.material), c, "((....))"), target)

    def test_fresh_snapshots_resist_changed_process_defaults(self):
        # Fresh model instances prevent the Python LRU from hiding native
        # global-default contamination. Alternate DNA/RNA loads at fixed RT.
        c = Conditions(celsius=25, sodium=.12)
        fixtures = [compatible_structure(structure, material)
                    for material in ("dna", "rna")
                    for structure in ("(....)", "(.(...).)", "((...)(...))")]
        expected = []
        for current, state in fixtures:
            model = ViennaMotifModel(current.strands[0].material)
            expected.append((model.energy(current, c, state), model.thermal_energy(c), model.provenance))
        changed = {"pf_smooth": 0, "betaScale": 2.0, "special_hp": 0,
                   "noGUclosure": 1, "logML": 1, "energy_set": 1,
                   "saltDPXInit": 700}
        original = {field: getattr(RNA.cvar, field) for field in changed}
        try:
            for field, value in changed.items():
                setattr(RNA.cvar, field, value)
            self.assertEqual(RNA.md().pf_smooth, 0)
            for (current, state), (energy, thermal, provenance) in zip(fixtures, expected):
                with self.subTest(material=current.strands[0].material, state=state):
                    model = ViennaMotifModel(current.strands[0].material)
                    self.assertEqual(model.energy(current, c, state), energy)
                    self.assertEqual(model.thermal_energy(c), thermal)
                    self.assertEqual(model.provenance, provenance)
                    self.assertIn("pf_smooth=1", model.provenance.ensemble)
            self.assertEqual(RNA.md().pf_smooth, 0)
        finally:
            for field, value in original.items():
                setattr(RNA.cvar, field, value)

    def test_symmetry_is_one_uniform_labelled_state_factor(self):
        for material in ("rna", "dna"):
            model = ViennaMotifModel(material)
            for count in (2, 3, 4, 5):
                current = spec(("GC",)*count, material)
                c = Conditions(celsius=37, sodium=1.021)
                correction = model.thermal_energy(c) * math.log(count)
                self.assertEqual(current.symmetry_factor, count)
                for state in enumerate_structures(current, pair_allowed=model.pair_allowed,
                                                   min_hairpin_size=3, max_nucleotides=18):
                    total = model.energy(current, c, state)
                    native = self.native.partition(current, c, state).free_energy_kcal_mol
                    self.assertAlmostEqual(total, native, delta=4e-6)
                    self.assertGreater(correction, 0)
                self.assertIn("labelled-state", model.provenance.standard_state)

    def test_large_hairpins_extrapolate_and_native_salt_slot_bug_is_reported(self):
        for material in ("rna", "dna"):
            model = ViennaMotifModel(material)
            for count in (29, 30, 31, 40):
                current, state = compatible_structure("(" + "."*count + ")", material)
                for conditions in (Conditions(sodium=.12), Conditions(sodium=1.021)):
                    params = model._parameters(conditions)
                    has_slot_bug = (conditions.sodium != RNA.MODEL_DEFAULT_SALT
                                    and RNA.doubleP_getitem(params.expSaltLoop, 31) == 0)
                    if count == 30 and has_slot_bug:
                        with self.assertRaisesRegex(ArithmeticError, "31-backbone"):
                            model.energy(current, conditions, state)
                        with self.assertRaises(ArithmeticError):
                            self.native.partition(current, conditions, state)
                    else:
                        energy = model.energy(current, conditions, state)
                        self.assertTrue(math.isfinite(energy))
                        if current.n_nucleotides >= 32 and has_slot_bug:
                            # The complete native backend conservatively blocks
                            # every long request on affected stock builds: an
                            # unconstrained PF can silently omit other states.
                            with self.assertRaises(ArithmeticError):
                                self.native.partition(current, conditions, state)
                        else:
                            self.assertAlmostEqual(energy,
                                self.native.partition(current, conditions, state).free_energy_kcal_mol, delta=4e-6)

    def test_same_salt_table_boundary_affects_29_base_interiors_and_bulges(self):
        for structure in ("(" + "."*29 + "(...))", "(" + "."*14 + "(...)" + "."*15 + ")"):
            current, state = compatible_structure(structure)
            model = ViennaMotifModel()
            conditions = Conditions(sodium=.12)
            has_slot_bug = RNA.doubleP_getitem(model._parameters(conditions).expSaltLoop, 31) == 0
            if has_slot_bug:
                with self.assertRaisesRegex(ArithmeticError, "31-backbone"):
                    model.energy(current, conditions, state)
                with self.assertRaises(ArithmeticError):
                    self.native.partition(current, conditions, state)
            else:
                # Audited patched builds initialize this slot. The public
                # constrained partition must agree with independent loop sums.
                self.assertAlmostEqual(model.energy(current, conditions, state),
                    self.native.partition(current, conditions, state).free_energy_kcal_mol, delta=4e-6)
            baseline = Conditions(sodium=1.021)
            self.assertAlmostEqual(model.energy(current, baseline, state),
                self.native.partition(current, baseline, state).free_energy_kcal_mol, delta=4e-6)


if __name__ == "__main__":
    unittest.main()
