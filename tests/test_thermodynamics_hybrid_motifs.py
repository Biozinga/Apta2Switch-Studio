"""Public duplex tables plus independently composed mixed-motif checks."""
import itertools
import math
from dataclasses import replace
import unittest

from aptaswitch_core.thermodynamics.enumeration import ExhaustiveBackend, enumerate_structures
from aptaswitch_core.thermodynamics.hybrid_motifs import HybridMotifModel
from aptaswitch_core.thermodynamics.loops import decompose_structure
from aptaswitch_core.thermodynamics.structures import PairingConstraints, format_structure, parse_structure
from aptaswitch_core.thermodynamics.types import ComplexSpec, Conditions, Strand
from aptaswitch_core.thermodynamics.vienna_motifs import ViennaMotifModel

try:
    import RNA
except ImportError:
    RNA = None
try:
    from Bio.SeqUtils.MeltingTemp import R_DNA_NN1
except ImportError:
    R_DNA_NN1 = None
NATIVE = RNA is not None and RNA.__version__ == "2.7.2"


def mixed_fixture(structure):
    lengths = tuple(map(len, structure.split("+")))
    state = parse_structure(structure, lengths)
    sequence = ["A"] * sum(lengths)
    for i, j in enumerate(state.partners):
        if j > i:
            sequence[i], sequence[j] = "G", "C"
    strands, offset = [], 0
    for index, n in enumerate(lengths):
        strands.append(Strand("".join(sequence[offset:offset+n]), "rna" if index % 2 == 0 else "dna", str(index)))
        offset += n
    return ComplexSpec(tuple(strands)), state


@unittest.skipUnless(NATIVE, "optional ViennaRNA 2.7.2 bindings are unavailable")
class HybridMotifTests(unittest.TestCase):
    def setUp(self):
        self.model = HybridMotifModel()

    def native_parameters(self, material, temperature):
        return ViennaMotifModel(material)._parameters(Conditions(celsius=temperature, sodium=1.021))

    @staticmethod
    def native_energy(weight, params):
        return -params.kT/1000 * math.log(weight)

    def test_pair_admissibility_preserves_real_chemical_identity(self):
        strands = tuple(Strand(base, material, material+base) for material, alphabet in
                        (("rna", "ACGU"), ("dna", "ACGT")) for base in alphabet)
        spec = ComplexSpec(strands)
        for i, left in enumerate(strands):
            for j, right in enumerate(strands):
                pair = (left.sequence+right.sequence).replace("U", "T")
                wanted = i != j and (pair in ("AT", "TA", "CG", "GC") or
                    left.material == right.material == "rna" and pair in ("GT", "TG"))
                self.assertEqual(self.model.pair_allowed(spec, i, j), wanted)
        self.assertEqual(spec.sequence, "ACGUACGT")
        with self.assertRaises(ValueError):
            self.model.pair_allowed(spec, True, 1)

    def test_domain_and_structure_guards_are_explicit(self):
        spec = ComplexSpec((Strand("G", "rna", "r"), Strand("T", "dna", "d")))
        for c in (Conditions(), Conditions(sodium=1.021), Conditions(sodium=1, magnesium=.001),
                  Conditions(celsius=24, sodium=1), Conditions(celsius=56, sodium=1)):
            with self.assertRaises(ValueError):
                self.model.energy(spec, c, "(+)")
        c = Conditions(sodium=1)
        self.assertEqual(self.model.energy(spec, c, "(+)"), math.inf)  # no fabricated hybrid wobble
        self.assertEqual(self.model.energy(spec, c, ".+."), math.inf)
        with self.assertRaises(ValueError):
            self.model.energy(spec, c, PairingConstraints.unconstrained(2))
        for structure in ("(.)", "("+"."*31+"(+)"+")"):
            current, state = mixed_fixture(structure)
            self.assertEqual(self.model.energy(current, c, state), math.inf)
        self.assertIn("experimental", self.model.provenance.model)
        self.assertIn("Na=1M-reference", self.model.provenance.standard_state)

    def test_homogeneous_endpoint_uses_its_own_continuous_loop_weights(self):
        for material, temperature in itertools.product(("rna", "dna"), (25, 37, 55)):
            mono = ViennaMotifModel(material)
            for sequence, structure in (("GCGAAAGC", "((....))"), ("GAAACGAAAC", "(...)(...)")):
                spec = ComplexSpec((Strand(sequence, material, "s"),))
                actual = self.model.energy(spec, Conditions(celsius=temperature, sodium=1), structure)
                expected = mono.energy(spec, Conditions(celsius=temperature, sodium=1.021), structure)
                self.assertAlmostEqual(actual, expected, delta=1e-12)

    @unittest.skipUnless(R_DNA_NN1 is not None, "optional independent Biopython table unavailable")
    def test_all_sixteen_hybrid_stacks_and_both_strand_orders_match_biopython(self):
        for key, (enthalpy, entropy) in R_DNA_NN1.items():
            if len(key) != 5 or key[2] != "/":
                continue
            rna, dna = key.split("/")
            for temperature in (25, 37, 55):
                c = Conditions(celsius=temperature, sodium=1)
                expected = (enthalpy+R_DNA_NN1["init"][0]
                            - c.kelvin*(entropy+R_DNA_NN1["init"][1])/1000)
                strands = (Strand(rna.replace("T", "U"), "rna", "r"), Strand(dna[::-1], "dna", "d"))
                for order in (strands, strands[::-1]):
                    self.assertAlmostEqual(self.model.energy(ComplexSpec(order), c, "((+))"), expected, delta=1e-12)

    @unittest.skipUnless(R_DNA_NN1 is not None, "optional independent Biopython table unavailable")
    def test_long_duplex_energy_curves_and_rna_orientation(self):
        complement = str.maketrans("ACGT", "TGCA")
        for rna in ("ACGTAGTC", "AAAAAA", "TTTTTT", "CGGATCAGT"):
            bottom = rna.translate(complement)
            keys = [rna[i:i+2]+"/"+bottom[i:i+2] for i in range(len(rna)-1)]
            h = math.fsum(R_DNA_NN1[k][0] for k in ["init", *keys])
            s = math.fsum(R_DNA_NN1[k][1] for k in ["init", *keys])
            spec = ComplexSpec((Strand(rna.replace("T", "U"), "rna", "r"), Strand(bottom[::-1], "dna", "d")))
            structure = "("*len(rna) + "+" + ")"*len(rna)
            for t in (25, 37, 55):
                c = Conditions(celsius=t, sodium=1)
                self.assertAlmostEqual(self.model.energy(spec, c, structure), h-c.kelvin*s/1000, delta=1e-11)
        # Swapping just the chemistry is not a cyclic rotation of the same
        # physical duplex: rAG/dTC and rCU/dGA have different measured stacks.
        first = ComplexSpec((Strand("AG", "rna", "r"), Strand("CT", "dna", "d")))
        swapped = ComplexSpec((Strand("AG", "dna", "d"), Strand("CU", "rna", "r")))
        self.assertNotAlmostEqual(self.model.energy(first, Conditions(sodium=1), "((+))"),
                                  self.model.energy(swapped, Conditions(sodium=1), "((+))"))

    def test_mixed_small_interiors_use_complete_projection_mean_not_missing_generic_table(self):
        structures = ("(.(+).)", "(.(+)..)", "(..(+).)", "(..(+)..)",
                      "(.(+)....)", "(...(+)....)")
        for structure, temperature in itertools.product(structures, (25, 37, 55)):
            spec, _ = mixed_fixture(structure)
            loop = next(face for face in decompose_structure(spec, structure) if face.kind == "interior")
            left, right = (b-a for a, b in loop.unpaired_runs)
            reference = []
            for material in ("rna", "dna"):
                p = self.native_parameters(material, temperature)
                # GC outer, reverse CG inner, four neighboring A bases.
                reference.append(self.native_energy(RNA.exp_E_IntLoop(left, right, 2, 1, 1, 1, 1, 1, p), p))
            observed = self.model.loop_energy(spec, Conditions(celsius=temperature, sodium=1), loop)
            self.assertTrue(math.isfinite(observed))
            self.assertAlmostEqual(observed, math.fsum(reference)/2, delta=1e-12)

    def test_small_interior_projection_tables_have_no_missing_wc_entries(self):
        for material, temperature in itertools.product(("rna", "dna"), (25, 37, 55)):
            p = self.native_parameters(material, temperature)
            for left, right in ((1, 1), (1, 2), (2, 1), (2, 2)):
                for outer, inner in itertools.product((1, 2, 5, 6), repeat=2):
                    for a, b, c, d in itertools.product(range(1, 5), repeat=4):
                        if (left == 1 and a != c) or (right == 1 and b != d):
                            continue  # The same unpaired base occurs at both ends.
                        weight = RNA.exp_E_IntLoop(left, right, outer, inner, a, b, c, d, p)
                        self.assertTrue(math.isfinite(weight) and weight > 0,
                            (material, temperature, left, right, outer, inner, a, b, c, d))

    @unittest.skipUnless(R_DNA_NN1 is not None, "optional independent Biopython table unavailable")
    def test_one_base_bulge_uses_hybrid_stack_instead_of_mean_homogeneous_stacks(self):
        for structure, temperature in itertools.product(("(.(+))", "((+).)", "(..(+))"), (25, 37, 55)):
            spec, _ = mixed_fixture(structure)
            loop = next(face for face in decompose_structure(spec, structure) if face.kind == "bulge")
            reference = []
            for material in ("rna", "dna"):
                p = self.native_parameters(material, temperature)
                reference.append(self.native_energy(RNA.doubleP_getitem(p.expbulge, loop.unpaired_count), p))
            expected = math.fsum(reference)/2
            c = Conditions(celsius=temperature, sodium=1)
            if loop.unpaired_count == 1:
                h, s = R_DNA_NN1["GG/CC"]
                expected += h-c.kelvin*s/1000
            self.assertAlmostEqual(self.model.loop_energy(spec, c, loop), expected, delta=1e-12)

    def test_three_to_five_strand_junctions_have_independent_composition_and_loop_sums(self):
        for count, temperature in itertools.product((3, 4, 5), (25, 37, 55)):
            # Each strand contributes two free bases to one multibranch face.
            strands = tuple(Strand("GAAC", "rna" if i % 2 == 0 else "dna", str(i)) for i in range(count))
            spec, partners = ComplexSpec(strands), list(range(4*count))
            for i, j in [(0, 4*count-1), *((4*k-1, 4*k) for k in range(1, count))]:
                partners[i], partners[j] = j, i
            c = Conditions(celsius=temperature, sodium=1)
            p = {m: self.native_parameters(m, temperature) for m in ("rna", "dna")}
            rna_fraction = sum(s.material == "rna" for s in strands)/count
            hybrid_init = 1.9 + c.kelvin*.0039
            major_init = self.native_energy(p["rna"].expDuplexInit, p["rna"])
            association = (count-1)*(hybrid_init+(2*rna_fraction-1)*(major_init-hybrid_init))
            loop = math.fsum(self.native_energy(v.expMLclosing, v)
                + count*self.native_energy(RNA.doubleP_getitem(v.expMLintern, 1), v)
                + 2*count*self.native_energy(v.expMLbase, v) for v in p.values())/2
            symmetry = p["rna"].kT/1000 * math.log(spec.symmetry_factor)
            self.assertAlmostEqual(self.model.energy(spec, c, PairingConstraints(tuple(partners))),
                                   association+loop+symmetry, delta=1e-11)

    def test_mixed_junction_keeps_actual_rna_terminal_cost_and_its_hairpin(self):
        spec = ComplexSpec((Strand("GAAAAUG", "rna", "r"), Strand("CC", "dna", "d")))
        faces = decompose_structure(spec, "((...)(+))")
        mixed = next(face for face in faces if face.kind == "multiloop")
        hairpin = next(face for face in faces if face.kind == "hairpin")
        for temperature in (25, 37, 55):
            c = Conditions(celsius=temperature, sodium=1)
            p = {m: self.native_parameters(m, temperature) for m in ("rna", "dna")}
            expected = math.fsum(self.native_energy(v.expMLclosing, v)
                + 3*self.native_energy(RNA.doubleP_getitem(v.expMLintern, 1), v) for v in p.values())/2
            expected += self.native_energy(p["rna"].expTermAU, p["rna"])
            self.assertAlmostEqual(self.model.loop_energy(spec, c, mixed), expected, delta=1e-12)
            expected_hairpin = self.native_energy(RNA.exp_E_Hairpin(3, 5, 1, 1, "AAAAU", p["rna"]), p["rna"])
            self.assertAlmostEqual(self.model.loop_energy(spec, c, hairpin), expected_hairpin, delta=1e-12)

    def test_association_interpolates_by_nucleotide_composition_in_both_directions(self):
        for rna_count, dna_count in ((0, 10), (1, 9), (5, 5), (9, 1), (10, 0)):
            strands = []
            if rna_count:
                strands.append(Strand("A"*rna_count, "rna", "r"))
            if dna_count:
                strands.append(Strand("A"*dna_count, "dna", "d"))
            if len(strands) == 1:  # homogeneous association endpoint, still L=2
                only = strands[0]
                strands = [Strand("AAAAA", only.material, "a"), Strand("AAAAA", only.material, "b")]
            spec = ComplexSpec(tuple(strands))
            for temperature in (25, 37, 55):
                c = Conditions(celsius=temperature, sodium=1)
                material = "rna" if rna_count >= dna_count else "dna"
                p = self.native_parameters(material, temperature)
                homogeneous = self.native_energy(p.expDuplexInit, p)
                hybrid = 1.9+c.kelvin*.0039
                expected = hybrid + abs(rna_count-dna_count)/10*(homogeneous-hybrid)
                self.assertAlmostEqual(self.model.association_energy(spec, c), expected, delta=1e-12)
        self.assertIn("native-sha256=", self.model.provenance.version)

    def test_malformed_local_faces_cannot_access_native_arrays(self):
        spec, _ = mixed_fixture("(.(+))")
        face = next(loop for loop in decompose_structure(spec, "(.(+))") if loop.kind == "bulge")
        for invalid in (replace(face, unpaired_runs=((-1, 0),)),
                        replace(face, unpaired_runs=((1, 1),)),
                        replace(face, child_pairs=()), replace(face, kind="stack")):
            with self.assertRaises(ValueError):
                self.model.loop_energy(spec, Conditions(sodium=1), invalid)

    def test_cyclic_rotation_preserves_all_small_mixed_state_energies(self):
        for count in (2, 3, 4, 5):
            spec = ComplexSpec(tuple(Strand("GC", "rna" if i % 2 == 0 else "dna", str(i)) for i in range(count)))
            rotated = ComplexSpec(spec.strands[1:]+spec.strands[:1])
            for state in enumerate_structures(spec, pair_allowed=self.model.pair_allowed, min_hairpin_size=3):
                partners = [None] * spec.n_nucleotides
                for i, j in enumerate(state.partners):
                    partners[(i-2) % len(partners)] = (j-2) % len(partners)
                c = Conditions(sodium=1)
                self.assertAlmostEqual(self.model.energy(spec, c, state),
                    self.model.energy(rotated, c, PairingConstraints(tuple(partners))), delta=1e-11)

    def test_exhaustive_ensemble_observables_share_the_same_energy_model(self):
        exhaustive = ExhaustiveBackend(pair_allowed=self.model.pair_allowed, energy=self.model.energy,
            thermal_energy=self.model.thermal_energy, provenance=self.model.provenance, min_hairpin_size=3)
        spec = ComplexSpec((Strand("GCGC", "rna", "r"), Strand("GCGC", "dna", "d")))
        result = exhaustive.evaluate(spec, Conditions(sodium=1))
        self.assertGreater(result.n_structures, 1)
        for row in result.pairs().matrix:
            self.assertAlmostEqual(math.fsum(row), 1.0, delta=1e-12)
        self.assertEqual(result.mfe()[0], min(result.energies_kcal_mol))
        self.assertEqual(result.sample(20, seed=2), result.sample(20, seed=2))
        self.assertEqual(result.partition.provenance, self.model.provenance)

    def test_published_sodium_option_is_explicit_versioned_and_bounded(self):
        corrected = HybridMotifModel(sodium_correction="published-pair")
        self.assertEqual(self.model.sodium_correction, "reference")
        self.assertNotEqual(corrected.provenance.parameter_digest, self.model.provenance.parameter_digest)
        self.assertTrue(corrected.provenance.version.startswith("2;"))
        with self.assertRaises(ValueError):
            HybridMotifModel(sodium_correction="automatic")
        for sodium in (.11999, 1.00001):
            with self.assertRaises(ValueError):
                corrected.pair_salt_energy(Conditions(sodium=sodium))
        for temperature in (25, 37, 55):
            c = Conditions(celsius=temperature, sodium=.12)
            entropy = 114/310.15 * math.log(.12)
            self.assertAlmostEqual(corrected.pair_salt_energy(c), -c.kelvin*entropy/1000, delta=1e-15)
            self.assertEqual(corrected.pair_salt_energy(Conditions(celsius=temperature, sodium=1)), 0.0)

    def test_published_sodium_counts_each_pair_once_and_keeps_one_molar_identical(self):
        corrected = HybridMotifModel(sodium_correction="published-pair")
        for spec, state in (mixed_fixture("(.(+).)"), mixed_fixture("((...)(+))"),
                            mixed_fixture("(...)."), mixed_fixture("((+))")):
            lengths = tuple(len(s.sequence) for s in spec.strands)
            faces = decompose_structure(spec, format_structure(state, lengths))
            pair_count = sum(j > i for i, j in enumerate(state.partners))
            self.assertEqual(sum(face.closing_pair is not None for face in faces), pair_count)
            for temperature in (25, 37, 55):
                reference = Conditions(celsius=temperature, sodium=1)
                baseline = self.model.energy(spec, reference, state)
                self.assertEqual(corrected.energy(spec, reference, state), baseline)
                for sodium in (.12, .5):
                    c = Conditions(celsius=temperature, sodium=sodium)
                    delta = pair_count * -.114*math.log(sodium)*c.kelvin/310.15
                    self.assertAlmostEqual(corrected.energy(spec, c, state)-baseline, delta, delta=1e-12)
                    self.assertEqual(corrected.association_energy(spec, c), self.model.association_energy(spec, reference))

    def test_published_sodium_preserves_homogeneous_motif_and_reaction_pair_balance(self):
        corrected = HybridMotifModel(sodium_correction="published-pair")
        rna = Strand("GAAACG", "rna", "r")
        dna = Strand("C", "dna", "d")
        mono, trigger, mixed = ComplexSpec((rna,)), ComplexSpec((dna,)), ComplexSpec((rna, dna))
        mono_hairpin = next(face for face in decompose_structure(mono, "(...).") if face.kind == "hairpin")
        mixed_hairpin = next(face for face in decompose_structure(mixed, "(...)(+)") if face.kind == "hairpin")
        for sodium in (.12, .5, 1):
            c, reference = Conditions(sodium=sodium), Conditions(sodium=1)
            self.assertEqual(corrected.loop_energy(mono, c, mono_hairpin), corrected.loop_energy(mixed, c, mixed_hairpin))
            salted_reaction = corrected.energy(mixed, c, "(...)(+)")-corrected.energy(mono, c, "(...).")-corrected.energy(trigger, c, ".")
            reference_reaction = self.model.energy(mixed, reference, "(...)(+)")-self.model.energy(mono, reference, "(...).")-self.model.energy(trigger, reference, ".")
            self.assertAlmostEqual(salted_reaction-reference_reaction, corrected.pair_salt_energy(c), delta=1e-12)

    def test_published_sodium_reweighting_reconstructs_pf_pairs_and_monotonicity(self):
        corrected = HybridMotifModel(sodium_correction="published-pair")
        def evaluator(model):
            return ExhaustiveBackend(pair_allowed=model.pair_allowed, energy=model.energy,
                thermal_energy=model.thermal_energy, provenance=model.provenance, min_hairpin_size=3)
        specs = (ComplexSpec((Strand("GCGAAAGC", "rna", "r"),)),
                 ComplexSpec((Strand("GCGAAAGC", "dna", "d"),)),
                 ComplexSpec((Strand("GCGC", "rna", "r"), Strand("GCGC", "dna", "d"))))
        for spec, temperature in itertools.product(specs, (25, 37, 55)):
            baseline = evaluator(self.model).evaluate(spec, Conditions(celsius=temperature, sodium=1))
            counts = [sum(j > i for i, j in enumerate(state.partners)) for state in baseline.structures]
            thermal = baseline.partition.thermal_energy_kcal_mol
            previous_pairs = -1
            for sodium in (.12, .3, .6, 1):
                c = Conditions(celsius=temperature, sodium=sodium)
                observed = evaluator(corrected).evaluate(spec, c)
                factor = -.114*math.log(sodium)*c.kelvin/310.15
                weights = [probability*math.exp(-n*factor/thermal) for probability, n in zip(baseline.probabilities, counts)]
                normalization = math.fsum(weights)
                self.assertEqual(observed.structures, baseline.structures)
                self.assertAlmostEqual(observed.partition.free_energy_kcal_mol,
                    baseline.partition.free_energy_kcal_mol-thermal*math.log(normalization), delta=1e-12)
                for expected, actual in zip(weights, observed.probabilities):
                    self.assertAlmostEqual(actual, expected/normalization, delta=1e-12)
                expected_pairs = math.fsum(n*w for n, w in zip(counts, weights))/normalization
                observed_pairs = math.fsum(1-row[i] for i, row in enumerate(observed.pairs().matrix))/2
                self.assertAlmostEqual(observed_pairs, expected_pairs, delta=1e-12)
                self.assertGreaterEqual(observed_pairs+1e-12, previous_pairs)
                previous_pairs = observed_pairs


if __name__ == "__main__":
    unittest.main()
