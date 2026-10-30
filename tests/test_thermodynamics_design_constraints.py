"""Sequence admissibility and architecture parity without native engines."""
from dataclasses import replace
import random
import unittest

from aptaswitch_core.architecture import architecture_layout, default_architecture_for_trigger, target_switch_structures
from aptaswitch_core.models import InverseExtensionSettings, ReporterContext, RunRequest, SequenceInput, ThermoConditions
from aptaswitch_core.preview import build_switch_model
from aptaswitch_core.sequences import reverse_complement_rna
from aptaswitch_core.thermodynamics.design_constraints import (
    BaseSite, DesignAssignment, DesignProblem, LibraryConstraint, NoMutableSitesError,
    StrandPattern, UnsatisfiableDesignError, build_design_problem, design_bases,
)
from aptaswitch_core.thermodynamics.types import Strand


def request(length=30, **changes):
    return RunRequest(SequenceInput("ACGTACGT", ("ACGT" * 20)[:length]),
                      default_architecture_for_trigger(length), trials=1, **changes)


def changed_sites(before, after):
    return {BaseSite(a.name, i) for a, b in zip(before.strands, after.strands)
            for i, (x, y) in enumerate(zip(a.sequence, b.sequence)) if x != y}


class GenericConstraintTests(unittest.TestCase):
    def test_design_alphabets_are_explicit_and_chemistry_aware(self):
        expected = {"A": "A", "C": "C", "G": "G", "T": "T", "R": "AG", "Y": "CT",
                    "S": "CG", "W": "AT", "K": "GT", "M": "AC", "B": "CGT",
                    "D": "AGT", "H": "ACT", "V": "ACG", "X": "ACGT"}
        for code, bases in expected.items():
            self.assertEqual(design_bases(code, "dna"), tuple(bases))
            self.assertEqual(design_bases(code.replace("T", "U"), "rna"), tuple(bases.replace("T", "U")))
        for code, chemistry in (("T", "rna"), ("U", "dna"), ("N", "dna"), ("XX", "rna"), ("X", "mixed")):
            with self.assertRaises(ValueError):
                design_bases(code, chemistry)

    def test_library_words_respect_the_declared_site_order(self):
        sites = (BaseSite("A", 1), BaseSite("B", 0), BaseSite("A", 0))
        problem = DesignProblem((StrandPattern("A", "rna", "XX"), StrandPattern("B", "dna", "X")),
                                (LibraryConstraint("words", sites, ("UTA", "CGG")),))
        observed = set()
        for seed in range(20):
            assignment = problem.generate(random.Random(seed))
            problem.validate(assignment)
            observed.add((assignment.sequence("A"), assignment.sequence("B")))
        self.assertEqual(observed, {("AU", "T"), ("GC", "G")})
        mutated = problem.mutate(assignment, random.Random(3), site=sites[0])
        self.assertEqual(len(changed_sites(assignment, mutated)), 3)

    def test_fixed_patterns_and_effectively_fixed_couplings_do_not_loop(self):
        fixed = DesignProblem((StrandPattern("fixed", "dna", "ACGT"),))
        value = fixed.generate(random.Random(1))
        self.assertEqual(value.sequence("fixed"), "ACGT")
        with self.assertRaises(NoMutableSitesError):
            fixed.mutate(value, random.Random(1))
        coupled = DesignProblem((StrandPattern("A", "rna", "AX"),),
            (LibraryConstraint("paired", (BaseSite("A", 0), BaseSite("A", 1)), ("AU", "UA", "CG", "GC")),))
        value = coupled.generate(random.Random(1))
        self.assertEqual(value.sequence("A"), "AU")
        with self.assertRaises(NoMutableSitesError):
            coupled.mutate(value, random.Random(1), site=BaseSite("A", 1))

    def test_incompatible_fixed_bases_and_cyclic_constraints_fail_explicitly(self):
        problem = DesignProblem((StrandPattern("A", "rna", "AA"),),
            (LibraryConstraint("paired", (BaseSite("A", 0), BaseSite("A", 1)), ("AU", "UA")),))
        with self.assertRaises(UnsatisfiableDesignError):
            problem.generate(random.Random(1))
        # Arc consistency alone cannot solve this odd cycle of binary complements.
        triangle = DesignProblem((StrandPattern("A", "rna", "WWW"),), tuple(
            LibraryConstraint(f"pair{i}", (BaseSite("A", i), BaseSite("A", (i + 1) % 3)), ("AU", "UA"))
            for i in range(3)))
        with self.assertRaises(UnsatisfiableDesignError):
            triangle.generate(random.Random(1))

    def test_strict_validation_rejects_unknown_sites_chemistry_and_rng(self):
        pattern = StrandPattern("A", "rna", "X")
        for sites, words in (((BaseSite("A", 1),), ("A",)), ((BaseSite("A", 0),), ("T",))):
            with self.assertRaises(ValueError):
                DesignProblem((pattern,), (LibraryConstraint("bad", sites, words),))
        problem = DesignProblem((pattern,))
        with self.assertRaises(ValueError):
            problem.generate(123)
        with self.assertRaises(ValueError):
            problem.validate(DesignAssignment((Strand("A", "dna", "A"),)))
        with self.assertRaises(ValueError):
            problem.mutate(problem.generate(random.Random(1)), random.Random(1), site=BaseSite("B", 0))


class CurrentSwitchDesignTests(unittest.TestCase):
    def assert_request_parity(self, req):
        normalized = req.normalized()
        before = normalized.to_dict()
        problem = build_design_problem(req)
        preview = build_switch_model(normalized.architecture, normalized.design_trigger_pattern(),
                                     reporter_sequence=normalized.reporter.sequence_after_start_rna)
        expected = {"Trigger": preview.trigger_sequence, "Switch": preview.switch_sequence}
        for pattern in problem.strands:
            self.assertEqual(pattern.pattern, expected[pattern.name])
        self.assertEqual(dict(problem.layout), architecture_layout(req.architecture))
        self.assertEqual(tuple(t.structure for t in problem.targets[:2]),
                         target_switch_structures(req.architecture, len(req.reporter.sequence_after_start_rna)))
        trigger = normalized.thermo.resolved_trigger_concentration_m
        switch = normalized.thermo.resolved_switch_concentration_m
        duplex = min(trigger, switch)
        self.assertEqual(problem.tubes[0].name, "tube_off")
        self.assertEqual(problem.tubes[0].on_targets, (("switch_off", switch),))
        expected_on = {"switch_on": duplex}
        if trigger > duplex:
            expected_on["trigger_free"] = trigger - duplex
        if switch > duplex:
            expected_on["switch_off"] = switch - duplex
        self.assertEqual(problem.tubes[1].name, "tube_on")
        self.assertEqual(dict(problem.tubes[1].on_targets), expected_on)
        for tube in problem.tubes:
            self.assertEqual(tube.max_complex_size, normalized.thermo.screening_max_complex_size)
        self.assertEqual(normalized.to_dict(), before)
        return problem

    def test_classic_24_and_30_match_current_patterns_offsets_and_tubes(self):
        for length in (24, 30):
            with self.subTest(length=length):
                problem = self.assert_request_parity(request(length))
                value = problem.generate(random.Random(2026))
                self.assertEqual(value.sequence("Trigger"), request(length).sequences.effective_trigger())
                switch = value.sequence("Switch")
                start = problem.layout["len_t7"]
                self.assertEqual(switch[start:start + length], reverse_complement_rna(value.sequence("Trigger")))
                self.assertEqual(switch[problem.layout["aug_start"]:problem.layout["aug_start"] + 3], "AUG")
                self.assertEqual(switch[problem.layout["reporter_start"]:], request(length).reporter.sequence_after_start_rna)
                self.assertEqual(problem.reporter, request(length).reporter)
                with self.assertRaises(TypeError):
                    problem.layout["aug_start"] = 0

    def test_inverse_both_ends_and_custom_layout_preserve_original_core(self):
        custom = replace(default_architecture_for_trigger(30), toehold_length=9, stem_length=8,
                         upper_stem2_length=2, bulge_length=2, upper_stem_length=7,
                         rbs_prefix_length=0, t7_leader="", rbs_loop="AGGAGG", aug_spacer=2,
                         frame_linker="ACUXX")
        for arch in (default_architecture_for_trigger(24), default_architecture_for_trigger(30), custom):
            target = arch.trigger_region_length
            core = ("ACGT" * 10)[:target - 5]
            for side in ("5prime", "3prime", "both"):
                with self.subTest(target=target, side=side):
                    req = replace(request(target), sequences=SequenceInput("ACGTACGT", core), architecture=arch,
                                  inverse_extension=InverseExtensionSettings(target, side))
                    problem = self.assert_request_parity(req)
                    prefix, _ = req.inverse_extension.extension_lengths(len(core))
                    rng = random.Random(51)
                    value = problem.generate(rng)
                    for _ in range(30):
                        problem.validate(value)
                        trigger = value.sequence("Trigger")
                        self.assertEqual(trigger[prefix:prefix + len(core)], core)
                        start = problem.layout["len_t7"]
                        self.assertEqual(value.sequence("Switch")[start:start + target], reverse_complement_rna(trigger))
                        value = problem.mutate(value, rng)

    def test_specific_inverse_mutation_repairs_wc_and_off_wobble_together(self):
        req = replace(request(30), sequences=SequenceInput("ACGT", "ACGT"),
                      inverse_extension=InverseExtensionSettings(30, "5prime"))
        problem = build_design_problem(req)
        rng = random.Random(81)
        value = problem.generate(rng)
        # 5' trigger base zero pairs with the upper-left arm, itself paired in OFF.
        site = BaseSite("Trigger", 0)
        mutation = problem.mutate(value, rng, site=site)
        changes = changed_sites(value, mutation)
        self.assertIn(site, changes)
        self.assertIn(BaseSite("Switch", problem.layout["len_t7"] + 29), changes)
        component = next(c for c in problem._components if site in c.sites)
        self.assertEqual(len(component.sites), 3)
        self.assertLessEqual(changes, set(component.sites))
        problem.validate(mutation)
        # Fix trigger C -> RNA left G, allowing RNA right C *or* U.
        toy = DesignProblem((StrandPattern("D", "dna", "C"), StrandPattern("R", "rna", "XX")),
            (LibraryConstraint("on", (BaseSite("D", 0), BaseSite("R", 0)), ("AU", "CG", "GC", "TA")),
             LibraryConstraint("off", (BaseSite("R", 0), BaseSite("R", 1)), ("AU", "UA", "CG", "GC", "GU", "UG"))))
        value = toy.generate(random.Random(10))
        mutation = toy.mutate(value, random.Random(10), site=BaseSite("R", 1))
        self.assertEqual({value.sequence("R"), mutation.sequence("R")}, {"GC", "GU"})

    def test_seeded_generation_and_mutations_are_reproducible(self):
        req = replace(request(30), sequences=SequenceInput("ACGT", "ACGT"),
                      inverse_extension=InverseExtensionSettings(30, "both"))
        problem = build_design_problem(req)
        def trajectory(seed):
            rng = random.Random(seed)
            value = problem.generate(rng)
            values = [value]
            for _ in range(50):
                value = problem.mutate(value, rng)
                values.append(value)
            return values
        self.assertEqual(trajectory(7), trajectory(7))
        self.assertNotEqual(trajectory(7), trajectory(8))

    def test_provided_sequences_and_reporter_prefix_never_change(self):
        # Concrete sequence records as supplied by FASTA import; the full
        # reporter remains stored separately and is not extended into the design.
        trigger = "GGGCACGTATACGGACTGCACTGGGACATGC"[:30]
        reporter = ReporterContext("fixed", "Provided reporter", "GCT" * 10, "GCT" * 100)
        req = replace(request(30), sequences=SequenceInput("ACGT", trigger), reporter=reporter)
        problem = build_design_problem(req)
        rng = random.Random(5)
        value = problem.generate(rng)
        for _ in range(100):
            self.assertEqual(value.sequence("Trigger"), trigger)
            self.assertEqual(value.sequence("Switch")[problem.layout["reporter_start"]:], "GCU" * 10)
            for region in problem.regions:
                if region.name in ("rbs", "t7_leader", "aug", "frame_linker"):
                    original = next(s.pattern for s in problem.strands if s.name == "Switch")[region.start:region.stop]
                    self.assertEqual(value.sequence("Switch")[region.start:region.stop], original)
            value = problem.mutate(value, rng)
        self.assertEqual(problem.reporter.full_sequence_after_start_rna, "GCT" * 100)
        altered = DesignAssignment((Strand("A" + trigger[1:], "dna", "Trigger"), value.strands[1]))
        with self.assertRaisesRegex(ValueError, "fixed base"):
            problem.validate(altered)

    def test_unavoidable_stop_fails_before_generation_and_opt_out_is_preserved(self):
        req = replace(request(30), sequences=SequenceInput("ACGT", "G" * 30),
                      architecture=replace(default_architecture_for_trigger(30), frame_linker="UAA"))
        with self.assertRaisesRegex(UnsatisfiableDesignError, "STOP"):
            build_design_problem(req)
        unfiltered = build_design_problem(replace(req, exclude_stop_candidates=False, exclude_non_linear_candidates=False))
        value = unfiltered.generate(random.Random(1))
        self.assertTrue(unfiltered.reading_frame_report(value)["has_premature_stop"])
        self.assertTrue(unfiltered.passes_sequence_screen(value))
        self.assertFalse(unfiltered.screening.exclude_non_linear_candidates)
        # A STOP at the reporter boundary is not a premature pre-reporter STOP.
        req = replace(req, architecture=replace(req.architecture, frame_linker=""),
                      reporter=ReporterContext("stop", "Boundary", "UAAGCU"))
        problem = build_design_problem(req)
        self.assertTrue(problem.passes_sequence_screen(problem.generate(random.Random(1))))

    def test_asymmetric_target_concentrations_keep_free_strand_targets(self):
        for trigger, switch, excess in ((2e-6, 1e-6, "trigger_free"), (1e-6, 2e-6, "switch_off")):
            thermo = replace(ThermoConditions(), trigger_concentration_m=trigger, switch_concentration_m=switch,
                             screening_max_complex_size=3)
            problem = self.assert_request_parity(request(30, thermo=thermo))
            self.assertEqual(dict(problem.tubes[1].on_targets), {"switch_on": 1e-6, excess: 1e-6})
            self.assertEqual(problem.tubes[0].max_complex_size, 3)

    def test_new_designs_default_to_three_strand_validation_and_reject_larger_catalogs(self):
        normalized = request(30).normalized()
        self.assertEqual(normalized.thermo.validation_max_complex_size, 3)
        for key in ('screening_max_complex_size', 'validation_max_complex_size'):
            for size in (4, 5, 6):
                historical = replace(ThermoConditions(), **{key: size})
                self.assertEqual(getattr(historical, key), size)
                with self.subTest(key=key, size=size), self.assertRaisesRegex(ValueError, 'two to three'):
                    request(30, thermo=historical).normalized()

    def test_switch_rna_trigger_dna_and_pairing_rules_are_fixed(self):
        problem = self.assert_request_parity(request(30))
        value = problem.generate(random.Random(8))
        self.assertEqual(tuple(s.material for s in value.strands), ("dna", "rna"))
        for library in problem.libraries:
            if library.name.startswith("switch_off"):
                self.assertIn("GU", library.allowed_words)
            else:
                self.assertNotIn("GU", library.allowed_words)
                self.assertNotIn("UG", library.allowed_words)

    def test_other_design_chemistries_are_rejected_without_conversion(self):
        for switch, trigger in (("dna", "rna"), ("rna", "rna"), ("dna", "dna")):
            req = request(30, thermo=replace(ThermoConditions(), switch_material=switch, trigger_material=trigger))
            with self.subTest(switch=switch, trigger=trigger), self.assertRaises(ValueError):
                build_design_problem(req)

    def test_free_x_patterns_resolve_to_concrete_candidates_and_n_is_rejected(self):
        req = replace(request(30), architecture=replace(request(30).architecture, frame_linker="CXX"))
        problem = build_design_problem(req)
        self.assertTrue(any("X" in strand.pattern for strand in problem.strands))
        self.assertTrue(all("N" not in strand.pattern for strand in problem.strands))
        for strand in problem.generate(random.Random(4)).strands:
            self.assertLessEqual(set(strand.sequence), set("ACGU" if strand.material == "rna" else "ACGT"))
        with self.assertRaises(ValueError):
            replace(req, architecture=replace(req.architecture, frame_linker="CNN")).normalized()
        for ambiguous in ("ACNX", "ACGN", "ACGX"):
            with self.subTest(ambiguous=ambiguous), self.assertRaises(ValueError):
                replace(req, sequences=replace(req.sequences, trigger_dna=ambiguous)).normalized()

    def test_full_x_linker_can_be_optimized_without_changing_fixed_bases_or_frame(self):
        fixed = request(30)
        req = replace(fixed, sequences=SequenceInput("ACGT", "G" * 30),
                      architecture=replace(fixed.architecture,
                          frame_linker="X" * len(fixed.architecture.resolved_frame_linker())))
        problem = build_design_problem(req)
        linker = next(region for region in problem.regions if region.name == "frame_linker")
        for offset in ("aug_start", "reporter_start"):
            self.assertEqual(problem.layout[offset], architecture_layout(fixed.architecture)[offset])
        rng = random.Random(17)
        value = problem.generate(rng)
        for index in range(linker.start, linker.stop):
            site = BaseSite("Switch", index)
            mutation = problem.mutate_guided(value, rng, {site: 1.0},
                                             exploration_fraction=0, multi_fraction=0)
            component = next(c for c in problem._components if site in c.sites)
            self.assertTrue(changed_sites(value, mutation))
            self.assertLessEqual(changed_sites(value, mutation), set(component.sites))
            self.assertFalse(problem.reading_frame_report(mutation)["has_premature_stop"])
            for strand, pattern in zip(mutation.strands, problem.strands):
                alphabet = set("ACGU" if strand.material == "rna" else "ACGT")
                self.assertLessEqual(set(strand.sequence), alphabet)
                for actual, expected in zip(strand.sequence, pattern.pattern):
                    if expected != "X":
                        self.assertEqual(actual, expected)
            value = mutation
        # Disabling the constraint retains the exact reading-frame diagnostic.
        unfiltered = build_design_problem(replace(req, exclude_stop_candidates=False))
        for prefix, has_stop in (("CCC", False), ("UAA", True)):
            switch = value.sequence("Switch")
            switch = (switch[:linker.start] + prefix + "C" * (linker.stop-linker.start-3)
                      + switch[linker.stop:])
            candidate = DesignAssignment((value.strands[0], Strand(switch, "rna", "Switch")))
            self.assertEqual(unfiltered.reading_frame_report(candidate)["has_premature_stop"], has_stop)
            if has_stop:
                with self.assertRaisesRegex(ValueError, "coupled sequence library"):
                    problem.validate(candidate)
            else:
                self.assertTrue(problem.passes_sequence_screen(candidate))

    def test_inverse_extension_requires_dna_for_both_initial_strands(self):
        for aptamer, trigger in (("dna", "rna"), ("rna", "rna"), ("rna", "dna")):
            settings = InverseExtensionSettings(30, thermo=replace(ThermoConditions(),
                                               switch_material=aptamer, trigger_material=trigger))
            req = replace(request(24), architecture=default_architecture_for_trigger(30), inverse_extension=settings)
            with self.subTest(aptamer=aptamer, trigger=trigger), self.assertRaises(ValueError):
                req.normalized()


if __name__ == "__main__":
    unittest.main()
