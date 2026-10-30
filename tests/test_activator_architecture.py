"""The activator defaults must reach actual design targets and frame scanning."""

from __future__ import annotations

import random
import unittest
from dataclasses import replace
from unittest.mock import Mock

from aptaswitch_core.architecture import architecture_layout, default_architecture_for_trigger, target_switch_structures
from aptaswitch_core.models import DEFAULT_RBS_SEQUENCE, ReporterContext, RunRequest, SequenceInput, SwitchArchitecture
from aptaswitch_core.thermodynamics.design_constraints import StrandPattern, UnsatisfiableDesignError, build_design_problem
from aptaswitch_core.thermodynamics.design_engine import design_candidates, evaluate_candidate
from tests.test_thermodynamics_design_engine import AnalyticalBackend
from aptaswitch_core.preview import build_switch_model, candidate_display_structure, candidate_structure_svg
from aptaswitch_core.sequences import reading_frame_report, reverse_complement_rna


TRIGGER = "GTCCAGGCTGGTATAATTAGATCCACGTAC"
CLASSIC_LINKER = "AACCUGGCGGCAGCGCAAAAG"


def configured_architecture(trigger_length, **changes):
    """Resize the full left arm's lower segment for custom spacing fixtures."""
    architecture = replace(default_architecture_for_trigger(trigger_length), **changes)
    if architecture.binding_layout == "full_left_arm" and "stem_length" not in changes:
        architecture = replace(
            architecture,
            stem_length=architecture.stem_length + trigger_length - architecture.trigger_region_length,
        )
    return architecture


def make_request(**architecture_changes):
    return RunRequest(
        sequences=SequenceInput(aptamer_dna="ACGTACGTACGT", trigger_dna=TRIGGER),
        architecture=configured_architecture(len(TRIGGER), **architecture_changes),
        trials=1,
    ).normalized()


def make_legacy_request(**architecture_changes):
    architecture_changes.setdefault("frame_linker", None)
    return make_request(binding_layout="legacy_lower_stem", stem_length=18, **architecture_changes)




def switch_domains(problem):
    """Read the actual native strand pattern using its declared region offsets."""
    switch = next(strand for strand in problem.strands if strand.name == "Switch")
    return [StrandPattern(region.name, switch.material, switch.pattern[region.start:region.stop])
            for region in problem.regions if region.strand == "Switch"]


def switch_targets(problem):
    return {target.name: target for target in problem.targets}


def resolved_test_sequence(problem):
    """Resolve free symbols solely for layout checks, not predicted physics."""
    return "".join(domain.pattern.replace("X", "A") for domain in switch_domains(problem))


class ActivatorArchitectureTests(unittest.TestCase):
    def test_standard_activator_matches_trigger_length_with_six_bases_to_aug(self):
        request = make_request()
        arch = request.architecture
        for trigger_length, split in ((24, (15, 9)), (30, (12, 9)), (32, (12, 20))):
            with self.subTest(trigger_length=trigger_length):
                standard = default_architecture_for_trigger(trigger_length)
                self.assertEqual((standard.toehold_length, standard.stem_length), split)
                self.assertEqual(standard.trigger_region_length, trigger_length)
                self.assertEqual(standard.upper_stem2_length, 0 if trigger_length == 30 else 4)
                self.assertEqual(standard.binding_layout, "full_left_arm" if trigger_length == 30 else "legacy_lower_stem")
        self.assertEqual((arch.rbs_prefix_length, arch.rbs_to_aug_distance), (3, 6))
        self.assertEqual(arch.aug_spacer, 0)
        layout = architecture_layout(arch)
        self.assertEqual(layout["aug_start"] - layout["rbs_start"] - len(arch.rbs_loop), 6)
        self.assertEqual(request.to_dict()["architecture"]["rbs_prefix_length"], 3)

    def test_standard_30nt_places_exact_linker_after_nine_stem_bases_before_reporter(self):
        request = make_request()
        arch = request.architecture
        self.assertEqual(arch.rbs_loop, "AGAGGAGA")
        self.assertEqual(arch.rbs_loop, DEFAULT_RBS_SEQUENCE)
        self.assertEqual(arch.upper_stem2_length, 0)
        self.assertEqual(arch.frame_linker, CLASSIC_LINKER)
        self.assertEqual(arch.resolved_frame_linker(), CLASSIC_LINKER)
        self.assertEqual(len(CLASSIC_LINKER), 21)
        design = build_design_problem(request)
        domains = switch_domains(design)
        names = [domain.name for domain in domains]
        self.assertNotIn("upper2_l", names)
        self.assertNotIn("upper2_r", names)
        self.assertEqual(names[names.index("aug") + 1:], ["stem_r", "frame_linker", "reporter_after_start"])
        self.assertEqual(next(domain.pattern for domain in domains if domain.name == "stem_r"), "X" * 9)
        self.assertEqual(next(domain.pattern for domain in domains if domain.name == "rbs"), "AGAGGAGA")
        self.assertEqual(next(domain.pattern for domain in domains if domain.name == "frame_linker"), CLASSIC_LINKER)
        layout = design.layout
        sequence = resolved_test_sequence(design)
        reporter_start = sum(len(domain.pattern) for domain in domains[:-1])
        self.assertEqual(sequence[layout["aug_start"]:layout["aug_start"] + 3], "AUG")
        self.assertEqual(reporter_start - (layout["aug_start"] + 3), 30)
        self.assertEqual(sequence[layout["aug_start"] + 3:reporter_start], "A" * 9 + CLASSIC_LINKER)
        self.assertEqual(sequence[reporter_start:], request.reporter.sequence_after_start_rna)
        self.assertEqual(layout["reporter_start"], 83)
        self.assertEqual((reporter_start - layout["aug_start"]) % 3, 0)
        self.assertEqual(request.to_dict()["architecture"]["frame_linker"], CLASSIC_LINKER)
        for stop_position, is_premature in ((reporter_start - 3, True), (reporter_start, False)):
            with self.subTest(stop_position=stop_position):
                with_stop = sequence[:stop_position] + "UAA" + sequence[stop_position + 3:]
                report = reading_frame_report(with_stop, layout["aug_start"], reporter_start=layout["reporter_start"])
                self.assertEqual(report["has_premature_stop"], is_premature)
        for target in (switch_targets(design)["switch_off"], switch_targets(design)["switch_on"]):
            self.assertEqual(len(target.structure.split("+")[-1]), len(sequence))
            self.assertEqual(target.structure.count("("), target.structure.count(")"))
            self.assertEqual(target.structure.split("+")[-1][62:83], "." * 21)

    def test_custom_frame_linker_reaches_native_as_two_free_unpaired_bases(self):
        request = make_legacy_request(upper_stem2_length=4)
        self.assertIsNone(request.architecture.frame_linker)
        self.assertEqual(request.architecture.resolved_frame_linker(), "XX")
        design = build_design_problem(request)
        domains = switch_domains(design)
        linker_index = next(i for i, domain in enumerate(domains) if domain.name == "frame_linker")
        self.assertEqual(domains[linker_index].pattern, "XX")
        self.assertEqual(domains[linker_index + 1].name, "reporter_after_start")
        linker_start = sum(len(domain.pattern) for domain in domains[:linker_index])
        for target in (switch_targets(design)["switch_off"], switch_targets(design)["switch_on"]):
            self.assertEqual(target.structure.split("+")[-1][linker_start:linker_start + 2], "..")

    def test_30nt_off_target_has_six_pairs_before_unpaired_aug_and_nine_after(self):
        request = make_request()
        design = build_design_problem(request)
        off = switch_targets(design)["switch_off"].structure
        sequence = resolved_test_sequence(design)
        # Independently parse the structure sent to native constraints, rather than
        # trusting the layout helper that also constructs the target.
        pairs = {}
        stack = []
        for index, token in enumerate(off):
            if token == "(":
                stack.append(index)
            elif token == ")":
                partner = stack.pop()
                pairs[index] = partner
                pairs[partner] = index
        self.assertEqual(stack, [])
        self.assertEqual(len(pairs) // 2, 15)
        self.assertEqual(sequence[50:53], "AUG")
        self.assertEqual(off[50:53], "...")
        self.assertEqual([pairs[index] for index in range(44, 50)], list(range(32, 26, -1)))
        self.assertEqual([pairs[index] for index in range(53, 62)], list(range(23, 14, -1)))
        self.assertEqual(off[24:27], "...")
        self.assertEqual(off[33:44], "." * 11)
        self.assertEqual(sequence[62:83], CLASSIC_LINKER)
        self.assertEqual(off[62:83], "." * 21)
        self.assertEqual(sequence[83:], request.reporter.sequence_after_start_rna)

    def test_30nt_recognizes_fixed_bulge_and_upper_left_bases_of_trigger(self):
        request = make_request()
        design = build_design_problem(request)
        domains = {domain.name: domain.pattern for domain in switch_domains(design)}
        binding = reverse_complement_rna(TRIGGER)
        self.assertEqual(domains["toehold"], binding[:12])
        self.assertEqual(domains["stem_l"], binding[12:21])
        self.assertEqual(domains["bulge_l"], binding[21:24])
        self.assertEqual(domains["upper_l"], binding[24:30])
        self.assertEqual(domains["stem_r"], "X" * 9)
        self.assertEqual(domains["upper_r"], "X" * 6)
        self.assertEqual(domains["rbs_prefix"], "X" * 3)
        on_trigger, on_switch = switch_targets(design)["switch_on"].structure.split("+")
        self.assertEqual(on_trigger, "(" * 30)
        self.assertEqual(on_switch[3:33], ")" * 30)
        self.assertTrue(all(token == "." for token in on_switch[33:]))
        self.assertEqual(request.to_dict()["architecture"]["binding_layout"], "full_left_arm")

    def test_missing_binding_layout_preserves_historical_30nt_offsets_and_domains(self):
        # This fixture represents a saved run from before the layout field;
        # it deliberately has the original 18-bp lower stem.
        saved = {
            "toehold_length": 12, "stem_length": 18,
            "upper_stem_length": 6, "upper_stem2_length": 0,
            "bulge_length": 3, "rbs_prefix_length": 3,
        }
        architecture = SwitchArchitecture(**saved)
        self.assertEqual(architecture.binding_layout, "legacy_lower_stem")
        architecture.validate_for_trigger(30)
        request = replace(make_request(), architecture=architecture).normalized()
        design = build_design_problem(request)
        domains = {domain.name: domain.pattern for domain in switch_domains(design)}
        self.assertEqual(domains["stem_l"], reverse_complement_rna(TRIGGER)[12:])
        self.assertEqual(domains["bulge_l"], "XXX")
        self.assertEqual(domains["upper_l"], "X" * 6)
        self.assertEqual(domains["stem_r"], "X" * 18)
        self.assertEqual(design.layout["aug_start"], 59)
        self.assertEqual(design.layout["reporter_start"], 80)
        self.assertEqual(design.layout["reporter_start"] - design.layout["aug_start"], 21)

    def test_saved_full_arm_with_none_linker_keeps_its_original_reporter_offset(self):
        # A saved beta design explicitly stored None before the classic linker
        # was added. Reloading it must not inject 21 bases into its sequence.
        architecture = SwitchArchitecture(
            toehold_length=12, stem_length=9, upper_stem2_length=0,
            binding_layout="full_left_arm", frame_linker=None,
        )
        request = replace(make_request(), architecture=architecture).normalized()
        design = build_design_problem(request)
        self.assertEqual(architecture.resolved_frame_linker(), "")
        self.assertNotIn("frame_linker", [domain.name for domain in switch_domains(design)])
        self.assertEqual(design.layout["reporter_start"], 62)
        self.assertEqual(design.layout["reporter_start"] - design.layout["aug_start"], 12)
        self.assertEqual(resolved_test_sequence(design)[62:], request.reporter.sequence_after_start_rna)

    def test_explicit_empty_linker_removes_classic_domain_without_fallback(self):
        request = make_request(frame_linker="")
        design = build_design_problem(request)
        self.assertEqual(request.architecture.frame_linker, "")
        self.assertNotIn("frame_linker", [domain.name for domain in switch_domains(design)])
        self.assertEqual(design.layout["reporter_start"], 62)
        self.assertNotIn(CLASSIC_LINKER, resolved_test_sequence(design))

    def test_24nt_target_is_unchanged(self):
        architecture = default_architecture_for_trigger(24)
        off, on = target_switch_structures(architecture, 3)
        expected = (
            "." * 18 + "(" * 13 + "." * 3 + "(" * 6
            + "." * 11 + ")" * 6 + "." * 3 + ")" * 13 + "." * 5
        )
        self.assertEqual(off, expected)
        self.assertEqual(on, "(" * 24 + "+" + "." * 3 + ")" * 24 + "." * (len(expected) - 27))

    def test_full_arm_dimension_mismatch_and_unknown_layout_are_rejected(self):
        for architecture in (
            replace(default_architecture_for_trigger(30), stem_length=18),
            replace(default_architecture_for_trigger(30), binding_layout="unknown"),
        ):
            with self.subTest(architecture=architecture), self.assertRaises(ValueError):
                architecture.validate_for_trigger(30)

    def test_automatic_linker_padding_preserves_aug_and_reporter_frame(self):
        for upper_stem2_length, padding in ((0, ""), (3, ""), (4, "XX"), (5, "X")):
            with self.subTest(upper_stem2_length=upper_stem2_length):
                request = make_legacy_request(upper_stem2_length=upper_stem2_length)
                design = build_design_problem(request)
                domains = switch_domains(design)
                layout = design.layout
                self.assertEqual(layout["frame_linker"], padding)
                self.assertEqual("frame_linker" in [domain.name for domain in domains], bool(padding))
                sequence = resolved_test_sequence(design)
                self.assertNotIn("X", sequence)
                self.assertEqual(sequence[layout["aug_start"]:layout["aug_start"] + 3], "AUG")
                reporter_start = sum(len(domain.pattern) for domain in domains[:-1])
                self.assertEqual((reporter_start - layout["aug_start"]) % 3, 0)
                self.assertEqual(reporter_start, layout["reporter_start"])
                self.assertEqual(sequence[reporter_start:], request.reporter.sequence_after_start_rna)
                self.assertEqual(len(sequence), sum(len(domain.pattern) for domain in domains))

    def test_native_keeps_explicit_custom_linker_bases_and_free_positions(self):
        for pattern in ("AGC", "CXX"):
            with self.subTest(pattern=pattern):
                request = make_request(frame_linker=pattern)
                design = build_design_problem(request)
                linker_domain = next(domain for domain in switch_domains(design) if domain.name == "frame_linker")
                self.assertEqual(linker_domain.pattern, pattern)

    def test_native_receives_variable_unpaired_prefix_and_fixed_user_sequences(self):
        request = make_request(t7_leader="GA", rbs_loop="AGGAGA")
        design = build_design_problem(request)
        domains = switch_domains(design)
        names = [domain.name for domain in domains]
        by_name = {domain.name: domain for domain in domains}
        self.assertEqual(by_name["t7_leader"].pattern, "GA")
        self.assertEqual(by_name["rbs"].pattern, "AGGAGA")
        self.assertEqual(by_name["rbs_prefix"].pattern, "XXX")
        self.assertEqual(by_name["upper_r"].pattern, "XXXXXX")
        self.assertEqual(names[names.index("rbs") - 1], "rbs_prefix")
        self.assertEqual(names[names.index("rbs") + 1:names.index("aug")], ["upper_r"])
        layout = design.layout
        targets = switch_targets(design)
        for name, structure in (
            ("OFF", targets["switch_off"].structure),
            ("ON", targets["switch_on"].structure.split("+")[1]),
        ):
            with self.subTest(state=name):
                self.assertEqual(len(structure), sum(len(domain.pattern) for domain in domains))
                self.assertEqual(structure[layout["rbs_prefix_start"]:layout["rbs_start"]], "...")
                self.assertEqual(structure[layout["rbs_start"]:layout["rbs_start"] + 6], "......")
        self.assertEqual(targets["switch_on"].structure.split("+")[0], "(" * 30)
        preview = build_switch_model(
            request.architecture, TRIGGER,
            reporter_sequence=request.reporter.sequence_after_start_rna,
        )
        self.assertEqual(preview.off_structure, targets["switch_off"].structure)
        self.assertEqual(preview.on_structure, targets["switch_on"].structure)

    def test_custom_spacing_bulge_and_legacy_prefix_keep_valid_target_lengths(self):
        for changes in (
            {"upper_stem_length": 8, "frame_linker": None},
            {"upper_stem_length": 4, "aug_spacer": 2, "frame_linker": None},
            {"rbs_prefix_length": 0, "upper_stem_length": 5, "frame_linker": None},
            {"rbs_prefix_length": 5, "bulge_length": 4, "t7_leader": "", "rbs_loop": "AGGA", "frame_linker": None},
        ):
            with self.subTest(changes=changes):
                request = make_request(**changes)
                design = build_design_problem(request)
                domains = switch_domains(design)
                expected_length = sum(len(domain.pattern) for domain in domains)
                targets = switch_targets(design)
                for target in (targets["switch_off"], targets["switch_on"]):
                    self.assertEqual(len(target.structure.split("+")[-1]), expected_length)
                    self.assertEqual(target.structure.count("("), target.structure.count(")"))
                sequence = resolved_test_sequence(design)
                self.assertEqual(len(sequence), expected_length)
                layout = design.layout
                self.assertEqual(sequence[layout["aug_start"]:layout["aug_start"] + 3], "AUG")
                if changes.get("rbs_prefix_length") == 0:
                    self.assertNotIn("rbs_prefix", [domain.name for domain in domains])

    def test_native_domains_match_shared_offsets_and_reporter_frame(self):
        request = make_request(t7_leader="GA", rbs_loop="AGGA", upper_stem_length=8,
                               frame_linker=CLASSIC_LINKER + "XX")
        # This fixture tests offsets with the fixed linker shifted relative to
        # its standard codon frame, so retain STOPs explicitly for this check.
        request = replace(request, exclude_stop_candidates=False)
        arch = request.architecture
        layout = architecture_layout(arch)
        design = build_design_problem(request)
        sequence = resolved_test_sequence(design)
        binding_start = len(arch.t7_leader)
        self.assertEqual(sequence[binding_start:binding_start + 30], reverse_complement_rna(TRIGGER))
        self.assertEqual(sequence[layout["rbs_start"]:layout["rbs_start"] + 4], "AGGA")
        self.assertEqual(sequence[layout["aug_start"]:layout["aug_start"] + 3], "AUG")
        self.assertNotIn("X", sequence)
        self.assertEqual(design.layout["rbs_start"], layout["rbs_start"])
        self.assertEqual(design.layout["aug_start"], layout["aug_start"])
        self.assertEqual(design.layout["reporter_start"], layout["reporter_start"])
        reporter_start = layout["reporter_start"]
        self.assertEqual(sequence[reporter_start:], request.reporter.sequence_after_start_rna)

    def test_rbs_linker_bounds_follow_domains_from_loop_start_to_before_reporter(self):
        for trigger_length, changes in (
            (24, {}),
            (30, {}),
            (30, {"rbs_prefix_length": 0, "t7_leader": "", "rbs_loop": "AGGA"}),
            (30, {"rbs_prefix_length": 5, "t7_leader": "GAGAGA", "rbs_loop": "AGGAGG", "upper_stem_length": 8, "frame_linker": None}),
            (24, {"frame_linker": "G", "upper_stem2_length": 5, "aug_spacer": 2}),
            (30, {"frame_linker": "ACC", "upper_stem2_length": 4}),
        ):
            with self.subTest(trigger_length=trigger_length, changes=changes):
                trigger = TRIGGER[:trigger_length]
                request = replace(
                    make_request(),
                    sequences=SequenceInput(aptamer_dna="ACGTACGT", trigger_dna=trigger),
                    architecture=configured_architecture(trigger_length, **changes),
                ).normalized()
                design = build_design_problem(request)
                domains = switch_domains(design)
                names = [domain.name for domain in domains]
                loop_index = names.index("rbs_prefix" if request.architecture.rbs_prefix_length else "rbs")
                expected_start = sum(len(domain.pattern) for domain in domains[:loop_index])
                expected_end = sum(len(domain.pattern) for domain in domains[:-1])
                layout = design.layout
                self.assertEqual(layout["rbs_linker_start"], expected_start)
                self.assertEqual(layout["rbs_linker_end"], expected_end)
                self.assertEqual(layout["reporter_start"], expected_end)
                self.assertEqual(layout["rbs_linker_length"], expected_end - expected_start)
                if not changes:
                    self.assertEqual(layout["rbs_linker_length"], 35 if trigger_length == 24 else 50)
                sequence = resolved_test_sequence(design)
                self.assertEqual(sequence[expected_end:], request.reporter.sequence_after_start_rna)



    def test_negative_or_empty_optimized_domains_are_rejected(self):
        for changes in (
            {"upper_stem_length": 0}, {"upper_stem2_length": -1},
            {"bulge_length": 0}, {"rbs_prefix_length": -1}, {"aug_spacer": -1},
        ):
            with self.subTest(changes=changes), self.assertRaises(ValueError):
                make_request(**changes)


    def test_rbs_linker_mfe_receives_actual_prefix_through_linker_and_no_reporter(self):
        for length, changes, reporter in (
            (24, {}, "CGUAAAGGCGAGGAGCUGUUC"), (30, {}, "CCC"),
            (30, {}, "UGC" * 30),
            (30, {"rbs_prefix_length": 0, "t7_leader": "", "rbs_loop": "AGGA"}, "UGC" * 10),
            (30, {"rbs_prefix_length": 5, "t7_leader": "GAGA", "rbs_loop": "AGGA", "frame_linker": "ACC"}, "UGC" * 5),
        ):
            with self.subTest(length=length, changes=changes):
                request = replace(make_request(), sequences=SequenceInput("ACGTACGT", TRIGGER[:length]),
                    architecture=configured_architecture(length, **changes),
                    reporter=ReporterContext("custom", "Custom", reporter),
                    thermo=replace(make_request().thermo, validation_max_complex_size=2),
                    exclude_stop_candidates=False, exclude_non_linear_candidates=False)
                problem = build_design_problem(request)
                assignment = problem.generate(random.Random(9))
                backend = AnalyticalBackend(request)
                backend.minimum = Mock(wraps=backend.minimum)
                candidate = evaluate_candidate(request, problem, assignment, backend, trial=1)
                region = next(call.args[0].strands[0] for call in backend.minimum.call_args_list
                              if call.args[0].strands[0].name == "RBS_linker")
                switch = assignment.sequence("Switch")
                self.assertEqual(region.sequence, switch[problem.layout["rbs_linker_start"]:problem.layout["reporter_start"]])
                self.assertEqual(region.material, request.thermo.switch_material)
                self.assertFalse(candidate.excluded)
                self.assertEqual(candidate.delta_g_rbs_linker, -1.)

    def test_stop_constraints_reject_fixed_stop_after_new_prefix_before_native_analysis(self):
        request = make_request(t7_leader="G", rbs_loop="AGGAGA", frame_linker="UAA")
        backend = AnalyticalBackend(request)
        with self.assertRaisesRegex(UnsatisfiableDesignError, "STOP"):
            design_candidates(request, backend)
        self.assertEqual(backend.calls, [])

    def test_native_constraints_accept_updated_domains_and_target_structures(self):
        for upper_stem2_length in (0, 3, 4, 5):
            with self.subTest(upper_stem2_length=upper_stem2_length):
                problem = build_design_problem(make_request(upper_stem2_length=upper_stem2_length))
                problem.validate(problem.generate(random.Random(12)))


class CalculatedStructureValidationTests(unittest.TestCase):
    def setUp(self):
        from aptaswitch_core.localization import language_context
        self.enterContext(language_context("fr"))

    def test_missing_and_malformed_structures_are_rejected_without_fallback(self):
        for structure in (None, "", "   ", "(((+....", "...x", "+...", "...+", "..++.."):
            with self.subTest(structure=structure), self.assertRaisesRegex(ValueError, "Structure calculée"):
                candidate_display_structure(structure)
        self.assertEqual(candidate_display_structure(" ((..)) + .. "), "((..))+..")

    def test_candidate_rejects_structure_length_and_strand_mismatches(self):
        request = make_request()
        model = build_switch_model(request.architecture, TRIGGER)
        switch = model.switch_sequence.replace("X", "A")
        for state, structure in (
            ("off", "." * (len(switch) - 1)),
            ("off", ".+" + "." * (len(switch) - 1)),
            ("on", "." * (len(TRIGGER) + len(switch))),
            ("on", "." * (len(TRIGGER) - 1) + "+" + "." * (len(switch) + 1)),
        ):
            with self.subTest(state=state, structure=structure), self.assertRaisesRegex(ValueError, "longueurs"):
                candidate_structure_svg(request.architecture, TRIGGER, switch, structure, state=state)

    def test_candidate_rejects_unresolved_or_missing_calculated_bases(self):
        request = make_request()
        model = build_switch_model(request.architecture, TRIGGER)
        for switch in ("", model.switch_sequence):
            with self.subTest(switch=switch), self.assertRaisesRegex(ValueError, "séquences calculées"):
                candidate_structure_svg(request.architecture, TRIGGER, switch, model.off_structure)


if __name__ == "__main__":
    unittest.main()
