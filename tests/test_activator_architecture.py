"""The activator defaults must reach actual design targets and frame scanning."""

from __future__ import annotations

import os
import unittest
from dataclasses import replace
from types import SimpleNamespace
from unittest.mock import Mock, patch

from aptaswitch_core.architecture import architecture_layout, default_architecture_for_trigger
from aptaswitch_core.models import DEFAULT_RBS_SEQUENCE, ReporterContext, RunRequest, SequenceInput
from aptaswitch_core.nupack_engine import (
    _build_design_spec,
    _create_model,
    nupack_import_context,
    run_nupack_engine,
)
from aptaswitch_core.preview import build_switch_model, candidate_display_structure, candidate_structure_svg
from aptaswitch_core.sequences import reading_frame_report, reverse_complement_rna


TRIGGER = "GTCCAGGCTGGTATAATTAGATCCACGTAC"


def make_request(**architecture_changes):
    return RunRequest(
        sequences=SequenceInput(aptamer_dna="ACGTACGTACGT", trigger_dna=TRIGGER),
        architecture=replace(default_architecture_for_trigger(len(TRIGGER)), **architecture_changes),
        trials=1,
        nupack_path="/test-only/nupack",
    ).normalized()


class CapturingNupack:
    class Domain:
        def __init__(self, sequence, *, name, material):
            if len(sequence) <= 1:
                raise ValueError("Empty NUPACK domains must be omitted.")
            self.sequence = sequence
            self.name = name
            self.material = material

    class TargetStrand:
        def __init__(self, domains, *, name):
            self.domains = domains
            self.name = name

    class TargetComplex:
        def __init__(self, strands, structure, *, name):
            self.strands = strands
            self.structure = structure
            self.name = name

    class TargetTube:
        def __init__(self, **kwargs):
            self.__dict__.update(kwargs)

    def SetSpec(self, **kwargs):
        return kwargs

    def DesignOptions(self, **kwargs):
        return kwargs

    def tube_design(self, **kwargs):
        return SimpleNamespace(**kwargs)


def switch_targets(design):
    return {
        target.name: target
        for tube in design["spec"].tubes
        for target in tube.on_targets
    }


def resolved_test_sequence(design):
    """Materialize captured domains with arbitrary A bases for offset checks only."""
    return "".join(domain.sequence[1:].replace("N", "A") for domain in design["switch"].domains)


class ActivatorArchitectureTests(unittest.TestCase):
    def test_standard_activator_matches_trigger_length_with_six_bases_to_aug(self):
        request = make_request()
        arch = request.architecture
        for trigger_length, split in ((24, (15, 9)), (30, (12, 18)), (32, (12, 20))):
            with self.subTest(trigger_length=trigger_length):
                standard = default_architecture_for_trigger(trigger_length)
                self.assertEqual((standard.toehold_length, standard.stem_length), split)
                self.assertEqual(standard.trigger_region_length, trigger_length)
                self.assertEqual(standard.upper_stem2_length, 0 if trigger_length == 30 else 4)
        self.assertEqual((arch.rbs_prefix_length, arch.rbs_to_aug_distance), (3, 6))
        self.assertEqual(arch.aug_spacer, 0)
        layout = architecture_layout(arch)
        self.assertEqual(layout["aug_start"] - layout["rbs_start"] - len(arch.rbs_loop), 6)
        self.assertEqual(request.to_dict()["architecture"]["rbs_prefix_length"], 3)

    def test_standard_30nt_has_exactly_18_bases_after_aug_before_reporter(self):
        request = make_request()
        arch = request.architecture
        self.assertEqual(arch.rbs_loop, "AGAGGAGA")
        self.assertEqual(arch.rbs_loop, DEFAULT_RBS_SEQUENCE)
        self.assertEqual(arch.upper_stem2_length, 0)
        self.assertIsNone(arch.frame_linker)
        self.assertEqual(arch.resolved_frame_linker(), "")
        design = _build_design_spec(CapturingNupack(), request, object())
        domains = design["switch"].domains
        names = [domain.name for domain in domains]
        self.assertNotIn("upper2_l", names)
        self.assertNotIn("upper2_r", names)
        self.assertNotIn("frame_linker", names)
        self.assertEqual(names[names.index("aug") + 1:], ["stem_r", "reporter_after_start"])
        self.assertEqual(next(domain.sequence for domain in domains if domain.name == "stem_r"), "r" + "N" * 18)
        self.assertEqual(next(domain.sequence for domain in domains if domain.name == "rbs"), "rAGAGGAGA")
        layout = design["layout"]
        sequence = resolved_test_sequence(design)
        reporter_start = sum(len(domain.sequence) - 1 for domain in domains[:-1])
        self.assertEqual(sequence[layout["aug_start"]:layout["aug_start"] + 3], "AUG")
        self.assertEqual(reporter_start - (layout["aug_start"] + 3), 18)
        self.assertEqual(sequence[reporter_start:], request.reporter.sequence_after_start_rna)
        self.assertEqual(layout["reporter_start_codon"], 8)
        for stop_position, is_premature in ((reporter_start - 3, True), (reporter_start, False)):
            with self.subTest(stop_position=stop_position):
                with_stop = sequence[:stop_position] + "UAA" + sequence[stop_position + 3:]
                report = reading_frame_report(with_stop, layout["aug_start"], layout["reporter_start_codon"])
                self.assertEqual(report["has_premature_stop"], is_premature)
        for target in (switch_targets(design)["switch_off"], switch_targets(design)["switch_on"]):
            self.assertEqual(len(target.structure.split("+")[-1]), len(sequence))
            self.assertEqual(target.structure.count("("), target.structure.count(")"))

    def test_custom_frame_linker_reaches_nupack_as_two_free_unpaired_bases(self):
        request = make_request(upper_stem2_length=4)
        self.assertIsNone(request.architecture.frame_linker)
        self.assertEqual(request.architecture.resolved_frame_linker(), "NN")
        design = _build_design_spec(CapturingNupack(), request, object())
        domains = design["switch"].domains
        linker_index = next(i for i, domain in enumerate(domains) if domain.name == "frame_linker")
        self.assertEqual(domains[linker_index].sequence, "rNN")
        self.assertEqual(domains[linker_index + 1].name, "reporter_after_start")
        linker_start = sum(len(domain.sequence) - 1 for domain in domains[:linker_index])
        for target in (switch_targets(design)["switch_off"], switch_targets(design)["switch_on"]):
            self.assertEqual(target.structure.split("+")[-1][linker_start:linker_start + 2], "..")

    def test_automatic_linker_padding_preserves_aug_and_reporter_frame(self):
        for upper_stem2_length, padding in ((0, ""), (3, ""), (4, "NN"), (5, "N")):
            with self.subTest(upper_stem2_length=upper_stem2_length):
                request = make_request(upper_stem2_length=upper_stem2_length)
                design = _build_design_spec(CapturingNupack(), request, object())
                domains = design["switch"].domains
                layout = design["layout"]
                self.assertEqual(layout["frame_linker"], padding)
                self.assertEqual("frame_linker" in [domain.name for domain in domains], bool(padding))
                sequence = resolved_test_sequence(design)
                self.assertNotIn("N", sequence)
                self.assertEqual(sequence[layout["aug_start"]:layout["aug_start"] + 3], "AUG")
                reporter_start = sum(len(domain.sequence) - 1 for domain in domains[:-1])
                self.assertEqual((reporter_start - layout["aug_start"]) % 3, 0)
                self.assertEqual(reporter_start, layout["aug_start"] + (layout["reporter_start_codon"] - 1) * 3)
                self.assertEqual(sequence[reporter_start:], request.reporter.sequence_after_start_rna)
                self.assertEqual(len(sequence), sum(len(domain.sequence) - 1 for domain in domains))

    def test_nupack_keeps_explicit_custom_linker_bases_and_free_positions(self):
        for pattern in ("AG", "CN"):
            with self.subTest(pattern=pattern):
                request = make_request(frame_linker=pattern)
                design = _build_design_spec(CapturingNupack(), request, object())
                linker_domain = next(domain for domain in design["switch"].domains if domain.name == "frame_linker")
                self.assertEqual(linker_domain.sequence, "r" + pattern)

    def test_open_legacy_run_keeps_its_original_offsets(self):
        architecture = replace(default_architecture_for_trigger(24), upper_stem_length=5, rbs_prefix_length=0)
        legacy = SimpleNamespace(**{
            name: value for name, value in vars(architecture).items()
            if name != "rbs_prefix_length"
        })
        legacy.trigger_region_length = 24
        legacy.resolved_frame_linker = architecture.resolved_frame_linker
        self.assertEqual(architecture_layout(legacy), architecture_layout(architecture))

    def test_nupack_receives_variable_unpaired_prefix_and_fixed_user_sequences(self):
        request = make_request(t7_leader="GA", rbs_loop="AGGAGA")
        design = _build_design_spec(CapturingNupack(), request, object())
        domains = design["switch"].domains
        names = [domain.name for domain in domains]
        by_name = {domain.name: domain for domain in domains}
        self.assertEqual(by_name["t7_leader"].sequence, "rGA")
        self.assertEqual(by_name["rbs"].sequence, "rAGGAGA")
        self.assertEqual(by_name["rbs_prefix"].sequence, "rNNN")
        self.assertEqual(by_name["upper_r"].sequence, "rNNNNNN")
        self.assertEqual(names[names.index("rbs") - 1], "rbs_prefix")
        self.assertEqual(names[names.index("rbs") + 1:names.index("aug")], ["upper_r"])
        layout = design["layout"]
        targets = switch_targets(design)
        for name, structure in (
            ("OFF", targets["switch_off"].structure),
            ("ON", targets["switch_on"].structure.split("+")[1]),
        ):
            with self.subTest(state=name):
                self.assertEqual(len(structure), sum(len(domain.sequence) - 1 for domain in domains))
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
            {"upper_stem_length": 8},
            {"upper_stem_length": 4, "aug_spacer": 2},
            {"rbs_prefix_length": 0, "upper_stem_length": 5},
            {"rbs_prefix_length": 5, "bulge_length": 4, "t7_leader": "", "rbs_loop": "AGGA"},
        ):
            with self.subTest(changes=changes):
                request = make_request(**changes)
                design = _build_design_spec(CapturingNupack(), request, object())
                domains = design["switch"].domains
                expected_length = sum(len(domain.sequence) - 1 for domain in domains)
                targets = switch_targets(design)
                for target in (targets["switch_off"], targets["switch_on"]):
                    self.assertEqual(len(target.structure.split("+")[-1]), expected_length)
                    self.assertEqual(target.structure.count("("), target.structure.count(")"))
                sequence = resolved_test_sequence(design)
                self.assertEqual(len(sequence), expected_length)
                layout = design["layout"]
                self.assertEqual(sequence[layout["aug_start"]:layout["aug_start"] + 3], "AUG")
                if changes.get("rbs_prefix_length") == 0:
                    self.assertNotIn("rbs_prefix", [domain.name for domain in domains])

    def test_nupack_domains_match_shared_offsets_and_reporter_frame(self):
        request = make_request(t7_leader="GA", rbs_loop="AGGA", upper_stem_length=8)
        arch = request.architecture
        layout = architecture_layout(arch)
        design = _build_design_spec(CapturingNupack(), request, object())
        sequence = resolved_test_sequence(design)
        binding_start = len(arch.t7_leader)
        self.assertEqual(sequence[binding_start:binding_start + 30], reverse_complement_rna(TRIGGER))
        self.assertEqual(sequence[layout["rbs_start"]:layout["rbs_start"] + 4], "AGGA")
        self.assertEqual(sequence[layout["aug_start"]:layout["aug_start"] + 3], "AUG")
        self.assertNotIn("N", sequence)
        self.assertEqual(design["layout"]["rbs_start"], layout["rbs_start"])
        self.assertEqual(design["layout"]["aug_start"], layout["aug_start"])
        self.assertEqual(design["layout"]["reporter_start_codon"], layout["reporter_start_codon"])
        reporter_start = layout["aug_start"] + (layout["reporter_start_codon"] - 1) * 3
        self.assertEqual(sequence[reporter_start:], request.reporter.sequence_after_start_rna)

    def test_rbs_linker_bounds_follow_domains_from_loop_start_to_before_reporter(self):
        for trigger_length, changes in (
            (24, {}),
            (30, {}),
            (30, {"rbs_prefix_length": 0, "t7_leader": "", "rbs_loop": "AGGA"}),
            (30, {"rbs_prefix_length": 5, "t7_leader": "GAGAGA", "rbs_loop": "AGGAGG", "upper_stem_length": 8}),
            (24, {"frame_linker": "G", "upper_stem2_length": 5, "aug_spacer": 2}),
            (30, {"frame_linker": "AC", "upper_stem2_length": 4}),
        ):
            with self.subTest(trigger_length=trigger_length, changes=changes):
                trigger = TRIGGER[:trigger_length]
                request = replace(
                    make_request(),
                    sequences=SequenceInput(aptamer_dna="ACGTACGT", trigger_dna=trigger),
                    architecture=replace(default_architecture_for_trigger(trigger_length), **changes),
                ).normalized()
                design = _build_design_spec(CapturingNupack(), request, object())
                domains = design["switch"].domains
                names = [domain.name for domain in domains]
                loop_index = names.index("rbs_prefix" if request.architecture.rbs_prefix_length else "rbs")
                expected_start = sum(len(domain.sequence) - 1 for domain in domains[:loop_index])
                expected_end = sum(len(domain.sequence) - 1 for domain in domains[:-1])
                layout = design["layout"]
                self.assertEqual(layout["rbs_linker_start"], expected_start)
                self.assertEqual(layout["rbs_linker_end"], expected_end)
                self.assertEqual(layout["reporter_start"], expected_end)
                self.assertEqual(layout["rbs_linker_length"], expected_end - expected_start)
                if not changes:
                    self.assertEqual(layout["rbs_linker_length"], 35 if trigger_length == 24 else 38)
                sequence = resolved_test_sequence(design)
                self.assertEqual(sequence[expected_end:], request.reporter.sequence_after_start_rna)

    def test_rbs_linker_mfe_receives_actual_prefix_through_linker_and_no_reporter(self):
        for trigger_length, changes, reporter_sequence in (
            (24, {}, "CGUAAAGGCGAGGAGCUGUUC"),
            (30, {}, "CCC"),
            (30, {}, "UGC" * 30),
            (30, {"rbs_prefix_length": 0, "t7_leader": "", "rbs_loop": "AGGA"}, "UGC" * 10),
            (30, {"rbs_prefix_length": 5, "t7_leader": "GAGA", "rbs_loop": "AGGA", "frame_linker": "AC"}, "UGC" * 5),
        ):
            with self.subTest(trigger_length=trigger_length, changes=changes, reporter_length=len(reporter_sequence)):
                trigger = TRIGGER[:trigger_length]
                request = replace(
                    make_request(),
                    sequences=SequenceInput(aptamer_dna="ACGTACGT", trigger_dna=trigger),
                    architecture=replace(default_architecture_for_trigger(trigger_length), **changes),
                    reporter=ReporterContext(key="custom", label="Custom reporter", sequence_after_start_rna=reporter_sequence),
                ).normalized()
                design = _build_design_spec(CapturingNupack(), request, object())
                domains = design["switch"].domains
                names = [domain.name for domain in domains]
                # Distinct resolved bases make a missing prefix or upstream
                # stem inclusion visible in the actual MFE call.
                parts = [
                    domain.sequence[1:].replace("N", "G" if domain.name == "rbs_prefix" else "C" if domain.name == "upper_l" else "A")
                    for domain in domains
                ]
                sequence = "".join(parts)
                loop_index = names.index("rbs_prefix" if request.architecture.rbs_prefix_length else "rbs")
                expected_segment = "".join(parts[loop_index:-1])
                calculated = SimpleNamespace(
                    to_analysis=lambda strand: sequence if strand is design["switch"] else trigger,
                    defects=SimpleNamespace(ensemble_defect=.01),
                )
                design["spec"] = Mock()
                design["spec"].run.return_value = [calculated]
                nupack = SimpleNamespace(
                    Model=Mock(return_value=object()),
                    mfe=Mock(side_effect=lambda *, strands, model: [SimpleNamespace(
                        energy=-1.0, structure="+".join("." * (len(seq) - 1) for seq in strands),
                    )]),
                )
                with (
                    patch("aptaswitch_core.nupack_engine.nupack_import_context") as context,
                    patch("aptaswitch_core.nupack_engine._build_design_spec", return_value=design),
                    patch("aptaswitch_core.nupack_engine._complex_pct", return_value=10.0),
                ):
                    context.return_value.__enter__.return_value = nupack
                    result = run_nupack_engine(request)
                self.assertFalse(result.candidates[0].excluded)
                self.assertEqual(nupack.mfe.call_count, 4)
                self.assertEqual(nupack.mfe.call_args_list[-1].kwargs["strands"], ["r" + expected_segment])
                self.assertEqual(result.candidates[0].delta_g_rbs_linker, -1.0)

    def test_stop_filter_scans_after_new_prefix_before_costly_nupack_analysis(self):
        request = make_request(t7_leader="G", rbs_loop="AGGAGA")
        design = _build_design_spec(CapturingNupack(), request, object())
        layout = design["layout"]
        nominal = resolved_test_sequence(design)
        sequence = (
            nominal[:layout["aug_start"]] + "AUGUAA"
            + nominal[layout["aug_start"] + 6:]
        )
        self.assertTrue(reading_frame_report(sequence, layout["aug_start"], layout["reporter_start_codon"])["has_premature_stop"])
        result = SimpleNamespace(to_analysis=lambda strand: sequence if strand is design["switch"] else TRIGGER)
        design["spec"] = Mock()
        design["spec"].run.return_value = [result]
        nupack = SimpleNamespace(Model=Mock(return_value=object()), mfe=Mock())
        with (
            patch("aptaswitch_core.nupack_engine.nupack_import_context") as context,
            patch("aptaswitch_core.nupack_engine._build_design_spec", return_value=design),
            patch("aptaswitch_core.nupack_engine._complex_pct") as tube_analysis,
        ):
            context.return_value.__enter__.return_value = nupack
            run = run_nupack_engine(request)
        self.assertTrue(run.candidates[0].excluded)
        self.assertEqual(run.candidates[0].first_stop, 2)
        nupack.mfe.assert_not_called()
        tube_analysis.assert_not_called()

    def test_negative_or_empty_optimized_domains_are_rejected(self):
        for changes in (
            {"upper_stem_length": 0}, {"upper_stem2_length": -1},
            {"bulge_length": 0}, {"rbs_prefix_length": -1}, {"aug_spacer": -1},
        ):
            with self.subTest(changes=changes), self.assertRaises(ValueError):
                make_request(**changes)

    @unittest.skipUnless(os.environ.get("APTASWITCH_TEST_NUPACK"), "User-provided NUPACK not configured")
    def test_real_nupack_accepts_updated_domains_and_target_structures(self):
        with nupack_import_context(os.environ["APTASWITCH_TEST_NUPACK"]) as nupack:
            for upper_stem2_length in (0, 3, 4, 5):
                with self.subTest(upper_stem2_length=upper_stem2_length):
                    request = make_request(upper_stem2_length=upper_stem2_length)
                    design = _build_design_spec(nupack, request, _create_model(nupack, request))
                    self.assertIsNotNone(design["spec"])


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
        switch = model.switch_sequence.replace("N", "A")
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
