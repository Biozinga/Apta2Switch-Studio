"""Automatic reading annotations must come from the supplied sequences."""

from dataclasses import replace
from pathlib import Path
import tempfile
import unittest

from aptaswitch_core.architecture import architecture_layout, default_architecture_for_trigger
from aptaswitch_core.localization import get_language, set_language
from aptaswitch_core.models import REPORTERS, ReporterContext, RunRequest, SequenceInput
from aptaswitch_core.reading import prepare_reading_request
from aptaswitch_core.saved_switch_analysis import build_saved_switch_context
from aptaswitch_core.sequences import reverse_complement_rna
from aptaswitch_studio import saved_switch_library
from tests.design_fixtures import design_result_fixture
from tests.test_reading_request import supplied_switch


TRIGGER = "GTCCAGGCTGGTATAATTAGATCCACGTAC"


def supplied_request(length=30, *, linker=None, leader="GGG", reporter=None, suffix=None):
    trigger = TRIGGER[:length]
    architecture = default_architecture_for_trigger(length)
    architecture = replace(architecture, t7_leader=leader)
    if linker is not None:
        architecture = replace(architecture, frame_linker=linker)
    request = RunRequest(SequenceInput("AACCGGTTAACC", trigger), architecture,
                         reporter=reporter or REPORTERS["sfgfp"])
    sequence = list(supplied_switch(request, suffix=suffix))
    start = len(leader)
    sequence[start:start + length] = reverse_complement_rna(trigger)
    return request, "".join(sequence)


class ReadingInferenceTests(unittest.TestCase):
    def setUp(self):
        self.addCleanup(set_language, get_language())
        set_language("fr")

    def test_known_12_24_and_30nt_families_ignore_unrelated_form_architecture(self):
        for length in (12, 24, 30):
            for reporter in (REPORTERS["sfgfp"], REPORTERS["mng"]):
                with self.subTest(length=length, reporter=reporter.key):
                    source, switch = supplied_request(length, reporter=reporter)
                    unrelated = replace(default_architecture_for_trigger(18), t7_leader="AAA",
                                        rbs_loop="CCCC", frame_linker="GCA" * 15)
                    form = replace(source, architecture=unrelated)
                    reading = prepare_reading_request(form, switch)
                    self.assertEqual(reading.fixed_switch_sequence, switch)
                    self.assertEqual(reading.sequences.trigger_dna, source.sequences.trigger_dna)
                    self.assertEqual(architecture_layout(reading.architecture),
                                     architecture_layout(replace(source.architecture,
                                         frame_linker=source.architecture.resolved_frame_linker().replace("X", "A"))))
                    self.assertEqual(reading.reporter, source.reporter)

    def test_actual_leader_is_recognized_without_a_form_control(self):
        for leader in ("", "AAA", "GGGAAGCU"):
            with self.subTest(leader=leader):
                source, switch = supplied_request(leader=leader)
                reading = prepare_reading_request(source, switch.lower().replace("u", "t"))
                self.assertEqual(reading.architecture.t7_leader, leader)
                self.assertEqual(reading.fixed_switch_sequence, switch)

    def test_selected_reporter_identifies_arbitrary_optimized_linker_bases_and_lengths(self):
        for linker in ("GCA" * 7, "GUC" * 2, "", "GCA" * 8):
            with self.subTest(linker=linker):
                source, switch = supplied_request(linker=linker)
                reading = prepare_reading_request(source, switch)
                self.assertEqual(reading.architecture.frame_linker, linker)
                self.assertEqual(reading.reporter.sequence_after_start_rna,
                                 source.reporter.sequence_after_start_rna)

    def test_historical_21nt_reporter_and_complete_cds_are_both_recognized(self):
        reporter = REPORTERS["sfgfp"]
        for suffix in (reporter.full_sequence_after_start_rna[:21],
                       reporter.full_sequence_after_start_rna[:45],
                       reporter.full_sequence_after_start_rna):
            with self.subTest(suffix_length=len(suffix)):
                source, switch = supplied_request(suffix=suffix)
                reading = prepare_reading_request(source, switch)
                self.assertEqual(reading.reporter.sequence_after_start_rna, suffix)
                self.assertEqual(reading.reporter.full_sequence_after_start_rna,
                                 reporter.full_sequence_after_start_rna)

    def test_complete_unknown_reporter_uses_exact_classic_linker_marker(self):
        source, switch = supplied_request(suffix="CCCGGGAUCUAA")
        reading = prepare_reading_request(source, switch, reporter_is_complete=True)
        self.assertEqual(reading.reporter.key, "reading_reporter")
        self.assertEqual(reading.reporter.full_sequence_after_start_rna, "CCCGGGAUCUAA")
        self.assertIsNone(reading.reporter.source_accession)
        self.assertEqual(reading.fixed_switch_sequence, switch)

    def test_complete_known_custom_reporter_can_identify_an_optimized_linker(self):
        reporter = ReporterContext(key="custom", label="Custom", sequence_after_start_rna="CCCGGGAUCUAA",
                                   full_sequence_after_start_rna="CCCGGGAUCUAA")
        source, switch = supplied_request(linker="GCA" * 7, reporter=reporter)
        reading = prepare_reading_request(source, switch, reporter_is_complete=True)
        self.assertEqual(reading.architecture.frame_linker, "GCA" * 7)
        self.assertEqual(reading.reporter.full_sequence_after_start_rna, reporter.full_sequence_after_start_rna)

    def test_unknown_complete_reporter_with_optimized_linker_is_not_guessed(self):
        source, switch = supplied_request(linker="GCA" * 7, suffix="CCCGGGAUCUAA")
        with self.assertRaisesRegex(ValueError, "frontière.*linker.*rapporteur"):
            prepare_reading_request(source, switch, reporter_is_complete=True)

    def test_short_coincidental_reporter_match_is_not_used_as_an_annotation(self):
        source, switch = supplied_request(suffix=REPORTERS["sfgfp"].sequence_after_start_rna[:20])
        with self.assertRaisesRegex(ValueError, "au moins 21 bases"):
            prepare_reading_request(source, switch)

    def test_multiple_valid_trigger_anchors_are_rejected(self):
        source, switch = supplied_request()
        layout = architecture_layout(source.architecture)
        linker_start = layout["aug_start"] + 3 + source.architecture.stem_length
        # Both activator interpretations have the same final reporter. The
        # first interpretation has a long linker containing the second switch.
        ambiguous = switch[:linker_start] + "A" + switch
        with self.assertRaisesRegex(ValueError, "Plusieurs découpages"):
            prepare_reading_request(source, ambiguous)

    def test_multiple_reporter_boundaries_are_rejected_without_choosing_a_default(self):
        reporter = ReporterContext(key="repeat", label="Repeat", sequence_after_start_rna="C" * 30,
                                   full_sequence_after_start_rna="C" * 60)
        source, switch = supplied_request(linker="CCC" * 7, reporter=reporter)
        with self.assertRaisesRegex(ValueError, "Plusieurs découpages"):
            prepare_reading_request(source, switch)

    def test_incompatible_trigger_or_rbs_is_reported_without_repairing_sequence(self):
        source, switch = supplied_request()
        layout = architecture_layout(source.architecture)
        for position in (len(source.architecture.t7_leader) + 2, layout["rbs_start"], layout["aug_start"]):
            changed = switch[:position] + ("A" if switch[position] != "A" else "C") + switch[position + 1:]
            with self.subTest(position=position), self.assertRaisesRegex(ValueError, "pas reconnue automatiquement"):
                prepare_reading_request(source, changed)

    def test_inferred_annotation_survives_native_saved_switch_roundtrip(self):
        source, switch = supplied_request(linker="GCA" * 7)
        with tempfile.TemporaryDirectory() as folder:
            source = replace(source, output_dir=Path(folder).resolve() / "runs",
                             architecture=default_architecture_for_trigger(24))
            request = prepare_reading_request(source, switch)
            result = design_result_fixture(request)
            result.optimization = {
                "algorithm": "fixed_sequence_evaluation", "objective": "supplied_sequences",
                "sequence_optimization": False, "evaluations": 1, "trials": [],
            }
            candidate = result.candidates[0]
            candidate.switch_seq = switch
            snapshot = saved_switch_library.candidate_snapshot(result, candidate, "Auto")
            saved = saved_switch_library.save_switch(Path(folder).resolve() / "saved", "Auto", snapshot)
            reopened = saved_switch_library.load_switch(saved.path)
            context = build_saved_switch_context(reopened.snapshot, {})
            self.assertEqual(context["designed_switch_sequence"], switch)
            self.assertEqual(context["reporter_start"], 83)
            self.assertEqual(context["aug_start"], 50)
            self.assertEqual(context["trigger_sequence"], source.sequences.trigger_dna)
            self.assertEqual(reopened.snapshot["request"]["architecture"]["frame_linker"], "GCA" * 7)


if __name__ == "__main__":
    unittest.main()
