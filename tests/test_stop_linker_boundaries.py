"""STOP detection follows AUG across custom domains, never substring matches."""

import csv
from dataclasses import replace
import json
from pathlib import Path
import random
import tempfile
import unittest
from unittest.mock import patch

from aptaswitch_core.architecture import architecture_layout, default_architecture_for_trigger
from aptaswitch_core.exports import export_run
from aptaswitch_core.localization import language_context
from aptaswitch_core.models import ReporterContext, RunRequest, SequenceInput
from aptaswitch_core.preview import premature_stop_positions
from aptaswitch_core.sequences import reading_frame_report
from aptaswitch_core.thermodynamics.design_constraints import (
    DesignAssignment, DesignProblem, UnsatisfiableDesignError, build_design_problem,
)
from aptaswitch_core.thermodynamics.design_engine import evaluate_candidate
from aptaswitch_core.thermodynamics.types import Strand
from tests.design_fixtures import design_result_fixture
from tests.test_thermodynamics_design_engine import AnalyticalBackend


STOPS = ("UAA", "UAG", "UGA")


def request(length=30, *, linker=None, reporter="GCU" * 3, **kwargs):
    architecture = default_architecture_for_trigger(length)
    if linker is not None:
        architecture = replace(architecture, frame_linker=linker)
    options = RunRequest(
        SequenceInput("ACGTACGT", "G" * length), architecture,
        reporter=ReporterContext("test", "Test", reporter, reporter),
        trials=1, **kwargs,
    )
    return replace(options, scoring=replace(options.scoring, rbs_calculator=0.))


def fixed_switch(options, *, actual_linker=None):
    """Synthetic supplied switch; no assertion about its thermodynamics."""
    layout = architecture_layout(options.architecture)
    sequence = list("C" * layout["reporter_start"])
    for start, motif in (
        (0, options.architecture.t7_leader),
        (layout["rbs_start"], options.architecture.rbs_loop),
        (layout["aug_start"], "AUG"),
        (layout["reporter_start"] - len(layout["frame_linker"]),
         actual_linker if actual_linker is not None else layout["frame_linker"].replace("X", "C")),
    ):
        sequence[start:start + len(motif)] = motif
    return "".join(sequence) + options.reporter.sequence_after_start_rna


class StopLinkerBoundaryTests(unittest.TestCase):
    def setUp(self):
        self.enterContext(language_context("fr"))

    def test_all_three_stops_in_fixed_custom_linker_fail_before_search(self):
        for stop in STOPS:
            for bases in (stop, stop.replace("U", "T")):
                options = request(linker=bases + "XXX")
                with self.subTest(linker=bases), patch.object(
                    DesignProblem, "generate", side_effect=AssertionError("No search for an unavoidable STOP")
                ), self.assertRaisesRegex(UnsatisfiableDesignError, "STOP prématuré au codon 5"):
                    build_design_problem(options)
                self.assertEqual(options.architecture.frame_linker, bases + "XXX")

    def test_all_three_stops_across_stem_and_custom_linker_are_constrained(self):
        for stop in STOPS:
            options = replace(request(24, linker=stop[1:]),
                              sequences=SequenceInput("ACGTACGT", "T" * 24))
            with self.subTest(stop=stop), self.assertRaisesRegex(UnsatisfiableDesignError, "STOP"):
                build_design_problem(options)
            unfiltered = build_design_problem(replace(options, exclude_stop_candidates=False))
            value = unfiltered.generate(random.Random(17))
            linker = next(region for region in unfiltered.regions if region.name == "frame_linker")
            sequence = value.sequence("Switch")
            self.assertEqual(sequence[linker.start - 1:linker.stop], stop)
            report = unfiltered.reading_frame_report(value)
            self.assertTrue(report["has_premature_stop"])
            self.assertEqual(report["stops_found"][-1]["position"], linker.start - 1)

    def test_optimized_and_partly_fixed_linkers_generate_and_mutate_without_stops(self):
        for linker in ("X" * 21, "UXXCCC", "UAXCCC"):
            with self.subTest(linker=linker):
                options = request(linker=linker)
                problem = build_design_problem(options)
                rng = random.Random(31)
                layout = problem.layout
                for _ in range(30):
                    value = problem.generate(rng)
                    for _ in range(3):
                        coding = value.sequence("Switch")[layout["aug_start"]:layout["reporter_start"]]
                        self.assertFalse(any(coding[i:i + 3] in STOPS for i in range(0, len(coding), 3)))
                        self.assertTrue(problem.passes_sequence_screen(value))
                        value = problem.mutate_guided(value, rng, ())

    def test_offsets_follow_aug_and_only_in_frame_linker_stops_are_premature(self):
        for leader in ("", "G", "GGGA"):
            for stop in STOPS:
                for phase in (0, 1, 2):
                    options = request(linker="C" * 9)
                    options = replace(options, architecture=replace(options.architecture, t7_leader=leader))
                    layout = architecture_layout(options.architecture)
                    start = layout["reporter_start"] - 9 + phase
                    sequence = fixed_switch(options)
                    sequence = sequence[:start] + stop + sequence[start + 3:]
                    with self.subTest(leader=leader, stop=stop, phase=phase):
                        report = reading_frame_report(sequence.replace("U", "T"), layout["aug_start"],
                                                      reporter_start=layout["reporter_start"])
                        self.assertEqual(report["has_premature_stop"], phase == 0)
                        self.assertEqual(premature_stop_positions(options.architecture, sequence),
                                         list(range(start, start + 3)) if phase == 0 else [])
                        if phase == 0:
                            self.assertEqual(report["first_stop"], 5)

    def test_literal_stop_across_aligned_linker_reporter_boundary_is_out_of_frame(self):
        for stop in STOPS:
            reporter = stop[1:] + "CCUAUAA"
            options = request(linker="CCU", reporter=reporter)
            sequence = fixed_switch(options)
            layout = architecture_layout(options.architecture)
            with self.subTest(stop=stop):
                self.assertEqual(sequence[layout["reporter_start"] - 1:layout["reporter_start"] + 2], stop)
                report = reading_frame_report(sequence, layout["aug_start"], reporter_start=layout["reporter_start"])
                self.assertFalse(report["has_premature_stop"])
                self.assertEqual(report["stops_found"][-1]["codon"], "UAA")
                self.assertEqual(report["stops_found"][-1]["position"], len(sequence) - 3)
                self.assertEqual(premature_stop_positions(options.architecture, sequence), [])

    def test_unaligned_raw_reporter_boundary_still_flags_a_crossing_stop(self):
        for stop in STOPS:
            for reporter_start in (5, 6):
                # AUG starts at1; the following STOP starts at4 and overlaps
                # an independently declared, unaligned reporter boundary.
                with self.subTest(stop=stop, boundary=reporter_start):
                    report = reading_frame_report("CAUG" + stop + "CCC", 1, reporter_start=reporter_start)
                    self.assertTrue(report["has_premature_stop"])
                    self.assertEqual(report["first_stop"], 2)
                    self.assertEqual(report["stops_found"][0]["position"], 4)
        for length in (1, 2, 4, 5):
            options = request(linker="X" * length)
            with self.subTest(linker_length=length), self.assertRaisesRegex(ValueError, "hors du cadre de lecture"):
                options.normalized()

    def test_reading_keeps_all_three_linker_stops_and_exports_the_detected_status(self):
        for stop in STOPS:
            with self.subTest(stop=stop), tempfile.TemporaryDirectory() as folder:
                options = request(linker="XXX", reporter="GCU" * 2 + "UAA",
                                  output_dir=Path(folder).resolve())
                sequence = fixed_switch(options, actual_linker=stop)
                options = replace(options, fixed_switch_sequence=sequence).normalized()
                problem = build_design_problem(options)
                value = DesignAssignment((Strand(options.sequences.trigger_dna, "dna", "Trigger"),
                                          Strand(sequence, "rna", "Switch")))
                candidate = evaluate_candidate(options, problem, value, AnalyticalBackend(options), trial=1)
                self.assertTrue(candidate.has_stop_codon)
                self.assertEqual(candidate.first_stop, 5)
                self.assertFalse(candidate.excluded)
                self.assertIsNotNone(candidate.score)
                self.assertEqual(candidate.switch_seq, sequence)
                result = design_result_fixture(options)
                result.candidates = [candidate]
                exports = export_run(result)
                payload = json.loads(Path(exports["json"]).read_text())
                self.assertTrue(payload["results"][0]["has_stop_codon"])
                self.assertEqual(payload["results"][0]["first_stop"], 5)
                with Path(exports["csv"]).open(newline="") as stream:
                    row = next(csv.DictReader(stream))
                self.assertEqual((row["STOP_codon"], row["First_STOP"], row["Status"]),
                                 ("yes", "5", "evaluated"))

    def test_terminal_reporter_stop_is_retained_without_a_premature_stop_flag(self):
        for stop in STOPS:
            options = request(linker="CCC", reporter="GCU" * 2 + stop)
            sequence = fixed_switch(options)
            options = replace(options, fixed_switch_sequence=sequence).normalized()
            problem = build_design_problem(options)
            value = DesignAssignment((Strand(options.sequences.trigger_dna, "dna", "Trigger"),
                                      Strand(sequence, "rna", "Switch")))
            with self.subTest(terminal_stop=stop):
                candidate = evaluate_candidate(options, problem, value, AnalyticalBackend(options), trial=1)
                self.assertFalse(candidate.has_stop_codon)
                self.assertEqual(candidate.first_stop, 8)
                self.assertFalse(candidate.excluded)
                self.assertIsNotNone(candidate.score)
                self.assertEqual(candidate.switch_seq, sequence)


if __name__ == "__main__":
    unittest.main()
