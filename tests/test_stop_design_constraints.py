"""STOP-free proposals, exact coding boundaries, and impossible fixed inputs."""
from dataclasses import replace
import random
import unittest
from unittest.mock import patch

from aptaswitch_core.architecture import default_architecture_for_trigger
from aptaswitch_core.models import InverseExtensionSettings, ReporterContext, RunRequest, SequenceInput
from aptaswitch_core.sequences import reading_frame_report
from aptaswitch_core.thermodynamics.design_constraints import (
    DesignProblem, UnsatisfiableDesignError, build_design_problem,
)


def request(length=30, *, linker=None, **changes):
    arch = default_architecture_for_trigger(length)
    if linker is not None:
        arch = replace(arch, frame_linker=linker)
    return RunRequest(SequenceInput("ACGT", "G" * length), arch,
                      reporter=ReporterContext("test", "Test", "GCU" * 4), **changes)


class StopDesignConstraintTests(unittest.TestCase):
    def test_exact_boundary_counts_crossing_stop_but_not_reporter_stop(self):
        sequence = "AUGUAAUGG"
        for boundary, expected in ((3, False), (4, True), (5, True), (6, True)):
            with self.subTest(boundary=boundary):
                report = reading_frame_report(sequence, 0, reporter_start=boundary)
                self.assertEqual(report["has_premature_stop"], expected)
                self.assertEqual(report["first_stop"], 2)
        self.assertFalse(reading_frame_report("AUGCUACAA", 0, reporter_start=6)["has_premature_stop"])

    def test_misaligned_custom_linker_fails_without_changing_input(self):
        for linker in ("U", "UA", "XXXX"):
            for exclude in (True, False):
                req = request(linker=linker, exclude_stop_candidates=exclude)
                with self.subTest(linker=linker, exclude=exclude), self.assertRaisesRegex(ValueError, "multiple .*3"):
                    build_design_problem(req)
                self.assertEqual(req.architecture.frame_linker, linker)

    def test_impossible_fixed_stop_fails_before_generation_and_opt_out_keeps_it(self):
        req = request(linker="UAAXXX")
        with patch.object(DesignProblem, "generate", side_effect=AssertionError("Must fail before proposals")):
            with self.assertRaisesRegex(UnsatisfiableDesignError, "STOP"):
                build_design_problem(req)
        unfiltered = build_design_problem(replace(req, exclude_stop_candidates=False))
        value = unfiltered.generate(random.Random(17))
        self.assertTrue(unfiltered.reading_frame_report(value)["has_premature_stop"])
        self.assertTrue(unfiltered.passes_sequence_screen(value))
        self.assertFalse(any(lib.name.startswith("sense_codon_") for lib in unfiltered.libraries))

    def test_safe_fixed_codons_add_no_constraint_components(self):
        req = request()
        filtered = build_design_problem(req)
        unfiltered = build_design_problem(replace(req, exclude_stop_candidates=False))
        self.assertEqual(filtered.libraries, unfiltered.libraries)
        self.assertEqual(filtered.mutable_components(), unfiltered.mutable_components())

    def test_pairing_and_fixed_linker_can_prove_junction_stop_impossible(self):
        # The fixed DNA recognition arm forces U at the last stem nucleotide;
        # the two fixed linker As complete UAA across the domain boundary.
        req = replace(request(24, linker="AA"), sequences=SequenceInput("ACGT", "T" * 24))
        with self.assertRaisesRegex(UnsatisfiableDesignError, "STOP"):
            build_design_problem(req)
        unfiltered = build_design_problem(replace(req, exclude_stop_candidates=False))
        value = unfiltered.generate(random.Random(9))
        self.assertTrue(unfiltered.reading_frame_report(value)["has_premature_stop"])

    def test_variable_linker_and_inverse_mutations_preserve_fixed_bases_and_sense_codons(self):
        req = replace(request(linker="X" * 21), sequences=SequenceInput("ACGT", "G" * 24),
                      inverse_extension=InverseExtensionSettings(30, "both"))
        problem = build_design_problem(req)
        rng = random.Random(17)
        for _ in range(50):
            value = problem.generate(rng)
            for step in range(5):
                self.assertTrue(problem.passes_sequence_screen(value))
                problem.validate(value)
                for strand, pattern in zip(value.strands, problem.strands):
                    for actual, expected in zip(strand.sequence, pattern.pattern):
                        if expected != "X":
                            self.assertEqual(actual, expected)
                value = (problem.mutate(value, rng) if step % 2 else
                         problem.mutate_guided(value, rng, ()))

    def test_stop_constraint_spans_stem_linker_junction(self):
        req = replace(request(24), sequences=SequenceInput("ACGT", "T" * 24))
        problem = build_design_problem(req)
        linker = next(region for region in problem.regions if region.name == "frame_linker")
        crossing = next(lib for lib in problem.libraries if lib.name == f"sense_codon_{linker.start-1}")
        self.assertEqual(tuple(site.index for site in crossing.sites),
                         (linker.start - 1, linker.start, linker.start + 1))
        rng = random.Random(3)
        for _ in range(100):
            value = problem.generate(rng)
            self.assertFalse(problem.reading_frame_report(value)["has_premature_stop"])
            self.assertEqual(value.sequence("Switch")[linker.start-1], "U")


if __name__ == "__main__":
    unittest.main()
