"""Exact sequence and reporter provenance in imported switch analyses."""

from dataclasses import replace
from pathlib import Path
import tempfile
import unittest

from aptaswitch_core.architecture import architecture_layout, default_architecture_for_trigger
from aptaswitch_core.localization import get_language, set_language
from aptaswitch_core.models import (
    REPORTERS, InverseExtensionSettings, ReporterContext, RunRequest,
    SequenceInput, ThermoConditions,
)
from aptaswitch_core.reading import prepare_reading_request as _prepare_reading_request
from aptaswitch_core.saved_switch_analysis import build_saved_switch_context
from aptaswitch_studio import saved_switch_library
from tests.design_fixtures import design_result_fixture


def prepare_manually_annotated_reading(*args, **kwargs):
    return _prepare_reading_request(*args, infer_architecture=False, **kwargs)


def supplied_switch(request, suffix=None, linker=None):
    """An arbitrary sequence with the selected architecture's fixed motifs."""
    architecture = request.architecture
    layout = architecture_layout(architecture)
    sequence = list("C" * layout["reporter_start"])
    for start, motif in (
        (0, architecture.t7_leader),
        (layout["rbs_start"], architecture.rbs_loop),
        (layout["aug_start"], "AUG"),
        (layout["reporter_start"] - len(layout["frame_linker"]),
         linker if linker is not None else layout["frame_linker"].replace("X", "A")),
    ):
        sequence[start:start + len(motif)] = motif
    suffix = request.reporter.sequence_after_start_rna if suffix is None else suffix
    switch = "".join(sequence) + suffix
    return switch.replace("T", "U") if request.thermo.switch_material == "rna" else switch.replace("U", "T")


class ReadingRequestTests(unittest.TestCase):
    def setUp(self):
        self.addCleanup(set_language, get_language())
        set_language("fr")

    def request(self, **overrides):
        return RunRequest(
            sequences=SequenceInput("AACCGGTTAACC", "AACC" * 6),
            architecture=default_architecture_for_trigger(24),
            **overrides,
        )

    def test_known_reporter_keeps_full_cds_and_provenance_for_any_exact_prefix(self):
        for key, reporter in REPORTERS.items():
            for suffix in (reporter.sequence_after_start_rna,
                           reporter.full_sequence_after_start_rna[:45],
                           reporter.full_sequence_after_start_rna):
                with self.subTest(reporter=key, suffix_length=len(suffix)):
                    base = self.request(reporter=reporter)
                    switch = supplied_switch(base, suffix=suffix, linker="GU")
                    result = prepare_manually_annotated_reading(base, switch)
                    self.assertEqual(result.fixed_switch_sequence, switch)
                    self.assertEqual(result.reporter, replace(reporter, sequence_after_start_rna=suffix))
                    self.assertEqual(result.architecture.frame_linker, "GU")
                    self.assertEqual(result.normalized(), result)

    def test_normalizes_rna_chemistry_and_whitespace_without_redesign(self):
        base = self.request()
        switch = supplied_switch(base, linker="GU")
        text = " \n".join(switch.lower().replace("u", "t")[i:i + 10] for i in range(0, len(switch), 10))
        result = prepare_manually_annotated_reading(base, text)
        self.assertEqual(result.fixed_switch_sequence, switch)
        self.assertNotIn("T", result.fixed_switch_sequence)
        self.assertEqual(result.architecture.frame_linker, "GU")
        self.assertEqual(result.sequences.aptamer_dna, base.sequences.aptamer_dna)
        self.assertEqual(result.sequences.trigger_dna, base.sequences.trigger_dna)
        self.assertEqual(result.reporter.full_sequence_after_start_rna,
                         base.reporter.full_sequence_after_start_rna)

    def test_native_reading_rejects_incompatible_strand_chemistry(self):
        base = self.request(thermo=ThermoConditions(switch_material="dna", trigger_material="rna"))
        with self.assertRaises(ValueError):
            prepare_manually_annotated_reading(base, supplied_switch(base))

    def test_explicit_full_supplied_reporter_is_frozen_without_invented_source(self):
        base = self.request()
        suffix = "CCCGGGAUCUAA"
        switch = supplied_switch(base, suffix=suffix)
        result = prepare_manually_annotated_reading(base, switch, reporter_is_complete=True)
        self.assertEqual(result.fixed_switch_sequence, switch)
        self.assertEqual(result.reporter, ReporterContext(
            key="reading_reporter", label="Rapporteur fourni",
            sequence_after_start_rna=suffix, full_sequence_after_start_rna=suffix,
        ))
        self.assertIsNone(result.reporter.source_accession)
        self.assertIsNone(result.reporter.source_url)

    def test_mismatch_or_missing_full_reporter_requires_explicit_choice(self):
        for reporter in (REPORTERS["sfgfp"], ReporterContext()):
            with self.subTest(reporter=reporter), self.assertRaisesRegex(
                ValueError, "Sélectionnez le rapporteur correspondant"
            ):
                base = self.request(reporter=reporter)
                prepare_manually_annotated_reading(base, supplied_switch(base, suffix="CCCGGGAUCUAA"))

    def test_empty_suffix_and_invalid_bases_are_rejected(self):
        base = self.request()
        for sequence in (supplied_switch(base, suffix=""), supplied_switch(base) + "N"):
            with self.subTest(sequence=sequence), self.assertRaises(ValueError):
                prepare_manually_annotated_reading(base, sequence, reporter_is_complete=True)

    def test_reading_always_computes_metrics_without_exclusion_or_inverse_design(self):
        base = self.request(trials=99, inverse_extension=InverseExtensionSettings(),
                            optimization_options={"unknown_old_setting": "unused"})
        result = prepare_manually_annotated_reading(base, supplied_switch(base))
        self.assertEqual(result.trials, 1)
        self.assertIsNone(result.inverse_extension)
        self.assertNotIn("unknown_old_setting", result.optimization_options)
        self.assertFalse(result.exclude_stop_candidates)
        self.assertFalse(result.exclude_non_linear_candidates)
        self.assertEqual(result.thermo, base.thermo)
        self.assertEqual(result.scoring, base.scoring)
        self.assertEqual(base.trials, 99)
        self.assertIsNotNone(base.inverse_extension)

    def test_reading_preserves_whole_trigger_despite_a_previous_binding_window(self):
        base = self.request()
        base = replace(base, sequences=replace(base.sequences, binding_start=4, binding_length=12))
        result = prepare_manually_annotated_reading(base, supplied_switch(base))
        self.assertEqual(result.sequences.trigger_dna, base.sequences.trigger_dna)
        self.assertEqual(result.sequences.effective_trigger(), base.sequences.trigger_dna)
        self.assertEqual(result.sequences.binding_start, 0)
        self.assertEqual(result.sequences.binding_length, len(base.sequences.trigger_dna))
        # An architecture describing only the previous window is incompatible
        # with a fixed, whole-trigger analysis and must not crop it silently.
        cropped = replace(base, architecture=default_architecture_for_trigger(12))
        with self.assertRaisesRegex(ValueError, "longueur du trigger"):
            prepare_manually_annotated_reading(cropped, supplied_switch(cropped))

    def test_fixed_architecture_motifs_and_custom_linker_are_never_replaced(self):
        base = self.request()
        layout = architecture_layout(base.architecture)
        switch = supplied_switch(base)
        for start in (0, layout["rbs_start"], layout["aug_start"]):
            with self.subTest(position=start), self.assertRaisesRegex(ValueError, "motif"):
                changed = switch[:start] + "C" + switch[start + 1:]
                prepare_manually_annotated_reading(base, changed)
        base = replace(base, architecture=replace(base.architecture, frame_linker="GX"))
        result = prepare_manually_annotated_reading(base, supplied_switch(base, linker="GU"))
        self.assertEqual(result.architecture.frame_linker, "GU")
        with self.assertRaisesRegex(ValueError, "linker"):
            prepare_manually_annotated_reading(base, supplied_switch(base, linker="AU"))

    def test_standard_30nt_architecture_keeps_its_full_linker(self):
        base = replace(self.request(), sequences=SequenceInput("AACCGG", "AACC" * 7 + "AA"),
                       architecture=default_architecture_for_trigger(30))
        result = prepare_manually_annotated_reading(base, supplied_switch(base))
        self.assertEqual(result.architecture.frame_linker,
                         base.architecture.frame_linker.replace("T", "U"))
        self.assertEqual(architecture_layout(result.architecture)["reporter_start"],
                         architecture_layout(base.architecture)["reporter_start"])

    def test_saved_switch_roundtrip_keeps_exact_input_and_full_context(self):
        for supplied_complete in (False, True):
            with self.subTest(supplied_complete=supplied_complete), tempfile.TemporaryDirectory() as folder:
                base = self.request(output_dir=Path(folder).resolve() / "runs")
                suffix = "CCCGGGAUCUAA" if supplied_complete else None
                switch = supplied_switch(base, suffix=suffix, linker="GU")
                request = prepare_manually_annotated_reading(base, switch, reporter_is_complete=supplied_complete)
                result = design_result_fixture(request)
                result.optimization = {
                    "algorithm": "fixed_sequence_evaluation", "objective": "supplied_sequences",
                    "sequence_optimization": False,
                    "evaluations": 1, "trials": [],
                }
                candidate = result.candidates[0]
                candidate.switch_seq = switch
                snapshot = saved_switch_library.candidate_snapshot(result, candidate, "Lecture")
                saved = saved_switch_library.save_switch(Path(folder).resolve() / "saved", "Lecture", snapshot)
                loaded = saved_switch_library.load_switch(saved.path)
                context = build_saved_switch_context(loaded.snapshot, {})
                self.assertEqual(loaded.snapshot["request"]["fixed_switch_sequence"], switch)
                self.assertEqual(context["designed_switch_sequence"], switch)
                self.assertEqual(context["trigger_sequence"], candidate.trigger_seq)
                self.assertTrue(loaded.snapshot["sequences"]["reporter_is_complete"])
                self.assertEqual(loaded.snapshot["request"]["engine"], "native-hybrid")
                self.assertEqual(context["material"], request.thermo.material)
                if supplied_complete:
                    self.assertEqual(context["full_switch_sequence"], switch)
                else:
                    self.assertTrue(context["full_switch_sequence"].startswith(switch))
