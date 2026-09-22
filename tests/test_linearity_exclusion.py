"""The optional ON RBS–linker screen rejects pairing before downstream analyses."""

from __future__ import annotations

import unittest
from contextlib import contextmanager
from dataclasses import replace
from types import SimpleNamespace
from unittest.mock import Mock, patch

from aptaswitch_core.architecture import architecture_layout, default_architecture_for_trigger
from aptaswitch_core.models import RunRequest, SequenceInput
from aptaswitch_core.nupack_engine import _first_non_linear_on_mfe, run_nupack_engine


TRIGGER = "GTCCAGGCTGGTATAATTAGATCC"


def request(**kwargs):
    return RunRequest(
        sequences=SequenceInput(aptamer_dna="ACGTACGT", trigger_dna=TRIGGER),
        architecture=default_architecture_for_trigger(len(TRIGGER)),
        trials=kwargs.pop("trials", 1), nupack_path="/test-only/nupack", **kwargs,
    )


def switch_sequence(options):
    layout = architecture_layout(options.architecture)
    size = layout["reporter_start"] + len(options.reporter.sequence_after_start_rna)
    sequence = list("C" * size)
    sequence[layout["aug_start"]:layout["aug_start"] + 3] = "AUG"
    return "".join(sequence)


def on_structure(sequence, pairs=()):
    """Test pairs are indices in trigger + switch, excluding the strand break."""
    symbols = list("." * (len(TRIGGER) + len(sequence)))
    for left, right in pairs:
        symbols[left], symbols[right] = "(", ")"
    return "".join(symbols[:len(TRIGGER)]) + "+" + "".join(symbols[len(TRIGGER):])


def record(structure, energy=-10.0):
    return SimpleNamespace(structure=structure, energy=energy)


@contextmanager
def fake_design(options, structures, sequences=None, *, include_defect=True):
    sequences = sequences or [switch_sequence(options)] * len(structures)
    returned_on = iter(structures)

    def mfe(*, strands, model):
        if len(strands) == 2:
            return next(returned_on)
        return [record("." * (len(strands[0]) - 1), -1.0)]

    nupack = SimpleNamespace(Model=Mock(return_value=object()), mfe=Mock(side_effect=mfe))
    outputs = []
    for sequence in sequences:
        output = SimpleNamespace(
            to_analysis=lambda strand, seq=sequence: seq if strand == "switch" else TRIGGER,
        )
        if include_defect:
            output.defects = SimpleNamespace(ensemble_defect=.01)
        outputs.append([output])
    spec = Mock()
    spec.run.side_effect = outputs
    design = {
        "spec": spec, "switch": "switch", "trigger": "trigger",
        "layout": architecture_layout(options.architecture),
    }
    with (
        patch("aptaswitch_core.nupack_engine.nupack_import_context") as context,
        patch("aptaswitch_core.nupack_engine._build_design_spec", return_value=design),
        patch("aptaswitch_core.nupack_engine._complex_pct", return_value=20.0) as tubes,
        patch("aptaswitch_core.nupack_engine.compute_selection_score", return_value=-12.0) as score,
    ):
        context.return_value.__enter__.return_value = nupack
        yield nupack, tubes, score, spec


class LinearityRegionTests(unittest.TestCase):
    def setUp(self):
        from aptaswitch_core.localization import language_context
        self.enterContext(language_context("fr"))
        self.options = request()
        self.sequence = switch_sequence(self.options)
        self.layout = architecture_layout(self.options.architecture)
        self.start = self.layout["rbs_linker_start"]
        self.end = self.layout["rbs_linker_end"]

    def screen(self, records, **kwargs):
        return _first_non_linear_on_mfe(records, **{
            "trigger_length": len(TRIGGER), "switch_length": len(self.sequence),
            "start": self.start, "end": self.end, **kwargs,
        })

    def test_exact_loop_prefix_and_last_linker_base_are_included_reporter_is_excluded(self):
        offset = len(TRIGGER)
        cases = (
            ((offset + self.start, offset + self.end + 5), True),
            ((offset + self.end - 1, offset + self.end + 5), True),
            ((offset + self.start - 5, offset + self.start - 1), False),
            ((offset + self.end, offset + self.end + 5), False),
        )
        self.assertEqual(self.start, self.layout["rbs_start"] - self.options.architecture.rbs_prefix_length)
        self.assertEqual(self.end, self.layout["reporter_start"])
        for pair, excluded in cases:
            with self.subTest(pair=pair):
                mfe = record(on_structure(self.sequence, [pair]))
                self.assertEqual(self.screen([mfe]) is mfe, excluded)

    def test_trigger_pairs_into_region_fail_but_binding_upstream_is_allowed(self):
        for local_index, excluded in ((self.start, True), (self.start - 1, False)):
            with self.subTest(local_index=local_index):
                mfe = record(on_structure(self.sequence, [(0, len(TRIGGER) + local_index)]))
                self.assertEqual(self.screen([mfe]) is mfe, excluded)

    def test_every_returned_mfe_must_pass_and_first_failing_record_is_returned(self):
        good = record(on_structure(self.sequence))
        bad = record(on_structure(self.sequence, [(0, len(TRIGGER) + self.start)]))
        self.assertIsNone(self.screen([good, good]))
        self.assertIs(self.screen([good, bad, bad]), bad)

    def test_invalid_structure_shape_syntax_and_interval_fail_explicitly(self):
        good = on_structure(self.sequence)
        malformed = (
            good.replace("+", ""), good.replace("+", "++"),
            "." + good, "(" + good[1:], ")" + good[1:], "x" + good[1:],
        )
        for structure in malformed:
            with self.subTest(structure=structure), self.assertRaisesRegex(ValueError, "Criblage RBS–linker"):
                self.screen([record(structure)])
        for kwargs in ({"start": -1}, {"end": len(self.sequence) + 1}, {"end": self.start}):
            with self.subTest(kwargs=kwargs), self.assertRaisesRegex(ValueError, "bornes"):
                self.screen([record(good)], **kwargs)
        with self.assertRaisesRegex(ValueError, "aucune structure"):
            self.screen([])


class LinearityEngineTests(unittest.TestCase):
    def setUp(self):
        from aptaswitch_core.localization import language_context
        self.enterContext(language_context("fr"))

    def test_default_and_serialized_policy_survive_normalization(self):
        self.assertTrue(request().exclude_non_linear_candidates)
        for enabled in (False, True):
            options = request(exclude_non_linear_candidates=enabled)
            self.assertEqual(options.normalized().exclude_non_linear_candidates, enabled)
            self.assertEqual(options.to_dict()["exclude_non_linear_candidates"], enabled)

    def test_rejects_after_one_on_call_without_defect_other_mfe_tubes_or_score(self):
        options = request()
        sequence = switch_sequence(options)
        start = architecture_layout(options.architecture)["rbs_linker_start"]
        good = record(on_structure(sequence))
        bad = record(on_structure(sequence, [(0, len(TRIGGER) + start)]), -9.87654321)
        progress = []
        with fake_design(options, [[good, bad]], include_defect=False) as (nupack, tubes, score, _):
            result = run_nupack_engine(options, progress_callback=lambda *args: progress.append(args))
        self.assertEqual(nupack.mfe.call_count, 1)
        self.assertEqual(len(nupack.mfe.call_args.kwargs["strands"]), 2)
        tubes.assert_not_called()
        score.assert_not_called()
        candidate = result.candidates[0]
        self.assertTrue(candidate.excluded)
        self.assertFalse(candidate.has_stop_codon)
        self.assertTrue(candidate.exclusion_reason.startswith("RBS–linker non linéaire en ON"))
        self.assertEqual(candidate.structure_on, bad.structure)
        self.assertEqual(candidate.delta_g_on, bad.energy)
        self.assertEqual(candidate.structure_off, "")
        for field in ("score", "defect", "on_yield_pct", "leak_pct", "delta_g_off", "ddg_activation", "delta_g_rbs_linker", "has_bad_rbs_linker"):
            self.assertIsNone(getattr(candidate, field), field)
        self.assertEqual(progress[-1][0], 1)
        self.assertIn("analyses en tube ignorées", progress[-1][2])

    def test_valid_candidate_reuses_on_mfe_and_sorts_before_excluded(self):
        options = request(trials=2)
        sequence = switch_sequence(options)
        start = architecture_layout(options.architecture)["rbs_linker_start"]
        bad = record(on_structure(sequence, [(0, len(TRIGGER) + start)]))
        good = record(on_structure(sequence))
        with fake_design(options, [[bad], [good]]) as (nupack, tubes, score, _):
            result = run_nupack_engine(options)
        self.assertEqual(nupack.mfe.call_count, 5)  # One excluded + four retained, including ON once.
        self.assertEqual(sum(len(call.kwargs["strands"]) == 2 for call in nupack.mfe.call_args_list), 2)
        self.assertEqual(tubes.call_count, 2)
        score.assert_called_once()
        self.assertEqual([item.trial for item in result.candidates], [2, 1])
        self.assertEqual(result.candidates[0].score, -12.0)
        self.assertEqual(result.candidates[0].structure_on, good.structure)

    def test_disabled_filter_scores_paired_candidate_without_warning(self):
        options = request(exclude_non_linear_candidates=False)
        sequence = switch_sequence(options)
        start = architecture_layout(options.architecture)["rbs_linker_start"]
        bad = record(on_structure(sequence, [(0, len(TRIGGER) + start)]))
        with fake_design(options, [[bad]]) as (nupack, tubes, score, _):
            result = run_nupack_engine(options)
        self.assertEqual(nupack.mfe.call_count, 4)
        self.assertEqual(tubes.call_count, 2)
        score.assert_called_once()
        self.assertFalse(result.candidates[0].excluded)
        self.assertEqual(result.warnings, [])

    def test_stop_screen_runs_first_and_disabled_stop_retains_accurate_exclusion_reason(self):
        options = request()
        layout = architecture_layout(options.architecture)
        sequence = switch_sequence(options)
        start = layout["aug_start"] + 3
        sequence = sequence[:start] + "UAA" + sequence[start + 3:]
        bad = record(on_structure(sequence, [(0, len(TRIGGER) + layout["rbs_linker_start"])]))
        with fake_design(options, [[bad]], [sequence], include_defect=False) as (nupack, tubes, score, _):
            result = run_nupack_engine(options)
        nupack.mfe.assert_not_called()
        tubes.assert_not_called()
        score.assert_not_called()
        self.assertTrue(result.candidates[0].exclusion_reason.startswith("Codon STOP"))
        with fake_design(options, [[bad]], [sequence], include_defect=False) as (nupack, tubes, score, _):
            result = run_nupack_engine(replace(options, exclude_stop_candidates=False))
        self.assertEqual(nupack.mfe.call_count, 1)
        self.assertTrue(result.candidates[0].has_stop_codon)
        self.assertEqual(result.candidates[0].first_stop, 2)
        self.assertTrue(result.candidates[0].exclusion_reason.startswith("RBS–linker"))


if __name__ == "__main__":
    unittest.main()
