"""Inverse extension conservation stays separate from switch performance."""
from __future__ import annotations

import copy
from dataclasses import asdict, replace
import os
from pathlib import Path
from types import SimpleNamespace
import unittest
from unittest.mock import Mock, patch

import numpy as np

from aptaswitch_core import inverse_extension as inverse
from aptaswitch_core import extension_thermodynamics
from aptaswitch_core.models import (
    DesignCandidate, InverseExtensionSettings, RunRequest, RunResult, SequenceInput,
    SwitchArchitecture, ThermoConditions,
)
from aptaswitch_core.trigger_extension import ComplexAnalysisResult, parse_pairs


def candidate(trial, trigger="AAGT", score=100.0, excluded=False):
    return DesignCandidate(trial=trial, switch_seq="ACGU", trigger_seq=trigger, score=score,
        defect=.01, on_yield_pct=90, leak_pct=2, delta_g_off=-10, delta_g_on=-15,
        ddg_activation=-4, delta_g_rbs_linker=0, has_bad_rbs_linker=False,
        has_stop_codon=False, excluded=excluded)


def result(candidates=None, *, side="5prime", target=4, top=20):
    settings = InverseExtensionSettings(target_length=target, extension_side=side, top_candidates=top,
        thermo=ThermoConditions(switch_material="dna", trigger_material="dna", temperature_c=29.5, sodium_m=.15, magnesium_m=0))
    request = RunRequest(SequenceInput("ACGT", "GT"), SwitchArchitecture(12, 18), inverse_extension=settings)
    return RunResult("inverse-test", request, candidates if candidates is not None else [candidate(1)], Path("unused"), "start", "end")


def analysis(structure):
    size = len(structure.replace("+", ""))
    matrix = np.eye(size)
    for left, right in parse_pairs(structure):
        matrix[left, left] = matrix[right, right] = .1
        matrix[left, right] = matrix[right, left] = .9
    return ComplexAnalysisResult(structure, -3.5, -4.0, .6, .7, matrix)


class InverseAnalysisTests(unittest.TestCase):
    def setUp(self):
        from test_thermodynamics_context_analysis import FixedWeightBackend
        self.backend = FixedWeightBackend()
        self.enterContext(patch("aptaswitch_core.thermodynamics.setup.get_backend", return_value=self.backend))
        self.analyze = Mock(side_effect=self.native)
        self.fake = SimpleNamespace(metadata={"thermodynamics": asdict(self.backend.provenance)}, analyze=self.analyze)
        self.context = self.enterContext(patch.object(inverse, "open_extension_model"))
        self.context.return_value.__enter__.return_value = self.fake

    @staticmethod
    def native(aptamer, trigger, target=None):
        return analysis(target or "..((+))")

    def test_top_switches_reference_once_duplicate_triggers_once_scores_untouched(self):
        values = [candidate(4, "TTGT", 50), candidate(2, "AAGT", 90), candidate(1, "CCGT", 99, True),
                  candidate(3, "AAGT", 80), candidate(5, "CCGT", 70)]
        run = result(values, top=3)
        original = [(id(c), c.score, c.trigger_seq, c.switch_seq) for c in values]
        inverse.analyze_inverse_extensions(run, "/fake")
        stage = run.inverse_extension_analysis
        self.assertEqual(stage["status"], "completed")
        self.assertEqual(stage["selected_trials"], [2, 3, 5])
        self.assertEqual(stage["analyzed_count"], 3)
        self.assertEqual(stage["unique_trigger_count"], 2)
        self.assertEqual(self.analyze.call_count, 3)
        self.assertEqual(stage["model"]["temperature_c"], 29.5)
        self.assertEqual(stage["model"]["sodium_m"], .15)
        self.assertEqual([(id(c), c.score, c.trigger_seq, c.switch_seq) for c in run.candidates], original)
        self.assertEqual(run.candidates[0].inverse_extension_metrics, {})
        self.assertEqual(run.candidates[2].inverse_extension_metrics, {})
        self.assertEqual(stage["reference"]["trigger_sequence"], "GT")
        self.assertEqual(stage["thermodynamics"], asdict(self.backend.provenance))

        self.assertEqual([r["trial"] for r in stage["candidates"]], [2, 3, 5])
        first, duplicate = values[1].inverse_extension_metrics, values[3].inverse_extension_metrics
        self.assertEqual(first["final_score"], duplicate["final_score"])
        self.assertEqual(first["thermodynamics"], stage["thermodynamics"])
        self.assertEqual(first["missing_reference_mfe_pairs"], 0)
        self.assertEqual(first["core_pair_probability_rmse"], 0)
        self.assertEqual(first["reference_pairs"][0]["reference_probability"], .9)
        self.assertEqual(first["reference_pairs"][0]["candidate_probability"], .9)
        first["reference_pairs"][0]["candidate_probability"] = 0
        self.assertEqual(duplicate["reference_pairs"][0]["candidate_probability"], .9)

    def test_extension_mapping_for_both_ends_preserves_dna(self):
        for side, target, trigger, prefix, suffix, structure in [
            ("5prime", 4, "AAGT", "AA", "", "..((+..))"),
            ("3prime", 4, "GTAA", "", "AA", "..((+)).."),
            ("both", 5, "AAGTC", "AA", "C", "..((+..))."),
        ]:
            with self.subTest(side=side):
                run = result([candidate(1, trigger)], side=side, target=target)
                inverse.analyze_inverse_extensions(run, "/fake")
                self.assertEqual(run.inverse_extension_analysis["status"], "completed")
                metrics = run.candidates[0].inverse_extension_metrics
                self.assertEqual((metrics["extension_5prime"], metrics["extension_3prime"]), (prefix, suffix))
                self.assertEqual(metrics["mfe_structure"], structure)
                self.assertEqual(metrics["extended_trigger"], trigger)
                self.assertEqual(metrics["reference_mfe_structure"], "..((+))")
                self.assertEqual(metrics["missing_reference_mfe_pairs"], 0)
                self.assertEqual(self.analyze.call_args.args[0], "ACGT")

    def test_original_core_and_target_length_are_strict(self):
        for trigger in ("AAAA", "AGT"):
            with self.subTest(trigger=trigger):
                run = result([candidate(1, trigger)])
                inverse.analyze_inverse_extensions(run, "/fake")
                self.assertEqual(run.inverse_extension_analysis["status"], "failed")
                self.assertTrue(run.inverse_extension_analysis["error"])
                self.assertEqual(run.candidates[0].inverse_extension_metrics, {})
        self.context.assert_not_called()

    def test_cancel_before_native_call_and_after_reference_are_explicit(self):
        run = result()
        inverse.analyze_inverse_extensions(run, "/fake", cancel_check=lambda: True)
        self.assertEqual(run.inverse_extension_analysis["status"], "cancelled")
        self.analyze.assert_not_called()
        inverse.analyze_inverse_extensions(run, "/fake", cancel_check=lambda: self.analyze.call_count > 0)
        self.assertEqual(run.inverse_extension_analysis["status"], "cancelled")
        self.assertEqual(run.inverse_extension_analysis["analyzed_count"], 0)
        self.assertEqual(run.candidates[0].inverse_extension_metrics, {})

    def test_failed_second_trigger_retains_only_first_completed_analysis(self):
        run = result([candidate(1), candidate(2, "CCGT", 90)])
        self.analyze.side_effect = [analysis("..((+))"), analysis("..((+..))"), RuntimeError("native failure")]
        inverse.analyze_inverse_extensions(run, "/fake")
        self.assertEqual(run.inverse_extension_analysis["status"], "failed")
        self.assertEqual(run.inverse_extension_analysis["error"], "native failure")
        self.assertEqual(run.inverse_extension_analysis["analyzed_count"], 1)
        self.assertEqual(run.candidates[0].inverse_extension_metrics["status"], "completed")
        self.assertEqual(run.candidates[1].inverse_extension_metrics, {})
        self.assertEqual(run.status, "completed")

    def test_cancellation_between_triggers_retains_completed_metrics(self):
        run = result([candidate(1), candidate(2, "CCGT", 90)])
        cancel = False
        def progress(done, total, message):
            nonlocal cancel
            if done == 2:
                cancel = True
        inverse.analyze_inverse_extensions(run, "/fake", progress_callback=progress, cancel_check=lambda: cancel)
        self.assertEqual(run.inverse_extension_analysis["status"], "cancelled")
        self.assertEqual(run.inverse_extension_analysis["analyzed_count"], 1)
        self.assertEqual(self.analyze.call_count, 2)
        self.assertEqual(run.candidates[0].inverse_extension_metrics["status"], "completed")
        self.assertEqual(run.candidates[1].inverse_extension_metrics, {})

    def test_invalid_native_probabilities_fail_without_fabricating_metrics(self):
        invalid = analysis("..((+..))")
        invalid.pair_matrix[0, 0] = float("nan")
        self.analyze.side_effect = [analysis("..((+))"), invalid]
        run = result()
        inverse.analyze_inverse_extensions(run, "/fake")
        self.assertEqual(run.inverse_extension_analysis["status"], "failed")
        self.assertEqual(run.candidates[0].inverse_extension_metrics, {})

    def test_no_eligible_switches_does_not_start_a_reference_calculation(self):
        run = result([candidate(1, excluded=True)])
        inverse.analyze_inverse_extensions(run, "/fake")
        self.assertEqual(run.inverse_extension_analysis["status"], "completed")
        self.assertEqual(run.inverse_extension_analysis["selected_count"], 0)
        self.analyze.assert_not_called()


class InverseRankingTests(unittest.TestCase):
    def run_with_metrics(self, scores=(100, 90, 80), impacts=(10, 0, 20)):
        run = result([candidate(i + 1, score=value) for i, value in enumerate(scores)])
        for value, impact in zip(run.candidates, impacts):
            value.inverse_extension_metrics = {"status": "completed", "final_score": impact}
        return run

    def test_weights_select_performance_or_conservation_and_do_not_mutate(self):
        run = self.run_with_metrics()
        before = copy.deepcopy(run.to_dict())
        self.assertEqual([row["trial"] for row in inverse.ranked_inverse_candidates(run, 1)], [1, 2, 3])
        self.assertEqual([row["trial"] for row in inverse.ranked_inverse_candidates(run, 0)], [2, 1, 3])
        rows = inverse.ranked_inverse_candidates(run)
        self.assertEqual([row["compromise_score"] for row in rows], [.25, .25, 1])
        self.assertEqual([row["pareto"] for row in rows], [True, True, False])
        self.assertEqual(run.to_dict(), before)

    def test_ties_use_average_ranks_and_single_candidate_has_zero_rank_cost(self):
        rows = inverse.ranked_inverse_candidates(self.run_with_metrics((100, 100, 80), (10, 10, 20)))
        self.assertEqual([row["performance_rank"] for row in rows], [1.5, 1.5, 3])
        self.assertEqual([row["structure_rank"] for row in rows], [1.5, 1.5, 3])
        self.assertEqual([row["compromise_score"] for row in rows], [.25, .25, 1])
        row = inverse.ranked_inverse_candidates(self.run_with_metrics((100,), (10,)))[0]
        self.assertEqual(row["compromise_score"], 0)

    def test_bounded_public_scores_keep_raw_performance_ranks_and_pareto(self):
        run = self.run_with_metrics((99.99, 99.99, 50), (10, 10, 20))
        run.request = replace(run.request, score_scale="bounded_v1")
        run.candidates[0].score_raw = 10000.01
        run.candidates[1].score_raw = 10000.02
        run.candidates[2].score_raw = 0
        rows = inverse.ranked_inverse_candidates(run, 1)
        self.assertEqual([row["trial"] for row in rows], [2, 1, 3])
        self.assertEqual([row["performance_rank"] for row in rows], [1, 2, 3])
        self.assertEqual([row["pareto"] for row in rows], [True, False, False])
        self.assertTrue(all(0 <= row["switch_score"] <= 100 for row in rows))
        self.assertEqual([row["switch_score_raw"] for row in rows], [10000.02, 10000.01, 0])


    def test_unanalyzed_excluded_or_nonfinite_values_are_not_ranked(self):
        run = self.run_with_metrics()
        run.candidates[0].inverse_extension_metrics = {}
        run.candidates[1].excluded = True
        run.candidates[2].inverse_extension_metrics["final_score"] = float("nan")
        self.assertEqual(inverse.ranked_inverse_candidates(run), [])
        for weight in (-.1, 1.1, float("nan"), float("inf"), True, None):
            with self.subTest(weight=weight), self.assertRaises(ValueError):
                inverse.ranked_inverse_candidates(run, weight)


@unittest.skipUnless(os.environ.get("APTASWITCH_HYBRID_NATIVE_LIBRARY"), "Native hybrid library required")
class RealInverseNativeTests(unittest.TestCase):
    def test_real_original_reference_and_extended_complex(self):
        run = result([candidate(1, "AACACGT"), candidate(2, "AACACGT", 90)], target=7)
        run.request = RunRequest(SequenceInput("ACGTACGT", "CACGT"), SwitchArchitecture(12, 18),
            inverse_extension=InverseExtensionSettings(target_length=7, top_candidates=2))
        inverse.analyze_inverse_extensions(run)
        stage = run.inverse_extension_analysis
        self.assertEqual(stage["status"], "completed", stage)
        self.assertEqual(stage["reference"]["trigger_sequence"], "CACGT")
        self.assertEqual(stage["analyzed_count"], 2)
        self.assertEqual(stage["unique_trigger_count"], 1)
        self.assertEqual([len(s) for s in run.candidates[0].inverse_extension_metrics["mfe_structure"].split("+")], [8, 7])
        self.assertTrue(stage["thermodynamics"]["parameter_digest"])


if __name__ == "__main__":
    unittest.main()
