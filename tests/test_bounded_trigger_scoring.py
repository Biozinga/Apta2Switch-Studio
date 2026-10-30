"""Bounded ranking and ensemble folding of free, inversely designed triggers."""

from __future__ import annotations

import unittest
from dataclasses import replace
from types import SimpleNamespace
from unittest.mock import patch

import numpy as np

from aptaswitch_core.models import InverseExtensionSettings, ScoringWeights, SequenceInput
from aptaswitch_core.scoring import compute_selection_score, selection_score_100, sort_candidates
from tests.test_thermodynamics_design_engine import AnalyticalBackend, request as native_request


def request(**kwargs):
    return replace(native_request(), **kwargs)


def inverse_request(**kwargs):
    return replace(request(**kwargs), sequences=SequenceInput("ACGTACGT", "ACGT" * 3),
                   inverse_extension=InverseExtensionSettings(target_length=30))


class BoundedScoreTests(unittest.TestCase):
    def test_fixed_mapping_is_bounded_and_monotone_without_run_relative_scaling(self):
        raw = [-1e308, -1e6, -199., -100., 0., 100., 199., 250., 1e6, 1e308]
        scores = [selection_score_100(value) for value in raw]
        self.assertTrue(all(0 <= value <= 100 for value in scores))
        self.assertEqual(scores, sorted(set(scores)))
        self.assertEqual(selection_score_100(0), 50)
        self.assertEqual(selection_score_100(100), 75)
        self.assertEqual(selection_score_100(-100), 25)
        for invalid in (float("nan"), float("inf"), -float("inf"), True):
            with self.subTest(value=invalid), self.assertRaises(ValueError):
                selection_score_100(invalid)

    def test_optional_free_trigger_penalty_has_the_expected_direction_and_weight(self):
        args = dict(on_yield_pct=90, leak_pct=2, functional_defect=.1, ddg_activation=-12,
                    delta_g_rbs_linker=-2, weights=ScoringWeights())
        baseline = compute_selection_score(**args)
        folded = compute_selection_score(**args, trigger_linearity_defect=.35)
        self.assertAlmostEqual(baseline - folded, 105)
        self.assertEqual(compute_selection_score(**args, trigger_linearity_defect=0), baseline)
        self.assertEqual(compute_selection_score(**{**args, "weights": replace(args["weights"], trigger_linearity=0)},
                                                 trigger_linearity_defect=.9), baseline)
        bounded = selection_score_100(compute_selection_score(**args, trigger_linearity_defect=.35))
        self.assertEqual(bounded, selection_score_100(folded))
        self.assertLess(bounded, selection_score_100(compute_selection_score(**args)))
        for invalid in (True, -.01, 1.01, float("nan"), float("inf")):
            with self.subTest(value=invalid), self.assertRaises(ValueError):
                compute_selection_score(**args, trigger_linearity_defect=invalid)
        with self.assertRaises(ValueError):
            compute_selection_score(**{**args, "on_yield_pct": float("inf")})


    def test_raw_score_retains_order_when_bounded_display_rounds_to_a_tie(self):
        first = SimpleNamespace(trial=1, score=selection_score_100(1e8), score_raw=1e8,
                                defect=.01, excluded=False, regional_defects={"balanced": .01})
        second = SimpleNamespace(trial=2, score=selection_score_100(1e8 + 100), score_raw=1e8 + 100,
                                 defect=.9, excluded=False, regional_defects={"balanced": .9})
        self.assertEqual(first.score, second.score)
        self.assertEqual(sort_candidates([first, second]), [second, first])

    def test_new_requests_persist_score_scale_and_validate_weight(self):
        options = request().normalized()
        self.assertEqual(options.score_scale, "bounded_v1")
        self.assertEqual(options.to_dict()["score_scale"], "bounded_v1")
        self.assertEqual(options.to_dict()["scoring"]["trigger_linearity"], 3.)
        with self.assertRaises(ValueError):
            replace(options, score_scale="raw").normalized()
        for invalid in (True, -1, float("nan"), float("inf"), "1"):
            with self.subTest(value=invalid), self.assertRaises(ValueError):
                replace(options, scoring=replace(options.scoring, trigger_linearity=invalid)).normalized()
        with self.assertRaises(ValueError):
            replace(options, score_scale="typo").normalized()


class TriggerLinearityTests(unittest.TestCase):
    def evaluate(self, options=None, *, matrix=None, bad_on=False, stop=False):
        import random
        from aptaswitch_core.thermodynamics.design_constraints import build_design_problem
        from aptaswitch_core.thermodynamics.design_engine import evaluate_candidate
        from aptaswitch_core.thermodynamics.types import PairProbabilities
        options = options or inverse_request()
        backend = AnalyticalBackend(options, bad_on=bad_on)
        original = backend.analyze
        def analyze(spec, conditions, *, need_pairs):
            result = original(spec, conditions, need_pairs=need_pairs)
            if need_pairs and tuple(s.name for s in spec.strands) == ("Trigger",) and matrix is not None:
                return replace(result, pairs=PairProbabilities(result.partition, matrix))
            return result
        backend.analyze = analyze
        problem = build_design_problem(options)
        assignment = problem.generate(random.Random(17))
        with patch.object(type(problem), "reading_frame_report", return_value={"has_premature_stop": stop, "first_stop": 2 if stop else None}):
            candidate = evaluate_candidate(options, problem, assignment, backend, trial=1)
        return candidate, backend

    def test_ensemble_mean_counts_both_ends_and_does_not_binary_score_mfe(self):
        matrix = np.eye(30)
        matrix[0, 0] = matrix[-1, -1] = .6
        matrix[0, -1] = matrix[-1, 0] = .4
        candidate, backend = self.evaluate(matrix=matrix)
        data = candidate.trigger_linearity
        self.assertAlmostEqual(data["normalized_defect"], .8 / 30)
        self.assertAlmostEqual(data["unpaired_fraction"], 1 - .8 / 30)
        self.assertAlmostEqual(data["expected_paired_nt"], .8)
        self.assertEqual(data["per_base_unpaired"], [.6] + [1.] * 28 + [.6])
        self.assertEqual(data["structure_mfe"], "." * 30)
        self.assertEqual([call for call in backend.calls if call == ("analyze", ("Trigger",), True)],
                         [("analyze", ("Trigger",), True)])

    def test_fixed_trigger_has_no_optional_probability_calculation(self):
        candidate, backend = self.evaluate(request())
        self.assertEqual(candidate.trigger_linearity, {})
        self.assertNotIn(("analyze", ("Trigger",), True), backend.calls)

    def test_sparse_invalid_and_mismatched_probability_matrices_fail(self):
        bad = [np.eye(29), .9 * np.eye(30), np.full((30, 30), float("nan")), np.eye(30, dtype=bool)]
        asymmetric = np.eye(30)
        asymmetric[0, 0], asymmetric[0, -1] = .6, .4
        bad.append(asymmetric)
        for matrix in bad:
            with self.subTest(matrix=matrix), self.assertRaises(ValueError):
                self.evaluate(matrix=matrix)

    def test_inverse_candidate_keeps_measurement_when_weight_is_zero(self):
        zero, _ = self.evaluate(inverse_request(scoring=ScoringWeights(trigger_linearity=0.)))
        weighted, _ = self.evaluate(inverse_request(scoring=ScoringWeights(trigger_linearity=1.)))
        self.assertEqual(zero.trigger_linearity, weighted.trigger_linearity)
        self.assertAlmostEqual(weighted.trigger_linearity["normalized_defect"], .2)
        self.assertAlmostEqual(zero.score_raw - weighted.score_raw, 20., places=4)
        self.assertEqual(weighted.score, selection_score_100(weighted.score_raw))

    def test_excluded_inverse_has_no_optional_trigger_measurement(self):
        for stop in (False, True):
            candidate, backend = self.evaluate(inverse_request(exclude_stop_candidates=True, exclude_non_linear_candidates=True),
                                                bad_on=True, stop=stop)
            self.assertTrue(candidate.excluded)
            self.assertEqual(candidate.trigger_linearity, {})
            self.assertIsNone(candidate.score_raw)
            self.assertNotIn(("analyze", ("Trigger",), True), backend.calls)


if __name__ == "__main__":
    unittest.main()
