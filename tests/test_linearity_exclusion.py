"""Native ON screening precedes candidate physics and keeps its exact witness."""
from dataclasses import replace
import random
from types import SimpleNamespace
import unittest

from aptaswitch_core.models import SequenceInput
from aptaswitch_core.scoring import sort_candidates
from aptaswitch_core.thermodynamics.design_constraints import build_design_problem
from aptaswitch_core.thermodynamics.design_engine import evaluate_candidate
from tests.test_thermodynamics_design_engine import AnalyticalBackend, request


class LinearityRegionTests(unittest.TestCase):
    def test_screen_includes_entire_loop_prefix_and_linker_but_not_reporter(self):
        for prefix, leader, linker in ((0, '', ''), (3, 'GGG', 'GCU'), (5, 'GAGAGA', 'ACC')):
            options = request()
            options = replace(options, architecture=replace(options.architecture,
                              rbs_prefix_length=prefix, t7_leader=leader, frame_linker=linker))
            problem = build_design_problem(options)
            assignment = problem.generate(random.Random(17))
            backend = AnalyticalBackend(options, bad_on=True)
            candidate = evaluate_candidate(options, problem, assignment, backend, trial=1)
            offset = len(assignment.sequence('Trigger'))
            self.assertEqual(backend.calls, [('screen', offset+problem.layout['rbs_start']-prefix,
                                             offset+problem.layout['reporter_start'])])
            self.assertTrue(candidate.excluded)
            self.assertEqual(problem.layout['rbs_linker_end'], problem.layout['reporter_start'])

    def test_failed_all_minimum_screen_requires_a_structural_witness(self):
        options = request()
        problem = build_design_problem(options)
        backend = AnalyticalBackend(options)
        backend.screen_unpaired_region = lambda *args, **kwargs: SimpleNamespace(all_minima_unpaired=False, witness=None)
        with self.assertRaisesRegex(ValueError, 'witness'):
            evaluate_candidate(options, problem, problem.generate(random.Random(17)), backend, trial=1)


class LinearityEngineTests(unittest.TestCase):
    def evaluate(self, options, *, bad_on=False, trial=1):
        problem = build_design_problem(options)
        backend = AnalyticalBackend(options, bad_on=bad_on)
        assignment = problem.generate(random.Random(17))
        candidate = evaluate_candidate(options, problem, assignment, backend, trial=trial)
        return candidate, backend

    def test_default_and_serialized_policy_survive_normalization(self):
        self.assertTrue(request().exclude_non_linear_candidates)
        for enabled in (False, True):
            options = request(exclude_non_linear_candidates=enabled)
            self.assertEqual(options.normalized().exclude_non_linear_candidates, enabled)
            self.assertEqual(options.to_dict()['exclude_non_linear_candidates'], enabled)

    def test_failure_keeps_on_witness_and_skips_defects_tubes_score_and_other_minima(self):
        candidate, backend = self.evaluate(request(), bad_on=True)
        self.assertEqual(len(backend.calls), 1)
        self.assertEqual(backend.calls[0][0], 'screen')
        self.assertTrue(candidate.excluded)
        self.assertTrue(candidate.exclusion_reason.startswith('RBS–linker non linéaire en ON'))
        self.assertIn('+', candidate.structure_on)
        self.assertEqual(candidate.delta_g_on, -25.)
        self.assertEqual(candidate.structure_off, '')
        for name in ('score', 'defect', 'on_yield_pct', 'leak_pct', 'delta_g_off',
                     'ddg_activation', 'delta_g_rbs_linker', 'has_bad_rbs_linker'):
            self.assertIsNone(getattr(candidate, name), name)

    def test_valid_candidate_uses_one_on_minimum_and_sorts_before_excluded(self):
        excluded, _ = self.evaluate(request(), bad_on=True)
        valid, backend = self.evaluate(request(), trial=2)
        on_calls = [call for call in backend.calls if call == ('minimum', ('Trigger', 'Switch'))]
        self.assertEqual(len(on_calls), 1)
        self.assertFalse(valid.excluded)
        self.assertEqual([row.trial for row in sort_candidates([excluded, valid])], [2, 1])

    def test_disabled_filter_does_not_request_screening_and_continues_scoring(self):
        candidate, backend = self.evaluate(request(exclude_non_linear_candidates=False), bad_on=True)
        self.assertFalse(candidate.excluded)
        self.assertIsNotNone(candidate.score)
        self.assertFalse(any(call[0] == 'screen' for call in backend.calls))
        self.assertTrue(any(call[0] == 'analyze' for call in backend.calls))

    def test_stop_free_generation_precedes_on_screen_and_opt_out_keeps_stop_diagnostic(self):
        options = replace(request(), exclude_stop_candidates=True,
                          sequences=SequenceInput('AACCGGTTAACC', 'GCTGATTCTCAATGGCTCGCCCGATATGCA'))
        filtered, backend = self.evaluate(options, bad_on=True)
        self.assertFalse(filtered.has_stop_codon)
        self.assertTrue(filtered.exclusion_reason.startswith('RBS–linker'))
        self.assertEqual(len(backend.calls), 1)
        failed, backend = self.evaluate(replace(options, exclude_stop_candidates=False), bad_on=True)
        self.assertTrue(failed.has_stop_codon)
        self.assertIsNotNone(failed.first_stop)
        self.assertEqual(len(backend.calls), 1)
        self.assertTrue(failed.exclusion_reason.startswith('RBS–linker'))


if __name__ == '__main__':
    unittest.main()
