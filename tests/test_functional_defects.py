"""Functional OFF lock/accessibility and overlapping ON opening defects."""
from dataclasses import asdict, replace
import math
from types import SimpleNamespace
from unittest import TestCase
from unittest.mock import Mock, call, patch

import numpy as np

from aptaswitch_core.architecture import architecture_layout, default_architecture_for_trigger, target_switch_structures
from aptaswitch_core.defect_metrics import functional_defect_metrics
from aptaswitch_core.models import FunctionalDefectWeights, RunRequest, ScoringWeights, SequenceInput
from aptaswitch_core.scoring import compute_selection_score, sort_candidates
from tests.defect_fixtures import target_matrix, disrupt_pairs, disrupt_unpaired


class FunctionalDefectTests(TestCase):
    def setUp(self):
        self.arch = default_architecture_for_trigger(30)
        self.reporter_length = 6
        self.layout = architecture_layout(self.arch)
        self.off_target, self.on_target = target_switch_structures(self.arch, self.reporter_length)
        self.off, self.on = target_matrix(self.off_target), target_matrix(self.on_target)

    def metrics(self, weights=None):
        return functional_defect_metrics(self.arch, self.reporter_length, self.off, self.on, weights)

    def test_classic_regions_match_both_lock_arms_and_overlapping_on_regions(self):
        result = self.metrics()
        self.assertEqual(result['scheme'], 'functional_v2')
        self.assertEqual(result['balanced'], 0)
        self.assertEqual(result['normalized_weights'],
                         {'off_toehold': 4 / 15, 'off_stem': 4 / 15,
                          'on_rbs_aug': 4 / 15, 'on_rbs_linker': 1 / 5})
        regions, parts = result['regions'], result['off_stem_parts']
        self.assertEqual(regions['off_toehold']['indices_switch'], list(range(3, 15)))
        self.assertEqual(parts['before_aug']['indices_switch'], list(range(27, 33)) + list(range(44, 50)))
        self.assertEqual(parts['after_aug']['indices_switch'], list(range(15, 24)) + list(range(53, 62)))
        self.assertEqual(parts['before_aug']['pair_count'], 6)
        self.assertEqual(parts['after_aug']['pair_count'], 9)
        self.assertEqual(regions['off_stem']['nucleotide_count'], 30)
        self.assertEqual(regions['on_rbs_aug']['indices_switch'], list(range(33, 53)))
        self.assertEqual(regions['on_rbs_linker']['indices_switch'], list(range(33, 83)))
        self.assertNotIn('on_binding', result)

    def test_separate_stem_defects_aggregate_per_nucleotide_not_equal_stem_weight(self):
        disrupt_pairs(self.off, list(range(44, 50)), .2)
        disrupt_pairs(self.off, list(range(53, 62)), .8)
        result = self.metrics()
        self.assertAlmostEqual(result['off_stem_parts']['before_aug']['normalized_defect'], .2)
        self.assertAlmostEqual(result['off_stem_parts']['after_aug']['normalized_defect'], .8)
        self.assertAlmostEqual(result['off_stem_parts']['before_aug']['expected_incorrect_nt'], 12 * .2)
        self.assertAlmostEqual(result['off_stem'], (6 * .2 + 9 * .8) / 15)
        self.assertAlmostEqual(result['balanced'], result['off_stem'] * 4 / 15)

    def test_off_wrong_partner_counts_as_error_even_when_base_remains_paired(self):
        i, k = 44, 45
        j, l = int(np.argmax(self.off[i])), int(np.argmax(self.off[k]))
        for a, b in ((i, j), (k, l)):
            self.off[a, b] = self.off[b, a] = 0
        for a, b in ((i, l), (k, j)):
            self.off[a, b] = self.off[b, a] = 1
        result = self.metrics()
        self.assertAlmostEqual(result['off_stem_parts']['before_aug']['normalized_defect'], 4 / 12)
        self.assertAlmostEqual(result['off_stem'], 4 / 30)

    def test_off_aug_bulge_and_on_duplex_errors_do_not_enter_selected_metrics(self):
        # AUG and the opposing bulge are unpaired by design. The locking
        # criterion is their flanking stem's correctness, not AUG pairing.
        disrupt_unpaired(self.off, [24, 25, 50, 51], .8)
        offset = self.layout['len_trigger']
        disrupt_pairs(self.on, list(range(offset + 3, offset + 33)), .8)
        result = self.metrics()
        self.assertGreater(result['switch_structural']['off'], 0)
        self.assertGreater(result['switch_structural']['on'], 0)
        self.assertEqual(result['balanced'], 0)

    def test_on_post_aug_and_linker_errors_only_affect_extended_region(self):
        offset = self.layout['len_trigger']
        disrupt_unpaired(self.on, [offset + 53, offset + 54, offset + 81, offset + 82], .5)
        result = self.metrics()
        self.assertEqual(result['on_rbs_aug'], 0)
        self.assertAlmostEqual(result['on_rbs_linker'], 4 * .5 / 50)
        self.assertAlmostEqual(result['balanced'], 4 * .5 / 50 / 5)

    def test_on_aug_and_prefix_errors_are_counted_in_both_on_regions(self):
        offset = self.layout['len_trigger']
        disrupt_unpaired(self.on, [offset + 33, offset + 34, offset + 51, offset + 52], .6)
        result = self.metrics()
        self.assertAlmostEqual(result['on_rbs_aug'], 4 * .6 / 20)
        self.assertAlmostEqual(result['on_rbs_linker'], 4 * .6 / 50)
        self.assertAlmostEqual(result['on'], (4 * .6 / 20 + .75 * 4 * .6 / 50) / 1.75)

    def test_off_toehold_measures_accessibility_and_excludes_leader(self):
        disrupt_unpaired(self.off, list(range(3, 15)), .7)
        disrupt_unpaired(self.off, [0, 1], .9)
        result = self.metrics()
        self.assertAlmostEqual(result['off_toehold'], .7)
        self.assertEqual(result['off_stem'], 0)
        self.assertAlmostEqual(result['balanced'], .7 * 4 / 15)

    def test_custom_geometry_tracks_lengths_and_excludes_aug_spacer_from_lock(self):
        arch = replace(default_architecture_for_trigger(24), upper_stem_length=4,
                       upper_stem2_length=2, rbs_prefix_length=1, aug_spacer=3,
                       frame_linker='ACGU')
        layout = architecture_layout(arch)
        off, on = target_switch_structures(arch, 0)
        result = functional_defect_metrics(arch, 0, target_matrix(off), target_matrix(on))
        parts = result['off_stem_parts']
        self.assertEqual(parts['before_aug']['pair_count'], 4)
        self.assertEqual(parts['after_aug']['pair_count'], 11)
        spacer = set(range(layout['aug_start'] - 3, layout['aug_start']))
        self.assertFalse(spacer & set(result['regions']['off_stem']['indices_switch']))
        self.assertTrue(spacer <= set(result['regions']['on_rbs_aug']['indices_switch']))
        self.assertEqual(result['regions']['on_rbs_aug']['indices_switch'][-1], layout['aug_start'] + 2)
        self.assertEqual(result['regions']['on_rbs_linker']['indices_switch'][-1], layout['reporter_start'] - 1)

    def test_reporter_is_excluded_from_both_on_regions(self):
        offset = self.layout['len_trigger']
        reporter_start = self.layout['reporter_start']
        disrupt_unpaired(self.on, list(range(offset + reporter_start, len(self.on))), .9)
        self.assertEqual(self.metrics()['balanced'], 0)

    def test_functional_weights_change_aggregation_and_keep_inactive_diagnostics(self):
        disrupt_unpaired(self.off, list(range(3, 15)), .4)
        before = self.metrics()
        result = self.metrics(FunctionalDefectWeights(1, 0, 0, 0))
        self.assertAlmostEqual(result['balanced'], .4)
        self.assertAlmostEqual(result['off'], .4)
        self.assertIsNone(result['on'])
        self.assertEqual(result['regions'], before['regions'])
        self.assertEqual(result['off_stem_parts'], before['off_stem_parts'])
        self.assertEqual(self.metrics(FunctionalDefectWeights(0, 0, 1, 0))['balanced'], 0)
        self.assertIsNone(self.metrics(FunctionalDefectWeights(0, 0, 1, 0))['off'])

    def test_complete_probability_matrices_required_for_both_states(self):
        for bad in (None, self.off[:-1, :-1], self.off * .9, self.off.astype(bool)):
            with self.subTest(bad=type(bad)), self.assertRaises(ValueError):
                functional_defect_metrics(self.arch, 6, bad, self.on)
        with self.assertRaises(ValueError):
            functional_defect_metrics(self.arch, 6, self.off, self.on * .9)

    def test_native_candidate_dispatch_preserves_functional_weights_and_strand_order(self):
        import random
        from aptaswitch_core.models import ReporterContext
        from aptaswitch_core.thermodynamics.design_constraints import build_design_problem
        from aptaswitch_core.thermodynamics.design_engine import evaluate_candidate
        from test_thermodynamics_design_engine import AnalyticalBackend
        weights = FunctionalDefectWeights(0, 2, 1, 1)
        request = RunRequest(SequenceInput("AACCGGTT", "ACGT" * 7 + "AC"), self.arch,
                             reporter=ReporterContext("fixture", "Fixture", "ACGUAC"),
                             functional_defect_weights=weights, exclude_stop_candidates=False,
                             exclude_non_linear_candidates=False)
        backend = AnalyticalBackend(request)
        problem = build_design_problem(request)
        candidate = evaluate_candidate(request, problem, problem.generate(random.Random(17)), backend, trial=1)
        self.assertEqual(candidate.regional_defects, functional_defect_metrics(self.arch, 6, self.off, self.on, weights))
        pair_calls = [call[1] for call in backend.calls if call[0] == "analyze" and call[2]]
        self.assertIn(("Switch",), pair_calls)
        self.assertIn(("Trigger", "Switch"), pair_calls)



class FunctionalWeightAndScoringTests(TestCase):
    def test_weight_validation_and_overflow_safe_normalization(self):
        for value in (True, np.bool_(False), -1, float('nan'), float('inf'), None, '1', 1j, 10 ** 1000):
            with self.subTest(value=repr(value)), self.assertRaises(ValueError):
                FunctionalDefectWeights(off_toehold=value).fractions()
        with self.assertRaises(ValueError):
            FunctionalDefectWeights(0, 0, 0, 0).fractions()
        self.assertEqual(FunctionalDefectWeights(*(1e308 for _ in range(4))).fractions(),
                         {name: .25 for name in asdict(FunctionalDefectWeights())})

    def test_new_request_defaults_functional_and_preserves_weights(self):
        request = RunRequest(SequenceInput('ACGTACGT', 'ACGT' * 7 + 'AC'),
                             default_architecture_for_trigger(30),
                             functional_defect_weights=FunctionalDefectWeights(1, 2, 3, 4))
        self.assertEqual(request.normalized().defect_metric, 'functional')
        self.assertEqual(request.normalized().to_dict()['functional_defect_weights'],
                         asdict(request.functional_defect_weights))
        for invalid in (None, {}, FunctionalDefectWeights(0, 0, 0, 0)):
            with self.subTest(invalid=invalid), self.assertRaises(ValueError):
                replace(request, functional_defect_weights=invalid).normalized()
        for unsupported in ('regional', 'global', 'unknown'):
            with self.subTest(metric=unsupported), self.assertRaises(ValueError):
                replace(request, defect_metric=unsupported).normalized()

    def test_functional_penalty_uses_requested_weight(self):
        args = dict(on_yield_pct=0, leak_pct=0, functional_defect=.2,
                    ddg_activation=0, delta_g_rbs_linker=0, weights=ScoringWeights(defect=6))
        self.assertEqual(compute_selection_score(**args), -120)
        for invalid in (None, True, float('nan'), 1.1, -.1):
            with self.subTest(value=invalid), self.assertRaises(ValueError):
                compute_selection_score(**{**args, 'functional_defect': invalid})

    def test_functional_tie_breaking_and_zero_penalty(self):
        candidates = [SimpleNamespace(trial=1, excluded=False, score=5, defect=.1, regional_defects={'balanced': .7}),
                      SimpleNamespace(trial=2, excluded=False, score=5, defect=.9, regional_defects={'balanced': .2})]
        self.assertEqual(sort_candidates(candidates), candidates[::-1])
        self.assertEqual(sort_candidates(candidates, defect_weight=0), candidates)
