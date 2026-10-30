"""RBS scoring steers refinement and survives persistence without relabelling."""
from copy import deepcopy
from dataclasses import replace
import json
import math
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from aptaswitch_core.engines import run_design
from aptaswitch_core.models import InverseExtensionSettings
from aptaswitch_core.scoring import compute_selection_score, selection_score_100
from aptaswitch_core.thermodynamics.design_engine import design_candidates
from aptaswitch_core.thermodynamics.optimization import optimize, TrialOptions
from aptaswitch_studio.run_library import _load_design
from tests.test_thermodynamics_design_engine import request, AnalyticalBackend
from tests.test_thermodynamics_optimization import binary_problem, objective, options


class RBSDesignIntegrationTests(unittest.TestCase):
    def test_merit_and_linearity_steer_search_without_buying_local_violations(self):
        for field, best, worse in (("rbs_calculator_merit", .9, .1),
                                    ("trigger_linearity_defect", .1, .9)):
            weight = "rbs_calculator_weight" if field.startswith('rbs') else 'trigger_linearity_weight'
            def loss(assignment):
                return objective(**{field: best if assignment.sequence('S') == 'U' else worse, weight: 3.})
            result = optimize(binary_problem(), loss, options=options())
            self.assertEqual(result.best_assignment.sequence('S'), 'U')
            def guarded(assignment):
                bad_local = assignment.sequence('S') == 'U'
                return objective(local=.4 if bad_local else .05,
                    **{field: best if bad_local else worse, weight: 3.})
            result = optimize(binary_problem(), guarded, options=options())
            self.assertEqual(result.best_assignment.sequence('S'), 'A')

    def test_disabled_rbs_never_loads_ostir_and_keeps_original_physical_metrics(self):
        req = request()
        req = replace(req, scoring=replace(req.scoring, rbs_calculator=0.))
        with patch('aptaswitch_core.rbs_calculator.predict_rbs_switch', side_effect=AssertionError('disabled')):
            batch = design_candidates(req, AnalyticalBackend(req), options=TrialOptions(max_evaluations=1, max_proposals=1))
        self.assertEqual(batch.candidates[0].rbs_calculator, {})
        self.assertEqual(batch.optimization_trials[0].best_objective.rbs_calculator_weight, 0.)

    def test_weighted_candidate_round_trip_and_tampering(self):
        with tempfile.TemporaryDirectory() as directory:
            req = replace(request(), output_dir=Path(directory), optimization_options={'max_evaluations':1, 'max_proposals':1})
            backend = AnalyticalBackend(req)
            backend.provenance = replace(backend.provenance, model=req.thermo.material)
            with patch('aptaswitch_core.thermodynamics.setup.backend_for_request', return_value=backend):
                result = run_design(req)
            payload = json.loads(json.dumps(result.to_dict(), allow_nan=False))
            loaded = _load_design(payload, result.output_dir/'result.json')
            self.assertEqual(loaded.candidates[0].rbs_calculator, result.candidates[0].rbs_calculator)
            self.assertEqual(loaded.candidates[0].score_raw, result.candidates[0].score_raw)
            from aptaswitch_studio.saved_switch_library import candidate_snapshot
            snapshot = candidate_snapshot(result, result.candidates[0], 'RBS candidate')
            self.assertEqual(snapshot['candidate']['rbs_calculator'], result.candidates[0].rbs_calculator)
            for mutate in (lambda p:p['results'][0]['rbs_calculator'].__setitem__('merit', .123),
                           lambda p:p['results'][0].__setitem__('rbs_calculator', {}),
                           lambda p:p['optimization']['trials'][0]['best_objective'].__setitem__('rbs_calculator_weight', 0.)):
                bad = deepcopy(payload); mutate(bad)
                with self.assertRaises(ValueError):
                    _load_design(bad, result.output_dir/'result.json')

    def test_bonus_is_exact_and_score_stays_bounded(self):
        weights = request().scoring
        kwargs=dict(on_yield_pct=100.,leak_pct=50.,ddg_activation=-10.,delta_g_rbs_linker=-1.,
                    weights=weights,functional_defect=.1)
        base=compute_selection_score(**kwargs,rbs_calculator_merit=0.)
        for merit in (0., .1, .5, 1.):
            score=compute_selection_score(**kwargs,rbs_calculator_merit=merit)
            self.assertAlmostEqual(score-base,100*weights.rbs_calculator*merit)
            self.assertTrue(0 <= selection_score_100(score) <= 100)
        for bad in (True, -.1, 1.1, math.nan):
            with self.assertRaises(ValueError):compute_selection_score(**kwargs,rbs_calculator_merit=bad)
