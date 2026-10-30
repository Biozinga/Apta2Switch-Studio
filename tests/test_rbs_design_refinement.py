"""RBS stability guides search without replacing functional-region priorities."""
from dataclasses import replace
from copy import deepcopy
import json
import math
import random
from types import SimpleNamespace
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from aptaswitch_core.thermodynamics.design_constraints import build_design_problem
from aptaswitch_core.thermodynamics.design_engine import design_candidates
from aptaswitch_core.thermodynamics.design_objective import objective_from_evaluation
from aptaswitch_core.thermodynamics.evaluation import evaluate_design
from aptaswitch_core.thermodynamics.optimization import optimize
from tests.test_thermodynamics_optimization import binary_problem, objective, options
from tests.test_thermodynamics_design_engine import AnalyticalBackend, request


class RBSRefinementTests(unittest.TestCase):
    def test_invalid_request_weight_fails_before_any_backend_work(self):
        for weight in (-1., math.nan, math.inf, True):
            req = request()
            req = replace(req, scoring=replace(req.scoring, rbs=weight))
            backend = AnalyticalBackend(req)
            with self.subTest(weight=weight), self.assertRaises(ValueError):
                design_candidates(req, backend)
            self.assertEqual(backend.calls, [])

    def test_continuous_penalty_and_disabled_weight(self):
        value = objective(global_loss=.2, rbs_linker_energy=-8., rbs_weight=6.)
        self.assertAlmostEqual(value.refinement_loss, .68)
        self.assertAlmostEqual(value.rbs_stability_penalty, .48)
        self.assertEqual(replace(value, rbs_weight=0).refinement_loss, .2)
        self.assertEqual(replace(value, rbs_linker_energy=0).refinement_loss, .2)
        for changes in ({'rbs_linker_energy': math.nan}, {'rbs_weight': -1},
                        {'rbs_linker_energy': True}, {'rbs_weight': math.inf}):
            with self.subTest(changes=changes), self.assertRaises(ValueError):
                replace(value, **changes)

    def test_less_stable_fragment_is_selected_among_equivalent_local_candidates(self):
        def loss(assignment):
            return objective(rbs_linker_energy=-12. if assignment.sequence('S') == 'A' else -1., rbs_weight=6.)
        result = optimize(binary_problem(), loss, options=options())
        self.assertEqual(result.best_assignment.sequence('S'), 'U')
        self.assertEqual(result.best_objective.rbs_linker_energy, -1.)

    def test_stability_cannot_buy_a_violation_of_local_ceilings(self):
        def loss(assignment):
            good = assignment.sequence('S') == 'A'
            return objective(local=.1 if good else .3, global_loss=.3 if good else 0.,
                             rbs_linker_energy=-20. if good else 0., rbs_weight=6.)
        result = optimize(binary_problem(), loss, options=options())
        self.assertEqual(result.best_assignment.sequence('S'), 'A')
        self.assertTrue(all(elite.objective.local_defects[0] <= result.local_limits[0]
                            for elite in result.elites))

    def test_bad_stability_prevents_false_target_reached(self):
        result = optimize(binary_problem(), lambda _: objective(local=0., global_loss=0.,
            rbs_linker_energy=-10., rbs_weight=6.), options=options(target_local_loss=.005, target_global_loss=.005))
        self.assertNotEqual(result.stop_reason, 'target_reached')

    def test_fragment_minimum_is_computed_before_tube_analyses_and_recorded(self):
        req = request()
        backend = AnalyticalBackend(req)
        from aptaswitch_core.thermodynamics.optimization import TrialOptions
        batch = design_candidates(req, backend, options=TrialOptions(max_evaluations=1, max_proposals=1))
        self.assertEqual(backend.calls[0], ('minimum', ('RBS_linker',)))
        self.assertEqual(sum(call == ('minimum', ('RBS_linker',)) for call in backend.calls), 1)
        self.assertEqual(batch.optimization_trials[0].best_objective.rbs_weight, req.scoring.rbs)
        self.assertEqual(batch.optimization_trials[0].best_objective.rbs_linker_energy, -1.)
        self.assertEqual(batch.candidates[0].rbs_accessibility['probability'], 1.)

    def test_mfe_paired_sites_guide_refinement_but_local_defects_stay_identical(self):
        req = request()
        problem = build_design_problem(req)
        backend = AnalyticalBackend(req)
        evaluation = evaluate_design(problem, problem.generate(random.Random(8)), backend.analyze,
                                     standard_state=backend.standard_state)
        plain = objective_from_evaluation(req, evaluation)
        length = problem.layout['rbs_linker_end'] - problem.layout['rbs_linker_start']
        guided = objective_from_evaluation(req, evaluation, rbs_minimum=SimpleNamespace(
            structure='(' + '.' * (length-2) + ')', energy_kcal_mol=-5.))
        self.assertEqual(plain.local_defects, guided.local_defects)
        self.assertEqual(plain.site_weights, guided.site_weights)
        self.assertEqual(plain.guard_defects, guided.guard_defects)
        differences = {site.index: weight-dict(plain.global_site_weights).get(site, 0.)
                       for site, weight in guided.global_site_weights if weight != dict(plain.global_site_weights).get(site, 0.)}
        self.assertEqual(set(differences), {problem.layout['rbs_linker_start'], problem.layout['rbs_linker_end']-1})
        self.assertAlmostEqual(sum(differences.values()), .3)

    def test_archived_accessibility_is_preserved_and_wrong_context_rejected(self):
        from aptaswitch_core.engines import run_design
        from aptaswitch_studio.run_library import _load_design
        with tempfile.TemporaryDirectory() as directory:
            req = replace(request(), output_dir=Path(directory),
                          optimization_options={'max_evaluations': 1, 'max_proposals': 1})
            backend = AnalyticalBackend(req)
            backend.provenance = replace(backend.provenance, model=req.thermo.material)
            with patch('aptaswitch_core.thermodynamics.setup.backend_for_request', return_value=backend):
                result = run_design(req)
            payload = json.loads(json.dumps(result.to_dict(), allow_nan=False))
            source = result.output_dir / 'result.json'
            loaded = _load_design(payload, source)
            self.assertEqual(result.candidates[0].rbs_accessibility, loaded.candidates[0].rbs_accessibility)
            mutations = (
                lambda p: p['results'][0]['rbs_accessibility'].update(start=0),
                lambda p: p['results'][0]['rbs_accessibility']['conditions'].update(celsius=25),
                lambda p: p['results'][0]['rbs_accessibility']['complex_spec']['strands'][0].update(material='rna', sequence='A'*30),
                lambda p: p['results'][0]['rbs_accessibility']['provenance'].update(version='other'),
                lambda p: p['optimization']['trials'][0]['best_objective'].update(rbs_weight=0),
                lambda p: p['optimization']['rbs_stability'].update(weight=0),
            )
            for mutate in mutations:
                changed = deepcopy(payload)
                mutate(changed)
                with self.assertRaises(ValueError):
                    _load_design(changed, source)

    def test_checkpoint_cannot_silently_resume_a_different_search_method(self):
        from aptaswitch_core.engines import run_design
        from aptaswitch_core.design_checkpoint import save_design_checkpoint
        with tempfile.TemporaryDirectory() as directory:
            req = replace(request(), output_dir=Path(directory),
                          optimization_options={'max_evaluations': 1, 'max_proposals': 1})
            backend = AnalyticalBackend(req)
            backend.provenance = replace(backend.provenance, model=req.thermo.material)
            with patch('aptaswitch_core.thermodynamics.setup.backend_for_request', return_value=backend):
                paused = run_design(req, control_check=lambda: 'pause')
                paused.optimization.pop('rbs_stability')
                save_design_checkpoint(paused)
                backend.calls.clear()
                with self.assertRaises(ValueError):
                    run_design(req, resume_result=paused)
                self.assertEqual(backend.calls, [])
