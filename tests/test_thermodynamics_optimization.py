"""Population/lineage optimizer invariants with deterministic toy objectives."""
from dataclasses import FrozenInstanceError, replace
import math
import unittest
from unittest.mock import patch

from aptaswitch_core.thermodynamics.design_constraints import (
    BaseSite, DesignAssignment, DesignProblem, LibraryConstraint, ScreeningPolicy,
    StrandPattern, UnsatisfiableDesignError,
)
from aptaswitch_core.thermodynamics.optimization import ObjectiveValue, TrialOptions, assignment_key, optimize
from aptaswitch_core.thermodynamics.types import Strand


def binary_problem(length=1):
    return DesignProblem((StrandPattern("S", "rna", "W" * length),))


def assignment(sequence):
    return DesignAssignment((Strand(sequence, "rna", "S"),))


def objective(local=.2, global_loss=.3, guards=(.1, .1, .1), weights=(1, 1, 1, .75), **kw):
    values = (local,) * 4 if isinstance(local, (int, float)) else local
    return ObjectiveValue(values, weights, global_loss, guards, **kw)


def options(**changes):
    return replace(TrialOptions(seed=7, max_evaluations=100, max_proposals=100,
                   initial_temperature=0, final_temperature=0, target_local_loss=0,
                   target_global_loss=0), **changes)


class OptimizationTests(unittest.TestCase):
    def test_defaults_and_normalized_weighted_local_objective(self):
        config = TrialOptions()
        self.assertEqual((config.population_size, config.lineages), (64, 8))
        self.assertEqual((config.local_fraction, config.local_tolerance, config.guard_tolerance), (.7, .005, .01))
        self.assertEqual((config.target_local_loss, config.target_global_loss), (.005, .005))
        value = objective((.1, .2, .3, .4))
        self.assertAlmostEqual(value.local_loss, (.1 + .2 + .3 + .4 * .75) / 3.75)
        self.assertEqual(value.local_weights, (1, 1, 1, .75))
        self.assertAlmostEqual(objective((.1, .2, .3, .4), weights=(1e308,) * 4).local_loss, .25)
        with self.assertRaises(FrozenInstanceError):
            value.global_defect = .1

    def test_coupled_population_and_mutations_preserve_libraries(self):
        patterns = (StrandPattern("D", "dna", "XXX"), StrandPattern("R", "rna", "XXX"))
        libraries = tuple(LibraryConstraint(f"pair{i}", (BaseSite("D", i), BaseSite("R", 2-i)),
                                           ("AU", "CG", "GC", "TA")) for i in range(3))
        problem = DesignProblem(patterns, libraries)
        evaluated, snapshots = [], []
        def loss(value):
            problem.validate(value)
            evaluated.append(value)
            mismatch = sum(a != b for a, b in zip(value.sequence("D"), "ACG")) / 3
            return objective(mismatch, mismatch)
        result = optimize(problem, loss, options=options(max_proposals=400), progress_callback=snapshots.append)
        self.assertEqual(result.stop_reason, "target_reached")
        self.assertEqual(result.best_assignment.sequence("D"), "ACG")
        self.assertEqual(result.best_assignment.sequence("R"), "CGU")
        self.assertEqual(result.evaluations, len(evaluated))
        self.assertEqual(result.best_global_loss, result.best_objective.global_defect)
        for row in snapshots:
            problem.validate(row.current_assignment)
            problem.validate(row.best_assignment)

    def test_seed_replays_every_phase_and_progress(self):
        problem = binary_problem(7)
        loss = lambda value: objective(.1 + value.sequence("S").count("A") / 10)
        config = options(max_evaluations=35, max_proposals=60, initial_temperature=.2, final_temperature=.01)
        first, second = [], []
        a = optimize(problem, loss, options=config, progress_callback=first.append)
        b = optimize(problem, loss, options=config, progress_callback=second.append)
        self.assertEqual(a, b)
        self.assertEqual(first, second)
        self.assertEqual({row.phase for row in first}, {"initialization", "local", "global"})
        self.assertLessEqual(a.population_evaluated, config.max_evaluations // 5)
        self.assertGreater(a.lineages_started, 1)
        alternate = []
        optimize(problem, loss, options=replace(config, seed=24), progress_callback=alternate.append)
        self.assertNotEqual(first, alternate)

    def test_unique_small_pool_and_diverse_lineages_get_round_robin_turns(self):
        origins = []
        def mutate(problem, value, rng, weights, **kwargs):
            origins.append(value.sequence("S"))
            return value
        pool = [assignment(seq) for seq in ("AAAA", "UUUU", "AUUA", "UAAU", "AAUU", "UUAA", "AUAU", "UAUA")]
        with patch.object(DesignProblem, "mutate_guided", mutate), patch.object(DesignProblem, "generate", side_effect=pool):
            result = optimize(binary_problem(4), lambda _: objective(),
                              options=options(max_evaluations=40, max_proposals=40, population_size=8, lineages=4))
        self.assertEqual(result.population_evaluated, 8)
        self.assertEqual(result.evaluations, 8)
        self.assertEqual(result.lineages_started, 4)
        self.assertEqual(len(set(origins[:4])), 4)
        self.assertEqual(origins[0:4], origins[4:8])
        # Farthest-first picks the opposite corner as its second seed.
        self.assertTrue(all(a != b for a, b in zip(origins[0], origins[1])))
        self.assertEqual(result.proposals, 40)
        self.assertLessEqual(len(result.elites), 4)

    def test_evaluation_and_generation_budgets_are_shared_and_never_exceeded(self):
        for evaluation_budget, proposal_budget in ((1, 1), (2, 3), (3, 50), (10, 2), (17, 50)):
            calls = []
            result = optimize(binary_problem(8), lambda value: calls.append(value) or objective(),
                options=options(max_evaluations=evaluation_budget, max_proposals=proposal_budget, cache_size=0))
            with self.subTest(evaluations=evaluation_budget, proposals=proposal_budget):
                self.assertLessEqual(result.evaluations, evaluation_budget)
                self.assertLessEqual(result.proposals, proposal_budget)
                self.assertEqual(result.evaluations, len(calls))
                self.assertEqual(result.stop_reason, "budget_exhausted")
                self.assertEqual(result.evaluations if result.budget_kind == "evaluations" else result.proposals,
                                 evaluation_budget if result.budget_kind == "evaluations" else proposal_budget)

    def test_cache_cycles_spend_proposals_and_disabled_cache_recomputes(self):
        for capacity in (0, 1, 2):
            calls = []
            result = optimize(binary_problem(), lambda value: calls.append(value) or objective(),
                initial_assignment=assignment("A"),
                options=options(population_size=1, max_proposals=20, cache_size=capacity))
            self.assertEqual(result.stop_reason, "budget_exhausted")
            self.assertEqual(result.budget_kind, "proposals")
            self.assertEqual(result.proposals, 20)
            self.assertEqual(result.evaluations, len(calls))
            if capacity == 2:
                self.assertEqual((result.evaluations, result.cache_hits), (2, 18))
            elif capacity == 0:
                self.assertEqual((result.evaluations, result.cache_hits), (20, 0))
            else:
                self.assertGreater(result.evaluations, 2)
                self.assertEqual(result.evaluations + result.cache_hits, 20)

    def test_global_refinement_cannot_drift_local_or_guard_ceilings(self):
        problem = binary_problem(4)
        start = assignment("AAAA")
        proposals = iter(("UAAA", "UUAA", "UUUA", "UUUU"))
        losses = {"AAAA": objective(.1, .8), "UAAA": objective(.105, .6),
                  "UUAA": objective(.110, .4), "UUUA": objective(.115, .2),
                  "UUUU": objective(.1, .1, guards=(.111, .1, .1))}
        phase = "initialization"
        snapshots = []
        def progress(row):
            nonlocal phase
            snapshots.append(row)
            phase = row.phase
        def mutate(problem, value, rng, weights, **kwargs):
            return value if phase == "local" else assignment(next(proposals, "UUUU"))
        with patch.object(DesignProblem, "mutate_guided", mutate):
            result = optimize(problem, lambda value: losses[value.sequence("S")], initial_assignment=start,
                options=options(population_size=1, lineages=1, local_fraction=.5, max_proposals=11),
                progress_callback=progress)
        self.assertEqual(result.best_assignment.sequence("S"), "UAAA")
        self.assertAlmostEqual(result.best_loss, .105)
        self.assertEqual(result.best_global_loss, .6)
        global_rows = [row for row in snapshots if row.phase == "global"]
        self.assertTrue(global_rows)
        self.assertEqual(len({row.local_limits for row in global_rows}), 1)
        self.assertEqual(len({row.guard_limits for row in global_rows}), 1)
        self.assertEqual(result.local_limits, (.10500000000000001,) * 4)
        self.assertEqual(result.guard_limits, (.11,) * 3)
        for elite in result.elites:
            self.assertTrue(all(a <= b + 1e-12 for a, b in zip(elite.objective.local_defects, result.local_limits)))
            self.assertTrue(all(a <= b for a, b in zip(elite.objective.guard_defects, result.guard_limits)))

    def test_global_lineages_start_only_inside_the_final_local_neighborhood(self):
        phase, origins, snapshots = "initialization", [], []
        def progress(row):
            nonlocal phase
            phase = row.phase
            snapshots.append(row)
        def mutate(problem, value, rng, weights, **kwargs):
            if phase == "global":
                origins.append(value.sequence("S"))
            return value
        with patch.object(DesignProblem, "generate", side_effect=[assignment("AA"), assignment("UU")]), \
                patch.object(DesignProblem, "mutate_guided", mutate):
            result = optimize(binary_problem(2), lambda value: objective(.1 if value.sequence("S") == "AA" else .3),
                options=options(max_evaluations=10, max_proposals=20, population_size=2, lineages=2),
                progress_callback=progress)
        self.assertEqual(result.lineages_started, 2)
        self.assertTrue(origins)
        self.assertEqual(set(origins), {"AA"})
        for row in snapshots:
            if row.phase == "global":
                self.assertTrue(all(value <= ceiling + 1e-12
                    for value, ceiling in zip(row.current_objective.local_defects, row.local_limits)))

    def test_disabled_local_zone_does_not_constrain_global_refinement(self):
        rows = {"A": objective((.1, .1, .1, .1), .8, weights=(1, 1, 1, 0)),
                "U": objective((.1, .1, .1, .9), .1, weights=(1, 1, 1, 0))}
        result = optimize(binary_problem(), lambda value: rows[value.sequence("S")],
            initial_assignment=assignment("A"), options=options(population_size=1))
        self.assertEqual(result.local_limits[-1], 1)
        self.assertEqual(result.best_assignment.sequence("S"), "U")

    def test_guard_anchor_uses_best_maximum_shortage_not_lowest_local_loss(self):
        rows = {"AA": objective(.8, .8, guards=(.1, .1, .1)),
                "AU": objective(.1, .1, guards=(.4, .1, .1)),
                "UA": objective(.2, .2, guards=(.1, .4, .1)),
                "UU": objective(.3, .3, guards=(.1, .1, .4))}
        result = optimize(binary_problem(2), lambda value: rows[value.sequence("S")],
                          options=options(max_evaluations=20, max_proposals=20, population_size=4))
        self.assertEqual(result.guard_limits, (.11,) * 3)
        self.assertEqual(result.best_assignment.sequence("S"), "AA")
        self.assertEqual(len(result.elites), 1)

    def test_target_needs_global_local_and_guard_thresholds(self):
        for local, global_loss, reason in ((.004, .004, "target_reached"),
                                         (.004, .2, "budget_exhausted"),
                                         (.2, .004, "budget_exhausted")):
            result = optimize(binary_problem(), lambda _: objective(local, global_loss),
                              options=replace(TrialOptions(), max_evaluations=1))
            self.assertEqual(result.stop_reason, reason)
        # A local target cannot cause target_reached if the returned global
        # compromise leaves that target while staying inside local tolerance.
        rows = {"A": objective(.001, .005), "U": objective(.006, .001)}
        result = optimize(binary_problem(), lambda value: rows[value.sequence("S")],
                          options=replace(TrialOptions(), max_evaluations=10, max_proposals=10))
        self.assertEqual(result.stop_reason, "space_exhausted")
        self.assertEqual(result.best_assignment.sequence("S"), "U")
        self.assertGreater(result.best_loss, result.options.target_local_loss)

    def test_sequence_prefilter_avoids_all_folding_calls_for_impossible_stop_policy(self):
        screening = ScreeningPolicy(0, 6, True, False, 0, 9)
        problem = DesignProblem((StrandPattern("Switch", "rna", "AUGUAAUGG"),), screening=screening)
        result = optimize(problem, lambda _: self.fail("STOP sequence must not be folded"), options=options())
        self.assertEqual(result.stop_reason, "no_admissible_sequence")
        self.assertEqual(result.evaluations, 0)
        self.assertIsNone(result.best_assignment)
        variable = replace(problem, strands=(StrandPattern("Switch", "rna", "AUGUAAXXX"),))
        result = optimize(variable, lambda _: self.fail("STOP sequence must not be folded"), options=options())
        self.assertEqual(result.stop_reason, "no_admissible_sequence")
        self.assertEqual(result.evaluations, 0)
        self.assertEqual(result.proposals, result.options.max_proposals // 5)
        with patch.object(DesignProblem, "passes_sequence_screen", side_effect=AssertionError("No policy")):
            optimize(binary_problem(), lambda _: objective(), options=options(max_evaluations=1))

    def test_guidance_uses_local_then_global_errors_and_exploration(self):
        local = ((BaseSite("S", 0), 1.0),)
        global_weights = ((BaseSite("S", 1), 1.0),)
        calls = []
        original = DesignProblem.mutate_guided
        def record(problem, value, rng, weights, **kwargs):
            calls.append((weights, kwargs))
            return original(problem, value, rng, weights, **kwargs)
        with patch.object(DesignProblem, "mutate_guided", record):
            optimize(binary_problem(4), lambda _: objective(site_weights=local, global_site_weights=global_weights),
                     options=options(population_size=1, max_proposals=20, exploration_fraction=.35))
        self.assertIn(local, [row[0] for row in calls])
        self.assertIn(global_weights, [row[0] for row in calls])
        self.assertTrue(all(row[1]["exploration_fraction"] == .35 for row in calls))

    def test_diverse_final_elites_are_presented_in_quality_order_for_screening(self):
        globals_by_sequence = {"AA": .1, "AU": .2, "UA": .3, "UU": .4}
        result = optimize(binary_problem(2), lambda value: objective(.2, globals_by_sequence[value.sequence("S")]),
                          options=options())
        self.assertEqual(tuple(item.objective.global_defect for item in result.elites), (.1, .2, .3, .4))
        self.assertEqual(result.best_assignment, result.elites[0].assignment)
        self.assertEqual(result.best_objective, result.elites[0].objective)

    def test_fully_enumerated_space_stops_without_redundant_cache_cycles(self):
        calls = []
        result = optimize(binary_problem(2), lambda value: calls.append(value) or objective(), options=options())
        self.assertEqual(result.stop_reason, "space_exhausted")
        self.assertEqual((result.evaluations, result.proposals), (4, 4))
        self.assertEqual(result.lineages_started, 0)
        self.assertEqual(len({assignment_key(item) for item in calls}), 4)

    def test_no_mutable_sites_is_distinct(self):
        fixed = DesignProblem((StrandPattern("S", "rna", "A"),))
        result = optimize(fixed, lambda _: objective(), options=options())
        self.assertEqual(result.stop_reason, "no_mutable_sites")
        self.assertEqual((result.evaluations, result.proposals), (1, 1))
        self.assertIsNone(result.budget_kind)

    def test_cancel_before_work_during_evaluation_and_after_progress_keeps_completed_work(self):
        with patch.object(DesignProblem, "generate", side_effect=AssertionError("Must not generate")):
            result = optimize(binary_problem(), lambda _: self.fail("Must not fold"), cancel_check=lambda: True)
        self.assertEqual(result.stop_reason, "cancelled")
        self.assertEqual(result.evaluations, 0)
        self.assertIsNone(result.best_assignment)
        cancelled = False
        def loss(value):
            nonlocal cancelled
            cancelled = True
            return objective()
        result = optimize(binary_problem(), loss, cancel_check=lambda: cancelled)
        self.assertEqual(result.stop_reason, "cancelled")
        self.assertEqual(result.evaluations, 1)
        self.assertAlmostEqual(result.best_loss, .2)
        cancelled = False
        snapshots = []
        def progress(row):
            nonlocal cancelled
            snapshots.append(row)
            if row.phase == "local":
                cancelled = True
        result = optimize(binary_problem(8), lambda _: objective(), progress_callback=progress,
                          cancel_check=lambda: cancelled)
        self.assertEqual(result.stop_reason, "cancelled")
        self.assertIsNotNone(result.best_assignment)
        self.assertTrue(any(result.best_assignment == item.assignment for item in result.elites))

    def test_loss_cache_identity_and_invocation_isolation(self):
        a, b = Strand("A", "rna", "A"), Strand("A", "rna", "B")
        keys = {assignment_key(DesignAssignment(strands)) for strands in (
            (a, b), (b, a), (replace(a, material="dna"), b), (replace(a, name="C"), b),
            (replace(a, sequence="G"), b))}
        self.assertEqual(len(keys), 5)
        a = optimize(binary_problem(), lambda _: objective(.1), options=options(max_evaluations=1))
        b = optimize(binary_problem(), lambda _: objective(.8), options=options(max_evaluations=1))
        self.assertAlmostEqual(a.best_loss, .1)
        self.assertAlmostEqual(b.best_loss, .8)

    def test_objective_validation_and_errors_propagate(self):
        for value in (-1, math.nan, math.inf, -math.inf, True, "0", 1j, None, 1.01):
            with self.subTest(value=value), self.assertRaises((ValueError, TypeError)):
                objective(global_loss=value)
        for weights in ((0, 0, 0, 0), (1, 1, -1, 1), (1, 1, math.nan, 1)):
            with self.assertRaises(ValueError):
                objective(weights=weights)
        for value in (.1, None, "0"):
            with self.assertRaises(ValueError):
                optimize(binary_problem(), lambda _: value)
        for callback in ("loss", "progress_callback", "cancel_check"):
            def raises(*args):
                raise RuntimeError("intended failure")
            arguments = {"loss": lambda _: objective(), callback: raises}
            with self.subTest(callback=callback), self.assertRaisesRegex(RuntimeError, "intended failure"):
                optimize(binary_problem(), **arguments)
        impossible = DesignProblem((StrandPattern("S", "rna", "AA"),),
            (LibraryConstraint("pair", (BaseSite("S", 0), BaseSite("S", 1)), ("AU", "UA")),))
        with self.assertRaises(UnsatisfiableDesignError):
            optimize(impossible, lambda _: objective())
        with self.assertRaisesRegex(ValueError, "length"):
            optimize(binary_problem(), lambda _: objective(), initial_assignment=assignment("AA"))

    def test_options_reject_invalid_numbers_but_adapt_tiny_budgets(self):
        invalid = [dict(seed=True), dict(seed=1.2), dict(max_evaluations=0), dict(max_evaluations=True),
                   dict(max_proposals=0), dict(max_proposals=1.2), dict(population_size=0), dict(lineages=True),
                   dict(cache_size=-1), dict(cache_size=1.2), dict(local_fraction=0), dict(local_fraction=1),
                   dict(initial_temperature=-1), dict(final_temperature=-1), dict(initial_temperature=math.inf),
                   dict(final_temperature=math.nan), dict(initial_temperature=True), dict(target_local_loss=-1),
                   dict(target_global_loss=math.nan), dict(target_local_loss=True), dict(target_global_loss="0.1"),
                   dict(initial_temperature=0, final_temperature=.1), dict(local_tolerance=1.1),
                   dict(exploration_fraction=-.1), dict(guard_tolerance=-.1)]
        for changes in invalid:
            with self.subTest(changes=changes), self.assertRaises(ValueError):
                replace(TrialOptions(), **changes)
        TrialOptions(max_evaluations=1, max_proposals=1, population_size=1, lineages=8)


if __name__ == "__main__":
    unittest.main()
