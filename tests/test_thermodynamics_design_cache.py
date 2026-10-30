"""Exact run-local reuse, lazy finalist exclusion and retained matrix budgets."""
from collections import Counter
from dataclasses import replace
import random
from types import SimpleNamespace
import unittest
from unittest.mock import patch

from aptaswitch_core.thermodynamics import design_engine
from aptaswitch_core.thermodynamics.design_constraints import build_design_problem
from aptaswitch_core.thermodynamics.evaluation import ComplexAnalysisCache, evaluate_design
from aptaswitch_core.thermodynamics.optimization import TrialOptions
from aptaswitch_core.thermodynamics.structures import PairingConstraints
from aptaswitch_core.thermodynamics.types import ComplexSpec, Conditions, PartitionResult, Strand
from tests.test_thermodynamics_design_engine import AnalyticalBackend, request


class ObservableBackend:
    provenance = AnalyticalBackend.provenance
    standard_state = AnalyticalBackend.standard_state

    def __init__(self):
        self.calls = []

    def minimum(self, spec, conditions):
        self.calls.append(("minimum", spec, conditions))
        context = SimpleNamespace(complex_spec=spec, conditions=conditions,
                                  provenance=self.provenance, constraints=None)
        return SimpleNamespace(context=context, energy_kcal_mol=-1., structure="." * spec.n_nucleotides)

    def screen_unpaired_region(self, spec, conditions, *, start, end):
        self.calls.append(("screen", spec, conditions, start, end))
        context = SimpleNamespace(complex_spec=spec, conditions=conditions,
                                  provenance=self.provenance, constraints=None)
        return SimpleNamespace(context=context, all_minima_unpaired=True, witness=None)

    def partition(self, spec, conditions, constraints=None):
        self.calls.append(("partition", spec, conditions, constraints))
        return PartitionResult(spec, conditions, self.provenance, 1., -1., 1., constraints)


class RunCalculationCacheTests(unittest.TestCase):
    def test_invariant_final_validation_and_trigger_minimum_are_reused_between_trials(self):
        req = replace(request(), trials=2)
        backend = AnalyticalBackend(req)
        batch = design_engine.design_candidates(req, backend, options=TrialOptions(
            max_evaluations=1, max_proposals=1))
        self.assertEqual(len(batch.candidates), 2)
        analyze = Counter(call[1] for call in backend.calls if call[0] == "analyze")
        for names in (("Aptamer",), ("Aptamer", "Aptamer"), ("Aptamer", "Trigger"),
                      ("Trigger",), ("Trigger", "Trigger")):
            self.assertEqual(analyze[names], 1, names)
        minima = Counter(call[1] for call in backend.calls if call[0] == "minimum")
        self.assertEqual(minima[("Trigger",)], 1)
        self.assertEqual(batch.analysis_cache["backend_calls"], len(backend.calls))
        # The analytical fixture implements a constrained PF by filtering one
        # prescribed analysis; production records that top-level PF separately.
        self.assertEqual(batch.analysis_cache["analysis_calls"] + batch.analysis_cache["partition_calls"],
                         sum(analyze.values()))
        self.assertEqual(batch.analysis_cache["minimum_calls"], sum(minima.values()))
        self.assertEqual(batch.analysis_cache["screen_calls"], 2)

    def test_observable_keys_preserve_labels_chemistry_order_conditions_and_region(self):
        backend = ObservableBackend()
        cache = ComplexAnalysisCache(backend)
        dna = Strand("ACG", "dna", "A")
        other = Strand("GCA", "rna", "B")
        specs = (ComplexSpec((dna,)), ComplexSpec((replace(dna, material="rna"),)),
                 ComplexSpec((replace(dna, name="renamed"),)),
                 ComplexSpec((dna, other)), ComplexSpec((other, dna)))
        conditions = Conditions()
        for spec in specs:
            first = cache.minimum(spec, conditions)
            self.assertIs(cache.minimum(spec, conditions), first)
        cache.minimum(specs[0], replace(conditions, celsius=40))
        self.assertEqual(len(backend.calls), len(specs) + 1)
        first = cache.screen_unpaired_region(specs[0], conditions, start=0, end=1)
        self.assertIs(cache.screen_unpaired_region(specs[0], conditions, start=0, end=1), first)
        cache.screen_unpaired_region(specs[0], conditions, start=0, end=2)
        self.assertEqual(len(backend.calls), len(specs) + 3)
        with self.assertRaises(ValueError):
            cache.screen_unpaired_region(specs[0], conditions, start=False, end=1)
        backend.provenance = replace(backend.provenance, version="changed")
        with self.assertRaisesRegex(ValueError, "backend changed"):
            cache.minimum(specs[0], conditions)
        with self.assertRaisesRegex(ValueError, "backend changed"):
            cache.screen_unpaired_region(specs[0], conditions, start=0, end=1)

    def test_completed_observables_share_the_bounded_entry_budget(self):
        backend = ObservableBackend()
        cache = ComplexAnalysisCache(backend, max_entries=1)
        spec, conditions = ComplexSpec((Strand("ACG", "rna", "S"),)), Conditions()
        cache.minimum(spec, conditions)
        cache.screen_unpaired_region(spec, conditions, start=0, end=1)
        cache.minimum(spec, conditions)
        self.assertEqual(len(backend.calls), 3)
        self.assertEqual(cache.statistics()["cached_observables"], 1)
        self.assertEqual(cache.statistics()["cached_pair_cells"], 0)

    def test_observable_with_incompatible_context_is_not_cached(self):
        backend = ObservableBackend()
        cache = ComplexAnalysisCache(backend)
        spec, conditions = ComplexSpec((Strand("ACG", "rna", "S"),)), Conditions()
        result = backend.minimum(spec, conditions)
        result.context.conditions = replace(conditions, celsius=40)
        with patch.object(backend, "minimum", return_value=result):
            with self.assertRaisesRegex(ValueError, "exact unrestricted"):
                cache.minimum(spec, conditions)
        self.assertEqual(cache.statistics()["cached_observables"], 0)

    def test_partial_partitions_preserve_constraints_and_share_entry_budget(self):
        backend = ObservableBackend()
        cache = ComplexAnalysisCache(backend, max_entries=2)
        spec, conditions = ComplexSpec((Strand("ACG", "rna", "S"),)), Conditions()
        first = PairingConstraints((0, -1, -1))
        second = PairingConstraints((-1, 1, -1))
        result = cache.partition(spec, conditions, first)
        self.assertIs(cache.partition(spec, conditions, first), result)
        self.assertEqual(cache.partition(spec, conditions, second).constraints, second)
        self.assertEqual(cache.statistics()["partition_calls"], 2)
        cache.minimum(spec, conditions)
        cache.partition(spec, conditions, first)
        self.assertEqual(cache.statistics()["partition_calls"], 3)
        self.assertEqual(cache.statistics()["cached_observables"], 2)
        with patch.object(backend, "partition", return_value=replace(result, constraints=second)):
            with self.assertRaisesRegex(ValueError, "exact requested"):
                cache.partition(spec, replace(conditions, celsius=40), first)

    def test_unrestricted_partition_reuses_analysis_including_unresolved_constraints(self):
        req = request()
        backend = AnalyticalBackend(req)
        problem = build_design_problem(req)
        spec = ComplexSpec((problem.generate(random.Random(1)).strands[0],))
        cache = ComplexAnalysisCache(backend)
        analysis = cache.analyze(spec, problem.conditions, need_pairs=True)
        self.assertIs(cache.partition(spec, problem.conditions), analysis.partition)
        self.assertIs(cache.partition(spec, problem.conditions,
                      PairingConstraints.unconstrained(spec.n_nucleotides)), analysis.partition)
        self.assertEqual(cache.statistics()["backend_calls"], 1)


class FinalistEvaluationCacheTests(unittest.TestCase):
    def test_excluded_finalist_cache_miss_does_not_rebuild_full_evaluation(self):
        req = request()
        backend = AnalyticalBackend(req, bad_on=True)
        cache_type = design_engine._EvaluationCache
        with patch.object(design_engine, "_EvaluationCache", side_effect=lambda **kwargs:
                          cache_type(max_entries=0, max_pair_cells=0)), \
                patch.object(design_engine, "evaluate_design", wraps=evaluate_design) as evaluate:
            batch = design_engine.design_candidates(req, backend, options=TrialOptions(
                max_evaluations=1, max_proposals=1))
        self.assertEqual(evaluate.call_count, 1)
        self.assertEqual(len(batch.candidates), 1)
        self.assertTrue(batch.candidates[0].excluded)
        self.assertEqual(batch.optimization_trials[0].candidate_screened, 1)

    def test_full_evaluation_cache_evicts_by_cells_before_entry_budget(self):
        req = request()
        problem, backend = build_design_problem(req), AnalyticalBackend(req)
        first = problem.generate(random.Random(1))
        second = problem.mutate_guided(first, random.Random(2), (), multi_fraction=0)
        evaluations = tuple(evaluate_design(problem, assignment, backend.analyze,
                            standard_state=backend.standard_state) for assignment in (first, second))
        cells = sum(a.partition.complex_spec.n_nucleotides ** 2
                    for a in evaluations[0].analyses if a.pairs is not None)
        cache = design_engine._EvaluationCache(max_entries=16, max_pair_cells=cells)
        for evaluation in evaluations:
            cache.remember(evaluation)
        self.assertIsNone(cache.get(first))
        self.assertIs(cache.get(second), evaluations[1])
        self.assertEqual(cache.pair_cells, cells)
        cache.remember(evaluations[1])
        self.assertEqual(cache.pair_cells, cells)
        too_small = design_engine._EvaluationCache(max_entries=16, max_pair_cells=cells - 1)
        too_small.remember(evaluations[0])
        self.assertIsNone(too_small.get(first))
        self.assertEqual(too_small.pair_cells, 0)


if __name__ == "__main__":
    unittest.main()
