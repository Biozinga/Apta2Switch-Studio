"""Physical objective reuse and constrained mutations for local-first design."""
from dataclasses import replace
import random
import unittest

from aptaswitch_core.models import FunctionalDefectWeights
from aptaswitch_core.thermodynamics.design_constraints import BaseSite, DesignProblem, StrandPattern, LibraryConstraint, UnsatisfiableDesignError, build_design_problem
from aptaswitch_core.thermodynamics.design_engine import design_candidates
from aptaswitch_core.thermodynamics.design_objective import objective_from_evaluation
from aptaswitch_core.thermodynamics.evaluation import ComplexAnalysis, ComplexAnalysisCache, evaluate_design
from aptaswitch_core.thermodynamics.optimization import TrialOptions
from aptaswitch_core.thermodynamics.types import ComplexSpec, PairProbabilities, Strand
from tests.test_thermodynamics_design_engine import AnalyticalBackend, request


class PriorityObjectiveTests(unittest.TestCase):
    def test_local_and_global_remain_distinct_and_no_extra_folding_is_needed(self):
        req = request()
        problem = build_design_problem(req)
        backend = AnalyticalBackend(req)
        evaluation = evaluate_design(problem, problem.generate(random.Random(1)), backend.analyze,
                                     standard_state=backend.standard_state)
        before = len(backend.calls)
        objective = objective_from_evaluation(req, evaluation)
        self.assertEqual(len(backend.calls), before)
        self.assertEqual(objective.local_defects, (0., 0., 0., 0.))
        self.assertEqual(objective.local_loss, 0.)
        self.assertGreater(objective.global_defect, 0.)
        self.assertGreater(objective.guard_defects[1], 0.)
        self.assertEqual(objective.guard_defects[2], 0.)
        self.assertAlmostEqual(objective.global_defect, evaluation.defect.normalized)

    def test_fixed_defective_toehold_guides_mutation_of_its_wrong_partner(self):
        req = replace(request(), functional_defect_weights=FunctionalDefectWeights(1., 0., 0., 0.))
        problem = build_design_problem(req)
        i, j = problem.layout['len_t7'], problem.layout['rbs_prefix_start']
        backend = AnalyticalBackend(req)
        analyze = backend.analyze
        def mispaired(spec, conditions, *, need_pairs):
            result = analyze(spec, conditions, need_pairs=need_pairs)
            if result.pairs and tuple(s.name for s in spec.strands) == ('Switch',):
                matrix = [list(row) for row in result.pairs.matrix]
                matrix[i][i] = matrix[j][j] = .75
                matrix[i][j] = matrix[j][i] = .25
                return ComplexAnalysis(result.partition, PairProbabilities(result.partition, tuple(map(tuple, matrix))))
            return result
        initial = problem.generate(random.Random(2))
        evaluation = evaluate_design(problem, initial, mispaired, standard_state=backend.standard_state)
        objective = objective_from_evaluation(req, evaluation)
        self.assertGreater(dict(objective.site_weights)[BaseSite('Switch', j)], 0.)
        mutated = problem.mutate_guided(initial, random.Random(1), objective.site_weights,
                                       exploration_fraction=0, multi_fraction=0)
        self.assertEqual(mutated.sequence('Switch')[i], initial.sequence('Switch')[i])
        self.assertNotEqual(mutated.sequence('Switch')[j], initial.sequence('Switch')[j])
        problem.validate(mutated)

    def test_impossible_stop_filter_spends_no_thermodynamic_evaluations(self):
        req = request()
        req = replace(req, exclude_stop_candidates=True, architecture=replace(req.architecture, frame_linker='UAA'))
        backend = AnalyticalBackend(req)
        with self.assertRaisesRegex(UnsatisfiableDesignError, 'STOP'):
            design_candidates(req, backend, options=TrialOptions(max_evaluations=10, max_proposals=40))
        self.assertEqual(backend.calls, [])


class ComplexReuseTests(unittest.TestCase):
    def test_unchanged_trigger_is_reused_between_distinct_candidates(self):
        req = request()
        problem = build_design_problem(req)
        backend = AnalyticalBackend(req)
        cache = ComplexAnalysisCache(backend)
        first = problem.generate(random.Random(1))
        second = problem.mutate_guided(first, random.Random(2), (), multi_fraction=0)
        for value in (first, second):
            evaluate_design(problem, value, cache, standard_state=backend.standard_state)
        trigger_calls = [call for call in backend.calls if call[1] == ('Trigger',)]
        self.assertEqual(len(trigger_calls), 1)
        self.assertGreater(cache.statistics()['cache_hits'], 0)
        self.assertLess(cache.statistics()['backend_calls'], 10)

    def test_cache_respects_conditions_probability_upgrade_and_provenance(self):
        req = request()
        backend = AnalyticalBackend(req)
        cache = ComplexAnalysisCache(backend)
        problem = build_design_problem(req)
        strand = problem.generate(random.Random(1)).strands[0]
        spec = ComplexSpec((strand,))
        for conditions in (problem.conditions, replace(problem.conditions, celsius=40)):
            cache(spec, conditions, need_pairs=False)
        cache(spec, problem.conditions, need_pairs=True)
        self.assertEqual(cache.statistics()['backend_calls'], 3)
        cache(spec, problem.conditions, need_pairs=True)
        self.assertEqual(cache.statistics()['cache_hits'], 1)
        backend.provenance = replace(backend.provenance, version='changed')
        with self.assertRaises(ValueError):
            cache(spec, problem.conditions, need_pairs=False)

    def test_identical_letters_do_not_merge_chemistries_labels_or_strand_order(self):
        req = request()
        backend = AnalyticalBackend(req)
        cache = ComplexAnalysisCache(backend)
        conditions = build_design_problem(req).conditions
        dna = Strand('ACG', 'dna', 'A')
        rna = Strand('ACG', 'rna', 'A')
        other = Strand('GCA', 'rna', 'B')
        contexts = (ComplexSpec((dna,)), ComplexSpec((rna,)),
                    ComplexSpec((replace(dna, name='renamed'),)),
                    ComplexSpec((dna, other)), ComplexSpec((other, dna)))
        for spec in contexts:
            result = cache(spec, conditions, need_pairs=True)
            self.assertEqual(result.partition.complex_spec, spec)
        self.assertEqual(cache.statistics()['backend_calls'], len(contexts))
        for spec in contexts:
            cache(spec, conditions, need_pairs=True)
        self.assertEqual(cache.statistics()['backend_calls'], len(contexts))
        self.assertEqual(cache.statistics()['cache_hits'], len(contexts))

    def test_probability_memory_is_bounded(self):
        req = request()
        backend = AnalyticalBackend(req)
        problem = build_design_problem(req)
        spec = ComplexSpec((problem.generate(random.Random(1)).strands[0],))
        cache = ComplexAnalysisCache(backend, max_pair_cells=1)
        cache(spec, problem.conditions, need_pairs=True)
        self.assertEqual(cache.statistics()['cached_pair_cells'], 0)
        self.assertEqual(cache.statistics()['cached_complexes'], 0)


class CoupledPopulationTests(unittest.TestCase):
    def problem(self):
        return DesignProblem((StrandPattern('A', 'rna', 'X'), StrandPattern('B', 'dna', 'X')),
            (LibraryConstraint('hybrid', (BaseSite('A', 0), BaseSite('B', 0)), ('AT', 'GC', 'CG', 'UA')),))

    def test_enumeration_and_distance_count_coupled_groups_once(self):
        problem = self.problem()
        values = problem.enumerate_assignments(4)
        self.assertEqual(len(values), 4)
        self.assertIsNone(problem.enumerate_assignments(3))
        self.assertEqual(len(problem.mutable_components()), 1)
        self.assertEqual(problem.design_distance(values[0], values[1]), 1.)
        for value in values:
            other = problem.mutate_guided(value, random.Random(1), ())
            problem.validate(other)
            self.assertNotEqual(value, other)

    def test_two_coupled_groups_can_change_together_without_breaking_constraints(self):
        problem = DesignProblem((StrandPattern('S', 'rna', 'XX'),))
        first = problem.generate(random.Random(1))
        other = problem.mutate_guided(first, random.Random(2), (), multi_fraction=1)
        self.assertEqual(problem.design_distance(first, other), 1.)
        self.assertEqual(sum(a != b for a,b in zip(first.sequence('S'), other.sequence('S'))), 2)


if __name__ == '__main__':
    unittest.main()
