"""Analytical OFF/ON tubes, coordinate integrity and full catalogue coverage."""
from dataclasses import replace
import math
import random
import unittest

from aptaswitch_core.thermodynamics.design_constraints import (
    DesignProblem, DesignTarget, DesignTube, StrandPattern,
)
from aptaswitch_core.thermodynamics.evaluation import ComplexAnalysis, evaluate_design
from aptaswitch_core.thermodynamics.structures import PairingConstraints, parse_structure
from aptaswitch_core.thermodynamics.tubes import StandardStateConvention
from aptaswitch_core.thermodynamics.types import (
    ComplexSpec, Conditions, ModelProvenance, PairProbabilities, PartitionResult,
)


MODEL = ModelProvenance("finite", "1", "analytical-test", "fixture:analytical",
                        "connected finite states", "fixture-1M")
STANDARD = StandardStateConvention("fixture-1M", 1., True, True, True)


def problem(max_size=2):
    return DesignProblem(
        strands=(StrandPattern("A", "rna", "GG"), StrandPattern("B", "dna", "CC")),
        targets=(DesignTarget("OFF", ("A",), ".."),
                 DesignTarget("ON", ("A", "B"), "((+))"),
                 DesignTarget("free_B", ("B",), "..")),
        tubes=(DesignTube("OFF", (("OFF", 1e-6),), max_size),
               DesignTube("ON", (("ON", 1e-6), ("free_B", 1e-6)), max_size)),
        conditions=Conditions(),
    )


class FiniteBackend:
    def __init__(self, empty_on=False):
        self.calls = []
        self.empty_on = empty_on

    def __call__(self, spec, conditions, *, need_pairs):
        self.calls.append((spec, conditions, need_pairs))
        # A=.6uM, B=1.6uM, AB=.4uM is an independently fixed equilibrium.
        names = tuple(s.name for s in spec.strands)
        if len(names) == 1:
            log_z = 0.0
            states = [("..", 1.)]
        elif names in (("A", "B"), ("B", "A")) and not self.empty_on:
            log_z = math.log(.4e-6 / (.6e-6 * 1.6e-6))
            states = [("((+))", .75), ("(.+.)", .25)]
        else:
            log_z, states = -math.inf, []
        partition = PartitionResult(spec, conditions, MODEL, log_z, -log_z, 1.)
        pairs = None
        if need_pairs and states:
            n = spec.n_nucleotides
            matrix = [[0.] * n for _ in range(n)]
            for structure, p in states:
                parsed = parse_structure(structure, [len(s.sequence) for s in spec.strands])
                for i, j in enumerate(parsed.partners):
                    matrix[i][j] += p
            pairs = PairProbabilities(partition, tuple(tuple(row) for row in matrix))
        return ComplexAnalysis(partition, pairs)


class DesignEvaluationTests(unittest.TestCase):
    def run_evaluation(self, p=None, backend=None):
        p = p or problem()
        backend = backend or FiniteBackend()
        return evaluate_design(p, p.generate(random.Random(0)), backend, standard_state=STANDARD)

    def test_analytical_multitube_loss_preserves_shortage_and_structural_parts(self):
        backend = FiniteBackend()
        result = self.run_evaluation(backend=backend)
        self.assertAlmostEqual(result.loss, (.2 + 2.4) / 6 / 2)
        self.assertAlmostEqual(result.defect.structural, .2 / 6 / 2)
        self.assertAlmostEqual(result.defect.concentration, 2.4 / 6 / 2)
        self.assertEqual([t.name for t in result.tubes], ["OFF", "ON"])
        self.assertEqual(len(backend.calls), 5)  # A/B/AA/AB/BB, A and AA shared.
        self.assertEqual(sum(need for _, _, need in backend.calls), 3)
        on = result.tubes[1].equilibrium
        self.assertAlmostEqual(on.total_concentrations_molar[0], 1e-6)
        self.assertAlmostEqual(on.total_concentrations_molar[1], 2e-6)
        for total, reconstructed in zip(on.total_concentrations_molar,
                                        on.stoichiometry @ on.concentrations_molar):
            self.assertAlmostEqual(reconstructed / total, 1., places=10)

    def test_complete_catalogue_to_five_includes_all_parasites_once(self):
        backend = FiniteBackend()
        result = self.run_evaluation(problem(5), backend)
        self.assertEqual([len(t.equilibrium.partitions) for t in result.tubes], [5, 23])
        self.assertEqual(len(backend.calls), 23)
        keys = [a.partition.complex_spec.cyclic_key for a in result.analyses]
        self.assertEqual(len(keys), len(set(keys)))
        self.assertEqual(sum(need for _, _, need in backend.calls), 3)

    def test_empty_on_ensemble_is_full_shortage_without_fabricated_pairs(self):
        result = self.run_evaluation(backend=FiniteBackend(empty_on=True))
        self.assertAlmostEqual(result.loss, (4 / 6) / 2)
        self.assertEqual(result.defect.structural, 0)
        self.assertIsNone(result.tubes[1].defect.contributions[0].complex_defect)

    def test_target_order_overrides_catalogue_without_rotating_its_matrix(self):
        p = problem()
        p = replace(p, targets=(p.targets[0], replace(p.targets[1], strands=("B", "A")), p.targets[2]))
        result = self.run_evaluation(p)
        target = result.tubes[1].defect.contributions[0]
        self.assertEqual([s.name for s in target.complex_defect.complex_spec.strands], ["B", "A"])
        self.assertAlmostEqual(result.loss, 2.6 / 6 / 2)

    def test_larger_explicit_target_included_beyond_size_cutoff(self):
        result = self.run_evaluation(problem(1))
        self.assertEqual([len(t.equilibrium.partitions) for t in result.tubes], [1, 3])
        self.assertAlmostEqual(result.loss, 2.6 / 6 / 2)

    def test_different_target_rotations_in_separate_tubes_are_not_reused(self):
        p = problem()
        p = replace(p, targets=p.targets + (DesignTarget("ON_BA", ("B", "A"), "((+))"),),
                    tubes=p.tubes + (DesignTube("reverse", (("ON_BA", 1e-6), ("free_B", 1e-6)), 2),))
        backend = FiniteBackend()
        result = self.run_evaluation(p, backend)
        self.assertEqual(len(backend.calls), 6)
        self.assertEqual(result.tubes[1].defect.normalized, result.tubes[2].defect.normalized)

    def test_wrong_backend_context_constraints_or_missing_matrix_fail(self):
        finite = FiniteBackend()
        def altered(kind):
            def analyze(spec, conditions, *, need_pairs):
                result = finite(spec, conditions, need_pairs=need_pairs)
                p = result.partition
                if kind == "matrix":
                    return ComplexAnalysis(p)
                if kind == "conditions":
                    p = replace(p, conditions=Conditions(celsius=25))
                elif kind == "constraints":
                    p = replace(p, constraints=PairingConstraints(tuple(range(spec.n_nucleotides))))
                elif kind == "order":
                    p = replace(p, complex_spec=ComplexSpec(tuple(reversed(spec.strands))))
                return ComplexAnalysis(p)
            return analyze
        for kind in ("matrix", "conditions", "constraints", "order"):
            with self.subTest(kind=kind), self.assertRaises(ValueError):
                self.run_evaluation(backend=altered(kind))

    def test_incompatible_model_between_tubes_is_not_silently_averaged(self):
        finite = FiniteBackend()
        def analyze(spec, conditions, *, need_pairs):
            result = finite(spec, conditions, need_pairs=need_pairs)
            if any(s.name == "B" for s in spec.strands):
                p = replace(result.partition, provenance=replace(MODEL, model="other"))
                return ComplexAnalysis(p, replace(result.pairs, partition=p) if result.pairs else None)
            return result
        with self.assertRaisesRegex(ValueError, "provenance"):
            self.run_evaluation(backend=analyze)
        # Even disjoint tubes must not conceal different thermodynamic models.
        p = problem()
        p = replace(p, tubes=(DesignTube("A", (("OFF", 1e-6),), 1),
                              DesignTube("B", (("free_B", 1e-6),), 1)))
        with self.assertRaisesRegex(ValueError, "provenance"):
            self.run_evaluation(p, analyze)

    def test_incomplete_problem_and_ambiguous_duplicate_targets_fail_before_backend(self):
        p = problem()
        invalid = [replace(p, conditions=None), replace(p, tubes=()),
                   replace(p, targets=p.targets + (replace(p.targets[1], name="otherON"),),
                           tubes=(DesignTube("two-targets", (("ON", 1e-6), ("otherON", 1e-6)), 2),))]
        backend = FiniteBackend()
        for entry in invalid:
            with self.subTest(entry=entry), self.assertRaises(ValueError):
                self.run_evaluation(entry, backend)
        self.assertEqual(backend.calls, [])


if __name__ == "__main__":
    unittest.main()
