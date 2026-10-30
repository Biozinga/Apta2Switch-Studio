"""Existing candidate rules exercised with an explicit analytical test backend."""
from dataclasses import replace
import math
import os
import random
from types import SimpleNamespace
import unittest

from aptaswitch_core.architecture import default_architecture_for_trigger, target_switch_structures
from aptaswitch_core.models import InverseExtensionSettings, RunRequest, SequenceInput, ThermoConditions
from aptaswitch_core.thermodynamics.design_constraints import build_design_problem
from aptaswitch_core.thermodynamics.design_engine import design_candidates, evaluate_candidate
from aptaswitch_core.thermodynamics.evaluation import ComplexAnalysis, evaluate_design
from aptaswitch_core.thermodynamics.optimization import TrialOptions
from aptaswitch_core.thermodynamics.structures import format_structure, parse_structure
from aptaswitch_core.thermodynamics.tubes import StandardStateConvention
from aptaswitch_core.thermodynamics.types import ComplexSpec, ModelProvenance, PairProbabilities, PartitionResult


def request(**changes):
    return RunRequest(SequenceInput("AACCGGTTAACC", "ACGT"*7+"AC"),
                      default_architecture_for_trigger(30), trials=1,
                      thermo=ThermoConditions(validation_max_complex_size=2),
                      exclude_stop_candidates=False, **changes)


class AnalyticalBackend:
    """Contract fixture: prescribed connected weights, not biological predictions."""
    provenance = ModelProvenance("finite-test", "1", "fixture", "fixture:design",
                                 "connected test states", "fixture-1M")
    standard_state = StandardStateConvention("fixture-1M", 1., True, True, True)

    def __init__(self, req, *, bad_on=False):
        self.request, self.bad_on, self.calls = req, bad_on, []
        self.off, self.on = target_switch_structures(req.architecture, len(req.reporter.sequence_after_start_rna))

    def structure(self, spec):
        names = tuple(s.name for s in spec.strands)
        if names == ("Switch",):
            return self.off
        if names == ("Trigger", "Switch"):
            return self.on
        return "+".join("."*len(s.sequence) for s in spec.strands)

    def analyze(self, spec, conditions, *, need_pairs):
        self.calls.append(("analyze", tuple(s.name for s in spec.strands), need_pairs))
        names = tuple(s.name for s in spec.strands)
        if len(names) == 1:
            log_z = 0.
        elif sorted(names) == ["Switch", "Trigger"]:
            log_z = math.log(1e6)
        elif sorted(names) == ["Aptamer", "Trigger"]:
            log_z = math.log(1e9)
        else:
            log_z = -math.inf
        p = PartitionResult(spec, conditions, self.provenance, log_z, -log_z, 1.)
        pairs = None
        if need_pairs and not p.is_empty:
            n = spec.n_nucleotides
            matrix = [[0.]*n for _ in range(n)]
            target = parse_structure(self.structure(spec), [len(s.sequence) for s in spec.strands])
            for i, j in enumerate(target.partners):
                matrix[i][j] = 1.
            if names == ("Trigger",) and self.request.inverse_extension is not None:
                for i in range(n):
                    matrix[i][i], matrix[i][n-1-i] = .8, .2
            pairs = PairProbabilities(p, tuple(tuple(row) for row in matrix))
        return ComplexAnalysis(p, pairs)

    def partition(self, spec, conditions, constraints=None):
        partition = self.analyze(spec, conditions, need_pairs=False).partition
        if constraints is None:
            return partition
        if len(constraints.partners) != spec.n_nucleotides:
            raise ValueError("Constraints must describe the complete fixture complex")
        state = parse_structure(self.structure(spec), [len(s.sequence) for s in spec.strands])
        allowed = all(wanted < 0 or state.partners[i] == wanted
                      for i, wanted in enumerate(constraints.partners))
        return replace(partition, constraints=constraints,
                       log_partition=partition.log_partition if allowed else -math.inf,
                       free_energy_kcal_mol=partition.free_energy_kcal_mol if allowed else math.inf)

    def minimum(self, spec, conditions):
        names = tuple(s.name for s in spec.strands)
        self.calls.append(("minimum", names))
        energy = {("Switch",): -10., ("Trigger", "Switch"): -25.,
                  ("Trigger",): -2., ("RBS_linker",): -1.}[names]
        return SimpleNamespace(structure=self.structure(spec), energy_kcal_mol=energy)

    def screen_unpaired_region(self, spec, conditions, *, start, end):
        self.calls.append(("screen", start, end))
        witness = None
        if self.bad_on:
            lengths = [len(s.sequence) for s in spec.strands]
            state = list(parse_structure(self.on, lengths).partners)
            state[start], state[start+4] = start+4, start
            from aptaswitch_core.thermodynamics.structures import PairingConstraints
            witness = SimpleNamespace(structure=format_structure(PairingConstraints(tuple(state)), lengths),
                                      energy_kcal_mol=-25.)
        return SimpleNamespace(all_minima_unpaired=not self.bad_on, witness=witness)


class NativeDesignOrchestrationTests(unittest.TestCase):
    def test_direct_engine_entry_rejects_other_strand_roles_before_work(self):
        req = request()
        backend = AnalyticalBackend(req)
        for switch, trigger in (("dna", "rna"), ("rna", "rna"), ("dna", "dna")):
            invalid = replace(req, thermo=replace(req.thermo,
                              switch_material=switch, trigger_material=trigger))
            with self.subTest(switch=switch, trigger=trigger), self.assertRaises(ValueError):
                design_candidates(invalid, backend)
        self.assertEqual(backend.calls, [])

    def candidate(self, req, backend=None, seed=20):
        problem = build_design_problem(req)
        backend = backend or AnalyticalBackend(req)
        candidate = evaluate_candidate(req, problem, problem.generate(random.Random(seed)), backend, trial=1)
        return candidate, backend

    def test_keeps_existing_defect_score_and_chemical_order_semantics(self):
        req = request()
        candidate, backend = self.candidate(req)
        self.assertFalse(candidate.excluded)
        self.assertEqual(candidate.ddg_activation, -13.)
        self.assertEqual(candidate.delta_g_rbs_linker, -1.)
        self.assertEqual(candidate.regional_defects["balanced"], 0.)
        self.assertEqual(candidate.defect_structural, 0.)
        self.assertGreater(candidate.defect_concentration, 0.)
        self.assertGreater(candidate.on_yield_pct, candidate.leak_pct)
        self.assertTrue(0 <= candidate.score <= 100)
        weights = req.scoring
        raw = (weights.on_yield*candidate.on_yield_pct-weights.leak*candidate.leak_pct
               -weights.ddg*candidate.ddg_activation+weights.rbs*candidate.delta_g_rbs_linker
               +100*weights.rbs_calculator*candidate.rbs_calculator['merit'])
        self.assertAlmostEqual(candidate.score_raw, raw, places=4)
        p = build_design_problem(req)
        self.assertIn(("screen", 30+p.layout["rbs_linker_start"], 30+p.layout["rbs_linker_end"]), backend.calls)
        self.assertEqual(candidate.trigger_linearity, {})
        self.assertEqual(len(candidate.switch_seq.split("+")), 1)

    def test_bad_on_stops_before_other_analyses_and_keeps_witness(self):
        req = request()
        backend = AnalyticalBackend(req, bad_on=True)
        candidate, _ = self.candidate(req, backend)
        self.assertTrue(candidate.excluded)
        self.assertIsNone(candidate.score)
        self.assertIsNone(candidate.defect)
        self.assertEqual(candidate.delta_g_on, -25.)
        self.assertEqual(len(backend.calls), 1)
        self.assertEqual(backend.calls[0][0], "screen")
        self.assertIn("+", candidate.structure_on)

    def test_stop_constraints_generate_admissible_sequence_with_existing_fixed_trigger(self):
        req = replace(request(), exclude_stop_candidates=True,
                      sequences=SequenceInput("AACCGGTTAACC", "GCTGATTCTCAATGGCTCGCCCGATATGCA"))
        problem = build_design_problem(req)
        assignment = problem.generate(random.Random(3))
        self.assertFalse(problem.reading_frame_report(assignment)["has_premature_stop"])
        self.assertEqual(assignment.sequence("Trigger"), req.sequences.trigger_dna)
        problem.validate(assignment)

    def test_inverse_trigger_linearity_has_existing_weight_and_real_backbone(self):
        req = replace(request(), sequences=SequenceInput("AACCGGTT", "ACGT"*6),
                      inverse_extension=InverseExtensionSettings(30, "both"))
        candidate, _ = self.candidate(req)
        self.assertAlmostEqual(candidate.trigger_linearity["normalized_defect"], .2)
        self.assertEqual(candidate.trigger_linearity["nucleotide_count"], 30)
        self.assertNotIn("U", candidate.trigger_seq)
        self.assertNotIn("T", candidate.switch_seq)
        without = replace(req, scoring=replace(req.scoring, trigger_linearity=0))
        other, _ = self.candidate(without)
        self.assertAlmostEqual(other.score_raw-candidate.score_raw, req.scoring.trigger_linearity*20, places=4)

    def test_budget_outcome_and_independent_seeds_are_retained(self):
        req = replace(request(), trials=2)
        options = TrialOptions(seed=17, max_evaluations=1, max_proposals=1)
        result = design_candidates(req, AnalyticalBackend(req), options=options)
        self.assertEqual(result.status, "completed")
        self.assertEqual(len(result.candidates), 2)
        self.assertEqual([r.options.seed for r in result.optimization_trials], [17,18])
        self.assertEqual([r.stop_reason for r in result.optimization_trials], ["budget_exhausted"]*2)
        self.assertEqual(result.provenance, AnalyticalBackend.provenance)
        self.assertEqual(result.conditions.celsius, req.thermo.temperature_c)

    def test_cancel_before_start_is_explicit_not_a_completed_empty_run(self):
        req = request()
        backend = AnalyticalBackend(req)
        result = design_candidates(req, backend, cancel_check=lambda: True)
        self.assertEqual(result.status, "cancelled")
        self.assertEqual(result.optimization_trials, ())
        self.assertEqual(backend.calls, [])

    def test_cancel_during_final_candidate_analysis_preserves_complete_candidate(self):
        req = request()
        backend = AnalyticalBackend(req)
        cancelled = False
        original_minimum = backend.minimum
        def minimum(*args, **kwargs):
            nonlocal cancelled
            if tuple(s.name for s in args[0].strands) == ("Trigger", "Switch"):
                cancelled = True
            return original_minimum(*args, **kwargs)
        backend.minimum = minimum
        batch = design_candidates(req, backend, options=TrialOptions(
            max_evaluations=1, max_proposals=1), cancel_check=lambda: cancelled)
        self.assertEqual(batch.status, "cancelled")
        self.assertEqual(len(batch.candidates), 1)
        self.assertIsNotNone(batch.candidates[0].score)
        self.assertEqual(batch.optimization_trials[0].stop_reason, "budget_exhausted")

    def test_cached_loss_from_other_targets_is_rejected(self):
        req = request()
        problem = build_design_problem(req)
        assignment = problem.generate(random.Random(20))
        backend = AnalyticalBackend(req)
        other = replace(problem, tubes=problem.tubes[:1])
        evaluated = evaluate_design(other, assignment, backend.analyze, standard_state=backend.standard_state)
        with self.assertRaisesRegex(ValueError, "different candidate"):
            evaluate_candidate(req, problem, assignment, backend, trial=1, evaluated=evaluated)

    def test_chemically_identical_aptamer_and_trigger_do_not_duplicate_species(self):
        req = request()
        req = replace(req, sequences=replace(req.sequences, aptamer_dna=req.sequences.trigger_dna))
        candidate, backend = self.candidate(req)
        self.assertFalse(candidate.excluded)
        # DNA origins are equiprobable in the mixed pool; the reported target
        # counts trigger-origin duplexes. Extra aptamer-origin copies change
        # the equilibrium before this one-half labelling factor is applied.
        total = req.thermo.resolved_trigger_concentration_m
        switch = req.thermo.resolved_switch_concentration_m
        k = 1e6
        physical_duplex = 2*(2*total)*switch / (2*total+switch+1/k +
            math.sqrt((2*total+switch+1/k)**2-4*(2*total)*switch))
        self.assertAlmostEqual(candidate.leak_pct, 50*physical_duplex/switch, places=4)


@unittest.skipUnless(os.environ.get("APTASWITCH_HYBRID_NATIVE_LIBRARY"), "Explicit compiled hybrid library required")
class NativeClassicSwitchIntegrationTests(unittest.TestCase):
    def test_inverse_x_domains_resolve_before_real_native_calculations(self):
        from aptaswitch_core.thermodynamics.native_backend import NativeHybridBackend
        backend = NativeHybridBackend(os.environ["APTASWITCH_HYBRID_NATIVE_LIBRARY"])
        core = "ACGT" * 6
        req = replace(request(exclude_non_linear_candidates=False),
                      sequences=SequenceInput("AACCGGTT", core),
                      inverse_extension=InverseExtensionSettings(30, "both"))
        self.assertEqual(req.design_trigger_pattern(), "XXX" + core + "XXX")
        batch = design_candidates(req, backend, options=TrialOptions(
            seed=17, max_evaluations=1, max_proposals=1))
        candidate, = batch.candidates
        self.assertEqual(len(candidate.trigger_seq), 30)
        self.assertEqual(candidate.trigger_seq[3:-3], core)
        self.assertLessEqual(set(candidate.trigger_seq), set("ACGT"))
        self.assertLessEqual(set(candidate.switch_seq), set("ACGU"))
        self.assertEqual(candidate.trigger_linearity["nucleotide_count"], 30)
        self.assertTrue(math.isfinite(candidate.defect))

    def test_classic30_reports_real_metrics_without_claiming_optimization_convergence(self):
        from aptaswitch_core.thermodynamics.native_backend import NativeHybridBackend
        backend = NativeHybridBackend(os.environ["APTASWITCH_HYBRID_NATIVE_LIBRARY"])
        # Dimers keep this integration check short. Actual request defaults
        # remain max_size=3; screening is disabled here to exercise all metrics.
        req = request(exclude_non_linear_candidates=False)
        batch = design_candidates(req, backend, options=TrialOptions(
            seed=17, max_evaluations=1, max_proposals=1))
        candidate, = batch.candidates
        outcome, = batch.optimization_trials
        self.assertEqual(outcome.stop_reason, "budget_exhausted")
        self.assertGreater(max(outcome.best_loss, outcome.best_global_loss),
                           max(outcome.options.target_local_loss, outcome.options.target_global_loss))
        problem = build_design_problem(req)
        problem.validate(outcome.best_assignment)
        self.assertTrue(any("X" in strand.pattern for strand in problem.strands))
        strands = {s.name: s for s in outcome.best_assignment.strands}
        self.assertEqual(strands["Trigger"].material, "dna")
        self.assertEqual(strands["Switch"].material, "rna")
        self.assertEqual(candidate.trigger_seq, strands["Trigger"].sequence)
        self.assertEqual(candidate.switch_seq, strands["Switch"].sequence)
        self.assertLessEqual(set(candidate.switch_seq), set("ACGU"))
        self.assertLessEqual(set(candidate.trigger_seq), set("ACGT"))
        self.assertEqual(len(candidate.trigger_seq), 30)
        self.assertEqual(len(candidate.switch_seq), 113)
        self.assertFalse(candidate.excluded)
        for structure, lengths in ((candidate.structure_off, (113,)),
                                   (candidate.structure_on, (30, 113))):
            parse_structure(structure, lengths)
        self.assertTrue(0 <= candidate.score <= 100)
        self.assertAlmostEqual(candidate.defect, candidate.defect_structural + candidate.defect_concentration, places=5)
        self.assertTrue(0 <= candidate.regional_defects["balanced"] <= 1)
        self.assertTrue(0 <= candidate.on_yield_pct <= 100)
        self.assertTrue(0 <= candidate.leak_pct <= 100)
        self.assertAlmostEqual(candidate.ddg_activation,
            candidate.delta_g_on - candidate.delta_g_off -
            backend.minimum(ComplexSpec((strands["Trigger"],)),
                            problem.conditions).energy_kcal_mol, places=5)
        self.assertEqual(batch.provenance, backend.provenance)


if __name__ == "__main__":
    unittest.main()
