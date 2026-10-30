"""Inverse conservation orchestration with an explicitly selected backend."""
from copy import deepcopy
from dataclasses import asdict, replace
import json
import math
import os
from pathlib import Path
import subprocess
import sys
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import numpy as np

from aptaswitch_core.inverse_extension import analyze_inverse_extensions, ranked_inverse_candidates
from aptaswitch_core.models import (
    DesignCandidate, InverseExtensionSettings, RunRequest, RunResult, SequenceInput,
    SwitchArchitecture, ThermoConditions,
)
from aptaswitch_core.thermodynamics.types import (
    Conditions, ModelProvenance, PairProbabilities, PartitionResult,
)
from aptaswitch_core.trigger_extension import design_trigger_extensions, padded_reference_structure


def candidate(trial, trigger="AAGT", score=100., excluded=False):
    return DesignCandidate(trial=trial, switch_seq="ACGU", trigger_seq=trigger, score=score,
        defect=.01, on_yield_pct=90, leak_pct=2, delta_g_off=-10, delta_g_on=-15,
        ddg_activation=-4, delta_g_rbs_linker=0, has_bad_rbs_linker=False,
        has_stop_codon=False, excluded=excluded)


def inverse_run(candidates=None, *, side="5prime", target=4, top=20,
                aptamer="ACGT", trigger="GT"):
    thermo = ThermoConditions(temperature_c=29.5, sodium_m=.15, magnesium_m=0,
                              switch_material="dna", trigger_material="dna")
    settings = InverseExtensionSettings(target_length=target, extension_side=side,
                                         top_candidates=top, thermo=thermo)
    request = RunRequest(SequenceInput(aptamer, trigger), SwitchArchitecture(12, 18),
                         inverse_extension=settings)
    return RunResult("native-inverse", request, candidates or [candidate(1)],
                     Path("unused"), "start", "end")


class FixtureBackend:
    """Deterministic orchestration fixture, not a physical energy model."""
    provenance = ModelProvenance("fixture", "1", "inverse-fixture", "fixture:1",
                                  "connected-fixture", "1M-fixture")
    condition_limits = dict(temperature_min_c=25., temperature_max_c=55., sodium_min_m=.12,
                            sodium_max_m=1., magnesium_min_m=0., magnesium_max_m=0.)

    def __init__(self):
        self.calls = []
        self.probability_calls = []
        self.fail_trigger = None

    def analyze(self, spec, conditions, *, need_pairs):
        self.calls.append((spec, conditions))
        if spec.strands[1].sequence == self.fail_trigger:
            raise ArithmeticError("later trigger failure")
        partition = PartitionResult(spec, conditions, self.provenance, 4., -2.4, .6)
        matrix = np.eye(spec.n_nucleotides)
        matrix[0, 0] = matrix[-1, -1] = .2
        matrix[0, -1] = matrix[-1, 0] = .8
        return SimpleNamespace(partition=partition, pairs=PairProbabilities(partition, matrix))

    def minimum(self, spec, conditions):
        aptamer, trigger = spec.strands
        structure = "(" + "." * (len(aptamer.sequence) - 1) + "+" + "." * (len(trigger.sequence) - 1) + ")"
        context = SimpleNamespace(complex_spec=spec, conditions=conditions,
            provenance=self.provenance, thermal_energy_kcal_mol=.6, constraints=None)
        return SimpleNamespace(context=context, structure=structure, energy_kcal_mol=-2.)

    def structure_probability(self, spec, conditions, structure, *, partition):
        self.probability_calls.append((spec, structure))
        return .6


class NativeInverseOrchestrationTests(unittest.TestCase):
    def test_reference_deduplication_metadata_and_scores_use_explicit_backend(self):
        backend = FixtureBackend()
        values = [candidate(1, "AAGT", 100), candidate(2, "AAGT", 90),
                  candidate(3, "CCGT", 80), candidate(4, "TTGT", 70, True)]
        run = inverse_run(values, top=3)
        before = [(id(c), c.score, c.score_raw, c.defect, c.trigger_seq, c.switch_seq) for c in values]
        with patch("aptaswitch_core.thermodynamics.setup.get_backend",
                   side_effect=AssertionError("Explicit backend must be retained")):
            analyze_inverse_extensions(run, backend=backend)
        stage = run.inverse_extension_analysis
        self.assertEqual(stage["status"], "completed", stage.get("error"))
        self.assertEqual(stage["selected_trials"], [1, 2, 3])
        self.assertEqual(stage["analyzed_count"], 3)
        self.assertEqual(stage["unique_trigger_count"], 2)
        self.assertEqual([spec.strands[1].sequence for spec, _ in backend.calls], ["GT", "AAGT", "CCGT"])
        self.assertTrue(all(conditions == Conditions(29.5, .15, 0) for _, conditions in backend.calls))
        self.assertTrue(all([s.material for s in spec.strands] == ["dna", "dna"] for spec, _ in backend.calls))
        self.assertEqual(stage["model"]["material"], backend.provenance.model)
        self.assertEqual(stage["model"]["ensemble"], backend.provenance.ensemble)
        self.assertEqual(stage["thermodynamics"], asdict(backend.provenance))

        self.assertEqual([(id(c), c.score, c.score_raw, c.defect, c.trigger_seq, c.switch_seq) for c in values], before)
        self.assertEqual(run.candidates[3].inverse_extension_metrics, {})
        for value in run.candidates[:3]:
            self.assertEqual(value.inverse_extension_metrics["thermodynamics"], stage["thermodynamics"])

        first, duplicate = [value.inverse_extension_metrics for value in run.candidates[:2]]
        self.assertEqual(first["final_score"], duplicate["final_score"])
        first["thermodynamics"]["version"] = "changed"
        self.assertEqual(duplicate["thermodynamics"], asdict(backend.provenance))
        self.assertEqual(stage["thermodynamics"], asdict(backend.provenance))
        self.assertEqual(len(ranked_inverse_candidates(run)), 3)
        json.dumps(stage, allow_nan=False)

    def test_each_extension_end_keeps_dna_and_reference_coordinates(self):
        for side, trigger in (("5prime", "AAGT"), ("3prime", "GTAA"), ("both", "AGTA")):
            with self.subTest(side=side):
                backend, run = FixtureBackend(), inverse_run([candidate(1, trigger)], side=side)
                analyze_inverse_extensions(run, backend=backend)
                stage = run.inverse_extension_analysis
                self.assertEqual(stage["status"], "completed", stage.get("error"))
                reference_spec = backend.calls[0][0]
                self.assertEqual([s.sequence for s in reference_spec.strands], ["ACGT", "GT"])
                self.assertEqual([s.material for s in reference_spec.strands], ["dna", "dna"])
                self.assertEqual(backend.calls[1][0].strands[1].sequence, trigger)
                target = padded_reference_structure(stage["reference"]["mfe_structure"], 2, side)
                self.assertIn((backend.calls[1][0], target), backend.probability_calls)

    def test_native_domain_failure_does_not_call_a_backend_or_change_switch_scores(self):
        for changes, constraint in (
            (dict(temperature_c=60), "25–55 °C"),
            (dict(sodium_m=.1), "Na⁺ 0.12–1 M"),
            (dict(magnesium_m=.001), "Mg²⁺ = 0"),
        ):
            with self.subTest(changes=changes):
                backend, run = FixtureBackend(), inverse_run()
                settings = run.request.inverse_extension
                run.request = replace(run.request, inverse_extension=replace(settings,
                    thermo=replace(settings.thermo, **changes)))
                with patch("aptaswitch_core.thermodynamics.setup.get_backend",
                           side_effect=AssertionError("Invalid conditions must not load a runtime")) as load:
                    analyze_inverse_extensions(run, backend=backend)
                load.assert_not_called()
                self.assertEqual(run.inverse_extension_analysis["status"], "failed")
                self.assertIn(constraint, run.inverse_extension_analysis["error"])
                self.assertEqual(backend.calls, [])
                self.assertEqual(run.candidates[0].score, 100.)
                self.assertEqual(run.candidates[0].inverse_extension_metrics, {})

    def test_other_initial_strand_chemistries_fail_before_thermodynamic_work(self):
        for aptamer, trigger in (("rna", "rna"), ("rna", "dna"), ("dna", "rna")):
            backend, run = FixtureBackend(), inverse_run()
            settings = run.request.inverse_extension
            run.request = replace(run.request, inverse_extension=replace(settings,
                thermo=replace(settings.thermo, switch_material=aptamer, trigger_material=trigger)))
            with self.subTest(aptamer=aptamer, trigger=trigger):
                analyze_inverse_extensions(run, backend=backend)
                self.assertEqual(run.inverse_extension_analysis["status"], "failed")
                self.assertEqual(backend.calls, [])
                self.assertEqual(run.candidates[0].inverse_extension_metrics, {})

    def test_later_failure_and_cancellation_preserve_only_completed_candidates(self):
        for failure in (False, True):
            with self.subTest(failure=failure):
                backend = FixtureBackend()
                run = inverse_run([candidate(1), candidate(2, "CCGT", 90)])
                if failure:
                    backend.fail_trigger = "CCGT"
                cancelled = False
                def progress(done, total, message):
                    nonlocal cancelled
                    if done == 2 and not failure:
                        cancelled = True
                analyze_inverse_extensions(run, backend=backend, progress_callback=progress,
                                             cancel_check=lambda: cancelled)
                stage = run.inverse_extension_analysis
                self.assertEqual(stage["status"], "failed" if failure else "cancelled")
                self.assertEqual(stage["analyzed_count"], 1)
                self.assertEqual(stage["unique_trigger_count"], 1)
                self.assertEqual(stage["candidates"][0]["trial"], 1)
                self.assertEqual(stage["thermodynamics"], asdict(backend.provenance))
                self.assertEqual(run.candidates[0].inverse_extension_metrics["status"], "completed")
                self.assertEqual(run.candidates[1].inverse_extension_metrics, {})

    def test_native_inverse_runs_in_a_fresh_process(self):
        code = '''
import sys
from test_thermodynamics_inverse_extension import FixtureBackend, inverse_run
from aptaswitch_core.inverse_extension import analyze_inverse_extensions
run = inverse_run()
analyze_inverse_extensions(run, backend=FixtureBackend())
assert run.inverse_extension_analysis['status'] == 'completed', run.inverse_extension_analysis
print('isolated')
'''
        process = subprocess.run([sys.executable, "-c", code], text=True, capture_output=True, timeout=30)
        self.assertEqual(process.returncode, 0, process.stderr)
        self.assertEqual(process.stdout.strip(), "isolated")


@unittest.skipUnless(os.environ.get("APTASWITCH_HYBRID_NATIVE_LIBRARY"), "Opt-in compiled native backend")
class RealNativeInverseTests(unittest.TestCase):
    def setUp(self):
        from aptaswitch_core.thermodynamics.native_backend import NativeHybridBackend
        self.backend = NativeHybridBackend(os.environ["APTASWITCH_HYBRID_NATIVE_LIBRARY"])

    def test_short_dna_extensions_complete_with_native_probabilities(self):
        run = inverse_run([candidate(1, "AACACGT"), candidate(2, "AACACGT", 90)],
                          aptamer="ACGTACGT", trigger="CACGT", target=7)
        analyze_inverse_extensions(run, backend=self.backend)
        stage = run.inverse_extension_analysis
        self.assertEqual(stage["status"], "completed", stage.get("error"))
        self.assertEqual(stage["analyzed_count"], 2)
        self.assertEqual(stage["unique_trigger_count"], 1)
        self.assertEqual(stage["thermodynamics"], asdict(self.backend.provenance))
        metrics = run.candidates[0].inverse_extension_metrics
        self.assertEqual(metrics["thermodynamics"], asdict(self.backend.provenance))
        self.assertEqual([len(part) for part in metrics["mfe_structure"].split("+")], [8, 7])
        self.assertTrue(math.isfinite(metrics["final_score"]))
        self.assertTrue(0 < metrics["mfe_structure_probability"] <= 1)
        self.assertTrue(0 <= metrics["target_structure_probability"] <= 1)


        json.dumps(stage, allow_nan=False)

    def test_one_direct_extension_and_inverse_conservation_metrics_agree(self):
        run = inverse_run(aptamer="ACGTACGT", trigger="CACGT", target=6)
        direct = design_trigger_extensions(aptamer_dna="ACGTACGT", trigger_dna="CACGT",
            target_lengths=[6], thermo=run.request.inverse_extension.thermo,
            top=1, max_ensemble_candidates=4, backend=self.backend)
        self.assertEqual(direct["status"], "completed", direct.get("error"))
        selected = direct["candidates"][0]
        run.candidates = [candidate(1, selected["extended_trigger"])]
        analyze_inverse_extensions(run, backend=self.backend)
        stage = run.inverse_extension_analysis
        self.assertEqual(stage["status"], "completed", stage.get("error"))
        actual = run.candidates[0].inverse_extension_metrics
        for name in ("mfe_structure", "mfe_energy_kcal_mol", "ensemble_free_energy_kcal_mol",
                     "mfe_structure_probability", "target_structure_probability", "final_score",
                     "reference_pair_probability_loss", "core_pair_probability_rmse",
                     "extension_pair_probability_sum"):
            self.assertEqual(actual[name], selected[name], name)

    def test_empty_connected_reference_fails_without_conservation_values(self):
        run = inverse_run([candidate(1, "AAAA")], aptamer="AAAA", trigger="AA")
        before = deepcopy(run.candidates[0].to_dict()) if hasattr(run.candidates[0], "to_dict") else asdict(run.candidates[0])
        analyze_inverse_extensions(run, backend=self.backend)
        stage = run.inverse_extension_analysis
        self.assertEqual(stage["status"], "failed")
        self.assertIn("empty", stage["error"])
        self.assertEqual(stage["reference"], {})
        self.assertEqual(stage["analyzed_count"], 0)
        self.assertEqual(run.candidates[0].inverse_extension_metrics, {})
        self.assertEqual(asdict(run.candidates[0]), before)


if __name__ == "__main__":
    unittest.main()
