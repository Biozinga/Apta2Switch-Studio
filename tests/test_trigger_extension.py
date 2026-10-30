"""Algorithm parity, genuine enumeration, cancellation and native integration."""

from __future__ import annotations

import json
import os
from pathlib import Path
import sys
from types import SimpleNamespace
from dataclasses import asdict, replace
import unittest
from unittest.mock import patch

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from aptaswitch_core import trigger_extension as extension
from aptaswitch_core.models import ThermoConditions


def original_reference():
    """Golden outputs evaluated by the original program, not this app's engine."""
    return json.loads((ROOT / "tests/fixtures/extension_original_reference.json").read_text())


def pair_matrix(size, pairs):
    matrix = np.eye(size)
    for (i, j), probability in pairs.items():
        matrix[i, j] = matrix[j, i] = probability
        matrix[i, i] -= probability
        matrix[j, j] -= probability
    return matrix


class AlgorithmParityTests(unittest.TestCase):
    def setUp(self):
        self.original = original_reference()

    def test_mfe_keys_and_values_match_source_with_generic_sequences(self):
        reference = extension.MfeResult("((+))", -3.2)
        for case in self.original["mfe_cases"]:
            candidate_structure = case["structure"]
            with self.subTest(structure=candidate_structure):
                candidate_mfe = extension.MfeResult(candidate_structure, -3.1)
                mapped = extension.map_reference_pairs_to_extended({(0, 3), (1, 2)}, 2, 2)
                candidate = extension.score_mfe_candidate("AC", candidate_mfe, reference, mapped, 2, "GT")
                actual = extension.asdict(candidate)
                self.assertEqual(actual.pop("extension_3prime"), "")
                actual["extended_z1"] = actual.pop("extended_trigger")
                self.assertEqual(actual, case["candidate"])
                self.assertEqual(list(extension.mfe_key(candidate)), case["rank_key"])

    def test_ensemble_formula_and_tie_breaks_match_original(self):
        reference = extension.ComplexAnalysisResult("((+))", -3.2, -3.5, .6, .6,
                                                     pair_matrix(4, {(0, 3): .8, (1, 2): .7}))
        analysis = extension.ComplexAnalysisResult("((+.))", -3.1, -3.45, .55, .55,
                                                    pair_matrix(5, {(0, 4): .75, (1, 3): .65, (0, 2): .05}))
        mfe = extension.MfeResult(analysis.mfe_structure, analysis.mfe_energy_kcal_mol)
        reference_mfe = extension.MfeResult(reference.mfe_structure, reference.mfe_energy_kcal_mol)
        mapped = extension.map_reference_pairs_to_extended({(0, 3), (1, 2)}, 2, 1)
        candidate = extension.score_mfe_candidate("A", mfe, reference_mfe, mapped, 2, "GT")
        actual = extension.ensemble_metrics(candidate, reference, analysis, 2, 2)
        expected = self.original["ensemble_case"]
        renamed = dict(actual)
        self.assertEqual(renamed.pop("extension_3prime"), "")
        renamed["extended_z1"] = renamed.pop("extended_trigger")
        renamed["extension_to_z0_probability_sum"] = renamed.pop("extension_to_aptamer_probability_sum")
        renamed["extension_to_z1_core_probability_sum"] = renamed.pop("extension_to_trigger_core_probability_sum")
        self.assertEqual(renamed, expected["candidate"])
        self.assertEqual(list(extension.final_rank_key(actual)), expected["rank_key"])

    def test_search_estimate_is_full_space_and_validates_target_lengths(self):
        estimate = extension.estimate_extension_search(18, [30, 24, 24])
        self.assertEqual(estimate["extension_lengths"], [6, 12])
        self.assertEqual(estimate["total_candidates_expected"], 4**6 + 4**12)
        for lengths in ([], [18], [17], [True], [24.0]):
            with self.subTest(lengths=lengths), self.assertRaises(ValueError):
                extension.estimate_extension_search(18, lengths)


class SyntheticExtensionModel:
    """Synthetic complete model outputs for the search orchestration boundary."""
    def __init__(self, backend):
        self.metadata = {"thermodynamics": asdict(backend.provenance)}
        self.mfe_sequences, self.ensemble_sequences = [], []
        self.fail_mfe = False

    def minimum(self, aptamer, trigger):
        self.mfe_sequences.append((aptamer, trigger))
        if self.fail_mfe:
            raise RuntimeError("Native calculation failure")
        return extension.MfeResult("." * len(aptamer) + "+" + "." * len(trigger), -1.)

    def analyze(self, aptamer, trigger, *, target=None):
        self.ensemble_sequences.append((aptamer, trigger))
        structure = "." * len(aptamer) + "+" + "." * len(trigger)
        return extension.ComplexAnalysisResult(structure, -1., -1.2, 1., 1., np.eye(len(aptamer) + len(trigger)))


class SearchOrchestrationTests(unittest.TestCase):
    def setUp(self):
        from aptaswitch_core.localization import language_context
        from test_thermodynamics_context_analysis import FixedWeightBackend
        self.enterContext(language_context("fr"))
        self.backend = FixedWeightBackend()
        self.fake = SyntheticExtensionModel(self.backend)
        self.enterContext(patch("aptaswitch_core.thermodynamics.setup.get_backend", return_value=self.backend))
        self.import_context = self.enterContext(patch.object(extension, "open_extension_model"))
        self.import_context.return_value.__enter__.return_value = self.fake

    def run_search(self, **kwargs):
        options = dict(aptamer_dna="AC", trigger_dna="GT", target_lengths=[3, 4],
                       thermo=ThermoConditions(temperature_c=37, switch_material="dna"), top=2)
        options.update(kwargs)
        return extension.design_trigger_extensions(**options)

    def test_enumerates_every_sequence_and_preserves_top_for_every_length(self):
        result = self.run_search()
        self.assertEqual(result["status"], "completed")
        self.assertEqual(len(self.fake.mfe_sequences), 20)
        self.assertEqual(len(set(self.fake.mfe_sequences)), 20)
        self.assertTrue(all(aptamer == "AC" and trigger.endswith("GT") for aptamer, trigger in self.fake.mfe_sequences))
        self.assertEqual([candidate["trigger_length"] for candidate in result["candidates"]], [3, 3, 4, 4])
        self.assertEqual([candidate["rank_within_length"] for candidate in result["candidates"]], [1, 2, 1, 2])
        self.assertTrue(result["mfe_screen"]["exhaustive"])
        self.assertTrue(result["mfe_screen"]["ensemble_ranking_exhaustive"])
        self.assertIn("nucleotide_probabilities", result["reference"])
        json.dumps(result, allow_nan=False)

    def test_cap_limits_only_ensemble_and_is_declared(self):
        result = self.run_search(max_ensemble_candidates=2)
        self.assertEqual(len(self.fake.mfe_sequences), 20)
        self.assertEqual(len(self.fake.ensemble_sequences), 5)  # reference + 2/length
        self.assertTrue(result["mfe_screen"]["exhaustive"])
        self.assertTrue(result["mfe_screen"]["mfe_survivor_cap_hit"])
        self.assertFalse(result["mfe_screen"]["ensemble_ranking_exhaustive"])
        self.assertTrue(any("classement d'ensemble est partiel" in warning for warning in result["warnings"]))

    def test_cancel_during_screen_returns_no_pretend_candidates(self):
        result = self.run_search(cancel_check=lambda: len(self.fake.mfe_sequences) >= 3)
        self.assertEqual(result["status"], "cancelled")
        self.assertEqual(result["mfe_screen"]["total_candidates_screened"], 3)
        self.assertFalse(result["mfe_screen"]["exhaustive"])
        self.assertEqual(result["candidates"], [])

    def test_cancel_during_ensemble_keeps_explicitly_partial_results(self):
        result = self.run_search(cancel_check=lambda: len(self.fake.ensemble_sequences) >= 3)
        self.assertEqual(result["status"], "cancelled")
        self.assertTrue(result["mfe_screen"]["exhaustive"])
        self.assertEqual(len(result["candidates"]), 2)
        self.assertFalse(result["mfe_screen"]["ensemble_complete"])

    def test_advance_stops_screen_but_runs_ensemble_and_preserves_survivors(self):
        phases = []
        checkpoints = []
        result = self.run_search(advance_check=lambda: len(self.fake.mfe_sequences) >= 3,
                                 phase_callback=phases.append, checkpoint_callback=checkpoints.append)
        self.assertEqual(result["status"], "completed")
        self.assertEqual(result["mfe_screen"]["total_candidates_screened"], 3)
        self.assertTrue(result["mfe_screen"]["stopped_early"])
        self.assertTrue(result["mfe_screen"]["ensemble_complete"])
        self.assertFalse(result["mfe_screen"]["ensemble_ranking_exhaustive"])
        self.assertEqual(len(result["mfe_survivors"]), 3)
        self.assertEqual(len(self.fake.ensemble_sequences), 4)
        self.assertEqual(phases, ["reference", "mfe", "ensemble"])
        self.assertTrue(checkpoints)

    def test_cancelled_screen_resumes_ensemble_without_any_mfe_replay(self):
        saved = self.run_search(cancel_check=lambda: len(self.fake.mfe_sequences) >= 3)
        before = len(self.fake.mfe_sequences)
        result = self.run_search(resume_result=saved)
        self.assertEqual(len(self.fake.mfe_sequences), before)
        self.assertEqual(result["status"], "completed")
        self.assertEqual(result["mfe_screen"]["ensemble_candidates_analyzed"], 3)
        self.assertEqual(result["resume"]["method"], "saved_survivors")
        self.assertEqual(saved["status"], "cancelled")

    def test_legacy_recovery_only_replays_first_retained_ties(self):
        saved = self.run_search(max_ensemble_candidates=2)
        saved.pop("mfe_survivors")
        self.fake.mfe_sequences.clear()
        recovered = self.run_search(max_ensemble_candidates=2, resume_result=saved)
        self.assertEqual(len(self.fake.mfe_sequences), 4)  # 2 per length, not all 20
        self.assertEqual(recovered["candidates"], saved["candidates"])
        self.assertEqual(recovered["mfe_screen"]["total_candidates_screened"], 20)

    def test_resume_rejects_different_inputs_or_reference(self):
        saved = self.run_search()
        with self.assertRaisesRegex(ValueError, "mêmes séquences"):
            self.run_search(resume_result=saved, trigger_dna="AA")
        with self.assertRaisesRegex(ValueError, "mêmes séquences"):
            self.run_search(resume_result=saved, thermo=ThermoConditions(material="experimental-Sugimoto95-Vienna04-loopmean", temperature_c=25, switch_material="dna"))
        saved["reference"]["mfe_energy_kcal_mol"] = -42
        with self.assertRaisesRegex(ValueError, "référence thermodynamique diffère"):
            self.run_search(resume_result=saved)

    def test_cancel_takes_precedence_over_advance(self):
        result = self.run_search(cancel_check=lambda: len(self.fake.mfe_sequences) >= 3,
                                 advance_check=lambda: len(self.fake.mfe_sequences) >= 3)
        self.assertEqual(result["status"], "cancelled")
        self.assertEqual(result["mfe_screen"]["ensemble_candidates_analyzed"], 0)

    def test_cancelling_legacy_recovery_keeps_recoverable_recipe(self):
        saved = self.run_search(max_ensemble_candidates=2)
        saved.pop("mfe_survivors")
        self.fake.mfe_sequences.clear()
        stopped = self.run_search(max_ensemble_candidates=2, resume_result=saved,
                                  cancel_check=lambda: len(self.fake.mfe_sequences) >= 1)
        self.assertEqual(stopped["status"], "cancelled")
        self.assertNotIn("mfe_survivors", stopped)
        result = self.run_search(max_ensemble_candidates=2, resume_result=stopped)
        self.assertEqual(result["candidates"], saved["candidates"])

    def test_periodic_checkpoint_resumes_a_crashed_screen(self):
        import itertools
        from copy import deepcopy
        checkpoints = []
        def crash_after_checkpoint(done, total, message):
            if len(self.fake.mfe_sequences) >= 3:
                raise RuntimeError("interrupted process")
        with patch.object(extension.time, "monotonic", side_effect=itertools.count(0, 61)):
            with self.assertRaisesRegex(RuntimeError, "interrupted process"):
                self.run_search(checkpoint_callback=lambda data: checkpoints.append(deepcopy(data)),
                                progress_callback=crash_after_checkpoint)
        saved = checkpoints[-1]
        before = len(self.fake.mfe_sequences)
        result = self.run_search(resume_result=saved)
        self.assertEqual(len(self.fake.mfe_sequences), before)
        self.assertEqual(result["mfe_screen"]["ensemble_candidates_analyzed"], len(saved["mfe_survivors"]))
        self.assertFalse(result["mfe_screen"]["exhaustive"])

    def test_resume_rejects_incomplete_survivor_checkpoint(self):
        saved = self.run_search()
        saved["mfe_survivors"].pop()
        with self.assertRaisesRegex(ValueError, "Sauvegarde MFE incohérente"):
            self.run_search(resume_result=saved)

    def test_model_honors_exact_conditions_and_dna_chemistry(self):
        result = self.run_search(thermo=ThermoConditions(switch_material="dna", trigger_material="dna", temperature_c=25, sodium_m=.22, magnesium_m=0))
        self.assertEqual(result["model"]["temperature_c"], 25)
        self.assertEqual(result["model"]["sodium_m"], .22)
        self.assertEqual(result["model"]["material"], self.backend.provenance.model)
        self.assertTrue(all(aptamer == "AC" and trigger.endswith("GT") and "U" not in trigger for aptamer, trigger in self.fake.mfe_sequences))
        self.assertEqual(result["model"]["aptamer_material"], "dna")

    def test_unsupported_conditions_fail_before_calculation(self):
        for changes in ({"trigger_material": "invalid"}, {"trigger_material": "rna"}, {"switch_material": "rna"}, {"magnesium_m": .01},
                        {"sodium_m": .119}, {"sodium_m": 1.001}, {"temperature_c": -.01}):
            with self.subTest(changes=changes), self.assertRaises(ValueError):
                self.run_search(thermo=replace(ThermoConditions(switch_material="dna"), **changes))
        self.import_context.assert_not_called()

    def test_real_calculation_errors_propagate_without_fallback(self):
        self.fake.fail_mfe = True
        with self.assertRaisesRegex(RuntimeError, "calculation failure"):
            self.run_search()
        with patch("aptaswitch_core.thermodynamics.setup.get_backend", side_effect=RuntimeError("missing native runtime")):
            with self.assertRaisesRegex(RuntimeError, "missing native"):
                self.run_search()

    def test_phase_progress_and_final_count(self):
        messages = []
        self.run_search(progress_callback=lambda *args: messages.append(args))
        ensembles = [message for message in messages if message[2].startswith("Ensemble ·")]
        self.assertEqual(ensembles[0][:2], (0, 20))
        self.assertEqual(messages[-1][:2], (20, 20))


@unittest.skipUnless(os.environ.get("APTASWITCH_HYBRID_NATIVE_LIBRARY"), "Native hybrid library required")
class RealNativeParityTests(unittest.TestCase):
    def test_dna_extensions_match_direct_native_calls(self):
        from aptaswitch_core.thermodynamics.setup import get_backend
        from aptaswitch_core.thermodynamics.types import ComplexSpec, Conditions, Strand
        backend = get_backend()
        aptamer, trigger = "ACGTACGT", "ACGT"
        result = extension.design_trigger_extensions(
            aptamer_dna=aptamer, trigger_dna=trigger, target_lengths=[5],
            thermo=ThermoConditions(switch_material="dna", trigger_material="dna", temperature_c=37))
        self.assertEqual(result["status"], "completed")
        self.assertEqual(result["mfe_screen"]["total_candidates_screened"], 4)
        self.assertEqual(result["model"]["material"], backend.provenance.model)
        self.assertEqual(result["model"]["aptamer_material"], "dna")
        self.assertTrue(result["candidates"])
        for row, sequence in [(result["reference"], trigger), *((r, r["extended_trigger"]) for r in result["candidates"])]:
            spec = ComplexSpec((Strand(aptamer, "dna", "aptamer"), Strand(sequence, "dna", "trigger")))
            expected = backend.minimum(spec, Conditions(37., .12, 0.))
            self.assertEqual(row["mfe_structure"], expected.structure)
            self.assertAlmostEqual(row["mfe_energy_kcal_mol"], expected.energy_kcal_mol, delta=1e-6)


if __name__ == "__main__":
    unittest.main()
