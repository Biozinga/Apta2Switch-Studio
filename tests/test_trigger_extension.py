"""Algorithm parity, genuine enumeration, cancellation and NUPACK integration."""

from __future__ import annotations

import json
import os
from pathlib import Path
import sys
from types import SimpleNamespace
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


class FakeNupack:
    """An explicitly synthetic NUPACK test double, only for orchestration tests."""

    def __init__(self):
        self.mfe_sequences = []
        self.ensemble_sequences = []
        self.model_options = None
        self.fail_mfe = False

    def Model(self, **options):
        self.model_options = options
        return options

    def Strand(self, sequence, name):
        return SimpleNamespace(sequence=sequence, name=name)

    class Complex:
        def __init__(self, strands, name):
            self.strands = strands
            self.name = name

    def mfe(self, *, strands, model):
        self.mfe_sequences.append(tuple(strands))
        if self.fail_mfe:
            raise RuntimeError("NUPACK calculation failure")
        return [SimpleNamespace(structure="+".join("." * (len(seq) - 1) for seq in strands), energy=-1.0)]

    def complex_analysis(self, *, complexes, model, compute):
        complex_obj = complexes[0]
        sequences = [strand.sequence for strand in complex_obj.strands]
        self.ensemble_sequences.append(tuple(sequences))
        structure = "+".join("." * (len(seq) - 1) for seq in sequences)
        matrix = np.eye(sum(len(seq) - 1 for seq in sequences))
        return {complex_obj: SimpleNamespace(
            mfe=[SimpleNamespace(structure=structure, energy=-1.0)], free_energy=-1.2,
            pairs=SimpleNamespace(to_array=lambda: matrix),
        )}

    def structure_probability(self, **kwargs):
        return 1.0


class SearchOrchestrationTests(unittest.TestCase):
    def setUp(self):
        from aptaswitch_core.localization import language_context
        self.enterContext(language_context("fr"))
        self.fake = FakeNupack()
        self.context = patch.object(extension, "nupack_import_context")
        self.import_context = self.context.start()
        self.import_context.return_value.__enter__.return_value = self.fake
        self.addCleanup(self.context.stop)

    def run_search(self, **kwargs):
        options = dict(aptamer_dna="AC", trigger_dna="GT", target_lengths=[3, 4],
                       thermo=ThermoConditions(material="dna04", temperature_c=37), import_path="/test/nupack", top=2)
        options.update(kwargs)
        return extension.design_trigger_extensions(**options)

    def test_enumerates_every_sequence_and_preserves_top_for_every_length(self):
        result = self.run_search()
        self.assertEqual(result["status"], "completed")
        self.assertEqual(len(self.fake.mfe_sequences), 20)
        self.assertEqual(len(set(self.fake.mfe_sequences)), 20)
        self.assertTrue(all(aptamer == "dAC" and trigger.endswith("GT") for aptamer, trigger in self.fake.mfe_sequences))
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

    def test_model_honors_exact_conditions_and_rna_chemistry(self):
        result = self.run_search(thermo=ThermoConditions(material="rna06", temperature_c=25, sodium_m=.22, magnesium_m=.01))
        self.assertEqual(self.fake.model_options, {"material": "rna06", "celsius": 25, "sodium": .22, "magnesium": .01, "ensemble": "stacking"})
        self.assertTrue(all(aptamer == "rAC" and trigger.endswith("GU") and "T" not in trigger for aptamer, trigger in self.fake.mfe_sequences))
        self.assertEqual(result["model"]["aptamer_material"], "rna")

    def test_real_calculation_errors_propagate_without_fallback(self):
        self.fake.fail_mfe = True
        with self.assertRaisesRegex(RuntimeError, "calculation failure"):
            self.run_search()
        with self.assertRaisesRegex(RuntimeError, "obligatoire"):
            self.run_search(import_path="")

    def test_phase_progress_and_final_count(self):
        messages = []
        self.run_search(progress_callback=lambda *args: messages.append(args))
        ensembles = [message for message in messages if message[2].startswith("Ensemble NUPACK")]
        self.assertEqual(ensembles[0][:2], (0, 20))
        self.assertEqual(messages[-1][:2], (20, 20))


@unittest.skipUnless(os.environ.get("APTASWITCH_TEST_NUPACK"), "Set APTASWITCH_TEST_NUPACK for real NUPACK parity")
class RealNupackParityTests(unittest.TestCase):
    def test_fixed_length_search_matches_original_nupack_numerically(self):
        original = original_reference()["real_nupack_case"]
        import_path = os.environ["APTASWITCH_TEST_NUPACK"]
        result = extension.design_trigger_extensions(
            aptamer_dna=original["aptamer"], trigger_dna=original["trigger"],
            target_lengths=[original["target_length"]],
            thermo=ThermoConditions(
                material=original["material"], temperature_c=original["temperature_c"],
                sodium_m=original["sodium_m"], magnesium_m=original["magnesium_m"],
            ),
            import_path=import_path,
        )
        expected = original["candidates"]
        self.assertEqual(result["mfe_screen"]["total_candidates_screened"], 4)
        self.assertEqual(len(expected), len(result["candidates"]))
        self.assertEqual(result["reference"]["mfe_structure"], original["reference_mfe_structure"])
        for actual, source in zip(result["candidates"], expected):
            self.assertEqual(actual["extended_trigger"], source["extended_z1"])
            self.assertEqual(actual["mfe_structure"], source["mfe_structure"])
            for key in ("final_score", "mfe_energy_kcal_mol", "core_pair_probability_rmse",
                        "reference_pair_probability_loss", "target_structure_probability", "extension_pair_probability_sum"):
                # The reference was calculated with NUPACK 4.1.0.1; allow only
                # small numerical differences between native platform builds.
                self.assertAlmostEqual(actual[key], source[key], delta=1e-6, msg=key)


if __name__ == "__main__":
    unittest.main()
