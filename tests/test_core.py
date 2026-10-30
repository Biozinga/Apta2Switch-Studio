from __future__ import annotations

import tempfile
import unittest
from types import SimpleNamespace
from pathlib import Path
from unittest.mock import patch
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from aptaswitch_core.architecture import default_architecture_for_trigger
from aptaswitch_core.exports import export_run
from aptaswitch_core.engines import run_design
from tests.design_fixtures import design_result_fixture
from aptaswitch_core.models import REPORTERS, RunRequest, ScoringWeights, SequenceInput, ThermoConditions
import aptaswitch_core.settings as settings
import aptaswitch_core.reporter_library as reporter_library
from aptaswitch_core.scoring import compute_selection_score, selection_score_100
from aptaswitch_core.sequences import (
    SequenceValidationError,
    normalize_dna,
)


class SequenceTests(unittest.TestCase):
    def test_normalize_dna_accepts_rna_input(self):
        self.assertEqual(normalize_dna("acgu acgu"), "ACGTACGT")

    def test_normalize_dna_rejects_iupac_ambiguity(self):
        with self.assertRaises(SequenceValidationError):
            normalize_dna("ACGR")



class ArchitectureTests(unittest.TestCase):
    def test_default_architecture_covers_trigger(self):
        architecture = default_architecture_for_trigger(24)
        self.assertEqual(architecture.toehold_length + architecture.stem_length, 24)
        architecture.validate_for_trigger(24)


class ThermoConditionsTests(unittest.TestCase):
    def test_unsupported_models_and_salts_are_rejected_before_calculation(self):
        invalid = [{"material": name} for name in ("dna04", "rna", "rna06", "rna95", "rna99")]
        invalid += [{"magnesium_m": .001}, {"sodium_m": .119}, {"sodium_m": 1.01},
                    {"temperature_c": float("nan")}, {"temperature_c": 24.},
                    {"temperature_c": 56.}, {"switch_material": "unknown"}]
        for values in invalid:
            with self.subTest(values=values), self.assertRaises(ValueError):
                ThermoConditions(**values).validate_model()

    def test_native_model_preserves_explicit_chemistry_at_supported_limits(self):
        from aptaswitch_core.thermodynamics.setup import MODEL_ID
        for switch, trigger in (("rna", "dna"), ("dna", "rna"), ("dna", "dna"), ("rna", "rna")):
            for temperature, salt in ((25., .12), (55., 1.)):
                thermo = ThermoConditions(switch_material=switch, trigger_material=trigger,
                                         temperature_c=temperature, sodium_m=salt)
                thermo.validate_model()
                self.assertEqual(thermo.material, MODEL_ID)
                self.assertEqual((thermo.switch_material, thermo.trigger_material), (switch, trigger))

    def test_separate_species_concentrations_override_legacy_value(self):
        thermo = ThermoConditions(
            concentration_m=5e-6,
            trigger_concentration_m=12e-9,
            switch_concentration_m=250e-9,
        )
        self.assertEqual(thermo.resolved_trigger_concentration_m, 12e-9)
        self.assertEqual(thermo.resolved_switch_concentration_m, 250e-9)

    def test_legacy_shared_concentration_remains_supported(self):
        thermo = ThermoConditions(concentration_m=8e-6)
        self.assertEqual(thermo.resolved_trigger_concentration_m, 8e-6)
        self.assertEqual(thermo.resolved_switch_concentration_m, 8e-6)


class ConcentrationTests(unittest.TestCase):
    def test_complex_percentage_uses_each_species_concentration(self):
        from aptaswitch_core.saved_switch_analysis import _population
        from aptaswitch_core.thermodynamics.types import ComplexSpec, Strand
        trigger, switch = Strand("AAAA", "dna", "Trigger"), Strand("UUUU", "rna", "Switch")
        populations = {ComplexSpec((trigger, switch)): 2e-9, ComplexSpec((switch,)): 2e-9}
        result = _population(populations, {"Trigger": 2e-9, "Switch": 4e-9}, 2)
        self.assertEqual(result["target_switch_fraction_pct"], 50.)
        self.assertEqual(result["target_trigger_fraction_pct"], 100.)


class ScoringTests(unittest.TestCase):
    def test_default_weights_prioritize_rbs_and_defect_with_correct_directions(self):
        metrics = dict(on_yield_pct=90, leak_pct=2, functional_defect=0.02,
                       ddg_activation=-12, delta_g_rbs_linker=-2,
                       weights=ScoringWeights())
        baseline = compute_selection_score(**metrics)
        gains = []
        for name, value, expected_gain in (
            ("delta_g_rbs_linker", -1, 6),
            ("functional_defect", 0.01, 6),
            ("on_yield_pct", 91, 3),
            ("ddg_activation", -13, 0.25),
            ("leak_pct", 1, 0.0),
        ):
            with self.subTest(metric=name):
                gain = compute_selection_score(**{**metrics, name: value}) - baseline
                self.assertAlmostEqual(gain, expected_gain)
                gains.append(gain)
        self.assertEqual(gains, sorted(gains, reverse=True))

    def test_score_keeps_ranking_information_above_100_and_below_zero(self):
        metrics = dict(on_yield_pct=90, leak_pct=2, functional_defect=0.02,
                       ddg_activation=-12, delta_g_rbs_linker=-2,
                       weights=ScoringWeights())
        baseline = compute_selection_score(**metrics)
        higher_on = compute_selection_score(**{**metrics, "on_yield_pct": 91})
        self.assertGreater(baseline, 100)
        self.assertAlmostEqual(higher_on - baseline, 3)

        poor = {**metrics, "on_yield_pct": 5, "functional_defect": 0.8}
        baseline = compute_selection_score(**poor)
        higher_on = compute_selection_score(**{**poor, "on_yield_pct": 6})
        self.assertLess(baseline, 0)
        self.assertAlmostEqual(higher_on - baseline, 3)

    def test_association_with_aptamer_is_informative_by_default_but_can_be_weighted(self):
        from dataclasses import replace
        metrics = dict(on_yield_pct=90, leak_pct=0, functional_defect=0.02,
                       ddg_activation=-12, delta_g_rbs_linker=-2,
                       weights=ScoringWeights())
        baseline = compute_selection_score(**metrics)
        self.assertEqual(compute_selection_score(**{**metrics, "leak_pct": 100}), baseline)
        weighted = {**metrics, "weights": replace(metrics["weights"], leak=.2)}
        self.assertAlmostEqual(compute_selection_score(**weighted) -
                               compute_selection_score(**{**weighted, "leak_pct": 10}), 2.)


class NativeValidationTests(unittest.TestCase):
    def test_missing_path_fails_cleanly(self):
        from aptaswitch_core.thermodynamics.setup import validate_runtime
        result = validate_runtime("/definitely/missing/native/library")
        self.assertFalse(result.ok)
        self.assertTrue(result.message)


class ReporterLibraryTests(unittest.TestCase):
    def setUp(self):
        self._tmpdir = tempfile.TemporaryDirectory()
        tmp_dir = Path(self._tmpdir.name)
        patcher_dir = patch.object(settings, "SETTINGS_DIR", tmp_dir)
        patcher_path = patch.object(settings, "SETTINGS_PATH", tmp_dir / "settings.json")
        patcher_dir.start()
        patcher_path.start()
        self.addCleanup(patcher_dir.stop)
        self.addCleanup(patcher_path.stop)
        self.addCleanup(self._tmpdir.cleanup)

    def test_all_reporters_includes_builtins_only_by_default(self):
        reporters = reporter_library.all_reporters()
        self.assertEqual(set(reporters), set(REPORTERS))

    def test_add_custom_reporter_appears_in_all_reporters(self):
        reporter = reporter_library.add_custom_reporter("mCherry", "atg" + "acg" * 15 + "taa")
        self.assertEqual(reporter.sequence_after_start_rna, "ACG" * 10)
        self.assertEqual(reporter.full_sequence_after_start_rna, "ACG" * 15 + "UAA")
        reporters = reporter_library.all_reporters()
        self.assertIn(reporter.key, reporters)
        self.assertFalse(reporter_library.is_builtin_reporter(reporter.key))

    def test_add_custom_reporter_rejects_invalid_sequence(self):
        with self.assertRaises(SequenceValidationError):
            reporter_library.add_custom_reporter("Bad", "ACGX")

    def test_add_custom_reporter_rejects_empty_name(self):
        with self.assertRaises(ValueError):
            reporter_library.add_custom_reporter("  ", "ACGU")

    def test_duplicate_labels_get_unique_keys(self):
        first = reporter_library.add_custom_reporter("Dup", "ATG" + "ACG" * 15 + "TAA")
        second = reporter_library.add_custom_reporter("Dup", "ATG" + "ACG" * 16 + "TAA")
        self.assertNotEqual(first.key, second.key)

    def test_remove_custom_reporter(self):
        reporter = reporter_library.add_custom_reporter("Temp", "ATG" + "ACG" * 15 + "TAA")
        reporter_library.remove_custom_reporter(reporter.key)
        self.assertNotIn(reporter.key, reporter_library.all_reporters())

    def test_cannot_remove_builtin_reporter(self):
        with self.assertRaises(ValueError):
            reporter_library.remove_custom_reporter("sfgfp")


class EngineDispatchTests(unittest.TestCase):
    def setUp(self):
        from aptaswitch_core.localization import language_context
        self.enterContext(language_context("fr"))

    def test_default_native_engine_passes_cancellation_to_design_search(self):
        from dataclasses import asdict
        from aptaswitch_core.thermodynamics.types import Conditions, ModelProvenance
        from aptaswitch_core.thermodynamics.setup import ENGINE_ID
        request = RunRequest(SequenceInput("ACGTACGT", "ACGT" * 6), default_architecture_for_trigger(24))
        self.assertEqual(request.engine, ENGINE_ID)
        backend = SimpleNamespace(provenance=ModelProvenance("fixture", "1", "fixture", "fixture:1", "connected", "1M"))
        batch = SimpleNamespace(candidates=[], status="cancelled", started_at="2026-01-01T00:00:00+00:00",
                                completed_at="2026-01-01T00:00:01+00:00", optimization_trials=[],
                                provenance=backend.provenance, conditions=Conditions(37., .12, 0.), analysis_cache={})
        cancelled = lambda: True
        progress = []
        with patch("aptaswitch_core.thermodynamics.setup.backend_for_request", return_value=backend), \
             patch("aptaswitch_core.engines.design_candidates", return_value=batch) as run:
            result = run_design(request, progress_callback=lambda *event: progress.append(event), cancel_check=cancelled)
        self.assertEqual(result.status, "cancelled")
        self.assertEqual(run.call_args.args, (request.normalized(), backend))
        self.assertIs(run.call_args.kwargs["cancel_check"], cancelled)
        self.assertTrue(progress)
        self.assertEqual(result.thermodynamics["parameter_digest"], backend.provenance.parameter_digest)

    def test_removed_and_unknown_engines_fail_before_any_calculation(self):
        from dataclasses import replace
        request = RunRequest(SequenceInput("ACGTACGT", "ACGT" * 6), default_architecture_for_trigger(24))
        for engine in ("mock", "", "other", "external"):
            with self.subTest(engine=engine), patch("aptaswitch_core.engines.design_candidates") as run:
                with self.assertRaises(ValueError):
                    run_design(replace(request, engine=engine))
                run.assert_not_called()


class OutputTests(unittest.TestCase):
    def test_run_exports_expected_files(self):
        with tempfile.TemporaryDirectory() as tmp:
            request = RunRequest(
                sequences=SequenceInput(
                    aptamer_dna="CTACCAGCTTTGAGGCTCGATCCAGCTTATTCAATTATACCAGCTTATTCAATTATACCAGC",
                    trigger_dna="GTCCAGGCTGGTATAATTAGATCC",
                ),
                architecture=default_architecture_for_trigger(24),
                trials=3,
                output_dir=Path(tmp),
            )
            result = design_result_fixture(request)
            result.candidates[0].score_raw = 250.75
            result.candidates[0].score = selection_score_100(250.75)
            paths = export_run(result)
            self.assertEqual(len(result.candidates), 3)
            self.assertTrue(Path(paths["json"]).exists())
            self.assertTrue(Path(paths["csv"]).exists())
            self.assertTrue(Path(paths["svg"]).exists())
            self.assertTrue((Path(paths["project"]).name == ".aptaswitch.json"))
            csv_text = Path(paths["csv"]).read_text()
            self.assertIn("Score_selection,", csv_text)
            self.assertNotIn("Score_selection_sur_100", csv_text)
            self.assertIn("250.75", csv_text)
            svg = Path(paths["svg"]).read_text()
            self.assertIn(f"{selection_score_100(250.75):.4f}", svg)
            self.assertIn("/100", svg)


if __name__ == "__main__":
    unittest.main()
