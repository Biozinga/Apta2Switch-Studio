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
from aptaswitch_core.nupack_engine import _complex_pct
import aptaswitch_core.nupack_setup as nupack_setup
from aptaswitch_core.nupack_setup import validate_nupack_path
import aptaswitch_core.reporter_library as reporter_library
from aptaswitch_core.scoring import compute_selection_score
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


class NupackConcentrationTests(unittest.TestCase):
    def test_complex_percentage_uses_each_species_concentration(self):
        class FakeStrand:
            def __init__(self, sequence, name):
                self.sequence = sequence
                self.name = name

        class FakeTube:
            def __init__(self, *, strands, complexes, name):
                self.strands = strands
                self.complexes = complexes
                self.name = name

        class FakeNupack:
            Strand = FakeStrand
            Tube = FakeTube
            SetSpec = lambda self, **kwargs: kwargs

            def tube_analysis(self, *, tubes, model):
                tube = tubes[0]
                self.last_tube = tube
                class ExpectedComplex:
                    strands = [
                        SimpleNamespace(name="Trigger"),
                        SimpleNamespace(name="Switch"),
                    ]

                expected = ExpectedComplex()
                analysis = SimpleNamespace(complex_concentrations={expected: 2e-9})
                return SimpleNamespace(tubes={tube: analysis})

        nupack = FakeNupack()
        percentage = _complex_pct(
            nupack,
            strands=[("Trigger", "AAAA", 2e-9), ("Switch", "UUUU", 4e-9)],
            model=object(),
            denominator_m=4e-9,
            max_size=2,
            expected_names=["Trigger", "Switch"],
        )

        concentrations = {strand.name: value for strand, value in nupack.last_tube.strands.items()}
        self.assertEqual(concentrations, {"Trigger": 2e-9, "Switch": 4e-9})
        self.assertEqual(percentage, 50.0)


class ScoringTests(unittest.TestCase):
    def test_default_weights_prioritize_rbs_and_defect_with_correct_directions(self):
        metrics = dict(on_yield_pct=90, leak_pct=2, defect=0.02,
                       ddg_activation=-12, delta_g_rbs_linker=-2,
                       weights=ScoringWeights())
        baseline = compute_selection_score(**metrics)
        gains = []
        for name, value, expected_gain in (
            ("delta_g_rbs_linker", -1, 6),
            ("defect", 0.01, 6),
            ("on_yield_pct", 91, 3),
            ("ddg_activation", -13, 0.25),
            ("leak_pct", 1, 0.05),
        ):
            with self.subTest(metric=name):
                gain = compute_selection_score(**{**metrics, name: value}) - baseline
                self.assertAlmostEqual(gain, expected_gain)
                gains.append(gain)
        self.assertEqual(gains, sorted(gains, reverse=True))

    def test_score_keeps_ranking_information_above_100_and_below_zero(self):
        metrics = dict(on_yield_pct=90, leak_pct=2, defect=0.02,
                       ddg_activation=-12, delta_g_rbs_linker=-2,
                       weights=ScoringWeights())
        baseline = compute_selection_score(**metrics)
        higher_on = compute_selection_score(**{**metrics, "on_yield_pct": 91})
        self.assertGreater(baseline, 100)
        self.assertAlmostEqual(higher_on - baseline, 3)

        poor = {**metrics, "on_yield_pct": 5, "defect": 0.8}
        baseline = compute_selection_score(**poor)
        higher_on = compute_selection_score(**{**poor, "on_yield_pct": 6})
        self.assertLess(baseline, 0)
        self.assertAlmostEqual(higher_on - baseline, 3)

    def test_leak_only_reverses_close_candidates_with_default_weights(self):
        metrics = dict(on_yield_pct=90, leak_pct=0, defect=0.02,
                       ddg_activation=-12, delta_g_rbs_linker=-2,
                       weights=ScoringWeights())
        baseline = compute_selection_score(**metrics)
        # Even the entire 0–100% leak range cannot outweigh two ON points.
        clear_on_advantage = compute_selection_score(
            **{**metrics, "on_yield_pct": 92, "leak_pct": 100},
        )
        self.assertGreater(clear_on_advantage, baseline)
        # A small leak improvement still differentiates otherwise near-equals.
        near_equal = compute_selection_score(
            **{**metrics, "on_yield_pct": 90.01, "leak_pct": 1},
        )
        self.assertGreater(baseline, near_equal)


class NupackValidationTests(unittest.TestCase):
    def setUp(self):
        from aptaswitch_core.localization import language_context
        self.enterContext(language_context("en"))

    def test_missing_path_fails_cleanly(self):
        result = validate_nupack_path("/definitely/missing/nupack/path")
        self.assertFalse(result.ok)
        self.assertIn("does not exist", result.message)


class ReporterLibraryTests(unittest.TestCase):
    def setUp(self):
        self._tmpdir = tempfile.TemporaryDirectory()
        tmp_dir = Path(self._tmpdir.name)
        patcher_dir = patch.object(nupack_setup, "SETTINGS_DIR", tmp_dir)
        patcher_path = patch.object(nupack_setup, "SETTINGS_PATH", tmp_dir / "settings.json")
        patcher_dir.start()
        patcher_path.start()
        self.addCleanup(patcher_dir.stop)
        self.addCleanup(patcher_path.stop)
        self.addCleanup(self._tmpdir.cleanup)

    def test_all_reporters_includes_builtins_only_by_default(self):
        reporters = reporter_library.all_reporters()
        self.assertEqual(set(reporters), set(REPORTERS))

    def test_add_custom_reporter_appears_in_all_reporters(self):
        reporter = reporter_library.add_custom_reporter("mCherry", "acgu" * 5)
        self.assertEqual(reporter.sequence_after_start_rna, "ACGU" * 5)
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
        first = reporter_library.add_custom_reporter("Dup", "ACGU")
        second = reporter_library.add_custom_reporter("Dup", "ACGUACGU")
        self.assertNotEqual(first.key, second.key)

    def test_remove_custom_reporter(self):
        reporter = reporter_library.add_custom_reporter("Temp", "ACGU")
        reporter_library.remove_custom_reporter(reporter.key)
        self.assertNotIn(reporter.key, reporter_library.all_reporters())

    def test_cannot_remove_builtin_reporter(self):
        with self.assertRaises(ValueError):
            reporter_library.remove_custom_reporter("sfgfp")


class EngineDispatchTests(unittest.TestCase):
    def setUp(self):
        from aptaswitch_core.localization import language_context
        self.enterContext(language_context("fr"))

    def test_default_engine_is_nupack_and_callbacks_reach_it(self):
        request = RunRequest(
            sequences=SequenceInput(aptamer_dna="ACGTACGT", trigger_dna="ACGT" * 6),
            architecture=default_architecture_for_trigger(24),
        )
        self.assertEqual(request.engine, "nupack")
        progress, cancel = object(), object()
        sentinel = object()
        with patch("aptaswitch_core.engines.run_nupack_engine", return_value=sentinel) as run:
            self.assertIs(run_design(request, progress_callback=progress, cancel_check=cancel), sentinel)
        run.assert_called_once_with(request, progress_callback=progress, cancel_check=cancel)

    def test_removed_and_unknown_engines_fail_before_any_calculation(self):
        for engine in ("mock", "", "other"):
            with self.subTest(engine=engine), patch("aptaswitch_core.engines.run_nupack_engine") as run:
                with self.assertRaisesRegex(ValueError, "Seul le moteur NUPACK"):
                    run_design(SimpleNamespace(engine=engine))
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
            result.candidates[0].score = 250.75
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
            self.assertIn("250.7500", svg)
            self.assertNotIn("/100", svg)


if __name__ == "__main__":
    unittest.main()
