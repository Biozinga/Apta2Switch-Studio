"""Tube displays keep strand fractions distinct from complex particle counts."""

from copy import deepcopy
import unittest
from unittest.mock import patch

from tests.ui_helpers import FrenchAppTest as AppTest

from aptaswitch_studio.tube_analysis_ui import (
    _composition_rows, _format_number, _reference_delta, _reference_is_compatible,
)


def _tube(mixed=False):
    summary = {
        "trigger_free_pct": 30.0 if mixed else 60.0,
        "trigger_target_pct": 50.0 if mixed else 0.0,
        "trigger_parasite_pct": 20.0 if mixed else 40.0,
        "aptamer_free_pct": 50.0 if mixed else 0.0,
        "aptamer_target_pct": 50.0 if mixed else 0.0,
        "aptamer_parasite_pct": 0.0,
        "parasite_particle_fraction_pct": 8.33 if mixed else 25.0,
        "parasite_complex_molecule_count": 1e10 if mixed else 2e10,
        "total_complex_molecule_count": 12e10 if mixed else 8e10,
    }
    # These small fixtures preserve the distinction between strand and particle fractions.
    populations = [("T", "free", 0, 1, 30 if mixed else 60, 0, 25 if mixed else 75),
                   ("2 T", "parasite", 0, 2, 20 if mixed else 40, 0, 8.33 if mixed else 25)]
    if mixed:
        populations += [("A + T", "target", 1, 1, 50, 50, 41.67),
                        ("A", "free", 1, 0, 0, 50, 25)]
    complexes = [{
        "label": label, "category": category, "aptamer_count": n_a,
        "trigger_count": n_t, "strand_count": n_a + n_t,
        "concentration_m": particle * 1e-8,
        "complex_molecule_count": particle * 1e8,
        "trigger_fraction_pct": trigger, "aptamer_fraction_pct": aptamer,
        "particle_fraction_pct": particle,
    } for label, category, n_a, n_t, trigger, aptamer, particle in populations]
    return {"summary": summary, "complexes": complexes}


def analysis_fixture():
    return {
        "status": "completed",
        "model": {"material": "experimental-Sugimoto95-Vienna04-loopmean", "aptamer_material": "dna", "trigger_material": "dna",
                  "temperature_c": 37., "sodium_m": .12, "magnesium_m": 0., "ensemble": "stacking"},
        "conditions": {"aptamer_concentration_m": 5e-6, "trigger_concentration_m": 5e-6,
                       "volume_l": 20e-6, "max_complex_size": 5},
        "tubes": {"trigger_only": _tube(), "aptamer_trigger": _tube(mixed=True)},
    }


def tube_app(analysis, reference=None):
    from aptaswitch_studio.tube_analysis_ui import render_tube_analysis
    render_tube_analysis(analysis, key="tube_test", reference=reference)


def english_tube_app(analysis):
    from aptaswitch_core.localization import set_language
    from aptaswitch_studio.tube_analysis_ui import render_tube_analysis
    set_language("en")
    render_tube_analysis(analysis, key="tube_test_english")


class TubeAnalysisUiTests(unittest.TestCase):
    def test_composition_charts_use_input_strands_not_complex_particles(self):
        rows = _composition_rows(_tube(), mixed=False)
        self.assertEqual([row["fraction_pct"] for row in rows], [60.0, 40.0])
        mixed = _composition_rows(_tube(mixed=True), mixed=True)
        self.assertEqual([row["fraction_pct"] for row in mixed], [30, 50, 20, 50, 50, 0])

    def test_completed_view_shows_both_tubes_and_expected_counts_without_new_computation(self):
        with patch("aptaswitch_core.tube_analysis.analyze_parasite_tubes") as compute:
            app = AppTest.from_function(tube_app, args=(analysis_fixture(),)).run()
        self.assertFalse(app.exception)
        compute.assert_not_called()
        self.assertEqual([tab.label for tab in app.tabs], ["Trigger seul", "Aptamère + trigger"])
        self.assertEqual(app.metric[0].value, "40 %")
        self.assertEqual(app.metric[3].value, "20 %")
        self.assertEqual(app.metric[4].value, "50 %")
        self.assertEqual(len(app.get("vega_lite_chart")), 2)
        self.assertEqual(len(app.table), 2)
        self.assertTrue(any("20 µL" in caption.value and "5 brins" in caption.value for caption in app.caption))
        self.assertTrue(any("Aptamère (tube mixte)" in caption.value for caption in app.caption))
        self.assertTrue(any("sans ligand" in metric.help for metric in app.metric))
        self.assertTrue(any("ni un appariement correct" in metric.help for metric in app.metric if metric.help))

    def test_reference_comparison_reports_percentage_points(self):
        reference = analysis_fixture()
        candidate = deepcopy(reference)
        candidate["tubes"]["trigger_only"]["summary"]["trigger_parasite_pct"] = 35.0
        app = AppTest.from_function(tube_app, args=(candidate, reference)).run()
        self.assertFalse(app.exception)
        self.assertIn("-5.000", app.metric[0].delta)
        self.assertIn("original", app.metric[0].delta)
        self.assertIsNone(_reference_delta({"status": "failed"}, "trigger_only", 10))

    def test_missing_failed_and_cancelled_results_never_show_zero_metrics(self):
        for analysis in ({}, {"status": "failed", "error": "Native calculation failed"}, {"status": "cancelled", "tubes": {}}):
            with self.subTest(analysis=analysis):
                app = AppTest.from_function(tube_app, args=(analysis,)).run()
                self.assertFalse(app.exception)
                self.assertFalse(app.metric)
                self.assertFalse(app.tabs)
                self.assertTrue(app.caption or app.warning)

    def test_small_expected_populations_remain_visible(self):
        self.assertEqual(_format_number(0), "0")
        self.assertEqual(_format_number(1e-7), "1.00 × 10⁻⁷")
        self.assertEqual(_format_number(6.02214076e13), "6.02 × 10¹³")
        self.assertEqual(_format_number(float("nan")), "—")
        self.assertEqual(_format_number(500000), "500 000")
        self.assertEqual(_format_number(15362), "15 400")

    def test_reference_deltas_require_matching_equilibrium_conditions(self):
        candidate = analysis_fixture()
        for field, value in (("trigger_concentration_m", 1e-6),
                             ("aptamer_concentration_m", 2e-6), ("max_complex_size", 2)):
            with self.subTest(field=field):
                reference = deepcopy(candidate)
                reference["conditions"][field] = value
                self.assertFalse(_reference_is_compatible(candidate, reference))
        reference = deepcopy(candidate)
        reference["model"] = {**reference["model"], "trigger_material": "rna", "aptamer_material": "rna"}
        self.assertFalse(_reference_is_compatible(candidate, reference))
        app = AppTest.from_function(tube_app, args=(candidate, reference)).run()
        self.assertFalse(app.exception)
        self.assertFalse(app.metric[0].delta)
        self.assertTrue(any("conditions de calcul diffèrent" in caption.value for caption in app.caption))
        reference = deepcopy(candidate)
        reference["conditions"]["volume_l"] = 1e-6
        self.assertTrue(_reference_is_compatible(candidate, reference))

    def test_english_view_translates_both_tubes_and_the_scientific_explanation(self):
        app = AppTest.from_function(english_tube_app, args=(analysis_fixture(),)).run()
        self.assertFalse(app.exception)
        self.assertEqual([tab.label for tab in app.tabs], ["Trigger only", "Aptamer + trigger"])
        self.assertEqual(app.metric[0].label, "Trigger in parasitic complexes")
        self.assertEqual(app.metric[0].value, "40 %")
        self.assertEqual(app.expander[0].label, "Detailed tube composition")
        self.assertTrue(any("Equilibrium prediction without ligand" in metric.help for metric in app.metric))
        self.assertTrue(any("does not establish ligand binding or correct base pairing" in metric.help
                            for metric in app.metric if metric.help))


if __name__ == "__main__":
    unittest.main()
