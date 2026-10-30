"""Saved inverse-extension controls and comparison; no calculations or network calls."""

from __future__ import annotations

import copy
import csv
import io
from dataclasses import replace
from pathlib import Path
from types import SimpleNamespace
import unittest
from unittest.mock import patch

from tests.ui_helpers import FrenchAppTest as AppTest

from aptaswitch_studio import inverse_extension_ui as ui
from aptaswitch_core.scoring import selection_score_100


def controls_app(trigger="GCTT"):
    from aptaswitch_studio.inverse_extension_ui import render_inverse_extension_controls
    render_inverse_extension_controls(trigger)


def make_result():
    candidates = []
    for trial, score, impact, extended in (
        (1, 30.0, 9.0, "CCGCTTAA"),
        (2, 20.0, 4.0, "AAGCTTCC"),
        (3, 10.0, 1.0, "GGGCTTTT"),
    ):
        metrics = {
            "status": "completed", "trial": trial, "original_switch_rank": trial,
            "extended_trigger": extended, "extension_5prime": extended[:2],
            "extension_3prime": extended[-2:], "extension_length": 4, "trigger_length": 8,
            "mfe_structure": "((((+..))))..", "mfe_energy_kcal_mol": -3.0 - trial,
            "nucleotide_probabilities": [0.7] * 12, "final_score": impact,
        }
        candidates.append(SimpleNamespace(trial=trial, score=selection_score_100(score), score_raw=score, defect=0.1,
                                          regional_defects={"scheme": "functional_v2", "balanced": 0.1},
                                          trigger_linearity={}, excluded=False, inverse_extension_metrics=metrics))
    return SimpleNamespace(
        run_id="inverse_ui_fixture", candidates=candidates,
        request=SimpleNamespace(inverse_extension=SimpleNamespace(top_candidates=3),
                                score_scale="bounded_v1", defect_metric="functional",
                                scoring=SimpleNamespace(defect=6.0)),
        inverse_extension_analysis={
            "status": "completed", "selected_count": 3, "analyzed_count": 3,
            "model": {"material": "experimental-Sugimoto95-Vienna04-loopmean", "aptamer_material": "dna", "trigger_material": "dna",
                      "temperature_c": 37, "sodium_m": 0.12, "magnesium_m": 0.0},
            "reference": {
                "aptamer_sequence": "AAGC", "trigger_sequence": "GCTT",
                "mfe_structure": "((((+))))", "mfe_energy_kcal_mol": -4.0,
                "nucleotide_probabilities": [0.8] * 8,
            },
            "candidates": [copy.deepcopy(candidate.inverse_extension_metrics) for candidate in candidates],
        },
    )


def results_app():
    import streamlit as st
    from tests.test_inverse_extension_ui import make_result
    from aptaswitch_studio.inverse_extension_ui import render_inverse_extension_results
    if "saved_result" not in st.session_state:
        st.session_state.saved_result = make_result()
    render_inverse_extension_results(st.session_state.saved_result)


def reanalysis_app():
    import streamlit as st
    from tests.test_inverse_extension_ui import make_result
    from aptaswitch_studio.inverse_extension_ui import render_inverse_extension_results
    if "saved_result" not in st.session_state:
        st.session_state.saved_result = make_result()

    def record_request(result, requested):
        st.session_state.reanalysis_requested = (result.run_id, requested)

    render_inverse_extension_results(st.session_state.saved_result, on_reanalyze=record_request)


def saved_comparison_app():
    import streamlit as st
    from aptaswitch_studio.inverse_extension_ui import render_inverse_extension_candidate
    result = st.session_state.saved_result
    render_inverse_extension_candidate(result, result.candidates[0], key_prefix="retained")


def historical_reanalysis_app():
    import streamlit as st
    from aptaswitch_studio.inverse_extension_ui import render_inverse_extension_results
    from aptaswitch_studio.web_app import _start_inverse_reanalysis

    render_inverse_extension_results(st.session_state.saved_result, on_reanalyze=_start_inverse_reanalysis)
    if st.session_state.get("run_error"):
        st.error(st.session_state.run_error)


def historical_result(magnesium=0., sodium=.15):
    from aptaswitch_core.architecture import default_architecture_for_trigger
    from aptaswitch_core.models import InverseExtensionSettings, RunRequest, RunResult, SequenceInput, ThermoConditions

    request = RunRequest(
        sequences=SequenceInput("ACGTACGT", "ACGT" * 6),
        architecture=default_architecture_for_trigger(30),
        inverse_extension=InverseExtensionSettings(
            target_length=30, top_candidates=3,
            thermo=ThermoConditions(material="dna04", switch_material="dna", trigger_material="dna",
                                    temperature_c=29.5, sodium_m=sodium, magnesium_m=magnesium),
        ),
        output_dir=Path("unused"),
    )
    return RunResult("saved-inverse", request, [], Path("unused/saved-inverse"), "start", "end")


class RunningJob:
    def snapshot(self):
        return SimpleNamespace(state="running")


class InverseExtensionControlTests(unittest.TestCase):
    def setUp(self):
        from tests.native_ui_fixtures import patch_native_runtime
        patch_native_runtime(self)

    def test_selected_mode_controls_leave_existing_inputs_untouched(self):
        app = AppTest.from_function(controls_app)
        app.session_state.trigger_input = "GCTT"
        app.session_state.aptamer_input = "AAGC"
        app.run()
        self.assertFalse(app.exception)
        self.assertFalse(app.checkbox)
        self.assertEqual(app.number_input(key="inverse_extension_target_length").value, 30)
        self.assertEqual(app.number_input(key="inverse_extension_top").value, 20)
        self.assertEqual(app.code[0].value, "X" * 26 + "GCTT")
        self.assertEqual(app.session_state.trigger_input, "GCTT")
        self.assertEqual(app.session_state.aptamer_input, "AAGC")
        self.assertTrue(any("conditions du système initial" in item.proto.help for item in app.subheader))

    def test_two_ends_assign_odd_extra_base_to_five_prime_and_preserve_core(self):
        app = AppTest.from_function(controls_app, args=("gc uu",))
        app.session_state.inverse_extension_target_length = 9
        app.session_state.inverse_extension_side = "both"
        app.run()
        self.assertFalse(app.exception)
        self.assertEqual(app.code[0].value, "XXXGCTTXX")
        app.radio(key="inverse_extension_side").set_value("3prime").run()
        self.assertEqual(app.code[0].value, "GCTTXXXXX")

    def test_incomplete_or_too_long_original_has_no_invented_preview(self):
        for trigger in ("", "ACNX", "GCTT"):
            with self.subTest(trigger=trigger):
                app = AppTest.from_function(controls_app, args=(trigger,))
                app.session_state.inverse_extension_target_length = 4
                app.run()
                self.assertFalse(app.exception)
                self.assertFalse(app.code)
                self.assertEqual(app.session_state.inverse_extension_target_length, 4)
                if trigger == "GCTT":
                    self.assertIn("dépasser les 4 nt", app.error[0].value)
                else:
                    self.assertTrue(app.info)


class InverseExtensionResultsTests(unittest.TestCase):
    def test_completed_stopped_run_shows_saved_comparison_and_optional_controls(self):
        result = make_result()
        result.status = "stopped"
        app = AppTest.from_function(reanalysis_app)
        app.session_state.saved_result = result
        app.run()
        self.assertFalse(app.exception)
        self.assertTrue(any("Comparaison enregistrée : 3/3" in item.value for item in app.caption))
        optional = next(item for item in app.expander if item.label == "Nouvelle comparaison (facultatif)")
        self.assertFalse(optional.proto.expanded)
        self.assertEqual(optional.button[0].label, "Lancer une nouvelle comparaison")
        self.assertNotIn("reanalysis_requested", app.session_state)

    def test_retained_comparison_uses_original_rank_and_never_recalculates(self):
        result = make_result()
        result.candidates = [result.candidates[1]]
        result.candidates[0].inverse_extension_metrics["original_switch_rank"] = 14
        before = copy.deepcopy(result)
        app = AppTest.from_function(saved_comparison_app)
        app.session_state.saved_result = result
        with patch.object(ui, "ranked_inverse_candidates", side_effect=AssertionError("No single-candidate ranking")), \
                patch("aptaswitch_core.inverse_extension.analyze_inverse_extensions", side_effect=AssertionError("No recalculation")):
            app.run()
        self.assertFalse(app.exception)
        values = {metric.label: metric.value for metric in app.metric}
        self.assertEqual(values["Rang switch initial"], "14")
        self.assertEqual(values["Impact structurel"], "4")
        self.assertEqual(len(app.get("html")), 2)
        self.assertFalse(app.slider)
        self.assertFalse(app.table)
        self.assertFalse(app.number_input)
        self.assertFalse(app.button)
        self.assertEqual(app.session_state.saved_result, before)

    def test_uncompared_retained_candidate_has_no_invented_metrics_or_structures(self):
        result = make_result()
        result.candidates[0].inverse_extension_metrics = {}
        app = AppTest.from_function(saved_comparison_app)
        app.session_state.saved_result = result
        app.run()
        self.assertFalse(app.exception)
        self.assertIn("Aucune comparaison", app.info[0].value)
        self.assertIn("3 meilleurs candidats", app.caption[0].value)
        self.assertFalse(app.metric)
        self.assertFalse(app.get("html"))

    def setUp(self):
        from tests.native_ui_fixtures import patch_native_runtime
        patch_native_runtime(self)
        # The SVG renderer still runs; PNG encoding is unrelated to this mode.
        self.enterContext(patch("aptaswitch_studio.figure_ui.svg_to_png_bytes", return_value=b"PNG"))

    def test_slider_reorders_saved_metrics_without_changing_scores_or_selection(self):
        app = AppTest.from_function(results_app).run()
        self.assertFalse(app.exception)
        saved = copy.deepcopy(app.session_state.saved_result)
        self.assertEqual(list(app.table[0].value["Essai"]), [1, 2, 3])
        app.slider(key="inverse_extension_performance_weight").set_value(0).run()
        self.assertFalse(app.exception)
        self.assertEqual(list(app.table[0].value["Essai"]), [3, 2, 1])
        self.assertEqual(app.selectbox(key="inverse_extension_selected_trial").value, 1)
        self.assertEqual(app.session_state.saved_result, saved)
        app.slider(key="inverse_extension_performance_weight").set_value(100).run()
        self.assertEqual(list(app.table[0].value["Essai"]), [1, 2, 3])
        help_text = app.slider(key="inverse_extension_performance_weight").help
        self.assertIn("rangs normalisés", help_text)
        self.assertIn("ne mesure pas un rendement expérimental", help_text)

    def test_two_diagrams_use_saved_reference_and_exact_extension_annotations(self):
        with patch.object(ui, "render_extension_structure_svg", wraps=ui.render_extension_structure_svg) as renderer:
            app = AppTest.from_function(results_app)
            app.session_state.aptamer_input = "CCCC"
            app.session_state.trigger_input = "AAAA"
            app.session_state.extension_temperature = 85.0
            app.run()
            self.assertFalse(app.exception)
            self.assertEqual(len(app.get("html")), 2)
            self.assertEqual(len(app.get("download_button")), 5)
            calls = renderer.call_args_list
            self.assertEqual(calls[0].args, ("AAGC", "GCTT", "((((+))))"))
            self.assertEqual(calls[2].args, ("AAGC", "CCGCTTAA", "((((+..)))).."))
            self.assertEqual(calls[2].kwargs["extension_length"], 4)
            self.assertEqual(calls[2].kwargs["extension_3prime_length"], 2)
            app.selectbox(key="inverse_extension_selected_trial").set_value(3).run()
            self.assertEqual(renderer.call_args_list[-2].args[1], "GGGCTTTT")
        self.assertTrue(any("37 °C" in item.value for item in app.caption))
        self.assertFalse(any("85 °C" in item.value for item in app.caption))

    def test_failed_comparisons_remain_visible_and_are_not_ranked_as_conserved(self):
        result = make_result()
        result.inverse_extension_analysis.update(status="failed", analyzed_count=1, error="Native calculation interrupted",
                                                  warnings=["Référence sans paire entre les brins"])
        for candidate in result.candidates[1:]:
            candidate.inverse_extension_metrics = {"status": "failed"}
        app = AppTest.from_function(results_app)
        app.session_state.saved_result = result
        app.run()
        self.assertFalse(app.exception)
        self.assertEqual(list(app.table[0].value["Essai"]), [1])
        self.assertIn("1/3", app.caption[1].value)
        warnings = " ".join(item.value for item in app.warning)
        self.assertIn("ne signifie pas", warnings)
        self.assertIn("Native calculation interrupted", warnings)
        self.assertIn("Référence sans paire", warnings)

    def test_no_analysis_hides_module_and_empty_analysis_does_not_fabricate_structures(self):
        for stage, expected_message in (({}, False), ({"status": "cancelled"}, True)):
            with self.subTest(stage=stage):
                result = make_result()
                result.inverse_extension_analysis = stage
                for candidate in result.candidates:
                    candidate.inverse_extension_metrics = {}
                app = AppTest.from_function(results_app)
                app.session_state.saved_result = result
                app.run()
                self.assertFalse(app.exception)
                self.assertEqual(bool(app.info), expected_message)
                self.assertFalse(app.table)
                self.assertFalse(app.get("html"))

    def test_missing_saved_reference_and_malformed_dot_bracket_show_diagnostic(self):
        for remove_reference in (False, True):
            with self.subTest(remove_reference=remove_reference):
                result = make_result()
                if remove_reference:
                    result.inverse_extension_analysis["reference"] = {}
                else:
                    result.inverse_extension_analysis["reference"]["mfe_structure"] = "((+))"
                app = AppTest.from_function(results_app)
                app.session_state.saved_result = result
                app.session_state.inverse_extension_selected_trial = 999
                app.run()
                self.assertFalse(app.exception)
                self.assertTrue(app.warning)
                self.assertEqual(app.selectbox(key="inverse_extension_selected_trial").value, 1)
                self.assertEqual(len(app.get("html")), 0 if remove_reference else 1)

    def test_reanalysis_only_runs_from_explicit_button_with_requested_top(self):
        app = AppTest.from_function(reanalysis_app).run()
        self.assertFalse(app.exception)
        self.assertNotIn("reanalysis_requested", app.session_state)
        self.assertEqual(app.number_input(key="inverse_extension_reanalyze_top").value, 3)
        saved = copy.deepcopy(app.session_state.saved_result)
        app.number_input(key="inverse_extension_reanalyze_top").set_value(2).run()
        self.assertNotIn("reanalysis_requested", app.session_state)
        app.button(key="inverse_extension_reanalyze").click().run()
        self.assertFalse(app.exception)
        self.assertEqual(app.session_state.reanalysis_requested, ("inverse_ui_fixture", 2))
        self.assertEqual(app.session_state.saved_result, saved)

    def test_historical_reanalysis_click_uses_current_model_without_changing_source(self):
        from aptaswitch_studio import web_app

        result = historical_result()
        original = copy.deepcopy(result)
        app = AppTest.from_function(historical_reanalysis_app)
        app.session_state.saved_result = result
        with patch.object(web_app, "configured_library_path", return_value="/test/native"), patch.object(
            web_app, "DesignRunJob",
        ) as create_job:
            create_job.return_value.snapshot.return_value = SimpleNamespace(state="running")
            app.run()
            self.assertFalse(app.exception)
            captions = " ".join(item.value for item in app.caption)
            self.assertIn("Nouvelle comparaison : moteur hybride expérimental Vienna + Sugimoto", captions)
            self.assertIn("conserve son modèle dna04", captions)
            create_job.assert_not_called()
            app.number_input(key="inverse_extension_reanalyze_top").set_value(2).run()
            app.button(key="inverse_extension_reanalyze").click().run()
            self.assertFalse(app.exception)
            self.assertFalse(app.error)
            create_job.assert_called_once()
            request = create_job.call_args.args[0]
            self.assertEqual(request.inverse_extension.top_candidates, 2)
            self.assertEqual(request.inverse_extension.thermo,
                             replace(original.request.inverse_extension.thermo, material="experimental-Sugimoto95-Vienna04-loopmean"))
            self.assertIs(create_job.call_args.kwargs["analysis_source"], result)
            create_job.return_value.start.assert_called_once()
            self.assertEqual(app.session_state.saved_result, original)

    def test_historical_reanalysis_rejects_incompatible_salts_before_start(self):
        from aptaswitch_studio import web_app

        for magnesium, sodium, message in ((.01, .15, "Mg²⁺ = 0"), (0., .05, "Na⁺ 0.12–1 M")):
            with self.subTest(magnesium=magnesium, sodium=sodium):
                result = historical_result(magnesium, sodium)
                original = copy.deepcopy(result)
                app = AppTest.from_function(historical_reanalysis_app)
                app.session_state.saved_result = result
                with patch.object(web_app, "configured_library_path", return_value="/test/native"), patch.object(
                    web_app, "DesignRunJob",
                ) as create_job:
                    app.run()
                    app.button(key="inverse_extension_reanalyze").click().run()
                    self.assertFalse(app.exception)
                    self.assertIn(message, app.error[0].value)
                    create_job.assert_not_called()
                self.assertEqual(app.session_state.saved_result, original)

    def test_any_active_job_blocks_reanalysis_and_missing_analysis_can_be_requested(self):
        for key in ("design_job", "extension_job", "system_tube_job", "blast_job",
                    "ensemble_job", "saved_switch_analysis_job", "saved_switch_ffs_job"):
            with self.subTest(job=key):
                app = AppTest.from_function(reanalysis_app)
                app.session_state[key] = RunningJob()
                app.run()
                self.assertFalse(app.exception)
                self.assertTrue(app.button(key="inverse_extension_reanalyze").disabled)
        app = AppTest.from_function(reanalysis_app)
        result = make_result()
        result.inverse_extension_analysis = {}
        app.session_state.saved_result = result
        app.run()
        self.assertFalse(app.exception)
        self.assertFalse(app.button(key="inverse_extension_reanalyze").disabled)
        self.assertTrue(app.info)

    def test_csv_export_tracks_current_order_scores_and_weight(self):
        with patch.object(ui, "render_export_menu", wraps=ui.render_export_menu) as export:
            app = AppTest.from_function(results_app).run()
            app.slider(key="inverse_extension_performance_weight").set_value(0).run()
            option = export.call_args.args[0][0]
        self.assertFalse(app.exception)
        self.assertTrue(option.filename.endswith("inverse_extension_comparison.csv"))
        rows = list(csv.DictReader(io.StringIO(option.data.decode("utf-8-sig"))))
        self.assertEqual([row["Essai"] for row in rows], ["3", "2", "1"])
        for row, raw_score in zip(rows, (10, 20, 30)):
            self.assertAlmostEqual(float(row["Score du switch"]), selection_score_100(raw_score))
        self.assertEqual({row["Priorité performance (%)"] for row in rows}, {"0"})

    def test_saved_trigger_linearity_is_displayed_without_adding_a_compromise_term(self):
        result = make_result()
        result.candidates[0].trigger_seq = "CCGCUUAA"
        result.candidates[0].trigger_linearity = {
            "unpaired_fraction": .8, "normalized_defect": .2,
            "structure_mfe": "........", "per_base_unpaired": [.8] * 8,
        }
        saved = copy.deepcopy(result)
        app = AppTest.from_function(results_app)
        app.session_state.saved_result = result
        app.run()
        self.assertFalse(app.exception)
        table = app.table[0].value
        self.assertEqual(list(table["Compromis"]), [.5, .5, .5])
        self.assertEqual(table.iloc[0]["Bases libres du trigger (%)"], 80)
        self.assertEqual(table.iloc[0]["Défaut du trigger isolé (0–1)"], .2)
        self.assertTrue(table.iloc[1:].get("Défaut du trigger isolé (0–1)").isna().all())
        self.assertTrue(any("CCGCUUAA\n........" == item.value for item in app.code))
        self.assertTrue(any("au modèle et aux conditions du design du switch" in item.help for item in app.metric))
        self.assertEqual(app.session_state.saved_result, saved)



if __name__ == "__main__":
    unittest.main()
