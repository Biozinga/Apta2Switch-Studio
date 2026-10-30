"""Result summaries keep the saved ranking and candidate selection stays separate."""

from __future__ import annotations

from dataclasses import replace
from pathlib import Path
from types import SimpleNamespace
import unittest
from unittest.mock import patch

from aptaswitch_core.localization import get_language, set_language
from aptaswitch_core.models import DesignCandidate
from aptaswitch_studio.results_table_ui import (
    _PENDING_SELECTION, _apply_pending_selection, _open_full_table, _table_identity,
    best_candidates, candidate_status, result_rows,
)
from tests.ui_helpers import FrenchAppTest as AppTest


def sample_result(count=41, *, run_id="recorded-run"):
    candidate = DesignCandidate(
        trial=40, switch_seq="AUG", trigger_seq="CAT", score=99.0,
        defect=.2, on_yield_pct=10.0, leak_pct=2.0,
        delta_g_off=-9.0, delta_g_on=-10.0, ddg_activation=-1.0,
        delta_g_rbs_linker=-.1, has_bad_rbs_linker=False, has_stop_codon=False,
        defect_structural=.15, defect_concentration=.05,
        regional_defects={"off": .1, "on": .2, "balanced": .15},
    )
    candidates = [replace(candidate, trial=max(1000, count * 7) - 7 * index, score=max(0.0, 99.0 - index))
                  for index in range(count)]
    return SimpleNamespace(run_id=run_id, candidates=candidates, output_dir=Path("/tmp") / run_id,
                           request=SimpleNamespace(score_scale="bounded_v1"))


def table_app():
    import streamlit as st
    from aptaswitch_studio.results_table_ui import render_results_table
    from tests.test_results_table_ui import sample_result
    result = st.session_state.get("test_result", sample_result())
    render_results_table(result)
    choices = [candidate.trial for candidate in result.candidates]
    if choices:
        if st.session_state.get("selected_trial") not in choices:
            st.session_state.selected_trial = choices[0]
        st.selectbox("Candidat", choices, key="selected_trial")


def table_dialog_app():
    import streamlit as st
    from aptaswitch_studio.results_table_ui import (
        _PENDING_SELECTION, _open_full_table, _table_identity, render_results_table,
    )
    from tests.test_results_table_ui import sample_result
    result = st.session_state.get("test_result", sample_result())
    opening_candidate = _PENDING_SELECTION in st.session_state
    render_results_table(result)
    if not opening_candidate:
        # AppTest re-executes the complete script for dialog widget events.
        # Keep the modal open here to exercise its separate selector/action handling.
        _open_full_table(result, _table_identity(result))
    choices = [candidate.trial for candidate in result.candidates]
    if st.session_state.get("selected_trial") not in choices:
        st.session_state.selected_trial = choices[0]
    st.selectbox("Candidat", choices, key="selected_trial")




class ResultRowsTests(unittest.TestCase):
    def setUp(self):
        set_language("fr")

    def test_five_percent_rounds_up_and_keeps_recorded_order(self):
        for count, expected in ((0, 0), (1, 1), (19, 1), (20, 1), (21, 2), (40, 2), (41, 3)):
            with self.subTest(count=count):
                result = sample_result(count)
                self.assertEqual(best_candidates(result), result.candidates[:expected])

    def test_excluded_candidates_do_not_count_toward_top_five_percent(self):
        result = sample_result(41)
        result.candidates[0] = replace(result.candidates[0], excluded=True, score=None)
        result.candidates[20] = replace(result.candidates[20], excluded=True, score=None)
        self.assertEqual([candidate.trial for candidate in best_candidates(result)],
                         [candidate.trial for candidate in result.candidates[1:3]])
        rows = result_rows(result)
        self.assertEqual(len(rows), 41)
        self.assertIsNone(rows[0]["Rang"])
        self.assertIsNone(rows[0]["Score"])
        self.assertEqual(rows[1]["Rang"], 1)
        self.assertEqual(rows[21]["Rang"], 20)

    def test_all_excluded_has_no_best_but_preserves_full_table(self):
        result = sample_result(5)
        result.candidates = [replace(candidate, excluded=True, score=None) for candidate in result.candidates]
        self.assertEqual(result_rows(result, compact=True), [])
        self.assertEqual(len(result_rows(result)), 5)

    def test_metrics_stay_numeric_without_rescoring(self):
        result = sample_result(21)
        result.candidates[:2] = [replace(result.candidates[0], score=9.0),
                                 replace(result.candidates[1], score=80.0)]
        with patch("aptaswitch_core.scoring.compute_selection_score", side_effect=AssertionError("No rescoring")):
            rows = result_rows(result)
            compact = result_rows(result, compact=True)
        self.assertEqual([row["Score"] for row in compact], [9.0, 80.0])
        self.assertEqual(sorted(row["Score"] for row in compact), [9.0, 80.0])
        for name, expected in (("Defect global", .2), ("Défaut structure (global)", .15),
                               ("Défaut concentration (global)", .05), ("Défaut OFF/ON", .15)):
            self.assertEqual(rows[0][name], expected)
        self.assertEqual(list(compact[0]), ["Rang", "Essai", "Score", "ON %", "Trigger–switch avec aptamère (%)", "Défaut OFF/ON", "Defect global"])

    def test_uncomputed_is_not_a_zero(self):
        result = sample_result(1)
        result.candidates[0] = replace(result.candidates[0], defect=None, defect_structural=0.0,
                                       defect_concentration=None, regional_defects={})
        row = result_rows(result)[0]
        self.assertIsNone(row["Defect global"])
        self.assertEqual(row["Défaut structure (global)"], 0.0)
        self.assertIsNone(row["Défaut concentration (global)"])
        self.assertIsNone(row["Défaut OFF/ON"])

    def test_absent_opening_measurement_is_blank_without_changing_score(self):
        result = sample_result(1)
        del result.candidates[0].rbs_accessibility
        row = result_rows(result)[0]
        self.assertIsNone(row["RBS → AUG entièrement libre ON (%)"])
        self.assertIsNone(row["ΔG ouverture RBS → AUG ON"])
        self.assertEqual(row["Score"], 99.)
        self.assertEqual(result_rows(result, compact=True)[0]["Score"], 99.)

    def test_measured_zero_opening_probability_is_preserved(self):
        result = sample_result(1)
        result.candidates[0].rbs_accessibility = {
            "probability": 0., "opening_free_energy_kcal_mol": None,
        }
        self.assertEqual(result_rows(result)[0]["RBS → AUG entièrement libre ON (%)"], 0.)

    def test_pending_selection_from_a_different_run_is_discarded(self):
        result = sample_result()
        previous = sample_result(run_id="previous")
        state = {_PENDING_SELECTION: {"identity": _table_identity(previous), "trial": result.candidates[1].trial}}
        with patch("aptaswitch_studio.results_table_ui.st.session_state", state):
            _apply_pending_selection(result, _table_identity(result))
        self.assertEqual(state, {})

    def test_status_preserves_screening_reason(self):
        excluded = replace(sample_result(1).candidates[0], excluded=True,
                           exclusion_reason="Codon STOP prématuré — codon 12")
        self.assertEqual(candidate_status(excluded), "Exclu — STOP prématuré")

    def test_dialog_restores_session_language_on_its_fragment_rerun(self):
        result = sample_result()
        set_language("en")
        with patch("aptaswitch_studio.results_table_ui.st.dialog", side_effect=lambda *args, **kwargs: lambda function: function), \
                patch("aptaswitch_studio.results_table_ui.st.session_state", {"ui_language": "fr"}), \
                patch("aptaswitch_studio.results_table_ui._render_paginated_table"), \
                patch("aptaswitch_studio.results_table_ui.st.caption") as caption, \
                patch("aptaswitch_studio.results_table_ui.st.selectbox", return_value=result.candidates[0].trial) as selectbox, \
                patch("aptaswitch_studio.results_table_ui.st.button", return_value=False):
            _open_full_table(result, _table_identity(result))
        self.assertEqual(get_language(), "fr")
        self.assertEqual(caption.call_args.args[0], "Tous les candidats du classement enregistré, y compris les exclus.")
        self.assertEqual(selectbox.call_args.args[0], "Candidat à explorer")


class ResultsTableAppTests(unittest.TestCase):
    def app(self, result=None):
        app = AppTest.from_function(table_app)
        if result is not None:
            app.session_state["test_result"] = result
        app.run()
        self.assertFalse(app.exception)
        return app

    def test_compact_table_and_expand_complete_table(self):
        result = sample_result(41)
        result.candidates[-1] = replace(result.candidates[-1], excluded=True, score=None)
        app = self.app(result)
        self.assertEqual(len(app.table[0].value), 2)
        app.button(key=f"results_table_expand_{_table_identity(result)}").click().run()
        self.assertFalse(app.exception)
        full = app.table[-1].value
        self.assertEqual(len(full), 41)
        self.assertIn("Défaut structure (global)", full.columns)
        self.assertIn("Défaut concentration (global)", full.columns)
        self.assertTrue(full.iloc[-1]["Statut"].startswith("Exclu"))
        self.assertFalse(app.dataframe)
        self.assertEqual(app.selectbox(key=f"results_table_choose_candidate_{_table_identity(result)}").value, result.candidates[0].trial)

    def test_static_preview_leaves_candidate_navigation_to_separate_selector(self):
        result = sample_result()
        app = self.app(result)
        self.assertFalse(app.dataframe)
        self.assertEqual(len(app.table), 1)
        self.assertEqual(app.selectbox(key="selected_trial").value, result.candidates[0].trial)
        app.selectbox(key="selected_trial").set_value(result.candidates[6].trial).run()
        self.assertFalse(app.exception)
        self.assertEqual(app.selectbox(key="selected_trial").value, result.candidates[6].trial)
        self.assertNotIn(_PENDING_SELECTION, app.session_state)

    def test_new_run_does_not_reuse_previous_candidate(self):
        result = sample_result()
        app = self.app(result)
        app.selectbox(key="selected_trial").set_value(result.candidates[1].trial).run()
        new = sample_result(41, run_id="next-run")
        new.candidates = [replace(candidate, trial=candidate.trial + 10000) for candidate in new.candidates]
        app.session_state["test_result"] = new
        app.run()
        self.assertFalse(app.exception)
        self.assertEqual(app.selectbox(key="selected_trial").value, new.candidates[0].trial)
        self.assertFalse(app.dataframe)

    def test_full_table_uses_explicit_candidate_selector_and_closes_dialog(self):
        result = sample_result()
        result.candidates[-1] = replace(result.candidates[-1], excluded=True, score=None)
        app = AppTest.from_function(table_dialog_app)
        app.session_state["test_result"] = result
        app.run()
        self.assertFalse(app.exception)
        self.assertFalse(app.dataframe)
        app.selectbox(key=f"results_table_choose_candidate_{_table_identity(result)}").set_value(result.candidates[-1].trial).run()
        self.assertFalse(app.exception)
        app.button(key=f"results_table_open_candidate_{_table_identity(result)}").click().run()
        self.assertFalse(app.exception)
        self.assertEqual(app.selectbox(key="selected_trial").value, result.candidates[-1].trial)
        self.assertEqual(len(app.table), 1)
        self.assertNotIn(_PENDING_SELECTION, app.session_state)

    def test_full_table_pages_preserve_recorded_order_and_global_candidate_choice(self):
        result = sample_result(450)
        result.candidates[210] = replace(result.candidates[210], excluded=True, score=None)
        original_candidates = list(result.candidates)
        app = AppTest.from_function(table_dialog_app)
        app.session_state["test_result"] = result
        app.run()
        identity = _table_identity(result)
        self.assertFalse(app.exception)
        self.assertEqual(len(app.table[-1].value), 200)
        app.number_input(key=f"results_table_full_{identity}_page").set_value(2).run()
        self.assertFalse(app.exception)
        frame = app.table[-1].value
        self.assertEqual(frame["Essai"].tolist(), [candidate.trial for candidate in result.candidates[200:400]])
        self.assertEqual(frame.iloc[0]["Rang"], 201)
        self.assertTrue(frame.iloc[10]["Statut"].startswith("Exclu"))
        self.assertEqual(app.selectbox(key="selected_trial").value, result.candidates[0].trial)
        chooser = app.selectbox(key=f"results_table_choose_candidate_{identity}")
        self.assertEqual(len(chooser.options), 450)
        self.assertEqual(chooser.value, result.candidates[0].trial)
        self.assertEqual(result.candidates, original_candidates)
        self.assertFalse(app.dataframe)

    def test_large_top_five_percent_preview_is_paginated(self):
        result = sample_result(4021)
        app = self.app(result)
        self.assertEqual(len(app.table[0].value), 200)
        app.number_input(key=f"results_table_best_{_table_identity(result)}_page").set_value(2).run()
        self.assertFalse(app.exception)
        self.assertEqual(app.table[0].value["Essai"].tolist(),
                         [candidate.trial for candidate in result.candidates[200:202]])
        self.assertEqual(app.selectbox(key="selected_trial").value, result.candidates[0].trial)

    def test_empty_and_all_excluded_keep_interface_accessible(self):
        empty = self.app(sample_result(0))
        self.assertEqual(len(empty.table), 0)
        self.assertEqual(len(empty.button), 0)
        excluded = sample_result(2)
        excluded.candidates = [replace(candidate, excluded=True, score=None) for candidate in excluded.candidates]
        app = self.app(excluded)
        self.assertEqual(len(app.table), 0)
        app.button(key=f"results_table_expand_{_table_identity(excluded)}").click().run()
        self.assertFalse(app.exception)
        self.assertEqual(len(app.table[-1].value), 2)


if __name__ == "__main__":
    unittest.main()
