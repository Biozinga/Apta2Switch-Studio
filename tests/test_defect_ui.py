"""Regional diagnostics remain distinct from global defects and saved scores."""

from __future__ import annotations

import copy
from dataclasses import replace
import unittest

import pandas as pd
from unittest.mock import patch

from tests.ui_helpers import FrenchAppTest as AppTest
from tests.design_fixtures import design_result_fixture

from aptaswitch_core.architecture import default_architecture_for_trigger
from aptaswitch_core.localization import language_context
from aptaswitch_core.models import FunctionalDefectWeights, RunRequest, SequenceInput
from aptaswitch_studio.results_table_ui import result_rows


def make_result():
    request = RunRequest(
        sequences=SequenceInput("ACGTACGT", "G" * 30),
        architecture=default_architecture_for_trigger(30), trials=1,
        functional_defect_weights=FunctionalDefectWeights(1, 1, 1, 1),
    )
    result = design_result_fixture(request)
    candidate = result.candidates[0]
    candidate.defect = 0.35
    candidate.defect_structural = 0.20
    candidate.defect_concentration = 0.15
    candidate.regional_defects = {
        "off_stem": 0.2, "off_toehold": 0.0, "on_rbs_aug": 0.1, "on_rbs_linker": 0.3,
        "off": 0.1, "on": 0.2, "balanced": 0.15,
        "regions": {
            "off_stem": {"indices_switch": [1, 2, 6, 7], "nucleotide_count": 4,
                         "expected_incorrect_nt": 0.8, "normalized_defect": 0.2, "share_of_switch_defect": 0.5},
            "off_toehold": {"indices_switch": [0, 3, 4], "nucleotide_count": 3,
                            "expected_incorrect_nt": 0.0, "normalized_defect": 0.0, "share_of_switch_defect": 0.0},
            "on_rbs_aug": {"indices_switch": list(range(6)), "nucleotide_count": 6,
                           "expected_incorrect_nt": 0.6, "normalized_defect": 0.1, "share_of_switch_defect": 0.25},
            "on_rbs_linker": {"indices_switch": [8, 9], "nucleotide_count": 2,
                              "expected_incorrect_nt": 0.6, "normalized_defect": 0.3, "share_of_switch_defect": 0.25},
        },
    }
    candidate.regional_defects.update(scheme="functional_v2",
        weights={key: 1. for key in ("off_toehold", "off_stem", "on_rbs_aug", "on_rbs_linker")},
        normalized_weights=request.functional_defect_weights.fractions())
    return result


def result_app():
    import streamlit as st
    from tests.test_defect_ui import make_result
    from aptaswitch_studio.defect_ui import render_defect_method, render_candidate_defects
    from aptaswitch_studio.results_table_ui import result_rows
    if "saved_result" not in st.session_state:
        st.session_state.saved_result = make_result()
    result = st.session_state.saved_result
    render_defect_method(result)
    st.table(result_rows(result))
    render_candidate_defects(result.candidates[0])


def candidate_app():
    from tests.test_defect_ui import make_result
    from aptaswitch_studio.web_app import _render_switch_candidate
    result = make_result()
    _render_switch_candidate(result, result.candidates[0])


class DefectRowsTests(unittest.TestCase):
    def test_table_keeps_global_parts_and_regional_metrics_separate(self):
        result = make_result()
        saved = copy.deepcopy(result)
        with language_context("fr"):
            row = result_rows(result)[0]
        self.assertEqual(row["Defect global"], 0.35)
        self.assertEqual(row["Défaut structure (global)"], 0.2)
        self.assertEqual(row["Défaut concentration (global)"], 0.15)
        self.assertEqual(row["Défaut OFF"], 0.1)
        self.assertEqual(row["Défaut ON"], 0.2)
        self.assertEqual(row["Défaut OFF/ON"], 0.15)
        self.assertEqual(row["Score"], 50.0)
        self.assertEqual(result, saved)


    def test_excluded_candidates_keep_missing_values_without_a_rank_or_score(self):
        result = make_result()
        result.candidates[0] = replace(result.candidates[0], excluded=True, score=None, defect=None, regional_defects={},
                                       exclusion_reason="Codon STOP prématuré — codon 2", has_stop_codon=True, first_stop=2)
        with language_context("fr"):
            row = result_rows(result)[0]
        self.assertIsNone(row["Rang"])
        self.assertIsNone(row["Score"])
        self.assertTrue(all(row[key] is None for key in ("Defect global", "Défaut OFF/ON", "Défaut OFF", "Défaut ON")))


class DefectUITests(unittest.TestCase):
    def test_detail_displays_four_saved_regions_and_explicit_percentages(self):
        app = AppTest.from_function(result_app).run()
        self.assertFalse(app.exception)
        table = app.table[1].value
        self.assertEqual(list(table["État"]), ["OFF", "OFF", "ON", "ON"])
        self.assertEqual(list(table["Défaut d’ensemble (0–1)"]), ["0", "0.2", "0.1", "0.3"])
        self.assertEqual(list(table["Part du défaut du switch (%)"]), ["0 %", "50 %", "25 %", "25 %"])
        self.assertEqual(list(table["Part du compromis (%)"]), ["25 %"] * 4)
        self.assertEqual(table["Positions du switch (à partir de 1)"].iloc[1], "2–3, 7–8")
        self.assertEqual(table["Nucléotides incorrects attendus"].iloc[1], "0.8")
        captions = " ".join(item.value for item in app.caption)
        help_text = " ".join(item.proto.help for item in app.markdown)
        self.assertIn("pas uniquement sur la structure MFE", help_text)
        self.assertIn("Poids effectifs enregistrés : OFF 50 % · ON 50 %", captions)
        self.assertIn("s’additionnent", help_text)
        self.assertIn("mesure indépendante", help_text)
        self.assertIn("AUG et le bulge restent non appariés", help_text)

    def test_saved_weights_override_form_values_and_disabled_state_keeps_local_defects(self):
        result = make_result()
        result.request = replace(result.request, functional_defect_weights=FunctionalDefectWeights(0, 0, 1, 3))
        regional = result.candidates[0].regional_defects
        regional.update(weights={"off_stem": 0., "off_toehold": 0., "on_rbs_aug": 1., "on_rbs_linker": 3.},
                        normalized_weights={"off_stem": 0., "off_toehold": 0., "on_rbs_aug": .25, "on_rbs_linker": .75},
                        off=None, on=.25, balanced=.25)
        saved = copy.deepcopy(result)
        app = AppTest.from_function(result_app)
        app.session_state.saved_result = result
        for key in ("off_stem", "off_toehold", "on_rbs_aug", "on_rbs_linker"):
            app.session_state[f"weight_functional_{key}"] = 20.0
        app.run()
        self.assertFalse(app.exception)
        self.assertTrue(pd.isna(app.table[0].value["Défaut OFF"].iloc[0]))
        self.assertEqual(list(app.table[1].value["Défaut d’ensemble (0–1)"]), ["0", "0.2", "0.1", "0.3"])
        self.assertEqual(list(app.table[1].value["Part du compromis (%)"]), ["0 %", "0 %", "25 %", "75 %"])
        captions = " ".join(item.value for item in app.caption)
        self.assertIn("Poids effectifs enregistrés : OFF 0 % · ON 100 %", captions)
        self.assertIn("OFF est désactivé", captions)
        self.assertEqual(app.session_state.saved_result, saved)

    def test_large_saved_raw_weights_normalize_without_overflow(self):
        result = make_result()
        result.request = replace(result.request, functional_defect_weights=FunctionalDefectWeights(1e308, 1e308, 1e308, 1e308))
        result.candidates[0].regional_defects["weights"] = {
            key: 1e308 for key in ("off_stem", "off_toehold", "on_rbs_aug", "on_rbs_linker")
        }
        app = AppTest.from_function(result_app)
        app.session_state.saved_result = result
        app.run()
        self.assertFalse(app.exception)
        self.assertEqual(list(app.table[1].value["Part du compromis (%)"]), ["25 %"] * 4)
        self.assertFalse(any("est désactivé" in item.value for item in app.caption))

    def test_zero_switch_defect_has_undefined_share_not_zero_percent(self):
        result = make_result()
        regional = result.candidates[0].regional_defects
        for zone in ("on_rbs_aug", "on_rbs_linker"):
            regional[zone] = 0.0
            regional["regions"][zone].update(expected_incorrect_nt=0.0, normalized_defect=0.0, share_of_switch_defect=None)
        regional.update(on=0.0, balanced=0.05)
        app = AppTest.from_function(result_app)
        app.session_state.saved_result = result
        app.run()
        self.assertFalse(app.exception)
        table = app.table[1].value
        self.assertEqual(list(table["Part du défaut du switch (%)"])[2:], ["—", "—"])
        self.assertEqual(list(table["Défaut d’ensemble (0–1)"])[2:], ["0", "0"])


    def test_selected_candidate_view_includes_regional_details_without_recalculation(self):
        with patch("aptaswitch_studio.web_app.render_diagram"):
            app = AppTest.from_function(candidate_app).run()
        self.assertFalse(app.exception)
        self.assertIn("Détail des défauts OFF/ON", [item.label for item in app.expander])
        self.assertEqual(len(app.table[0].value), 4)


if __name__ == "__main__":
    unittest.main()
