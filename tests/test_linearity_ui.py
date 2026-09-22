"""The optional ON-linearity screen must remain distinct from STOP exclusion."""

from __future__ import annotations

import tempfile
import unittest
from contextlib import ExitStack
from dataclasses import replace
from pathlib import Path
from unittest.mock import patch

from tests.ui_helpers import FrenchAppTest as AppTest

from aptaswitch_core.architecture import architecture_layout, default_architecture_for_trigger
from aptaswitch_core.models import RunRequest, SequenceInput
from aptaswitch_core.preview import candidate_structure_svg
from tests.design_fixtures import design_result_fixture


APP_PATH = Path(__file__).resolve().parents[1] / "src" / "aptaswitch_studio" / "web_app.py"
TRIGGER = "GTCCAGGCTGGTATAATTAGATCC"
REASON = "RBS–linker non linéaire en ON"


class LinearityUITests(unittest.TestCase):
    def setUp(self):
        self.context = ExitStack()
        self.addCleanup(self.context.close)
        self.output = Path(self.context.enter_context(tempfile.TemporaryDirectory()))
        for target, value in (
            ("aptaswitch_core.nupack_setup.load_settings", {}),
            ("aptaswitch_core.nupack_setup.detect_nupack", None),
            ("aptaswitch_core.nupack_setup.configured_nupack_import_path", "/test/nupack"),
            ("aptaswitch_studio.extension_ui.configured_nupack_import_path", None),
        ):
            self.context.enter_context(patch(target, return_value=value))
        for target in ("save_settings", "save_nupack_settings"):
            self.context.enter_context(patch(
                f"aptaswitch_core.nupack_setup.{target}",
                side_effect=AssertionError("UI tests must not write settings"),
            ))
        self.start = self.context.enter_context(patch(
            "aptaswitch_studio.web_runtime.DesignRunJob.start", autospec=True,
        ))

    def app(self):
        app = AppTest.from_file(str(APP_PATH), default_timeout=30)
        app.session_state["dependency_detection_done"] = True
        app.session_state["output_dir"] = str(self.output)
        app.session_state["aptamer_input"] = "ACGTACGT"
        app.session_state["trigger_input"] = TRIGGER
        app.run()
        self.assertFalse(app.exception)
        return app

    def launch(self, app):
        count = self.start.call_count
        next(button for button in app.button if button.label == "Lancer la conception").click().run()
        self.assertFalse(app.exception)
        self.assertEqual(self.start.call_count, count + 1)
        return self.start.call_args.args[0].request

    def result(self):
        request = RunRequest(
            sequences=SequenceInput("ACGTACGT", TRIGGER),
            architecture=default_architecture_for_trigger(len(TRIGGER)),
            trials=3, output_dir=self.output,
        )
        result = design_result_fixture(request)
        valid, non_linear, stop = result.candidates
        # Give the recorded ON structure a hairpin inside the RBS-loop region.
        # These are synthetic rendering fixtures, not a thermodynamic claim.
        layout = architecture_layout(request.architecture)
        structure = list(non_linear.structure_on)
        start = len(TRIGGER) + 1 + int(layout["rbs_linker_start"])
        structure[start], structure[start + 4] = "(", ")"
        non_linear = replace(
            non_linear, trial=1, score=None, defect=None, on_yield_pct=None,
            leak_pct=None, delta_g_off=None, ddg_activation=None,
            delta_g_rbs_linker=None, has_bad_rbs_linker=None,
            structure_off="", structure_on="".join(structure),
            excluded=True, exclusion_reason=REASON,
        )
        stop = replace(
            stop, trial=2, score=None, defect=None, on_yield_pct=None,
            leak_pct=None, delta_g_off=None, delta_g_on=None,
            ddg_activation=None, delta_g_rbs_linker=None,
            has_bad_rbs_linker=None, has_stop_codon=True, first_stop=2,
            structure_off="", structure_on="", excluded=True,
            exclusion_reason="Codon STOP prématuré — codon 2",
        )
        return replace(result, candidates=[replace(valid, trial=3), non_linear, stop], warnings=[])

    def test_default_filter_is_in_scoring_and_reaches_run_request(self):
        app = self.app()
        scoring = next(tab for tab in app.tabs if tab.label == "Scoring")
        self.assertTrue(scoring.checkbox(key="exclude_non_linear_candidates").value)
        request = self.launch(app)
        self.assertTrue(request.exclude_non_linear_candidates)
        self.assertTrue(request.normalized().exclude_non_linear_candidates)
        self.assertTrue(request.to_dict()["exclude_non_linear_candidates"])

    def test_opt_out_is_silent_and_survives_navigation_and_submission(self):
        app = self.app()
        warnings_before = [warning.value for warning in app.warning]
        app.checkbox(key="exclude_non_linear_candidates").uncheck().run()
        self.assertFalse(app.exception)
        self.assertEqual([warning.value for warning in app.warning], warnings_before)
        self.assertTrue(app.checkbox(key="exclude_stop_candidates").value)
        app.radio(key="nav").set_value("Résultats").run()
        app.radio(key="nav").set_value("Conception").run()
        self.assertFalse(app.exception)
        self.assertFalse(app.checkbox(key="exclude_non_linear_candidates").value)
        self.assertEqual([warning.value for warning in app.warning], warnings_before)
        request = self.launch(app)
        self.assertFalse(request.exclude_non_linear_candidates)
        self.assertFalse(request.normalized().exclude_non_linear_candidates)
        self.assertFalse(request.to_dict()["exclude_non_linear_candidates"])

    def test_mixed_exclusions_keep_reasons_and_reuse_existing_on_structure(self):
        result = self.result()
        app = self.app()
        app.session_state["last_result"] = result
        app.session_state["selected_trial"] = 1
        with patch("aptaswitch_core.preview.candidate_structure_svg", wraps=candidate_structure_svg) as render, patch(
            "aptaswitch_core.nupack_engine.run_nupack_engine",
            side_effect=AssertionError("Inspecting a candidate must not recompute thermodynamics"),
        ):
            app.radio(key="nav").set_value("Résultats").run()
        self.assertFalse(app.exception)
        rows = app.get("table")[0].value
        self.assertEqual(list(rows["Essai"]), [3, 1, 2])
        row = rows[rows["Essai"] == 1].iloc[0]
        self.assertIn(REASON, row["Statut"])
        self.assertNotIn("STOP", row["Statut"])
        self.assertEqual(row["Rang"], "—")
        self.assertEqual(row["Score"], "—")
        self.assertEqual(row["STOP"], "non")
        self.assertEqual(row["ΔG ON"], str(result.candidates[1].delta_g_on))
        self.assertTrue(any(REASON in option for option in app.selectbox(key="selected_trial").options))
        self.assertTrue(any(REASON in item.value for item in app.error))
        self.assertTrue(render.called)
        self.assertTrue(all(call.kwargs["state"] == "on" for call in render.call_args_list))
        self.assertTrue(all(call.args[3] == result.candidates[1].structure_on for call in render.call_args_list))
        self.assertTrue(any(button.label == "Télécharger la structure PNG" for button in app.get("download_button")))
        self.start.assert_not_called()

        # The earlier STOP screen has no calculated structure to display.
        with patch("aptaswitch_core.preview.candidate_structure_svg") as render:
            app.selectbox(key="selected_trial").set_value(2).run()
        self.assertFalse(app.exception)
        render.assert_not_called()
        self.assertTrue(any("STOP" in item.value for item in app.error))
        self.assertFalse(any(button.label == "Télécharger la structure PNG" for button in app.get("download_button")))

    def test_all_non_linear_exclusions_are_not_mislabeled_as_stops(self):
        result = self.result()
        app = self.app()
        app.session_state["last_result"] = replace(result, candidates=[result.candidates[1]])
        app.radio(key="nav").set_value("Résultats").run()
        self.assertFalse(app.exception)
        notices = [item.value for item in (*app.info, *app.warning, *app.error)]
        self.assertTrue(any("Tous les candidats sont exclus" in text for text in notices))
        self.assertTrue(any(REASON in text for text in notices))
        self.assertFalse(any("STOP" in text for text in notices))
        self.assertEqual(app.get("table")[0].value.iloc[0]["Score"], "—")


if __name__ == "__main__":
    unittest.main()
