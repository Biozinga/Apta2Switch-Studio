"""Inspect ranked results without recomputing or losing screened candidates."""

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


APP_PATH = Path(__file__).resolve().parents[1] / "src/aptaswitch_studio/web_app.py"
TRIGGER = "GTCCAGGCTGGTATAATTAGATCC"


class CandidateNavigationTests(unittest.TestCase):
    def setUp(self):
        self.context = ExitStack()
        self.addCleanup(self.context.close)
        self.output = Path(self.context.enter_context(tempfile.TemporaryDirectory()))
        for target, value in (
            ("aptaswitch_core.nupack_setup.load_settings", {}),
            ("aptaswitch_core.nupack_setup.detect_nupack", None),
            ("aptaswitch_core.nupack_setup.configured_nupack_import_path", None),
            ("aptaswitch_studio.extension_ui.configured_nupack_import_path", None),
        ):
            self.context.enter_context(patch(target, return_value=value))
        for target in (
            "aptaswitch_core.nupack_setup.save_settings",
            "aptaswitch_core.nupack_setup.save_nupack_settings",
            "aptaswitch_core.nupack_engine.run_nupack_engine",
            "aptaswitch_core.nupack_engine.analyze_aptamer_trigger",
            "aptaswitch_studio.web_runtime.DesignRunJob.start",
        ):
            self.context.enter_context(patch(
                target, side_effect=AssertionError("Result navigation must not run an engine or write settings"),
            ))
        self.render = self.context.enter_context(patch(
            "aptaswitch_core.preview.candidate_structure_svg", wraps=candidate_structure_svg,
        ))

    def result(self):
        request = RunRequest(
            sequences=SequenceInput("ACGTACGT", TRIGGER),
            architecture=default_architecture_for_trigger(len(TRIGGER)),
            trials=3, output_dir=self.output,
        )
        result = design_result_fixture(request)
        # Deliberately unsorted trial IDs and distinguishable recorded data:
        # navigation follows the persisted ranking, not the order of trials.
        candidates = [
            replace(candidate, trial=trial, score=score,
                    switch_seq=base + candidate.switch_seq[1:])
            for candidate, trial, score, base in zip(
                result.candidates, (40, 7, 19), (90.0, 70.0, 60.0), "ACG",
            )
        ]
        candidates[1] = replace(
            candidates[1],
            structure_off="." * len(candidates[1].structure_off),
            structure_on="".join("+" if char == "+" else "." for char in candidates[1].structure_on),
        )
        return replace(result, candidates=candidates, warnings=[])

    def app(self, result):
        app = AppTest.from_file(str(APP_PATH), default_timeout=30)
        for key, value in {
            "dependency_detection_done": True,
            "output_dir": str(self.output),
            "nav": "Résultats",
            "last_result": result,
            "last_exports": {},
        }.items():
            app.session_state[key] = value
        app.run()
        self.assertFalse(app.exception)
        return app

    def assert_candidate(self, app, candidate, *, state="on"):
        self.assertFalse(app.exception)
        self.assertEqual(app.selectbox(key="selected_trial").value, candidate.trial)
        self.assertEqual(app.session_state["selected_trial"], candidate.trial)
        self.assertIn(candidate.switch_seq, [block.value for block in app.code])
        self.assertIn(candidate.trigger_seq, [block.value for block in app.code])
        call = self.render.call_args
        self.assertEqual(call.args[1:3], (candidate.trigger_seq, candidate.switch_seq))
        self.assertEqual(call.args[3], getattr(candidate, f"structure_{state}"))
        self.assertEqual(call.kwargs["state"], state)

    def test_arrows_follow_rank_and_stay_synchronized_with_selector(self):
        result = self.result()
        app = self.app(result)
        self.assertTrue(app.button(key="previous_candidate").disabled)
        self.assertFalse(app.button(key="next_candidate").disabled)
        app.radio(key="candidate_state").set_value("ON + trigger").run()
        self.assert_candidate(app, result.candidates[0])

        app.button(key="next_candidate").click().run()
        self.assert_candidate(app, result.candidates[1])
        self.assertFalse(app.button(key="previous_candidate").disabled)
        self.assertFalse(app.button(key="next_candidate").disabled)

        app.button(key="next_candidate").click().run()
        self.assert_candidate(app, result.candidates[2])
        self.assertTrue(app.button(key="next_candidate").disabled)
        app.button(key="previous_candidate").click().run()
        self.assert_candidate(app, result.candidates[1])

        app.selectbox(key="selected_trial").set_value(40).run()
        self.assertTrue(app.button(key="previous_candidate").disabled)
        app.button(key="next_candidate").click().run()
        self.assert_candidate(app, result.candidates[1])

    def test_arrows_reach_both_exclusions_and_reuse_only_recorded_structures(self):
        result = self.result()
        valid, nonlinear, stop = result.candidates
        structure = list(nonlinear.structure_on)
        loop_start = len(TRIGGER) + 1 + int(architecture_layout(result.request.architecture)["rbs_linker_start"])
        structure[loop_start], structure[loop_start + 4] = "(", ")"
        screened = dict(score=None, defect=None, on_yield_pct=None, leak_pct=None,
                        delta_g_off=None, ddg_activation=None, delta_g_rbs_linker=None,
                        has_bad_rbs_linker=None, structure_off="", excluded=True)
        nonlinear = replace(nonlinear, **screened, structure_on="".join(structure),
                            exclusion_reason="RBS–linker non linéaire en ON")
        stop = replace(stop, **screened, delta_g_on=None, structure_on="",
                       has_stop_codon=True, first_stop=2,
                       exclusion_reason="Codon STOP prématuré — codon 2")
        result = replace(result, candidates=[valid, nonlinear, stop])
        app = self.app(result)

        app.button(key="next_candidate").click().run()
        self.assert_candidate(app, nonlinear)
        self.assertTrue(any(nonlinear.exclusion_reason in item.value for item in app.error))

        self.render.reset_mock()
        app.button(key="next_candidate").click().run()
        self.assertFalse(app.exception)
        self.render.assert_not_called()
        self.assertEqual(app.selectbox(key="selected_trial").value, stop.trial)
        self.assertIn(stop.switch_seq, [block.value for block in app.code])
        self.assertTrue(any(stop.exclusion_reason in item.value for item in app.error))
        self.assertTrue(app.button(key="next_candidate").disabled)
        app.button(key="previous_candidate").click().run()
        self.assert_candidate(app, nonlinear)

    def test_empty_and_single_candidate_bounds(self):
        result = self.result()
        empty = self.app(replace(result, candidates=[]))
        self.assertNotIn("selected_trial", {item.key for item in empty.selectbox})
        self.assertFalse({"previous_candidate", "next_candidate"} & {item.key for item in empty.button})
        single = self.app(replace(result, candidates=result.candidates[:1]))
        self.assertTrue(single.button(key="previous_candidate").disabled)
        self.assertTrue(single.button(key="next_candidate").disabled)

    def test_no_stop_success_notice_but_actual_stop_remains_visible(self):
        result = self.result()
        app = self.app(result)
        self.assertFalse(app.success)
        self.assertFalse(any("STOP" in item.value or "Cadre de lecture valide" in item.value
                             for item in (*app.info, *app.caption, *app.warning, *app.error)))
        stop = replace(result.candidates[1], has_stop_codon=True, first_stop=2)
        app.session_state["last_result"] = replace(result, candidates=[result.candidates[0], stop])
        app.button(key="next_candidate").click().run()
        self.assert_candidate(app, stop, state="off")
        self.assertTrue(any("Codon STOP prématuré — codon 2" in item.value for item in app.error))


if __name__ == "__main__":
    unittest.main()
