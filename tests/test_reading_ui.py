"""Full-app reading journeys use exact inputs and the existing saved switches."""

from __future__ import annotations

import json
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import Mock, patch

from streamlit.testing.v1 import AppTest

from aptaswitch_core.architecture import architecture_layout, default_architecture_for_trigger
from aptaswitch_core.models import REPORTERS
from aptaswitch_core.preview import build_switch_model
from aptaswitch_core.saved_switch_analysis import build_saved_switch_context
from aptaswitch_core.thermodynamics.setup import ENGINE_ID, MODEL_ID
from aptaswitch_studio.design_library import create_design, load_design
from aptaswitch_studio.saved_switch_library import list_saved_switches, saved_switch_root
from aptaswitch_studio.web_runtime import JobSnapshot
from tests.design_fixtures import design_result_fixture
from tests.native_ui_fixtures import patch_native_runtime
from tests.ui_helpers import prepare_test_design


APP_PATH = Path(__file__).resolve().parents[1] / "src/aptaswitch_studio/web_app.py"
APTAMER = "AACCGGTTAACCGG"
TRIGGER = "GTCCAGGCTGGTATAATTAGATCCACGTAC"


def supplied_switch() -> str:
    architecture = default_architecture_for_trigger(len(TRIGGER))
    model = build_switch_model(
        architecture, TRIGGER,
        reporter_sequence=REPORTERS["sfgfp"].sequence_after_start_rna,
    )
    sequence = list(model.switch_sequence.replace("X", "A"))
    # A supplied free base deliberately differs from the rendering fixture.
    sequence[architecture_layout(architecture)["rbs_prefix_start"]] = "C"
    return "".join(sequence)


SWITCH = supplied_switch()


class ReadingUITests(unittest.TestCase):
    def setUp(self):
        patch_native_runtime(self)
        self.root = Path(self.enterContext(tempfile.TemporaryDirectory())).resolve()
        self.project = create_design(self.root, "Reading UI")
        for target, value in (
            ("aptaswitch_core.settings.load_settings", {}),
            ("aptaswitch_core.thermodynamics.setup.configured_library_path", "/test/native"),
            ("aptaswitch_studio.extension_ui.configured_library_path", None),
        ):
            self.enterContext(patch(target, return_value=value))
        for target in ("save_settings",):
            self.enterContext(patch(
                "aptaswitch_core.settings." + target,
                side_effect=AssertionError("Reading UI tests must not modify user settings"),
            ))
        self.native = self.enterContext(patch(
            "aptaswitch_core.thermodynamics.design_engine.optimize",
            side_effect=AssertionError("Viewing supplied sequences must not calculate"),
        ))
        self.start = self.enterContext(patch(
            "aptaswitch_studio.web_runtime.DesignRunJob.start", autospec=True,
        ))

    def app(self, *, language="fr", project=None, **state):
        app = AppTest.from_file(str(APP_PATH), default_timeout=30)
        prepare_test_design(app, project=project or self.project)
        app.session_state["dependency_detection_done"] = True
        app.session_state["ui_language"] = language
        app.session_state["design_mode"] = "reading"
        for key, value in state.items():
            app.session_state[key] = value
        app.run()
        self.assertFalse(app.exception)
        self.assertFalse(app.error)
        return app

    def valid_app(self, **state):
        return self.app(aptamer_input=APTAMER, trigger_input=TRIGGER,
                        reading_switch_input=SWITCH, **state)

    def button(self, app, label="Calculer les métriques"):
        return next(button for button in app.button if button.label == label)

    def launch(self, app):
        self.assertFalse(self.button(app).disabled)
        self.button(app).click().run()
        self.assertFalse(app.exception)
        self.assertFalse(app.error)
        self.start.assert_called_once()
        return self.start.call_args.args[0]

    def test_mode_switch_and_navigation_preserve_three_sequences_then_reopen_design(self):
        app = self.valid_app()
        self.assertEqual(app.radio(key="design_mode").options, ["Conception", "Lecture de séquences"])
        app.radio(key="design_mode").set_value("design").run()
        self.assertFalse(app.exception)
        self.assertEqual(app.text_area(key="aptamer_input").value, APTAMER)
        app.radio(key="design_mode").set_value("reading").run()
        self.assertEqual(app.text_area(key="reading_switch_input").value, SWITCH)
        app.radio(key="nav").set_value("Réglages").run()
        app.radio(key="nav").set_value("Conception").run()
        self.assertEqual(app.text_area(key="reading_switch_input").value, SWITCH)
        app.button(key="design_home").click().run()
        saved = load_design(self.project.path)
        self.assertEqual(saved.state["form"]["design_mode"], "reading")
        self.assertEqual(saved.state["form"]["reading_switch_input"], SWITCH)
        other = create_design(self.root, "Other project")
        app.run()
        app.selectbox(key="design_to_open").set_value(str(other.path))
        app.button(key="open_design").click().run()
        self.assertEqual(app.radio(key="design_mode").value, "design")
        self.assertEqual(app.text_area(key="aptamer_input").value, "")
        app.button(key="design_home").click().run()
        app.selectbox(key="design_to_open").set_value(str(self.project.path))
        app.button(key="open_design").click().run()
        self.assertFalse(app.exception)
        self.assertEqual(app.radio(key="design_mode").value, "reading")
        for key, expected in (("aptamer_input", APTAMER), ("trigger_input", TRIGGER),
                              ("reading_switch_input", SWITCH)):
            self.assertEqual(app.text_area(key=key).value, expected)
        self.native.assert_not_called()
        self.start.assert_not_called()

    def test_dispatches_exact_fixed_switch_once_without_inverse_or_exclusion_filters(self):
        raw_switch = "  " + SWITCH[:40].lower() + "\n" + SWITCH[40:].lower() + "  "
        app = self.app(
            aptamer_input=APTAMER.lower(), trigger_input=TRIGGER.lower(),
            reading_switch_input=raw_switch, extension_enabled=True,
            inverse_extension_enabled=True, inverse_extension_target_length=60,
            customize_tswitch=True, toehold_length=12, binding_start=4, binding_end=12,
            trials=100, exclude_stop_candidates=True, exclude_non_linear_candidates=True,
        )
        self.assertFalse(any(widget.key == "trials" for widget in app.number_input))
        self.assertFalse(any(widget.key in {
            "optimization_max_evaluations", "optimization_population_size",
            "optimization_lineages", "optimization_local_tolerance",
        } for widget in app.number_input))
        self.assertFalse(any(widget.key == "long_run_confirm" for widget in app.checkbox))
        job = self.launch(app)
        request = job.request
        self.assertEqual(request.fixed_switch_sequence, SWITCH)
        self.assertEqual(request.engine, ENGINE_ID)
        self.assertEqual(request.thermo.material, MODEL_ID)
        self.assertEqual(request.thermodynamics_library_path, "/test/native")
        self.assertEqual(request.sequences.aptamer_dna, APTAMER)
        self.assertEqual(request.sequences.trigger_dna, TRIGGER)
        self.assertEqual(request.sequences.effective_trigger(), TRIGGER)
        self.assertEqual(request.sequences.binding_start, 0)
        self.assertEqual(request.sequences.binding_length, len(TRIGGER))
        self.assertEqual(request.trials, 1)
        self.assertIsNone(request.inverse_extension)
        self.assertFalse(request.exclude_stop_candidates)
        self.assertFalse(request.exclude_non_linear_candidates)
        self.assertEqual(request.output_dir, self.project.path / "calculations")
        self.assertEqual(app.session_state["reading_switch_input"], raw_switch)
        self.assertEqual(app.radio(key="nav").value, "Résultats")
        self.native.assert_not_called()

    def test_invalid_switch_disables_analysis_and_reports_validation(self):
        for sequence in (SWITCH + "?", SWITCH.replace("AGAGGAGA", "CCCCCCCC", 1)):
            with self.subTest(sequence=sequence[-15:]):
                app = self.app(aptamer_input=APTAMER, trigger_input=TRIGGER,
                               reading_switch_input=sequence)
                self.assertTrue(self.button(app).disabled)
                self.assertTrue(app.warning)
                self.assertFalse(app.exception)
        self.native.assert_not_called()
        self.start.assert_not_called()

    def test_architecture_is_automatic_and_hidden_design_settings_are_ignored(self):
        app = self.valid_app(customize_tswitch=True, toehold_length=7,
                             t7_leader="AAA", rbs_loop="CCCC", rbs_aug_distance=11,
                             linker_auto=False, frame_linker="ACGACG")
        hidden = {"t7_leader", "rbs_loop", "toehold_length", "rbs_aug_distance",
                  "customize_tswitch", "frame_linker", "linker_auto", "binding_start", "binding_end"}
        rendered = {widget.key for kind in (app.text_input, app.number_input, app.slider,
                                           app.checkbox, app.toggle) for widget in kind}
        self.assertFalse(hidden & rendered)
        self.assertTrue(any("Découpage reconnu automatiquement" in item.value for item in app.caption))
        request = self.launch(app).request
        self.assertEqual(request.architecture.t7_leader, "GGG")
        self.assertEqual(request.architecture.rbs_loop, "AGAGGAGA")
        self.assertEqual(request.architecture.toehold_length, 12)
        self.assertEqual(request.architecture.stem_length, 9)
        self.assertEqual(len(request.architecture.resolved_frame_linker()), 21)
        self.assertEqual(request.fixed_switch_sequence, SWITCH)

    def test_finished_reading_result_saves_and_reopens_with_full_reporter_actions(self):
        app = self.valid_app()
        job = self.launch(app)
        result = design_result_fixture(job.request)
        result.candidates[0].switch_seq = job.request.fixed_switch_sequence
        result.optimization = {
            "algorithm": "fixed_sequence_evaluation", "objective": "supplied_sequences",
            "sequence_optimization": False, "evaluations": 1, "trials": [],
        }
        result.output_dir.mkdir(parents=True)
        result_file = result.output_dir / (result.run_id + "_results.json")
        result_file.write_text(json.dumps(result.to_dict()), encoding="utf-8")
        job._state = "done"
        job._result = result
        job._exports = {"json": str(result_file)}
        app.run()
        self.assertFalse(app.exception)
        self.assertIs(app.session_state["last_result"], result)
        self.button(app, "Enregistrer ce toehold switch").click().run()
        self.assertFalse(app.exception)
        self.assertFalse(app.error)
        saved = list_saved_switches(saved_switch_root(self.project.path))[0]
        self.assertEqual(saved.snapshot["sequences"]["switch"], SWITCH)
        self.assertEqual(saved.snapshot["sequences"]["trigger"], TRIGGER)
        self.assertEqual(saved.snapshot["sequences"]["aptamer"], APTAMER)
        self.assertEqual(saved.snapshot["request"]["fixed_switch_sequence"], SWITCH)
        self.assertEqual(saved.snapshot["candidate"]["on_yield_pct"], 60.)
        context = build_saved_switch_context(saved.snapshot, {})
        self.assertTrue(context["switch_sequence"].endswith(REPORTERS["sfgfp"].full_sequence_after_start_rna))
        self.button(app, "Ouvrir les toehold switches enregistrés").click().run()
        self.assertFalse(app.exception)
        self.assertFalse(app.error)
        self.assertEqual(app.radio(key="nav").value, "Toehold switches enregistrés")
        self.assertIn("Contexte du toehold switch avec rapporteur complet",
                      app.radio(key="saved_switch_section_selector").options)
        self.assertIn("Simulation cinétique et film", app.radio(key="saved_switch_section_selector").options)
        self.native.assert_not_called()

    def test_english_reading_page_translates_mode_and_main_controls(self):
        app = self.valid_app(language="en")
        self.assertEqual(app.radio(key="design_mode").options, ["Design", "Sequence reading"])
        self.assertEqual(app.text_area(key="reading_switch_input").label, "Switch sequence")
        self.assertFalse(self.button(app, "Compute metrics").disabled)
        self.assertTrue(any(item.value == "One analysis of the supplied sequences, without optimization."
                            for item in app.caption))
        self.native.assert_not_called()


class ReadingContinuationUITests(unittest.TestCase):
    def test_reading_progress_uses_evaluations_with_resume_or_pause_and_no_finish(self):
        from aptaswitch_core.design_checkpoint import CHECKPOINT_NAME

        script = """
import streamlit as st
from aptaswitch_core.localization import set_language
from aptaswitch_studio.design_run_ui import render_design_continuation, render_design_run_controls
set_language('fr')
job = st.session_state.get('design_job')
if job is not None and job.snapshot().state == 'running':
    render_design_run_controls(job, job.snapshot())
render_design_continuation(st.session_state['result'])
"""
        with tempfile.TemporaryDirectory() as folder:
            output = Path(folder)
            (output / CHECKPOINT_NAME).touch()
            request = SimpleNamespace(fixed_switch_sequence=SWITCH, trials=1)
            for state in ("paused", "running"):
                for evaluated in (0, 1):
                    with self.subTest(state=state, evaluated=evaluated):
                        result = SimpleNamespace(
                            request=request, output_dir=output, status=state,
                            optimization={"evaluations": evaluated,
                                          "checkpoint": {"phase": "design", "completed_trials": 19},
                                          "trials": [None] * 19},
                        )
                        snapshot = JobSnapshot(state, evaluated, 1, "", can_request_control=True)
                        job = SimpleNamespace(request=request, supports_continuation=True,
                                              snapshot=lambda: snapshot, pause=Mock())
                        app = AppTest.from_string(script, default_timeout=15)
                        app.session_state["ui_language"] = "fr"
                        app.session_state["result"] = result
                        app.session_state["design_job"] = job
                        app.run()
                        self.assertFalse(app.exception)
                        self.assertFalse(app.error)
                        text = [item.value for item in (app.markdown if state == "paused" else app.caption)]
                        self.assertTrue(any(f"Lecture de séquences · {evaluated}/1 analyse terminée" in value
                                            for value in text))
                        self.assertFalse(any("Terminer" in button.label for button in app.button))
                        if state == "paused":
                            self.assertEqual([button.key for button in app.button], ["resume_design_run"])
                            self.assertFalse(app.button(key="resume_design_run").disabled)
                        else:
                            self.assertEqual([button.key for button in app.button], ["pause_design_run"])
                            app.button(key="pause_design_run").click().run()
                            self.assertFalse(app.exception)
                            job.pause.assert_called_once()


if __name__ == "__main__":
    unittest.main()
