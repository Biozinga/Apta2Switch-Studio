"""The standard activator and opt-in customization must reach the run unchanged."""

from __future__ import annotations

import tempfile
import unittest
from contextlib import ExitStack
from html.parser import HTMLParser
from pathlib import Path
from unittest.mock import patch

from tests.ui_helpers import FrenchAppTest as AppTest

from aptaswitch_core.architecture import architecture_layout


APP_PATH = Path(__file__).resolve().parents[1] / "src" / "aptaswitch_studio" / "web_app.py"
TRIGGER = "GTCCAGGCTGGTATAATTAGATCCACGTAC"


class ArchitectureOptions(HTMLParser):
    def __init__(self):
        super().__init__()
        self.options = []
        self.current = None

    def handle_starttag(self, tag, attrs):
        if tag == "option":
            self.current = [dict(attrs), ""]
            self.options.append(self.current)

    def handle_data(self, data):
        if self.current is not None:
            self.current[1] += data

    def handle_endtag(self, tag):
        if tag == "option":
            self.current = None


class ActivatorUITests(unittest.TestCase):
    def setUp(self):
        self.context = ExitStack()
        self.addCleanup(self.context.close)
        self.output = self.context.enter_context(tempfile.TemporaryDirectory())
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
                side_effect=AssertionError("UI tests must not modify user settings"),
            ))
        self.start = self.context.enter_context(patch(
            "aptaswitch_studio.web_runtime.DesignRunJob.start", autospec=True,
        ))

    def app(self, trigger=TRIGGER, **existing_state):
        app = AppTest.from_file(str(APP_PATH), default_timeout=20)
        app.session_state["dependency_detection_done"] = True
        app.session_state["output_dir"] = self.output
        app.session_state["aptamer_input"] = "ACGTACGT"
        app.session_state["trigger_input"] = trigger
        for key, value in existing_state.items():
            app.session_state[key] = value
        app.run()
        self.assertFalse(app.exception)
        return app

    def launch(self, app):
        count = self.start.call_count
        next(button for button in app.button if button.label == "Lancer la conception").click().run()
        self.assertFalse(app.exception)
        self.assertEqual(self.start.call_count, count + 1)
        return self.start.call_args.args[0].request

    def resume_design(self, app):
        app.session_state["design_job"] = None
        app.radio(key="nav").set_value("Conception").run()
        self.assertFalse(app.exception)

    def assert_standard(self, request, trigger=TRIGGER):
        architecture = request.architecture
        expected_split = {24: (15, 9), 30: (12, 18), 32: (12, 20)}[len(trigger)]
        self.assertEqual((architecture.toehold_length, architecture.stem_length), expected_split)
        self.assertEqual(architecture.rbs_prefix_length, 3)
        self.assertEqual(architecture.upper_stem_length, 6)
        self.assertEqual(architecture.upper_stem2_length, 0 if len(trigger) == 30 else 4)
        self.assertEqual(architecture.rbs_to_aug_distance, 6)
        self.assertEqual(architecture.aug_spacer, 0)
        self.assertIsNone(architecture.frame_linker)
        self.assertEqual(architecture.resolved_frame_linker(), {24: "NN", 30: "", 32: ""}[len(trigger)])
        if len(trigger) == 30:
            layout = architecture_layout(architecture)
            # AUG is codon 1, followed by six codons (18 nt), then the reporter.
            self.assertEqual(layout["reporter_start_codon"], 8)
        self.assertEqual(request.sequences.binding_start, 0)
        self.assertEqual(request.sequences.binding_length, len(trigger))
        self.assertEqual(request.sequences.effective_trigger(), trigger)

    def test_standard_form_has_only_supported_architecture_and_editable_sequences(self):
        app = self.app()
        self.assertFalse(app.checkbox(key="customize_tswitch").value)
        for collection, keys in (
            (app.number_input, {"binding_start", "binding_end", "rbs_aug_distance", "aug_spacer"}),
            (app.slider, {"toehold_length"}),
            (app.text_input, {"frame_linker"}),
            (app.toggle, {"linker_auto"}),
        ):
            self.assertTrue(keys.isdisjoint({item.key for item in collection}))
        self.assertEqual(app.text_input(key="t7_leader").label, "Leader 5′")
        self.assertEqual(app.text_input(key="rbs_loop").label, "RBS")
        self.assertEqual(app.text_input(key="rbs_loop").value, "AGAGGAGA")
        self.assertFalse(app.text_input(key="t7_leader").disabled)
        self.assertFalse(app.text_input(key="rbs_loop").disabled)

        selector = next(item.proto.body for item in app.get("html")
                        if 'id="switch-architecture-choice"' in item.proto.body)
        options = ArchitectureOptions()
        options.feed(selector)
        self.assertEqual([label for _, label in options.options], [
            "Toehold activateur", "Toehold répresseur", "Toehold répresseur 3WJ",
        ])
        self.assertIn("selected", options.options[0][0])
        self.assertNotIn("disabled", options.options[0][0])
        self.assertTrue(all("disabled" in attrs for attrs, _ in options.options[1:]))

        app.text_input(key="t7_leader").set_value("AAGGG")
        app.text_input(key="rbs_loop").set_value("AGGAGG").run()
        request = self.launch(app)
        self.assert_standard(request)
        self.assertEqual(request.architecture.t7_leader, "AAGGG")
        self.assertEqual(request.architecture.rbs_loop, "AGGAGG")
        self.assertEqual(request.sequences.trigger_dna, TRIGGER)

        # Updating a previously configured trigger must update the standard
        # binding region and preview/run split without truncating any bases.
        for trigger in (TRIGGER[:24], TRIGGER + "GT"):
            with self.subTest(trigger_length=len(trigger)):
                self.resume_design(app)
                app.text_area(key="trigger_input").set_value(trigger).run()
                self.assert_standard(self.launch(app), trigger)

    def test_rbs_default_migrates_once_and_preserves_custom_sequences(self):
        app = self.app(rbs_loop="UAGAGGAGAAC")
        self.assertEqual(app.text_input(key="rbs_loop").value, "AGAGGAGA")
        request = self.launch(app)
        self.assert_standard(request)
        self.assertEqual(request.architecture.rbs_loop, "AGAGGAGA")

        # A deliberately entered sequence remains editable even when it happens
        # to match the old default. Migration only affects a preexisting session.
        self.resume_design(app)
        app.text_input(key="rbs_loop").set_value("UAGAGGAGAAC").run()
        self.assertEqual(self.launch(app).architecture.rbs_loop, "UAGAGGAGAAC")

        custom_app = self.app(rbs_loop="AGGAGG")
        self.assertEqual(custom_app.text_input(key="rbs_loop").value, "AGGAGG")
        self.assertEqual(self.launch(custom_app).architecture.rbs_loop, "AGGAGG")

    def test_empty_rbs_uses_new_default_in_run(self):
        app = self.app()
        app.text_input(key="rbs_loop").set_value("").run()
        request = self.launch(app)
        self.assert_standard(request)
        self.assertEqual(request.architecture.rbs_loop, "AGAGGAGA")

    def test_short_trigger_requires_extension_or_explicit_customization(self):
        app = self.app(trigger="ACGTACGTACGT")
        button = next(button for button in app.button if button.label == "Lancer la conception")
        self.assertTrue(button.disabled)
        self.assertFalse(app.exception)
        self.start.assert_not_called()
        self.assertIsNone(app.session_state["design_job"])
        errors = " ".join(item.value for collection in (app.error, app.info, app.warning) for item in collection)
        self.assertIn("24 nt", errors)
        self.assertIn("Étendez", errors)
        self.assertIn("personnalisez le tSwitch", errors)

        app.checkbox(key="customize_tswitch").check().run()
        request = self.launch(app)
        self.assertEqual(request.architecture.trigger_region_length, 12)
        self.assertEqual(request.sequences.effective_trigger(), "ACGTACGTACGT")

    def test_custom_values_reach_run_and_survive_standard_toggle_and_navigation(self):
        app = self.app()
        app.checkbox(key="customize_tswitch").check().run()
        self.assertTrue(any("déconseillée" in warning.value for warning in app.warning))
        app.number_input(key="binding_start").set_value(4)
        app.number_input(key="binding_end").set_value(25).run()
        app.slider(key="toehold_length").set_value(10)
        app.number_input(key="rbs_aug_distance").set_value(8)
        app.toggle(key="linker_auto").set_value(False).run()
        app.text_input(key="frame_linker").set_value("CN").run()
        custom = self.launch(app)
        self.assertEqual(custom.sequences.binding_start, 3)
        self.assertEqual(custom.sequences.binding_length, 22)
        self.assertEqual(custom.sequences.effective_trigger(), TRIGGER[3:25])
        self.assertEqual((custom.architecture.toehold_length, custom.architecture.stem_length), (10, 12))
        self.assertEqual(custom.architecture.rbs_to_aug_distance, 8)
        self.assertEqual(custom.architecture.rbs_prefix_length, 3)
        self.assertEqual(custom.architecture.frame_linker, "CN")

        self.resume_design(app)
        app.checkbox(key="customize_tswitch").uncheck().run()
        self.assert_standard(self.launch(app))
        self.resume_design(app)
        app.radio(key="nav").set_value("Réglages").run()
        app.radio(key="nav").set_value("Conception").run()
        app.checkbox(key="customize_tswitch").check().run()
        self.assertFalse(app.exception)
        self.assertEqual(app.number_input(key="binding_start").value, 4)
        self.assertEqual(app.number_input(key="binding_end").value, 25)
        self.assertEqual(app.slider(key="toehold_length").value, 10)
        self.assertEqual(app.number_input(key="rbs_aug_distance").value, 8)
        self.assertFalse(app.toggle(key="linker_auto").value)
        self.assertEqual(app.text_input(key="frame_linker").value, "CN")
        self.assertEqual(self.launch(app).architecture, custom.architecture)

        self.resume_design(app)
        extended_trigger = "ACGT" + TRIGGER
        app.text_area(key="trigger_input").set_value(extended_trigger).run()
        self.assertFalse(app.exception)
        self.assertEqual(app.number_input(key="binding_start").value, 4)
        self.assertEqual(app.number_input(key="binding_end").value, 25)
        self.assertEqual(app.slider(key="toehold_length").value, 10)
        self.assertEqual(app.number_input(key="rbs_aug_distance").value, 8)
        extended = self.launch(app)
        self.assertEqual(extended.architecture, custom.architecture)
        self.assertEqual(extended.sequences.binding_start, 3)
        self.assertEqual(extended.sequences.binding_length, 22)
        self.assertEqual(extended.sequences.effective_trigger(), extended_trigger[3:25])


if __name__ == "__main__":
    unittest.main()
