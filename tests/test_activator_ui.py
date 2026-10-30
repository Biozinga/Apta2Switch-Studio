"""Standard geometry and editable leader/RBS/linker reach the run unchanged."""

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
CLASSIC_LINKER = "AACCUGGCGGCAGCGCAAAAG"


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
        from tests.native_ui_fixtures import patch_native_runtime
        patch_native_runtime(self)
        self.context = ExitStack()
        self.addCleanup(self.context.close)
        self.output = self.context.enter_context(tempfile.TemporaryDirectory())
        for target, value in (
            ("aptaswitch_core.settings.load_settings", {}),
            ("aptaswitch_core.thermodynamics.setup.configured_library_path", "/test/native"),
            ("aptaswitch_studio.extension_ui.configured_library_path", None),
        ):
            self.context.enter_context(patch(target, return_value=value))
        for target in ("save_settings",):
            self.context.enter_context(patch(
                f"aptaswitch_core.settings.{target}",
                side_effect=AssertionError("UI tests must not modify user settings"),
            ))
        self.start = self.context.enter_context(patch(
            "aptaswitch_studio.web_runtime.DesignRunJob.start", autospec=True,
        ))

    def app(self, trigger=TRIGGER, **existing_state):
        app = AppTest.from_file(str(APP_PATH), default_timeout=20)
        app.session_state["dependency_detection_done"] = True
        self.output = app.session_state["output_dir"]
        app.session_state["aptamer_input"] = "ACGTACGT"
        app.session_state["trigger_input"] = trigger
        for key, value in existing_state.items():
            app.session_state[key] = value
        app.run()
        self.assertEqual(app.radio(key="design_section").value, "system")
        app.radio(key="design_section").set_value("switch").run()
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

    def assert_standard(self, request, trigger=TRIGGER, manual_linker=None):
        architecture = request.architecture
        expected_split = {24: (15, 9), 30: (12, 9), 32: (12, 20)}[len(trigger)]
        self.assertEqual((architecture.toehold_length, architecture.stem_length), expected_split)
        self.assertEqual(architecture.rbs_prefix_length, 3)
        self.assertEqual(architecture.upper_stem_length, 6)
        self.assertEqual(architecture.upper_stem2_length, 0 if len(trigger) == 30 else 4)
        self.assertEqual(architecture.rbs_to_aug_distance, 6)
        self.assertEqual(architecture.aug_spacer, 0)
        self.assertEqual(architecture.binding_layout, "full_left_arm" if len(trigger) == 30 else "legacy_lower_stem")
        expected_linker = (manual_linker if manual_linker is not None
                           else CLASSIC_LINKER if len(trigger) == 30 else None)
        resolved_linker = (manual_linker if manual_linker is not None
                           else {24: "XX", 30: CLASSIC_LINKER, 32: ""}[len(trigger)])
        self.assertEqual(architecture.frame_linker, expected_linker)
        self.assertEqual(architecture.resolved_frame_linker(), resolved_linker)
        if len(trigger) == 30:
            layout = architecture_layout(architecture)
            # AUG is codon 1, then 9 stem bases and the selected linker.
            self.assertEqual(layout["reporter_start"] - layout["aug_start"], 12 + len(resolved_linker))
            self.assertEqual(layout["reporter_start"] - layout["aug_start"] - 3, 9 + len(resolved_linker))
        self.assertEqual(request.sequences.binding_start, 0)
        self.assertIn(request.sequences.binding_length, (None, len(trigger)))
        self.assertEqual(request.sequences.effective_trigger(), trigger)

    def test_stale_custom_geometry_and_binding_never_change_standard_request(self):
        for trigger in (TRIGGER[:24], TRIGGER, TRIGGER + "GT"):
            with self.subTest(trigger_length=len(trigger)):
                app = self.app(trigger=trigger, trigger_snapshot=trigger,
                               customize_tswitch=True, toehold_length=1,
                               rbs_aug_distance=25, binding_start=4, binding_end=12)
                self.assertNotIn("customize_tswitch", {item.key for item in app.checkbox})
                self.assert_standard(self.launch(app), trigger)

    def test_fixed_chemistry_forces_zero_magnesium_in_design_and_request(self):
        app = self.app(
            magnesium_concentration_value=10.0,
            magnesium_concentration_unit="mM",
            extension_magnesium=7.0,
            extension_magnesium_unit="mM",
        )
        value = app.number_input(key="magnesium_concentration_value")
        self.assertEqual(value.value, 0.0)
        self.assertTrue(value.disabled)
        self.assertTrue(app.selectbox(key="magnesium_concentration_unit").disabled)
        app.radio(key="design_section").set_value("system").run()
        self.assertEqual(app.number_input(key="extension_magnesium").value, 0.0)
        self.assertTrue(app.number_input(key="extension_magnesium").disabled)
        app.radio(key="design_section").set_value("switch").run()
        self.assertEqual(self.launch(app).thermo.magnesium_m, 0.0)

    def test_switch_is_rna_trigger_is_dna_without_a_chemistry_choice(self):
        app = self.app()
        self.assertNotIn("material_label", app.session_state)
        self.assertNotIn("extension_material", app.session_state)
        self.assertFalse(any(widget.label == "Chimie des brins" for widget in app.selectbox))
        sodium = app.number_input(key="sodium_concentration_value")
        self.assertEqual((sodium.min, sodium.max), (120., 1000.))
        self.assertEqual(app.number_input(key="temperature_c").min, 25.)
        self.assertEqual(app.number_input(key="temperature_c").max, 55.)
        self.assertTrue(app.number_input(key="magnesium_concentration_value").disabled)
        request = self.launch(app)
        self.assertEqual((request.thermo.material, request.thermo.switch_material, request.thermo.trigger_material),
                         ("experimental-Sugimoto95-Vienna04-loopmean", "rna", "dna"))
        self.assertEqual(request.thermo.magnesium_m, 0.)

    def test_removed_spacing_never_pads_default_linker(self):
        for spacing in (6, 7, 9):
            with self.subTest(spacing=spacing):
                app = self.app(customize_tswitch=True, rbs_aug_distance=spacing)
                self.assertEqual(app.text_input(key="frame_linker").value, CLASSIC_LINKER)
                app.text_input(key="frame_linker").set_value("CCX").run()
                app.button(key="restore_architecture_defaults").click().run()
                self.assertTrue(app.session_state["linker_auto"])
                self.assertEqual(app.text_input(key="frame_linker").value, CLASSIC_LINKER)
                request = self.launch(app)
                self.assert_standard(request)
                layout = architecture_layout(request.architecture)
                self.assertEqual((layout["reporter_start"] - layout["aug_start"]) % 3, 0)

    def test_optimized_linker_keeps_length_and_reaches_design_as_free_positions(self):
        app = self.app()
        for pattern, length in ((CLASSIC_LINKER, 21), (" gct aac ", 6)):
            with self.subTest(pattern=pattern):
                app.text_input(key="frame_linker").set_value(pattern).run()
                app.button(key="optimize_linker").click().run()
                self.assertFalse(app.exception)
                self.assertFalse(app.session_state["linker_auto"])
                self.assertEqual(app.text_input(key="frame_linker").value, "X" * length)
                self.assertTrue(app.button(key="optimize_linker").disabled)
                self.assert_standard(self.launch(app), manual_linker="X" * length)
                self.resume_design(app)
        app.text_input(key="frame_linker").set_value("").run()
        self.assertTrue(app.button(key="optimize_linker").disabled)

    def test_restore_defaults_resets_editable_sequences_only(self):
        app = self.app(temperature_c=31., weight_defect=3.,
                       customize_tswitch=True, toehold_length=10,
                       rbs_aug_distance=8, binding_start=4, binding_end=25)
        app.text_input(key="t7_leader").set_value("AAA")
        app.text_input(key="rbs_loop").set_value("AGGAGG")
        app.text_input(key="frame_linker").set_value("XXX").run()
        app.button(key="restore_architecture_defaults").click().run()
        self.assertFalse(app.exception)
        self.assertEqual(app.text_input(key="t7_leader").value, "GGG")
        self.assertEqual(app.text_input(key="rbs_loop").value, "AGAGGAGA")
        self.assertEqual(app.text_input(key="frame_linker").value, CLASSIC_LINKER)
        self.assertEqual(app.session_state["temperature_c"], 31.)
        self.assertEqual(app.session_state["weight_defect"], 3.)
        self.assertEqual(app.session_state["aptamer_input"], "ACGTACGT")
        self.assertEqual(app.session_state["trigger_input"], TRIGGER)
        self.assert_standard(self.launch(app))

    def test_manual_linker_survives_navigation_and_can_restore_default(self):
        for manual, normalized in (("CCX", "CCX"), ("", ""), (" ac txxx ", "ACUXXX")):
            with self.subTest(manual=manual):
                app = self.app()
                app.text_input(key="frame_linker").set_value(manual).run()
                self.assertFalse(app.session_state["linker_auto"])
                request = self.launch(app)
                self.assert_standard(request, manual_linker=normalized)
                self.assertEqual(architecture_layout(request.architecture)["reporter_start"], 62 + len(normalized))
                self.resume_design(app)
                app.radio(key="nav").set_value("Réglages").run()
                app.radio(key="nav").set_value("Conception").run()
                self.assertFalse(app.session_state["linker_auto"])
                self.assertEqual(app.text_input(key="frame_linker").value, manual)
                self.assert_standard(self.launch(app), manual_linker=normalized)
                self.resume_design(app)
                app.button(key="restore_architecture_defaults").click().run()
                self.assertTrue(app.session_state["linker_auto"])
                self.assertEqual(app.text_input(key="frame_linker").value, CLASSIC_LINKER)
                self.assert_standard(self.launch(app))

    def test_obsolete_n_linker_is_not_silently_removed_or_used_for_design(self):
        app = self.app()
        app.text_input(key="frame_linker").set_value("CCN").run()
        self.assertFalse(app.exception)
        self.assertEqual(app.text_input(key="frame_linker").value, "CCN")
        self.assertTrue(next(button for button in app.button if button.label == "Lancer la conception").disabled)
        self.assertTrue(any("Le linker" in warning.value for warning in app.warning))
        self.start.assert_not_called()
        app.text_input(key="frame_linker").set_value("CCX").run()
        self.assertFalse(app.exception)
        self.assertEqual(self.launch(app).architecture.frame_linker, "CCX")

    def test_standard_form_has_only_supported_architecture_and_editable_sequences(self):
        app = self.app()
        self.assertNotIn("customize_tswitch", {item.key for item in app.checkbox})
        for collection, keys in (
            (app.number_input, {"binding_start", "binding_end", "rbs_aug_distance", "aug_spacer"}),
            (app.slider, {"toehold_length"}),
            (app.toggle, {"linker_auto"}),
        ):
            self.assertTrue(keys.isdisjoint({item.key for item in collection}))
        self.assertEqual(app.text_input(key="t7_leader").label, "Leader 5′")
        self.assertEqual(app.text_input(key="rbs_loop").label, "RBS")
        self.assertEqual(app.text_input(key="rbs_loop").value, "AGAGGAGA")
        self.assertFalse(app.text_input(key="t7_leader").disabled)
        self.assertFalse(app.text_input(key="rbs_loop").disabled)
        self.assertEqual(app.text_input(key="frame_linker").label, "Linker")
        self.assertEqual(app.text_input(key="frame_linker").value, CLASSIC_LINKER)
        self.assertFalse(app.text_input(key="frame_linker").disabled)
        self.assertTrue(app.button(key="restore_architecture_defaults").disabled)
        self.assertTrue(any("Linker : 21 nt" in item.value for item in app.caption))
        self.assertTrue(any("après les 9 nt de tige" in item.value for item in app.caption))

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
                app.radio(key="design_section").set_value("system").run()
                app.text_area(key="trigger_input").set_value(trigger).run()
                app.radio(key="design_section").set_value("switch").run()
                self.assertEqual(app.text_input(key="frame_linker").value, {24: "XX", 32: ""}[len(trigger)])
                self.assertFalse(app.text_input(key="frame_linker").disabled)
                self.assertTrue(app.session_state["linker_auto"])
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

    def test_short_trigger_requires_extension_even_with_stale_customization(self):
        app = self.app(trigger="ACGTACGTACGT", customize_tswitch=True,
                       toehold_length=8, binding_start=1, binding_end=12)
        button = next(button for button in app.button if button.label == "Lancer la conception")
        self.assertTrue(button.disabled)
        self.assertFalse(app.exception)
        self.start.assert_not_called()
        self.assertIsNone(app.session_state["design_job"])
        errors = " ".join(item.value for collection in (app.error, app.info, app.warning) for item in collection)
        self.assertIn("24 nt", errors)
        self.assertIn("Activez l’extension du trigger", errors)
        self.assertNotIn("personnalisez le toehold switch", errors)
        self.assertNotIn("customize_tswitch", {item.key for item in app.checkbox})

    def test_manual_linker_is_never_padded_and_editable_sequences_survive_trigger_change(self):
        app = self.app(customize_tswitch=True, toehold_length=10,
                       rbs_aug_distance=8, binding_start=4, binding_end=25)
        app.text_input(key="t7_leader").set_value("AAA")
        app.text_input(key="rbs_loop").set_value("AGGAGG")
        app.text_input(key="frame_linker").set_value("CX").run()
        # An old customized architecture accepted this length, but the standard
        # one does not. Keep the explicit sequence and report the frame error.
        self.assertEqual(app.text_input(key="frame_linker").value, "CX")
        self.assertTrue(next(button for button in app.button if button.label == "Lancer la conception").disabled)
        self.assertTrue(any("cadre de lecture" in item.value for item in app.warning))
        self.start.assert_not_called()
        app.text_input(key="frame_linker").set_value("CCX").run()
        request = self.launch(app)
        self.assert_standard(request, manual_linker="CCX")
        self.assertEqual(request.architecture.t7_leader, "AAA")
        self.assertEqual(request.architecture.rbs_loop, "AGGAGG")
        self.resume_design(app)
        app.radio(key="nav").set_value("Réglages").run()
        app.radio(key="nav").set_value("Conception").run()
        self.assertFalse(app.exception)
        self.assertFalse(app.session_state["linker_auto"])
        self.assertEqual(app.text_input(key="frame_linker").value, "CCX")
        self.assertEqual(self.launch(app).architecture, request.architecture)

        self.resume_design(app)
        extended_trigger = TRIGGER + "GT"
        app.radio(key="design_section").set_value("system").run()
        app.text_area(key="trigger_input").set_value(extended_trigger).run()
        app.radio(key="design_section").set_value("switch").run()
        self.assertFalse(app.exception)
        extended = self.launch(app)
        self.assert_standard(extended, extended_trigger, manual_linker="CCX")
        self.assertEqual(extended.architecture.t7_leader, "AAA")
        self.assertEqual(extended.architecture.rbs_loop, "AGGAGG")


if __name__ == "__main__":
    unittest.main()
