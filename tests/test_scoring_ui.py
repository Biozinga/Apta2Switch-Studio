"""Scoring controls must describe and submit the same physical calculation."""

from __future__ import annotations

import tempfile
import unittest
from contextlib import ExitStack
from pathlib import Path
from unittest.mock import patch

from tests.ui_helpers import FrenchAppTest as AppTest


APP_PATH = Path(__file__).resolve().parents[1] / "src" / "aptaswitch_studio" / "web_app.py"
DEFAULT_WEIGHTS = {
    "weight_rbs": 6.0,
    "weight_defect": 6.0,
    "weight_on": 3.0,
    "weight_ddg": 0.25,
    "weight_leak": 0.05,
}
PREVIOUS_WEIGHTS = {
    "weight_rbs": 6.0,
    "weight_defect": 6.0,
    "weight_leak": 3.0,
    "weight_on": 1.0,
    "weight_ddg": 0.25,
}
OLD_WEIGHTS = {
    "weight_on": 1.0,
    "weight_leak": 3.0,
    "weight_defect": 3.0,
    "weight_ddg": 0.5,
    "weight_rbs": 4.0,
}


class ScoringUITests(unittest.TestCase):
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

    def app(self, **existing_state):
        app = AppTest.from_file(str(APP_PATH), default_timeout=20)
        app.session_state["dependency_detection_done"] = True
        app.session_state["output_dir"] = self.output
        app.session_state["aptamer_input"] = "ACGTACGT"
        app.session_state["trigger_input"] = "GTCCAGGCTGGTATAATTAGATCCACGTAC"
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

    def assert_weights(self, app, expected):
        self.assertEqual(
            {key: app.number_input(key=key).value for key in expected}, expected,
        )

    def assert_request_weights(self, request, expected):
        self.assertEqual({
            "weight_rbs": request.scoring.rbs,
            "weight_defect": request.scoring.defect,
            "weight_leak": request.scoring.leak,
            "weight_on": request.scoring.on_yield,
            "weight_ddg": request.scoring.ddg,
        }, expected)

    def test_default_priority_order_and_weights_reach_actual_run_request(self):
        app = self.app()
        controls = [item for item in app.number_input if item.key in DEFAULT_WEIGHTS]
        self.assertEqual([item.key for item in controls], list(DEFAULT_WEIGHTS))
        self.assertTrue(all(item.proto.help for item in controls))
        self.assert_weights(app, DEFAULT_WEIGHTS)
        self.assert_request_weights(self.launch(app), DEFAULT_WEIGHTS)

    def test_migration_is_once_and_preserves_custom_weights_across_navigation(self):
        app = self.app(**OLD_WEIGHTS)
        self.assert_weights(app, DEFAULT_WEIGHTS)
        self.assertEqual(app.session_state["scoring_defaults_revision"], 2)

        # Re-entering the old numbers deliberately must not trigger migration
        # again, including when Streamlit removes the widgets on another page.
        for key, value in OLD_WEIGHTS.items():
            app.number_input(key=key).set_value(value)
        app.run()
        app.radio(key="nav").set_value("Réglages").run()
        app.radio(key="nav").set_value("Conception").run()
        self.assertFalse(app.exception)
        self.assert_weights(app, OLD_WEIGHTS)
        self.assert_request_weights(self.launch(app), OLD_WEIGHTS)

        # Even a single deliberate change means a preexisting set is custom.
        custom = {**OLD_WEIGHTS, "weight_on": 2.0}
        custom_app = self.app(**custom)
        self.assert_weights(custom_app, custom)
        self.assert_request_weights(self.launch(custom_app), custom)

    def test_previous_default_migrates_but_customized_revision_one_is_preserved(self):
        for revision in (0, 1):
            with self.subTest(previous_revision=revision):
                app = self.app(scoring_defaults_revision=revision, **PREVIOUS_WEIGHTS)
                self.assert_weights(app, DEFAULT_WEIGHTS)
                self.assertEqual(app.session_state["scoring_defaults_revision"], 2)
                self.assert_request_weights(self.launch(app), DEFAULT_WEIGHTS)

        custom = {**PREVIOUS_WEIGHTS, "weight_leak": 0.0}
        custom_app = self.app(scoring_defaults_revision=1, **custom)
        self.assert_weights(custom_app, custom)
        self.assertEqual(custom_app.session_state["scoring_defaults_revision"], 2)
        self.assert_request_weights(self.launch(custom_app), custom)

        # Restoring earlier defaults after migration is a deliberate user choice.
        custom_app.radio(key="nav").set_value("Conception").run()
        for key, value in PREVIOUS_WEIGHTS.items():
            custom_app.number_input(key=key).set_value(value)
        custom_app.run()
        custom_app.radio(key="nav").set_value("Réglages").run()
        custom_app.radio(key="nav").set_value("Conception").run()
        self.assertFalse(custom_app.exception)
        self.assert_weights(custom_app, PREVIOUS_WEIGHTS)

        # At revision 1, the original vector can only be a deliberate reset.
        original_app = self.app(scoring_defaults_revision=1, **OLD_WEIGHTS)
        self.assert_weights(original_app, OLD_WEIGHTS)
        self.assert_request_weights(self.launch(original_app), OLD_WEIGHTS)

    def test_help_describes_boundaries_metric_directions_and_model_limits(self):
        app = self.app()
        help_text = {key: app.number_input(key=key).proto.help for key in DEFAULT_WEIGHTS}
        rbs = help_text["weight_rbs"]
        self.assertIn("première base de la boucle", rbs)
        self.assertIn("base précédant le premier nucléotide du rapporteur", rbs)
        self.assertIn("avant le RBS (NNN dans la cible) sont incluses", rbs)
        self.assertIn("tige située avant la boucle et le rapporteur sont exclus", rbs)
        self.assertIn("moins négatif, proche de 0", rbs)

        defect = help_text["weight_defect"]
        self.assertIn("cibles OFF et ON", defect)
        self.assertIn("structures et aux concentrations cibles", defect)
        self.assertIn("Plus faible est meilleur", defect)
        self.assertIn("0 est idéal", defect)

        leak = help_text["weight_leak"]
        self.assertIn("Plus faible est meilleur", leak)
        self.assertIn("aptamère:trigger:switch de 1:1:1", leak)
        self.assertIn("concentration d’aptamère est actuellement égale à celle du trigger", leak)
        self.assertIn("[complexe trigger–switch] / [switch] initial", leak)
        self.assertIn("pas une mesure de traduction", leak)
        self.assertIn("ligand n’est pas simulé", leak)
        self.assertIn("aptamères sont dans un autre tube", leak)
        self.assertIn("scénario de comparaison", leak)
        self.assertIn("pas une simulation de la fuite ou du transfert entre vos tubes", leak)
        self.assertIn("poids de 0,05", leak)
        self.assertIn("candidats très proches", leak)

        on = help_text["weight_on"]
        self.assertIn("Plus élevé est meilleur", on)
        self.assertIn("trigger:switch de 1:1", on)
        self.assertIn("sans aptamère ni ligand", on)
        self.assertIn("[complexe trigger–switch] / [switch] initial", on)
        self.assertIn("pas directement l’activation de la traduction", on)
        self.assertIn("prioritaire sur la fuite", on)
        self.assertIn("60 fois plus", on)

        ddg = help_text["weight_ddg"]
        self.assertIn("− ΔG MFE du switch seul − ΔG MFE du trigger seul", ddg)
        self.assertIn("plus négative", ddg)
        self.assertIn("ni une barrière d’activation ni une vitesse de réaction", ddg)

    def test_tube_help_tracks_physical_concentrations_and_ratio_not_display_units(self):
        app = self.app()
        self.assertIn("trigger:switch de 1:1", app.number_input(key="weight_on").proto.help)

        app.number_input(key="trigger_concentration_value").set_value(10.0).run()
        self.assertFalse(app.exception)
        on = app.number_input(key="weight_on").proto.help
        leak = app.number_input(key="weight_leak").proto.help
        self.assertIn("trigger libre (10 µM) et le switch (5 µM)", on)
        self.assertIn("trigger:switch de 2:1", on)
        self.assertIn("aptamère (10 µM), trigger (10 µM) et switch (5 µM)", leak)
        self.assertIn("aptamère:trigger:switch de 2:2:1", leak)

        # 10 µM and 10000 nM describe the same tube; changing the display unit
        # must not falsely change the ratio or the wording of the experiment.
        app.selectbox(key="trigger_concentration_unit").set_value("nM").run()
        self.assertFalse(app.exception)
        self.assertAlmostEqual(app.number_input(key="trigger_concentration_value").value, 10000.0)
        self.assertEqual(app.number_input(key="weight_on").proto.help, on)
        self.assertEqual(app.number_input(key="weight_leak").proto.help, leak)
        request = self.launch(app)
        self.assertAlmostEqual(request.thermo.resolved_trigger_concentration_m, 1e-5, delta=1e-15)
        self.assertAlmostEqual(request.thermo.resolved_switch_concentration_m, 5e-6, delta=1e-15)


if __name__ == "__main__":
    unittest.main()
