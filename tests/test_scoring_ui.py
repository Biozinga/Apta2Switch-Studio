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
    "weight_leak": 0.0,
    "weight_rbs_calculator": 3.0,
}
FUNCTIONAL_KEYS = ("weight_functional_off_stem", "weight_functional_off_toehold",
                 "weight_functional_on_rbs_aug", "weight_functional_on_rbs_linker")
PREVIOUS_WEIGHTS = {
    "weight_rbs": 6.0,
    "weight_defect": 6.0,
    "weight_leak": 3.0,
    "weight_on": 1.0,
    "weight_ddg": 0.25,
    "weight_rbs_calculator": 3.0,
}
OLD_WEIGHTS = {
    "weight_on": 1.0,
    "weight_leak": 3.0,
    "weight_defect": 3.0,
    "weight_ddg": 0.5,
    "weight_rbs": 4.0,
    "weight_rbs_calculator": 3.0,
}


class ScoringUITests(unittest.TestCase):
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

    def app(self, **existing_state):
        app = AppTest.from_file(str(APP_PATH), default_timeout=20)
        app.session_state["dependency_detection_done"] = True
        self.output = app.session_state["output_dir"]
        app.session_state["aptamer_input"] = "ACGTACGT"
        app.session_state["trigger_input"] = "GTCCAGGCTGGTATAATTAGATCCACGTAC"
        for key, value in existing_state.items():
            app.session_state[key] = value
        app.run()
        self.assertFalse(app.exception)
        app.radio(key="design_section").set_value("switch").run()
        self.assertFalse(app.exception)
        return app

    def launch(self, app):
        count = self.start.call_count
        next(button for button in app.button if button.label == "Lancer la conception").click().run()
        self.assertFalse(app.exception)
        self.assertEqual(self.start.call_count, count + 1)
        return self.start.call_args.args[0].request

    def test_optimization_budget_survives_navigation_and_reaches_native_request(self):
        app = self.app()
        self.assertEqual(app.number_input(key="optimization_max_evaluations").value, 1000)
        app.number_input(key="optimization_max_evaluations").set_value(37).run()
        app.radio(key="nav").set_value("Résultats").run()
        app.radio(key="nav").set_value("Conception").run()
        self.assertEqual(app.number_input(key="optimization_max_evaluations").value, 37)
        request = self.launch(app)
        from dataclasses import asdict
        from aptaswitch_core.thermodynamics.optimization import TrialOptions
        self.assertEqual(request.optimization_options, asdict(TrialOptions(max_evaluations=37)))
        self.assertEqual(request.engine, "native-hybrid")

    def test_raising_evaluation_budget_also_allows_more_proposals(self):
        from dataclasses import asdict
        from aptaswitch_core.thermodynamics.optimization import TrialOptions
        app = self.app(optimization_max_evaluations=100000)
        request = self.launch(app)
        self.assertEqual(request.optimization_options,
                         asdict(TrialOptions(max_evaluations=100000, max_proposals=500000)))

    def test_population_lineages_and_local_tolerance_reach_the_saved_request(self):
        app = self.app()
        self.assertEqual(app.number_input(key="optimization_population_size").value, 64)
        self.assertEqual(app.number_input(key="optimization_lineages").value, 8)
        self.assertEqual(app.number_input(key="optimization_local_tolerance").value, .005)
        app.number_input(key="optimization_population_size").set_value(32)
        app.number_input(key="optimization_lineages").set_value(4)
        app.number_input(key="optimization_local_tolerance").set_value(.01).run()
        app.radio(key="nav").set_value("Résultats").run()
        app.radio(key="nav").set_value("Conception").run()
        request = self.launch(app)
        self.assertEqual(request.optimization_options["population_size"], 32)
        self.assertEqual(request.optimization_options["lineages"], 4)
        self.assertEqual(request.optimization_options["local_tolerance"], .01)

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
            "weight_rbs_calculator": request.scoring.rbs_calculator,
        }, expected)

    def test_default_priority_order_and_weights_reach_actual_run_request(self):
        app = self.app()
        controls = [item for item in app.number_input if item.key in DEFAULT_WEIGHTS]
        self.assertEqual([item.key for item in controls], list(DEFAULT_WEIGHTS))
        self.assertTrue(all(item.proto.help for item in controls))
        self.assert_weights(app, DEFAULT_WEIGHTS)
        self.assertEqual(app.number_input(key="weight_defect").label, "Défaut OFF/ON · malus")
        self.assertTrue(any("D_{\\mathrm{OFF/ON}}" in item.value for item in app.latex))
        request = self.launch(app)
        self.assert_request_weights(request, DEFAULT_WEIGHTS)
        self.assertEqual(request.defect_metric, "functional")
        self.assertEqual(request.functional_defect_weights.fractions(), {
            "off_toehold": 1/3.75, "off_stem": 1/3.75, "on_rbs_aug": 1/3.75, "on_rbs_linker": .2,
        })

    def test_rbs_calculator_weight_survives_navigation_and_reaches_request(self):
        app = self.app()
        control = app.number_input(key="weight_rbs_calculator")
        self.assertEqual(control.value, 3.)
        self.assertEqual(control.label, "RBS ON/OFF · bonus")
        self.assertIn("OSTIR 1.1.3", control.proto.help)
        self.assertIn("aucune énergie de liaison ADN–ARN", control.proto.help)
        self.assertIn("normalisation", control.proto.help)
        self.assertIn("pas des paramètres optimaux", control.proto.help)
        control.set_value(.75).run()
        app.radio(key="nav").set_value("Réglages").run()
        app.radio(key="nav").set_value("Conception").run()
        self.assertEqual(app.number_input(key="weight_rbs_calculator").value, .75)
        self.assertEqual(self.launch(app).scoring.rbs_calculator, .75)

    def test_relative_region_weights_survive_navigation_and_reach_request(self):
        app = self.app()
        self.assertEqual({key: app.number_input(key=key).value for key in FUNCTIONAL_KEYS},
                         {key: .75 if key.endswith("on_rbs_linker") else 1.0 for key in FUNCTIONAL_KEYS})
        custom = dict(zip(FUNCTIONAL_KEYS, (0.0, 1.0, 2.0, 7.0)))
        for key, value in custom.items():
            app.number_input(key=key).set_value(value)
        app.run()
        self.assertFalse(app.exception)
        captions = " ".join(item.value for item in app.caption)
        self.assertIn("Poids effectifs : OFF 10 % · ON 90 %", captions)
        self.assertIn("Part du compromis : 70 %", captions)
        app.radio(key="nav").set_value("Réglages").run()
        app.radio(key="nav").set_value("Conception").run()
        self.assertFalse(app.exception)
        self.assertEqual({key: app.number_input(key=key).value for key in FUNCTIONAL_KEYS}, custom)
        request = self.launch(app)
        self.assertEqual({key: getattr(request.functional_defect_weights, key.removeprefix("weight_functional_"))
                          for key in FUNCTIONAL_KEYS}, custom)
        self.assertEqual(request.scoring.defect, 6.0)
        self.assert_request_weights(request, DEFAULT_WEIGHTS)

    def test_zero_region_sum_blocks_launch_and_global_coefficient_is_separate(self):
        app = self.app()
        for key in FUNCTIONAL_KEYS:
            app.number_input(key=key).set_value(0.0)
        app.run()
        self.assertFalse(app.exception)
        launch = next(button for button in app.button if button.label == "Lancer la conception")
        self.assertTrue(launch.disabled)
        self.assertTrue(any("Au moins un poids de zone" in item.value for item in app.error))
        self.assertEqual(self.start.call_count, 0)
        app.number_input(key="weight_defect").set_value(0.0).run()
        self.assertTrue(next(button for button in app.button if button.label == "Lancer la conception").disabled)
        app.number_input(key="weight_functional_on_rbs_aug").set_value(1.0).run()
        self.assertFalse(next(button for button in app.button if button.label == "Lancer la conception").disabled)
        self.assertTrue(any("OFF est désactivé" in item.value for item in app.caption))
        request = self.launch(app)
        self.assertEqual(request.scoring.defect, 0.0)
        self.assertEqual(request.functional_defect_weights.on_rbs_aug, 1.0)

    def test_existing_weights_are_never_replaced_by_new_defaults(self):
        for weights in (OLD_WEIGHTS, PREVIOUS_WEIGHTS,
                        {**DEFAULT_WEIGHTS, "weight_leak": .05},
                        {**DEFAULT_WEIGHTS, "weight_leak": .7}):
            with self.subTest(weights=weights):
                app = self.app(weight_trigger_linearity=1., **weights)
                self.assert_weights(app, weights)
                app.radio(key="nav").set_value("Réglages").run()
                app.radio(key="nav").set_value("Conception").run()
                self.assert_weights(app, weights)
                self.assert_request_weights(self.launch(app), weights)
                self.assertEqual(app.session_state["weight_trigger_linearity"], 1.)

    def test_help_describes_boundaries_metric_directions_and_model_limits(self):
        app = self.app()
        help_text = {key: app.number_input(key=key).proto.help for key in DEFAULT_WEIGHTS}
        rbs = help_text["weight_rbs"]
        self.assertIn("ViennaRNA", rbs)
        self.assertIn("fragment ARN isolé", rbs)
        self.assertIn("du préfixe de boucle avant le RBS à la fin du linker", rbs)
        self.assertIn("Un ΔG moins négatif", rbs)
        self.assertIn("après les défauts locaux", rbs)
        self.assertIn("Un poids nul désactive", rbs)
        self.assertIn("aucune exclusion automatique", rbs)

        defect = help_text["weight_defect"]
        self.assertIn("cibles OFF et ON", defect)
        self.assertIn("RBS → fin du linker pèse 20 %", defect)
        self.assertIn("26,7 %", defect)
        self.assertIn("Ces zones se recouvrent", defect)
        self.assertIn("moyenne pondérée des quatre zones", defect)
        self.assertIn("séparément des poids relatifs", defect)
        self.assertIn("pas seulement sur la structure MFE", defect)
        self.assertIn("ne contient pas le défaut de concentration", defect)
        self.assertIn("fractions de 0 à 1", defect)
        self.assertIn("RBS et AUG peuvent être libres en OFF", defect)
        self.assertIn("Ces quatre poids relatifs guident aussi la première phase d’optimisation", defect)
        self.assertIn("Aucun seuil d’exclusion n’est ajouté", defect)
        self.assertIn("Plus faible est meilleur", defect)
        self.assertIn("0 est idéal", defect)

        leak = help_text["weight_leak"]
        self.assertIn("Poids 0", leak)
        self.assertIn("aptamère:trigger:switch de 1:1:1", leak)
        self.assertIn("concentration d’aptamère est égale à celle du trigger", leak)
        self.assertIn("[complexe trigger–switch] / [switch] initial", leak)
        self.assertIn("pas une fuite de traduction en OFF", leak)
        self.assertIn("ligand n’est pas simulé", leak)
        self.assertIn("aptamères sont dans un autre tube", leak)
        self.assertIn("scénario de comparaison", leak)
        self.assertIn("pas une simulation du transfert entre les tubes", leak)
        self.assertIn("n’intervient pas dans le score", leak)
        self.assertIn("exactement 1 trigger et 1 switch", leak)
        self.assertIn("tous les switches introduits", leak)
        self.assertIn("complexes concurrents jusqu’à 3 brins", leak)

        on = help_text["weight_on"]
        self.assertIn("Plus élevé est meilleur", on)
        self.assertIn("trigger:switch de 1:1", on)
        self.assertIn("sans aptamère ni ligand", on)
        self.assertIn("[complexe trigger–switch] / [switch] initial", on)
        self.assertIn("pas directement l’activation de la traduction", on)
        self.assertIn("reste informatif (poids 0)", on)

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
