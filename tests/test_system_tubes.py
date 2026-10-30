"""Original-system tubes run explicitly and remain tied to their physical inputs."""

import copy
from dataclasses import asdict
import json
from pathlib import Path
import threading
import unittest
from unittest.mock import patch

from aptaswitch_studio.system_tube_runtime import SystemTubeJob
from aptaswitch_core.models import ThermoConditions
from aptaswitch_core.tube_analysis import TubeAnalysisConditions
from tests import test_system_sections as system_sections
from tests.test_tube_analysis_ui import analysis_fixture
from aptaswitch_studio.web_runtime import JobSnapshot

button = system_sections.button


class SystemTubeUITests(unittest.TestCase):
    setUp = system_sections.SystemSectionsTests.setUp
    app = system_sections.SystemSectionsTests.app
    fill_sequences = system_sections.SystemSectionsTests.fill_sequences
    open_extension = system_sections.SystemSectionsTests.open_extension
    analysis_result = staticmethod(system_sections.SystemSectionsTests.analysis_result)

    @staticmethod
    def population_result(path, aptamer, trigger, thermo, conditions, **kwargs):
        result = analysis_fixture()
        result.update(aptamer_sequence=aptamer, trigger_sequence=trigger,
                      model=asdict(thermo), conditions=asdict(conditions))
        if kwargs.get("progress_callback"):
            kwargs["progress_callback"](2, 2, "Tubes calculés")
        return result

    def start_and_finish(self, app):
        app.button(key="analyze_system_tubes").click().run()
        job = app.session_state["system_tube_job"]
        job._thread.join(timeout=5)
        self.assertFalse(job._thread.is_alive())
        app.run()
        self.assertFalse(app.exception)
        return job

    def test_explicit_button_preserves_results_across_ligand_changes_and_invalidates_concentrations(self):
        with patch("aptaswitch_studio.system_tube_runtime.analyze_parasite_tubes", side_effect=self.population_result) as compute:
            app = self.app()
            self.fill_sequences(app)
            compute.assert_not_called()
            self.open_extension(app)
            self.assertEqual(app.number_input(key="extension_parasite_top").value, 5)
            app.radio(key="design_section").set_value("system").run()
            self.assertFalse(app.exception)
            job = self.start_and_finish(app)
            compute.assert_called_once()
            conditions = compute.call_args.args[4]
            # New native calculations support complexes of at most three strands.
            self.assertEqual(conditions.max_complex_size, 3)
            self.assertAlmostEqual(conditions.volume_l, 20e-6)
            self.assertEqual(job.snapshot().state, "done")
            self.assertTrue(any(metric.label == "Trigger dans des complexes parasites" for metric in app.metric))
            app.checkbox(key="extension_ligand_base_3").check().run()
            compute.assert_called_once()
            self.assertTrue(any(metric.label == "Trigger dans des complexes parasites" for metric in app.metric))
            self.assertTrue(any(download.key == "system_tubes_json" for download in app.get("download_button")))
            app.number_input(key="system_trigger_concentration").set_value(7.0).run()
            compute.assert_called_once()
            self.assertFalse(any(metric.label == "Trigger dans des complexes parasites" for metric in app.metric))
            self.start_and_finish(app)
            self.assertEqual(compute.call_count, 2)
            self.assertAlmostEqual(compute.call_args.args[4].trigger_concentration_m, 7e-6)

    def test_physical_conditions_and_screening_count_are_sent_to_extension_request(self):
        with patch("aptaswitch_studio.extension_ui.ExtensionRunJob") as job_class:
            job_class.return_value.snapshot.return_value = JobSnapshot("running", 0, 1, "Préparation…")
            app = self.app()
            self.fill_sequences(app)
            app.number_input(key="system_aptamer_concentration").set_value(3.0)
            app.number_input(key="system_trigger_concentration").set_value(8.0)
            app.number_input(key="system_tube_volume_ul").set_value(50.0).run()
            self.open_extension(app)
            app.number_input(key="extension_parasite_top").set_value(7).run()
            button(app, "Étendre le trigger").click().run()
            request = job_class.call_args.args[0]
            self.assertEqual(request.parasite_top, 7)
            self.assertAlmostEqual(request.tube_conditions.aptamer_concentration_m, 3e-6)
            self.assertAlmostEqual(request.tube_conditions.trigger_concentration_m, 8e-6)
            self.assertAlmostEqual(request.tube_conditions.volume_l, 50e-6)
            self.assertEqual(request.tube_conditions.max_complex_size, 3)
            job_class.return_value.start.assert_called_once()

    def test_failed_original_analysis_displays_error_without_fake_metrics(self):
        with patch("aptaswitch_studio.system_tube_runtime.analyze_parasite_tubes", side_effect=RuntimeError("Engine failed")):
            app = self.app()
            self.fill_sequences(app)
            self.start_and_finish(app)
            self.assertTrue(any("Engine failed" in error.value for error in app.error))
            self.assertFalse(any(metric.label == "Trigger dans des complexes parasites" for metric in app.metric))

    def test_changing_input_cancels_pending_job_and_keeps_expensive_jobs_disabled(self):
        release = threading.Event()
        started = threading.Event()

        def long_analysis(*args, **kwargs):
            started.set()
            release.wait(10)
            return {"status": "cancelled", "tubes": {}}

        with patch("aptaswitch_studio.system_tube_runtime.analyze_parasite_tubes", side_effect=long_analysis):
            app = self.app()
            self.fill_sequences(app)
            app.text_area(key="trigger_input").set_value("ACGT" * 6).run()
            # A valid standard switch can start before the tube worker runs.
            app.radio(key="design_section").set_value("switch").run()
            self.assertFalse(app.exception)
            self.assertFalse(button(app, "Lancer la conception").disabled)
            app.radio(key="design_section").set_value("system").run()
            app.button(key="analyze_system_tubes").click().run()
            job = app.session_state["system_tube_job"]
            try:
                self.assertTrue(started.wait(2))
                self.open_extension(app)
                self.assertTrue(button(app, "Étendre le trigger").disabled)
                app.radio(key="design_section").set_value("switch").run()
                self.assertFalse(app.exception)
                self.assertTrue(button(app, "Lancer la conception").disabled)
                app.radio(key="design_section").set_value("system").run()
                app.number_input(key="extension_temperature").set_value(25.0).run()
                self.assertTrue(job._cancel.is_set())
                self.assertTrue(app.button(key="analyze_system_tubes").disabled)
            finally:
                release.set()
                job._thread.join(5)
            app.run()
            self.assertFalse(app.exception)
            self.assertFalse(app.button(key="analyze_system_tubes").disabled)

    def test_extended_results_display_saved_tubes_and_unscreened_candidates(self):
        result = json.loads((Path(__file__).parent / "fixtures" / "extension_native_fixture.json").read_text())
        original = self.app()
        candidates = result["candidates"]
        result["reference"]["tube_analysis"] = analysis_fixture()
        result["parasite_analysis"] = {"status": "completed", "completed_count": 1, "selected_count": 1, "top_per_length": 1}
        candidates[0]["tube_analysis"] = analysis_fixture()
        candidates[0].update(trigger_only_parasite_pct=40.0, mixed_trigger_parasite_pct=20.0, mixed_trigger_target_pct=50.0)
        expected = copy.deepcopy(result)
        with patch("aptaswitch_studio.system_tube_runtime.analyze_parasite_tubes") as compute:
            original.session_state["extension_result"] = result
            original.session_state["nav"] = "Résultats"
            original.session_state["results_kind"] = "Extension de triggers"
            original.run()
            self.assertFalse(original.exception)
            compute.assert_not_called()
            self.assertTrue(any("Parasites · trigger seul (%)" in table.value.columns for table in original.table))
            self.assertTrue(any(metric.label == "Trigger dans des complexes parasites" for metric in original.metric))
            self.assertEqual(original.session_state["extension_result"], expected)


class SystemTubeJobTests(unittest.TestCase):
    def test_cancel_before_start_does_not_call_native_partition(self):
        from tests.native_ui_fixtures import patch_native_runtime
        from unittest.mock import Mock
        backend = patch_native_runtime(self)
        backend.partition = Mock(side_effect=AssertionError("Cancelled job must not calculate"))
        signature = ("ACGT", "TGCA", "/fixture/native", ThermoConditions(switch_material="dna"), TubeAnalysisConditions())
        job = SystemTubeJob(signature, backend=backend)
        job.cancel()
        job.start()
        job._thread.join(5)
        self.assertFalse(job._thread.is_alive())
        backend.partition.assert_not_called()
        self.assertEqual(job.snapshot().state, "done", job.snapshot().error)
        self.assertEqual(job.snapshot().result["status"], "cancelled")
        self.assertTrue(job.consume_once())
        self.assertFalse(job.consume_once())


if __name__ == "__main__":
    unittest.main()
