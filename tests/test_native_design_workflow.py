"""Native app routing and portable provenance exercised with the real library."""
from dataclasses import asdict, replace
from copy import deepcopy
import json
import os
from pathlib import Path
import shutil
import tempfile
import unittest
from unittest.mock import patch

from aptaswitch_core.engines import run_design
from aptaswitch_core.models import ReporterContext, RunRequest, SequenceInput, SwitchArchitecture, ThermoConditions
from aptaswitch_core.thermodynamics.optimization import TrialOptions
from aptaswitch_core.thermodynamics.setup import ENGINE_ID, MODEL_ID
from aptaswitch_core.workflow import run_and_export
from aptaswitch_studio.run_library import load_run
from aptaswitch_studio.saved_switch_design import load_saved_design
from aptaswitch_studio.saved_switch_library import candidate_snapshot, save_switch


LIBRARY = os.environ.get("APTASWITCH_HYBRID_NATIVE_LIBRARY", "")


class NativeRequestTests(unittest.TestCase):
    def test_native_defaults_and_model_domain_are_explicit(self):
        thermo = ThermoConditions()
        self.assertEqual(thermo.material, MODEL_ID)
        thermo.validate_model()
        for temperature in (24.9, 55.1):
            with self.subTest(temperature=temperature), self.assertRaises(ValueError):
                replace(thermo, temperature_c=temperature).validate_model()
        with self.assertRaises(ValueError):
            replace(thermo, material="external-archive", temperature_c=80).validate_model()

    def test_legacy_engine_cannot_silently_compute_under_new_identity(self):
        request = RunRequest(SequenceInput("GGTT", "AACC"), SwitchArchitecture(2, 2),
                             engine="external-archive", thermo=ThermoConditions(material="experimental-Sugimoto95-Vienna04-loopmean"))
        with self.assertRaisesRegex(ValueError, "native-hybrid"):
            run_design(request)


class NativeArchiveValidationTests(unittest.TestCase):
    def test_target_status_requires_both_consistent_objectives_to_reach_their_targets(self):
        from tests.design_fixtures import design_result_fixture
        from aptaswitch_studio.run_library import _load_design
        request = RunRequest(SequenceInput("GGTT", "AACC"), SwitchArchitecture(2, 2), trials=1)
        original = design_result_fixture(request).to_dict()
        for kind in ("local", "global"):
            payload = deepcopy(original)
            report = payload["optimization"]["trials"][0]
            if kind == "local":
                report["best_loss"] = .2
                report["best_objective"]["local_defects"] = [.2] * 4
            else:
                report["best_global_loss"] = .2
                report["best_objective"]["global_defect"] = .2
            with self.subTest(kind=kind), self.assertRaisesRegex(ValueError, "target was not reached"):
                _load_design(payload, Path("/archive/results.json"))

    def test_native_archives_require_model_conditions_and_real_optimizer_outcomes(self):
        from tests.design_fixtures import design_result_fixture
        from aptaswitch_studio.run_library import _load_design
        request = RunRequest(SequenceInput("GGTT", "AACC"), SwitchArchitecture(2, 2), trials=2)
        payload = design_result_fixture(request).to_dict()
        mutations = [
            lambda p: p.pop("thermodynamics"),
            lambda p: p.update(thermodynamics={}),
            lambda p: p["thermodynamics"].pop("parameter_digest"),
            lambda p: p["thermodynamics"].pop("engine"),
            lambda p: p["thermodynamics"]["conditions"].update(celsius=40),
            lambda p: p["thermodynamics"]["conditions"].update(sodium=True),
            lambda p: p.pop("optimization"),
            lambda p: p.update(optimization={}),
            lambda p: p["optimization"].update(algorithm="unknown"),
            lambda p: p["optimization"].update(objective="unknown"),
            lambda p: p["optimization"].update(seed_rule="random"),
            lambda p: p["optimization"]["functional_weights"].update(off_toehold=2),
            lambda p: p["optimization"]["options"].update(max_evaluations=75),
            lambda p: p["optimization"].update(trials=[]),
            lambda p: p["optimization"]["trials"][0].pop("options"),
            lambda p: p["optimization"]["trials"][1].update(trial=1),
            lambda p: p["optimization"]["trials"][1]["options"].update(seed=0),
            lambda p: p["optimization"]["trials"][0].update(evaluations=10**6),
            lambda p: p["optimization"]["trials"][0].update(proposals=True),
            lambda p: p["optimization"]["trials"][0].update(best_loss=-.1),
            lambda p: p["optimization"]["trials"][0].update(best_loss=.2),
            lambda p: p["optimization"]["trials"][0].update(best_global_loss=.2),
            lambda p: p["optimization"]["trials"][0]["best_objective"].update(global_defect=.2),
            lambda p: p["optimization"]["trials"][0]["best_objective"].update(local_defects=[.2]*4),
            lambda p: p["optimization"]["trials"][0].update(phase="unknown"),
            lambda p: p["optimization"]["trials"][0].update(phase="global"),
            lambda p: p["optimization"]["trials"][0].update(population_evaluated=2),
            lambda p: p["optimization"]["trials"][0].update(evaluations=0),
            lambda p: p["optimization"]["trials"][0].update(stop_reason="budget_exhausted", budget_kind="evaluations"),
            lambda p: p["optimization"]["trials"][0].update(budget_kind="proposals"),
            lambda p: p["optimization"]["trials"][0].update(budget_kind=[]),
            lambda p: p["optimization"]["trials"][0].update(stop_reason=[]),
            lambda p: p["results"][0].update(trial=3),
        ]
        with patch("aptaswitch_core.thermodynamics.setup.get_backend", side_effect=AssertionError("Offline archive")):
            _load_design(deepcopy(payload), Path("/archive/results.json"))
            for index, mutate in enumerate(mutations):
                changed = deepcopy(payload)
                mutate(changed)
                with self.subTest(mutation=index), self.assertRaises(ValueError):
                    _load_design(changed, Path("/archive/results.json"))

    def test_unsupported_archive_is_not_relabelled_as_native(self):
        from tests.design_fixtures import design_result_fixture
        from aptaswitch_studio.run_library import _load_design
        request = RunRequest(SequenceInput("GGTT", "AACC"), SwitchArchitecture(2, 2), trials=1)
        payload = design_result_fixture(request).to_dict()
        payload["request"]["engine"] = "external-archive"
        payload["request"]["thermo"]["material"] = "external-archive"
        payload.pop("optimization")
        payload.pop("thermodynamics")
        with patch("aptaswitch_core.thermodynamics.setup.get_backend", side_effect=AssertionError("Offline archive")):
            with self.assertRaises(ValueError):
                _load_design(payload, Path("/archive/results.json"))


@unittest.skipUnless(LIBRARY, "Explicit compiled hybrid library required")
class NativeDesignWorkflowTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.folder = Path(temporary.name).resolve()
        self.request = RunRequest(
            SequenceInput("GGTT", "AACC"),
            SwitchArchitecture(2, 2, upper_stem_length=1, upper_stem2_length=0,
                               rbs_loop="AGGA", t7_leader="", rbs_prefix_length=0),
            reporter=ReporterContext(key="short", label="Short test reporter", sequence_after_start_rna="GCU"),
            thermo=ThermoConditions(temperature_c=37, validation_max_complex_size=2),
            trials=2, output_dir=self.folder / "designs", thermodynamics_library_path=LIBRARY,
            exclude_stop_candidates=False, exclude_non_linear_candidates=False,
            optimization_options=asdict(TrialOptions(seed=17, max_evaluations=1, max_proposals=1)),
        )

    def test_real_design_export_snapshot_and_native_cache_preserve_identity(self):
        progress = []
        result, exports = run_and_export(self.request, progress_callback=lambda *event: progress.append(event))
        self.assertEqual(result.status, "completed")
        self.assertEqual(result.request.engine, ENGINE_ID)
        self.assertEqual(result.thermodynamics["model"], MODEL_ID)
        self.assertTrue(result.thermodynamics["parameter_digest"].startswith("sha256:"))
        self.assertTrue(all(isinstance(event[2], str) for event in progress))
        trials = result.optimization["trials"]
        self.assertEqual([row["options"]["seed"] for row in trials], [17, 18])
        self.assertEqual([row["stop_reason"] for row in trials], ["budget_exhausted"] * 2)
        self.assertTrue(all(row["evaluations"] == 1 for row in trials))
        loaded = load_run(Path(exports["json"])).result
        self.assertEqual(loaded.to_dict(), result.to_dict())
        self.assertEqual(json.loads(Path(exports["json"]).read_text())["optimization"], result.optimization)
        self.assertIn(MODEL_ID, Path(exports["csv"]).read_text())

        from aptaswitch_core.structure_ensemble import analyze_structure_ensemble, load_ensemble, save_ensemble
        candidate = result.candidates[0]
        families = analyze_structure_ensemble(result.request, candidate, num_sample=100, seed=23)
        save_ensemble(result.output_dir, families)
        snapshot = candidate_snapshot(result, candidate, "Native retained")
        saved = save_switch(self.folder / "selections", snapshot["name"], snapshot,
                            source_output_dir=result.output_dir)
        shutil.rmtree(result.output_dir)
        with patch("aptaswitch_core.thermodynamics.setup.get_backend", side_effect=AssertionError("Read-only restoration")):
            restored, retained = load_saved_design(saved)
        self.assertEqual(retained.to_dict(), candidate.to_dict())
        self.assertEqual(restored.thermodynamics, result.thermodynamics)
        self.assertEqual(restored.optimization["trials"],
                         [row for row in result.optimization["trials"] if row["trial"] == candidate.trial])
        self.assertEqual(restored.optimization["source_trial_count"], result.request.trials)
        self.assertEqual(restored.optimization["source_run_id"], result.run_id)
        self.assertEqual(restored.optimization["scope"], "retained_candidate")
        self.assertEqual(load_ensemble(restored.output_dir, families["inputs"]), families)

    def test_cancellation_exports_no_invented_metrics_or_completed_trials(self):
        result, exports = run_and_export(self.request, cancel_check=lambda: True)
        self.assertEqual(result.status, "cancelled")
        self.assertEqual(result.candidates, [])
        self.assertEqual(result.optimization["trials"], [])
        self.assertEqual(result.thermodynamics["model"], MODEL_ID)
        self.assertEqual(load_run(Path(exports["json"])).result.status, "cancelled")

    def test_global_phase_preserves_guards_and_reopens_with_separate_objectives(self):
        request = replace(self.request, trials=1, optimization_options=asdict(TrialOptions(
            seed=17, max_evaluations=12, max_proposals=120, population_size=4, lineages=2,
            target_local_loss=0., target_global_loss=0.,
        )))
        result, exports = run_and_export(request)
        report, = result.optimization["trials"]
        self.assertEqual(report["phase"], "global")
        self.assertEqual(len(report["local_limits"]), 4)
        self.assertEqual(len(report["guard_limits"]), 3)
        self.assertEqual(report["best_global_loss"], report["best_objective"]["global_defect"])
        self.assertLessEqual(report["evaluations"], 12)
        self.assertLessEqual(report["proposals"], 120)
        self.assertEqual(load_run(Path(exports["json"])).result.optimization, result.optimization)

    def test_model_mismatch_in_recorded_provenance_is_rejected(self):
        result, exports = run_and_export(replace(self.request, trials=1), cancel_check=lambda: True)
        path = Path(exports["json"])
        payload = json.loads(path.read_text())
        payload["thermodynamics"]["model"] = "external-archive"
        path.write_text(json.dumps(payload))
        with self.assertRaisesRegex(ValueError, "differs"):
            load_run(path)


if __name__ == "__main__":
    unittest.main()
