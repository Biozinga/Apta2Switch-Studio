from __future__ import annotations

import csv
import json
import sys
import tempfile
import threading
import unittest
import zipfile
from copy import deepcopy
from dataclasses import replace
from pathlib import Path
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from aptaswitch_core.models import ThermoConditions
from aptaswitch_core.tube_analysis import TubeAnalysisConditions
from aptaswitch_studio.extension_runtime import (
    ExtensionRunJob,
    ExtensionRunRequest,
    run_extension_and_export,
    request_from_extension_result,
)
from tests.export_fixtures import export_extension_result


def request(output_dir: Path, *, trigger: str = "ACGT") -> ExtensionRunRequest:
    return ExtensionRunRequest(
        aptamer_dna="ACGTACGT", trigger_dna=trigger, target_lengths=(5, 6),
        thermo=ThermoConditions(material="experimental-Sugimoto95-Vienna04-loopmean", temperature_c=37, switch_material="dna"),
        import_path="/configured/native", output_dir=output_dir,
        molecule_name="=custom system", top=12, max_ensemble_candidates=40, parasite_top=0,
    )


def core_result(status: str = "completed") -> dict:
    return {
        "status": status,
        "input": {"aptamer_dna": "ACGTACGT", "trigger_dna": "ACGT", "target_lengths": [5, 6]},
        "reference": {"mfe_structure": "........+....", "mfe_energy_kcal_mol": -1.2},
        "mfe_screen": {"total_candidates_expected": 20, "total_candidates_screened": 20},
        "candidates": [{
            "rank": 1, "rank_within_length": 1, "final_score": 0.25,
            "trigger_length": 5, "extension_length": 1, "extension_5prime": "T",
            "extended_trigger": "TACGT", "mfe_structure": "........+.....",
            "mfe_energy_kcal_mol": -1.21,
            "reference_pairs": [{"i": 0, "j": 9, "reference_probability": 0.8, "candidate_probability": 0.79}],
        }],
        "curves": {"reference_pairs": []},
        "warnings": [],
    }


def write_visual_fixture(result: dict, output_dir: Path) -> dict[str, str]:
    path = output_dir / f"{result['run_id']}_reference.svg"
    path.write_text('<svg xmlns="http://www.w3.org/2000/svg"/>', encoding="utf-8")
    return {"svg_reference": str(path)}


def tube_result(trigger: str, conditions: TubeAnalysisConditions) -> dict:
    return {
        "status": "completed", "aptamer_sequence": "ACGTACGT", "trigger_sequence": trigger,
        "model": {"material": "experimental-Sugimoto95-Vienna04-loopmean", "aptamer_material": "dna", "trigger_material": "dna", "temperature_c": 37, "sodium_m": 0.12, "magnesium_m": 0},
        "conditions": {
            "aptamer_concentration_m": conditions.aptamer_concentration_m,
            "trigger_concentration_m": conditions.trigger_concentration_m,
            "volume_l": conditions.volume_l, "max_complex_size": conditions.max_complex_size,
        },
        "tubes": {
            "trigger_only": {
                "summary": {"trigger_parasite_pct": 20, "trigger_target_pct": 0},
                "complexes": [{
                    "label": "2 T", "aptamer_count": 0, "trigger_count": 2, "strand_count": 2,
                    "category": "parasite", "concentration_m": 5e-7,
                    "complex_molecule_count": 6.02214076e12, "trigger_fraction_pct": 20,
                }],
            },
            "aptamer_trigger": {
                "summary": {"trigger_parasite_pct": 10, "trigger_target_pct": 75},
                "complexes": [{
                    "label": "A + T", "aptamer_count": 1, "trigger_count": 1, "strand_count": 2,
                    "category": "target", "concentration_m": 3.75e-6,
                    "complex_molecule_count": 4.51660557e13, "trigger_fraction_pct": 75,
                }],
            },
        },
    }


class ExtensionRuntimeTests(unittest.TestCase):
    def setUp(self):
        from aptaswitch_core.localization import language_context
        self.enterContext(language_context("fr"))
        from test_thermodynamics_context_analysis import FixedWeightBackend
        self.backend = FixedWeightBackend()
        self.enterContext(patch("aptaswitch_studio.extension_runtime.get_backend", return_value=self.backend))

    def test_advance_event_is_distinct_from_cancellation_and_only_for_screen(self):
        job = ExtensionRunJob(request(Path("/unused")))
        job.advance_to_ensemble()
        self.assertFalse(job._advance.is_set())
        job._set_phase("mfe")
        job._progress(3, 20, "MFE")
        self.assertTrue(job.can_advance)
        job.advance_to_ensemble()
        self.assertTrue(job._advance.is_set())
        self.assertFalse(job._cancel.is_set())
        self.assertFalse(job.can_advance)
        job._set_phase("ensemble")
        self.assertFalse(job.can_advance)

    def test_resume_uses_saved_request_and_writes_atomic_checkpoint(self):
        with tempfile.TemporaryDirectory() as tmp:
            original = request(Path("/old/location"))
            saved = core_result("cancelled")
            saved["request"] = original.to_dict()
            resumed = request_from_extension_result(saved, import_path="/new/native", output_dir=Path(tmp))
            self.assertEqual(resumed.trigger_dna, original.trigger_dna)
            self.assertEqual(resumed.thermo, original.thermo)
            self.assertEqual(resumed.import_path, "")
            from dataclasses import asdict
            from aptaswitch_core.extension_thermodynamics import extension_model_settings
            saved["thermodynamics"] = asdict(self.backend.provenance)
            saved["model"] = extension_model_settings(resumed.thermo, self.backend)

            def search(**kwargs):
                self.assertIs(kwargs["resume_result"], saved)
                kwargs["checkpoint_callback"]({**core_result("running"), "mfe_survivors": []})
                return core_result()

            with (patch("aptaswitch_studio.extension_runtime.design_trigger_extensions", side_effect=search),
                  patch("aptaswitch_studio.extension_runtime.export_extension_visuals", return_value={})):
                result, _ = run_extension_and_export(resumed, resume_result=saved)
            checkpoint = Path(result["output_dir"]) / "mfe_checkpoint.json"
            payload = json.loads(checkpoint.read_text())
            self.assertEqual(payload["request"]["trigger_dna"], original.trigger_dna)
            self.assertIn("mfe_survivors", payload)
            self.assertFalse(checkpoint.with_suffix(".json.tmp").exists())

    def test_resume_caps_legacy_catalogs_at_three_without_editing_saved_results(self):
        saved = core_result("cancelled")
        saved["request"] = request(Path("/old/location")).to_dict()
        saved["request"]["tube_conditions"]["max_complex_size"] = 5
        saved["request"]["thermo"]["screening_max_complex_size"] = 5
        saved["request"]["thermo"]["validation_max_complex_size"] = 5
        original = deepcopy(saved)
        resumed = request_from_extension_result(saved, output_dir=Path("/new/location"))
        self.assertEqual(resumed.tube_conditions.max_complex_size, 3)
        self.assertEqual(resumed.thermo.screening_max_complex_size, 3)
        self.assertEqual(resumed.thermo.validation_max_complex_size, 3)
        self.assertEqual(saved, original)

    def test_resume_rejects_unknown_or_missing_engine_without_relabeling_saved_results(self):
        for engine in ("external-archive", None):
            with self.subTest(engine=engine):
                saved = core_result("cancelled")
                saved["request"] = request(Path("/old/location")).to_dict()
                if engine is None:
                    saved["request"].pop("engine")
                else:
                    saved["request"]["engine"] = engine
                original = deepcopy(saved)
                with self.assertRaisesRegex(ValueError, "native-hybrid"):
                    request_from_extension_result(saved, output_dir=Path("/new/location"))
                self.assertEqual(saved, original)

    def test_ligand_annotation_is_exported_without_changing_candidate_scores(self):
        with tempfile.TemporaryDirectory() as tmp:
            run_request = replace(request(Path(tmp)), ligand_aptamer_positions=(1, 4, 8))
            with (
                patch("aptaswitch_studio.extension_runtime.design_trigger_extensions", return_value=core_result()),
                patch("aptaswitch_studio.extension_runtime.export_extension_visuals", side_effect=write_visual_fixture),
            ):
                result, exports = run_extension_and_export(run_request)
            saved = json.loads(Path(exports["json"]).read_text())
            self.assertEqual(saved["input"]["ligand_aptamer_positions"], [1, 4, 8])
            self.assertEqual(saved["request"]["ligand_aptamer_positions"], [1, 4, 8])
            self.assertEqual(result["candidates"], core_result()["candidates"])

    def test_search_results_and_all_exports_retain_provenance_and_structures(self):
        import openpyxl

        with tempfile.TemporaryDirectory() as tmp:
            run_request = request(Path(tmp))
            progress = []

            def search(**kwargs):
                kwargs["progress_callback"](12, 20, "Criblage · 12 extensions évaluées.")
                kwargs["progress_callback"](20, 20, "Validation thermodynamique terminée.")
                return core_result()

            with (
                patch("aptaswitch_studio.extension_runtime.design_trigger_extensions", side_effect=search) as engine,
                patch("aptaswitch_studio.extension_runtime.export_extension_visuals", side_effect=write_visual_fixture),
            ):
                result, exports = run_extension_and_export(
                    run_request, progress_callback=lambda *args: progress.append(args),
                )

            self.assertEqual(engine.call_count, 1)
            self.assertEqual(engine.call_args.kwargs["target_lengths"], (5, 6))
            self.assertEqual(engine.call_args.kwargs["thermo"], run_request.thermo)
            self.assertEqual(engine.call_args.kwargs["import_path"], "/configured/native")
            self.assertEqual(engine.call_args.kwargs["top"], 12)
            self.assertEqual(engine.call_args.kwargs["max_ensemble_candidates"], 40)
            self.assertEqual(set(exports), {"json", "csv", "xlsx", "log", "zip", "svg_reference"})
            self.assertEqual(Path(result["output_dir"]).parent, Path(tmp) / "extensions")
            self.assertRegex(result["run_id"], r"^\d{8}_\d{6}_extension_[a-f0-9]{6}$")
            saved = json.loads(Path(exports["json"]).read_text())
            self.assertEqual(saved["candidates"], core_result()["candidates"])
            self.assertEqual(saved["reference"], core_result()["reference"])
            self.assertEqual(saved["request"]["molecule_name"], "=custom system")
            self.assertEqual(saved["request"]["output_dir"], tmp)
            with Path(exports["csv"]).open(encoding="utf-8-sig", newline="") as handle:
                rows = list(csv.DictReader(handle))
            self.assertEqual(rows[0]["extended_trigger"], "TACGT")
            self.assertEqual(rows[0]["mfe_structure"], "........+.....")
            self.assertEqual(json.loads(rows[0]["reference_pairs"]), core_result()["candidates"][0]["reference_pairs"])
            workbook = openpyxl.load_workbook(exports["xlsx"])
            self.assertIn("Référence", workbook.sheetnames)
            name_cells = [row[1] for row in workbook["Paramètres"] if row[0].value == "molecule_name"]
            self.assertEqual(name_cells[0].value, "=custom system")
            self.assertEqual(name_cells[0].data_type, "s")
            workbook.close()
            journal = Path(exports["log"]).read_text()
            self.assertRegex(journal, r"\[\d{2}:\d{2}:\d{2}\] Criblage · 12 extensions évaluées\.")
            self.assertIn("Extension terminée — résultats exportés.", journal)
            with zipfile.ZipFile(exports["zip"]) as archive:
                self.assertEqual(set(archive.namelist()), {Path(path).name for key, path in exports.items() if key != "zip"})
            self.assertTrue(any("structures 2D" in item[2] for item in progress))

    def test_progress_reports_elapsed_eta_and_normalizes_input_length(self):
        with patch("aptaswitch_studio.extension_runtime.time.monotonic", side_effect=[100.0, 110.0]):
            job = ExtensionRunJob(request(Path("/unused"), trigger="AC GT\n"))
            job._progress(5, 20, "Criblage MFE")
            snapshot = job.snapshot()
        self.assertEqual(snapshot.total, 20)
        self.assertEqual(snapshot.current, 5)
        self.assertEqual(snapshot.elapsed_seconds, 10)
        self.assertEqual(snapshot.eta_seconds, 30)
        self.assertRegex(snapshot.log[-1], r"^\[\d{2}:\d{2}:\d{2}\] Criblage MFE$")
        self.assertFalse(job.consume_once())

    def test_cancellation_reaches_engine_and_exports_partial_results(self):
        entered = threading.Event()
        released = threading.Event()
        with tempfile.TemporaryDirectory() as tmp:
            job = ExtensionRunJob(request(Path(tmp)))

            def search(**kwargs):
                entered.set()
                self.assertTrue(released.wait(timeout=5))
                self.assertTrue(kwargs["cancel_check"]())
                kwargs["progress_callback"](4, 20, "Annulation · 4 extensions évaluées.")
                result = core_result("cancelled")
                result["mfe_screen"]["total_candidates_screened"] = 4
                return result

            with (
                patch("aptaswitch_studio.extension_runtime.design_trigger_extensions", side_effect=search),
                patch("aptaswitch_studio.extension_runtime.export_extension_visuals", side_effect=write_visual_fixture),
            ):
                job.start()
                try:
                    self.assertTrue(entered.wait(timeout=5))
                    job.cancel()
                finally:
                    released.set()
                job._thread.join(timeout=5)
                self.assertFalse(job._thread.is_alive())
            snapshot = job.snapshot()
            self.assertEqual(snapshot.state, "done")
            self.assertEqual(snapshot.result["status"], "cancelled")
            self.assertEqual(snapshot.current, 4)
            self.assertEqual(snapshot.eta_seconds, 0)
            self.assertIn("partiels", snapshot.message)
            self.assertTrue(Path(snapshot.exports["zip"]).is_file())
            self.assertTrue(job.consume_once())
            self.assertFalse(job.consume_once())

    def test_eta_resets_when_search_moves_to_ensemble_validation(self):
        with patch("aptaswitch_studio.extension_runtime.time.monotonic", side_effect=[100.0, 120.0, 130.0]):
            job = ExtensionRunJob(request(Path("/unused")))
            job._progress(20, 20, "Criblage terminé")
            job._progress(0, 4, "Validation d’ensemble")
            job._progress(1, 4, "Validation d’ensemble 1/4")
            snapshot = job.snapshot()
        self.assertEqual(snapshot.elapsed_seconds, 30)
        self.assertEqual(snapshot.eta_seconds, 30)

    def test_engine_failure_is_reported_without_replacing_it_with_fake_results(self):
        job = ExtensionRunJob(request(Path("/unused")))
        with patch("aptaswitch_studio.extension_runtime.design_trigger_extensions", side_effect=RuntimeError("Moteur natif indisponible")):
            job._run()
        snapshot = job.snapshot()
        self.assertEqual(snapshot.state, "failed")
        self.assertEqual(snapshot.error, "Moteur natif indisponible")
        self.assertIsNone(snapshot.result)
        self.assertIsNone(snapshot.exports)
        self.assertTrue(any("ERREUR" in line for line in snapshot.log))
        self.assertTrue(job.consume_once())
        self.assertFalse(job.consume_once())

    def test_empty_cancelled_run_can_be_exported(self):
        with tempfile.TemporaryDirectory() as tmp:
            result = core_result("cancelled")
            result.update({"run_id": "cancelled", "started_at": "now", "completed_at": "now", "output_dir": tmp})
            result["candidates"] = []
            result["reference"] = {}
            exports = export_extension_result(result)
            with Path(exports["csv"]).open(encoding="utf-8-sig", newline="") as handle:
                self.assertEqual(list(csv.DictReader(handle)), [])
            self.assertTrue(Path(exports["zip"]).is_file())

    def test_xlsx_preserves_long_curves_and_probability_data_without_cell_truncation(self):
        import openpyxl

        with tempfile.TemporaryDirectory() as tmp:
            result = core_result()
            result.update({"run_id": "long_data", "started_at": "now", "completed_at": "now", "output_dir": tmp})
            points = [{"analyzed": index, "final_score": index / 1000} for index in range(1500)]
            result["curves"] = {"ensemble_progress": points}
            probabilities = [{"index": index, "probability": 0.123456789} for index in range(1500)]
            result["candidates"][0]["probabilities"] = probabilities
            with patch("aptaswitch_studio.extension_runtime.export_extension_visuals", return_value={}):
                exports = export_extension_result(result)
            workbook = openpyxl.load_workbook(exports["xlsx"])
            self.assertEqual(workbook["Progression ensemble"].max_row, 1501)
            self.assertEqual(workbook["Progression ensemble"].cell(1501, 1).value, 1499)
            chunks = [row[3].value for row in workbook["Données longues"].iter_rows(min_row=2)]
            self.assertEqual(json.loads("".join(chunks)), probabilities)
            self.assertTrue(all(len(chunk) <= 30000 for chunk in chunks))
            workbook.close()

    def test_tube_analysis_default_and_request_validation(self):
        self.assertEqual(ExtensionRunRequest.__dataclass_fields__["parasite_top"].default, 5)
        base = request(Path("/unused"))
        for value in (-1, 1.5, True, "5"):
            with self.subTest(value=value), self.assertRaises(ValueError):
                replace(base, parasite_top=value)
        custom = replace(base, parasite_top=3, tube_conditions=TubeAnalysisConditions(volume_l=10e-6))
        self.assertEqual(custom.to_dict()["tube_conditions"]["volume_l"], 10e-6)

    def test_best_candidates_per_length_and_original_are_analyzed_and_exported(self):
        import openpyxl

        with tempfile.TemporaryDirectory() as tmp:
            run_request = replace(request(Path(tmp)), parasite_top=2,
                                  tube_conditions=TubeAnalysisConditions(volume_l=10e-6))
            scientific = core_result()
            original_candidate = scientific["candidates"][0]
            # Deliberately shuffled global order: selection follows each length's rank.
            scientific["candidates"] = [
                {**original_candidate, "rank": rank, "rank_within_length": within,
                 "trigger_length": length, "extended_trigger": sequence, "final_score": rank / 10}
                for rank, within, length, sequence in (
                    (4, 3, 5, "CACGT"), (5, 2, 6, "TTACGT"), (2, 1, 6, "AAACGT"),
                    (1, 1, 5, "TACGT"), (3, 2, 5, "AACGT"),
                )
            ]
            original_order = [row["extended_trigger"] for row in scientific["candidates"]]
            original_scores = [row["final_score"] for row in scientific["candidates"]]

            def analyze(backend, aptamer, trigger, thermo, conditions, **kwargs):
                self.assertIs(backend, self.backend)
                self.assertEqual(aptamer, run_request.aptamer_dna)
                self.assertEqual(thermo, run_request.thermo)
                self.assertEqual(conditions, run_request.tube_conditions)
                kwargs["progress_callback"](2, 2, "Conservation vérifiée")
                return tube_result(trigger, conditions)

            progress = []
            with (
                patch("aptaswitch_studio.extension_runtime.design_trigger_extensions", return_value=scientific),
                patch("aptaswitch_core.thermodynamics.application_analysis.analyze_parasite_tubes", side_effect=analyze) as engine,
                patch("aptaswitch_studio.extension_runtime.export_extension_visuals", side_effect=write_visual_fixture),
            ):
                result, exports = run_extension_and_export(run_request, progress_callback=lambda *args: progress.append(args))
            self.assertEqual([call.args[2] for call in engine.call_args_list], ["ACGT", "TACGT", "AACGT", "AAACGT", "TTACGT"])
            self.assertEqual([row["extended_trigger"] for row in result["candidates"]], original_order)
            self.assertEqual([row["final_score"] for row in result["candidates"]], original_scores)
            self.assertNotIn("tube_analysis", result["candidates"][0])
            self.assertEqual(result["candidates"][1]["mixed_trigger_target_pct"], 75)
            self.assertEqual(result["reference"]["trigger_only_parasite_pct"], 20)
            self.assertEqual(result["parasite_analysis"]["selected_count"], 4)
            self.assertEqual(result["parasite_analysis"]["completed_count"], 4)
            self.assertEqual(result["parasite_analysis"]["status"], "completed")
            self.assertTrue(any(item[0:2] == (10, 10) for item in progress))
            saved = json.loads(Path(exports["json"]).read_text())
            self.assertEqual(saved["reference"]["tube_analysis"]["conditions"]["volume_l"], 10e-6)
            with Path(exports["tube_csv"]).open(encoding="utf-8-sig", newline="") as handle:
                populations = list(csv.DictReader(handle))
            self.assertEqual(len(populations), 10)
            self.assertEqual(populations[0]["sample"], "original")
            self.assertEqual(float(populations[0]["complex_molecule_count"]), 6.02214076e12)
            self.assertEqual(float(populations[0]["volume_l"]), 10e-6)
            self.assertEqual({row["tube"] for row in populations}, {"trigger_only", "aptamer_trigger"})
            workbook = openpyxl.load_workbook(exports["xlsx"])
            self.assertEqual(workbook["Composition des tubes"].max_row, 11)
            workbook.close()
            with zipfile.ZipFile(exports["zip"]) as archive:
                self.assertIn(Path(exports["tube_csv"]).name, archive.namelist())

    def test_tube_failure_retains_candidate_and_continues_without_zero_metrics(self):
        with tempfile.TemporaryDirectory() as tmp:
            scientific = core_result()
            candidate = scientific["candidates"][0]
            scientific["candidates"].append({**candidate, "rank": 2, "rank_within_length": 2, "extended_trigger": "AACGT"})
            run_request = replace(request(Path(tmp)), parasite_top=5)
            with (
                patch("aptaswitch_studio.extension_runtime.design_trigger_extensions", return_value=scientific),
                patch("aptaswitch_core.thermodynamics.application_analysis.analyze_parasite_tubes", side_effect=[
                    tube_result("ACGT", run_request.tube_conditions), RuntimeError("solver failure"),
                    tube_result("AACGT", run_request.tube_conditions),
                ]) as engine,
                patch("aptaswitch_studio.extension_runtime.export_extension_visuals", return_value={}),
            ):
                result, exports = run_extension_and_export(run_request)
            self.assertEqual(engine.call_count, 3)
            self.assertEqual(result["status"], "completed")
            self.assertEqual(result["parasite_analysis"]["status"], "completed_with_errors")
            self.assertEqual(result["parasite_analysis"]["completed_count"], 1)
            failed = result["candidates"][0]
            self.assertEqual(failed["tube_analysis"]["status"], "failed")
            self.assertNotIn("mixed_trigger_parasite_pct", failed)
            self.assertEqual(failed["final_score"], 0.25)
            self.assertIn("solver failure", Path(exports["log"]).read_text())
            self.assertEqual(result["candidates"][1]["mixed_trigger_parasite_pct"], 10)

    def test_cancelling_tube_analysis_keeps_baseline_and_partial_extensions(self):
        with tempfile.TemporaryDirectory() as tmp:
            run_request = replace(request(Path(tmp)), parasite_top=5)
            with (
                patch("aptaswitch_studio.extension_runtime.design_trigger_extensions", return_value=core_result()),
                patch("aptaswitch_core.thermodynamics.application_analysis.analyze_parasite_tubes", side_effect=[
                    tube_result("ACGT", run_request.tube_conditions), {"status": "cancelled", "tubes": {}},
                ]) as engine,
                patch("aptaswitch_studio.extension_runtime.export_extension_visuals", return_value={}),
            ):
                result, exports = run_extension_and_export(run_request)
            self.assertEqual(engine.call_count, 2)
            self.assertEqual(result["status"], "cancelled")
            self.assertEqual(result["parasite_analysis"]["status"], "cancelled")
            self.assertEqual(result["parasite_analysis"]["completed_count"], 0)
            self.assertEqual(result["reference"]["tube_analysis"]["status"], "completed")
            self.assertEqual(result["candidates"][0]["final_score"], 0.25)
            self.assertNotIn("trigger_only_parasite_pct", result["candidates"][0])
            self.assertTrue(Path(exports["zip"]).is_file())

    def test_already_cancelled_search_does_not_start_tube_analysis(self):
        with tempfile.TemporaryDirectory() as tmp:
            run_request = replace(request(Path(tmp)), parasite_top=5)
            with (
                patch("aptaswitch_studio.extension_runtime.design_trigger_extensions", return_value=core_result("cancelled")),
                patch("aptaswitch_core.thermodynamics.application_analysis.analyze_parasite_tubes") as engine,
                patch("aptaswitch_studio.extension_runtime.export_extension_visuals", return_value={}),
            ):
                result, _ = run_extension_and_export(run_request)
            engine.assert_not_called()
            self.assertEqual(result["parasite_analysis"]["status"], "cancelled")


if __name__ == "__main__":
    unittest.main()
