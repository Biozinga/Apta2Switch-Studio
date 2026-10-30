"""Global defect and its structural/concentration contributions remain informative."""
from __future__ import annotations

import copy
import csv
import json
from dataclasses import fields, make_dataclass
from pathlib import Path
import tempfile
import unittest

from aptaswitch_core.architecture import default_architecture_for_trigger, target_switch_structures
from aptaswitch_core.defect_metrics import functional_defect_metrics
from aptaswitch_core.exports import CSV_HEADERS, candidate_row, export_run, write_csv, write_xlsx
from aptaswitch_core.models import DesignCandidate, RunRequest, SequenceInput
from aptaswitch_core.scoring import selection_score_100
from aptaswitch_studio.run_library import load_run
from tests.design_fixtures import design_result_fixture
from tests.defect_fixtures import target_matrix


class DefectExportsTests(unittest.TestCase):
    def setUp(self):
        self.root = Path(self.enterContext(tempfile.TemporaryDirectory())).resolve()
        request = RunRequest(SequenceInput("ACGTACGT", "ACGT" * 6), default_architecture_for_trigger(24),
                             output_dir=self.root, trials=1, defect_metric="functional", score_scale="bounded_v1")
        self.result = design_result_fixture(request)
        self.result.output_dir.mkdir(parents=True)
        value = self.result.candidates[0]
        value.score_raw = -27.321
        value.score = selection_score_100(value.score_raw)
        value.defect = .1
        value.defect_structural = .07
        value.defect_concentration = .03
        value.defect_breakdown = {
            "global": .1, "structural": .07, "concentration": .03,
            "tubes": [{"name": "tube_off", "global": .1, "structural": .07, "concentration": .03,
                       "dimensional_defect": 1e-7, "dimensional_structural": 7e-8,
                       "dimensional_concentration": 3e-8, "target_nucleotide_concentration_m": 1e-6,
                       "complexes": [{"name": "switch_off", "structural": 7e-8, "concentration": 3e-8, "total": 1e-7}]}],
        }
        prefix_length = len(request.reporter.sequence_after_start_rna)
        off, on = target_switch_structures(request.architecture, prefix_length)
        value.regional_defects = functional_defect_metrics(request.architecture, prefix_length,
                                                          target_matrix(off), target_matrix(on),
                                                          weights=request.functional_defect_weights)

    def test_json_roundtrip_preserves_informative_defect_and_current_score(self):
        exports = export_run(self.result)
        path = Path(exports["json"])
        original = path.read_bytes()
        loaded = load_run(path).result
        self.assertEqual(loaded.to_dict(), self.result.to_dict())
        value = loaded.candidates[0]
        self.assertEqual((value.defect, value.defect_structural, value.defect_concentration), (.1, .07, .03))
        self.assertAlmostEqual(value.defect_structural + value.defect_concentration, value.defect)
        self.assertEqual(value.score, selection_score_100(-27.321))
        self.assertEqual(loaded.request.defect_metric, "functional")
        self.assertEqual(loaded.request.score_scale, "bounded_v1")
        self.assertEqual(path.read_bytes(), original)

    def test_csv_retains_global_and_both_contributions(self):
        path = write_csv(self.result, self.result.output_dir)
        with path.open(newline="") as handle:
            rows = list(csv.reader(handle))
        self.assertEqual(rows[0], CSV_HEADERS)
        self.assertEqual(len(rows[0]), len(rows[1]))
        row = dict(zip(rows[0], rows[1]))
        for key, expected in {"Defect": "0.1", "Defect_structural": "0.07",
                              "Defect_concentration": "0.03", "Defect_metric_for_score": "functional",
                              "Score_scale": "bounded_v1"}.items():
            self.assertEqual(row[key], expected)

    def test_xlsx_retains_contributions_and_their_definitions(self):
        path = write_xlsx(self.result, self.result.output_dir)
        if path is None:
            self.skipTest("openpyxl is unavailable")
        import openpyxl
        workbook = openpyxl.load_workbook(path)
        self.addCleanup(workbook.close)
        sheet = workbook.worksheets[0]
        self.assertEqual([cell.value for cell in sheet[1]], CSV_HEADERS)
        values = dict(zip(CSV_HEADERS, [cell.value for cell in sheet[2]]))
        self.assertEqual((values["Defect"], values["Defect_structural"], values["Defect_concentration"]), (.1, .07, .03))
        definitions = {row[0]: row[1] for row in workbook["Defect metrics"].iter_rows(values_only=True)}
        self.assertIn("Structural contribution", definitions["Defect_structural"])
        self.assertIn("concentration shortfall", definitions["Defect_concentration"])
        self.assertIn("never mean zero", definitions["Missing metrics"])

    def test_excluded_candidate_has_no_fabricated_global_defect(self):
        value = self.result.candidates[0]
        excluded = DesignCandidate.excluded_for_stop(trial=2, switch_seq=value.switch_seq,
            trigger_seq=value.trigger_seq, first_stop=2, reporter_label="sfGFP")
        row = dict(zip(CSV_HEADERS, candidate_row(2, excluded)))
        for key in ("Defect", "Defect_structural", "Defect_concentration"):
            self.assertIsNone(row[key])

    def test_absent_opening_measurement_exports_blank_cells(self):
        # Retained pre-reload objects keep their original dataclass definition.
        retained_type = make_dataclass("DesignCandidate", [
            (field.name, field.type) for field in fields(DesignCandidate)
            if field.name != "rbs_accessibility"
        ], namespace={"to_dict": DesignCandidate.to_dict})
        values = self.result.candidates[0].to_dict()
        values.pop("rbs_accessibility")
        candidate = retained_type(**values)
        self.result.candidates[0] = candidate
        row = dict(zip(CSV_HEADERS, candidate_row(1, candidate)))
        for key in ("ON_RBS_AUG_joint_unpaired_probability",
                    "ON_RBS_AUG_opening_kcal_mol", "ON_RBS_AUG_opening_impossible"):
            self.assertIsNone(row[key])
        self.assertEqual(row["Score_selection"], candidate.score)
        exports = export_run(self.result)
        with Path(exports["csv"]).open(newline="") as handle:
            exported = next(csv.DictReader(handle))
        self.assertEqual(exported["ON_RBS_AUG_joint_unpaired_probability"], "")
        loaded = load_run(Path(exports["json"])).result.candidates[0]
        self.assertEqual(loaded.rbs_accessibility, {})
        self.assertEqual(loaded.score, candidate.score)
        self.assertFalse(hasattr(candidate, "rbs_accessibility"))

    def test_zero_defect_can_keep_an_unknown_normalization_denominator(self):
        payload = self.result.to_dict()
        candidate = payload["results"][0]
        candidate.update(defect=0., defect_structural=0., defect_concentration=0.)
        breakdown = candidate["defect_breakdown"]
        for key in ("global", "structural", "concentration"):
            breakdown[key] = 0.
            breakdown["tubes"][0][key] = 0.
        tube = breakdown["tubes"][0]
        tube.update(dimensional_defect=0., dimensional_structural=0., dimensional_concentration=0.,
                    target_nucleotide_concentration_m=None)
        tube["complexes"][0].update(structural=0., concentration=0., total=0.)
        path = self.result.output_dir / "zero_results.json"
        path.write_text(json.dumps(payload))
        loaded = load_run(path).result.candidates[0]
        self.assertEqual(loaded.defect, 0.)
        self.assertIsNone(loaded.defect_breakdown["tubes"][0]["target_nucleotide_concentration_m"])

    def test_corrupted_global_defect_and_concentration_measurements_are_rejected(self):
        base = self.result.to_dict()
        mutations = [
            lambda p: p["results"][0].update(defect_structural=float("nan")),
            lambda p: p["results"][0].update(defect_concentration=1.01),
            lambda p: p["results"][0].update(defect_structural=True),
            lambda p: p["results"][0].update(defect_breakdown=[]),
            lambda p: p["results"][0]["defect_breakdown"].update({"global": float("inf")}),
            lambda p: p["results"][0]["defect_breakdown"].update(structural=-.01),
            lambda p: p["results"][0]["defect_breakdown"].update(tubes={}),
            lambda p: p["results"][0]["defect_breakdown"]["tubes"][0].update(dimensional_defect=float("inf")),
            lambda p: p["results"][0]["defect_breakdown"]["tubes"][0].update(dimensional_structural=-1),
            lambda p: p["results"][0]["defect_breakdown"]["tubes"][0].update(target_nucleotide_concentration_m=True),
            lambda p: p["results"][0]["defect_breakdown"]["tubes"][0].update(complexes={}),
            lambda p: p["results"][0]["defect_breakdown"]["tubes"][0]["complexes"][0].update(concentration=-1),
        ]
        path = self.result.output_dir / "invalid_results.json"
        for index, mutate in enumerate(mutations):
            with self.subTest(index=index):
                payload = copy.deepcopy(base)
                mutate(payload)
                path.write_text(json.dumps(payload))
                with self.assertRaises(ValueError):
                    load_run(path)


if __name__ == "__main__":
    unittest.main()
