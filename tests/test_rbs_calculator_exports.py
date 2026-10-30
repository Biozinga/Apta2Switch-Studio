"""Exports retain recorded OSTIR values and their limits without recalculation."""
from dataclasses import asdict, replace
import csv
import json
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

from aptaswitch_core.architecture import default_architecture_for_trigger
from aptaswitch_core.exports import CSV_HEADERS, candidate_row, export_run, write_csv
from aptaswitch_core.models import RunRequest, SequenceInput
from tests.design_fixtures import design_result_fixture


OSTIR_COLUMNS = ("OSTIR_OFF_reference_rate", "OSTIR_ON_constrained_rate",
                 "OSTIR_log10_ON_OFF", "OSTIR_exploratory_merit")


class RBSCalculatorExportsTests(unittest.TestCase):
    def setUp(self):
        self.root = Path(self.enterContext(tempfile.TemporaryDirectory())).resolve()
        request = RunRequest(SequenceInput("ACGTACGT", "ACGT" * 6),
                             default_architecture_for_trigger(24), output_dir=self.root, trials=1)
        self.result = design_result_fixture(request)
        self.result.request = replace(self.result.request,
            scoring=replace(self.result.request.scoring, rbs_calculator=.75))
        self.candidate = self.result.candidates[0]
        self.candidate.score_raw = 0.
        # Arbitrary rendering fixtures; no thermodynamic prediction is made.
        self.measurement = {
            "method": "synthetic-export-fixture",
            "off_rate": 12.5, "on_rate": 125., "log10_on_off": 1., "merit": .0625,
            "reference": {"ostir_version": "1.1.3", "vienna_version": "2.7.2",
                          "celsius": 37., "sodium_m": 1.021,
                          "anti_shine_dalgarno": "ACCUCCUUA",
                          "parameter_sha256": "synthetic-parameter-fixture"},
            "window_start": 12, "window_end": 82, "window_aug_start": 35,
            "window_sequence": "A" * 35 + "AUG" + "C" * 32,
            "window_occupied_indices": [0, 1, 2], "on_constraints": "xxx" + "." * 67,
        }
        self.candidate.rbs_calculator = self.measurement
        self.enterContext(patch("aptaswitch_core.rbs_calculator.predict_rbs_switch",
                                side_effect=AssertionError("Export must not recalculate OSTIR")))

    def test_csv_json_and_log_retain_the_measured_values_weight_and_raw_score(self):
        before = self.candidate.to_dict()
        exported = export_run(self.result)
        with Path(exported["csv"]).open(newline="") as handle:
            row = next(csv.DictReader(handle))
        self.assertEqual(tuple(float(row[key]) for key in OSTIR_COLUMNS), (12.5, 125., 1., .0625))
        self.assertEqual(float(row["Weight_rbs_calculator"]), .75)
        self.assertEqual(float(row["Score_raw"]), self.candidate.score_raw)
        payload = json.loads(Path(exported["json"]).read_text())
        self.assertEqual(payload["results"][0]["rbs_calculator"], self.measurement)
        self.assertEqual(payload["request"]["scoring"]["rbs_calculator"], .75)
        log = Path(exported["log"]).read_text()
        self.assertIn("100*w_rbs_calculator*OSTIR_exploratory_merit", log)
        self.assertIn("Weight_rbs_calculator: 0.75", log)
        self.assertEqual(before, self.candidate.to_dict())

    def test_xlsx_defines_the_rna_proxy_and_preserves_recorded_reference_parameters(self):
        exported = export_run(self.result)
        if "xlsx" not in exported:
            self.skipTest("openpyxl is unavailable")
        import openpyxl
        workbook = openpyxl.load_workbook(exported["xlsx"])
        self.addCleanup(workbook.close)
        sheet = workbook.worksheets[0]
        row = dict(zip(CSV_HEADERS, (cell.value for cell in sheet[2])))
        self.assertEqual(tuple(row[key] for key in OSTIR_COLUMNS), (12.5, 125., 1., .0625))
        self.assertEqual(row["Weight_rbs_calculator"], .75)
        definitions = {row[0]: row[1:] for row in workbook["Defect metrics"].iter_rows(values_only=True)}
        self.assertIn("100*w_rbs_calculator*OSTIR_exploratory_merit", definitions["Score_raw"][0])
        self.assertIn("not OFF translation leakage", definitions["Leak_pct"][1])
        self.assertIn("no DNA binding energy", definitions["OSTIR_ON_constrained_rate"][1])
        self.assertIn("not an experimentally validated induction ratio", definitions["OSTIR_log10_ON_OFF"][1])
        self.assertIn("not a validated ON target", definitions["OSTIR_exploratory_merit"][1])
        reference_text = " ".join(definitions["OSTIR reference model"])
        for fragment in ("OSTIR 1.1.3", "ViennaRNA 2.7.2", "E. coli", "37 C", "1.021 M", "+/-35 nt"):
            self.assertIn(fragment, reference_text)
        references = list(workbook["RBS Calculator model"].iter_rows(values_only=True))
        ref = next(row for row in references if row[1] == "reference")
        self.assertEqual(ref[0], self.candidate.trial)
        self.assertEqual(json.loads(ref[2]), self.measurement["reference"])

    def test_unmeasured_and_retained_candidates_have_blank_measurements_not_zero(self):
        self.candidate.rbs_calculator = {}
        retained = SimpleNamespace(**{key: value for key, value in asdict(self.candidate).items()
                                      if key != "rbs_calculator"})
        for candidate in (self.candidate, retained):
            with self.subTest(retained=candidate is retained):
                values = candidate_row(1, candidate)
                self.assertEqual(len(values), len(CSV_HEADERS))
                row = dict(zip(CSV_HEADERS, values))
                for column in (*OSTIR_COLUMNS, "Weight_rbs_calculator"):
                    self.assertIsNone(row[column])
                self.assertEqual(row["Score_raw"], self.candidate.score_raw)
        exported = export_run(self.result)
        with Path(exported["csv"]).open(newline="") as handle:
            row = next(csv.DictReader(handle))
        self.assertTrue(all(row[key] == "" for key in OSTIR_COLUMNS))
        self.assertEqual(float(row["Weight_rbs_calculator"]), .75)
        if "xlsx" in exported:
            import openpyxl
            workbook = openpyxl.load_workbook(exported["xlsx"])
            self.addCleanup(workbook.close)
            self.assertNotIn("RBS Calculator model", workbook.sheetnames)

    def test_unrecorded_weight_is_exported_blank_without_inventing_a_default(self):
        old_weights = SimpleNamespace(**{key: value for key, value in asdict(self.result.request.scoring).items()
                                        if key != "rbs_calculator"})
        self.result.request = replace(self.result.request, scoring=old_weights)
        self.candidate.rbs_calculator = {}
        self.result.output_dir.mkdir(parents=True, exist_ok=True)
        path = write_csv(self.result, self.result.output_dir)
        with path.open(newline="") as handle:
            row = next(csv.DictReader(handle))
        self.assertEqual(row["Weight_rbs_calculator"], "")
        self.assertTrue(all(row[key] == "" for key in OSTIR_COLUMNS))
