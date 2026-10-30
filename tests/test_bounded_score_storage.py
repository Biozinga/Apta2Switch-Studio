"""Saved score scales and trigger accessibility remain explicit and reproducible."""
from __future__ import annotations

import copy
import csv
from dataclasses import replace
import json
from pathlib import Path
import tempfile
import unittest

from aptaswitch_core.architecture import default_architecture_for_trigger
from aptaswitch_core.exports import export_run
from aptaswitch_core.models import RunRequest, SequenceInput
from aptaswitch_core.scoring import selection_score_100
from aptaswitch_studio.run_library import load_run
from tests.design_fixtures import design_result_fixture


class BoundedScoreStorageTests(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.root = Path(self.directory.name).resolve()
        request = RunRequest(SequenceInput("ACGTACGT", "ACGT" * 6),
                             default_architecture_for_trigger(24), output_dir=self.root, trials=1)
        self.result = design_result_fixture(request)
        self.result.output_dir.mkdir(parents=True)
        self.candidate = self.result.candidates[0]
        self.candidate.score_raw = 199.0
        self.candidate.score = selection_score_100(199.0)
        length = len(self.candidate.trigger_seq)
        self.candidate.trigger_linearity = {
            "normalized_defect": .2, "unpaired_fraction": .8,
            "per_base_unpaired": [.8] * length, "structure_mfe": "." * length,
            "nucleotide_count": length, "expected_paired_nt": length * .2,
            "model": {"material": "experimental-Sugimoto95-Vienna04-loopmean", "trigger_material": "dna",
                      "temperature_c": 37, "sodium_m": .12, "magnesium_m": 0.},
        }

    def write_payload(self, payload):
        path = self.result.output_dir / "test_fixture_results.json"
        path.write_text(json.dumps(payload))
        return path

    def test_bounded_json_roundtrip_retains_raw_score_and_full_trigger_diagnostics(self):
        exported = export_run(self.result)
        path = Path(exported["json"])
        original_bytes = path.read_bytes()
        loaded = load_run(path).result
        self.assertEqual(loaded.to_dict(), self.result.to_dict())
        self.assertEqual(loaded.request.score_scale, "bounded_v1")
        self.assertEqual(loaded.candidates[0].score_raw, 199.0)
        self.assertEqual(loaded.request.scoring.trigger_linearity, 3.0)
        self.assertEqual(path.read_bytes(), original_bytes)


    def test_csv_xlsx_svg_and_log_describe_the_same_scale_and_penalty(self):
        self.result.request = replace(self.result.request,
                                      scoring=replace(self.result.request.scoring, trigger_linearity=2.5))
        exported = export_run(self.result)
        with Path(exported["csv"]).open(newline="") as handle:
            row = next(csv.DictReader(handle))
        self.assertEqual(float(row["Score_selection"]), self.candidate.score)
        self.assertEqual(float(row["Score_raw"]), 199)
        self.assertEqual(float(row["Trigger_linearity_defect"]), .2)
        self.assertEqual(float(row["Trigger_unpaired_fraction"]), .8)
        self.assertEqual(float(row["Weight_trigger_linearity"]), 2.5)
        self.assertIn(f"{self.candidate.score:.4f}", Path(exported["svg"]).read_text())
        self.assertNotIn("199.0000", Path(exported["svg"]).read_text())
        log = Path(exported["log"]).read_text()
        self.assertIn("atan(Score_raw/100)", log)
        self.assertIn("Weight_trigger_linearity: 2.5", log)
        self.assertIn("last term only for inverse extension", log)
        if "xlsx" in exported:
            import openpyxl
            workbook = openpyxl.load_workbook(exported["xlsx"])
            self.addCleanup(workbook.close)
            sheet = workbook.worksheets[0]
            values = dict(zip((cell.value for cell in sheet[1]), (cell.value for cell in sheet[2])))
            self.assertAlmostEqual(values["Score_selection"], self.candidate.score)
            definitions = {row[0]: row[1:] for row in workbook["Defect metrics"].iter_rows(values_only=True)}
            self.assertIn("atan", definitions["Score_selection"][0])
            self.assertIn("only", definitions["Score_raw"][1])
            self.assertIn("not the probability", definitions["Trigger_unpaired_fraction"][1])

    def test_invalid_score_scale_and_trigger_measurements_are_rejected(self):
        base = self.result.to_dict()
        mutations = [
            lambda p: p["request"].pop("score_scale"),
            lambda p: p["request"].update(score_scale="unknown"),
            lambda p: p["request"].update(score_scale=[]),
            lambda p: p["request"]["scoring"].update(trigger_linearity=-1),
            lambda p: p["request"]["scoring"].update(trigger_linearity=True),
            lambda p: p["results"][0].update(score=199),
            lambda p: p["results"][0].update(score=-1),
            lambda p: p["results"][0].update(score_raw=float("inf")),
            lambda p: p["results"][0].update(score_raw=True),
            lambda p: p["results"][0].update(score_raw=0),
            lambda p: p["results"][0].update(trigger_linearity=[]),
            lambda p: p["results"][0]["trigger_linearity"].update(normalized_defect=1.1),
            lambda p: p["results"][0]["trigger_linearity"].update(normalized_defect=.5),
            lambda p: p["results"][0]["trigger_linearity"].update(unpaired_fraction=.4),
            lambda p: p["results"][0]["trigger_linearity"].update(per_base_unpaired=[.8]),
            lambda p: p["results"][0]["trigger_linearity"].update(per_base_unpaired=[float("nan")] * 24),
            lambda p: p["results"][0]["trigger_linearity"].update(structure_mfe="(" * 24),
            lambda p: p["results"][0]["trigger_linearity"].update(structure_mfe="." * 23 + "+"),
            lambda p: p["results"][0]["trigger_linearity"].update(nucleotide_count=23),
            lambda p: p["results"][0]["trigger_linearity"].update(expected_paired_nt=100),
            lambda p: p["results"][0]["trigger_linearity"].update(model=[]),
        ]
        for index, mutate in enumerate(mutations):
            with self.subTest(index=index):
                payload = copy.deepcopy(base)
                mutate(payload)
                with self.assertRaises(ValueError):
                    load_run(self.write_payload(payload))


if __name__ == "__main__":
    unittest.main()
