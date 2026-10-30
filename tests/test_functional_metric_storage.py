"""Functional metrics retain their exact definitions and weights in archives."""
from __future__ import annotations

import copy
import csv
from dataclasses import replace
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from aptaswitch_core.architecture import default_architecture_for_trigger, target_switch_structures
from aptaswitch_core.defect_metrics import functional_defect_metrics
from aptaswitch_core.exports import CSV_HEADERS, candidate_row, export_run
from aptaswitch_core.models import DesignCandidate, FunctionalDefectWeights, RunRequest, SequenceInput
from aptaswitch_studio.run_library import load_run
from aptaswitch_core.scoring import selection_score_100
from tests.design_fixtures import design_result_fixture, synthetic_native_metadata
from tests.defect_fixtures import disrupt_pairs, disrupt_unpaired, target_matrix


class FunctionalMetricStorageTests(unittest.TestCase):
    def setUp(self):
        directory = tempfile.TemporaryDirectory()
        self.addCleanup(directory.cleanup)
        self.root = Path(directory.name).resolve()
        request = RunRequest(SequenceInput("ACGTACGT", "ACGT" * 7 + "AC"),
                             default_architecture_for_trigger(30), output_dir=self.root, trials=2,
                             defect_metric="functional")
        self.result = design_result_fixture(request)
        self.result.output_dir.mkdir(parents=True)
        for candidate, raw in zip(self.result.candidates, [-20.5, 100]):
            candidate.score_raw = raw
            candidate.score = selection_score_100(raw)
        self.remeasure()

    def remeasure(self):
        request = self.result.request
        # These fixtures change the request directly, without rerunning a real
        # design. Keep their recorded optimization weights in the same context.
        self.result.optimization = synthetic_native_metadata(request)["optimization"]
        reporter_length = len(request.reporter.sequence_after_start_rna)
        off_target, on_target = target_switch_structures(request.architecture, reporter_length)
        off, on = target_matrix(off_target), target_matrix(on_target)
        regions = functional_defect_metrics(request.architecture, reporter_length, off, on)["regions"]
        disrupt_pairs(off, regions["off_stem"]["indices_switch"], .2)
        disrupt_unpaired(off, regions["off_toehold"]["indices_switch"], .4)
        offset = request.architecture.trigger_region_length
        disrupt_unpaired(on, [offset + index for index in regions["on_rbs_aug"]["indices_switch"][:2]], .6)
        metrics = functional_defect_metrics(request.architecture, reporter_length, off, on,
                                           weights=request.functional_defect_weights)
        for candidate in self.result.candidates:
            candidate.regional_defects = copy.deepcopy(metrics)

    def write_payload(self, payload):
        path = self.result.output_dir / "functional_results.json"
        path.write_text(json.dumps(payload))
        return path

    def test_roundtrip_preserves_conditions_metrics_order_scores_and_source(self):
        exports = export_run(self.result)
        path = Path(exports["json"])
        original = path.read_bytes()
        with patch("aptaswitch_core.thermodynamics.setup.get_backend", side_effect=AssertionError("No recalculation")), \
             patch("aptaswitch_core.engines.design_candidates", side_effect=AssertionError("No recalculation")):
            loaded = load_run(path).result
        self.assertEqual(loaded.to_dict(), self.result.to_dict())
        self.assertEqual([value.score for value in loaded.candidates], [selection_score_100(raw) for raw in (-20.5, 100)])
        self.assertEqual(loaded.request.defect_metric, "functional")
        self.assertEqual(loaded.candidates[0].regional_defects["scheme"], "functional_v2")
        self.assertEqual(path.read_bytes(), original)

    def test_exports_contain_only_current_functional_fields(self):
        exports = export_run(self.result)
        with Path(exports["csv"]).open(newline="") as handle:
            row = next(csv.DictReader(handle))
        self.assertNotIn("Defect_ON_binding", CSV_HEADERS)
        self.assertFalse(any(name.startswith("Regional_weight_") for name in CSV_HEADERS))
        self.assertEqual(row["Defect_metric_for_score"], "functional")
        self.assertEqual(row["Defect_metric_scheme"], "functional_v2")
        self.assertTrue(all(row[name] == "" for name in CSV_HEADERS if name.startswith("Regional_weight_")))
        self.assertAlmostEqual(float(row["Defect_OFF_stem_before_AUG"]), .2)
        self.assertAlmostEqual(float(row["Defect_OFF_stem_after_AUG"]), .2)
        self.assertGreater(float(row["Defect_ON_RBS_AUG"]), float(row["Defect_ON_RBS_linker"]))
        self.assertEqual(float(row["Functional_weight_ON_RBS_linker"]), self.result.request.functional_defect_weights.on_rbs_linker)
        journal = Path(exports["log"]).read_text()
        self.assertIn("Functional_defect_weights:", journal)
        self.assertIn("functional_v2", journal)
        self.assertNotIn("Regional_defect_weights:", journal)

    def test_workbook_documents_saved_weights_and_on_overlap(self):
        self.result.request = replace(self.result.request, functional_defect_weights=FunctionalDefectWeights(0, 2, 3, 1))
        self.remeasure()
        exports = export_run(self.result)
        if "xlsx" not in exports:
            self.skipTest("openpyxl is unavailable")
        import openpyxl
        workbook = openpyxl.load_workbook(exports["xlsx"])
        self.addCleanup(workbook.close)
        definitions = {row[0]: (row[1], row[2]) for row in workbook["Defect metrics"].iter_rows(values_only=True)}
        self.assertIn("/ sum(functional_weights)", definitions["Defect_regional_balanced"][0])
        self.assertIn("on_rbs_aug and on_rbs_linker", definitions["Defect_ON"][0])
        self.assertEqual(definitions["Functional weight: on_rbs_aug"][0], 3)
        self.assertIn("50.000000%", definitions["Functional weight: on_rbs_aug"][1])
        self.assertIn("contained", definitions["Overlapping ON regions"][0])
        self.assertNotIn("Regional weight: on_binding", definitions)

    def test_disabled_state_keeps_local_measurements_and_absent_aggregate(self):
        self.result.request = replace(self.result.request, functional_defect_weights=FunctionalDefectWeights(0, 0, 3, 0))
        self.remeasure()
        exports = export_run(self.result)
        loaded = load_run(Path(exports["json"])).result
        metrics = loaded.candidates[0].regional_defects
        self.assertIsNone(metrics["off"])
        self.assertGreater(metrics["off_stem"], 0)
        self.assertEqual(metrics["balanced"], metrics["on_rbs_aug"])
        with Path(exports["csv"]).open(newline="") as handle:
            row = next(csv.DictReader(handle))
        self.assertEqual(row["Defect_OFF"], "")
        self.assertEqual(row["Functional_weight_OFF_stem"], "0")
        self.assertNotEqual(row["Defect_OFF_stem"], "")

    def test_invalid_schema_weights_and_stem_diagnostics_are_rejected(self):
        base = self.result.to_dict()
        mutations = [
            lambda p: p["request"].pop("defect_metric"),
            lambda p: p["request"].pop("functional_defect_weights"),
            lambda p: p["request"].update(defect_metric="regional"),
            lambda p: p["request"]["functional_defect_weights"].update(on_binding=1),
            lambda p: p["request"]["functional_defect_weights"].update(on_rbs_aug=2),
            lambda p: p["results"][0]["regional_defects"].pop("scheme"),
            lambda p: p["results"][0]["regional_defects"].update(scheme="future_unknown"),
            lambda p: p["results"][0]["regional_defects"]["weights"].pop("on_rbs_aug"),
            lambda p: p["results"][0]["regional_defects"]["weights"].update(on_binding=1),
            lambda p: p["results"][0]["regional_defects"]["normalized_weights"].update(on_rbs_aug=.99),
            lambda p: p["results"][0]["regional_defects"].update(off=None),
            lambda p: p["results"][0]["regional_defects"].update(on=None),
            lambda p: p["results"][0]["regional_defects"].update(on_binding=.3),
            lambda p: p["results"][0]["regional_defects"]["regions"].pop("on_rbs_aug"),
            lambda p: p["results"][0]["regional_defects"].pop("on_rbs_aug"),
            lambda p: p["results"][0]["regional_defects"]["off_stem_parts"]["before_aug"].update(pair_count=999),
            lambda p: p["results"][0]["regional_defects"]["off_stem_parts"]["after_aug"].update(indices_switch=[99999]),
            lambda p: p["results"][0]["regional_defects"]["off_stem_parts"].pop("after_aug"),
        ]
        for index, mutate in enumerate(mutations):
            payload = copy.deepcopy(base)
            mutate(payload)
            with self.subTest(index=index), self.assertRaises(ValueError):
                load_run(self.write_payload(payload))


    def test_excluded_candidate_has_saved_policy_but_no_fabricated_measurements(self):
        original = self.result.candidates[0]
        candidate = DesignCandidate.excluded_for_stop(trial=3, switch_seq=original.switch_seq,
            trigger_seq=original.trigger_seq, first_stop=2, reporter_label="sfGFP")
        row = dict(zip(CSV_HEADERS, candidate_row(1, candidate)))
        self.assertIsNone(row["Defect_regional_balanced"])
        self.assertIsNone(row["Defect_ON_RBS_AUG"])
        self.assertIsNone(row["Defect_OFF_stem_before_AUG"])
        self.assertEqual(row["Functional_weight_ON_RBS_linker"], FunctionalDefectWeights().on_rbs_linker)


if __name__ == "__main__":
    unittest.main()
