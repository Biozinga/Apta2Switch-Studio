"""Native exhaustive extensions: enumeration, recovery and exported identity."""
from copy import deepcopy
from dataclasses import asdict, replace
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch
import zipfile

from aptaswitch_core.models import ThermoConditions
from aptaswitch_core.trigger_extension import design_trigger_extensions
from aptaswitch_core.thermodynamics.native_backend import NativeHybridBackend
from aptaswitch_core.thermodynamics.types import EmptyEnsembleError

LIBRARY = os.environ.get("APTASWITCH_HYBRID_NATIVE_LIBRARY")


@unittest.skipUnless(LIBRARY, "explicit native hybrid library required")
class NativeExtensionTests(unittest.TestCase):
    def setUp(self):
        self.backend = NativeHybridBackend(LIBRARY)
        self.thermo = ThermoConditions(temperature_c=37, switch_material="dna")

    def search(self, **options):
        args = dict(aptamer_dna="AC", trigger_dna="GT", target_lengths=[3, 4],
                    thermo=self.thermo, backend=self.backend, top=20, max_ensemble_candidates=20)
        args.update(options)
        return design_trigger_extensions(**args)

    def test_all_ends_and_lengths_keep_exact_dna_sequence_space(self):
        for side in ("5prime", "3prime", "both"):
            with self.subTest(side=side):
                result = self.search(extension_side=side)
                self.assertEqual(result["status"], "completed")
                self.assertEqual(result["mfe_screen"]["total_candidates_screened"], 20)
                self.assertEqual([p["total_candidates_screened"] for p in result["mfe_screen"]["by_length"]], [4, 16])
                self.assertTrue(result["mfe_screen"]["ensemble_ranking_exhaustive"])
                self.assertEqual(result["thermodynamics"], asdict(self.backend.provenance))
                self.assertEqual(result["model"]["material"], self.backend.provenance.model)
                self.assertEqual(result["model"]["aptamer_material"], "dna")

                self.assertEqual({r["trigger_length"] for r in result["candidates"]}, {3, 4})
                core = "GT"
                for row in result["candidates"]:
                    self.assertEqual(row["extended_trigger"], row["extension_5prime"] + core + row["extension_3prime"])
                    self.assertEqual(row["thermodynamics"], result["thermodynamics"])
                    self.assertTrue(0 <= row["target_structure_probability"] <= 1)
                json.dumps(result, allow_nan=False)

    def test_cancel_then_resume_reuses_survivors_without_reenumeration(self):
        with patch.object(self.backend, "minimum", wraps=self.backend.minimum) as minimum:
            saved = self.search(cancel_check=lambda: minimum.call_count >= 4)
        self.assertEqual(saved["status"], "cancelled")
        self.assertEqual(saved["mfe_screen"]["total_candidates_screened"], 3)
        self.assertEqual(saved["candidates"], [])
        with patch.object(self.backend, "minimum", wraps=self.backend.minimum) as minimum:
            resumed = self.search(resume_result=json.loads(json.dumps(saved)))
        self.assertEqual(minimum.call_count, 1 + len(saved["mfe_survivors"]))
        self.assertEqual(resumed["status"], "completed")
        self.assertEqual(resumed["resume"]["method"], "saved_survivors")
        self.assertFalse(resumed["mfe_screen"]["ensemble_ranking_exhaustive"])
        self.assertEqual(resumed["mfe_screen"]["total_candidates_screened"], 3)
        self.assertEqual(saved["status"], "cancelled")

    def test_advance_keeps_partial_label_and_cap_does_not_cap_search(self):
        with patch.object(self.backend, "minimum", wraps=self.backend.minimum) as minimum:
            result = self.search(advance_check=lambda: minimum.call_count >= 4)
        self.assertEqual(result["status"], "completed")
        self.assertTrue(result["mfe_screen"]["stopped_early"])
        self.assertEqual(result["mfe_screen"]["total_candidates_screened"], 3)
        self.assertTrue(result["mfe_screen"]["ensemble_complete"])
        capped = self.search(max_ensemble_candidates=1)
        self.assertEqual(capped["mfe_screen"]["total_candidates_screened"], 20)
        self.assertTrue(capped["mfe_screen"]["exhaustive"])
        self.assertEqual(len(capped["candidates"]), 2)

    def test_resume_rejects_missing_or_changed_model_provenance_before_computation(self):
        saved = self.search()
        for field in (None, "version", "model", "parameter_digest", "ensemble", "standard_state"):
            bad = deepcopy(saved)
            if field is None:
                bad.pop("thermodynamics")
            else:
                bad["thermodynamics"][field] += "-changed"
            with self.subTest(field=field), patch.object(self.backend, "analyze") as analyze:
                with self.assertRaises(ValueError):
                    self.search(resume_result=bad)
                analyze.assert_not_called()

    def test_invalid_conditions_or_backend_output_and_empty_reference_do_not_fall_back(self):
        with patch("aptaswitch_core.thermodynamics.setup.get_backend",
                   side_effect=AssertionError("Explicit backend must be retained")):
            for changes in (dict(trigger_material="rna"), dict(switch_material="rna"), dict(temperature_c=24), dict(sodium_m=.1), dict(magnesium_m=.001)):
                with self.subTest(changes=changes), self.assertRaises(ValueError):
                    self.search(thermo=replace(self.thermo, **changes))
            with self.assertRaises(EmptyEnsembleError):
                self.search(aptamer_dna="AA", trigger_dna="AA")
            original = self.backend.minimum
            def wrong(spec, conditions):
                result = original(spec, conditions)
                return replace(result, context=replace(result.context, conditions=replace(conditions, celsius=38)))
            with patch.object(self.backend, "minimum", side_effect=wrong), self.assertRaisesRegex(ValueError, "changed"):
                self.search()

    def test_clean_process_never_imports_legacy_engine_or_python_vienna(self):
        code = '''
import sys, os
class Reject:
    def find_spec(self, name, path=None, target=None):
        if name == 'RNA':
            raise AssertionError('legacy engine imported: ' + name)
sys.meta_path.insert(0, Reject())
from aptaswitch_core.trigger_extension import design_trigger_extensions
from aptaswitch_core.models import ThermoConditions
from aptaswitch_core.thermodynamics.native_backend import NativeHybridBackend
result = design_trigger_extensions(aptamer_dna='AC', trigger_dna='GT', target_lengths=[3],
    thermo=ThermoConditions(switch_material="dna"), backend=NativeHybridBackend(os.environ['APTASWITCH_HYBRID_NATIVE_LIBRARY']))
assert result['status'] == 'completed'
assert result['mfe_screen']['total_candidates_screened'] == 4
'''
        completed = subprocess.run([sys.executable, "-c", code], text=True, capture_output=True, timeout=30)
        self.assertEqual(completed.returncode, 0, completed.stderr)

    def test_export_and_checkpoint_roundtrip_include_real_model_and_native_tubes(self):
        from aptaswitch_studio.extension_runtime import (
            ExtensionRunRequest, run_extension_and_export, request_from_extension_result,
        )
        from aptaswitch_core.thermodynamics.tube_reporting import TubeAnalysisConditions
        from openpyxl import load_workbook
        with tempfile.TemporaryDirectory() as root:
            request = ExtensionRunRequest("AC", "GT", (3,), self.thermo, "", Path(root),
                top=1, parasite_top=1, tube_conditions=TubeAnalysisConditions(max_complex_size=3))
            with patch("aptaswitch_studio.extension_runtime.analyze_parasite_tubes",
                       side_effect=AssertionError("legacy tube engine called")):
                result, paths = run_extension_and_export(request, backend=self.backend)
            self.assertEqual(result["status"], "completed")
            self.assertEqual(result["parasite_analysis"]["status"], "completed")
            identity = asdict(self.backend.provenance)
            for row in [result["reference"], *result["candidates"]]:
                self.assertEqual(row["tube_analysis"]["thermodynamics"], identity)
                self.assertEqual(row["tube_analysis"]["conditions"]["max_complex_size"], 3)
                for tube in row["tube_analysis"]["tubes"].values():
                    self.assertIn("equilibrium", tube)
            stored = json.loads(Path(paths["json"]).read_text())
            from aptaswitch_studio.run_library import load_run
            reopened = load_run(Path(paths["json"]))
            self.assertEqual(reopened.kind, "extension")
            self.assertEqual(reopened.result["thermodynamics"], identity)
            self.assertEqual(reopened.result["candidates"], result["candidates"])
            self.assertEqual(stored["thermodynamics"], identity)
            self.assertEqual(stored["request"]["thermo"]["material"], identity["model"])
            checkpoint = json.loads((Path(result["output_dir"]) / "mfe_checkpoint.json").read_text())
            self.assertEqual(checkpoint["thermodynamics"], identity)
            self.assertEqual(checkpoint["request"]["thermo"]["material"], identity["model"])
            workbook = load_workbook(paths["xlsx"], read_only=True)
            self.assertTrue(any('thermodynamics' == row[0] and identity["parameter_digest"] in str(row[1])
                                for sheet in workbook for row in sheet.values if len(row) > 1))
            workbook.close()
            self.assertIn(identity["model"], Path(paths["log"]).read_text())
            with zipfile.ZipFile(paths["zip"]) as archive:
                self.assertTrue(any(name.endswith("_results.json") for name in archive.namelist()))
            resume = request_from_extension_result(stored, import_path="", output_dir=Path(root))
            repeated, _ = run_extension_and_export(resume, backend=self.backend, resume_result=stored)
            self.assertEqual(repeated["candidates"], result["candidates"])
