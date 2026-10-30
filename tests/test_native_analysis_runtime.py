"""The application workers use native calculations and retain exact model identity."""
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

from aptaswitch_core.models import ReporterContext, RunRequest, SequenceInput, SwitchArchitecture, ThermoConditions
from aptaswitch_core.structure_ensemble import ensemble_inputs, ensemble_key, load_ensemble
from aptaswitch_core.system_analysis import analyze_system
from aptaswitch_core.tube_analysis import analyze_parasite_tubes, TubeAnalysisConditions
from aptaswitch_core.thermodynamics.setup import get_backend
from aptaswitch_studio.ensemble_runtime import EnsembleJob
from aptaswitch_studio.system_tube_runtime import SystemTubeJob
from aptaswitch_studio import saved_switch_runtime, saved_switch_ensemble_runtime, saved_switch_library
from aptaswitch_studio.extension_runtime import ExtensionRunJob, ExtensionRunRequest, run_extension_and_export
from design_fixtures import design_result_fixture

LIBRARY = os.environ.get("APTASWITCH_HYBRID_NATIVE_LIBRARY")


@unittest.skipUnless(LIBRARY, "Explicit native hybrid library required")
class NativeAnalysisRuntimeTests(unittest.TestCase):
    def setUp(self):
        self.root = Path(self.enterContext(tempfile.TemporaryDirectory())).resolve()
        self.backend = get_backend(LIBRARY)
        self.provenance = asdict(self.backend.provenance)
        self.request = RunRequest(
            sequences=SequenceInput("AACCGGTT", "AACC"),
            architecture=SwitchArchitecture(toehold_length=2, stem_length=2),
            reporter=ReporterContext("test", "Short reporter", "GCU", "GCUUAA"),
            thermo=ThermoConditions(temperature_c=37., sodium_m=.12),
            trials=1, output_dir=self.root,
            thermodynamics_library_path=LIBRARY,
        )
        self.system_thermo = replace(self.request.thermo, switch_material="dna")
        self.result = design_result_fixture(self.request)
        self.candidate = self.result.candidates[0]
        self.snapshot = saved_switch_library.candidate_snapshot(self.result, self.candidate, "Native candidate")
        self.saved = saved_switch_library.save_switch(self.root / "saved", self.snapshot["name"], self.snapshot)
        self.config = dict(temperature_c=37., sodium_m=.12, magnesium_m=0., max_complex_size=2)

    def finish(self, job):
        job.start()
        job._thread.join(20)
        self.assertFalse(job._thread.is_alive())
        state = job.snapshot()
        self.assertEqual(state.state, "done", state.error)
        return state.result

    def test_public_initial_system_routes_ignore_historical_import_path(self):
        mfe = analyze_system("/obsolete/native", "AC", "GT", self.system_thermo)
        tubes = analyze_parasite_tubes("/obsolete/native", "AC", "GT", self.system_thermo,
                                      TubeAnalysisConditions(max_complex_size=2))
        for data in (mfe, tubes):
            self.assertEqual(data["thermodynamics"], self.provenance)

        self.assertAlmostEqual(sum(row["trigger_fraction_pct"] for row in tubes["tubes"]["aptamer_trigger"]["complexes"]), 100., places=6)
        with patch("aptaswitch_core.thermodynamics.setup.get_backend", side_effect=RuntimeError("missing native runtime")):
            with self.assertRaisesRegex(RuntimeError, "missing native"):
                analyze_system("/native", "AC", "GT", self.system_thermo)

    def test_system_and_ensemble_workers_use_native_provenance_for_results_and_cache(self):
        signature = ("AC", "GT", LIBRARY, self.system_thermo,
                     TubeAnalysisConditions(max_complex_size=2), json.dumps(self.provenance, sort_keys=True))
        tubes = self.finish(SystemTubeJob(signature))
        self.assertEqual(tubes["thermodynamics"], self.provenance)
        analysis = self.finish(EnsembleJob(self.result, self.candidate, LIBRARY, num_sample=100))
        inputs = ensemble_inputs(self.request, self.candidate, 100, backend=self.backend)
        self.assertEqual(analysis["inputs"], inputs)
        self.assertEqual(load_ensemble(self.result.output_dir, inputs), analysis)
        old_inputs = ensemble_inputs(self.request, self.candidate, 100)
        self.assertNotEqual(ensemble_key(old_inputs), analysis["key"])
        self.assertIsNone(load_ensemble(self.result.output_dir, old_inputs))

    def test_saved_context_and_families_keep_baseline_identity_and_portable_history(self):
        original = deepcopy(self.saved.snapshot)
        job = saved_switch_runtime.SavedSwitchAnalysisJob(self.saved, self.config, LIBRARY)
        analysis = self.finish(job)
        self.assertEqual(self.saved.snapshot, original)
        self.assertEqual(analysis["thermodynamics"], self.provenance)
        self.assertEqual(analysis["baseline"]["thermodynamics"], self.provenance)
        with patch("aptaswitch_core.thermodynamics.setup.get_backend", side_effect=RuntimeError("uninstalled")):
            source, candidate = saved_switch_ensemble_runtime.load_context_ensemble_source(self.saved, analysis)
        self.assertEqual(candidate.switch_seq, self.candidate.switch_seq + "UAA")
        family_job = saved_switch_ensemble_runtime.create_context_ensemble_job(
            self.saved, analysis, source, candidate, LIBRARY, num_sample=100)
        families = self.finish(family_job)
        self.assertEqual(families["model"], self.provenance)
        self.assertEqual(families["inputs"]["thermo"]["temperature_c"], 37.)
        self.assertEqual(load_ensemble(source.output_dir, family_job.inputs)["key"], families["key"])
        altered = deepcopy(analysis)
        altered["thermodynamics"]["version"] += "-changed"
        Path(analysis["result_path"]).write_text(json.dumps(altered))
        with self.assertRaisesRegex(ValueError, "modèle exact"):
            saved_switch_ensemble_runtime.create_context_ensemble_job(
                self.saved, altered, source, candidate, LIBRARY, num_sample=100)

    def test_extension_worker_native_checkpoint_and_changed_model_resume_rejection(self):
        request = ExtensionRunRequest("AC", "GT", (3,), self.system_thermo, "/obsolete/native", self.root / "extensions",
            top=1, parasite_top=1, tube_conditions=TubeAnalysisConditions(max_complex_size=2),
            thermodynamics_library_path=LIBRARY)
        result = self.finish(ExtensionRunJob(request))
        self.assertEqual(result["thermodynamics"], self.provenance)
        self.assertEqual(result["request"]["engine"], "native-hybrid")
        self.assertEqual(result["request"]["import_path"], "")
        self.assertEqual(result["parasite_analysis"]["status"], "completed")
        checkpoint = json.loads((Path(result["output_dir"]) / "mfe_checkpoint.json").read_text())
        self.assertEqual(checkpoint["thermodynamics"], self.provenance)
        bad = deepcopy(checkpoint)
        bad["thermodynamics"]["parameter_digest"] += "changed"
        destination = self.root / "must-remain-absent"
        with self.assertRaisesRegex(ValueError, "different thermodynamic"):
            run_extension_and_export(replace(request, output_dir=destination), resume_result=bad)
        self.assertFalse(destination.exists())

    def test_clean_process_routes_do_not_import_python_reference_engine(self):
        code = '''
import sys
class Reject:
    def find_spec(self, fullname, path=None, target=None):
        if fullname.split('.')[0] == 'RNA':
            raise AssertionError('Legacy calculation imported: ' + fullname)
sys.meta_path.insert(0, Reject())
from aptaswitch_core.system_analysis import analyze_system
from aptaswitch_core.tube_analysis import analyze_parasite_tubes, TubeAnalysisConditions
from aptaswitch_core.models import ThermoConditions
from aptaswitch_studio.system_tube_runtime import SystemTubeJob
from aptaswitch_studio.ensemble_runtime import EnsembleJob
from aptaswitch_studio.saved_switch_runtime import SavedSwitchAnalysisJob
from aptaswitch_studio.saved_switch_ensemble_runtime import SavedSwitchEnsembleJob
from aptaswitch_studio.extension_runtime import ExtensionRunJob
assert analyze_system('/native', 'AC', 'GT', ThermoConditions(switch_material='dna'))['thermodynamics']
assert analyze_parasite_tubes('/native', 'AC', 'GT', ThermoConditions(switch_material='dna'), TubeAnalysisConditions(max_complex_size=2))['status'] == 'completed'
'''
        completed = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, timeout=20)
        self.assertEqual(completed.returncode, 0, completed.stderr)


if __name__ == "__main__":
    unittest.main()
