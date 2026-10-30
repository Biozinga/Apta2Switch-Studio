"""Supplied sequences use the existing native measurements, without optimization."""
from copy import deepcopy
from dataclasses import replace
import json
import os
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from aptaswitch_core.design_checkpoint import request_digest
from aptaswitch_core.engines import run_design
from aptaswitch_core.models import ReporterContext, RunRequest, SequenceInput, SwitchArchitecture, ThermoConditions
from aptaswitch_core.preview import build_switch_model
from aptaswitch_core.thermodynamics.design_constraints import DesignAssignment, build_design_problem
from aptaswitch_core.thermodynamics.types import Strand
from aptaswitch_core.workflow import run_and_export
from aptaswitch_studio.run_library import _load_design, load_run
from aptaswitch_studio.design_job_registry import completed_trials
from aptaswitch_studio.saved_switch_design import load_saved_design
from aptaswitch_studio.saved_switch_library import candidate_snapshot, save_switch
from tests.test_thermodynamics_design_engine import AnalyticalBackend


LIBRARY = os.environ.get("APTASWITCH_HYBRID_NATIVE_LIBRARY", "")


def request(folder, *, rbs=False, stop=False):
    options = RunRequest(
        SequenceInput("GGTT", "AACC"),
        SwitchArchitecture(2, 2, upper_stem_length=1, upper_stem2_length=0,
                           rbs_loop="AGGA", t7_leader="", rbs_prefix_length=0),
        reporter=ReporterContext(key="supplied", label="Supplied reporter",
                                 sequence_after_start_rna="GCU", full_sequence_after_start_rna="GCU"),
        thermo=ThermoConditions(temperature_c=37, validation_max_complex_size=2),
        trials=1, output_dir=Path(folder) / "designs", thermodynamics_library_path=LIBRARY or None,
    )
    sequence = build_switch_model(options.architecture, options.sequences.trigger_dna,
                                  reporter_sequence="GCU").switch_sequence.replace("X", "A")
    if stop:
        from aptaswitch_core.architecture import architecture_layout
        start = architecture_layout(options.architecture)["aug_start"] + 3
        sequence = sequence[:start] + "UAA" + sequence[start + 3:]
    return replace(options, fixed_switch_sequence=sequence,
                   scoring=replace(options.scoring, rbs_calculator=options.scoring.rbs_calculator if rbs else 0.)).normalized()


def analytical_backend(options):
    backend = AnalyticalBackend(options, bad_on=True)
    backend.provenance = replace(backend.provenance, model=options.thermo.material)
    return backend


class NativeReadingTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.folder = Path(temporary.name).resolve()

    def test_fixed_patterns_have_no_pair_or_stop_constraints_and_never_change_inputs(self):
        options = request(self.folder, stop=True)
        self.assertFalse(options.exclude_stop_candidates)
        self.assertFalse(options.exclude_non_linear_candidates)
        # A mismatched left recognition arm is measured, never repaired.
        sequence = "AA" + options.fixed_switch_sequence[2:]
        options = replace(options, fixed_switch_sequence=sequence).normalized()
        problem = build_design_problem(options)
        self.assertEqual(problem.libraries, ())
        self.assertEqual(tuple(strand.pattern for strand in problem.strands), ("AACC", sequence))
        self.assertTrue(all(len(values) == 1 for values in problem._domains.values()))
        assignment = DesignAssignment((Strand("AACC", "dna", "Trigger"), Strand(sequence, "rna", "Switch")))
        problem.validate(assignment)
        self.assertTrue(problem.reading_frame_report(assignment)["has_premature_stop"])
        with self.assertRaises(ValueError):
            problem.validate(DesignAssignment((Strand("AACC", "dna", "Trigger"),
                                               Strand("C" + sequence[1:], "rna", "Switch"))))

    def test_ordinary_request_serialization_and_digest_do_not_gain_a_new_null_field(self):
        options = replace(request(self.folder), fixed_switch_sequence=None)
        self.assertNotIn("fixed_switch_sequence", options.to_dict())
        self.assertEqual(request_digest(options), request_digest(options.normalized()))
        fixed = request(self.folder)
        self.assertEqual(fixed.to_dict()["fixed_switch_sequence"], fixed.fixed_switch_sequence)

    def test_complete_native_candidate_metrics_survive_stop_without_optimizer_records(self):
        options = request(self.folder, stop=True)
        backend = analytical_backend(options)
        with (patch("aptaswitch_core.thermodynamics.setup.backend_for_request", return_value=backend),
              patch("aptaswitch_core.engines.design_candidates", side_effect=AssertionError("No sequence design")),
              patch("aptaswitch_core.thermodynamics.optimization.optimize", side_effect=AssertionError("No optimizer"))):
            result, exports = run_and_export(options)
        candidate, = result.candidates
        self.assertEqual(result.status, "completed")
        self.assertEqual(result.optimization["algorithm"], "fixed_sequence_evaluation")
        self.assertFalse(result.optimization["sequence_optimization"])
        self.assertEqual(result.optimization["evaluations"], 1)
        self.assertEqual(completed_trials(result), 1)
        self.assertEqual(result.optimization["trials"], [])
        self.assertTrue(candidate.has_stop_codon)
        self.assertFalse(candidate.excluded)
        self.assertEqual(candidate.switch_seq, options.fixed_switch_sequence)
        self.assertEqual(candidate.trigger_seq, options.sequences.trigger_dna)
        self.assertTrue(candidate.rbs_accessibility)
        self.assertIsNotNone(candidate.score)
        self.assertIsNotNone(candidate.on_yield_pct)
        self.assertEqual(load_run(exports["json"]).result.to_dict(), result.to_dict())

    def test_cancel_or_finish_before_analysis_never_calls_the_candidate_calculator(self):
        options = request(self.folder)
        backend = analytical_backend(options)
        with (patch("aptaswitch_core.thermodynamics.setup.backend_for_request", return_value=backend),
              patch("aptaswitch_core.thermodynamics.design_engine.evaluate_candidate", side_effect=AssertionError("No requested analysis"))):
            cancelled, exports = run_and_export(options, cancel_check=lambda: True)
            paused = run_design(options, control_check=lambda: "finish")
        self.assertEqual(cancelled.status, "cancelled")
        self.assertEqual(cancelled.candidates, [])
        self.assertEqual(cancelled.optimization["evaluations"], 0)
        self.assertEqual(completed_trials(cancelled), 0)
        self.assertEqual(load_run(exports["json"]).result.status, "cancelled")
        self.assertEqual(paused.status, "paused")
        self.assertEqual(paused.candidates, [])
        self.assertEqual(backend.calls, [])

    def test_pause_after_evaluation_resumes_exports_without_repeating_measurements(self):
        options = request(self.folder)
        controls = iter((None, "pause"))
        with patch("aptaswitch_core.thermodynamics.setup.backend_for_request", return_value=analytical_backend(options)):
            paused, exports = run_and_export(options, control_check=lambda: next(controls))
            self.assertEqual(paused.status, "paused")
            self.assertEqual(paused.optimization["evaluations"], 1)
            self.assertEqual(paused.optimization["checkpoint"]["phase"], "exports")
            loaded = load_run(exports["checkpoint"]).result
            with patch("aptaswitch_core.thermodynamics.design_engine.evaluate_candidate", side_effect=AssertionError("Never repeat a completed evaluation")):
                completed, exports = run_and_export(options, resume_result=loaded)
        self.assertEqual(completed.status, "completed")
        self.assertEqual(completed.candidates[0].to_dict(), paused.candidates[0].to_dict())
        self.assertEqual(completed.optimization["checkpoint"]["phase"], "done")
        self.assertEqual(load_run(exports["json"]).result.to_dict(), completed.to_dict())

    def test_loader_rejects_modified_fixed_sequences_and_fake_optimizer_metadata(self):
        options = request(self.folder)
        with patch("aptaswitch_core.thermodynamics.setup.backend_for_request", return_value=analytical_backend(options)):
            result, _ = run_and_export(options)
        payload = result.to_dict()
        mutations = (
            lambda data: data["optimization"].update(sequence_optimization=True),
            lambda data: data["optimization"].update(evaluations=0),
            lambda data: data["optimization"].update(trials=[{"trial": 1}]),
            lambda data: data["results"][0].update(trigger_seq="AAAA"),
            lambda data: data["results"][0].update(switch_seq="A" + options.fixed_switch_sequence[1:]),
            lambda data: data["request"].update(exclude_stop_candidates=True),
        )
        for index, mutate in enumerate(mutations):
            changed = deepcopy(payload)
            mutate(changed)
            with self.subTest(mutation=index), self.assertRaises(ValueError):
                _load_design(changed, result.output_dir / "result.json")


@unittest.skipUnless(LIBRARY, "Explicit compiled hybrid library required")
class RealNativeReadingTests(unittest.TestCase):
    def test_real_native_ostir_export_and_saved_snapshot_never_optimize(self):
        with tempfile.TemporaryDirectory() as folder:
            options = request(Path(folder).resolve(), rbs=True, stop=True)
            with (patch("aptaswitch_core.engines.design_candidates", side_effect=AssertionError("No design")),
                  patch("aptaswitch_core.thermodynamics.optimization.optimize", side_effect=AssertionError("No optimizer"))):
                result, exports = run_and_export(options)
            candidate, = result.candidates
            self.assertEqual(result.status, "completed")
            self.assertEqual(result.thermodynamics["engine"], "native-hybrid")
            self.assertEqual(result.thermodynamics["model"], options.thermo.material)
            self.assertEqual(candidate.switch_seq, options.fixed_switch_sequence)
            self.assertEqual(candidate.trigger_seq, options.sequences.trigger_dna)
            self.assertEqual(result.request.sequences.aptamer_dna, "GGTT")
            self.assertTrue(candidate.has_stop_codon)
            self.assertFalse(candidate.excluded)
            self.assertTrue(candidate.rbs_accessibility)
            self.assertEqual(candidate.rbs_calculator["method"], "ostir_ideal_occupancy_v1")
            self.assertEqual(candidate.rbs_calculator["reference"]["ostir_version"], "1.1.3")
            for name in ("score", "defect", "on_yield_pct", "leak_pct", "delta_g_off", "delta_g_on",
                         "ddg_activation", "delta_g_rbs_linker"):
                self.assertIsNotNone(getattr(candidate, name))
            self.assertEqual(load_run(exports["json"]).result.to_dict(), result.to_dict())
            snapshot = candidate_snapshot(result, candidate, "Native reading")
            saved = save_switch(Path(folder).resolve() / "retained", snapshot["name"], snapshot,
                                source_output_dir=result.output_dir)
            with patch("aptaswitch_core.thermodynamics.setup.get_backend", side_effect=AssertionError("Offline snapshot restoration")):
                restored, retained = load_saved_design(saved)
            self.assertEqual(retained.to_dict(), candidate.to_dict())
            self.assertEqual(restored.request.fixed_switch_sequence, options.fixed_switch_sequence)
            self.assertEqual(restored.optimization["trials"], [])
