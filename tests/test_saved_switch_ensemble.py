"""Full-context ensembles reuse complete strands and the selected reference."""
from contextlib import contextmanager
from copy import deepcopy
from dataclasses import replace
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import Mock, patch

from aptaswitch_core.architecture import default_architecture_for_trigger, target_switch_structures
from aptaswitch_core.models import RunRequest, SequenceInput
from aptaswitch_core.saved_switch_analysis import build_saved_switch_context, validate_analysis_config
from aptaswitch_core.saved_switch_ensemble import context_ensemble_source
from aptaswitch_core.structure_ensemble import (
    analyze_structure_ensemble, cluster_structures, ensemble_inputs, ensemble_key,
    load_ensemble, save_ensemble,
)
from aptaswitch_studio.saved_switch_library import candidate_snapshot
from design_fixtures import design_result_fixture
from test_structure_ensemble import matrix_for


def context_fixture(*, switch_material="rna", trigger_material="dna", include_aptamer=False,
                    temperature_c=42., sodium_m=.25):
    request = RunRequest(
        sequences=SequenceInput("AACCGGTTAACC", "ACGT" * 7 + "AC"),
        architecture=default_architecture_for_trigger(30), trials=1,
    )
    request = replace(request, thermo=replace(request.thermo,
        switch_material=switch_material, trigger_material=trigger_material))
    result = design_result_fixture(request)
    snapshot = candidate_snapshot(result, result.candidates[0], "Full reporter")
    config = validate_analysis_config(snapshot, {
        "temperature_c": temperature_c, "sodium_m": sodium_m,
        "include_aptamer": include_aptamer, "switch_concentration_m": 2e-7,
        "trigger_concentration_m": 8e-7, "aptamer_concentration_m": 1e-6,
    })
    context = build_saved_switch_context(snapshot, config)
    baseline = {"conditions": {key: config[key] for key in ("temperature_c", "sodium_m", "magnesium_m")}}
    for state in ("off", "on"):
        seqs = ([context["trigger_sequence"]] if state == "on" else []) + [context["switch_sequence"]]
        baseline[state] = {"mfe": {
            "strand_order": ["Trigger", "Switch"] if state == "on" else ["Switch"],
            "sequences": seqs, "structure": "+".join("." * len(s) for s in seqs),
        }}
    return result, snapshot, {"status": "completed", "config": config, "context": context, "baseline": baseline}


class SavedSwitchEnsembleTests(unittest.TestCase):
    def test_complete_sequence_and_reference_conditions_are_detached_from_design(self):
        design, snapshot, analysis = context_fixture()
        before = deepcopy((snapshot, analysis))
        source, candidate = context_ensemble_source(snapshot, analysis, Path("/trusted/context"))
        self.assertEqual((snapshot, analysis), before)
        self.assertEqual(candidate.switch_seq, analysis["context"]["switch_sequence"])
        self.assertTrue(candidate.switch_seq.startswith(design.candidates[0].switch_seq))
        self.assertGreater(len(candidate.switch_seq), len(design.candidates[0].switch_seq) + 600)
        self.assertEqual(source.request.reporter.sequence_after_start_rna,
                         design.request.reporter.full_sequence_after_start_rna)
        self.assertEqual(source.request.thermo.temperature_c, 42.)
        self.assertEqual(source.request.thermo.sodium_m, .25)
        self.assertEqual(source.request.thermo.magnesium_m, 0.)
        self.assertEqual(source.request.thermo.resolved_switch_concentration_m, 2e-7)
        self.assertEqual(source.request.thermo.resolved_trigger_concentration_m, 8e-7)
        self.assertEqual(source.output_dir, Path("/trusted/context"))
        snapshot["candidate"]["switch_seq"] = "ACGT"
        analysis["config"]["temperature_c"] = 90
        self.assertEqual(source.request.thermo.temperature_c, 42.)
        self.assertEqual(candidate.switch_seq, before[1]["context"]["switch_sequence"])

    def test_candidate_trigger_not_original_design_input_is_used(self):
        _, snapshot, analysis = context_fixture()
        snapshot["request"]["sequences"]["trigger_dna"] = "A" * 24
        _, candidate = context_ensemble_source(snapshot, analysis, Path("/context"))
        self.assertEqual(candidate.trigger_seq, snapshot["candidate"]["trigger_seq"])
        self.assertEqual(len(candidate.trigger_seq), 30)

    def test_reporter_acquisition_provenance_does_not_change_scientific_identity(self):
        _, snapshot, analysis = context_fixture()
        source, candidate = context_ensemble_source(snapshot, analysis, Path("/context"))
        analysis["context"]["reporter_source"] = "user_supplied_after_start"
        analysis["context"]["candidate_name"] = "Original display name"
        before = deepcopy(analysis)
        same_source, same_candidate = context_ensemble_source(snapshot, analysis, Path("/context"))
        self.assertEqual(analysis, before)
        self.assertEqual((same_source, same_candidate), (source, candidate))
        reporter = analysis["context"]["reporter_sequence_after_start"]
        analysis["context"]["reporter_sequence_after_start"] = reporter[:-1] + ("A" if reporter[-1] != "A" else "C")
        with self.assertRaises(ValueError):
            context_ensemble_source(snapshot, analysis, Path("/context"))

    def test_reference_survives_cancelled_or_failed_sweeps(self):
        _, snapshot, analysis = context_fixture()
        for status in ("cancelled", "failed"):
            with self.subTest(status=status):
                analysis["status"] = status
                source, candidate = context_ensemble_source(snapshot, analysis, Path("/context"))
                self.assertEqual(len(candidate.switch_seq), analysis["context"]["switch_length"])

    def test_reversed_context_chemistry_is_rejected_before_sampling(self):
        _, snapshot, analysis = context_fixture()
        snapshot["request"]["thermo"].update(switch_material="dna", trigger_material="rna")
        with self.assertRaises(ValueError):
            context_ensemble_source(snapshot, analysis, Path("/context"))

    def test_inconsistent_or_missing_baseline_is_rejected_before_sampling(self):
        _, snapshot, analysis = context_fixture()
        mutations = [
            lambda a: a.update(baseline=None),
            lambda a: a.update(context={}),
            lambda a: a["baseline"]["conditions"].update(temperature_c=22),
            lambda a: a["config"].update(magnesium_m=.01),
            lambda a: a["context"].update(switch_sequence=snapshot["candidate"]["switch_seq"]),
            lambda a: a["baseline"]["on"]["mfe"].update(strand_order=["Switch", "Trigger"]),
            lambda a: a["baseline"]["on"]["mfe"]["sequences"].reverse(),
            lambda a: a["baseline"]["off"]["mfe"].update(structure="."),
            lambda a: a["baseline"]["off"]["mfe"].update(structure="(" + a["baseline"]["off"]["mfe"]["structure"][1:]),
        ]
        for index, mutate in enumerate(mutations):
            with self.subTest(index=index):
                damaged = deepcopy(analysis)
                mutate(damaged)
                with self.assertRaises(ValueError):
                    context_ensemble_source(snapshot, damaged, Path("/context"))

    def test_sampling_uses_full_strands_four_regions_and_reference_native_model(self):
        _, snapshot, analysis = context_fixture(include_aptamer=True)
        source, candidate = context_ensemble_source(snapshot, analysis, Path("/context"))
        off, on = target_switch_structures(source.request.architecture,
                                          len(source.request.reporter.sequence_after_start_rna))
        from test_thermodynamics_ensemble_analysis import EnsembleFixture
        backend = EnsembleFixture(source.request)
        with patch.object(backend, "ensemble", wraps=backend.ensemble) as compute:
            ensemble = analyze_structure_ensemble(source.request, candidate, backend=backend, num_sample=100)
        self.assertEqual(ensemble["inputs"]["thermodynamics"], ensemble["model"])
        self.assertEqual([call.args[1].celsius for call in compute.call_args_list], [42., 42.])
        self.assertEqual([call.args[1].sodium for call in compute.call_args_list], [.25, .25])
        self.assertEqual([[s.name for s in call.args[0].strands] for call in compute.call_args_list],
                         [["Switch"], ["Trigger", "Switch"]])
        self.assertEqual(ensemble["states"]["off"]["strand_lengths"], [len(candidate.switch_seq)])
        self.assertEqual(ensemble["states"]["on"]["strand_lengths"], [30, len(candidate.switch_seq)])
        self.assertEqual(set(ensemble["states"]["off"]["regions"]), {"off_toehold", "off_stem"})
        self.assertEqual(set(ensemble["states"]["on"]["regions"]), {"on_rbs_aug", "on_rbs_linker"})
        for state in ensemble["states"].values():
            self.assertEqual(cluster_structures(state)["sample_count"], 100)
        with tempfile.TemporaryDirectory() as folder:
            save_ensemble(Path(folder), ensemble)
            self.assertEqual(load_ensemble(Path(folder), ensemble["inputs"]), ensemble)
            # Dense full-reporter matrices legitimately exceed 32 MiB. Check
            # the file-size allowance without constructing a huge test file.
            with patch.object(Path, "stat", return_value=SimpleNamespace(st_size=33 * 1024 * 1024)):
                self.assertEqual(load_ensemble(Path(folder), ensemble["inputs"]), ensemble)
            with patch.object(Path, "stat", return_value=SimpleNamespace(st_size=513 * 1024 * 1024)):
                with self.assertRaisesRegex(ValueError, "too large"):
                    load_ensemble(Path(folder), ensemble["inputs"])

    def test_cache_separates_full_context_and_reference_temperature_or_salt(self):
        design, snapshot, analysis = context_fixture()
        source, candidate = context_ensemble_source(snapshot, analysis, Path("/context"))
        full_key = ensemble_key(ensemble_inputs(source.request, candidate))
        self.assertNotEqual(full_key, ensemble_key(ensemble_inputs(design.request, design.candidates[0])))
        for conditions in ({"temperature_c": 43.}, {"sodium_m": .35}):
            _, other_snapshot, other_analysis = context_fixture(**conditions)
            other_source, other_candidate = context_ensemble_source(other_snapshot, other_analysis, Path("/context"))
            self.assertNotEqual(full_key, ensemble_key(ensemble_inputs(other_source.request, other_candidate)))
