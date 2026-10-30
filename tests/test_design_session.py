"""Named projects retain scientific inputs and artifacts across switches/moves."""

from __future__ import annotations

import json
import shutil
import tempfile
import unittest
from dataclasses import asdict, replace
from pathlib import Path
from types import MappingProxyType
from unittest.mock import Mock, patch

from aptaswitch_core.architecture import default_architecture_for_trigger
from aptaswitch_core.models import REPORTERS, ReporterContext, RunRequest, SequenceInput, ThermoConditions
from aptaswitch_core.tube_analysis import TubeAnalysisConditions
from aptaswitch_studio.design_library import create_design, import_design, load_design
from aptaswitch_studio.design_session import (
    active_job_running, attach_saved_switch_jobs, open_design_into_session,
    project_reporters, remember_original_sequences, reset_design_session,
    save_active_design, validate_design_state,
)
from aptaswitch_studio.web_runtime import JobSnapshot
from tests.design_fixtures import design_result_fixture


class DesignSessionTests(unittest.TestCase):
    def setUp(self):
        from tests.native_ui_fixtures import patch_native_runtime
        patch_native_runtime(self)
        self.root = Path(self.enterContext(tempfile.TemporaryDirectory())).resolve()
        self.enterContext(patch("aptaswitch_studio.design_session.configured_library_path", return_value="/current/native"))
        self.enterContext(patch("aptaswitch_studio.design_session.all_reporters", return_value={"sfgfp": REPORTERS["sfgfp"]}))
        self.project = create_design(self.root, "Design A")
        self.session = {"ui_language": "fr", "settings_thermodynamics_library_path": "/user/native",
                        "dependency_detection_done": True}
        open_design_into_session(self.project.path, self.session)

    def design_result(self, run_id="design_one", *, architecture=None):
        request = RunRequest(sequences=SequenceInput("ACGTACGT", "ACGT" * 6, "System"),
                             architecture=architecture or default_architecture_for_trigger(24),
                             output_dir=self.project.path / "calculations", trials=1)
        result = design_result_fixture(request)
        result.run_id = run_id
        result.output_dir = self.project.path / "calculations" / "designs" / run_id
        result.output_dir.mkdir(parents=True)
        path = result.output_dir / f"{run_id}_results.json"
        path.write_text(json.dumps(result.to_dict()))
        return result, path

    def extension_result(self):
        source = Path(__file__).parent / "fixtures" / "extension_native_fixture.json"
        data = json.loads(source.read_text())
        output = self.project.path / "calculations" / "extensions" / data["run_id"]
        output.mkdir(parents=True)
        path = output / f"{data['run_id']}_results.json"
        data["output_dir"] = str(output)
        path.write_text(json.dumps(data))
        return data, path


    def test_form_roundtrip_preserves_sequence_edits_and_discards_retired_geometry(self):
        retired = dict(binding_start=4, binding_end=17, customize_tswitch=True,
                       toehold_length=8, rbs_aug_distance=11)
        self.session.update(aptamer_input="ACGT", trigger_input=" acgu\n" * 6,
                            **retired, t7_leader="AAA", rbs_loop="AGGAGGA",
                            frame_linker="GCUAAC", linker_auto=False,
                            extension_ligand_positions=[2, 4], extension_ligand_snapshot="ACGT",
                            settings_oxdna_path="/private/engine", unknown_secret="not portable")
        saved = save_active_design(self.session)
        self.assertNotIn("settings_oxdna_path", saved.state["form"])
        self.assertNotIn("unknown_secret", saved.state["form"])
        self.assertTrue(set(retired).isdisjoint(saved.state["form"]))
        self.assertTrue(set(retired).isdisjoint(self.session))
        reopened = {"ui_language": "en", "settings_oxdna_path": "/new/engine"}
        open_design_into_session(self.project.path, reopened)
        self.assertTrue(set(retired).isdisjoint(reopened))
        self.assertEqual(reopened["trigger_snapshot"], "ACGT" * 6)
        self.assertEqual(reopened["t7_leader"], "AAA")
        self.assertEqual(reopened["rbs_loop"], "AGGAGGA")
        self.assertEqual(reopened["frame_linker"], "GCUAAC")
        self.assertEqual(reopened["extension_ligand_positions"], [2, 4])
        self.assertEqual(reopened["settings_oxdna_path"], "/new/engine")
        self.assertEqual(reopened["ui_language"], "en")

    def test_historical_geometry_form_is_accepted_but_not_restored_or_recaptured(self):
        retired = dict(binding_start=4, binding_end=17, customize_tswitch=True,
                       toehold_length=8, rbs_aug_distance=11)
        editable = dict(t7_leader="AAA", rbs_loop="AGGAGGA", frame_linker="CCX", linker_auto=False)
        historical = create_design(self.root, "Historical geometry", state={"form": {**retired, **editable}})
        imported = import_design(historical.path, self.root, name="Imported historical geometry")
        self.assertEqual(imported.state['form'], {**retired, **editable})
        manifest = historical.path / "design.json"
        before = manifest.read_bytes()
        self.session.update(retired)
        validate_design_state(historical.state)
        open_design_into_session(historical.path, self.session)
        self.assertTrue(set(retired).isdisjoint(self.session))
        self.assertEqual({key: self.session[key] for key in editable}, editable)
        self.assertEqual(manifest.read_bytes(), before)
        saved = save_active_design(self.session)
        self.assertTrue(set(retired).isdisjoint(saved.state["form"]))
        self.assertEqual({key: saved.state["form"][key] for key in editable}, editable)
        self.session.update(retired)
        reset_design_session(self.session)
        self.assertTrue(set(retired).isdisjoint(self.session))

    def test_historical_custom_result_architecture_is_not_normalized_by_form_migration(self):
        historical_architecture = replace(default_architecture_for_trigger(24),
                                          toehold_length=10, stem_length=14,
                                          architecture_name="historical_custom_geometry")
        result, path = self.design_result(architecture=historical_architecture)
        before = path.read_bytes()
        self.session.update(last_result=result, last_exports={"json": str(path)}, customize_tswitch=True)
        save_active_design(self.session)
        reopened = {}
        open_design_into_session(self.project.path, reopened)
        self.assertEqual(reopened["last_result"].request.architecture, historical_architecture)
        self.assertNotIn("customize_tswitch", reopened)
        self.assertEqual(path.read_bytes(), before)

    def test_originals_follow_initial_edits_and_survive_extension_and_reopen(self):
        self.session.update(aptamer_input="ACGT", trigger_input="AC")
        save_active_design(self.session)
        self.session["trigger_input"] = "ACGTAC"
        save_active_design(self.session)
        remember_original_sequences(self.session)
        self.session.update(design_using_extension=True, trigger_input="TTACGTACGG")
        save_active_design(self.session)
        reopened = {}
        open_design_into_session(self.project.path, reopened)
        self.assertEqual(reopened["design_original_trigger"], "ACGTAC")
        self.assertEqual(reopened["trigger_input"], "TTACGTACGG")
        self.assertTrue(reopened["design_using_extension"])
        original = json.loads((self.project.path / "sequences" / "original.json").read_text())
        current = json.loads((self.project.path / "sequences" / "current.json").read_text())
        self.assertEqual(original["trigger"], "ACGTAC")
        self.assertEqual(current["trigger"], "TTACGTACGG")

    def test_extension_opt_in_and_method_are_scoped_to_each_saved_design(self):
        self.session.update(extension_enabled=True, inverse_extension_enabled=True,
                            inverse_extension_target_length=30, inverse_extension_top=37,
                            extension_top=42, design_section="extension")
        saved = save_active_design(self.session)
        self.assertTrue(saved.state["form"]["extension_enabled"])
        self.assertNotIn("design_section", saved.state["form"])
        other = create_design(self.root, "Without extension")
        open_design_into_session(other.path, self.session)
        self.assertNotIn("extension_enabled", self.session)
        self.assertNotIn("design_section", self.session)
        open_design_into_session(self.project.path, self.session)
        for key, value in {"extension_enabled": True, "inverse_extension_enabled": True,
                           "inverse_extension_target_length": 30, "inverse_extension_top": 37,
                           "extension_top": 42}.items():
            self.assertEqual(self.session[key], value)
        self.session.update(extension_enabled=False, inverse_extension_enabled=False)
        save_active_design(self.session)
        reopened = {}
        open_design_into_session(self.project.path, reopened)
        self.assertFalse(reopened["extension_enabled"])
        self.assertFalse(reopened["inverse_extension_enabled"])
        self.assertEqual(reopened["inverse_extension_top"], 37)
        with self.assertRaises(ValueError):
            validate_design_state({"form": {"extension_enabled": "yes"}})

    def test_regional_weights_roundtrip_keep_raw_values_and_reset_between_designs(self):
        weights = {"weight_functional_off_stem": 0.0, "weight_functional_off_toehold": 0.1,
                   "weight_functional_on_rbs_aug": 4.5, "weight_functional_on_rbs_linker": 20.0}
        self.session.update(weights, weight_defect=6.0)
        saved = save_active_design(self.session)
        self.assertEqual({key: saved.state["form"][key] for key in weights}, weights)
        reopened = {}
        open_design_into_session(self.project.path, reopened)
        self.assertEqual({key: reopened[key] for key in weights}, weights)
        self.assertEqual(reopened["weight_defect"], 6.0)
        other = create_design(self.root, "Different priorities")
        open_design_into_session(other.path, reopened)
        self.assertTrue(all(key not in reopened for key in weights))
        open_design_into_session(self.project.path, reopened)
        self.assertEqual({key: reopened[key] for key in weights}, weights)
        reset_design_session(reopened)
        self.assertTrue(all(key not in reopened for key in weights))

    def test_regional_weight_drafts_allow_zero_sum_but_reject_invalid_widget_values(self):
        keys = ("weight_functional_off_stem", "weight_functional_off_toehold",
                "weight_functional_on_rbs_aug", "weight_functional_on_rbs_linker")
        zero_draft = dict.fromkeys(keys, 0.0)
        validate_design_state({"form": zero_draft})
        self.session.update(zero_draft)
        save_active_design(self.session)
        reopened = {}
        open_design_into_session(self.project.path, reopened)
        self.assertEqual({key: reopened[key] for key in keys}, zero_draft)
        for key in keys:
            for value in (-0.1, 20.1, True, "1.0", float("nan"), float("inf")):
                with self.subTest(key=key, value=value), self.assertRaises(ValueError):
                    validate_design_state({"form": {key: value}})


    def test_switch_saves_previous_and_removes_other_project_results_widgets_and_originals(self):
        self.session.update(aptamer_input="AAAA", trigger_input="TTTT", last_result="old",
                            extension_result={"old": 1}, extension_ligand_base_4=True,
                            selected_trial=8, extension_selected_30=2, design_reporter=asdict(ReporterContext()))
        other = create_design(self.root, "Design B", state={"form": {"aptamer_input": "GGGG"}})
        open_design_into_session(other.path, self.session)
        self.assertEqual(load_design(self.project.path).state["form"]["trigger_input"], "TTTT")
        self.assertEqual(self.session["aptamer_input"], "GGGG")
        for key in ("last_result", "extension_result", "extension_ligand_base_4",
                    "selected_trial", "extension_selected_30", "design_reporter"):
            self.assertNotIn(key, self.session)
        self.assertEqual(self.session["design_original_trigger"], "")
        self.assertEqual(self.session["settings_thermodynamics_library_path"], "/user/native")

    def test_all_running_workers_prevent_switch_and_reset_without_mutation(self):
        other = create_design(self.root, "Design B")
        for key in ("design_job", "extension_job", "system_tube_job", "blast_job"):
            with self.subTest(job=key):
                job = Mock()
                job.snapshot.return_value = JobSnapshot("running", 0, 1, "")
                self.session[key] = job
                self.assertTrue(active_job_running(self.session))
                before = dict(self.session)
                with self.assertRaises(ValueError):
                    open_design_into_session(other.path, self.session)
                with self.assertRaises(ValueError):
                    reset_design_session(self.session)
                self.assertEqual(self.session, before)
                del self.session[key]

    def test_ffs_can_be_excluded_from_busy_state_without_losing_its_worker(self):
        ffs = Mock()
        ffs.snapshot.return_value = JobSnapshot("running", 1, 4, "FFS in progress")
        self.enterContext(patch("aptaswitch_studio.saved_switch_runtime.active_saved_analysis_for_design", return_value=None))
        find_ffs = self.enterContext(patch("aptaswitch_studio.saved_switch_ffs_runtime.active_ffs_for_design", return_value=ffs))
        self.enterContext(patch("aptaswitch_studio.saved_switch_free_md_runtime.active_free_md_for_design", return_value=None))
        self.enterContext(patch("aptaswitch_studio.saved_switch_ensemble_runtime.active_context_ensemble_for_design", return_value=None))
        self.enterContext(patch("aptaswitch_studio.web_runtime.find_design_job", return_value=None))

        self.assertFalse(attach_saved_switch_jobs(self.session, ignore_ffs=True))
        self.assertIs(self.session["saved_switch_ffs_job"], ffs)
        find_ffs.assert_called_once_with(self.session["active_design_path"])
        self.assertTrue(attach_saved_switch_jobs(self.session))
        self.assertFalse(active_job_running(self.session, ignore_ffs=True))
        self.assertTrue(active_job_running(self.session))
        self.assertIs(self.session["saved_switch_ffs_job"], ffs)
        ffs.cancel.assert_not_called()

        readonly = MappingProxyType({"active_design_path": self.session["active_design_path"]})
        self.assertFalse(attach_saved_switch_jobs(readonly, ignore_ffs=True))
        self.assertTrue(attach_saved_switch_jobs(readonly))
        self.assertFalse(active_job_running(readonly, ignore_ffs=True))
        self.assertTrue(active_job_running(readonly))
        self.assertNotIn("saved_switch_ffs_job", readonly)

    def test_ignoring_ffs_still_blocks_for_free_md_and_every_other_running_job(self):
        ffs, other = Mock(), Mock()
        ffs.snapshot.return_value = JobSnapshot("running", 1, 4, "FFS in progress")
        other.snapshot.return_value = JobSnapshot("running", 0, 1, "Other calculation")
        self.enterContext(patch("aptaswitch_studio.saved_switch_runtime.active_saved_analysis_for_design", return_value=None))
        self.enterContext(patch("aptaswitch_studio.saved_switch_ffs_runtime.active_ffs_for_design", return_value=ffs))
        find_free_md = self.enterContext(patch("aptaswitch_studio.saved_switch_free_md_runtime.active_free_md_for_design", return_value=other))
        self.enterContext(patch("aptaswitch_studio.saved_switch_ensemble_runtime.active_context_ensemble_for_design", return_value=None))
        self.enterContext(patch("aptaswitch_studio.web_runtime.find_design_job", return_value=None))

        self.assertTrue(attach_saved_switch_jobs(self.session, ignore_ffs=True))
        self.assertIs(self.session["saved_switch_free_md_job"], other)
        self.assertIs(self.session["saved_switch_ffs_job"], ffs)
        self.assertTrue(active_job_running(self.session, ignore_ffs=True))
        find_free_md.return_value = None
        self.session.pop("saved_switch_free_md_job")
        for key in ("design_job", "extension_job", "system_tube_job", "blast_job",
                    "ensemble_job", "saved_switch_analysis_job", "saved_switch_free_md_job"):
            with self.subTest(job=key):
                session = {**self.session, key: other}
                other.snapshot.return_value = JobSnapshot("running", 0, 1, "Other calculation")
                self.assertTrue(active_job_running(session, ignore_ffs=True))
                self.assertIs(session["saved_switch_ffs_job"], ffs)
                other.snapshot.return_value = JobSnapshot("done", 1, 1, "Finished")
                self.assertFalse(active_job_running(session, ignore_ffs=True))

    def test_blast_preferences_reports_and_relative_exports_survive_import_without_consent_or_network(self):
        from aptaswitch_studio.blast_runtime import run_blast_and_export
        from aptaswitch_studio.design_library import export_design
        from tests.test_blast_runtime import client_result, screen_request
        request = screen_request()
        result, exports = run_blast_and_export(request, self.project.path / "calculations",
                                               client=Mock(run=Mock(return_value=client_result())))
        job = Mock()
        job.snapshot.return_value = JobSnapshot("done", 2, 2, "", result=result, exports=exports)
        self.session.update(blast_taxa=[{"taxid": 9606, "name": "Homo sapiens"}],
                            blast_scopes=["genomic", "transcript"], blast_max_hits=100,
                            blast_job=job, blast_consent=True, blast_consent_signature=request.fingerprint(),
                            blast_lookup_query="human", blast_lookup_results=[{"name": "temporary"}])
        saved = save_active_design(self.session)
        self.assertEqual(saved.state["blast"], {"taxa": [{"taxid": 9606, "name": "Homo sapiens"}],
                                               "scopes": ["genomic", "transcript"], "max_hits": 100})
        self.assertNotIn("consent", json.dumps(saved.state))
        self.assertNotIn("lookup", json.dumps(saved.state))
        self.assertEqual(saved.state["selected_runs"]["blast"], Path(exports["json"]).relative_to(self.project.path).as_posix())
        imported = import_design(export_design(self.project.path), self.root / "imported")
        reopened = {"blast_consent": True, "blast_consent_signature": "previous-project",
                    "blast_lookup_query": "previous", "blast_result": {"old": True}}
        with patch("aptaswitch_studio.blast_runtime.BlastClient", side_effect=AssertionError("Never resume imported RID")):
            open_design_into_session(imported.path, reopened)
        self.assertEqual(reopened["blast_taxa"], saved.state["blast"]["taxa"])
        self.assertEqual(reopened["blast_result"]["fingerprint"], result["fingerprint"])
        self.assertEqual(reopened["blast_result"]["runs"], result["runs"])
        for key in ("blast_job", "blast_consent", "blast_consent_signature", "blast_lookup_query"):
            self.assertNotIn(key, reopened)
        self.assertTrue(all(Path(path).is_relative_to(imported.path) and Path(path).is_file()
                            for path in reopened["blast_exports"].values()))

    def test_blast_failed_result_is_captured_before_fragment_consumes_it(self):
        from aptaswitch_studio.blast_runtime import run_blast_and_export
        from tests.test_blast_runtime import screen_request
        result, exports = run_blast_and_export(screen_request(), self.project.path / "calculations",
                                               client=Mock(run=Mock(side_effect=TimeoutError("Offline"))))
        worker = Mock()
        worker.snapshot.return_value = JobSnapshot("failed", 0, 2, "", result=result, exports=exports, error="Offline")
        self.session["blast_job"] = worker
        save_active_design(self.session)
        self.assertEqual(self.session["blast_result"]["status"], "failed")
        reopened = {}
        open_design_into_session(self.project.path, reopened)
        self.assertEqual(reopened["blast_result"]["status"], "failed")
        self.assertEqual(reopened["blast_result"]["error"], "Offline")

    def test_blast_restores_latest_completed_report_if_fragment_did_not_save_its_index(self):
        from aptaswitch_studio.blast_runtime import run_blast_and_export
        from tests.test_blast_runtime import client_result, screen_request
        result, exports = run_blast_and_export(screen_request(), self.project.path / "calculations",
                                               client=Mock(run=Mock(return_value=client_result())))
        later, later_exports = run_blast_and_export(screen_request(), self.project.path / "calculations",
                                                   client=Mock(run=Mock(return_value=client_result(status="partial"))))
        first_path, second_path = Path(exports["json"]), Path(later_exports["json"])
        # Filename ordering is intentionally the opposite of completion order.
        first_path.parent.rename(first_path.parent.parent / "zzz_old")
        second_path.parent.rename(second_path.parent.parent / "aaa_new")
        reopened = {}
        open_design_into_session(self.project.path, reopened)
        self.assertEqual(reopened["blast_result"]["completed_at"], later["completed_at"])
        self.assertEqual(reopened["blast_result"]["status"], "partial")

    def test_blast_draft_preferences_allow_empty_selection_and_reject_bad_import_parameters(self):
        validate_design_state({"blast": {"taxa": [], "scopes": [], "max_hits": 500}})
        for blast in ({"taxa": "human"}, {"taxa": [{"taxid": True, "name": "Human"}]},
                      {"taxa": [{"taxid": 9606, "name": ""}]},
                      {"taxa": [{"taxid": 9606, "name": "Homo sapiens"}] * 2},
                      {"scopes": ["everything"]}, {"scopes": [{}]},
                      {"scopes": ["genomic", "genomic"]}, {"max_hits": 501}):
            with self.subTest(blast=blast), self.assertRaises(ValueError):
                validate_design_state({"blast": blast})

    def test_invalid_target_preserves_active_state(self):
        target = create_design(self.root, "Bad state", state={"form": {"trials": "oops"}})
        before = dict(self.session)
        with self.assertRaises(ValueError):
            open_design_into_session(target.path, self.session)
        self.assertEqual(self.session, before)

    def test_unchanged_save_leaves_manifest_timestamp_unchanged(self):
        self.session.update(aptamer_input="AAAA", trigger_input="TTTT")
        saved = save_active_design(self.session)
        stamp = (self.project.path / "design.json").stat().st_mtime_ns
        self.assertEqual(save_active_design(self.session).updated_at, saved.updated_at)
        self.assertEqual((self.project.path / "design.json").stat().st_mtime_ns, stamp)

    def test_finished_workers_are_captured_before_their_fragment_consumes(self):
        result, path = self.design_result()
        job = Mock()
        job.snapshot.return_value = JobSnapshot("done", 1, 1, "", result=result, exports={"json": str(path)})
        self.session["design_job"] = job
        saved = save_active_design(self.session)
        self.assertIs(self.session["last_result"], result)
        self.assertEqual(saved.state["selected_runs"]["design"], path.relative_to(self.project.path).as_posix())
        job.consume_once.assert_not_called()

    def test_paused_design_is_selected_and_restored_from_disk(self):
        result, old_path = self.design_result()
        old_path.unlink()
        result.status = "paused"
        result.optimization["checkpoint"] = {"phase": "design", "completed_trials": 1}
        path = result.output_dir / "design_checkpoint.json"
        path.write_text(json.dumps(result.to_dict()))
        job = Mock()
        job.snapshot.return_value = JobSnapshot("paused", 1, 1, "", result=result,
                                                exports={"checkpoint": str(path)})
        self.session["design_job"] = job
        saved = save_active_design(self.session)
        self.assertEqual(saved.state["selected_runs"]["design"], path.relative_to(self.project.path).as_posix())
        reopened = {}
        with patch("aptaswitch_core.thermodynamics.setup.get_backend", side_effect=AssertionError("Do not recompute")):
            open_design_into_session(self.project.path, reopened)
        self.assertEqual(reopened["last_result"].status, "paused")
        self.assertEqual(reopened["last_result"].run_id, result.run_id)
        self.assertEqual(reopened["opened_design_source"], str(path))

    def test_active_design_worker_is_reattached_but_finished_job_does_not_replace_selected_run(self):
        job = Mock()
        job.snapshot.return_value = JobSnapshot("running", 1, 3, "")
        with patch("aptaswitch_studio.web_runtime.find_design_job", return_value=job) as find:
            self.assertTrue(active_job_running(self.session))
        find.assert_called_once_with(self.project.path / "calculations")
        self.assertIs(self.session["design_job"], job)
        self.session.pop("design_job")
        job.snapshot.return_value = JobSnapshot("paused", 1, 3, "")
        with patch("aptaswitch_studio.web_runtime.find_design_job", return_value=job):
            self.assertFalse(active_job_running(self.session))
        self.assertNotIn("design_job", self.session)

    def test_result_reopen_after_folder_move_relocates_design_and_extension(self):
        result, design_path = self.design_result()
        extension, extension_path = self.extension_result()
        self.session.update(last_result=result, last_exports={"json": str(design_path)},
                            extension_result=extension, extension_exports={"json": str(extension_path)})
        save_active_design(self.session)
        destination = self.root / "Moved project"
        shutil.move(str(self.project.path), destination)
        reopened = {}
        with patch("aptaswitch_core.thermodynamics.setup.get_backend", side_effect=AssertionError("Do not recompute")):
            open_design_into_session(destination, reopened)
        self.assertEqual(reopened["last_result"].candidates[0].score, result.candidates[0].score)
        self.assertEqual(reopened["extension_result"]["candidates"], extension["candidates"])

    def test_saved_custom_reporter_overrides_same_global_key_only_for_this_project(self):
        reporter = ReporterContext("sfgfp", "Project custom", "ACG" * 10, "ACG" * 15 + "UAA")
        self.session.update(reporter_key="sfgfp", design_reporter=asdict(reporter))
        save_active_design(self.session)
        reopened = {}
        open_design_into_session(self.project.path, reopened)
        self.assertEqual(project_reporters(reopened)["sfgfp"], reporter)
        self.assertEqual(project_reporters({})["sfgfp"], REPORTERS["sfgfp"])


    def test_initial_structure_and_tube_signatures_restore_without_local_dependency_paths(self):
        thermo, conditions = ThermoConditions(material="experimental-Sugimoto95-Vienna04-loopmean", switch_material="dna", trigger_material="dna"), TubeAnalysisConditions()
        structure = {"mfe_structure": "....+....", "paired_trigger_indices": []}
        tube_result = {"status": "completed", "tubes": {}}
        tube = Mock(signature=("AAAA", "TTTT", "/old/native", thermo, conditions))
        tube.snapshot.return_value = JobSnapshot("done", 2, 2, "", result=tube_result, log=("done",))
        self.session.update(system_analysis_signature=("AAAA", "TTTT", "/old/native", "experimental-Sugimoto95-Vienna04-loopmean", "dna", 37., .12, 0.),
                            system_analysis_result=structure, system_tube_job=tube)
        saved = save_active_design(self.session)
        self.assertNotIn("/old/native", json.dumps(saved.state))
        reopened = {}
        open_design_into_session(self.project.path, reopened)
        self.assertEqual(reopened["system_analysis_signature"], ("AAAA", "TTTT", "/current/native", "experimental-Sugimoto95-Vienna04-loopmean", "dna", 37., .12, 0., None))
        self.assertEqual(reopened["system_analysis_result"], structure)
        self.assertTrue(reopened["system_analysis_requires_refresh"])
        self.assertEqual(reopened["system_tube_job"].signature, ("AAAA", "TTTT", "/current/native", thermo, conditions, None))
        self.assertEqual(reopened["system_tube_job"].snapshot().result, tube_result)
        self.assertTrue(reopened["system_tube_job"].restored)
        self.assertFalse(reopened["system_tube_job"].consume_once())

    def test_historical_five_strand_tubes_restore_only_as_a_read_only_snapshot(self):
        from aptaswitch_studio.design_session import _restore_initial_analysis
        thermo = asdict(ThermoConditions(material="experimental-Sugimoto95-Vienna04-loopmean", switch_material="dna", trigger_material="dna"))
        thermo.update(screening_max_complex_size=5, validation_max_complex_size=5)
        conditions = {**asdict(TubeAnalysisConditions()), "max_complex_size": 5}
        result = {"status": "completed", "conditions": dict(conditions), "tubes": {}}
        saved = {"tubes": {"aptamer": "AAAA", "trigger": "TTTT", "thermo": thermo,
                           "conditions": conditions, "result": result}}
        reopened = {}
        _restore_initial_analysis(saved, reopened)
        job = reopened["system_tube_job"]
        self.assertTrue(job.restored)
        self.assertFalse(job.consume_once())
        self.assertEqual(job.signature[4].max_complex_size, 5)
        self.assertNotIsInstance(job.signature[4], TubeAnalysisConditions)
        self.assertNotEqual(job.signature[4], TubeAnalysisConditions())
        self.assertEqual(job.snapshot().result["conditions"]["max_complex_size"], 5)
        self.assertEqual(saved["tubes"]["conditions"]["max_complex_size"], 5)

    def test_native_initial_analysis_keeps_recorded_identity_across_runtime_changes(self):
        from aptaswitch_core.thermodynamics.setup import MODEL_ID
        identity = json.dumps({"backend": "native-hybrid", "parameter_digest": "sha256:recorded"})
        structure = {"mfe_structure": "....+....", "paired_trigger_indices": []}
        self.session.update(system_analysis_signature=("AAAA", "TTTT", "/old/native", MODEL_ID,
                            "dna", 37., .12, 0., identity), system_analysis_result=structure)
        saved = save_active_design(self.session)
        self.assertNotIn("/old/native", json.dumps(saved.state))
        reopened = {}
        open_design_into_session(self.project.path, reopened)
        self.assertEqual(reopened["system_analysis_signature"][-1], identity)
        self.assertEqual(reopened["system_analysis_result"], structure)
        self.assertTrue(reopened["system_analysis_requires_refresh"])


    def test_rejects_malformed_portable_state_before_import(self):
        for state in ({"form": []}, {"original_sequences": []},
                      {"form": {"temperature_c": float("nan")}},
                      {"form": {"linker_auto": "yes"}},
                      {"selected_runs": {"design": "/native/run.json"}},
                      {"selected_runs": {"design": "../run.json"}},
                      {"reporter": {"key": "x"}},
                      {"initial_analysis": {"structure": {"signature": None}}}):
            with self.subTest(state=state), self.assertRaises(ValueError):
                validate_design_state(state)

    def test_invalid_widget_values_cannot_replace_existing_design_or_active_session(self):
        self.session.update(aptamer_input="ACGT", trigger_input="ACGTAC", molecule_name="Keep me")
        save_active_design(self.session)
        previous_manifest = (self.project.path / "design.json").read_bytes()
        previous_session = dict(self.session)
        invalid_states = [
            {"form": {"trials": 9999999}}, {"form": {"temperature_c": -1}},
            {"form": {"optimization_population_size": 513}},
            {"form": {"optimization_lineages": 65}},
            {"form": {"optimization_local_tolerance": -0.001}},
            {"form": {"optimization_lineages": 8, "optimization_population_size": 4}},
            {"form": {"temperature_c": 101}}, {"form": {"extension_temperature": -1}},
            {"form": {"extension_top": 1001}}, {"form": {"extension_parasite_top": 1001}},
            {"form": {"weight_on": 21}},
            {"form": {"weight_defect": -0.1}}, {"form": {"extension_target_length": 1}},
            {"form": {"extension_target_length": 31, "extension_allow_long": False}},
            {"form": {"system_tube_volume_ul": 0.0009}},
            {"form": {"trigger_concentration_value": 10001, "trigger_concentration_unit": "µM"}},
            {"form": {"switch_concentration_value": 0.5, "switch_concentration_unit": "nM"}},
            {"form": {"sodium_concentration_value": 1001, "sodium_concentration_unit": "mM"}},
            {"form": {"magnesium_concentration_value": .00001, "magnesium_concentration_unit": "M"}},
            {"form": {"system_aptamer_concentration": 0.0001, "system_aptamer_concentration_unit": "nM"}},
            {"form": {"extension_sodium": 119}}, {"form": {"extension_magnesium": .01}},
            {"form": {"trigger_concentration_unit": "M"}},
            {"form": {"extension_ensemble_cap": 2 ** 53}},
            {"reporter": {"key": "bad", "label": [], "sequence_after_start_rna": "ACGU"}},
            {"form": {"design_reporter": {"key": "bad", "label": "Bad", "sequence_after_start_rna": "ACGX"}}},
        ]
        for index, state in enumerate(invalid_states):
            with self.subTest(state=state):
                incoming = create_design(self.root / "incoming", f"Invalid {index}", state=state)
                with self.assertRaises(ValueError):
                    open_design_into_session(incoming.path, self.session)
                self.assertEqual(self.session, previous_session)
                with self.assertRaises(ValueError):
                    import_design(incoming.path, self.root, name=self.project.name, overwrite=True)
                self.assertEqual((self.project.path / "design.json").read_bytes(), previous_manifest)
                self.assertFalse((self.root / ".backups").exists())

    def test_widget_bounds_accept_endpoints_units_and_unfinished_sequence_edits(self):
        form = {
            "aptamer_input": "", "trigger_input": " acg ", "frame_linker": "",
            "trials": 100000, "temperature_c": 0., "extension_temperature": 100.,
            "extension_top": 1000, "extension_parasite_top": 0,
            "weight_on": 20., "weight_leak": 0., "t7_leader": "AAA", "rbs_loop": "AGGAGGA",
            "system_tube_volume_ul": 0.001, "extension_target_length": 31,
            "extension_allow_long": True,
            "trigger_concentration_value": 1., "trigger_concentration_unit": "nM",
            "switch_concentration_value": 10., "switch_concentration_unit": "mM",
            "sodium_concentration_value": 1., "sodium_concentration_unit": "M",
            "magnesium_concentration_value": 0., "magnesium_concentration_unit": "mM",
            "system_aptamer_concentration": .001, "system_aptamer_concentration_unit": "nM",
            "system_trigger_concentration": 1e-6, "system_trigger_concentration_unit": "µM",
            "extension_sodium": 120., "extension_magnesium": 0.,
        }
        validate_design_state({"form": form})
        form.update(aptamer_input="ACG?", trigger_input="", frame_linker="in progress")
        validate_design_state({"form": form})



if __name__ == "__main__":
    unittest.main()
