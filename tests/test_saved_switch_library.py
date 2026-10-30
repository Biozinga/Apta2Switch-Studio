"""A retained candidate must remain exact, recoverable and self-contained."""

from __future__ import annotations

import json
from copy import deepcopy
import tempfile
import unittest
from concurrent.futures import ThreadPoolExecutor
from dataclasses import replace
from pathlib import Path
from unittest.mock import patch

from aptaswitch_core.models import REPORTERS, ReporterContext, InverseExtensionSettings, RunRequest, SequenceInput, SwitchArchitecture
from aptaswitch_studio import saved_switch_library as library
from tests.design_fixtures import design_result_fixture


class SavedSwitchLibraryTests(unittest.TestCase):
    def test_design_archive_carries_saved_candidate_full_reporter_and_context(self):
        from aptaswitch_studio.design_library import create_design, export_design, import_design
        from aptaswitch_studio.saved_switch_runtime import save_context, load_context

        project = create_design(self.folder / "projects", "Experiment")
        result = design_result_fixture(replace(self.request, reporter=REPORTERS["sfgfp"]))
        snapshot = library.candidate_snapshot(result, result.candidates[0], "Retained")
        saved = library.save_switch(library.saved_switch_root(project.path), snapshot["name"], snapshot)
        save_context(saved, {})
        marker = (saved.path / "calculations" / "structures3d") / "retained_model.pdb"
        marker.parent.mkdir()
        marker.write_text("REMARK test artifact\n")
        copied = import_design(export_design(project.path), self.folder / "imported")
        restored = library.list_saved_switches(library.saved_switch_root(copied.path))[0]
        self.assertEqual(restored.snapshot, saved.snapshot)
        self.assertEqual(load_context(restored), load_context(saved))
        self.assertEqual(((restored.path / "calculations" / "structures3d") / marker.name).read_bytes(), marker.read_bytes())
        self.assertFalse((restored.path.parent / ".saved-switch-library.lock").exists())

    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.folder = Path(temporary.name).resolve()
        self.root = self.folder / "design" / library.LIBRARY_FOLDER
        self.request = RunRequest(
            sequences=SequenceInput("AACCGGTTAACC", "AACC"),
            architecture=SwitchArchitecture(toehold_length=2, stem_length=2),
            reporter=ReporterContext(),
            trials=2, output_dir=self.folder / "original-runs",
        )
        self.result = design_result_fixture(self.request)
        self.candidate = self.result.candidates[0]

    def snapshot(self, name="Candidat β", **kwargs):
        return library.candidate_snapshot(self.result, self.candidate, name, **kwargs)

    def save(self, name="Candidat β", **kwargs):
        return library.save_switch(self.root, name, self.snapshot(name), **kwargs)

    def test_selection_snapshot_is_detached_complete_and_has_exact_candidate_sequences(self):
        self.candidate.inverse_extension_metrics = {"comparison": {"impact": 0.15}}
        snapshot = self.snapshot(original_sequences={"aptamer": "AACCGG", "trigger": "GGTT"})
        self.assertEqual(snapshot["request"], self.result.request.to_dict())
        self.assertEqual(snapshot["candidate"], self.candidate.to_dict())
        self.assertEqual(snapshot["sequences"]["switch"], self.candidate.switch_seq)
        self.assertEqual(snapshot["sequences"]["trigger"], self.candidate.trigger_seq)
        self.assertEqual(snapshot["sequences"]["original_trigger"], "GGTT")
        self.assertEqual(snapshot["sequences"]["original_trigger_source"], "design_original_sequences")
        self.assertFalse(snapshot["sequences"]["reporter_is_complete"])
        self.candidate.inverse_extension_metrics["comparison"]["impact"] = 0.8
        self.assertEqual(snapshot["candidate"]["inverse_extension_metrics"]["comparison"]["impact"], 0.15)

    def test_candidate_must_belong_to_selected_result(self):
        stranger = replace(self.candidate, trial=99)
        with self.assertRaises(ValueError):
            library.candidate_snapshot(self.result, stranger, "Absent")

    def test_candidate_snapshot_keeps_optimizer_settings_without_run_resume_checkpoint(self):
        from aptaswitch_studio.saved_switch_design import load_saved_design
        checkpoint = {"version": 1, "phase": "design", "completed_trials": 2,
                      "final_status": "completed", "payload_digest": "original-full-run-only"}
        self.result.optimization["checkpoint"] = checkpoint
        original = deepcopy(self.result.optimization)
        snapshot = self.snapshot()
        optimization = snapshot["source"]["optimization"]
        self.assertNotIn("checkpoint", optimization)
        self.assertEqual(optimization["options"], original["options"])
        self.assertEqual(optimization["functional_weights"], original["functional_weights"])
        self.assertEqual(optimization["trials"], [original["trials"][0]])
        saved = library.save_switch(self.root, snapshot["name"], snapshot)
        restored, candidate = load_saved_design(saved)
        self.assertNotIn("checkpoint", restored.optimization)
        self.assertEqual(candidate.to_dict(), self.candidate.to_dict())
        self.assertEqual(self.result.optimization, original)
        self.assertIs(self.result.optimization["checkpoint"], checkpoint)

    def test_retained_native_trial_920_of_100000_is_small_and_reopens_offline(self):
        from aptaswitch_studio.saved_switch_design import load_saved_design
        self.result.request = replace(self.result.request, trials=100000)
        self.candidate.trial = 920
        outcome = self.result.optimization["trials"][0]
        options = outcome["options"]
        reports = [{**outcome, "trial": trial, "options": {**options, "seed": options["seed"] + trial - 1}}
                   for trial in range(1, 100001)]
        self.result.optimization["trials"] = reports
        selected_before = deepcopy(reports[919])
        snapshot = self.snapshot()
        compact = snapshot["source"]["optimization"]
        self.assertEqual(compact["scope"], "retained_candidate")
        self.assertEqual(compact["source_run_id"], self.result.run_id)
        self.assertEqual(compact["source_trial_count"], 100000)
        self.assertEqual(compact["retained_trial"], 920)
        self.assertEqual(compact["trials"], [selected_before])
        original_size = len(json.dumps(snapshot))
        self.assertLess(original_size, 25000)
        # Large unrelated diagnostics cannot make one selected candidate grow.
        reports[0]["unrelated_trace"] = "x" * (library.MAX_SNAPSHOT_BYTES + 1)
        self.assertEqual(len(json.dumps(self.snapshot())), original_size)
        with patch("aptaswitch_core.thermodynamics.setup.get_backend", side_effect=AssertionError("Offline selection")):
            saved = library.save_switch(self.root, snapshot["name"], snapshot)
            restored, selected = load_saved_design(library.load_switch(saved.path))
        self.assertEqual(selected.trial, 920)
        self.assertEqual(restored.request.trials, 100000)
        self.assertEqual(restored.thermodynamics, self.result.thermodynamics)
        self.assertEqual(restored.optimization, compact)
        self.assertIs(self.result.optimization["trials"], reports)
        self.assertEqual(len(reports), 100000)
        self.assertEqual(reports[919], selected_before)

    def test_retained_optimizer_scope_cannot_disguise_an_incomplete_full_run(self):
        from aptaswitch_studio.run_library import _load_design
        snapshot = self.snapshot()
        payload = {**snapshot["source"], "request": snapshot["request"], "results": [snapshot["candidate"]]}
        with self.assertRaisesRegex(ValueError, "not a complete design run"):
            _load_design(payload, Path("/archive/run.json"))
        mutations = [
            lambda p: p["source"]["optimization"].update(source_run_id="another-run"),
            lambda p: p["source"]["optimization"].update(source_trial_count=100),
            lambda p: p["source"]["optimization"].update(retained_trial=2),
            lambda p: p["source"]["optimization"]["trials"][0].update(trial=2),
            lambda p: p["source"]["optimization"].update(scope="full_run"),
            lambda p: p["source"]["optimization"]["trials"][0]["options"].update(seed=100),
        ]
        for index, mutate in enumerate(mutations):
            changed = deepcopy(snapshot)
            mutate(changed)
            with self.subTest(mutation=index), self.assertRaises(ValueError):
                library._validate_snapshot(changed)

    def test_actual_inverse_extension_trigger_is_preserved_over_request_core(self):
        self.result.request = replace(self.result.request, inverse_extension=InverseExtensionSettings(
            target_length=8, extension_side="5prime", top_candidates=1,
        ))
        self.candidate.trigger_seq = "GCGCAACC"
        snapshot = self.snapshot()
        saved = library.save_switch(self.root, snapshot["name"], snapshot)
        reopened = library.load_switch(saved.path)
        self.assertEqual(reopened.snapshot["candidate"]["trigger_seq"], "GCGCAACC")
        self.assertEqual(reopened.snapshot["sequences"]["trigger"], "GCGCAACC")
        self.assertEqual(reopened.snapshot["sequences"]["input_trigger"], "AACC")
        self.assertEqual(reopened.snapshot["sequences"]["original_trigger"], "AACC")

    def test_saved_candidate_with_recorded_inverse_model_stays_visible_and_exact(self):
        self.result.request = replace(self.result.request, inverse_extension=InverseExtensionSettings(
            target_length=8, extension_side="5prime", top_candidates=1,
        ))
        self.candidate.trigger_seq = "GCGCAACC"
        saved = self.save()
        snapshot = json.loads(saved.manifest_path.read_text())
        snapshot["request"]["inverse_extension"]["thermo"].update(material="dna04", magnesium_m=.01)
        saved.manifest_path.write_text(json.dumps(snapshot))
        before = saved.manifest_path.read_bytes()
        with patch("aptaswitch_core.thermodynamics.setup.get_backend", side_effect=AssertionError("Must not compute")):
            visible = library.list_saved_switches(self.root)
        self.assertEqual(len(visible), 1)
        self.assertEqual(visible[0].name, saved.name)
        self.assertEqual(visible[0].snapshot, snapshot)
        self.assertEqual(visible[0].snapshot["request"]["thermo"]["material"], self.result.request.thermo.material)
        self.assertEqual(saved.manifest_path.read_bytes(), before)

    def test_complete_reporter_and_transcript_survive_save_without_duplicate_prefix(self):
        result = design_result_fixture(replace(self.request, reporter=REPORTERS["sfgfp"]))
        candidate = result.candidates[0]
        snapshot = library.candidate_snapshot(result, candidate, "Complet")
        saved = library.save_switch(self.root, "Complet", snapshot)
        reopened = library.load_switch(saved.path)
        sequences = reopened.snapshot["sequences"]
        self.assertTrue(sequences["reporter_is_complete"])
        self.assertEqual(sequences["reporter_full_after_start"], REPORTERS["sfgfp"].full_sequence_after_start_rna)
        self.assertEqual(sequences["switch_full_transcript"], candidate.switch_seq[:-30] + REPORTERS["sfgfp"].full_sequence_after_start_rna)
        self.assertIn(">switch_full_transcript\n", reopened.fasta_path.read_text())
        self.assertEqual(reopened.snapshot["request"]["reporter"]["source_accession"], "HQ873313.1")


    def test_roundtrip_moves_without_source_run_or_current_reporter_library(self):
        snapshot = self.snapshot()
        saved = library.save_switch(self.root, snapshot["name"], snapshot)
        snapshot["candidate"]["score"] = 99
        self.assertEqual(saved.snapshot["candidate"]["score"], self.candidate.score)
        self.assertEqual(library.load_switch(saved.manifest_path), saved)
        moved = self.folder / "imported-design"
        self.root.parent.rename(moved)
        with patch("aptaswitch_core.reporter_library.all_reporters", side_effect=AssertionError("must not resolve current library")):
            reopened = library.load_switch(moved / library.LIBRARY_FOLDER / saved.name)
        self.assertEqual(reopened.snapshot, saved.snapshot)
        self.assertEqual(reopened.path.parent.parent, moved)
        self.assertEqual(library.list_saved_switches(reopened.path.parent), [reopened])
        self.assertTrue(reopened.analyses_dir.is_dir())
        self.assertFalse((reopened.path / "calculations" / "structures3d").exists())
        self.assertEqual(reopened.context_path, reopened.path / "context.json")

    def test_fasta_does_not_append_reporter_again_or_claim_full_context(self):
        saved = self.save()
        records = {}
        for part in saved.fasta_path.read_text().split(">"):
            if part:
                header, *lines = part.splitlines()
                records[header] = "".join(lines)
        self.assertEqual(records["switch_designed_including_reporter_segment"], self.candidate.switch_seq)
        self.assertEqual(records["trigger_candidate"], self.candidate.trigger_seq)
        self.assertEqual(records["reporter_segment_after_start_completeness_unspecified"], self.request.reporter.sequence_after_start_rna)

    def test_duplicate_unicode_casefold_requires_explicit_overwrite_and_backs_up_everything(self):
        first = self.save("Élu")
        artifact = (first.path / "calculations" / "structures3d") / "valuable.pdb"
        artifact.parent.mkdir()
        artifact.write_text("EXPERIMENT")
        for name in ("élu", "E\u0301LU", " Élu "):
            with self.subTest(name=name), self.assertRaises(library.SavedSwitchExistsError):
                self.save(name)
        self.assertEqual(artifact.read_text(), "EXPERIMENT")
        self.candidate.score = 75.0
        replacement = self.save("ÉLU", overwrite=True)
        self.assertEqual(replacement.snapshot["candidate"]["score"], 75.0)
        backups = list((self.root / ".backups").iterdir())
        self.assertEqual(len(backups), 1)
        old = library.load_switch(backups[0])
        self.assertEqual(old.snapshot, first.snapshot)
        self.assertEqual(((old.path / "calculations" / "structures3d") / "valuable.pdb").read_text(), "EXPERIMENT")
        self.assertEqual(library.list_saved_switches(self.root), [replacement])

    def test_publication_failure_restores_existing_selection_and_artifacts(self):
        first = self.save()
        marker = first.analyses_dir / "result.json"
        marker.write_text("old data")
        rename = Path.rename

        def fail_staging(path, destination):
            if path.name.startswith(".switch-staging-"):
                raise OSError("simulated full disk")
            return rename(path, destination)

        with patch.object(Path, "rename", fail_staging), self.assertRaises(OSError):
            self.save(overwrite=True)
        self.assertEqual(library.load_switch(first.path), first)
        self.assertEqual(marker.read_text(), "old data")
        self.assertFalse(list(self.root.glob(".switch-staging-*")))

    def test_concurrent_same_name_publishes_only_one_selection(self):
        def create():
            try:
                return self.save()
            except library.SavedSwitchExistsError:
                return None

        with ThreadPoolExecutor(max_workers=4) as pool:
            results = list(pool.map(lambda _: create(), range(4)))
        self.assertEqual(sum(item is not None for item in results), 1)
        self.assertEqual(len(library.list_saved_switches(self.root)), 1)

    def test_unsafe_names_and_nonfinite_or_inconsistent_metadata_do_not_touch_library(self):
        for name in ("", "../outside", "/tmp/outside", "a/b", "a\\b", "CON", "..", "a\x00b", "x" * 181):
            with self.subTest(name=name), self.assertRaises(ValueError):
                self.snapshot(name)
        for mutate in (
            lambda d: d["candidate"].update(score=float("nan")),
            lambda d: d["candidate"].update(inverse_extension_metrics={"unsafe": float("inf")}),
            lambda d: d["sequences"].update(trigger="AAAA"),
            lambda d: d["sequences"].update(reporter_is_complete=True),
            lambda d: d["source"].update(original_sequences={"trigger_dna": "ACGT"}),
            lambda d: d["request"]["reporter"].update(sequence_after_start_rna="GGGGGG"),
        ):
            data = self.snapshot()
            mutate(data)
            with self.subTest(data=data), self.assertRaises(ValueError):
                library.save_switch(self.root, data["name"], data)
        self.assertFalse(self.root.exists())

    def test_rejects_symlink_root_ancestor_candidate_manifest_and_backup(self):
        saved = self.save()
        link = self.folder / "linked"
        link.symlink_to(self.root, target_is_directory=True)
        with self.assertRaises(ValueError):
            library.load_switch(link / saved.name)
        with self.assertRaises(ValueError):
            library.save_switch(link, "Elsewhere", self.snapshot("Elsewhere"))
        alias = self.root / "alias"
        alias.symlink_to(saved.path, target_is_directory=True)
        self.assertEqual(library.list_saved_switches(self.root), [saved])
        with self.assertRaises(ValueError):
            library.load_switch(alias)
        real = saved.path / "real.json"
        saved.manifest_path.rename(real)
        saved.manifest_path.symlink_to(real)
        with self.assertRaises(ValueError):
            library.load_switch(saved.path)
        saved.manifest_path.unlink()
        real.rename(saved.manifest_path)
        backup_target = self.folder / "outside"
        backup_target.mkdir()
        (self.root / ".backups").symlink_to(backup_target, target_is_directory=True)
        with self.assertRaises(ValueError):
            self.save(overwrite=True)
        self.assertEqual(library.load_switch(saved.path), saved)
        self.assertFalse(list(backup_target.iterdir()))

    def test_malformed_and_duplicate_json_fields_are_ignored_by_listing_and_rejected_by_open(self):
        saved = self.save()
        original = saved.manifest_path.read_text()
        corruptions = ("{}", "{", original.replace('"version": 1,', '"version": 1, "version": 2,'))
        for payload in corruptions:
            saved.manifest_path.write_text(payload)
            with self.subTest(payload=payload[:50]), self.assertRaises(ValueError):
                library.load_switch(saved.path)
            self.assertEqual(library.list_saved_switches(self.root), [])

    def test_roots_are_design_local_or_export_local_and_missing_listing_is_readonly(self):
        self.assertEqual(library.saved_switch_root(self.folder / "active"), self.folder / "active" / library.LIBRARY_FOLDER)
        self.assertEqual(library.saved_switch_root(export_root=self.folder / "runs"), self.folder / "runs" / library.LIBRARY_FOLDER)
        self.assertEqual(library.saved_switch_root(export_root=self.root), self.root)
        self.assertEqual(library.list_saved_switches(self.root), [])
        self.assertFalse(self.root.exists())


if __name__ == "__main__":
    unittest.main()
