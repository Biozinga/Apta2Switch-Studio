"""Original design structures stay exact and portable with a saved candidate."""

from copy import deepcopy
from dataclasses import replace
import json
import os
from pathlib import Path
import shutil
import tempfile
import unittest
from unittest.mock import patch

from aptaswitch_core.structure_ensemble import ensemble_inputs, ensemble_key, ensemble_path, load_ensemble, save_ensemble
from aptaswitch_studio import saved_switch_design as design
from aptaswitch_studio import saved_switch_library as library
from test_structure_ensemble import fake_analysis, fixture


class SavedSwitchDesignTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.folder = Path(temporary.name).resolve()
        self.root = self.folder / "design" / library.LIBRARY_FOLDER
        self.result = fixture(self.folder / "source")
        self.candidate = self.result.candidates[0]
        self.snapshot = library.candidate_snapshot(self.result, self.candidate, "Retained")
        self.default = self.folder / "run-library"
        self.legacy = self.folder / "previous-library"
        self.enterContext(patch.object(design, "default_export_dir", return_value=self.default))
        self.enterContext(patch.object(design, "legacy_export_dir", return_value=self.legacy))

    def artifact(self, count=1000, *, result=None, directory=None):
        result = result or self.result
        analysis = fake_analysis(result, count)
        path = save_ensemble(directory or result.output_dir, analysis)
        return analysis, path

    def save(self, **kwargs):
        return library.save_switch(self.root, self.snapshot["name"], self.snapshot, **kwargs)

    def test_save_detaches_all_sample_variants_and_source_timestamps(self):
        originals = {}
        for count, timestamp in ((100, 1), (1000, 10), (3000, 30), (5000, 20)):
            analysis, path = self.artifact(count)
            os.utime(path, ns=(timestamp * 1_000_000_000, timestamp * 1_000_000_000))
            originals[count] = analysis, path.stat().st_mtime_ns
        before = deepcopy(self.snapshot)
        saved = self.save(source_output_dir=self.result.output_dir)
        manifest = saved.manifest_path.read_bytes()
        shutil.rmtree(self.result.output_dir)
        moved_path = self.folder / "portable-candidate"
        saved.path.rename(moved_path)
        moved = library.load_switch(moved_path)
        with patch("aptaswitch_core.thermodynamics.setup.get_backend", side_effect=AssertionError("No calculation")):
            result, candidate = design.load_saved_design(moved)
        self.assertEqual(self.snapshot, before)
        self.assertEqual(moved.manifest_path.read_bytes(), manifest)
        self.assertEqual(result.output_dir, moved.path / "calculations" / "design")
        self.assertEqual(candidate.to_dict(), self.candidate.to_dict())
        self.assertEqual(result.request.thermo, self.result.request.thermo)
        self.assertEqual(result.request.reporter, self.result.request.reporter)
        self.assertEqual(result.warnings, self.result.warnings)
        for count, (analysis, timestamp) in originals.items():
            inputs = ensemble_inputs(result.request, candidate, count)
            self.assertEqual(load_ensemble(result.output_dir, inputs), analysis)
            self.assertEqual(ensemble_path(result.output_dir, analysis["key"]).stat().st_mtime_ns, timestamp)

    def test_existing_selection_recovers_only_its_own_design_run(self):
        source = self.root.parent / "calculations" / "designs" / self.result.run_id
        analysis, _ = self.artifact(directory=source)
        saved = self.save()
        before = saved.manifest_path.read_bytes()
        self.assertFalse((saved.path / "calculations" / "design").exists())
        result, candidate = design.load_saved_design(saved)
        self.assertEqual(load_ensemble(result.output_dir, ensemble_inputs(result.request, candidate)), analysis)
        self.assertEqual(saved.manifest_path.read_bytes(), before)
        shutil.rmtree(source)
        restored, candidate = design.load_saved_design(saved)
        self.assertEqual(load_ensemble(restored.output_dir, ensemble_inputs(restored.request, candidate)), analysis)

    def test_conventional_shared_run_libraries_can_restore_existing_selection(self):
        for root in (self.default, self.legacy):
            with self.subTest(root=root.name):
                analysis, path = self.artifact(directory=root / "designs" / self.result.run_id)
                self.snapshot = library.candidate_snapshot(self.result, self.candidate, root.name)
                saved = self.save()
                result, candidate = design.load_saved_design(saved)
                self.assertEqual(load_ensemble(result.output_dir, ensemble_inputs(result.request, candidate)), analysis)
                path.unlink()

    def test_missing_source_keeps_exact_mfe_and_detached_snapshot(self):
        self.candidate.inverse_extension_metrics = {"recorded": {"value": 1}}
        self.snapshot = library.candidate_snapshot(self.result, self.candidate, "Retained")
        saved = self.save(source_output_dir=self.folder / "absent")
        before = deepcopy(saved.snapshot)
        result, candidate = design.load_saved_design(saved)
        self.assertEqual(candidate.structure_off, self.candidate.structure_off)
        self.assertEqual(candidate.structure_on, self.candidate.structure_on)
        self.assertEqual(candidate.delta_g_off, self.candidate.delta_g_off)
        self.assertIsNone(load_ensemble(result.output_dir, ensemble_inputs(result.request, candidate)))
        candidate.inverse_extension_metrics["recorded"]["value"] = 9
        self.assertEqual(saved.snapshot, before)

    def test_paths_recorded_in_json_and_other_run_caches_are_never_followed(self):
        self.artifact()
        self.artifact(directory=self.root.parent / "calculations" / "designs" / "another_run")
        self.snapshot["source"]["output_dir"] = str(self.result.output_dir)
        self.snapshot["request"]["output_dir"] = str(self.result.output_dir)
        saved = self.save()
        result, candidate = design.load_saved_design(saved)
        self.assertIsNone(load_ensemble(result.output_dir, ensemble_inputs(result.request, candidate)))
        self.assertFalse((result.output_dir / "structure_ensembles").exists())

    def test_changed_sequence_conditions_or_broken_matrix_do_not_match(self):
        valid, _ = self.artifact(1000)
        foreign = deepcopy(self.result)
        foreign.candidates[0].switch_seq = "C" + foreign.candidates[0].switch_seq[1:]
        self.artifact(3000, result=foreign, directory=self.result.output_dir)
        foreign = deepcopy(self.result)
        foreign.request = replace(foreign.request, thermo=replace(foreign.request.thermo, temperature_c=42))
        foreign_analysis, _ = self.artifact(5000, result=foreign, directory=self.result.output_dir)
        forged = ensemble_path(self.result.output_dir, ensemble_key(ensemble_inputs(self.result.request, self.candidate, 5000)))
        forged.write_text(json.dumps(foreign_analysis))
        invalid, bad_path = self.artifact(100)
        invalid["states"]["off"]["pairs"] = [[0]]
        bad_path.write_text(json.dumps(invalid))
        saved = self.save(source_output_dir=self.result.output_dir)
        result, candidate = design.load_saved_design(saved)
        files = list((result.output_dir / "structure_ensembles").glob("*.json"))
        self.assertEqual([path.stem for path in files], [valid["key"]])
        self.assertEqual(load_ensemble(result.output_dir, ensemble_inputs(result.request, candidate)), valid)

    def test_valid_saved_families_are_not_overwritten_by_new_source_samples(self):
        source = self.root.parent / "calculations" / "designs" / self.result.run_id
        original, _ = self.artifact(directory=source)
        saved = self.save(source_output_dir=source)
        target = ensemble_path(saved.path / "calculations" / "design", original["key"])
        original_bytes, original_time = target.read_bytes(), target.stat().st_mtime_ns
        replacement = deepcopy(original)
        replacement["model"]["version"] = "later recomputation"
        save_ensemble(source, replacement)
        with patch.object(design, "load_ensemble", side_effect=AssertionError("Already detached")) as reread:
            design.load_saved_design(saved)
        reread.assert_not_called()
        self.assertEqual(target.read_bytes(), original_bytes)
        self.assertEqual(target.stat().st_mtime_ns, original_time)

    def test_foreign_candidate_caches_are_skipped_without_parsing(self):
        source = self.root.parent / "calculations" / "designs" / self.result.run_id
        foreign = deepcopy(self.result)
        foreign.request = replace(foreign.request, thermo=replace(foreign.request.thermo, temperature_c=42))
        self.artifact(result=foreign, directory=source)
        saved = self.save()
        with patch.object(design, "load_ensemble", side_effect=AssertionError("Foreign cache")) as read:
            result, _ = design.load_saved_design(saved)
        read.assert_not_called()
        self.assertFalse((result.output_dir / "structure_ensembles").exists())

    def test_failed_cache_copy_never_replaces_existing_candidate(self):
        saved = self.save()
        original = saved.manifest_path.read_bytes()
        self.artifact()
        with patch.object(design, "save_ensemble", side_effect=OSError("disk full")):
            with self.assertRaisesRegex(OSError, "disk full"):
                self.save(overwrite=True, source_output_dir=self.result.output_dir)
        self.assertEqual(saved.manifest_path.read_bytes(), original)
        self.assertFalse((self.root / ".backups").exists())
        self.assertFalse(list(self.root.glob(".switch-staging-*")))

    def test_linked_source_cache_is_not_followed(self):
        analysis, path = self.artifact()
        source = self.root.parent / "calculations" / "designs" / self.result.run_id
        source.mkdir(parents=True)
        (source / "structure_ensembles").symlink_to(path.parent, target_is_directory=True)
        saved = self.save()
        result, candidate = design.load_saved_design(saved)
        self.assertIsNone(load_ensemble(result.output_dir, ensemble_inputs(result.request, candidate)))
        with self.assertRaises(ValueError):
            self.save(overwrite=True, source_output_dir=source)
        self.assertEqual(json.loads(path.read_text()), analysis)

    def test_linked_destination_and_replaced_snapshot_are_rejected(self):
        analysis, source = self.artifact()
        saved = self.save()
        target = ensemble_path(saved.path / "calculations" / "design", analysis["key"])
        target.parent.mkdir(parents=True)
        target.symlink_to(source)
        with self.assertRaises(ValueError):
            design.load_saved_design(saved)
        target.unlink()
        current = json.loads(saved.manifest_path.read_text())
        current["candidate"]["structure_off"] = "edited"
        saved.manifest_path.write_text(json.dumps(current))
        with self.assertRaisesRegex(ValueError, "replaced"):
            design.load_saved_design(saved)


if __name__ == "__main__":
    unittest.main()
