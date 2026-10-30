"""A named design must remain portable and survive rejected imports or saves."""

from __future__ import annotations

import io
import json
import stat
import tempfile
import unittest
import zipfile
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from unittest.mock import patch

from aptaswitch_studio import design_library as library


class DesignLibraryTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name).resolve()
        self.library = self.root / "designs"

    def create(self, name="Mon aptamère", **kwargs):
        return library.create_design(self.library, name, **kwargs)

    def archive(self, members):
        output = io.BytesIO()
        with zipfile.ZipFile(output, "w", compression=zipfile.ZIP_DEFLATED) as archive:
            for name, content in members:
                archive.writestr(name, content)
        return output.getvalue()

    def manifest(self, name="Mon aptamère"):
        return json.dumps({
            "format": library.FORMAT, "version": 1, "name": name,
            "created_at": "2026-09-23T10:00:00+00:00",
            "updated_at": "2026-09-23T10:00:00+00:00", "state": {},
        })

    def test_create_and_save_preserve_identity_and_do_not_mutate_input(self):
        state = {"aptamer_dna": "ACGT", "trigger_dna": "TGCA", "selected": [1, 2]}
        project = self.create(state=state)
        state["selected"].append(3)
        self.assertEqual(project.name, "Mon aptamère")
        self.assertEqual(project.path, self.library / project.name)
        self.assertEqual(project.state["selected"], [1, 2])
        self.assertTrue((project.path / "sequences").is_dir())
        self.assertTrue(library.design_output_dir(project).is_dir())
        self.assertEqual(library.design_output_dir(project.path), project.path / "calculations")
        saved = library.save_design(project.path, {"aptamer_dna": "AACCGGTT"})
        self.assertEqual(saved.created_at, project.created_at)
        self.assertGreaterEqual(saved.updated_at, project.updated_at)
        self.assertEqual(library.load_design(project.path / "design.json"), saved)
        self.assertEqual(library.list_designs(self.library), [saved])

    def test_default_root_keeps_old_run_library_separate(self):
        with patch.object(library, "load_settings", return_value={}), patch.object(library, "default_export_dir", return_value=self.root / "runs"):
            self.assertEqual(library.default_design_root(), self.root / "designs")
        with patch.object(library, "load_settings", return_value={"design_library_dir": str(self.root / "custom")}):
            self.assertEqual(library.default_design_root(), self.root / "custom")

    def test_unchanged_autosave_preserves_timestamp_and_manifest_bytes(self):
        project = self.create(state={"aptamer_dna": "ACGT", "selected": [2, 4]})
        before = (project.path / "design.json").read_bytes()
        with patch.object(library, "_write_manifest", side_effect=AssertionError("unchanged state must not write")):
            saved = library.save_design(project.path, {"selected": [2, 4], "aptamer_dna": "ACGT"})
        self.assertEqual(saved, project)
        self.assertEqual((project.path / "design.json").read_bytes(), before)

    def test_casefold_and_unicode_normalization_collisions_need_explicit_overwrite(self):
        original = self.create("ÉTUDE", state={"keep": True})
        marker = original.path / "calculations" / "valuable.pdb"
        marker.write_text("ATOMS")
        for name in ("étude", "e\u0301tude", " ÉTUDE "):
            with self.subTest(name=name), self.assertRaises(library.DesignExistsError):
                self.create(name)
            self.assertEqual(library.load_design(original.path).state, {"keep": True})
            self.assertEqual(marker.read_text(), "ATOMS")
        self.assertFalse((self.library / ".backups").exists())

    def test_overwrite_keeps_old_folder_as_hidden_backup(self):
        original = self.create(state={"old": True})
        (original.path / "calculations" / "model.pdb").write_text("original")
        fresh = self.create("MON APTAMÈRE", overwrite=True, state={"new": True})
        self.assertEqual(fresh.state, {"new": True})
        self.assertFalse((fresh.path / "calculations" / "model.pdb").exists())
        backups = list((self.library / ".backups").iterdir())
        self.assertEqual(len(backups), 1)
        self.assertEqual(library.load_design(backups[0]).state, {"old": True})
        self.assertEqual((backups[0] / "calculations" / "model.pdb").read_text(), "original")
        self.assertEqual(library.list_designs(self.library), [fresh])

    def test_failed_publication_restores_old_folder(self):
        original = self.create(state={"old": True})
        real_rename = Path.rename

        def fail_staging_rename(path, target):
            if path.name.startswith(".design-staging-"):
                raise OSError("simulated full disk")
            return real_rename(path, target)

        with patch.object(Path, "rename", fail_staging_rename), self.assertRaises(OSError):
            self.create(overwrite=True, state={"new": True})
        self.assertEqual(library.load_design(original.path).state, {"old": True})
        self.assertFalse(list(self.library.glob(".design-staging-*")))
        self.assertEqual(list((self.library / ".backups").iterdir()), [])

    def test_invalid_names_cannot_escape_library_or_use_hidden_names(self):
        for name in ("", " ", ".", "..", "../outside", "a/b", "a\\b", "/absolute", ".hidden", "name\x00", "name\nnext", "C:drive", "ending.", "NUL"):
            with self.subTest(name=name), self.assertRaises(ValueError):
                self.create(name)
        self.assertFalse(self.library.exists())

    def test_invalid_state_and_atomic_write_failure_preserve_last_manifest(self):
        original = self.create(state={"saved": True})
        before = (original.path / "design.json").read_bytes()
        for state in ([], {"invalid": object()}, {"nan": float("nan")}):
            with self.subTest(state=state), self.assertRaises(ValueError):
                library.save_design(original.path, state)
            self.assertEqual((original.path / "design.json").read_bytes(), before)
        with patch.object(library.os, "replace", side_effect=OSError("simulated full disk")), self.assertRaises(OSError):
            library.save_design(original.path, {"new": True})
        self.assertEqual((original.path / "design.json").read_bytes(), before)
        self.assertFalse(list(original.path.glob(".design-*.json")))

    def test_list_skips_broken_manifests_hidden_folders_and_symlinks(self):
        project = self.create()
        broken = self.library / "broken"
        broken.mkdir()
        (broken / "design.json").write_text("not JSON")
        (self.library / "linked").symlink_to(project.path, target_is_directory=True)
        hidden = self.library / ".hidden"
        hidden.mkdir()
        (hidden / "design.json").write_text(self.manifest())
        self.assertEqual(library.list_designs(self.library), [project])
        self.assertEqual(library.list_designs(self.root / "missing"), [])

    def test_zip_roundtrip_preserves_inputs_results_and_md_files_at_new_location(self):
        project = self.create(state={"inputs": {"aptamer": "ACGT", "trigger": "GCGC"}, "result": "calculations/extensions/run/results.json"})
        sequences = project.path / "sequences" / "original.fasta"
        sequences.write_text(">aptamer\nACGT\n>trigger\nGCGC\n")
        results = project.path / "calculations" / "extensions" / "run"
        results.mkdir(parents=True)
        (results / "results.json").write_text('{"extended": "AAGCGCCC"}')
        models = results / "molecular_dynamics"
        models.mkdir()
        (models / "last_conf.dat").write_text("actual coarse-grained coordinates\n")
        (models / "viewer.html").write_text("<html>offline viewer</html>")
        archive = library.export_design(project.path)
        with zipfile.ZipFile(io.BytesIO(archive)) as zipped:
            self.assertIn("design.json", zipped.namelist())
            self.assertTrue(all(not name.startswith("/") for name in zipped.namelist()))
        restored = library.import_design(archive, self.root / "other computer")
        self.assertEqual(restored.state, project.state)
        self.assertEqual(restored.created_at, project.created_at)
        self.assertEqual((restored.path / "sequences" / "original.fasta").read_bytes(), sequences.read_bytes())
        self.assertEqual((restored.path / "calculations/extensions/run/molecular_dynamics/last_conf.dat").read_text(), "actual coarse-grained coordinates\n")
        self.assertNotIn(str(project.path), (restored.path / "design.json").read_text())

    def test_local_folder_manifest_and_wrapped_zip_import(self):
        original = self.create()
        for source in (original.path, original.path / "design.json"):
            imported = library.import_design(source, self.root / "copies", name=source.name)
            self.assertEqual(imported.state, original.state)
        wrapped = self.archive([
            ("Export/design.json", self.manifest("Renamed import")),
            ("Export/calculations/model.pdb", "ATOM"),
            ("__MACOSX/ignored.txt", "metadata"),
        ])
        imported = library.import_design(wrapped, self.library)
        self.assertEqual(imported.name, "Renamed import")
        self.assertEqual((imported.path / "calculations/model.pdb").read_text(), "ATOM")

    def test_export_and_local_copy_skip_symlinks_and_external_engines(self):
        original = self.create()
        private = self.root / "external-private-file"
        private.write_text("must not ship")
        (original.path / "calculations" / "linked.pdb").symlink_to(private)
        (original.path / "external_user_site").mkdir()
        (original.path / "external_user_site" / "engine").write_text("restricted")
        (original.path / "external").symlink_to(private.parent, target_is_directory=True)
        for copied in (library.import_design(original.path, self.root / "local-copy"), library.import_design(library.export_design(original.path), self.root / "zip-copy")):
            self.assertFalse((copied.path / "external_user_site").exists())
            self.assertFalse((copied.path / "external").exists())
            self.assertFalse((copied.path / "calculations/linked.pdb").exists())

    def test_import_collision_is_non_destructive_until_confirmed(self):
        original = self.create(state={"original": True})
        archive = self.archive([("design.json", self.manifest("MON APTAMÈRE")), ("calculations/new.pdb", "new")])
        with self.assertRaises(library.DesignExistsError):
            library.import_design(archive, self.library)
        self.assertEqual(library.load_design(original.path).state, {"original": True})
        imported = library.import_design(archive, self.library, overwrite=True)
        self.assertEqual((imported.path / "calculations/new.pdb").read_text(), "new")
        self.assertEqual(len(list((self.library / ".backups").iterdir())), 1)

    def test_malformed_import_with_overwrite_cannot_touch_existing_design(self):
        original = self.create(state={"valuable": True})
        before = (original.path / "design.json").read_bytes()
        malformed = [
            b"not a zip",
            self.archive([("design.json", "not JSON")]),
            self.archive([("design.json", self.manifest()), ("../escaped", "bad")]),
            self.archive([("design.json", self.manifest()), ("/absolute", "bad")]),
            self.archive([("design.json", self.manifest()), ("C:/drive", "bad")]),
            self.archive([("design.json", self.manifest()), ("sub\\escape", "bad")]),
            self.archive([("design.json", self.manifest()), ("data/file", "first"), ("DATA/FILE", "second")]),
            self.archive([("design.json", self.manifest()), ("data", "file"), ("data/file", "nested")]),
            self.archive([("one/design.json", self.manifest()), ("two/design.json", self.manifest())]),
            self.archive([("one/design.json", self.manifest()), ("outside", "bad")]),
        ]
        for archive in malformed:
            with self.subTest(size=len(archive)), self.assertRaises(ValueError):
                library.import_design(archive, self.library, overwrite=True)
            self.assertEqual((original.path / "design.json").read_bytes(), before)
            self.assertFalse(list(self.library.glob(".design-staging-*")))
        self.assertFalse((self.root / "escaped").exists())
        self.assertFalse((self.library / ".backups").exists())

    def test_symlink_zip_member_is_rejected(self):
        link = zipfile.ZipInfo("calculations/link")
        link.create_system = 3
        link.external_attr = (stat.S_IFLNK | 0o777) << 16
        archive = self.archive([("design.json", self.manifest()), (link, "/private/native")])
        with self.assertRaises(ValueError):
            library.import_design(archive, self.library)
        self.assertEqual(library.list_designs(self.library), [])

    def test_archive_file_count_and_expanded_size_limits_apply_before_publish(self):
        archive = self.archive([("design.json", self.manifest()), ("data/large", "X" * 3000)])
        for variable, limit in (("MAX_FILES", 1), ("MAX_FILE_BYTES", 2000), ("MAX_TOTAL_BYTES", 2000)):
            with self.subTest(variable=variable), patch.object(library, variable, limit), self.assertRaises(ValueError):
                library.import_design(archive, self.library)
        self.assertEqual(library.list_designs(self.library), [])

    def test_unknown_schema_and_invalid_dates_are_rejected_before_overwrite(self):
        original = self.create()
        for field, value in (("format", "other"), ("version", 2), ("version", True), ("created_at", "yesterday"), ("updated_at", "2026-09-23"), ("state", [])):
            data = json.loads(self.manifest())
            data[field] = value
            archive = self.archive([("design.json", json.dumps(data))])
            with self.subTest(field=field, value=value), self.assertRaises(ValueError):
                library.import_design(archive, self.library, overwrite=True)
        self.assertEqual(library.load_design(original.path), original)

    def test_invalid_known_session_state_is_rejected_before_overwrite(self):
        original = self.create(state={"valuable": True})
        data = json.loads(self.manifest())
        data["state"] = {"form": "invalid form shape"}
        archive = self.archive([("design.json", json.dumps(data))])
        with self.assertRaises(ValueError):
            library.import_design(archive, self.library, overwrite=True)
        self.assertEqual(library.load_design(original.path), original)
        self.assertFalse((self.library / ".backups").exists())

    def test_concurrent_creation_does_not_silently_overwrite(self):
        def create(index):
            try:
                return self.create(state={"thread": index})
            except library.DesignExistsError:
                return None

        with ThreadPoolExecutor(max_workers=4) as pool:
            results = list(pool.map(create, range(4)))
        winners = [project for project in results if project is not None]
        self.assertEqual(len(winners), 1)
        self.assertEqual(library.list_designs(self.library), winners)


if __name__ == "__main__":
    unittest.main()
