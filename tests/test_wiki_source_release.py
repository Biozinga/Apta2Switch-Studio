from __future__ import annotations

import hashlib
import importlib.util
import io
from pathlib import Path
import tarfile
import tempfile
import unittest
from unittest.mock import patch


ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("wiki_source_packager", ROOT / "scripts/package_wiki_source.py")
packager = importlib.util.module_from_spec(spec)
spec.loader.exec_module(packager)
AUDITED_UPSTREAM = packager.UPSTREAM


class WikiSourceReleaseTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.work = Path(temporary.name)
        self.root = self.work / "source"
        self.root.mkdir()
        self.output = self.work / "release.tar.gz"
        self.upstream = []
        for original in AUDITED_UPSTREAM:
            entry = dict(original)
            path = self.root / "vendor/sources" / entry["filename"]
            path.parent.mkdir(parents=True, exist_ok=True)
            with tarfile.open(path, "w:gz") as archive:
                notice = (entry["name"] + " fixture license\n").encode()
                member = tarfile.TarInfo(entry["notice"])
                member.size = len(notice)
                archive.addfile(member, io.BytesIO(notice))
            entry["sha256"] = hashlib.sha256(path.read_bytes()).hexdigest()
            self.upstream.append(entry)
        override = patch.object(packager, "UPSTREAM", tuple(self.upstream))
        override.start()
        self.addCleanup(override.stop)
        for name in sorted(packager.REQUIRED_FILES):
            self.write(name, "public source or notice\n")
        self.write("pyproject.toml", '[project]\nversion = "0.2.1b2"\n')
        self.write("scripts/prepare_thermodynamics_runtime.py",
                   f'VERSION = {self.upstream[0]["version"]!r}\nARCHIVE_SHA256 = {self.upstream[0]["sha256"]!r}\n')
        self.write("scripts/prepare_bundled_science.py",
                   f'OXDNA_REVISION = {self.upstream[1]["version"]!r}\nOXDNA_ARCHIVE_SHA256 = {self.upstream[1]["sha256"]!r}\n')
        self.write("scripts/prepare_ostir_sources.py",
                   f'VERSION = {self.upstream[2]["version"]!r}\nSHA256 = {self.upstream[2]["sha256"]!r}\n')

    def write(self, name, text):
        path = self.root / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(text)
        return path

    def build(self):
        with patch.object(packager.urllib.request, "urlopen", side_effect=AssertionError("Unexpected network access")):
            return packager.build(self.root, self.output, offline=True)

    def rewrite(self, *, change=None, extra=None):
        changed = self.work / "changed.tar.gz"
        with tarfile.open(self.output) as original, tarfile.open(changed, "w:gz") as archive:
            for member in original.getmembers():
                data = original.extractfile(member).read()
                if change and member.name.endswith("/" + change):
                    data += b"modified"
                    member.size = len(data)
                archive.addfile(member, io.BytesIO(data))
            if extra is not None:
                archive.addfile(extra, io.BytesIO(b"x") if extra.isfile() else None)
        return changed

    def test_roundtrip_contains_sources_not_local_data_and_keeps_executable_scripts(self):
        self.write("designs/private/design.json", "private design")
        self.write("runs/private.json", "private run")
        self.write("src/designs/private.json", "nested private design")
        self.write("src/.env.local", "private token")
        self.write("src/.aws/credentials", "private credential")
        self.write("src/__pycache__/cache.pyc", "derived cache")
        self.write("tests/private/results/private.json", "private calculation")
        self.write("tests/native_reference/results/private.json", "private calculation")
        self.write("tests/validation-results/private.json", "private calculation")
        self.write("src/module.py", "public module")
        (self.root / "scripts/build_linux.sh").chmod(0o755)
        manifest = self.build()
        self.assertEqual(packager.check(self.output), manifest)
        self.assertIn("src/module.py", manifest["files"])
        self.assertFalse(any("private" in name or ".env" in name or "__pycache__" in name
                             or "credentials" in name for name in manifest["files"]))
        self.assertEqual(len(manifest["upstream"]), 3)
        with tarfile.open(self.output) as archive:
            self.assertEqual(archive.getmember("Apta2Switch-Studio-0.2.1b2/scripts/build_linux.sh").mode, 0o755)
        self.assertEqual(self.output.with_name(self.output.name + ".sha256").read_text().split()[0],
                         hashlib.sha256(self.output.read_bytes()).hexdigest())

    def test_source_selection_never_follows_top_level_or_nested_symlinks(self):
        private = self.work / "private"
        private.mkdir()
        (private / "secret.py").write_text("private source")
        (self.root / "README.md").unlink()
        (self.root / "README.md").symlink_to(private / "secret.py")
        (self.root / "src/nested.py").symlink_to(private / "secret.py")
        (self.root / "packaging").symlink_to(private, target_is_directory=True)
        selected = {path.relative_to(self.root).as_posix() for path in packager.selected_files(self.root)}
        self.assertNotIn("README.md", selected)
        self.assertNotIn("src/nested.py", selected)
        self.assertNotIn("packaging/secret.py", selected)

    def test_tampered_source_is_rejected_even_when_archive_structure_is_valid(self):
        self.build()
        with self.assertRaisesRegex(RuntimeError, "differs from manifest"):
            packager.check(self.rewrite(change="src/aptaswitch_core/thermodynamics/native/connected_dp.cpp"))

    def test_unsafe_archive_paths_and_links_are_rejected(self):
        self.build()
        for name, kind in (("Apta2Switch-Studio-0.2.1b2/../outside.py", tarfile.REGTYPE),
                           ("Apta2Switch-Studio-0.2.1b2/src/link.py", tarfile.SYMTYPE),
                           ("Apta2Switch-Studio-0.2.1b2/src/./ambiguous.py", tarfile.REGTYPE),
                           ("Apta2Switch-Studio-0.2.1b2/src//ambiguous.py", tarfile.REGTYPE),
                           ("/absolute/file.py", tarfile.REGTYPE)):
            member = tarfile.TarInfo(name)
            member.type = kind
            if member.isfile():
                member.size = 1
            else:
                member.linkname = "../../private"
            with self.subTest(name=name), self.assertRaisesRegex(RuntimeError, "Unexpected source archive entry"):
                packager.check(self.rewrite(extra=member))

    def test_duplicate_archive_entries_are_rejected(self):
        self.build()
        member = tarfile.TarInfo("Apta2Switch-Studio-0.2.1b2/README.md")
        member.size = 1
        with self.assertRaisesRegex(RuntimeError, "Duplicate source archive path"):
            packager.check(self.rewrite(extra=member))

    def test_missing_source_in_offline_mode_never_downloads_or_replaces_output(self):
        (self.root / "vendor/sources" / self.upstream[0]["filename"]).unlink()
        self.output.write_bytes(b"keep previous release")
        with self.assertRaisesRegex(RuntimeError, "Missing upstream source in offline mode"):
            self.build()
        self.assertEqual(self.output.read_bytes(), b"keep previous release")

    def test_changed_builder_pins_fail_before_an_archive_is_written(self):
        self.write("scripts/prepare_bundled_science.py", "OXDNA_REVISION = 'different'\nOXDNA_ARCHIVE_SHA256 = 'different'\n")
        with self.assertRaisesRegex(RuntimeError, "pins differ from the actual builder"):
            self.build()
        self.assertFalse(self.output.exists())

    def test_actual_release_pins_match_the_three_production_builders(self):
        with patch.object(packager, "UPSTREAM", AUDITED_UPSTREAM):
            packager.check_build_pins(ROOT)

    def test_checkout_contains_required_build_inputs_and_reference_fixtures(self):
        selected = {path.relative_to(ROOT).as_posix() for path in packager.selected_files(ROOT)}
        packager.check_required_files(selected)

    def test_missing_native_builder_fails_before_upstream_downloads(self):
        (self.root / "scripts/build_hybrid_native.py").unlink()
        with patch.object(packager, "upstream_archive") as upstream:
            with self.assertRaisesRegex(RuntimeError, "scripts/build_hybrid_native.py"):
                self.build()
        upstream.assert_not_called()
        self.assertFalse(self.output.exists())


if __name__ == "__main__":
    unittest.main()
