from __future__ import annotations

import importlib.util
import hashlib
import io
import os
from pathlib import Path
import subprocess
import sys
import tarfile
import tempfile
import unittest
from unittest.mock import patch
from types import SimpleNamespace


ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("release_science_builder", ROOT / "scripts/prepare_bundled_science.py")
science = importlib.util.module_from_spec(spec)
spec.loader.exec_module(science)


def load_builder(name, filename):
    specification = importlib.util.spec_from_file_location(name, ROOT / "scripts" / filename)
    module = importlib.util.module_from_spec(specification)
    specification.loader.exec_module(module)
    return module


vienna = load_builder("release_vienna_builder", "prepare_thermodynamics_runtime.py")
ostir = load_builder("release_ostir_builder", "prepare_ostir_sources.py")


class ReleaseBuildTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.work = Path(self.temporary.name)

    def environment_check(self, **overrides):
        environment = dict(os.environ)
        environment.pop("PYTHON_BIN", None)
        environment.pop("APTASWITCH_PYTHON", None)
        environment.update(overrides)
        return subprocess.run(
            ["bash", "-c", 'source "$1"; aptaswitch_prepare_build_python 1',
             "build-test", str(ROOT / "scripts/build_environment.sh")],
            cwd=self.work, env=environment, capture_output=True, text=True, timeout=10)

    @unittest.skipUnless(sys.version_info[:2] == (3, 11), "Native build requires CPython 3.11")
    def test_check_environment_uses_existing_python_without_creating_environment(self):
        result = self.environment_check(APTASWITCH_PYTHON=sys.executable)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(list(self.work.iterdir()), [])

    @unittest.skipUnless(sys.version_info[:2] == (3, 11), "Native build requires CPython 3.11")
    def test_explicit_python_takes_precedence_over_bootstrap(self):
        result = self.environment_check(PYTHON_BIN=sys.executable,
                                        APTASWITCH_PYTHON=str(self.work / "missing"))
        self.assertEqual(result.returncode, 0, result.stderr)

    def test_missing_explicit_environment_is_reported_without_creating_default(self):
        result = self.environment_check(PYTHON_BIN=str(self.work / "other-env/bin/python"))
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("PYTHON_BIN does not name an executable", result.stderr)
        self.assertEqual(list(self.work.iterdir()), [])

    def test_missing_bootstrap_is_reported_without_installing_anything(self):
        result = self.environment_check(APTASWITCH_PYTHON=str(self.work / "missing"))
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("Install CPython 3.11", result.stderr)
        self.assertEqual(list(self.work.iterdir()), [])

    def test_non_python_executable_is_rejected(self):
        result = self.environment_check(PYTHON_BIN="false")
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual(list(self.work.iterdir()), [])

    def test_tampered_vendored_scientific_source_cannot_trigger_download_or_build(self):
        source = self.work / "vendor/sources" / f"oxDNA-{science.OXDNA_REVISION}.tar.gz"
        source.parent.mkdir(parents=True)
        source.write_bytes(b"tampered source")
        with patch.object(science, "ROOT", self.work), patch.object(science, "download") as download, \
                patch.object(science.subprocess, "run") as build:
            with self.assertRaisesRegex(RuntimeError, "Vendored oxDNA source archive checksum"):
                science.prepare_oxdna(self.work / "cache", self.work / "output", 2)
        download.assert_not_called()
        build.assert_not_called()
        self.assertFalse((self.work / ".local_deps").exists())

    def test_authenticated_vendored_scientific_source_does_not_access_the_network(self):
        source = self.work / "vendor/sources" / f"oxDNA-{science.OXDNA_REVISION}.tar.gz"
        source.parent.mkdir(parents=True)
        source.write_bytes(b"authenticated source fixture")
        with patch.object(science, "ROOT", self.work), patch.object(science, "download") as download, \
                patch.object(science, "OXDNA_ARCHIVE_SHA256", hashlib.sha256(source.read_bytes()).hexdigest()):
            self.assertEqual(science.oxdna_source_archive(self.work / "cache"), source)
        download.assert_not_called()

    def fixture_archive(self, filename, member_name, data):
        archive = self.work / "vendor/sources" / filename
        archive.parent.mkdir(parents=True, exist_ok=True)
        with tarfile.open(archive, "w:gz") as output:
            member = tarfile.TarInfo(member_name)
            member.size = len(data)
            output.addfile(member, io.BytesIO(data))
        return archive

    def test_vienna_vendor_source_is_authenticated_and_cached_before_building(self):
        params = b"audited parameter source fixture\n"
        archive = self.fixture_archive(f"ViennaRNA-{vienna.VERSION}.tar.gz",
                                       f"ViennaRNA-{vienna.VERSION}/src/ViennaRNA/params/params.c", params)
        authenticated = archive.read_bytes()
        cache = self.work / "build-cache"
        with patch.object(vienna, "ROOT", self.work), \
                patch.object(vienna, "ARCHIVE_SHA256", hashlib.sha256(authenticated).hexdigest()), \
                patch.object(vienna, "PATCHED_PARAMS", hashlib.sha256(params).hexdigest()), \
                patch.object(vienna.urllib.request, "urlopen") as network, \
                patch.object(vienna.subprocess, "run") as build:
            source = vienna.prepare_source(cache)
            self.assertEqual((source / "src/ViennaRNA/params/params.c").read_bytes(), params)
            self.assertEqual((cache / archive.name).read_bytes(), authenticated)
            self.assertTrue(build.called)
            network.assert_not_called()
            archive.write_bytes(b"modified vendor source")
            build.reset_mock()
            with self.assertRaisesRegex(RuntimeError, "ViennaRNA source archive checksum"):
                vienna.prepare_source(cache)
            build.assert_not_called()
            self.assertEqual((cache / archive.name).read_bytes(), authenticated)

    def test_ostir_vendor_source_is_authenticated_without_network_or_installed_source_changes(self):
        archive = self.fixture_archive(f"ostir-{ostir.VERSION}.tar.gz",
                                       f"ostir-{ostir.VERSION}/LICENSE", b"OSTIR fixture license\n")
        authenticated = archive.read_bytes()
        vienna_cache = self.work / ".local_deps/thermodynamics-build/ViennaRNA-2.7.2.tar.gz"
        vienna_cache.parent.mkdir(parents=True)
        vienna_cache.write_bytes(b"authenticated Vienna source fixture")
        with patch.object(ostir, "SHA256", hashlib.sha256(authenticated).hexdigest()), \
                patch.object(ostir.importlib.metadata, "version", side_effect=lambda name: {"OSTIR": "1.1.3", "ViennaRNA": "2.7.2"}[name]), \
                patch.object(ostir.importlib.metadata, "distribution", return_value=SimpleNamespace(files=())), \
                patch.object(ostir.urllib.request, "urlopen") as network, \
                patch("runpy.run_path", return_value={"ARCHIVE_SHA256": hashlib.sha256(vienna_cache.read_bytes()).hexdigest()}):
            output = ostir.prepare_sources(self.work)
            self.assertEqual((output / archive.name).read_bytes(), authenticated)
            self.assertEqual((output / "OSTIR-LICENSE.txt").read_bytes(), b"OSTIR fixture license\n")
            self.assertTrue((output / "Apta2Switch-source.tar.gz").is_file())
            network.assert_not_called()
            archive.write_bytes(b"modified vendor source")
            with self.assertRaisesRegex(RuntimeError, "OSTIR source checksum"):
                ostir.prepare_sources(self.work)
            self.assertEqual((output / archive.name).read_bytes(), authenticated)

    def test_python_license_is_found_in_debian_documentation(self):
        version = f"python{sys.version_info.major}.{sys.version_info.minor}"
        expected = Path("/usr/share/doc") / version / "copyright"
        with patch.object(Path, "is_file", lambda path: path == expected):
            self.assertEqual(science.python_license_file(), expected)

    def test_missing_python_notice_is_not_silently_omitted(self):
        with patch.object(Path, "is_file", return_value=False):
            with self.assertRaisesRegex(RuntimeError, "complete license"):
                science.python_license_file()

    def test_nonpositive_parallelism_is_rejected_before_building(self):
        result = subprocess.run(
            [sys.executable, str(ROOT / "scripts/prepare_bundled_science.py"), "--jobs", "0"],
            cwd=self.work, capture_output=True, text=True, timeout=10)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("--jobs must be positive", result.stderr)


if __name__ == "__main__":
    unittest.main()
