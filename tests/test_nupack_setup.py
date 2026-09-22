from __future__ import annotations

import json
import subprocess
import sys
import tempfile
import unittest
import zipfile
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

from packaging.tags import Tag

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from aptaswitch_core import nupack_setup as setup


class PortableNupackTests(unittest.TestCase):
    tag = Tag("cp311", "cp311", "macosx_11_0_arm64")

    def make_wheel(self, folder: Path, version: str = "4.1.0.1", extra: dict | None = None) -> Path:
        folder.mkdir(parents=True, exist_ok=True)
        wheel = folder / f"nupack-{version}-{self.tag}.whl"
        with zipfile.ZipFile(wheel, "w") as archive:
            archive.writestr("nupack/__init__.py", f"__version__ = '{version}'")
            archive.writestr("nupack/parameters/rna.json", "{}")
            for path, content in (extra or {}).items():
                archive.writestr(path, content)
        return wheel

    def test_official_download_selects_newest_compatible_wheel(self):
        with tempfile.TemporaryDirectory() as temporary:
            folder = Path(temporary)
            current = self.make_wheel(folder / "nupack-4.1.0.1" / "package")
            older = self.make_wheel(folder / "package", "4.0.2.1")
            (folder / "nupack-9.0-cp310-cp310-macosx_11_0_arm64.whl").touch()
            (folder / "nupack-9.0-cp311-cp311-linux_x86_64.whl").touch()
            (folder / "other-9.0-cp311-cp311-macosx_11_0_arm64.whl").touch()
            with patch.object(setup, "sys_tags", return_value=[self.tag]):
                self.assertEqual(setup.compatible_nupack_wheels(folder), [current, older])

    def test_official_folder_installs_without_pip_or_host_python(self):
        with tempfile.TemporaryDirectory() as temporary:
            folder = Path(temporary)
            download = folder / "download"
            self.make_wheel(download / "package")
            with (
                patch.object(setup, "SETTINGS_DIR", folder / "settings"),
                patch.object(setup, "sys_tags", return_value=[self.tag]),
                patch.object(setup, "_python_import_check", return_value=setup.NupackValidationResult(True, "", "OK", "4.1.0.1")),
                patch.object(setup.subprocess, "run") as run,
            ):
                result = setup.validate_nupack_path(str(download))
                self.assertTrue(result.ok, result.message)
                self.assertEqual(result.path, str(download.resolve()))
                self.assertTrue((Path(result.import_path) / "nupack" / "parameters" / "rna.json").is_file())
                run.assert_not_called()

    def test_incompatible_wheel_is_rejected_before_import(self):
        with tempfile.TemporaryDirectory() as temporary:
            wheel = self.make_wheel(Path(temporary))
            with (
                patch.object(setup, "sys_tags", return_value=[Tag("cp312", "cp312", "linux_x86_64")]),
                patch.object(setup, "_python_import_check") as check,
            ):
                result = setup.validate_nupack_wheel(wheel)
                self.assertFalse(result.ok)
                self.assertIn("compatible", result.message)
                check.assert_not_called()

    def test_wheel_cannot_extract_outside_its_cache(self):
        with tempfile.TemporaryDirectory() as temporary:
            folder = Path(temporary)
            wheel = self.make_wheel(folder, extra={"../../outside.txt": "bad"})
            with (
                patch.object(setup, "SETTINGS_DIR", folder / "settings"),
                patch.object(setup, "sys_tags", return_value=[self.tag]),
                patch.object(setup, "_python_import_check") as check,
            ):
                result = setup.validate_nupack_wheel(wheel)
                self.assertFalse(result.ok)
                check.assert_not_called()
                self.assertEqual(list(folder.rglob("outside.txt")), [])

    def test_failed_new_wheel_keeps_previous_working_cache(self):
        with tempfile.TemporaryDirectory() as temporary:
            folder = Path(temporary)
            wheel = self.make_wheel(folder)
            with (
                patch.object(setup, "SETTINGS_DIR", folder / "settings"),
                patch.object(setup, "sys_tags", return_value=[self.tag]),
                patch.object(setup, "_python_import_check", return_value=setup.NupackValidationResult(True, "", "OK")) as check,
            ):
                original = setup.validate_nupack_wheel(wheel)
                self.make_wheel(folder, extra={"nupack/broken.py": "fail"})
                check.return_value = setup.NupackValidationResult(False, "", "Bad binary")
                replacement = setup.validate_nupack_wheel(wheel)
                self.assertFalse(replacement.ok)
                self.assertTrue(Path(original.import_path).is_dir())

    def test_frozen_validation_dispatches_the_embedded_probe(self):
        def run(command, **kwargs):
            self.assertEqual(command[1], setup.PROBE_FLAG)
            self.assertNotIn("-m", command)
            self.assertNotIn("-c", command)
            Path(command[-1]).write_text(json.dumps({"ok": True, "version": "4.1"}))
            return subprocess.CompletedProcess(command, 0, "", "")

        with (
            patch.object(setup.sys, "frozen", True, create=True),
            patch.object(setup.subprocess, "run", side_effect=run),
        ):
            result = setup._python_import_check(Path("/user/nupack"))
            self.assertTrue(result.ok)
            self.assertEqual(result.version, "4.1")

    def test_probe_checks_actual_structure_calculation(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            calls = []
            nupack = SimpleNamespace(
                __file__=str(root / "nupack" / "__init__.py"), __version__="4.1",
                Model=lambda **kwargs: kwargs, config=SimpleNamespace(threads=0),
                mfe=lambda **kwargs: calls.append(kwargs) or ["(((...)))"],
            )
            result_file = root / "result.json"
            with patch.dict(sys.modules, {"nupack": nupack}), patch.object(sys, "path", list(sys.path)):
                self.assertTrue(setup.handle_nupack_probe([setup.PROBE_FLAG, str(root), str(result_file)]))
            self.assertTrue(json.loads(result_file.read_text())["ok"])
            self.assertEqual(calls[0]["strands"], ["GGGAAACCC"])

    def test_wrong_installed_package_cannot_make_empty_folder_validate(self):
        with tempfile.TemporaryDirectory() as temporary:
            with patch.object(setup, "_python_import_check") as check:
                self.assertFalse(setup.validate_nupack_path(temporary).ok)
                check.assert_not_called()

    def test_timed_out_binary_returns_validation_error(self):
        with patch.object(setup.subprocess, "run", side_effect=subprocess.TimeoutExpired("probe", 45)):
            self.assertFalse(setup._python_import_check(Path("/user/nupack")).ok)

    def test_drop_folder_uses_application_location_not_launch_directory(self):
        with tempfile.TemporaryDirectory() as temporary:
            folder = Path(temporary)
            with patch.object(setup, "application_directory", return_value=folder):
                self.assertEqual(setup.nupack_drop_directory(create=True), folder / "nupack")
                self.assertIn(folder / "nupack", setup.default_nupack_search_paths())

    def test_read_only_application_uses_user_folder(self):
        with tempfile.TemporaryDirectory() as temporary:
            folder = Path(temporary)
            with (
                patch.object(setup, "application_directory", return_value=folder / "readonly"),
                patch.object(setup, "SETTINGS_DIR", folder / "settings"),
                patch.object(setup.os, "access", return_value=False),
            ):
                self.assertEqual(setup.nupack_drop_directory(), folder / "settings" / "nupack")


if __name__ == "__main__":
    unittest.main()
