from __future__ import annotations

import sys
import tempfile
import unittest
from datetime import datetime, timezone
from pathlib import Path
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from aptaswitch_core.models import RunRequest, SequenceInput, SwitchArchitecture
from aptaswitch_core.run_storage import (
    application_directory,
    default_export_dir,
    legacy_export_dir,
    new_run_directory,
    run_output_directory,
)


class RunStorageTests(unittest.TestCase):
    def test_source_library_is_in_project_even_when_launched_elsewhere(self):
        with (
            patch("aptaswitch_core.run_storage.sys.frozen", False, create=True),
            patch("aptaswitch_core.run_storage.Path.cwd", return_value=Path("/unrelated")),
        ):
            self.assertEqual(application_directory(), ROOT)
            self.assertEqual(default_export_dir(), ROOT / "runs")
            request = RunRequest(
                sequences=SequenceInput(aptamer_dna="ACGT", trigger_dna="ACGT"),
                architecture=SwitchArchitecture(toehold_length=2, stem_length=2),
            )
            self.assertEqual(request.output_dir, ROOT / "runs")

    def test_mac_application_exports_beside_signed_bundle(self):
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp).resolve()
            executable = directory / "Apta2Switch-Studio.app" / "Contents" / "MacOS" / "Apta2Switch-Studio"
            with (
                patch("aptaswitch_core.run_storage.sys.frozen", True, create=True),
                patch("aptaswitch_core.run_storage.sys.executable", str(executable)),
            ):
                self.assertEqual(application_directory(), directory)
                self.assertEqual(default_export_dir(), directory / "runs")

    def test_source_and_built_application_in_checkout_share_one_library(self):
        executable = ROOT / "dist" / "Apta2Switch-Studio.app" / "Contents" / "MacOS" / "Apta2Switch-Studio"
        with (
            patch("aptaswitch_core.run_storage.sys.frozen", True, create=True),
            patch("aptaswitch_core.run_storage.sys.executable", str(executable)),
        ):
            self.assertEqual(application_directory(), ROOT)
            self.assertEqual(default_export_dir(), ROOT / "runs")

    def test_external_app_does_not_use_an_unrelated_python_project(self):
        with tempfile.TemporaryDirectory() as tmp:
            project = Path(tmp).resolve()
            (project / "pyproject.toml").write_text("[project]\nname = 'unrelated'\n", encoding="utf-8")
            app_directory = project / "apps"
            executable = app_directory / "Apta2Switch-Studio.app" / "Contents" / "MacOS" / "Apta2Switch-Studio"
            with (
                patch("aptaswitch_core.run_storage.sys.frozen", True, create=True),
                patch("aptaswitch_core.run_storage.sys.executable", str(executable)),
            ):
                self.assertEqual(application_directory(), app_directory)

    def test_other_packaged_executable_exports_to_its_directory(self):
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp).resolve()
            with (
                patch("aptaswitch_core.run_storage.sys.frozen", True, create=True),
                patch("aptaswitch_core.run_storage.sys.executable", str(directory / "Apta2Switch-Studio")),
            ):
                self.assertEqual(application_directory(), directory)

    def test_design_and_extension_share_a_root_with_separate_categories(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            started = datetime(2026, 9, 22, 14, 5, 6, tzinfo=timezone.utc)
            design_id, design_dir = new_run_directory(root, "design", started)
            extension_id, extension_dir = new_run_directory(root, "extension", started)
            self.assertRegex(design_id, r"^20260922_140506_design_[a-f0-9]{6}$")
            self.assertRegex(extension_id, r"^20260922_140506_extension_[a-f0-9]{6}$")
            self.assertEqual(design_dir, root / "designs" / design_id)
            self.assertEqual(extension_dir, root / "extensions" / extension_id)
            self.assertFalse(design_dir.exists())
            self.assertNotEqual(new_run_directory(root, "design", started)[0], design_id)

    def test_category_selection_does_not_nest_category_folders(self):
        root = Path("/custom/library")
        self.assertEqual(run_output_directory(root / "designs", "design"), root / "designs")
        self.assertEqual(run_output_directory(root / "extensions", "extension"), root / "extensions")
        self.assertEqual(run_output_directory(root / "designs", "extension"), root / "extensions")
        self.assertEqual(run_output_directory(root / "extensions", "design"), root / "designs")
        with self.assertRaises(ValueError):
            new_run_directory(root, "unknown")

    def test_previous_library_remains_addressable(self):
        self.assertEqual(legacy_export_dir(Path("/user/home")), Path("/user/home/Apta2Switch-Studio Runs"))


if __name__ == "__main__":
    unittest.main()
