from __future__ import annotations

import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import Mock, patch

from aptaswitch_core import platform_paths as paths
from aptaswitch_core import run_storage
from aptaswitch_studio import desktop


class PlatformPathTests(unittest.TestCase):
    def test_wsl_detection_requires_linux_and_accepts_kernel_fallback(self):
        with patch.dict(os.environ, {"WSL_DISTRO_NAME": "Apta2Switch"}, clear=True):
            with patch.object(sys, "platform", "darwin"):
                self.assertFalse(paths.is_wsl())
            with patch.object(sys, "platform", "linux"):
                self.assertTrue(paths.is_wsl())
        with patch.dict(os.environ, {}, clear=True), patch.object(sys, "platform", "linux"):
            with patch.object(Path, "read_text", return_value="6.6.87.2-microsoft-standard-WSL2"):
                self.assertTrue(paths.is_wsl())
            with patch.object(Path, "read_text", return_value="6.8.0-generic"):
                self.assertFalse(paths.is_wsl())

    def test_posix_paths_are_not_reinterpreted(self):
        with patch.object(paths, "is_wsl", return_value=False), patch.object(paths.subprocess, "run") as run:
            for value in ("relative folder", "/tmp/has spaces", r"C:\a literal POSIX name"):
                self.assertEqual(paths.normalize_user_path(value), Path(value))
            run.assert_not_called()

    def test_drive_paths_are_converted_as_one_argument_without_shell(self):
        value = r'C:\Users\Some Name\a $(touch sentinel); folder'
        with patch.object(paths, "is_wsl", return_value=True), patch.object(
            paths.subprocess, "run", return_value=Mock(stdout="/mnt/c/Users/Some Name/tools\n")
        ) as run:
            self.assertEqual(paths.normalize_user_path('"' + value + '"'), Path("/mnt/c/Users/Some Name/tools"))
        self.assertEqual(run.call_args.args[0], ["wslpath", "-u", value])
        self.assertNotIn("shell", run.call_args.kwargs)
        self.assertEqual(run.call_args.kwargs["timeout"], 5)

    def test_only_current_distribution_unc_is_accepted(self):
        with patch.object(paths, "is_wsl", return_value=True), patch.dict(
            os.environ, {"WSL_DISTRO_NAME": "Ubuntu-22.04"}
        ), patch.object(paths.subprocess, "run", return_value=Mock(stdout="/home/aptaswitch-studio/models\n")) as run:
            for host in ("wsl.localhost", "wsl$"):
                raw = "\\\\" + host + r"\Ubuntu-22.04\home\aptaswitch-studio\models"
                self.assertEqual(paths.normalize_user_path(raw), Path("/home/aptaswitch-studio/models"))
            paths.normalize_user_path(r"\\WSL.LOCALHOST\UBUNTU-22.04\home\aptaswitch-studio\models")
            self.assertEqual(run.call_args.args[0][-1], r"\\wsl.localhost\Ubuntu-22.04\home\aptaswitch-studio\models")
            run.reset_mock()
            for raw in (r"\\server\share\engine", r"\\wsl.localhost\AnotherDistro\home\model", "C:relative"):
                with self.subTest(raw=raw), self.assertRaises(ValueError):
                    paths.normalize_user_path(raw)
            run.assert_not_called()

    def test_wsl_system_converter_does_not_load_bundled_libraries(self):
        with patch.dict(os.environ, {"LD_LIBRARY_PATH": "/app/bundled", "LD_LIBRARY_PATH_ORIG": "/system/lib"}), patch.object(
            sys, "frozen", True, create=True
        ), patch.object(paths.subprocess, "run", return_value=Mock(stdout="/mnt/c/models")) as run:
            paths.wsl_path(r"C:\models")
            self.assertEqual(run.call_args.kwargs["env"]["LD_LIBRARY_PATH"], "/system/lib")
            self.assertEqual(os.environ["LD_LIBRARY_PATH"], "/app/bundled")

    def test_failed_conversion_never_becomes_a_relative_project_path(self):
        with patch.object(paths, "is_wsl", return_value=True):
            for outcome in (Mock(stdout=""), Mock(stdout="C:\\still-windows"),
                            Mock(stdout="/mnt/c/one\n/mnt/c/two")):
                with patch.object(paths.subprocess, "run", return_value=outcome), self.assertRaises(ValueError):
                    paths.normalize_user_path(r"C:\models")
            with patch.object(paths.subprocess, "run", side_effect=subprocess.TimeoutExpired("wslpath", 5)):
                with self.assertRaises(ValueError):
                    paths.normalize_user_path(r"C:\models")



class WslDesktopTests(unittest.TestCase):
    def test_folder_opens_windows_explorer_with_converted_path(self):
        unc = r"\\wsl.localhost\Ubuntu-22.04\home\user\Saved designs"
        with patch.object(desktop, "is_wsl", return_value=True), patch.object(
            desktop, "wsl_path", return_value=unc
        ) as convert, patch.object(desktop.subprocess, "Popen") as launch:
            desktop.open_directory("/home/user/Saved designs")
        convert.assert_called_once_with("/home/user/Saved designs", to_windows=True)
        self.assertEqual(launch.call_args.args[0], ["explorer.exe", unc])
        self.assertNotIn("shell", launch.call_args.kwargs)

    def test_browser_uses_windows_handler_and_reports_missing_interop(self):
        with patch.object(desktop, "is_wsl", return_value=True), patch.object(desktop.subprocess, "Popen") as launch:
            self.assertTrue(desktop.open_browser("http://127.0.0.1:8771"))
        self.assertEqual(launch.call_args.args[0], ["explorer.exe", "http://127.0.0.1:8771"])
        with patch.object(desktop, "is_wsl", return_value=True), patch.object(
            desktop.subprocess, "Popen", side_effect=FileNotFoundError
        ):
            self.assertFalse(desktop.open_browser("http://127.0.0.1:8771"))


class PersistentDataDirectoryTests(unittest.TestCase):
    def test_data_directory_survives_versioned_executable_location(self):
        from aptaswitch_core.oxdna_setup import default_oxdna_directory
        from aptaswitch_studio.design_library import default_design_root

        with tempfile.TemporaryDirectory() as directory, patch.dict(os.environ, {"APTASWITCH_DATA_DIR": directory}), patch.object(
            sys, "frozen", True, create=True
        ), patch("aptaswitch_studio.design_library.load_settings", return_value={}):
            for version in ("release-one", "release-two"):
                with patch.object(sys, "executable", "/opt/" + version + "/Apta2Switch-Studio"):
                    self.assertEqual(run_storage.application_directory(), Path(directory).resolve())
                    self.assertEqual(run_storage.default_export_dir(), Path(directory).resolve() / "runs")
                    self.assertEqual(default_design_root(), Path(directory).resolve() / "designs")
                    self.assertEqual(default_oxdna_directory(), Path("/opt/" + version + "/oxdna_runtime"))

    def test_relative_override_is_rejected_instead_of_depending_on_launch_directory(self):
        with patch.dict(os.environ, {"APTASWITCH_DATA_DIR": "relative/data"}), self.assertRaisesRegex(ValueError, "absolute"):
            run_storage.application_directory()


if __name__ == "__main__":
    unittest.main()
