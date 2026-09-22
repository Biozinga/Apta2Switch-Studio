from __future__ import annotations

import tempfile
import unittest
import zipfile
import struct
from io import BytesIO
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from aptaswitch_studio.web_runtime import (
    concentration_from_m,
    concentration_to_m,
    contiguous_ranges,
    default_export_dir,
    DesignRunJob,
    ensure_writable_export_dir,
    exports_zip_bytes,
    largest_unbound_window,
    normalized_export_dir,
    preferred_concentration_unit,
    svg_to_png_bytes,
)


class WebRuntimeTests(unittest.TestCase):
    def test_svg_is_exported_as_high_resolution_png(self):
        svg = (
            '<svg xmlns="http://www.w3.org/2000/svg" width="20" height="10">'
            '<rect width="20" height="10" fill="#fff"/></svg>'
        )
        png = svg_to_png_bytes(svg)
        width, height = struct.unpack(">II", png[16:24])

        self.assertTrue(png.startswith(b"\x89PNG\r\n\x1a\n"))
        self.assertEqual((width, height), (40, 20))

    def test_design_job_snapshot_exposes_log_elapsed_time_and_eta(self):
        with patch("aptaswitch_studio.web_runtime.time.monotonic", side_effect=[100.0, 110.0]):
            job = DesignRunJob(SimpleNamespace(trials=4))
            job._progress(1, 4, "Moteur · candidat 1 terminé.")
            snapshot = job.snapshot()

        self.assertEqual(snapshot.elapsed_seconds, 10.0)
        self.assertEqual(snapshot.eta_seconds, 30.0)
        self.assertEqual(snapshot.current, 1)
        self.assertIn("Moteur · candidat 1 terminé.", snapshot.log[-1])

    def test_concentration_units_convert_to_molar_and_back(self):
        for value, unit, expected_m in (
            (5.0, "nM", 5e-9),
            (5.0, "µM", 5e-6),
            (120.0, "mM", 0.12),
            (0.5, "M", 0.5),
        ):
            self.assertAlmostEqual(concentration_to_m(value, unit), expected_m)
            self.assertAlmostEqual(concentration_from_m(expected_m, unit), value)

    def test_ascii_and_greek_micro_aliases_are_supported(self):
        self.assertAlmostEqual(concentration_to_m(5, "uM"), 5e-6)
        self.assertAlmostEqual(concentration_to_m(5, "μM"), 5e-6)

    def test_preferred_concentration_unit_is_readable(self):
        self.assertEqual(preferred_concentration_unit(5e-9), "nM")
        self.assertEqual(preferred_concentration_unit(5e-6), "µM")
        self.assertEqual(preferred_concentration_unit(0.12), "mM")
        self.assertEqual(preferred_concentration_unit(2.0), "M")

    def test_blank_and_relative_export_paths_are_moved_under_app_dir(self):
        with tempfile.TemporaryDirectory() as tmp:
            app_dir = Path(tmp)
            self.assertEqual(default_export_dir(app_dir), app_dir / "runs")
            self.assertEqual(normalized_export_dir("", app_dir), default_export_dir(app_dir))
            self.assertEqual(
                normalized_export_dir("my-run", app_dir),
                default_export_dir(app_dir) / "my-run",
            )
            self.assertEqual(normalized_export_dir(app_dir / "custom", app_dir), app_dir / "custom")

    def test_read_only_style_export_failure_uses_fallback(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            invalid = root / "not-a-directory"
            invalid.write_text("occupied", encoding="utf-8")
            fallback = root / "fallback"
            selected, used_fallback = ensure_writable_export_dir(invalid, fallback=fallback)
            self.assertEqual(selected, fallback)
            self.assertTrue(used_fallback)

    def test_read_only_app_directory_falls_back_to_previous_user_library(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            invalid_app = root / "read-only-app"
            invalid_app.write_text("occupied", encoding="utf-8")
            default = invalid_app / "runs"
            user_library = root / "user" / "Apta2Switch-Studio Runs"
            with (
                patch("aptaswitch_studio.web_runtime.default_export_dir", return_value=default),
                patch("aptaswitch_studio.web_runtime.legacy_export_dir", return_value=user_library),
            ):
                selected, used_fallback = ensure_writable_export_dir(default)
            self.assertEqual(selected, user_library)
            self.assertTrue(used_fallback)
            self.assertTrue(selected.is_dir())

    def test_writable_custom_library_is_preserved(self):
        with tempfile.TemporaryDirectory() as tmp:
            custom = Path(tmp) / "custom"
            selected, used_fallback = ensure_writable_export_dir(custom)
            self.assertEqual(selected, custom)
            self.assertFalse(used_fallback)

    def test_contiguous_ranges_are_collapsed(self):
        self.assertEqual(contiguous_ranges([5, 1, 2, 2, 4]), [(1, 2), (4, 5)])

    def test_largest_unbound_window(self):
        self.assertEqual(largest_unbound_window(10, [0, 1, 5, 6]), (2, 3))
        self.assertIsNone(largest_unbound_window(3, [0, 2]))

    def test_exports_are_bundled_without_directory_paths(self):
        with tempfile.TemporaryDirectory() as tmp:
            first = Path(tmp) / "result.csv"
            second = Path(tmp) / "result.json"
            first.write_text("rank,score\n1,95\n", encoding="utf-8")
            second.write_text('{"score": 95}', encoding="utf-8")
            payload = exports_zip_bytes({"csv": str(first), "json": str(second)})
            with zipfile.ZipFile(BytesIO(payload)) as archive:
                self.assertEqual(set(archive.namelist()), {"result.csv", "result.json"})


if __name__ == "__main__":
    unittest.main()
