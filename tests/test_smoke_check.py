"""Release UI checks isolate settings even after consumers have been imported."""
from pathlib import Path
from importlib.metadata import PackageNotFoundError, version
import tempfile
import unittest
from unittest.mock import patch

from aptaswitch_core import settings
from aptaswitch_core import reporter_library
from aptaswitch_studio import design_library, design_library_ui
from aptaswitch_studio.smoke_check import _check_rbs_runtime, _isolated_ui_settings


class SmokeSettingsTests(unittest.TestCase):
    def test_preimported_consumers_share_disposable_settings_and_restore_user_paths(self):
        root = Path(self.enterContext(tempfile.TemporaryDirectory())).resolve()
        remembered_path = root / "remembered" / "settings.json"
        remembered_path.parent.mkdir()
        remembered_path.write_text('{"private": "untouched", "design_library_dir": "/old/designs"}')
        original = remembered_path.read_bytes()
        with (
            patch.object(settings, "SETTINGS_DIR", remembered_path.parent),
            patch.object(settings, "SETTINGS_PATH", remembered_path),
        ):
            with _isolated_ui_settings(root / "smoke"):
                isolated = settings.load_settings()
                self.assertNotIn("private", isolated)
                for consumer in (reporter_library, design_library):
                    self.assertEqual(consumer.load_settings(), isolated)
                self.assertEqual(design_library_ui.default_design_root(), root / "smoke" / "designs")
                reporter_library.save_settings({**isolated, "smoke_write": True})
                self.assertTrue(settings.load_settings()["smoke_write"])
                self.assertEqual(remembered_path.read_bytes(), original)
            self.assertEqual(settings.SETTINGS_PATH, remembered_path)
            self.assertEqual(settings.SETTINGS_DIR, remembered_path.parent)
            self.assertEqual(remembered_path.read_bytes(), original)

    def test_settings_paths_are_restored_when_the_ui_check_fails(self):
        root = Path(self.enterContext(tempfile.TemporaryDirectory())).resolve()
        original = settings.SETTINGS_DIR, settings.SETTINGS_PATH
        with self.assertRaisesRegex(RuntimeError, "page failed"):
            with _isolated_ui_settings(root):
                raise RuntimeError("page failed")
        self.assertEqual((settings.SETTINGS_DIR, settings.SETTINGS_PATH), original)


class SmokeRBSRuntimeTests(unittest.TestCase):
    def setUp(self):
        try:
            available = version("OSTIR") == "1.1.3" and version("ViennaRNA") == "2.7.2"
        except PackageNotFoundError:
            available = False
        if not available:
            self.skipTest("Pinned OSTIR/ViennaRNA runtime unavailable")

    def test_real_runtime_and_corrupted_constraint_provenance(self):
        from aptaswitch_core import rbs_calculator

        rbs_calculator._prediction.cache_clear()
        checked = _check_rbs_runtime()
        self.assertGreater(checked["on_rate"], checked["off_rate"])
        self.assertEqual(checked["on_constraints"].count("x"), 6)
        predict = rbs_calculator.predict_rbs_switch

        def drop_mask(*args):
            result = predict(*args)
            result["on_constraints"] = "." * len(result["on_constraints"])
            return result

        with patch.object(rbs_calculator, "predict_rbs_switch", side_effect=drop_mask):
            with self.assertRaisesRegex(ValueError, "different on_constraints"):
                _check_rbs_runtime()


if __name__ == "__main__":
    unittest.main()
