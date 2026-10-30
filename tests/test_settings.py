"""Shared settings remain independent of any scientific engine or installation."""
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from aptaswitch_core import settings


class SettingsTests(unittest.TestCase):
    def setUp(self):
        self.root = Path(self.enterContext(tempfile.TemporaryDirectory())) / "settings"
        self.enterContext(patch.object(settings, "SETTINGS_DIR", self.root))
        self.enterContext(patch.object(settings, "SETTINGS_PATH", self.root / "settings.json"))

    def test_missing_and_malformed_settings_have_no_implicit_runtime(self):
        self.assertEqual(settings.load_settings(), {})
        self.assertFalse(self.root.exists())
        self.root.mkdir()
        settings.SETTINGS_PATH.write_text("not json")
        self.assertEqual(settings.load_settings(), {})

    def test_existing_sections_and_custom_metadata_round_trip_unchanged(self):
        original = {"language": "fr", "reporters": [{"label": "custom"}],
                    "custom_plugin": {"path": "/custom/tool", "enabled": True},
                    "thermodynamics": {"library_path": "/native/lib"}}
        settings.save_settings(original)
        self.assertEqual(settings.load_settings(), original)
        changed = settings.load_settings()
        changed["language"] = "en"
        settings.save_settings(changed)
        self.assertEqual(settings.load_settings(), {**original, "language": "en"})


if __name__ == "__main__":
    unittest.main()
