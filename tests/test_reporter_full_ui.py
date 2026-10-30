"""The reporter library displays complete CDS and preserves custom additions."""

import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from streamlit.testing.v1 import AppTest
from aptaswitch_core import settings
from aptaswitch_core.models import REPORTERS


class ReporterFullUiTests(unittest.TestCase):
    def setUp(self):
        self.root = Path(self.enterContext(tempfile.TemporaryDirectory())).resolve()
        self.enterContext(patch.object(settings, "SETTINGS_DIR", self.root))
        self.enterContext(patch.object(settings, "SETTINGS_PATH", self.root / "settings.json"))


    def test_add_complete_reporter_persists_full_sequence_and_30_base_prefix(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary).resolve()
            self.enterContext(patch.object(settings, "SETTINGS_DIR", root))
            self.enterContext(patch.object(settings, "SETTINGS_PATH", root / "settings.json"))
            app = AppTest.from_string(f'''
import streamlit as st
from pathlib import Path
from aptaswitch_core import settings
from aptaswitch_core.localization import set_language
from aptaswitch_studio.web_app import _render_reporter_settings
set_language("fr")
settings.SETTINGS_DIR = Path({str(root)!r})
settings.SETTINGS_PATH = Path({str(root / 'settings.json')!r})
st.session_state.setdefault("reporter_return_to_design", False)
st.session_state.setdefault("reporter_key", "sfgfp")
st.session_state.setdefault("last_reporter_key", "sfgfp")
_render_reporter_settings()
''').run()
            self.assertFalse(app.exception)
            displayed = [element.value for element in app.code]
            self.assertEqual(len(displayed), 5)
            self.assertCountEqual(displayed, ["AUG" + reporter.full_sequence_after_start_rna
                                             for reporter in REPORTERS.values()])
            self.assertFalse(any(len(sequence) == 30 for sequence in displayed))
            self.assertFalse(any("bases utilisées pour le design" in item.value for item in app.caption))
            full = "ATG" + "ACG" * 15 + "TAA"
            app.text_input(key="new_reporter_label").input("Contexte complet")
            app.text_area(key="new_reporter_sequence").input(full)
            next(button for button in app.button if button.label == "Ajouter").click().run()
            self.assertFalse(app.exception)
            data = json.loads((root / "settings.json").read_text())["reporters"][0]
            self.assertEqual(data["sequence"], "ACG" * 10)
            self.assertEqual(data["full_sequence_after_start_rna"], "ACG" * 15 + "UAA")
            self.assertEqual(app.session_state["reporter_key"], data["key"])
            self.assertIn(full.replace("T", "U"), [element.value for element in app.code])
            self.assertFalse(any(len(element.value) == 30 for element in app.code))


if __name__ == "__main__":
    unittest.main()
