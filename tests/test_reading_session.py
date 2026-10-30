"""Reading projects retain their own sequences and named-library mode."""

from __future__ import annotations

from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from streamlit.testing.v1 import AppTest

from aptaswitch_studio.design_library import create_design, load_design
from aptaswitch_studio.design_session import (
    open_design_into_session, reset_design_session, save_active_design,
    validate_design_state,
)


LIBRARY_SCRIPT = """
import streamlit as st
from aptaswitch_core.localization import set_language
from aptaswitch_studio.design_library_ui import render_design_library
set_language('fr')
render_design_library()
"""


class ReadingSessionTests(unittest.TestCase):
    def setUp(self):
        self.root = Path(self.enterContext(tempfile.TemporaryDirectory())).resolve()
        self.project = create_design(self.root, "Reading project")
        self.enterContext(patch("aptaswitch_studio.design_session.all_reporters", return_value={}))

    def test_reading_form_roundtrip_keeps_all_three_sequences_and_reporter_choice(self):
        session = {}
        open_design_into_session(self.project.path, session)
        expected = {
            "design_mode": "reading", "aptamer_input": "AACCGG",
            "trigger_input": "CCGGTT", "reading_switch_input": "gggACGU\nACGU",
            "reading_use_full_reporter": True, "reporter_key": "custom_reporter",
        }
        session.update(expected)
        saved = save_active_design(session)
        self.assertEqual({key: saved.state["form"][key] for key in expected}, expected)
        reopened = {}
        open_design_into_session(self.project.path, reopened)
        self.assertEqual({key: reopened[key] for key in expected}, expected)

    def test_reading_fields_reset_when_switching_to_another_design(self):
        session = {}
        open_design_into_session(self.project.path, session)
        session.update(design_mode="reading", reading_switch_input="AUGACGU",
                       reading_use_full_reporter=False)
        other = create_design(self.root, "Design project", state={"form": {"design_mode": "design"}})
        open_design_into_session(other.path, session)
        self.assertEqual(session["design_mode"], "design")
        self.assertNotIn("reading_switch_input", session)
        self.assertNotIn("reading_use_full_reporter", session)
        open_design_into_session(self.project.path, session)
        self.assertEqual(session["design_mode"], "reading")
        self.assertEqual(session["reading_switch_input"], "AUGACGU")
        self.assertIs(session["reading_use_full_reporter"], False)
        reset_design_session(session)
        self.assertTrue(all(key not in session for key in (
            "design_mode", "reading_switch_input", "reading_use_full_reporter",
        )))

    def test_reading_drafts_allow_partial_sequences_and_validate_field_types(self):
        validate_design_state({"form": {"design_mode": "reading", "reading_switch_input": "",
                                        "reading_use_full_reporter": False}})
        validate_design_state({"form": {"design_mode": "design"}})
        validate_design_state({"form": {}})  # Existing projects need no migration.
        for field, value in (
            ("design_mode", "unknown"), ("design_mode", True),
            ("reading_switch_input", ["AUG"]),
            ("reading_use_full_reporter", "yes"), ("reading_use_full_reporter", 1),
        ):
            with self.subTest(field=field, value=value), self.assertRaises(ValueError):
                validate_design_state({"form": {field: value}})


class ReadingLibraryUITests(unittest.TestCase):
    def setUp(self):
        self.root = Path(self.enterContext(tempfile.TemporaryDirectory())).resolve()
        self.enterContext(patch("aptaswitch_studio.design_library_ui.default_design_root", return_value=self.root))
        self.enterContext(patch("aptaswitch_studio.design_session.all_reporters", return_value={}))

    def app(self):
        app = AppTest.from_string(LIBRARY_SCRIPT, default_timeout=30).run()
        self.assertFalse(app.exception)
        return app

    def test_creation_chooses_reading_mode_and_displays_it_in_library(self):
        app = self.app()
        selector = app.radio(key="design_new_mode")
        self.assertEqual(selector.value, "design")
        self.assertEqual(selector.options, ["Conception", "Lecture de séquences"])
        selector.set_value("reading")
        app.text_input(key="design_new_name").set_value("Imported switch")
        app.button(key="create_design").click().run()
        self.assertFalse(app.exception)
        self.assertFalse(app.error)
        saved = load_design(Path(app.session_state["active_design_path"]))
        self.assertEqual(saved.name, "Imported switch")
        self.assertEqual(saved.state["form"]["design_mode"], "reading")
        self.assertEqual(app.session_state["design_mode"], "reading")
        self.assertEqual(app.session_state["nav"], "Conception")
        self.assertTrue(any(item.value == "Lecture de séquences · Switch : 0 nt" for item in app.caption))

    def test_default_creation_retains_design_mode(self):
        app = self.app()
        app.text_input(key="design_new_name").set_value("Normal design")
        app.button(key="create_design").click().run()
        self.assertFalse(app.exception)
        saved = load_design(Path(app.session_state["active_design_path"]))
        self.assertEqual(saved.state["form"]["design_mode"], "design")


if __name__ == "__main__":
    unittest.main()
