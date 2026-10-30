"""Removing the sidebar monitor does not affect simulations across navigation."""
from __future__ import annotations

from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import Mock, patch

from aptaswitch_studio import saved_switch_ffs_runtime as runtime
from aptaswitch_studio import web_app
from aptaswitch_studio.web_runtime import JobSnapshot
from tests.ui_helpers import FrenchAppTest as AppTest


SCRIPT = """
import streamlit as st
from aptaswitch_studio.web_app import _initialize_session, _render_sidebar
from aptaswitch_studio.design_session import attach_saved_switch_jobs
_initialize_session()
attach_saved_switch_jobs(st.session_state)
page = _render_sidebar()
st.write(page)
if page == 'Toehold switches enregistrés':
    st.code(st.session_state.get('saved_switch_selected', ''))
"""


class FFSNavigationTests(unittest.TestCase):
    def setUp(self):
        self.folder = Path(self.enterContext(tempfile.TemporaryDirectory())).resolve()
        self.saved = SimpleNamespace(name='Retained second candidate',
                                     path=self.folder / 'tswitch_enregistres' / 'Second')
        self.state = {'steps': 100}
        self.job = SimpleNamespace(saved=self.saved, cancel=Mock(), snapshot=lambda: JobSnapshot(
            state='running', current=3, total=9,
            message=f"FFS · contact_2 · {self.state['steps']} pas", elapsed_seconds=120, log=()))
        self.enterContext(patch.dict(runtime._JOBS, {str(self.saved.path): self.job}, clear=True))
        self.enterContext(patch.object(web_app, 'persist_active_design', return_value=True))

    def app(self):
        app = AppTest.from_string(SCRIPT)
        app.session_state['active_design_path'] = str(self.folder)
        app.session_state['active_design_name'] = 'Test design'
        app.session_state['nav'] = 'Toehold switches enregistrés'
        app.session_state['saved_switch_selected'] = str(self.saved.path)
        app.run()
        self.assertFalse(app.exception)
        return app

    def test_navigation_preserves_job_and_candidate_without_sidebar_monitor(self):
        app = self.app()
        self.assertIs(app.session_state['saved_switch_ffs_job'], self.job)
        for page in ('Réglages', 'Conception', 'Mes designs'):
            self.state['steps'] += 100
            app.radio(key='nav').set_value(page).run()
            self.assertFalse(app.exception)
            self.assertIs(app.session_state['saved_switch_ffs_job'], self.job)
            self.assertEqual(app.session_state['saved_switch_selected'], str(self.saved.path))
            self.assertFalse(any(f"{self.state['steps']} pas" in caption.value for caption in app.sidebar.caption))
            self.assertFalse(any('Le calcul continue' in caption.value for caption in app.sidebar.caption))
            self.assertNotIn('sidebar_ffs_open', {button.key for button in app.sidebar.button})
        self.job.cancel.assert_not_called()
        app.radio(key='nav').set_value('Toehold switches enregistrés').run()
        self.assertFalse(app.exception)
        self.assertEqual(app.session_state['nav'], 'Toehold switches enregistrés')
        self.assertEqual(app.code[0].value, str(self.saved.path))
        self.job.cancel.assert_not_called()

    def test_new_browser_session_reattaches_without_a_sidebar_shortcut(self):
        app = self.app()
        del app.session_state['saved_switch_ffs_job']
        app.radio(key='nav').set_value('Réglages').run()
        self.assertFalse(app.exception)
        self.assertIs(app.session_state['saved_switch_ffs_job'], self.job)
        self.assertNotIn('sidebar_ffs_open', {button.key for button in app.sidebar.button})
        app.radio(key='nav').set_value('Toehold switches enregistrés').run()
        self.assertFalse(app.exception)
        self.assertEqual(app.session_state['saved_switch_selected'], str(self.saved.path))
        self.assertEqual(app.code[0].value, str(self.saved.path))
        self.job.cancel.assert_not_called()


if __name__ == '__main__':
    unittest.main()
