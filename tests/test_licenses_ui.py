"""The ready-to-use settings expose credits, never engine installation forms."""
import unittest
from unittest.mock import patch

from streamlit.testing.v1 import AppTest

from tests.ui_helpers import FrenchAppTest
from aptaswitch_studio.licenses_ui import CREDITS


class LicensesUITests(unittest.TestCase):
    def test_bilingual_credits_have_primary_references_and_actual_licences(self):
        for language, title in (('fr', 'Licences et remerciements'),
                                ('en', 'Licences and acknowledgements')):
            with self.subTest(title=title):
                app = AppTest.from_string(
                    'from aptaswitch_core.localization import set_language\n'
                    f'set_language({language!r})\n'
                    'from aptaswitch_studio.licenses_ui import render_licenses\n'
                    'render_licenses()'
                ).run()
                self.assertFalse(app.exception)
                self.assertEqual(app.subheader[0].value, title)
                self.assertFalse(app.text_input)
                self.assertFalse(app.button)
                text = '\n'.join(x.value for x in app.caption)
                self.assertIn('GPLv3', text)
                self.assertIn('OSTIR 1.1.3', text)
                self.assertIn('MIT', text)
                self.assertIn('BSD', text)
                for _, _, _, links in CREDITS:
                    for _, url in links:
                        self.assertTrue(any(url in item.value for item in app.markdown), url)

    def test_clean_settings_need_no_engine_registration_and_keep_custom_reporters(self):
        with patch('aptaswitch_studio.web_app.load_settings', return_value={}), \
             patch('aptaswitch_core.reporter_library.load_settings', return_value={}):
            app = FrenchAppTest.from_string('''
from aptaswitch_studio.web_app import _initialize_session, _render_settings_page
_initialize_session()
_render_settings_page()
''').run()
        self.assertFalse(app.exception)
        self.assertEqual([tab.label for tab in app.tabs], ['Licences et remerciements', 'Reporters', 'Général'])
        self.assertTrue(any(area.key == 'new_reporter_sequence' for area in app.text_area))
        self.assertTrue(any(button.label == 'Ajouter' for button in app.button))
        self.assertFalse(any('path' in (field.key or '') for field in app.text_input))
        self.assertFalse(any('Valider' in b.label or 'Vérifier' in b.label for b in app.button))


if __name__ == '__main__':
    unittest.main()
