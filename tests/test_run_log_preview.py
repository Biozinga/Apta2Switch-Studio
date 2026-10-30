"""Large optimizer records must not block saved-run display or hide later messages."""

from __future__ import annotations

import io
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from aptaswitch_core.localization import language_context, tr
from aptaswitch_studio.run_library_ui import (
    _LOG_PREVIEW_LIMIT,
    _bounded_log_lines,
    _log_preview,
)


class RunLogPreviewTests(unittest.TestCase):
    def test_ten_megabyte_optimization_keeps_following_messages_and_raw_file(self):
        payload = 'Optimization: {"trace": [' + '"x: y",' * 1_500_000 + 'null]}\n'
        raw = 'Début du calcul\n' + payload + 'Exports terminés.\n'
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "run.log"
            path.write_text(raw, encoding="utf-8")
            with patch("aptaswitch_studio.run_library_ui.tr", wraps=tr) as translate:
                with language_context("en"), path.open(encoding="utf-8") as stream:
                    preview, shortened = _log_preview(_bounded_log_lines(stream))
            self.assertTrue(shortened)
            self.assertLessEqual(len(preview), _LOG_PREVIEW_LIMIT)
            self.assertIn("Optimization:", preview)
            self.assertTrue(preview.endswith("Exports complete."))
            self.assertFalse(any("Optimization:" in call.args[0] for call in translate.call_args_list))
            self.assertEqual(path.read_text(encoding="utf-8"), raw)

    def test_short_json_is_opaque_but_messages_are_translated(self):
        raw = 'Thermodynamics: {"backend": "test"}\nExports terminés.\n'
        with patch("aptaswitch_studio.run_library_ui.tr", wraps=tr) as translate:
            with language_context("en"):
                preview, shortened = _log_preview(_bounded_log_lines(io.StringIO(raw)))
        self.assertFalse(shortened)
        self.assertEqual(preview, 'Thermodynamics: {"backend": "test"}\nExports complete.')
        translate.assert_called_once_with("Exports terminés.")

    def test_many_lines_keep_the_final_messages_with_bounded_output(self):
        raw = "\n".join(f"Message {index}: " + "x" * 100 for index in range(1000))
        raw += "\nExports terminés.\n"
        with language_context("fr"):
            preview, shortened = _log_preview(_bounded_log_lines(io.StringIO(raw)))
        self.assertTrue(shortened)
        self.assertLessEqual(len(preview), _LOG_PREVIEW_LIMIT)
        self.assertNotIn("Message 0:", preview)
        self.assertTrue(preview.endswith("Exports terminés."))

    def test_long_unterminated_line_and_multibyte_text(self):
        raw = "Résumé : " + "é" * 5000
        preview, shortened = _log_preview(_bounded_log_lines(io.StringIO(raw)))
        self.assertTrue(shortened)
        self.assertTrue(preview.startswith("Résumé : "))
        self.assertTrue(preview.endswith(" …"))
        self.assertNotIn("\ufffd", preview)


if __name__ == "__main__":
    unittest.main()
