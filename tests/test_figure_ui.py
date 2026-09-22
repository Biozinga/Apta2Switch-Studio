"""Live molecular drawings stay clean while downloads keep their annotations."""

import base64
import re
import unittest
from unittest.mock import patch
import xml.etree.ElementTree as ET

from tests.ui_helpers import FrenchAppTest as AppTest

from aptaswitch_studio import figure_ui


STRUCTURE = "((((+))))"


def figure_app(linear=False):
    from functools import partial
    from aptaswitch_core.extension_preview import aptamer_trigger_linear_svg
    from aptaswitch_core.extension_visuals import render_extension_structure_svg
    from aptaswitch_studio.figure_ui import render_diagram

    common = dict(title="Système annoté", subtitle="ΔG = -4.25 kcal/mol", ligand_aptamer_positions=[2])
    diagram = (partial(aptamer_trigger_linear_svg, "AAGC", "GCTT", **common) if linear else
               partial(render_extension_structure_svg, "AAGC", "GCTT", "((((+))))", **common))
    render_diagram(diagram, structure="((((+))))", filename="molecule", key="molecule")


class FigureUiTests(unittest.TestCase):
    def test_both_views_keep_annotations_only_in_downloads_and_copy_the_exact_structure(self):
        for linear in (False, True):
            with self.subTest(linear=linear), patch.object(figure_ui, "svg_to_png_bytes", wraps=figure_ui.svg_to_png_bytes) as png:
                app = AppTest.from_function(figure_app, args=(linear,)).run()
            self.assertFalse(app.exception)
            html = app.get("html")[0].proto.body
            svg = base64.b64decode(re.search(r'data:image/svg\+xml;base64,([^" ]+)', html).group(1)).decode()
            live_text = " ".join(ET.fromstring(svg).itertext())
            self.assertNotIn("ΔG", live_text)
            self.assertNotIn("Système annoté", live_text)
            self.assertNotIn("((((", live_text)
            self.assertNotIn("))))", live_text)
            self.assertIn("5", live_text)
            exported = png.call_args.args[0]
            self.assertIn("ΔG = -4.25 kcal/mol", exported)
            self.assertIn("Système annoté", exported)
            if not linear:
                self.assertIn("(((( ))))", exported)
                self.assertIn('class="base-annotation" data-position="2"', svg)
            popovers = app.get("popover")
            self.assertEqual([p.proto.popover.label for p in popovers], ["Structure", "Exporter"])
            self.assertEqual([c.value for c in popovers[0].code], [STRUCTURE])
            self.assertEqual(len(popovers[1].get("download_button")), 2)
            self.assertEqual(len(app.get("download_button")), 2)


if __name__ == "__main__":
    unittest.main()
