from __future__ import annotations

import base64
import sys
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

from tests.ui_helpers import FrenchAppTest as AppTest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from aptaswitch_studio.web_app import (
    SIDEBAR_LOGO_PATH,
    _convert_concentration_display_unit,
    _inject_styles,
    _render_svg,
)
from aptaswitch_studio.web_runtime import DesignRunJob
from aptaswitch_core.architecture import default_architecture_for_trigger
from tests.design_fixtures import design_result_fixture
from aptaswitch_core.models import RunRequest, SequenceInput
from aptaswitch_core.preview import parse_structure


class WebAppRenderingTests(unittest.TestCase):
    def test_results_use_static_table_and_offer_png_structure(self):
        request = RunRequest(
            sequences=SequenceInput(
                aptamer_dna="ACGTACGT",
                trigger_dna="GTCCAGGCTGGTATAATTAGATCC",
            ),
            architecture=default_architecture_for_trigger(24),
            trials=1,
        )
        app = AppTest.from_file(
            str(ROOT / "src" / "aptaswitch_studio" / "web_app.py"),
            default_timeout=20,
        ).run()
        app.session_state["last_result"] = design_result_fixture(request)

        app.radio(key="nav").set_value("Résultats").run()

        self.assertEqual(len(app.get("table")), 1)
        self.assertEqual(len(app.get("dataframe")), 0)
        download_labels = {button.label for button in app.get("download_button")}
        self.assertIn("Télécharger la structure PNG", download_labels)
        self.assertIn("Séquence linéaire PNG", download_labels)
        self.assertFalse(any("linéaire indisponible" in error.value for error in app.error))
        self.assertFalse(app.exception)

        candidate = app.session_state["last_result"].candidates[0]
        for state, expected in (("OFF", candidate.structure_off), ("ON + trigger", candidate.structure_on)):
            app.radio(key="candidate_state").set_value(state).run()
            self.assertFalse(app.exception)
            popovers = app.get("popover")
            structures = [p for p in popovers if p.proto.popover.label == "Structure"]
            self.assertEqual(len(structures), 2)
            self.assertTrue(all([c.value for c in p.code] == [expected] for p in structures))
            parse_structure(expected)
            self.assertEqual(sum(len(p.get("download_button")) for p in popovers), len(app.get("download_button")))

    def test_design_preview_offers_png_download(self):
        app = AppTest.from_file(
            str(ROOT / "src" / "aptaswitch_studio" / "web_app.py"),
            default_timeout=20,
        ).run()
        app.text_area(key="aptamer_input").set_value("ACGTACGT")
        app.text_area(key="trigger_input").set_value("GTCCAGGCTGGTATAATTAGATCC").run()

        download_labels = {button.label for button in app.get("download_button")}
        self.assertIn("Télécharger cet aperçu PNG", download_labels)
        self.assertFalse(app.exception)

    def test_results_show_live_engine_console_and_timing(self):
        app = AppTest.from_file(
            str(ROOT / "src" / "aptaswitch_studio" / "web_app.py"),
            default_timeout=20,
        ).run()
        job = DesignRunJob(SimpleNamespace(trials=10, engine="nupack"))
        job._progress(2, 10, "NUPACK · essai 2/10 terminé.")
        app.session_state["design_job"] = job

        app.radio(key="nav").set_value("Résultats").run()

        metrics = {metric.label: metric.value for metric in app.metric}
        self.assertEqual(metrics["Moteur"], "NUPACK")
        self.assertEqual(metrics["Progression"], "2 / 10")
        self.assertIn("Temps écoulé", metrics)
        self.assertIn("Temps restant estimé", metrics)
        self.assertTrue(
            any("NUPACK · essai 2/10 terminé." in code.value for code in app.code)
        )
        self.assertFalse(app.exception)

    def test_design_values_survive_page_navigation(self):
        app = AppTest.from_file(
            str(ROOT / "src" / "aptaswitch_studio" / "web_app.py"),
            default_timeout=20,
        ).run()

        app.text_area(key="aptamer_input").set_value("ACGTACGT")
        app.text_area(key="trigger_input").set_value(
            "CCGATGCTCTCCTTACGCCACCCACACCCG"
        )
        app.text_input(key="molecule_name").set_value("test persistance")
        app.number_input(key="trigger_concentration_value").set_value(12.5)
        app.number_input(key="switch_concentration_value").set_value(7.25)
        app.run()

        app.radio(key="nav").set_value("Résultats").run()
        app.radio(key="nav").set_value("Conception").run()

        self.assertEqual(app.text_area(key="aptamer_input").value, "ACGTACGT")
        self.assertEqual(
            app.text_area(key="trigger_input").value,
            "CCGATGCTCTCCTTACGCCACCCACACCCG",
        )
        self.assertEqual(
            app.text_input(key="molecule_name").value,
            "test persistance",
        )
        self.assertEqual(
            app.number_input(key="trigger_concentration_value").value,
            12.5,
        )
        self.assertEqual(
            app.number_input(key="switch_concentration_value").value,
            7.25,
        )
        self.assertFalse(app.exception)

    def test_changing_concentration_unit_preserves_molar_value(self):
        state = {
            "species_value": 5.0,
            "species_unit": "nM",
            "species_previous_unit": "µM",
        }
        with patch("aptaswitch_studio.web_app.st.session_state", state):
            _convert_concentration_display_unit(
                "species_value",
                "species_unit",
                "species_previous_unit",
            )

        self.assertEqual(state["species_value"], 5000.0)
        self.assertEqual(state["species_previous_unit"], "nM")

    def test_sidebar_igem_logo_is_available(self):
        self.assertEqual(SIDEBAR_LOGO_PATH.name, "igem_SU_logo.png")
        self.assertTrue(SIDEBAR_LOGO_PATH.is_file())

    def test_sidebar_logo_has_a_white_backdrop(self):
        with patch("aptaswitch_studio.web_app.st.markdown") as render_markdown:
            _inject_styles()

        stylesheet = render_markdown.call_args.args[0]
        self.assertIn('[data-testid="stSidebar"] [data-testid="stImage"]', stylesheet)
        self.assertIn('background:#ffffff', stylesheet)

    def test_svg_preview_is_embedded_as_an_image(self):
        svg = '<svg xmlns="http://www.w3.org/2000/svg" width="20" height="10"><circle cx="5" cy="5" r="3"/></svg>'

        with patch("aptaswitch_studio.web_app.st.html") as render_html:
            _render_svg(svg, height=120, natural_width=640, scroll=True)

        document = render_html.call_args.args[0]
        encoded = base64.b64encode(svg.encode("utf-8")).decode("ascii")
        self.assertIn(f'data:image/svg+xml;base64,{encoded}', document)
        self.assertIn('<img ', document)
        self.assertIn('width:640px', document)
        self.assertIn('overflow-x:auto', document)
        self.assertIn('max-height:100%', document)
        self.assertIn('object-fit:contain', document)
        self.assertNotIn('<circle', document)

    def test_regular_preview_fits_both_dimensions(self):
        svg = '<svg xmlns="http://www.w3.org/2000/svg" width="1000" height="500" />'

        with patch("aptaswitch_studio.web_app.st.html") as render_html:
            _render_svg(svg, height=510)

        document = render_html.call_args.args[0]
        self.assertIn('width:100%;height:100%', document)
        self.assertIn('object-fit:contain', document)


if __name__ == "__main__":
    unittest.main()
