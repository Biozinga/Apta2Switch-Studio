"""Screening never sends a trigger without specific, current permission."""

import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from tests.ui_helpers import FrenchAppTest as AppTest
from aptaswitch_core.blast_specificity import BlastScreenRequest, Taxon
from aptaswitch_studio.design_library import create_design, load_design
from aptaswitch_studio.web_runtime import JobSnapshot


TRIGGER = "ACGTACGTACGTACGTACGTACGTACGTAC"
TAXON = {"taxid": 4565, "name": "Triticum aestivum"}


class BlastUiTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.project = create_design(Path(self.temp.name), "Screening test")
        self.factory = self.enterContext(patch("aptaswitch_studio.blast_ui.BlastScreenJob"))
        self.job = self.factory.return_value
        self.job.snapshot.return_value = JobSnapshot(state="running", current=0, total=2, message="Pending")
        self.job.consume_once.return_value = False
        self.taxonomy = self.enterContext(patch("aptaswitch_studio.blast_ui.search_taxa", return_value=[Taxon(**TAXON)]))

    def app(self, **state):
        app = AppTest.from_string(
            "import streamlit as st\n"
            "from aptaswitch_studio.blast_ui import render_blast_screen\n"
            "render_blast_screen(st.session_state.get('trigger_input', ''))\n",
            default_timeout=20,
        )
        for key, value in {
            "active_design_path": str(self.project.path), "active_design_name": self.project.name,
            "trigger_input": TRIGGER, "blast_taxa": [TAXON], **state,
        }.items():
            app.session_state[key] = value
        app.run()
        self.assertFalse(app.exception)
        return app

    def test_render_and_taxonomy_lookup_never_submit_a_trigger(self):
        app = self.app(blast_taxa=[])
        self.assertTrue(app.button(key="blast_start").disabled)
        self.taxonomy.assert_not_called()
        app.text_input(key="blast_lookup_query").set_value("blé")
        app.button(key="blast_find_taxa").click().run()
        self.taxonomy.assert_called_once_with("blé")
        app.button(key="blast_add_taxon").click().run()
        self.assertEqual(app.session_state["blast_taxa"], [TAXON])
        self.assertFalse(app.checkbox(key="blast_consent").value)
        self.assertTrue(app.button(key="blast_start").disabled)
        self.factory.assert_not_called()

    def test_permission_consumed_for_one_specific_submission(self):
        app = self.app()
        self.assertTrue(app.button(key="blast_start").disabled)
        app.checkbox(key="blast_consent").check().run()
        self.assertFalse(app.button(key="blast_start").disabled)
        app.button(key="blast_start").click().run()
        self.assertFalse(app.exception)
        self.job.start.assert_called_once()
        request, directory = self.factory.call_args.args
        self.assertEqual(request.sequence, TRIGGER)
        self.assertEqual(request.taxa, (Taxon(**TAXON),))
        self.assertEqual(request.scopes, ("genomic", "transcript"))
        self.assertEqual(directory, self.project.path / "calculations")
        self.assertFalse(app.session_state["blast_consent"])
        self.assertNotIn("blast_consent", str(load_design(self.project.path).state))

    def test_editing_trigger_or_databases_requires_fresh_permission(self):
        app = self.app()
        app.checkbox(key="blast_consent").check().run()
        app.session_state["trigger_input"] = "T" + TRIGGER[1:]
        app.run()
        self.assertFalse(app.checkbox(key="blast_consent").value)
        self.assertTrue(app.button(key="blast_start").disabled)
        app.checkbox(key="blast_consent").check().run()
        app.multiselect(key="blast_scopes").set_value(["genomic"]).run()
        self.assertFalse(app.checkbox(key="blast_consent").value)
        self.assertTrue(app.button(key="blast_start").disabled)
        self.factory.assert_not_called()

    def test_changed_organism_does_not_reuse_permission(self):
        app = self.app()
        app.checkbox(key="blast_consent").check().run()
        app.button(key="blast_remove_4565").click().run()
        self.assertFalse(app.checkbox(key="blast_consent").value)
        self.assertTrue(app.button(key="blast_start").disabled)
        self.factory.assert_not_called()

    def test_full_app_navigation_keeps_selected_organisms_and_search_preferences(self):
        app = AppTest.from_file(str(Path(__file__).resolve().parents[1] / "src/aptaswitch_studio/web_app.py"), default_timeout=30)
        app.session_state["dependency_detection_done"] = True
        app.session_state["trigger_input"] = TRIGGER
        app.session_state["blast_taxa"] = [TAXON]
        app.run()
        self.assertFalse(app.exception)
        app.multiselect(key="blast_scopes").set_value(["transcript"]).run()
        app.number_input(key="blast_max_hits").set_value(25).run()
        app.radio(key="nav").set_value("Mes designs").run()
        app.radio(key="nav").set_value("Conception").run()
        self.assertFalse(app.exception)
        self.assertEqual(app.session_state["blast_taxa"], [TAXON])
        self.assertEqual(app.multiselect(key="blast_scopes").value, ["transcript"])
        self.assertEqual(app.number_input(key="blast_max_hits").value, 25)
        self.factory.assert_not_called()
        self.taxonomy.assert_not_called()

    def test_old_report_is_visible_but_marked_stale_without_resubmission(self):
        request = BlastScreenRequest(TRIGGER, (Taxon(**TAXON),)).normalized()
        result = {"request": request.to_dict(), "fingerprint": request.fingerprint(),
                  "status": "partial", "runs": [
                      {"taxon": TAXON, "scope": "genomic", "status": "completed", "hits": [], "database": "refseq_genomes"},
                      {"taxon": TAXON, "scope": "transcript", "status": "failed", "error": "Network unavailable", "hits": []},
                  ]}
        app = self.app(trigger_input="T" + TRIGGER[1:], blast_result=result)
        warnings = "\n".join(w.value for w in app.warning)
        self.assertIn("ancien trigger", warnings)
        self.assertIn("incomplet", warnings)
        self.assertIn("Network unavailable", warnings)
        self.assertTrue(any("ne valide pas l’absence" in item.value for item in app.info))
        self.factory.assert_not_called()
        self.taxonomy.assert_not_called()

    def test_complementary_predicted_transcript_is_not_labeled_available_in_sample(self):
        request = BlastScreenRequest(TRIGGER, (Taxon(**TAXON),)).normalized()
        hit = {"accession": "XM_123456.1", "title": "Predicted transcript", "hit_strand": "Minus",
               "identity_pct": 100.0, "query_coverage_pct": 100.0, "alignment_length": len(TRIGGER),
               "mismatches": 0, "gaps": 0, "evalue": 0.02, "qseq": TRIGGER, "hseq": TRIGGER}
        result = {"request": request.to_dict(), "fingerprint": request.fingerprint(), "status": "completed",
                  "runs": [{"taxon": TAXON, "scope": "transcript", "status": "completed", "hits": [hit],
                            "database": "refseq_rna", "scope_warning": "Scope not independently verified"}]}
        app = self.app(blast_result=result)
        captions = "\n".join(c.value for c in app.caption)
        self.assertIn("Transcrit prédit", captions)
        self.assertIn("pas automatiquement équivalent", captions)
        self.assertIn("présence dans votre échantillon n’est pas démontrée", captions)
        self.assertEqual(len(app.get("table")), 1)
        self.assertIn("Non renseigné par NCBI", app.get("table")[0].value.to_string())
        self.factory.assert_not_called()


if __name__ == "__main__":
    unittest.main()
