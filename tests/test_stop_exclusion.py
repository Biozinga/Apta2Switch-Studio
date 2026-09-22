"""STOP exclusion must precede thermodynamics, survive exports and remain optional."""

from __future__ import annotations

import csv
import json
import sys
import tempfile
import unittest
from contextlib import ExitStack, contextmanager
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import Mock, patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from aptaswitch_core.architecture import default_architecture_for_trigger
from aptaswitch_core.exports import export_run
from aptaswitch_core.models import RunRequest, ScoringWeights, SequenceInput
from aptaswitch_core.nupack_engine import run_nupack_engine
from aptaswitch_core.scoring import sort_candidates


def request(**kwargs):
    return RunRequest(
        sequences=SequenceInput(aptamer_dna="ACGTACGT", trigger_dna="GTCCAGGCTGGTATAATTAGATCC"),
        architecture=default_architecture_for_trigger(24), trials=kwargs.pop("trials", 2),
        nupack_path="/test-only/nupack", **kwargs,
    )


@contextmanager
def fake_design(sequences):
    """Real frame scanner, with only the external NUPACK boundary stubbed."""
    nupack = SimpleNamespace(
        Model=Mock(return_value=object()),
        mfe=Mock(side_effect=lambda *, strands, model: [SimpleNamespace(
            energy=-1.0, structure="+".join("." * (len(seq) - 1) for seq in strands),
        )]),
    )
    results = []
    for sequence in sequences:
        results.append([SimpleNamespace(
            to_analysis=lambda strand, seq=sequence: seq if strand == "switch" else "GTCCAGGCTGGTATAATTAGATCC",
            defects=SimpleNamespace(ensemble_defect=.01),
        )])
    spec = Mock()
    spec.run.side_effect = results
    design = {"spec": spec, "switch": "switch", "trigger": "trigger", "layout": {
        "aug_start": 0, "reporter_start_codon": 3,
        "rbs_linker_start": 0, "rbs_linker_end": 6,
    }}

    def tube_fraction(module, **kwargs):
        strands = kwargs["strands"]
        if len(strands) == 3:
            return 0.0
        return 100.0 if "UAA" in strands[1][1] else 20.0

    with (
        patch("aptaswitch_core.nupack_engine.nupack_import_context") as context,
        patch("aptaswitch_core.nupack_engine._build_design_spec", return_value=design),
        patch("aptaswitch_core.nupack_engine._complex_pct", side_effect=tube_fraction) as tubes,
    ):
        context.return_value.__enter__.return_value = nupack
        yield nupack, tubes, spec


class StopExclusionTests(unittest.TestCase):
    def setUp(self):
        from aptaswitch_core.localization import language_context
        self.enterContext(language_context("fr"))

    def test_default_skips_mfe_tubes_and_score_for_stop_but_keeps_candidate(self):
        progress = []
        with fake_design(["AUGUAACCC", "AUGCAACCC"]) as (nupack, tubes, spec):
            result = run_nupack_engine(request(), progress_callback=lambda *args: progress.append(args))
        self.assertEqual(spec.run.call_count, 2)  # The sequence must exist before it can be checked.
        self.assertEqual(nupack.mfe.call_count, 4)  # Four MFE calculations for the valid candidate only.
        self.assertEqual(tubes.call_count, 2)  # ON + leak for the valid candidate only.
        valid, excluded = result.candidates
        self.assertEqual((valid.trial, excluded.trial), (2, 1))
        self.assertTrue(excluded.excluded)
        self.assertEqual(excluded.first_stop, 2)
        self.assertEqual(excluded.switch_seq, "AUGUAACCC")
        for field in ("score", "defect", "on_yield_pct", "leak_pct", "delta_g_off", "delta_g_on", "ddg_activation", "delta_g_rbs_linker", "has_bad_rbs_linker"):
            self.assertIsNone(getattr(excluded, field), field)
        self.assertEqual(excluded.structure_off, "")
        self.assertEqual(excluded.structure_on, "")
        self.assertTrue(any(done == 1 and "exclu" in message for done, _, message in progress))
        self.assertEqual(progress[-1][0], 2)
        self.assertEqual(result.status, "completed")

    def test_disabled_filter_scores_stops_normally_and_they_can_win(self):
        options = request(exclude_stop_candidates=False,
                          scoring=ScoringWeights(on_yield=1, leak=0, defect=0, ddg=0, rbs=0))
        self.assertFalse(options.normalized().exclude_stop_candidates)
        self.assertFalse(options.to_dict()["exclude_stop_candidates"])
        with fake_design(["AUGUAACCC", "AUGCAACCC"]) as (nupack, tubes, _):
            result = run_nupack_engine(options)
        self.assertEqual(nupack.mfe.call_count, 8)
        self.assertEqual(tubes.call_count, 4)
        self.assertTrue(result.candidates[0].has_stop_codon)
        self.assertFalse(any(c.excluded for c in result.candidates))
        self.assertEqual(result.candidates[0].score, 100)
        self.assertEqual(result.candidates[1].score, 20)
        self.assertTrue(any("bien classés" in warning for warning in result.warnings))

    def test_all_stops_skip_scoring_and_export_without_a_false_best_candidate(self):
        with tempfile.TemporaryDirectory() as tmp:
            with fake_design(["AUGUAACCC", "AUGUGACCC"]) as (nupack, tubes, _), patch(
                "aptaswitch_core.nupack_engine.compute_selection_score", side_effect=AssertionError("Excluded candidates must not be scored")
            ):
                result = run_nupack_engine(request(output_dir=Path(tmp)))
            nupack.mfe.assert_not_called()
            tubes.assert_not_called()
            paths = export_run(result)
            self.assertEqual([c.trial for c in result.candidates], [1, 2])
            self.assertEqual(result.status, "completed")
            self.assertNotIn("svg", paths)
            saved = json.loads(Path(paths["json"]).read_text())
            self.assertTrue(saved["request"]["exclude_stop_candidates"])
            self.assertTrue(all(c["excluded"] and c["score"] is None for c in saved["results"]))
            with Path(paths["csv"]).open(newline="") as handle:
                rows = list(csv.DictReader(handle))
            self.assertTrue(all(row["Rank"] == "" and row["Score_selection"] == "" and row["Status"] == "excluded" for row in rows))
            self.assertIsNone(json.loads(Path(paths["project"]).read_text())["top_svg"])
            import openpyxl
            workbook = openpyxl.load_workbook(paths["xlsx"])
            self.assertIsNone(workbook.active["C2"].value)
            self.assertEqual(workbook.active["Q2"].value, "excluded")
            workbook.close()

    def test_exclusion_orders_after_even_negative_score_and_does_not_use_penalty(self):
        with fake_design(["AUGUAACCC", "AUGCAACCC"]):
            result = run_nupack_engine(request())
        valid, excluded = result.candidates
        valid.score = -100
        excluded.score = 100  # Even stale metrics cannot override exclusion status.
        self.assertEqual(sort_candidates([excluded, valid]), [valid, excluded])

    def test_frame_detector_accepts_dna_stops_and_ignores_out_of_frame_or_reporter_stops(self):
        for sequence, excluded in (("ATGTAACCC", True), ("AUGCUACAA", False), ("AUGCAAUAA", False)):
            with self.subTest(sequence=sequence), fake_design([sequence]):
                result = run_nupack_engine(request(trials=1))
            self.assertEqual(result.candidates[0].excluded, excluded)



class StopExclusionUITests(unittest.TestCase):
    def setUp(self):
        self.stack = ExitStack()
        self.addCleanup(self.stack.close)
        self.output = self.stack.enter_context(tempfile.TemporaryDirectory())
        self.stack.enter_context(patch("aptaswitch_core.nupack_setup.load_settings", return_value={}))
        self.stack.enter_context(patch("aptaswitch_core.nupack_setup.configured_nupack_import_path", return_value="/test-only/nupack"))

    def app(self):
        from tests.ui_helpers import FrenchAppTest as AppTest
        app = AppTest.from_file(str(ROOT / "src/aptaswitch_studio/web_app.py"), default_timeout=20)
        app.session_state["dependency_detection_done"] = True
        app.session_state["output_dir"] = self.output
        app.run()
        self.assertFalse(app.exception)
        return app

    def test_default_checkbox_opt_out_warning_and_request_persist_across_navigation(self):
        app = self.app()
        self.assertTrue(app.checkbox(key="exclude_stop_candidates").value)
        app.checkbox(key="exclude_stop_candidates").uncheck().run()
        self.assertTrue(any("bien classés" in item.value for item in app.warning))
        app.radio(key="nav").set_value("Résultats").run()
        app.radio(key="nav").set_value("Conception").run()
        self.assertFalse(app.checkbox(key="exclude_stop_candidates").value)
        app.text_area(key="aptamer_input").set_value("ACGTACGT")
        app.text_area(key="trigger_input").set_value("GTCCAGGCTGGTATAATTAGATCC").run()
        with patch("aptaswitch_studio.web_runtime.DesignRunJob.start"):
            next(b for b in app.button if b.label == "Lancer la conception").click().run()
        self.assertFalse(app.exception)
        self.assertFalse(app.session_state["design_job"].request.exclude_stop_candidates)

    def test_excluded_rows_are_last_without_rank_or_metrics_and_selection_skips_analysis(self):
        with fake_design(["AUGUAACCC", "AUGCAACCC"]):
            result = run_nupack_engine(request())
        app = self.app()
        app.session_state["last_result"] = result
        app.session_state["selected_trial"] = 1
        with patch("aptaswitch_studio.web_app.candidate_structure_svg") as structures:
            app.radio(key="nav").set_value("Résultats").run()
        self.assertFalse(app.exception)
        structures.assert_not_called()
        rows = app.get("table")[0].value
        self.assertEqual(rows.iloc[-1]["Statut"], "Exclu — STOP prématuré")
        self.assertEqual(rows.iloc[-1]["Score"], "—")
        self.assertEqual(rows.iloc[-1]["Rang"], "—")
        self.assertEqual(list(rows["Essai"]), [2, 1])
        self.assertTrue(any("Candidat exclu" in item.value for item in app.error))
        self.assertTrue(any(item.value == "AUGUAACCC" for item in app.code))

    def test_all_excluded_results_render_without_fabricated_structure(self):
        with fake_design(["AUGUAACCC"]):
            result = run_nupack_engine(request(trials=1))
        app = self.app()
        app.session_state["last_result"] = result
        app.radio(key="nav").set_value("Résultats").run()
        self.assertFalse(app.exception)
        self.assertTrue(any("Tous les candidats sont exclus" in item.value for item in app.warning))
        self.assertEqual(len(app.get("table")), 1)
        self.assertFalse(any(b.label == "Télécharger la structure PNG" for b in app.get("download_button")))


if __name__ == "__main__":
    unittest.main()
