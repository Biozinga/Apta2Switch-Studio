"""STOP screening precedes candidate physics and survives exports and UI."""
from __future__ import annotations
import csv
import json
import random
import tempfile
import unittest
from contextlib import ExitStack
from dataclasses import replace
from pathlib import Path
from unittest.mock import patch
import pandas as pd

from aptaswitch_core.exports import export_run
from aptaswitch_core.models import DesignCandidate, RunResult
from aptaswitch_core.scoring import sort_candidates
from aptaswitch_core.sequences import reading_frame_report
from aptaswitch_core.thermodynamics.design_constraints import UnsatisfiableDesignError, build_design_problem
from aptaswitch_core.thermodynamics.design_engine import design_candidates, evaluate_candidate
from tests.test_thermodynamics_design_engine import AnalyticalBackend, request as base_request
from tests.design_fixtures import design_result_fixture

ROOT = Path(__file__).resolve().parents[1]


def request(**changes):
    base = base_request(exclude_non_linear_candidates=False)
    return replace(base, exclude_stop_candidates=True, **changes)


def stopped_request(**changes):
    base = request()
    return replace(base, architecture=replace(base.architecture, frame_linker="UAA"), **changes)


def evaluate(options, trial=1):
    problem = build_design_problem(options)
    assignment = problem.generate(random.Random(17))
    backend = AnalyticalBackend(options)
    return evaluate_candidate(options, problem, assignment, backend, trial=trial), backend


def rendering_result(all_stops=False):
    """Synthetic UI record; target drawings are not thermodynamic predictions."""
    result = design_result_fixture(request())
    valid = result.candidates[0]
    valid.trial = 2
    start = build_design_problem(result.request).layout['aug_start'] + 3
    sequence = valid.switch_seq[:start] + 'UAA' + valid.switch_seq[start+3:]
    excluded = DesignCandidate.excluded_for_stop(trial=1, switch_seq=sequence,
        trigger_seq=valid.trigger_seq, first_stop=2, reporter_label=valid.reporter_label)
    result.candidates = [excluded] if all_stops else [valid, excluded]
    return result


class StopExclusionTests(unittest.TestCase):
    def setUp(self):
        from aptaswitch_core.localization import language_context
        self.enterContext(language_context('fr'))

    def test_default_rejects_unavoidable_stop_before_candidate_physics(self):
        options = stopped_request()
        self.assertTrue(options.exclude_stop_candidates)
        backend = AnalyticalBackend(options)
        with self.assertRaisesRegex(UnsatisfiableDesignError, "STOP"):
            design_candidates(options, backend)
        self.assertEqual(backend.calls, [])

    def test_stop_exclusion_record_keeps_sequence_without_fabricated_metrics(self):
        evaluated, _ = evaluate(stopped_request(exclude_stop_candidates=False))
        candidate = DesignCandidate.excluded_for_stop(trial=1, switch_seq=evaluated.switch_seq,
            trigger_seq=evaluated.trigger_seq, first_stop=evaluated.first_stop,
            reporter_label=evaluated.reporter_label)
        self.assertTrue(candidate.excluded)
        self.assertTrue(candidate.has_stop_codon)
        self.assertGreaterEqual(candidate.first_stop, 2)
        for field in ('score', 'defect', 'on_yield_pct', 'leak_pct', 'delta_g_off',
                      'delta_g_on', 'ddg_activation', 'delta_g_rbs_linker', 'has_bad_rbs_linker'):
            self.assertIsNone(getattr(candidate, field), field)
        self.assertEqual((candidate.structure_off, candidate.structure_on), ('', ''))

    def test_disabled_filter_scores_stops_and_preserves_optional_policy(self):
        options = stopped_request(exclude_stop_candidates=False)
        self.assertFalse(options.normalized().exclude_stop_candidates)
        self.assertFalse(options.to_dict()['exclude_stop_candidates'])
        candidate, backend = evaluate(options)
        self.assertTrue(candidate.has_stop_codon)
        self.assertFalse(candidate.excluded)
        self.assertIsNotNone(candidate.score)
        self.assertTrue(any(call[0] == 'analyze' for call in backend.calls))
        # Retained STOP candidates participate in ordinary score ranking.
        lower = replace(candidate, trial=2, has_stop_codon=False,
                        score=candidate.score-1, score_raw=candidate.score_raw-100)
        self.assertEqual(sort_candidates([lower, candidate])[0], candidate)

    def test_all_stops_export_without_false_best_or_fabricated_metrics(self):
        with tempfile.TemporaryDirectory() as tmp:
            options = stopped_request(output_dir=Path(tmp), trials=2)
            evaluated, _ = evaluate(replace(options, exclude_stop_candidates=False))
            candidates = [DesignCandidate.excluded_for_stop(trial=trial,
                switch_seq=evaluated.switch_seq, trigger_seq=evaluated.trigger_seq,
                first_stop=evaluated.first_stop, reporter_label=evaluated.reporter_label)
                for trial in (1, 2)]
            result = RunResult('stop_test', options, candidates, Path(tmp)/'run', 'start', 'end')
            paths = export_run(result)
            self.assertNotIn('svg', paths)
            saved = json.loads(Path(paths['json']).read_text())
            self.assertTrue(all(row['excluded'] and row['score'] is None for row in saved['results']))
            with Path(paths['csv']).open(newline='') as stream:
                rows = list(csv.DictReader(stream))
            self.assertTrue(all(row['Rank'] == row['Score_selection'] == '' and row['Status'] == 'excluded' for row in rows))
            self.assertIsNone(json.loads(Path(paths['project']).read_text())['top_svg'])
            import openpyxl
            book = openpyxl.load_workbook(paths['xlsx'])
            self.assertIsNone(book.active['C2'].value)
            self.assertEqual(book.active['Q2'].value, 'excluded')
            book.close()

    def test_excluded_candidates_remain_last_even_with_stale_high_scores(self):
        result = rendering_result()
        valid, excluded = result.candidates
        valid.score = -100
        excluded.score = 100
        self.assertEqual(sort_candidates([excluded, valid]), [valid, excluded])

    def test_frame_detector_accepts_dna_and_ignores_out_of_frame_and_reporter_stops(self):
        for sequence, excluded in (('ATGTAACCC', True), ('AUGCUACAA', False), ('AUGCAAUAA', False)):
            with self.subTest(sequence=sequence):
                self.assertEqual(reading_frame_report(sequence, 0, reporter_start=6)['has_premature_stop'], excluded)


class StopExclusionUITests(unittest.TestCase):
    def setUp(self):
        from tests.native_ui_fixtures import patch_native_runtime
        patch_native_runtime(self)
        self.stack = ExitStack()
        self.addCleanup(self.stack.close)
        self.output = self.stack.enter_context(tempfile.TemporaryDirectory())
        self.stack.enter_context(patch("aptaswitch_core.settings.load_settings", return_value={}))
        self.stack.enter_context(patch("aptaswitch_core.thermodynamics.setup.configured_library_path", return_value="/test-only/native"))

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
        app.radio(key="design_section").set_value("switch").run()
        self.assertTrue(app.checkbox(key="exclude_stop_candidates").value)
        app.checkbox(key="exclude_stop_candidates").uncheck().run()
        self.assertTrue(any("bien classés" in item.value for item in app.warning))
        app.radio(key="nav").set_value("Résultats").run()
        app.radio(key="nav").set_value("Conception").run()
        self.assertFalse(app.checkbox(key="exclude_stop_candidates").value)
        app.radio(key="design_section").set_value("system").run()
        app.text_area(key="aptamer_input").set_value("ACGTACGT")
        app.text_area(key="trigger_input").set_value("GTCCAGGCTGGTATAATTAGATCC").run()
        app.radio(key="design_section").set_value("switch").run()
        with patch("aptaswitch_studio.web_runtime.DesignRunJob.start"):
            next(b for b in app.button if b.label == "Lancer la conception").click().run()
        self.assertFalse(app.exception)
        self.assertFalse(app.session_state["design_job"].request.exclude_stop_candidates)

    def test_excluded_rows_are_last_without_rank_or_metrics_and_selection_skips_analysis(self):
        from tests.ui_helpers import candidate_tables
        result = rendering_result()
        app = self.app()
        app.session_state["last_result"] = result
        app.session_state["selected_trial"] = 1
        with patch("aptaswitch_studio.web_app.candidate_structure_svg") as structures:
            app.radio(key="nav").set_value("Résultats").run()
        self.assertFalse(app.exception)
        structures.assert_not_called()
        self.assertEqual(list(candidate_tables(app)[0].value["Essai"]), [2])
        next(button for button in app.button if button.label == "Agrandir le tableau complet").click().run()
        self.assertFalse(app.exception)
        rows = app.get("table")[-1].value
        self.assertEqual(rows.iloc[-1]["Statut"], "Exclu — STOP prématuré")
        self.assertTrue(pd.isna(rows.iloc[-1]["Score"]))
        self.assertTrue(pd.isna(rows.iloc[-1]["Rang"]))
        self.assertEqual(list(rows["Essai"]), [2, 1])
        self.assertTrue(any("Candidat exclu" in item.value for item in app.error))
        self.assertTrue(any(item.value == result.candidates[-1].switch_seq for item in app.code))

    def test_all_excluded_results_render_without_fabricated_structure(self):
        from tests.ui_helpers import candidate_tables
        result = rendering_result(all_stops=True)
        app = self.app()
        app.session_state["last_result"] = result
        app.radio(key="nav").set_value("Résultats").run()
        self.assertFalse(app.exception)
        self.assertTrue(any("Tous les candidats sont exclus" in item.value for item in app.warning))
        self.assertEqual(len(candidate_tables(app)), 0)
        next(button for button in app.button if button.label == "Agrandir le tableau complet").click().run()
        self.assertFalse(app.exception)
        self.assertEqual(len(app.get("table")[-1].value), 1)
        self.assertTrue(pd.isna(app.get("table")[-1].value.iloc[0]["Score"]))
        self.assertFalse(any(b.label == "Télécharger la structure PNG" for b in app.get("download_button")))


if __name__ == "__main__":
    unittest.main()
