from __future__ import annotations

import copy
import json
import tempfile
import unittest
from unittest.mock import patch

from aptaswitch_core.structure_ensemble import EnsembleCancelled, load_ensemble, save_ensemble
from aptaswitch_studio.design_session import active_job_running
from aptaswitch_studio.ensemble_runtime import EnsembleJob
from tests.test_structure_ensemble import fixture as original_fixture, fake_analysis
from dataclasses import asdict, replace

def fixture(*args, backend=None, **kwargs):
    result = original_fixture(*args, **kwargs)
    if backend is not None:
        result.thermodynamics = {**asdict(backend.provenance), "engine": "native-hybrid",
                                 "conditions": {"celsius": result.request.thermo.temperature_c}}
    return result

from tests.ui_helpers import FrenchAppTest as AppTest


def ensemble_app():
    import streamlit as st
    from aptaswitch_core.localization import set_language
    from aptaswitch_studio.ensemble_ui import render_structure_ensemble
    set_language(st.session_state.get('ui_language', 'fr'))
    result = st.session_state.saved_result
    render_structure_ensemble(result, result.candidates[0])


class InlineJob(EnsembleJob):
    def start(self):
        self._run()


class RuntimeTests(unittest.TestCase):
    def setUp(self):
        from tests.native_ui_fixtures import patch_native_runtime
        self.backend = patch_native_runtime(self)

    def test_worker_saves_without_waiting_for_page_and_reports_write_failures(self):
        with tempfile.TemporaryDirectory() as directory:
            result = fixture(directory, backend=self.backend)
            data = fake_analysis(result, backend=self.backend)
            with patch('aptaswitch_studio.ensemble_runtime.analyze_structure_ensemble', return_value=data):
                job = EnsembleJob(result, result.candidates[0], '/fake')
                self.assertTrue(active_job_running({'ensemble_job': job}))
                job.start()
                job._thread.join(5)
                self.assertEqual(job.snapshot().state, 'done')
                self.assertEqual(load_ensemble(result.output_dir, data['inputs']), data)
                with patch('aptaswitch_studio.ensemble_runtime.save_ensemble', side_effect=OSError('disk full')):
                    job = InlineJob(result, result.candidates[0], '/fake')
                    job.start()
                    self.assertEqual(job.snapshot().state, 'failed')
                    self.assertIn('disk full', job.snapshot().error)

    def test_cancelled_worker_does_not_write_partial_results(self):
        result = fixture(backend=self.backend)
        with patch('aptaswitch_studio.ensemble_runtime.analyze_structure_ensemble', side_effect=EnsembleCancelled), \
             patch('aptaswitch_studio.ensemble_runtime.save_ensemble') as save:
            job = InlineJob(result, result.candidates[0], '/fake')
            job.cancel()
            job.start()
            self.assertEqual(job.snapshot().state, 'cancelled')
            self.assertFalse(active_job_running({'ensemble_job': job}))
            save.assert_not_called()


class EnsembleUITests(unittest.TestCase):
    def setUp(self):
        from tests.native_ui_fixtures import patch_native_runtime
        self.backend = patch_native_runtime(self)
        self.directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.result = fixture(self.directory.name, backend=self.backend)
        self.data = fake_analysis(self.result, backend=self.backend)

    def app(self):
        app = AppTest.from_function(ensemble_app, default_timeout=30)
        app.session_state.saved_result = self.result
        return app

    def test_loading_results_never_starts_native_calculation(self):
        with patch('aptaswitch_studio.ensemble_ui.EnsembleJob') as job:
            app = self.app().run()
            self.assertFalse(app.exception)
            self.assertEqual(app.button(key='ensemble_start').label, 'Analyser les sous-structures et familles')
            job.assert_not_called()

    def test_start_persist_reload_off_on_matrix_and_adjust_grouping_without_rescoring(self):
        saved = copy.deepcopy(self.result)
        with patch('aptaswitch_studio.ensemble_ui.EnsembleJob', InlineJob), \
             patch('aptaswitch_studio.ensemble_ui.configured_library_path', return_value='/fake'), \
             patch('aptaswitch_studio.ensemble_runtime.analyze_structure_ensemble', return_value=self.data) as analyze:
            app = self.app().run()
            app.button(key='ensemble_start').click().run()
            self.assertFalse(app.exception)
            self.assertEqual(analyze.call_count, 1)
            self.assertEqual(analyze.call_args.args[0].thermo.temperature_c, 31)
            self.assertEqual(app.table[0].value['Fréquence estimée (%)'].tolist(), [100.])
            state_control = next(r for r in app.radio if r.label == 'Ensemble à explorer')
            state_control.set_value('on').run()
            self.assertFalse(app.exception)
            self.assertIn('Préfixe RBS → fin AUG ON · écart (%)', app.table[0].value.columns)
            matrix = next(c for c in app.checkbox if c.label == 'Afficher la matrice de probabilités')
            matrix.check().run()
            self.assertFalse(app.exception)
            chart = app.get('vega_lite_chart')[0]
            spec = json.loads(chart.proto.spec)
            self.assertEqual(spec['mark']['type'], 'rect')
            self.assertEqual(spec['encoding']['color']['scale']['domain'], [0, 1])
            next(s for s in app.slider if s.label == 'Tolérance du regroupement (%)').set_value(0).run()
            self.assertFalse(app.exception)
            self.assertEqual(analyze.call_count, 1)
            self.assertEqual(app.session_state.saved_result, saved)
        # A completely new UI session reads the sidecar; it does not need the worker.
        with patch('aptaswitch_studio.ensemble_ui.EnsembleJob') as job:
            reopened = self.app().run()
            self.assertFalse(reopened.exception)
            self.assertEqual(len(reopened.table), 2)
            job.assert_not_called()

    def test_english_text_and_imported_analysis_validation(self):
        path = save_ensemble(self.result.output_dir, self.data)
        app = self.app()
        app.session_state.ui_language = 'en'
        app.run()
        self.assertFalse(app.exception)
        self.assertEqual(app.expander[0].label, 'OFF/ON substructures and families')
        self.assertEqual(app.radio[0].label, 'Ensemble to explore')
        path.write_text('{broken json')
        app.run()
        self.assertFalse(app.exception)
        self.assertIn('Unable to read saved analysis', app.warning[0].value)
        self.assertTrue(app.button(key='ensemble_start'))

    def test_reopen_finds_saved_larger_sample(self):
        data = fake_analysis(self.result, num_sample=3000, backend=self.backend)
        save_ensemble(self.result.output_dir, data)
        app = self.app().run()
        self.assertFalse(app.exception)
        self.assertEqual(app.selectbox[0].value, 3000)
        self.assertEqual(app.table[0].value['Tirages'].tolist(), [3000])

    def test_other_native_job_disables_analysis(self):
        job = EnsembleJob(self.result, self.result.candidates[0], '/fake')
        app = self.app()
        app.session_state.design_job = job
        app.run()
        self.assertFalse(app.exception)
        self.assertTrue(app.button(key='ensemble_start').disabled)
