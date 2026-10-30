"""Full-reporter family exploration uses the selected persisted context."""
from pathlib import Path
from dataclasses import asdict
from types import SimpleNamespace
import tempfile
import unittest
from unittest.mock import patch

from aptaswitch_core.models import ReporterContext, RunRequest, SequenceInput, SwitchArchitecture
from aptaswitch_core.saved_switch_analysis import build_saved_switch_context, validate_analysis_config
from aptaswitch_core.structure_ensemble import load_ensemble
from aptaswitch_studio import saved_switch_library as library
from aptaswitch_studio import saved_switch_runtime as runtime
from aptaswitch_studio.ensemble_runtime import EnsembleJob
from aptaswitch_studio.saved_switch_ensemble_runtime import load_context_ensemble_source
from tests.design_fixtures import design_result_fixture
from tests.test_structure_ensemble import fake_analysis
from tests.ui_helpers import FrenchAppTest as AppTest
from tests.native_ui_fixtures import patch_native_runtime


SCRIPT = """
import streamlit as st
from aptaswitch_core.localization import set_language
from aptaswitch_studio.saved_switch_ui import _render_thermodynamics
set_language(st.session_state.get('ui_language', 'fr'))
_render_thermodynamics(st.session_state.saved_candidate)
"""


class FullContextEnsembleUITests(unittest.TestCase):
    def setUp(self):
        self.backend = patch_native_runtime(self)
        self.enterContext(patch("aptaswitch_studio.saved_switch_ensemble_runtime.get_backend", return_value=self.backend))
        self.folder = Path(self.enterContext(tempfile.TemporaryDirectory())).resolve()
        request = RunRequest(sequences=SequenceInput('AACCGGTT', 'AACC'),
            architecture=SwitchArchitecture(toehold_length=2, stem_length=2), trials=1,
            reporter=ReporterContext(key='test', label='Test reporter',
                sequence_after_start_rna='ACG' * 10, full_sequence_after_start_rna='ACG' * 20))
        design = design_result_fixture(request)
        snapshot = library.candidate_snapshot(design, design.candidates[0], 'Full candidate')
        self.saved = library.save_switch(library.saved_switch_root(self.folder), snapshot['name'], snapshot)
        self.analysis = self.save_analysis(29.)
        self.source, self.candidate = load_context_ensemble_source(self.saved, self.analysis)
        self.ensemble = fake_analysis(SimpleNamespace(request=self.source.request,
            candidates=[self.candidate], output_dir=self.source.output_dir), backend=self.backend)
        self.enterContext(patch('aptaswitch_studio.saved_switch_ui.configured_library_path', return_value='/fake'))
        self.enterContext(patch('aptaswitch_studio.ensemble_ui.configured_library_path', return_value='/fake'))

    def save_analysis(self, temperature):
        config = validate_analysis_config(self.saved.snapshot, {'temperature_c': temperature}, backend=self.backend)
        context = build_saved_switch_context(self.saved.snapshot, config, backend=self.backend)
        point = {'conditions': {field: config[field] for field in ('temperature_c', 'sodium_m', 'magnesium_m')},
                 'thermodynamics': asdict(self.backend.provenance), 'metrics': {'on_yield_pct': 65., 'on_trigger_bound_pct': 60., 'off_switch_monomer_pct': 75.}}
        for state in ('off', 'on'):
            seqs = ([context['trigger_sequence']] if state == 'on' else []) + [context['switch_sequence']]
            point[state] = {'population': {'complexes': []},
                'mfe': {'structure': '+'.join('.' * len(s) for s in seqs), 'sequences': seqs,
                        'strand_order': ['Trigger', 'Switch'] if state == 'on' else ['Switch'], 'energy_kcal_mol': -1.}}
        job = runtime.SavedSwitchAnalysisJob(self.saved, config, '/fake')
        with patch.object(runtime, 'analyze_saved_switch', return_value={'status': 'completed', 'baseline': point,
                  'sweeps': {}, 'config': config, 'context': context,
                  'model_limits': self.backend.condition_limits, 'thermodynamics': asdict(self.backend.provenance)}):
            job._run()
        return job.snapshot().result

    def app(self):
        app = AppTest.from_string(SCRIPT, default_timeout=30)
        app.session_state.saved_candidate = self.saved
        app.session_state.active_design_path = str(self.folder)
        return app.run()

    def family_button(self, app):
        return next(b for b in app.button if b.label == 'Analyser les sous-structures et familles')

    def test_full_context_families_sample_saved_conditions_and_reopen_without_recalculation(self):
        from aptaswitch_core.preview import candidate_structure_svg
        with patch('aptaswitch_studio.ensemble_runtime.analyze_structure_ensemble', return_value=self.ensemble) as compute, \
             patch.object(EnsembleJob, 'start', lambda job: job._run()), \
             patch('aptaswitch_studio.ensemble_ui.candidate_structure_svg', wraps=candidate_structure_svg) as diagram:
            app = self.app()
            self.assertFalse(app.exception)
            self.assertFalse(app.error)
            compute.assert_not_called()
            next(n for n in app.number_input if n.label == 'Température (°C)').set_value(42.).run()
            self.family_button(app).click().run()
            self.assertFalse(app.exception)
            self.assertFalse(app.error)
            self.assertEqual(compute.call_count, 1)
            request, candidate, _ = compute.call_args.args
            self.assertEqual(request.thermo.temperature_c, 29.)
            self.assertEqual(len(request.reporter.sequence_after_start_rna), 60)
            self.assertEqual(candidate.switch_seq, self.analysis['context']['switch_sequence'])
            self.assertGreater(len(candidate.switch_seq), len(self.saved.snapshot['candidate']['switch_seq']))
            self.assertEqual(diagram.call_args.args[2], candidate.switch_seq)
            self.assertEqual(len(diagram.call_args.kwargs['reporter_sequence']), 60)
            table = next(t.value for t in app.get("table") if 'Fréquence estimée (%)' in t.value.columns)
            self.assertEqual(table['Fréquence estimée (%)'].tolist(), [100.])
            next(r for r in app.radio if r.label == 'Ensemble à explorer').set_value('on').run()
            self.assertFalse(app.exception)
            self.assertEqual(diagram.call_args.kwargs['state'], 'on')
            self.assertIn('+', diagram.call_args.args[3])
            next(s for s in app.slider if s.label == 'Tolérance du regroupement (%)').set_value(0).run()
            self.assertEqual(compute.call_count, 1)
            stored = load_ensemble(self.source.output_dir, self.ensemble['inputs'])
            self.assertEqual({key: stored[key] for key in self.ensemble}, self.ensemble)
            self.assertEqual(stored['context_source']['run_id'], self.analysis['run_id'])
            fresh = self.app()
            self.assertFalse(fresh.exception)
            self.assertTrue(any(r.label == 'Ensemble à explorer' for r in fresh.radio))
            self.assertEqual(compute.call_count, 1)

    def test_different_selected_context_has_its_own_family_analysis(self):
        from aptaswitch_core.structure_ensemble import save_ensemble
        save_ensemble(self.source.output_dir, self.ensemble)
        other = self.save_analysis(35.)
        with patch('aptaswitch_studio.ensemble_runtime.analyze_structure_ensemble') as compute:
            app = self.app()
            self.assertFalse(app.exception)
            self.assertTrue(self.family_button(app))
            self.assertFalse(any(r.label == 'Ensemble à explorer' for r in app.radio))
            selector = next(s for s in app.selectbox if s.label == 'Analyse thermodynamique enregistrée')
            selector.set_value(1).run()
            self.assertFalse(app.exception)
            self.assertTrue(any(r.label == 'Ensemble à explorer' for r in app.radio))
            self.assertNotEqual(other['run_id'], self.analysis['run_id'])
            compute.assert_not_called()


if __name__ == '__main__':
    unittest.main()
