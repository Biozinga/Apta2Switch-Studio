"""Unexpected full-context complexes can be inspected without changing conditions."""
from copy import deepcopy
from pathlib import Path
from types import SimpleNamespace
import tempfile
import unittest
from unittest.mock import patch

from aptaswitch_core.models import ReporterContext, RunRequest, SequenceInput, SwitchArchitecture
from aptaswitch_core.saved_switch_analysis import MAX_COMPLEX_SIZE, build_saved_switch_context, validate_analysis_config
from aptaswitch_studio.saved_switch_library import candidate_snapshot
from tests.design_fixtures import design_result_fixture
from tests.ui_helpers import FrenchAppTest as AppTest
from tests.native_ui_fixtures import patch_native_runtime


SCRIPT = """
import streamlit as st
from aptaswitch_core.localization import set_language
from aptaswitch_studio.saved_switch_ui import _render_analysis
set_language(st.session_state.get('ui_language', 'fr'))
_render_analysis(st.session_state.analysis, st.session_state.saved_candidate)
"""


class SavedComplexUITests(unittest.TestCase):
    def setUp(self):
        self.backend = patch_native_runtime(self)
        request = RunRequest(sequences=SequenceInput('AACCGGTT', 'AACC'),
            architecture=SwitchArchitecture(toehold_length=2, stem_length=2), trials=1,
            reporter=ReporterContext(key='test', label='Test reporter',
                sequence_after_start_rna='ACG'*10, full_sequence_after_start_rna='ACG'*20))
        result = design_result_fixture(request)
        snapshot = candidate_snapshot(result, result.candidates[0], 'Candidate')
        folder = Path(self.enterContext(tempfile.TemporaryDirectory()))
        self.saved = SimpleNamespace(snapshot=snapshot, path=folder, name='Candidate')
        config = validate_analysis_config(snapshot, {'temperature_c': 29., 'include_aptamer': True,
                                                     'temperature_values_c': [25., 30.]})
        self.context = build_saved_switch_context(snapshot, config)
        point = {'conditions': {k: config[k] for k in ('temperature_c', 'sodium_m', 'magnesium_m')},
                 'metrics': {'on_yield_pct': 65., 'on_trigger_bound_pct': 60., 'off_switch_monomer_pct': 75.}}
        self.extra_mfe = self.mfe(['Switch', 'Switch'])
        for state in ('off', 'on'):
            names = ['Switch'] if state == 'off' else ['Trigger', 'Switch']
            populations = [self.row(names, self.mfe(names)), self.row(['Switch', 'Switch'], self.extra_mfe),
                           self.row(['Aptamer', 'Aptamer'], self.mfe(['Aptamer', 'Aptamer']))]
            point[state] = {'mfe': self.mfe(names), 'population': {'complexes': populations}}
        self.analysis = {'baseline': point, 'config': config, 'context': self.context,
                         'run_id': 'selected_analysis', 'status': 'completed', 'sweeps': {}}
        self.enterContext(patch('aptaswitch_studio.saved_switch_ui.active_job_running', return_value=False))

    def mfe(self, names):
        sequences = [self.context[{'Switch':'switch_sequence', 'Trigger':'trigger_sequence', 'Aptamer':'aptamer_sequence'}[n]] for n in names]
        return {'strand_order': names, 'sequences': sequences, 'energy_kcal_mol': -4.2,
                'structure': '+'.join('.'*len(s) for s in sequences)}

    def row(self, names, mfe):
        return {'label': ' + '.join(names), 'complex_key': '|'.join(names), 'strand_order': list(names),
                'concentration_m': 1e-7,
                'switch_count': names.count('Switch'), 'trigger_count': names.count('Trigger'),
                'aptamer_count': names.count('Aptamer'), 'mfe': mfe}

    def app(self):
        app = AppTest.from_string(SCRIPT, default_timeout=30)
        app.session_state.analysis = self.analysis
        app.session_state.saved_candidate = self.saved
        return app.run()

    def test_selects_homodimers_and_aptamers_with_their_exact_full_structure(self):
        from aptaswitch_core.saved_switch_diagrams import context_complex_diagram
        before = deepcopy(self.analysis)
        with patch('aptaswitch_studio.saved_switch_ui.context_complex_diagram', wraps=context_complex_diagram) as render, \
             patch('aptaswitch_studio.saved_switch_ui.start_analysis') as start:
            app = self.app()
            self.assertFalse(app.exception)
            self.assertFalse(app.error)
            self.assertTrue(any('jusqu’à 3 brins' in c.value for c in app.caption))
            selectors = [s for s in app.selectbox if s.label == 'Complexe à visualiser']
            self.assertEqual(len(selectors), 2)
            selectors[0].set_value('Switch|Switch').run()
            self.assertFalse(app.exception)
            self.assertFalse(app.error)
            self.assertEqual(render.call_args_list[-2].args[2], self.extra_mfe)
            self.assertEqual(len(self.extra_mfe['sequences'][0]), self.context['switch_length'])
            self.assertEqual([t.label for t in app.tabs], ['Structure 2D', 'Séquence linéaire']*2)
            next(s for s in app.selectbox if s.key.endswith('_on_complex')).set_value('Aptamer|Aptamer').run()
            self.assertFalse(app.exception)
            self.assertEqual(render.call_args.args[2]['strand_order'], ['Aptamer', 'Aptamer'])
            start.assert_not_called()
            self.assertEqual(self.analysis, before)

    def test_missing_structure_explicitly_recalculates_saved_baseline_without_sweeps(self):
        del self.analysis['baseline']['off']['population']['complexes'][1]['mfe']
        self.analysis['config']['max_complex_size'] = 2
        self.analysis['config']['temperature_c'] = 41.  # the recorded baseline, not form state, is authoritative
        with patch('aptaswitch_studio.saved_switch_ui.start_analysis', return_value=(SimpleNamespace(run_id='new_context'), True)) as start, \
             patch('aptaswitch_studio.saved_switch_ui.configured_library_path', return_value='/native'):
            app = self.app()
            next(s for s in app.selectbox if s.key.endswith('_off_complex')).set_value('Switch|Switch').run()
            self.assertFalse(app.exception)
            self.assertTrue(any('pas encore été calculée' in i.value for i in app.info))
            start.assert_not_called()
            next(b for b in app.button if b.label == 'Calculer les structures de ces complexes').click().run()
            self.assertFalse(app.exception)
            start.assert_called_once()
            self.assertEqual(start.call_args.args[1]['temperature_c'], 29.)
            self.assertEqual(start.call_args.args[1]['max_complex_size'], 3)
            self.assertEqual(self.analysis['config']['max_complex_size'], 2)
            self.assertTrue(start.call_args.args[1]['include_aptamer'])
            for axis in ('temperature_values_c', 'sodium_values_m', 'magnesium_values_m'):
                self.assertEqual(start.call_args.args[1][axis], [])
            self.assertEqual(self.analysis['config']['temperature_values_c'], [25., 30.])
            self.assertEqual(app.session_state['saved_switch_analysis_pending_' + self.saved.snapshot['snapshot_id']], 'new_context')

    def test_pentamers_keep_repeated_copies_and_distinct_cyclic_orders_selectable(self):
        from aptaswitch_core.saved_switch_diagrams import context_complex_diagram
        self.analysis['config']['max_complex_size'] = 5  # historical result remains inspectable
        orders = [['Trigger', 'Switch', 'Switch', 'Aptamer', 'Aptamer'],
                  ['Trigger', 'Switch', 'Aptamer', 'Switch', 'Aptamer']]
        added = []
        for names in orders:
            row = self.row(names, self.mfe(names))
            row['label'] = 'Trigger + 2 Switch + 2 Aptamer [' + ' → '.join(names) + ']'
            added.append(row)
        self.analysis['baseline']['on']['population']['complexes'].extend(added)
        repeated = self.row(['Switch'] * 5, self.mfe(['Switch'] * 5))
        self.analysis['baseline']['off']['population']['complexes'].append(repeated)
        before = deepcopy(self.analysis)
        with patch('aptaswitch_studio.saved_switch_ui.context_complex_diagram', wraps=context_complex_diagram) as render, \
             patch('aptaswitch_studio.saved_switch_ui.start_analysis') as start:
            app = self.app()
            self.assertFalse(app.exception)
            self.assertTrue(any('dans ce résultat : jusqu’à 5 brins' in item.value for item in app.caption))
            on = next(widget for widget in app.selectbox if widget.key.endswith('_on_complex'))
            self.assertIn(added[0]['label'], on.options)
            self.assertIn(added[1]['label'], on.options)
            for row in added:
                next(widget for widget in app.selectbox if widget.key.endswith('_on_complex')).set_value(row['complex_key']).run()
                self.assertFalse(app.exception)
                self.assertFalse(app.error)
                self.assertEqual(render.call_args.args[2]['strand_order'], row['strand_order'])
                self.assertEqual(len(render.call_args.args[2]['sequences']), 5)
                self.assertEqual(render.call_args.args[2]['sequences'].count(self.context['switch_sequence']), 2)
            next(widget for widget in app.selectbox if widget.key.endswith('_off_complex')).set_value(repeated['complex_key']).run()
            self.assertFalse(app.exception)
            self.assertEqual(render.call_args_list[-2].args[2]['strand_order'], ['Switch'] * 5)
            self.assertTrue(any('Traits bleus' in item.value and 'Arcs gris' in item.value and 'structure MFE' in item.value
                                for item in app.caption))
            start.assert_not_called()
        self.assertEqual(self.analysis, before)

    def test_historical_result_caption_keeps_recorded_two_strand_limit(self):
        self.analysis['config']['max_complex_size'] = 2
        app = self.app()
        self.assertFalse(app.exception)
        self.assertTrue(any('dans ce résultat : jusqu’à 2 brins' in item.value for item in app.caption))
        provenance = next(item for item in app.expander if item.label == 'Conditions, portée du calcul et provenance')
        self.assertFalse(provenance.proto.expanded)
        for prefix in ('Conditions de ce résultat : 29 °C', 'Complexes explorés dans ce résultat : jusqu’à 2 brins'):
            self.assertEqual(sum(item.value.startswith(prefix) for item in app.caption), 1)
            self.assertTrue(any(item.value.startswith(prefix) for item in provenance.caption))
        self.assertIn('Télécharger le contexte FASTA', [item.proto.label for item in provenance.get('download_button')])
        self.assertFalse(any('jusqu’à 5 brins' in item.value for item in app.caption))
        self.assertEqual(self.analysis['config']['max_complex_size'], 2)
        del self.analysis['config']['max_complex_size']
        app.run()
        self.assertFalse(app.exception)
        self.assertTrue(any('taille maximale des complexes n’est pas enregistrée' in item.value for item in app.caption))
        self.assertFalse(any('jusqu’à 5 brins' in item.value for item in app.caption))

    def test_new_context_uses_three_even_when_saved_settings_and_draft_used_five(self):
        self.assertEqual(MAX_COMPLEX_SIZE, 3)
        old = {**self.analysis['config'], 'max_complex_size': 5}
        script = """
import streamlit as st
from aptaswitch_studio.saved_switch_ui import _context_controls
config, context = _context_controls(st.session_state.saved_candidate)
st.session_state.shown_config = config
"""
        with patch('aptaswitch_studio.saved_switch_ui.load_context', return_value={'config': old}), \
             patch('aptaswitch_studio.saved_switch_ui.start_analysis', return_value=(SimpleNamespace(run_id='max3'), True)) as start, \
             patch('aptaswitch_studio.saved_switch_ui.configured_library_path', return_value='/native'):
            app = AppTest.from_string(script)
            app.session_state.saved_candidate = self.saved
            app.session_state.saved_switch_context_drafts = {self.saved.snapshot['snapshot_id']: old}
            app.run()
            self.assertFalse(app.exception)
            self.assertEqual(app.session_state.shown_config['max_complex_size'], 3)
            launch = next(button for button in app.button if button.label == 'Lancer l’analyse thermodynamique')
            self.assertIn('Nouveau calcul : complexes de un à 3 brins', launch.proto.help)
            next(button for button in app.button if button.label == 'Lancer l’analyse thermodynamique').click().run()
            self.assertFalse(app.exception)
            self.assertEqual(start.call_args.args[1]['max_complex_size'], 3)
            self.assertEqual(old['max_complex_size'], 5)

    def test_english_scope_and_mfe_pairing_caption_are_translated(self):
        app = self.app()
        app.session_state['ui_language'] = 'en'
        app.run()
        self.assertFalse(app.exception)
        captions = ' '.join(item.value for item in app.caption)
        self.assertIn('up to 3 strands', captions)
        self.assertIn('Blue lines', captions)
        self.assertIn('Gray arcs', captions)
        self.assertIn('MFE structure', captions)


if __name__ == '__main__':
    unittest.main()
