"""Native FFS controls, aggregate live progress and model-rate reporting."""
from __future__ import annotations

import json
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import Mock, patch

from aptaswitch_core.models import RunRequest, SequenceInput, SwitchArchitecture
from aptaswitch_studio import saved_switch_ffs_ui as ui
from aptaswitch_studio import saved_switch_library as library
from aptaswitch_studio.web_runtime import JobSnapshot
from tests.ui_helpers import FrenchAppTest as AppTest
from design_fixtures import design_result_fixture

SCRIPT = '''
import streamlit as st
from aptaswitch_studio.saved_switch_ffs_ui import render_saved_switch_ffs
render_saved_switch_ffs(st.session_state['saved'])
'''


class SavedSwitchKineticsUITests(unittest.TestCase):
    def setUp(self):
        from tests.native_ui_fixtures import patch_native_runtime
        patch_native_runtime(self)
        self.folder = Path(self.enterContext(tempfile.TemporaryDirectory())).resolve()
        request = RunRequest(sequences=SequenceInput('AACCGGTTAACC', 'ACGT' * 7 + 'AC'),
                             architecture=SwitchArchitecture(toehold_length=12, stem_length=18),
                             output_dir=self.folder, trials=1)
        result = design_result_fixture(request)
        snapshot = library.candidate_snapshot(result, result.candidates[0], 'Kinetic candidate')
        self.saved = library.save_switch(self.folder / library.LIBRARY_FOLDER, 'Kinetic candidate', snapshot)
        self.external = self.enterContext(patch.object(ui, 'configured_library_path', return_value='/native'))
        self.runtime = self.enterContext(patch.object(ui, 'detect_oxdna_runtime', return_value={'home': '/oxdna'}))
        self.active = self.enterContext(patch.object(ui, 'active_ffs', return_value=None))
        self.busy = self.enterContext(patch.object(ui, 'active_job_running', return_value=False))
        self.results = self.enterContext(patch.object(ui, 'saved_ffs_results', return_value=[]))

    def app(self):
        app = AppTest.from_string(SCRIPT)
        app.session_state['saved'] = self.saved
        app.run()
        self.assertFalse(app.exception)
        return app

    @staticmethod
    def button(app, label='Lancer la simulation cinétique'):
        return next(button for button in app.button if button.label == label)

    def result(self, *, mode='opening', status='model_estimate', value=2., interval=(1., 3.), unit='1/oxDNA_time'):
        run = self.folder / 'kinetic_run'
        run.mkdir(exist_ok=True)
        result = {'status': 'completed', 'run_id': 'kinetic_run', 'run_dir': str(run),
                  'method': 'NATIVE_FFS', 'kinetics': {'mode': mode, 'phase': 'report',
                  'replica': 4, 'replicas': 4, 'completed_replicas': 4,
                  'rate': {'status': status, 'value_model': value, 'interval95_model': interval,
                           'unit': unit, 'value_m_inv_model_time': 1234., 'experimental_rate': None}}}
        self.results.return_value = [result]
        return result

    def test_adjusted_counts_are_forwarded_to_one_complete_simulation(self):
        with patch.object(ui, 'start_ffs', return_value=(SimpleNamespace(saved=self.saved), True)) as start:
            app = self.app()
            controls = {widget.label: widget for widget in app.number_input}
            controls['Franchissements initiaux par répétition'].set_value(200)
            controls['Essais de production par interface'].set_value(1000)
            controls['Répliques indépendantes'].set_value(6)
            workers = min(2, ui.os.cpu_count() or 1)
            controls['Calculs simultanés'].set_value(workers).run()
            mode = next(widget for widget in app.radio if widget.label == 'Mode de simulation')
            self.assertEqual(mode.options, ['Cinétique complète', 'Départ proche du toehold'])
            self.assertEqual(mode.value, 'opening')
            self.button(app).click().run()
            start.assert_called_once()
            options = start.call_args.args[2]
            self.assertEqual(options['kinetics_mode'], 'opening')
            self.assertEqual([options[key] for key in ('kinetics_source_samples', 'kinetics_trials',
                                                      'kinetics_replicas', 'kinetics_workers')],
                             [200, 1000, 6, workers])
            self.assertNotIn('kinetics_pilot', options)
            self.assertNotIn('simulation_mode', options)
        self.assertFalse(app.exception)

    def test_start_disabled_for_unsupported_chemistry_or_length(self):
        self.saved.snapshot['candidate']['trigger_seq'] = 'ACGT'
        self.assertTrue(self.button(self.app()).disabled)
        self.saved.snapshot['candidate']['trigger_seq'] = 'ACGT' * 7 + 'AC'
        self.saved.snapshot['request']['thermo']['trigger_material'] = 'rna'
        self.assertTrue(self.button(self.app()).disabled)

    def test_rate_scope_units_interval_and_finite_box_concentration(self):
        self.result()
        app = self.app()
        self.assertEqual([(metric.label, metric.value) for metric in app.metric], [('Taux d’ouverture dans la boîte', '2 / τ')])
        self.assertTrue(any('Intervalle à 95 %' in item.value for item in app.caption))
        self.assertTrue(any('1234 M⁻¹ / temps oxDNA' in item.value for item in app.caption))
        self.assertIn('Unité : 1 / temps oxDNA', app.metric[0].proto.help)
        self.assertIn('pas d’une cinétique expérimentale', app.metric[0].proto.help)
        self.assertFalse(app.success)

    def test_missing_or_malformed_interval_is_not_invented(self):
        for interval in (None, [3., 4.], [-1., 3.], [float('nan'), 3.]):
            self.result(interval=interval)
            app = self.app()
            self.assertEqual(len(app.metric), 1)
            self.assertTrue(any('incertitude du taux n’est pas encore estimable' in item.value for item in app.caption))
            self.assertFalse(any('Intervalle à 95 %' in item.value for item in app.caption))

    def test_no_success_incomplete_and_invalid_rate_never_become_valid_metrics(self):
        for status, value, unit in [('no_success', 0., '1/oxDNA_time'), ('incomplete', 2., '1/oxDNA_time'),
                                   ('exploratory_modified_protocol', 2., '1/oxDNA_time'),
                                   ('model_estimate', 2., 'per_second'), ('model_estimate', float('nan'), '1/oxDNA_time'),
                                   ('model_estimate', True, '1/oxDNA_time')]:
            self.result(status=status, value=value, unit=unit)
            app = self.app()
            self.assertFalse(app.metric)
            if status == 'no_success':
                self.assertTrue(any('ne démontre pas que le taux est nul' in item.value for item in app.info))
            if status == 'incomplete':
                self.assertTrue(any('calcul est incomplet' in item.value for item in app.info))

    def test_flux_progress_shows_real_crossings_workers_and_completed_replicas(self):
        result = self.result(status='incomplete')
        result['kinetics'].update(phase='flux', workers=4, flux_crossings=59, flux_target=800,
                                  completed_replicas=0)
        app = self.app()
        captions = ' '.join(item.value for item in app.caption)
        self.assertIn('59/800 franchissements', captions)
        self.assertIn('4 simulations en parallèle', captions)
        self.assertIn('0/4 répétitions terminées', captions)
        crossing_control = next(item for item in app.number_input if item.label == 'Franchissements initiaux par répétition')
        self.assertIn('retour à l’état séparé', crossing_control.proto.help)
        self.assertIn('mesurer le flux initial', crossing_control.proto.help)
        self.assertFalse(any('pilote' in item.value.lower() or 'maintenus séparés' in item.value
                             for item in [*app.caption, *app.info]))

    def test_live_counts_refresh_per_replica_without_changing_active_run(self):
        run = self.folder / 'live'
        run.mkdir()
        request = {'kinetics': {'mode': 'opening', 'replicas': 4, 'trials': 1000, 'source_samples': 200},
                   'plan': {'interfaces': [{'name': 'contact_8', 'distance_max': 8., 'min_binding': 0},
                                           {'name': 'contact_7', 'distance_max': 7., 'min_binding': 0}]}}
        progress = {'method': 'NATIVE_FFS', 'status': 'running',
                    'kinetics': {**request['kinetics'], 'phase': 'production', 'workers': 4, 'completed_replicas': 0}}
        (run / 'progress.json').write_text(json.dumps(progress))
        folder = run / 'replica_001'
        folder.mkdir()
        state = {'index': 0, 'phase': 'production', 'production_index': 0, 'flux': {'crossings': 200, 'complete': True},
                 'stages': [{'name': 'contact_7', 'attempted': 12, 'successes': 2, 'failures': 9,
                             'pending': 1, 'target_trials': 1000, 'completed': False}]}
        (folder / 'progress.json').write_text(json.dumps(state))
        job = SimpleNamespace(saved=self.saved, folder=run, run_id='live', request=request,
            snapshot=lambda: JobSnapshot(state='running', current=0, total=4, message='FFS',
                                         elapsed_seconds=120, log=()), cancel=Mock())
        self.active.return_value = job
        self.busy.return_value = True
        app = self.app()
        self.assertEqual(len(app.get('vega_lite_chart')), 1)
        self.assertFalse(app.dataframe)
        table = app.table[0].value
        self.assertEqual(len(table), 4)
        self.assertIn('11', str(table.loc[0, 'Avancement étape']))
        self.assertIn('2', str(table.loc[0, 'Succès interface']))
        self.assertEqual(table.loc[1, 'État'], 'En attente')
        state['stages'][0].update(successes=3, pending=0)
        (folder / 'progress.json').write_text(json.dumps(state))
        app.run()
        self.assertFalse(app.exception)
        self.assertIn('3', str(app.table[0].value.loc[0, 'Succès interface']))
        self.assertIn('12', str(app.table[0].value.loc[0, 'Avancement étape']))
        job.cancel.assert_not_called()

        # A changed test protocol remains explicit across live rerenders.
        request['protocol_amendments'] = [{'reason': 'Test distance correction'}]
        request['plan']['distance_scope'] = 'trigger_to_toehold'
        app.run()
        self.assertFalse(app.exception)
        self.assertTrue(any('aucun taux cinétique' in message.value for message in app.info))
        self.assertIn('Distance toehold (nm)', app.table[0].value.columns)
        job.cancel.assert_not_called()

    def test_actual_native_aggregate_and_report_render_with_per_replica_details(self):
        from aptaswitch_core.oxdna_ffs import prepare_ffs_request
        from aptaswitch_core.oxdna_kinetics import direct_ffs_report
        from aptaswitch_core.oxdna_kinetics_worker import KineticSampler

        request = prepare_ffs_request(self.saved.snapshot, {}, {
            'kinetics_mode': 'opening', 'kinetics_replicas': 4,
            'kinetics_source_samples': 200, 'kinetics_trials': 1000,
        })
        output = self.folder / 'aggregate_controller'
        output.mkdir()
        worker = KineticSampler(request, output, started=1., clock=lambda: 16.)
        self.assertFalse(worker.proximity)
        states = {}
        for index, successes in enumerate((400, 450, 500, 550)):
            states[index] = {'replica_id': f'replica_{index+1}', 'index': index, 'completed': True,
                'phase': 'completed', 'steps': 100000, 'observation': {'binding': 1, 'distance': 1.},
                'flux': {'bank': list(range(200)), 'complete': True, 'exposure_steps': 100000},
                'stages': [{'name': interface['name'], 'successes': successes, 'failures': 1000-successes,
                            'attempted': 1000, 'target_trials': 1000, 'completed': True}
                           for interface in worker.interfaces[1:]]}
        result = worker.aggregate(states, 0, status='completed')
        self.assertEqual(result['simulation_mode'], 'opening')
        self.assertEqual(result['kinetics']['completed_replicas'], 4)
        self.assertEqual(result['kinetics']['flux_crossings'], 800)
        self.assertEqual(result['kinetics']['flux_target'], 800)
        self.assertEqual(result['elapsed_seconds'], 15.)
        result['kinetics'].update(phase='report', rate=direct_ffs_report(list(states.values()),
            dt=worker.request['sampling']['dt'], concentration_m=worker.request['concentration_each_m']))
        result.update(run_id='actual_native_report', run_dir=str(self.folder), repetitions=list(states.values()))
        self.results.return_value = [result]
        app = self.app()
        self.assertEqual([metric.label for metric in app.metric], ['Taux d’ouverture dans la boîte'])
        self.assertFalse(app.dataframe)
        table = next(exp for exp in app.expander if exp.label == 'Détails du calcul').table[0].value
        self.assertEqual(set(table['Répétition']), {1, 2, 3, 4})
        self.assertEqual(table.groupby('Répétition')['Succès'].first().tolist(), [400, 450, 500, 550])

    def test_native_controls_and_report_are_translated(self):
        from streamlit.testing.v1 import AppTest as RawAppTest

        self.result()
        app = RawAppTest.from_string('from aptaswitch_core.localization import set_language\nset_language("en")\n' + SCRIPT)
        app.session_state['saved'] = self.saved
        app.run()
        self.assertFalse(app.exception)
        self.assertIn('Simultaneous calculations', [widget.label for widget in app.number_input])
        self.assertIn('Initial crossings per replicate', [widget.label for widget in app.number_input])
        self.assertEqual([metric.label for metric in app.metric], ['Opening rate in the box'])
        self.assertTrue(any('Approximate 95% interval' in item.value for item in app.caption))


if __name__ == '__main__':
    unittest.main()
