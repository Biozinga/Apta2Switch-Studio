"""Proximity equilibration counts and safe snapshots remain visible in the UI."""
from copy import deepcopy
import json
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import Mock, patch

from aptaswitch_core.localization import language_context
from aptaswitch_studio import saved_switch_ffs_ui as ui
from aptaswitch_studio.ffs_progress_charts import replica_rows, selected_replica_chart
from aptaswitch_studio.web_runtime import JobSnapshot
from tests.ui_helpers import FrenchAppTest as AppTest


REQUEST = {'kinetics': {'mode': 'proximity', 'version': 4, 'replicas': 2,
                       'equilibration_steps': 200000,
                       'equilibration_method': 'isolated_monomers_before_close_placement_v1'},
           'plan': {}}


def monomers(switch, trigger):
    return {'target_steps_each': 200000, 'total_steps': 400000,
            'completed_steps': switch + trigger, 'switch_steps': switch,
            'trigger_steps': trigger, 'completed': switch == trigger == 200000}


COMPONENT = '''
import streamlit as st
from aptaswitch_studio import saved_switch_ffs_ui as ui
if st.session_state['action'] == 'monitor':
    ui._monitor(st.session_state['job'])
elif st.session_state['action'] == 'prepared':
    ui._render_prepared_view(st.session_state['folder'], st.session_state['request'],
                             st.session_state['data'], 'prepared_test')
elif st.session_state['action'] == 'result':
    ui._render_result(st.session_state['data'], 'result_test')
else:
    ui.render_saved_switch_ffs(st.session_state['saved'])
'''


class ProximityProgressTests(unittest.TestCase):
    def setUp(self):
        self.enterContext(language_context('fr'))

    def test_isolated_strand_and_total_counts_are_exact_without_ffs_trials(self):
        states = {0: {'phase': 'equilibration_switch', 'monomer_equilibration': monomers(45000, 0)},
                  1: {'phase': 'equilibration_trigger', 'monomer_equilibration': monomers(200000, 5000)}}
        before = deepcopy(states)
        rows = replica_rows(states, REQUEST)
        self.assertEqual(rows[0]['État'], 'Équilibrage du switch isolé')
        self.assertEqual(rows[0]['Avancement étape'], '45000/200000 pas')
        self.assertEqual(rows[0]['Équilibrage switch'], '45000/200000 pas')
        self.assertEqual(rows[0]['Équilibrage trigger'], '0/200000 pas')
        self.assertEqual(rows[0]['Équilibrage total'], '45000/400000 pas')
        self.assertEqual(rows[1]['État'], 'Équilibrage du trigger isolé')
        self.assertEqual(rows[1]['Avancement étape'], '5000/200000 pas')
        self.assertEqual(rows[1]['Équilibrage total'], '205000/400000 pas')
        self.assertEqual([row['Succès interface'] for row in rows], ['—', '—'])
        self.assertNotIn('essais', str(rows))
        self.assertNotIn('franchissement', str(rows))
        self.assertIsNone(selected_replica_chart(states[0], REQUEST))
        self.assertEqual(states, before)

    def test_missing_or_inconsistent_monomer_metadata_remains_missing(self):
        for value in ({}, 'invalid', {**monomers(10, 0), 'completed_steps': 11},
                      {**monomers(10, 0), 'switch_steps': 200001}):
            with self.subTest(value=value):
                rows = replica_rows({0: {'phase': 'equilibration_switch', 'monomer_equilibration': value}}, REQUEST)
                self.assertEqual(rows[0]['Avancement étape'], '—')
                self.assertEqual(rows[0]['Équilibrage total'], '—')
                self.assertEqual(rows[1]['Équilibrage trigger'], '—')

    def test_completed_equilibration_counts_survive_pause_and_translate(self):
        state = {'phase': 'initialization', 'status': 'paused', 'monomer_equilibration': monomers(200000, 200000)}
        with language_context('en'):
            row = replica_rows({0: state}, REQUEST, status='cancelled')[0]
            self.assertEqual(row['State'], 'Paused')
            self.assertEqual(row['Switch equilibration'], '200000/200000 steps')
            self.assertEqual(row['Total equilibration'], '400000/400000 steps')


class ProximitySnapshotTests(unittest.TestCase):
    def setUp(self):
        self.folder = Path(self.enterContext(tempfile.TemporaryDirectory())).resolve()
        (self.folder / 'preparation').mkdir()
        (self.folder / 'preparation/system.top').write_text('topology fixture')
        (self.folder / 'preparation/validated.dat').write_text('validated configuration fixture')
        (self.folder / 'preparation/preparation.json').write_text(json.dumps({
            'configuration_file': str(self.folder / 'preparation/validated.dat'),
            'topology_file': str(self.folder / 'preparation/system.top'), 'status': 'completed'}))
        self.movie = self.enterContext(patch.object(ui, 'render_oxdna_movie'))
        self.request = deepcopy(REQUEST)
        self.request['kinetics']['replicas'] = 1
        self.data = {'status': 'cancelled', 'run_dir': str(self.folder),
                     'requested_parameters': self.request, 'simulation_mode': 'proximity'}

    def app(self, action='result', **state):
        app = AppTest.from_string(COMPONENT, default_timeout=15)
        app.session_state['action'] = action
        app.session_state['folder'] = self.folder
        app.session_state['request'] = self.request
        app.session_state['data'] = self.data
        for key, value in state.items():
            app.session_state[key] = value
        app.run()
        self.assertFalse(app.exception)
        return app

    def replica(self):
        replica = self.folder / 'replica_001'
        replica.mkdir()
        (replica / 'progress.json').write_text(json.dumps({
            'index': 0, 'phase': 'equilibration_trigger', 'monomer_equilibration': monomers(200000, 5000)}))
        return replica

    def test_paused_replica_without_success_movie_shows_last_configuration(self):
        replica = self.replica()
        (replica / 'live.dat').write_text('last real MD configuration fixture')
        app = self.app()
        self.movie.assert_called_once()
        movie = self.movie.call_args.args[0]
        self.assertEqual(movie['movie_file'], 'replica_001/live.dat')
        self.assertEqual(movie['movie_kind'], 'prepared')
        self.assertIs(movie['final_interface_reached'], False)
        self.assertTrue(any('Dernière configuration 3D' in item.value for item in app.caption))
        self.assertFalse(any('Trajectoire FFS réussie assemblée' in item.value for item in app.caption))

    def test_successful_replica_movie_takes_priority_over_snapshot_and_preparation(self):
        replica = self.replica()
        (replica / 'live.dat').write_text('last configuration fixture')
        (replica / 'movie.dat').write_text('successful trajectory fixture')
        (replica / 'movie_metadata.json').write_text(json.dumps({'final_interface_reached': True}))
        self.app()
        self.movie.assert_called_once()
        self.assertEqual(self.movie.call_args.args[0]['movie_file'], 'replica_001/movie.dat')
        self.assertIs(self.movie.call_args.args[0]['final_interface_reached'], True)

    def test_prepared_configuration_remains_visible_when_replica_has_only_progress(self):
        self.replica()
        app = self.app()
        self.movie.assert_called_once()
        self.assertEqual(self.movie.call_args.args[0]['movie_file'], 'preparation/validated.dat')
        self.assertTrue(any('Configuration 3D préparée' in item.value for item in app.caption))

    def test_monitor_displays_validated_preparation_before_replica_progress_exists(self):
        saved = SimpleNamespace(path=self.folder / 'saved', snapshot={})
        job = SimpleNamespace(saved=saved, folder=self.folder, request=self.request, run_id='test_run', cancel=Mock(),
                              snapshot=lambda: JobSnapshot('running', 0, 1, 'Préparation 3D guidée', elapsed_seconds=4))
        with patch.object(ui, 'ffs_visualization', return_value=None):
            app = self.app('monitor', job=job)
        self.movie.assert_called_once()
        self.assertEqual(self.movie.call_args.args[0]['movie_file'], 'preparation/validated.dat')
        self.assertFalse(app.exception)
        job.cancel.assert_not_called()

    def test_global_trajectory_is_not_replaced_by_prepared_fallback(self):
        (self.folder / 'movie.dat').write_text('actual global trajectory fixture')
        self.data.update(movie_file='movie.dat', movie_kind='ffs')
        self.app()
        self.movie.assert_called_once()
        self.assertEqual(self.movie.call_args.args[0]['movie_file'], 'movie.dat')

    def test_fallback_rejects_outside_paths_and_symlinked_replica(self):
        outside = Path(self.enterContext(tempfile.TemporaryDirectory())).resolve()
        (outside / 'live.dat').write_text('outside fixture')
        (self.folder / 'replica_001').symlink_to(outside, target_is_directory=True)
        self.data['preparation'] = {'configuration_file': str(outside / 'live.dat')}
        (self.folder / 'preparation/preparation.json').unlink()
        for value in ('../outside.dat', str(outside / 'live.dat'), 'replica_001/live.dat'):
            self.assertIsNone(ui._contained_movie_file(self.folder, value))
        self.app('prepared')
        self.movie.assert_not_called()

    def test_new_proximity_form_describes_3d_isolated_equilibration_and_no_flux(self):
        saved = SimpleNamespace(path=self.folder / 'saved', snapshot={
            'request': {'thermo': {'switch_material': 'rna', 'trigger_material': 'dna',
                                   'temperature_c': 37., 'sodium_m': .12, 'magnesium_m': 0.}},
            'candidate': {'trigger_seq': 'ACGT' * 7 + 'AC'}})
        with patch.object(ui, 'active_ffs', return_value=None), \
             patch.object(ui, 'active_job_running', return_value=False), \
             patch.object(ui, 'saved_ffs_results', return_value=[]), \
             patch.object(ui, 'detect_oxdna_runtime', return_value={'home': 'fixture'}), \
             patch.object(ui, 'native_backend_identity', return_value='fixture'):
            app = self.app('form', saved=saved)
            app.radio(key=ui._key(saved) + '_kinetics_mode').set_value('proximity').run()
        self.assertFalse(app.exception)
        text = ' '.join(item.value for item in app.caption)
        self.assertIn('même préparation 3D guidée', text)
        self.assertIn('200 000 pas chacun', text)
        self.assertIn('Aucun flux initial n’est mesuré', text)
        self.assertFalse(any('Franchissements initiaux' in widget.label for widget in app.number_input))


if __name__ == '__main__':
    unittest.main()
