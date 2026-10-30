"""The free-MD flow exposes only saved extensions and descriptive comparisons."""
from __future__ import annotations
import json
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import Mock, patch
from aptaswitch_core.models import RunRequest, SequenceInput, SwitchArchitecture
from aptaswitch_studio import saved_switch_free_md_ui as ui
from aptaswitch_studio import saved_switch_library as library
from aptaswitch_studio.web_runtime import JobSnapshot
from tests.ui_helpers import FrenchAppTest as AppTest
from design_fixtures import design_result_fixture

SCRIPT = '''
import streamlit as st
from aptaswitch_studio.saved_switch_free_md_ui import render_saved_switch_free_md
render_saved_switch_free_md(st.session_state.saved)
'''


class FreeMDUITests(unittest.TestCase):
    def setUp(self):
        from tests.native_ui_fixtures import patch_native_runtime
        self.native_backend_fixture = patch_native_runtime(self)
        self.root = Path(self.enterContext(tempfile.TemporaryDirectory())).resolve()
        request = RunRequest(sequences=SequenceInput('AACCGGTTAACC', 'AACC'),
            architecture=SwitchArchitecture(toehold_length=2, stem_length=2), trials=1, output_dir=self.root)
        self.result = design_result_fixture(request)
        self.saved = self.save(True)
        self.enterContext(patch.object(ui, 'detect_oxdna_runtime', return_value={'home': '/oxdna'}))
        self.enterContext(patch.object(ui, 'configured_library_path', return_value='/integrated-backend'))
        self.native_backend = self.enterContext(patch.object(ui, 'native_backend_identity', return_value='integrated-test-backend'))
        self.active = self.enterContext(patch.object(ui, 'active_free_md', return_value=None))
        self.busy = self.enterContext(patch.object(ui, 'active_job_running', return_value=False))
        self.results = self.enterContext(patch.object(ui, 'saved_free_md_results', return_value=[]))

    def save(self, extended):
        original = {'aptamer_dna': 'AACCGGTTAACC', 'trigger_dna': 'CC'} if extended else None
        name = 'Extended' if extended else 'No extension'
        snap = library.candidate_snapshot(self.result, self.result.candidates[0], name, original_sequences=original)
        return library.save_switch(library.saved_switch_root(self.root), name, snap)

    def app(self, script=SCRIPT):
        app = AppTest.from_string(script, default_timeout=15)
        app.session_state.saved = self.saved
        app.run()
        self.assertFalse(app.exception)
        return app

    def button(self, app, label='Lancer la MD libre ADN–ADN'):
        return next(widget for widget in app.button if widget.label == label)

    def test_absent_extension_hides_feature_without_manual_sequence_entry(self):
        self.saved = self.save(False)
        app = self.app()
        self.assertFalse(app.button)
        self.assertFalse(app.text_input)
        self.assertFalse(app.number_input)

    def test_saved_candidate_section_exists_only_for_an_actual_extension(self):
        plain = self.save(False)
        script = '''
import streamlit as st
from aptaswitch_studio.saved_switch_ui import render_saved_switches_page
render_saved_switches_page()
'''
        with patch('aptaswitch_studio.saved_switch_ui._render_candidate'):
            app = AppTest.from_string(script)
            app.session_state.active_design_path = str(self.root)
            app.session_state.saved_switch_selected = str(self.saved.path)
            app.run()
            self.assertFalse(app.exception)
            self.assertIn('MD libre ADN–ADN', app.radio(key='saved_switch_section_selector').options)
            app.session_state.saved_switch_selected = str(plain.path)
            app.run()
            self.assertFalse(app.exception)
            self.assertNotIn('MD libre ADN–ADN', app.radio(key='saved_switch_section_selector').options)

    def test_monitored_region_even_allocation_and_time_help_forwarded(self):
        job = SimpleNamespace(saved=self.saved)
        with patch.object(ui, 'start_free_md', return_value=(job, True)) as start:
            app = self.app()
            self.assertTrue(self.button(app).disabled)
            app.text_input[0].set_value('2-4, 9').run()
            controls = {widget.label: widget for widget in app.number_input}
            self.assertIn('propre échelle de temps', controls['Durée approximative par réplique (µs)'].proto.help)
            controls['Nombre total de répliques (pair)'].set_value(6).run()
            self.assertTrue(any('3 répliques avec le trigger original et 3' in item.value for item in app.caption))
            self.button(app).click().run()
        self.assertFalse(app.exception)
        options = start.call_args.args[2]
        self.assertEqual(options['monitored_bases'], [2, 3, 4, 9])
        self.assertEqual(options['replicas'], 6)
        self.assertEqual(options['duration_us'], .1)
        self.assertNotIn('trigger_sequence', options)
        self.assertEqual(options['thermodynamics_library_path'], '/integrated-backend')
        self.assertNotIn('import_path', options)
        self.assertIs(app.session_state.saved_switch_free_md_job, job)

    def test_native_backend_required_and_physical_controls_match_integrated_model(self):
        self.native_backend.return_value = None
        app = self.app()
        app.text_input[0].set_value('2-4').run()
        self.assertTrue(self.button(app).disabled)
        self.assertTrue(any('moteur thermodynamique intégré est indisponible' in item.value for item in app.info))
        self.native_backend.return_value = 'integrated-test-backend'
        app.run()
        controls = {widget.label: widget for widget in app.number_input}
        temperature = controls['Température de simulation (°C)'].proto
        salt = controls['Sel monovalent Na⁺ (M)'].proto
        self.assertEqual((temperature.min, temperature.max), (25., 55.))
        self.assertEqual((salt.min, salt.max), (.12, 1.))
        self.assertFalse(self.button(app).disabled)

    def test_odd_count_and_out_of_range_base_cannot_launch(self):
        app = self.app()
        app.text_input[0].set_value('500').run()
        self.assertTrue(self.button(app).disabled)
        self.assertTrue(app.error)
        app.text_input[0].set_value('2').run()
        next(w for w in app.number_input if w.label == 'Nombre total de répliques (pair)').set_value(5).run()
        self.assertTrue(self.button(app).disabled)
        self.assertTrue(any('nombre pair' in item.value for item in app.error))

    def test_pause_keeps_controls_disabled_until_worker_stops_and_resume_uses_history(self):
        folder = self.root / 'run'
        folder.mkdir()
        job = SimpleNamespace(saved=self.saved, folder=folder, run_id='run', cancel=Mock(),
                              snapshot=lambda: JobSnapshot(state='running', current=1, total=4, message='Sampling', elapsed_seconds=1, log=()))
        self.active.return_value = job
        self.busy.return_value = True
        app = self.app()
        self.assertTrue(app.text_input[0].disabled)
        self.button(app, 'Mettre en pause').click().run()
        job.cancel.assert_called_once()
        self.active.return_value = None
        self.busy.return_value = False
        self.results.return_value = [{'run_id': 'run', 'run_dir': str(folder), 'status': 'cancelled', 'resumable': True}]
        with patch.object(ui, 'resume_free_md', return_value=(job, True)) as resume:
            app.run()
            self.button(app, 'Reprendre la simulation').click().run()
        resume.assert_called_once_with(self.saved, folder)
        self.assertFalse(app.exception)

    def test_partial_curves_show_group_counts_and_raw_exports(self):
        folder = self.root / 'partial'
        folder.mkdir()
        comparison = {'status': 'partial', 'expected_per_group': 3, 'groups': {
            variant: {'completed': 1, 'available': 2, 'summary': {},
                      'site': [{'base': 2, 'contact_fraction': .2, 'rmsf_nm': .3}]}
            for variant in ('original', 'extended')},
            'curves': [{'variant': variant, 'step': step, 'time_us': step / 100., 'metric': 'site_contact_fraction',
                        'mean': .3, 'sd': .1, 'n': 2} for variant in ('original', 'extended') for step in (0, 1)]}
        self.results.return_value = [{'run_id': 'partial', 'run_dir': str(folder), 'status': 'cancelled',
                                     'resumable': True, 'comparison': comparison, 'replicas': []}]
        app = self.app()
        self.assertEqual(len(app.get('vega_lite_chart')), 3)
        self.assertTrue(any('Comparaison partielle' in item.value for item in app.info))
        self.assertTrue(any('pas un intervalle de confiance' in item.value for item in app.caption))
        with patch.object(ui, 'make_free_md_archive', return_value=b'zip') as archive:
            self.button(app, 'Préparer l’archive des données brutes').click().run()
        archive.assert_called_once_with(folder)
        self.assertFalse(app.exception)

    def test_monitor_form_survives_navigation(self):
        script = '''
import streamlit as st
from aptaswitch_studio.saved_switch_free_md_ui import render_saved_switch_free_md
if st.radio('Page', ['MD', 'Autre']) == 'MD':
    render_saved_switch_free_md(st.session_state.saved)
'''
        app = self.app(script)
        app.text_input[0].set_value('3-5').run()
        app.radio[0].set_value('Autre').run()
        app.radio[0].set_value('MD').run()
        self.assertEqual(app.text_input[0].value, '3-5')
        self.assertFalse(self.button(app).disabled)

    def test_resumed_run_invalidates_previously_prepared_raw_archive(self):
        folder = self.root / 'resumed_archive'
        folder.mkdir()
        result = {'run_id': 'same_run', 'run_dir': str(folder), 'status': 'cancelled',
                  'completed_at': '2026-10-14T01:00:00+00:00', 'elapsed_seconds': 10.,
                  'replicas': [{'replica_index': 0, 'variant': 'original', 'sample_count': 2,
                                'production_steps': 20, 'status': 'cancelled'}]}
        self.results.return_value = [result]
        app = self.app()
        cache_key = ui._key(self.saved) + '_same_run_archive'
        with patch.object(ui, 'make_free_md_archive', return_value=b'old data'):
            self.button(app, 'Préparer l’archive des données brutes').click().run()
        self.assertEqual(app.session_state[cache_key], b'old data')
        # A further pause after resuming retains the same run ID and status.
        result.update(completed_at='2026-10-14T02:00:00+00:00', elapsed_seconds=20.)
        result['replicas'][0].update(sample_count=4, production_steps=40)
        app.run()
        self.assertNotIn(cache_key, app.session_state)
        with patch.object(ui, 'make_free_md_archive', return_value=b'new data'):
            self.button(app, 'Préparer l’archive des données brutes').click().run()
        self.assertEqual(app.session_state[cache_key], b'new data')
        self.assertFalse(app.exception)

    def test_selected_replica_uses_its_own_movie_topology_and_free_md_title(self):
        folder = self.root / 'movies'
        folder.mkdir()
        replicas = [{'replica_index': index, 'variant': variant,
                     'trajectory_file': f'replica_{index:03d}/trajectory.dat',
                     'topology_file': f'replica_{index:03d}/system.top',
                     'visualization': {'domains': [variant]}}
                    for index, variant in enumerate(('original', 'extended'))]
        self.results.return_value = [{'run_id': 'movies', 'run_dir': str(folder), 'status': 'completed',
                                     'replicas': replicas}]
        with patch.object(ui, 'render_oxdna_movie') as movie:
            app = self.app()
            self.assertEqual(movie.call_args.args[0]['movie_file'], 'replica_000/trajectory.dat')
            next(w for w in app.selectbox if w.label == 'Réplique à visualiser').set_value(1).run()
            shown = movie.call_args.args[0]
            self.assertEqual(shown['movie_file'], 'replica_001/trajectory.dat')
            self.assertEqual(shown['topology_file'], 'replica_001/system.top')
            self.assertEqual(shown['visualization']['domains'], ['extended'])
            self.assertEqual(shown['trajectory_endpoint'], 'free_md')
        self.assertFalse(app.exception)

    def test_live_native_replica_states_and_production_counter_are_shown(self):
        folder = self.root / 'live'
        folder.mkdir()
        (folder / 'progress.json').write_text(json.dumps({'fraction': .25, 'production_steps': 20,
            'target_steps': 80, 'replicas': 2, 'replica_states': [
                {'replica_index': 0, 'variant': 'original', 'phase': 'production', 'status': 'running',
                 'production_steps': 20, 'samples_count': 3},
                {'replica_index': None, 'variant': None, 'status': None}]}))
        job = SimpleNamespace(saved=self.saved, folder=folder, run_id='live', cancel=Mock(),
            snapshot=lambda: JobSnapshot(state='running', current=20, total=80, message='Sampling', elapsed_seconds=1, log=()))
        self.active.return_value = job
        self.busy.return_value = True
        app = self.app()
        self.assertTrue(any('20 / 80 pas MD' in item.value for item in app.caption))
        self.assertEqual(len(app.dataframe), 1)
        self.assertEqual(app.dataframe[0].value['Observations'].tolist(), [3, 0])


    def _live_app(self, *, fraction=.25, steps=100, total=400, replicas=4):
        folder = self.root / 'live_dashboard'
        folder.mkdir()
        request = {'replicas': replicas, 'sampling': {'production_steps': total // replicas, 'dt': .003},
                   'monitored_bases': [2], 'variants': [
                       {'name': name, 'aptamer_sequence': 'AACCGGTTAACC', 'trigger_sequence': seq,
                        'topology_file': f'preparation/{name}/system.top'}
                       for name, seq in (('original', 'CC'), ('extended', 'AACC'))]}
        (folder / 'worker_request.json').write_text(json.dumps(request))
        progress = {'production_steps': steps, 'target_steps': total, 'fraction': fraction,
                    'completed_replicas': 0, 'replicas': replicas, 'replica_states': []}
        (folder / 'progress.json').write_text(json.dumps(progress))
        # Detached jobs historically count completed replicas in this snapshot.
        # The progress bar must instead use the production counters on disk.
        self.active.return_value = SimpleNamespace(saved=self.saved, folder=folder, run_id='live_dashboard',
            request=request, cancel=Mock(), snapshot=lambda: JobSnapshot(state='running', current=0,
                total=replicas, message='oxDNA DNA2 · 0/4 réplicats terminés', elapsed_seconds=1, log=()))
        self.busy.return_value = True
        return folder, request

    def test_live_percent_uses_steps_instead_of_completed_replicas_or_stale_fraction(self):
        self._live_app(fraction=.9)
        app = self.app()
        progress = app.get('progress')[0].proto
        self.assertEqual(progress.value, 25)
        self.assertIn('25.00', progress.text)
        self.assertIn('0/4', progress.text)
        self.assertEqual(len(app.dataframe[0].value), 4)
        selector = next(w for w in app.selectbox if w.label == 'Réplique à visualiser')
        self.assertEqual(len(selector.options), 4)
        selector.set_value(3).run()
        self.assertFalse(app.exception)
        self.assertEqual(selector.value, 3)
        self.assertTrue(any('première sauvegarde' in c.value for c in app.caption))

    def test_live_selected_checkpoint_refreshes_and_uses_correct_group_topology(self):
        folder, request = self._live_app(replicas=2, total=200)
        for index, variant in enumerate(request['variants']):
            replica = folder / f'replica_{index:03d}'
            replica.mkdir()
            topology = folder / variant['topology_file']
            topology.parent.mkdir(parents=True)
            topology.write_text('16 2 5->3\nAACCGGTTAACC type=DNA\nAACC type=DNA\n')
            (replica / 'checkpoint_0.dat').write_text('test checkpoint')
            state = {'replica_index': index, 'variant': variant['name'], 'generation': 2,
                     'status': 'running', 'phase': 'production', 'production_steps': 25 + 10 * index,
                     'samples_count': 3, 'checkpoint_file': 'checkpoint_0.dat',
                     'topology_file': variant['topology_file']}
            (replica / 'state.json').write_text(json.dumps(state))
        with patch.object(ui, 'render_oxdna_movie') as movie:
            app = self.app()
            shown = movie.call_args.args[0]
            self.assertEqual(shown['movie_file'], 'replica_000/checkpoint_0.dat')
            self.assertEqual(shown['movie_kind'], 'live')
            self.assertEqual(shown['trajectory_endpoint'], 'free_md')
            self.assertEqual(shown['replica_number'], 1)
            next(w for w in app.selectbox if w.label == 'Réplique à visualiser').set_value(1).run()
            shown = movie.call_args.args[0]
            self.assertEqual(shown['movie_file'], 'replica_001/checkpoint_0.dat')
            self.assertEqual(shown['topology_file'], 'preparation/extended/system.top')
            self.assertEqual(shown['replica_number'], 2)
            self.assertEqual(app.get('progress')[1].proto.value, 35)
            self.assertIsInstance(shown['visualization'], dict)
            state_path = folder / 'replica_001/state.json'
            state = json.loads(state_path.read_text())
            state.update(generation=3, checkpoint_file='checkpoint_1.dat', production_steps=60)
            (folder / 'replica_001/checkpoint_1.dat').write_text('new checkpoint')
            state_path.write_text(json.dumps(state))
            app.run()
            self.assertEqual(movie.call_args.args[0]['movie_file'], 'replica_001/checkpoint_1.dat')
            self.assertEqual(app.get('progress')[1].proto.value, 60)
            self.assertEqual(next(w for w in app.selectbox if w.label == 'Réplique à visualiser').value, 1)
        self.assertFalse(app.exception)

    def test_live_preparation_and_small_percent_do_not_show_fake_completion(self):
        folder, _ = self._live_app(steps=27, total=10000)
        replica = folder / 'replica_000'
        replica.mkdir()
        (replica / 'state.json').write_text(json.dumps({'replica_index': 0, 'variant': 'original',
            'phase': 'equilibration', 'phase_steps': 2000, 'production_steps': 0, 'status': 'running'}))
        app = self.app()
        self.assertIn('0.27', app.get('progress')[0].proto.text)
        self.assertEqual(app.get('progress')[1].proto.value, 0)
        self.assertFalse(app.exception)

    def test_live_rejects_checkpoint_outside_run(self):
        folder, _ = self._live_app()
        replica = folder / 'replica_000'
        replica.mkdir()
        outside = self.root / 'outside.dat'
        outside.write_text('private')
        (replica / 'state.json').write_text(json.dumps({'replica_index': 0, 'variant': 'original',
            'phase': 'production', 'production_steps': 1, 'status': 'running',
            'checkpoint_file': str(outside), 'topology_file': 'preparation/original/system.top'}))
        with patch.object(ui, 'render_oxdna_movie') as movie:
            app = self.app()
            movie.assert_not_called()
        self.assertFalse(app.exception)
