"""Single kinetic simulation controls retain active jobs, settings and saved results."""
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


class SavedSwitchFFSUITests(unittest.TestCase):
    def setUp(self):
        from tests.native_ui_fixtures import patch_native_runtime
        patch_native_runtime(self)
        self.folder = Path(self.enterContext(tempfile.TemporaryDirectory())).resolve()
        request = RunRequest(sequences=SequenceInput('AACCGGTTAACC', 'ACGT' * 7 + 'AC'),
                             architecture=SwitchArchitecture(toehold_length=12, stem_length=18),
                             output_dir=self.folder, trials=1)
        result = design_result_fixture(request)
        snapshot = library.candidate_snapshot(result, result.candidates[0], 'Retained')
        self.saved = library.save_switch(self.folder / library.LIBRARY_FOLDER, 'Retained', snapshot)
        self.external = self.enterContext(patch.object(ui, 'configured_library_path', return_value='/native'))
        self.detect = self.enterContext(patch.object(ui, 'detect_oxdna_runtime', return_value={'home': '/local/oxdna'}))
        self.active = self.enterContext(patch.object(ui, 'active_ffs', return_value=None))
        self.busy = self.enterContext(patch.object(ui, 'active_job_running', return_value=False))
        self.results = self.enterContext(patch.object(ui, 'saved_ffs_results', return_value=[]))

    def app(self, script=SCRIPT):
        app = AppTest.from_string(script)
        app.session_state['saved'] = self.saved
        app.run()
        self.assertFalse(app.exception)
        return app

    @staticmethod
    def button(app, label='Lancer la simulation cinétique'):
        return next(button for button in app.button if button.label == label)

    def interrupted_result(self, resumable=True):
        run = self.folder / 'interrupted'
        run.mkdir(exist_ok=True)
        self.results.return_value = [{'status': 'interrupted', 'method': 'NATIVE_FFS', 'run_id': 'interrupted',
            'run_dir': str(run), 'resumable': resumable, 'kinetics': {'status': 'sampling', 'phase': 'production'}}]
        return run

    def running_job(self):
        folder = self.folder / 'running'
        folder.mkdir()
        job = SimpleNamespace(saved=self.saved, folder=folder, run_id='running', request={},
            snapshot=lambda: JobSnapshot(state='running', current=1, total=4, message='FFS',
                                         elapsed_seconds=30, log=('Sampling',)), cancel=Mock())
        self.active.return_value = job
        self.busy.return_value = True
        return job

    def test_opening_default_exposes_flux_and_parallel_controls(self):
        app = self.app()
        self.assertEqual([button.label for button in app.button], ['Lancer la simulation cinétique'])
        self.assertEqual([title.value for title in app.subheader], ['Cinétique · rapprochement, hybridation et ouverture'])
        self.assertEqual(app.radio[0].value, "opening")
        self.assertEqual(len(app.number_input), 6)
        self.assertEqual([expander.label for expander in app.expander], ['Paramètres avancés', 'Modèle et critères'])
        self.assertEqual(len(app.expander[0].number_input), 4)
        self.assertFalse(app.expander[0].proto.expanded)
        self.assertFalse(app.expander[1].proto.expanded)
        self.assertFalse(app.warning)
        self.assertFalse(app.info)
        labels = [widget.label for widget in app.number_input]
        self.assertIn('Franchissements initiaux par répétition', labels)
        self.assertIn('Calculs simultanés', labels)
        self.assertFalse(any('pilote' in label.lower() for label in labels))
        self.assertFalse(any('Budget' in label or 'Mg' in label for label in labels))
        controls = {widget.label: widget for widget in app.number_input}
        for label, default, minimum, maximum in (
            ('Franchissements initiaux par répétition', 75, 75, 200),
            ('Essais de production par interface', 250, 250, 1000),
        ):
            widget = controls[label]
            self.assertEqual(widget.value, default)
            self.assertEqual(widget.min, minimum)
            self.assertEqual(widget.max, maximum)

    def test_opening_start_sends_flux_and_worker_settings_with_conditions(self):
        job = SimpleNamespace(saved=self.saved)
        with patch.object(ui, 'start_ffs', return_value=(job, True)) as start:
            app = self.app()
            controls = {widget.label: widget for widget in app.number_input}
            controls['Température de simulation (°C)'].set_value(42.)
            app.selectbox(key=ui._key(self.saved) + '_salt_unit').set_value('M').run()
            app.number_input(key=ui._key(self.saved) + '_salt_value').set_value(.25).run()
            self.button(app).click().run()
        self.assertFalse(app.exception)
        self.assertEqual(start.call_args.args, (self.saved, {'home': '/local/oxdna'}, {
            'kinetics_mode': 'opening', 'kinetics_trials': 250, 'kinetics_replicas': 10,
            'kinetics_source_samples': 75, 'kinetics_workers': min(10, ui.os.cpu_count() or 1), 'temperature_c': 42., 'salt_m': .25, 'thermodynamics_library_path': '/native'}))
        self.assertIs(app.session_state['saved_switch_ffs_job'], job)

    def test_salt_units_preserve_molar_options_without_starting_a_run(self):
        key = ui._key(self.saved)
        with patch.object(ui, 'start_ffs', return_value=(SimpleNamespace(saved=self.saved), True)) as start:
            app = self.app()
            app.number_input(key=key + '_salt_value').set_value(175.).run()
            for unit, value in (('M', .175), ('µM', 175000.), ('nM', 175000000.), ('mM', 175.)):
                app.selectbox(key=key + '_salt_unit').set_value(unit).run()
                self.assertFalse(app.exception)
                self.assertEqual(app.number_input(key=key + '_salt_value').value, value)
                self.assertEqual(app.session_state[key + '_draft']['salt_m'], .175)
                start.assert_not_called()
            self.button(app).click().run()
            self.assertFalse(app.exception)
            self.assertEqual(start.call_args.args[2]['salt_m'], .175)

    def test_salt_value_unit_and_physical_draft_survive_navigation(self):
        script = SCRIPT.replace("render_saved_switch_ffs(st.session_state['saved'])", """
page = st.radio('Page', ['Simulation', 'Réglages'], key='page')
if page == 'Simulation':
    render_saved_switch_ffs(st.session_state['saved'])
""")
        key = ui._key(self.saved)
        app = self.app(script)
        app.selectbox(key=key + '_salt_unit').set_value('M').run()
        app.number_input(key=key + '_salt_value').set_value(.225).run()
        app.radio(key='page').set_value('Réglages').run()
        self.assertFalse(app.number_input)
        app.radio(key='page').set_value('Simulation').run()
        self.assertFalse(app.exception)
        self.assertEqual(app.selectbox(key=key + '_salt_unit').value, 'M')
        self.assertEqual(app.number_input(key=key + '_salt_value').value, .225)
        self.assertEqual(app.session_state[key + '_draft']['salt_m'], .225)

    def test_active_run_disables_both_salt_value_and_unit_without_changing_it(self):
        key = ui._key(self.saved)
        app = self.app()
        app.selectbox(key=key + '_salt_unit').set_value('M').run()
        app.number_input(key=key + '_salt_value').set_value(.175).run()
        job = self.running_job()
        with patch.object(ui, 'start_ffs') as start:
            app.run()
            self.assertFalse(app.exception)
            self.assertTrue(app.number_input(key=key + '_salt_value').disabled)
            self.assertTrue(app.selectbox(key=key + '_salt_unit').disabled)
            self.assertEqual(app.session_state[key + '_draft']['salt_m'], .175)
            start.assert_not_called()
            job.cancel.assert_not_called()

    def test_explains_complete_route_from_random_separated_strands(self):
        app = self.app()
        criteria = next(item for item in app.expander if item.label == 'Modèle et critères')
        self.assertFalse(criteria.proto.expanded)
        self.assertTrue(any('positions et orientations aléatoires' in item.value for item in criteria.markdown))
        self.assertTrue(any('la dernière base du trigger côté 3′' in item.value
                            and 'la première base du toehold après le leader' in item.value
                            and 'ne contraint pas leur orientation' in item.value for item in criteria.markdown))
        self.assertIn('Choisissez le départ', app.subheader[0].proto.help)
        self.assertFalse(any('réutilisé' in item.value for item in app.caption))
        self.assertEqual(app.radio[0].value, "opening")

    def test_running_monitor_requires_pause_before_new_launch(self):
        job = self.running_job()
        (job.folder / 'progress.json').write_text(json.dumps({'kinetics': {'mode': 'opening', 'phase': 'flux'}}))
        with patch.object(ui, 'start_ffs') as start:
            app = self.app()
            self.assertTrue(self.button(app).disabled)
            self.assertFalse(self.button(app, 'Mettre en pause').disabled)
            self.assertTrue(any('Simulation en cours' in item.value for item in app.caption))
            self.assertFalse(any('pause avant' in item.value for item in app.info))
            start.assert_not_called()
            job.cancel.assert_not_called()

    def test_isolated_live_refresh_restores_session_language_without_touching_job(self):
        from streamlit.testing.v1 import AppTest as RawAppTest
        job = self.running_job()
        script = '''
import streamlit as st
from aptaswitch_core.localization import set_language
from aptaswitch_studio.saved_switch_ffs_ui import _monitor
set_language('en' if st.session_state['ui_language'] == 'fr' else 'fr')
_monitor(st.session_state['job'])
'''
        for language, label in (('fr', 'Mettre en pause'), ('en', 'Pause')):
            with self.subTest(language=language):
                app = RawAppTest.from_string(script)
                app.session_state['ui_language'] = language
                app.session_state['job'] = job
                app.run()
                self.assertFalse(app.exception)
                self.assertEqual([button.label for button in app.button], [label])
        job.cancel.assert_not_called()

    def test_default_ten_replicas_workers_fit_available_cpu(self):
        for cores, workers in ((4, 4), (16, 10)):
            with self.subTest(cores=cores), patch.object(ui.os, 'cpu_count', return_value=cores):
                app = self.app()
                values = {widget.label: widget.value for widget in app.number_input}
                self.assertEqual(values['Répliques indépendantes'], 10)
                self.assertEqual(values['Calculs simultanés'], workers)

    def test_missing_dependencies_point_to_settings_and_disable_start(self):
        for engine in (self.external, self.detect):
            with self.subTest(engine=engine):
                original = engine.return_value
                engine.return_value = None
                app = self.app()
                self.assertTrue(self.button(app).disabled)
                self.assertEqual(len(app.info), 1)
                self.assertIn('paquet complet', app.info[0].value)
                self.assertFalse(app.text_input)
                engine.return_value = original

    def test_start_failures_are_visible(self):
        with patch.object(ui, 'start_ffs', side_effect=ValueError('Guide incompatible')):
            app = self.app()
            self.button(app).click().run()
        self.assertFalse(app.exception)
        self.assertTrue(any('Guide incompatible' in item.value for item in app.error))

    def test_running_job_can_pause_with_missing_or_malformed_progress(self):
        job = self.running_job()
        for content in ('{', '[]', '{}'):
            (job.folder / 'progress.json').write_text(content)
            app = self.app()
            self.assertTrue(self.button(app).disabled)
            self.assertTrue(all(widget.disabled for widget in app.number_input))
            self.assertFalse(self.button(app, 'Mettre en pause').disabled)
        self.button(app, 'Mettre en pause').click().run()
        self.assertFalse(app.exception)
        job.cancel.assert_called_once()

    def test_live_movie_is_bound_to_the_selected_replica(self):
        job = self.running_job()
        job.request = {'kinetics': {'replicas': 4, 'source_samples': 200}}
        for number in (1, 2):
            folder = job.folder / f'replica_{number:03d}'
            folder.mkdir()
            (folder / 'live.dat').write_text('native frame')
            (folder / 'progress.json').write_text(json.dumps({'index': number-1, 'phase': 'flux',
                'flux': {'crossings': number}, 'observation': {'binding': number-1, 'distance': 8.}}))
        with patch.object(ui, 'render_oxdna_movie') as movie:
            app = self.app()
            self.assertEqual(movie.call_args.args[0]['movie_file'], 'replica_001/live.dat')
            self.assertEqual(movie.call_args.args[0]['replica_number'], 1)
            self.assertEqual(movie.call_args.args[0]['visualization'], ui.ffs_visualization(self.saved.snapshot))
            next(widget for widget in app.selectbox if widget.label == 'Réplique à visualiser').set_value(1).run()
            self.assertEqual(movie.call_args.args[0]['movie_file'], 'replica_002/live.dat')
            self.assertEqual(movie.call_args.args[0]['replica_number'], 2)
            self.assertEqual(next(w for w in app.selectbox if w.label == 'Réplique à visualiser').value, 1)
            self.assertFalse(any('Vue live · réplique' in item.value for item in app.markdown))
            dashboard = next(item for item in app.markdown if 'Suivi par réplique' in item.value)
            self.assertIn('Une ligne = une simulation indépendante', dashboard.proto.help)
            self.assertFalse(app.dataframe)
            self.assertFalse(app.get('vega_lite_chart'))  # No empty production charts during flux.
        job.cancel.assert_not_called()

    def test_modified_test_run_stays_exploratory_in_unified_dashboard(self):
        job = self.running_job()
        job.request = {'kinetics': {'replicas': 10},
                       'protocol_amendments': [{'reason': 'Changed distance criteria'}]}
        (job.folder / 'progress.json').write_text(json.dumps({'status': 'running',
            'kinetics': {'phase': 'flux', 'replicas': 10},
            'protocol_amendments': job.request['protocol_amendments']}))
        app = self.app()
        self.assertTrue(any('aucun taux cinétique ne sera estimé' in item.value for item in app.info))
        self.assertFalse(app.metric)
        self.assertEqual(len(app.table[0].value), 10)
        self.assertFalse(app.dataframe)
        self.assertFalse(self.button(app, 'Mettre en pause').disabled)

    def test_live_terminal_distance_uses_published_scope_and_keeps_separation_distinct(self):
        job = self.running_job()
        job.request = {'kinetics': {'replicas': 4}, 'plan': {'distance_scope': 'trigger_to_toehold'}}
        (job.folder / 'progress.json').write_text(json.dumps({
            'plan': {'distance_scope': 'trigger_3prime_to_toehold_start'},
            'observation': {'distance': 7., 'separation': 2.},
        }))
        folder = job.folder / 'replica_001'
        folder.mkdir()
        (folder / 'progress.json').write_text(json.dumps({
            'index': 0, 'phase': 'flux', 'observation': {'distance': 7.},
        }))
        app = self.app()
        self.assertIn('Distance 3′–toehold (nm)', app.table[0].value.columns)
        dashboard = next(item for item in app.markdown if 'Suivi par réplique' in item.value)
        self.assertIn('la dernière base du trigger côté 3′', dashboard.proto.help)
        self.assertIn('ne contraint pas leur orientation', dashboard.proto.help)
        self.assertIn('tous les nucléotides des deux brins sont éloignés', dashboard.proto.help)
        captions = [item.value for item in app.caption]
        self.assertTrue(any('Distance trigger 3′–début du toehold' in caption for caption in captions))
        self.assertTrue(any('Distance minimale entre les deux brins' in caption for caption in captions))
        job.cancel.assert_not_called()

    def test_previous_failure_is_hidden_while_another_run_is_active(self):
        job = self.running_job()
        self.results.return_value = [{'status': 'failed', 'run_id': 'old', 'run_dir': str(self.folder),
                                     'error': 'Old preparation failed', 'simulation_mode': 'opening'}]
        app = self.app()
        self.assertFalse(app.error)
        next(widget for widget in app.toggle if widget.label == 'Afficher les simulations précédentes').set_value(True).run()
        self.assertTrue(any('Old preparation failed' in item.value for item in app.error))
        self.assertTrue(any('ne décrivent pas le calcul en cours' in item.value for item in app.caption))
        job.cancel.assert_not_called()

    def test_recorded_movie_follows_replica_selection(self):
        run = self.interrupted_result()
        result = self.results.return_value[0]
        result.update(status='completed', resumable=False, movie_file='movie.dat',
                      requested_parameters={'kinetics': {'replicas': 4}})
        for number in (1, 2):
            folder = run / f'replica_{number:03d}'
            folder.mkdir()
            (folder / 'progress.json').write_text(json.dumps({'index': number-1, 'completed': True}))
            (folder / 'movie.dat').write_text('movie')
            (folder / 'movie_metadata.json').write_text(json.dumps({'final_interface_reached': True}))
        with patch.object(ui, 'render_oxdna_movie') as movie:
            app = self.app()
            self.assertEqual(movie.call_args.args[0]['movie_file'], 'replica_001/movie.dat')
            next(w for w in app.selectbox if w.label == 'Réplique à visualiser').set_value(1).run()
            self.assertEqual(movie.call_args.args[0]['movie_file'], 'replica_002/movie.dat')
            self.assertTrue(movie.call_args.args[0]['final_interface_reached'])

    def test_replica_failure_reads_terminal_result_but_ignores_old_result_after_resume(self):
        run = self.folder / 'failure'
        folder = run / 'replica_001'
        folder.mkdir(parents=True)
        (folder / 'progress.json').write_text(json.dumps({'index': 0, 'phase': 'production'}))
        (folder / 'result.json').write_text(json.dumps({'status': 'failed', 'error': 'Engine failed'}))
        request = {'kinetics': {'replicas': 4}}
        states = ui._replica_progress(run, request, {'status': 'failed'})
        self.assertEqual(states[0]['status'], 'failed')
        self.assertEqual(states[0]['error'], 'Engine failed')
        resumed = ui._replica_progress(run, request, {'status': 'running'})
        self.assertNotIn('status', resumed[0])
        self.assertNotIn('error', resumed[0])

    def test_more_than_four_replicas_and_workers_are_visible_and_bounded(self):
        with patch.object(ui.os, 'cpu_count', return_value=10):
            app = self.app()
            controls = {widget.label: widget for widget in app.number_input}
            controls['Répliques indépendantes'].set_value(8).run()
            next(w for w in app.number_input if w.label == 'Calculs simultanés').set_value(8).run()
            self.assertFalse(app.exception)
            self.assertEqual(next(w.value for w in app.number_input if w.label == 'Calculs simultanés'), 8)
            self.assertIn('8 répliques seront calculées', next(w for w in app.number_input if w.label == 'Calculs simultanés').proto.help)
            next(w for w in app.number_input if w.label == 'Répliques indépendantes').set_value(4).run()
            self.assertFalse(app.exception)
            self.assertLessEqual(next(w.value for w in app.number_input if w.label == 'Calculs simultanés'), 4)

    def test_navigation_and_lost_session_handle_do_not_cancel_the_job(self):
        job = self.running_job()
        script = '''
import streamlit as st
from aptaswitch_studio.saved_switch_ffs_ui import render_saved_switch_ffs
page = st.radio('Page', ['Simulation', 'Réglages'])
if page == 'Simulation':
    render_saved_switch_ffs(st.session_state['saved'])
else:
    st.write('Réglages')
'''
        app = self.app(script)
        app.radio[0].set_value('Réglages').run()
        self.assertFalse(app.exception)
        del app.session_state['saved_switch_ffs_job']
        app.radio[0].set_value('Simulation').run()
        self.assertFalse(app.exception)
        self.assertIs(app.session_state['saved_switch_ffs_job'], job)
        job.cancel.assert_not_called()

    def test_unlaunched_options_survive_leaving_the_lazy_section(self):
        script = '''
import streamlit as st
from aptaswitch_studio.saved_switch_ffs_ui import render_saved_switch_ffs
page = st.radio('Page', ['Simulation', 'Réglages'])
if page == 'Simulation':
    render_saved_switch_ffs(st.session_state['saved'])
'''
        app = self.app(script)
        values = {'Température de simulation (°C)': 42., 'Sel monovalent Na⁺': 300.,
                  'Franchissements initiaux par répétition': 150, 'Essais de production par interface': 750,
                  'Répliques indépendantes': 6, 'Calculs simultanés': min(2, ui.os.cpu_count() or 1)}
        for widget in app.number_input:
            widget.set_value(values[widget.label])
        app.run()
        self.assertFalse(app.exception)
        app.radio[0].set_value('Réglages').run()
        self.assertFalse(app.number_input)
        app.radio[0].set_value('Simulation').run()
        self.assertFalse(app.exception)
        self.assertEqual({widget.label: widget.value for widget in app.number_input}, values)

    def test_progress_counts_completed_replicas_instead_of_the_running_replica(self):
        job = self.running_job()
        for completed in (0, 1, 3):
            (job.folder / 'progress.json').write_text(json.dumps({'kinetics': {
                'phase': 'flux', 'completed_replicas': completed, 'replica': 4, 'replicas': 4}}))
            app = self.app()
            self.assertFalse(app.get('progress'))
            self.assertTrue(any(f'{completed}/4 répliques terminées' in c.value for c in app.caption))

    def test_resume_keeps_frozen_run_without_a_current_runtime(self):
        run = self.interrupted_result()
        self.external.return_value = ''
        self.detect.return_value = None
        job = SimpleNamespace(saved=self.saved)
        with patch.object(ui, 'resume_ffs', return_value=(job, True)) as resume, patch.object(ui, 'start_ffs') as start:
            app = self.app()
            next(widget for widget in app.number_input if widget.label == 'Température de simulation (°C)').set_value(40.).run()
            self.assertFalse(self.button(app, 'Reprendre la simulation').disabled)
            self.button(app, 'Reprendre la simulation').click().run()
        self.assertFalse(app.exception)
        resume.assert_called_once_with(self.saved, run)
        start.assert_not_called()

    def test_resume_unavailable_without_checkpoint_and_disabled_when_busy(self):
        self.interrupted_result(resumable=False)
        app = self.app()
        self.assertFalse(any(button.label == 'Reprendre la simulation' for button in app.button))
        self.results.return_value[0]['resumable'] = True
        self.busy.return_value = True
        app = self.app()
        self.assertTrue(self.button(app, 'Reprendre la simulation').disabled)

    def test_resume_race_and_errors_are_handled(self):
        self.interrupted_result()
        app = self.app()
        self.busy.return_value = True
        with patch.object(ui, 'resume_ffs') as resume:
            self.button(app, 'Reprendre la simulation').click().run()
        resume.assert_not_called()
        self.busy.return_value = False
        app.run()
        with patch.object(ui, 'resume_ffs', side_effect=ValueError('Sauvegarde illisible')):
            self.button(app, 'Reprendre la simulation').click().run()
        self.assertFalse(app.exception)
        self.assertTrue(any('Sauvegarde illisible' in item.value for item in app.error))

    def test_saved_movie_comes_before_rate_and_table_is_in_collapsed_details(self):
        self.interrupted_result()
        result = self.results.return_value[0]
        result.update(status='completed', resumable=False, movie_file='movie.dat',
                      stages=[{'name': 'opening', 'successes': 4, 'failures': 6, 'attempted': 10,
                               'target_trials': 10, 'completed': True}])
        events = []
        original = ui._render_kinetics
        def metrics(data):
            events.append('metrics')
            return original(data)
        with patch.object(ui, 'render_oxdna_movie', side_effect=lambda *args: events.append('film')), \
             patch.object(ui, '_render_kinetics', side_effect=metrics):
            app = self.app()
        self.assertEqual(events, ['film', 'metrics'])
        details = next(expander for expander in app.expander if expander.label == 'Détails du calcul')
        self.assertFalse(details.proto.expanded)
        self.assertEqual(len(details.table), 1)
        self.assertEqual(details.table[0].value.loc[0, 'Essais prévus'], 10)
        self.assertFalse(app.dataframe)

    def test_english_entry_controls_are_translated(self):
        from streamlit.testing.v1 import AppTest as RawAppTest
        app = RawAppTest.from_string('from aptaswitch_core.localization import set_language\nset_language("en")\n' + SCRIPT)
        app.session_state['saved'] = self.saved
        app.run()
        self.assertFalse(app.exception)
        self.assertEqual(app.subheader[0].value, 'Kinetics · approach, hybridization and opening')
        self.assertEqual([button.label for button in app.button], ['Start kinetic simulation'])
        self.assertEqual(app.radio[0].value, "opening")
        self.assertEqual([expander.label for expander in app.expander], ['Advanced parameters', 'Model and criteria'])

    def test_proximity_start_omits_flux_and_explains_zero_bonds(self):
        with patch.object(ui, 'start_ffs', return_value=(SimpleNamespace(saved=self.saved), True)) as start:
            app = self.app()
            app.radio[0].set_value('proximity').run()
            self.assertFalse(app.exception)
            self.assertEqual(len(app.number_input), 5)
            self.assertNotIn('Franchissements initiaux par répétition', [w.label for w in app.number_input])
            captions = ' '.join(item.value for item in app.caption)
            self.assertIn('sans aucune liaison interbrin', captions)
            self.assertIn('1,7–2,1 nm', captions)
            self.button(app, 'Lancer depuis le toehold').click().run()
        self.assertFalse(app.exception)
        options = start.call_args.args[2]
        self.assertEqual(options['kinetics_mode'], 'proximity')
        self.assertNotIn('kinetics_source_samples', options)
        self.assertEqual(options['kinetics_trials'], 250)

    def test_proximity_mode_survives_navigation_and_restores_complete_flux_setting(self):
        script = '''
import streamlit as st
from aptaswitch_studio.saved_switch_ffs_ui import render_saved_switch_ffs
page = st.radio('Page', ['Simulation', 'Réglages'])
if page == 'Simulation':
    render_saved_switch_ffs(st.session_state['saved'])
'''
        app = self.app(script)
        next(w for w in app.number_input if w.label == 'Franchissements initiaux par répétition').set_value(150).run()
        app.radio[1].set_value('proximity').run()
        app.radio[0].set_value('Réglages').run()
        app.radio[0].set_value('Simulation').run()
        self.assertFalse(app.exception)
        self.assertEqual(app.radio[1].value, 'proximity')
        app.radio[1].set_value('opening').run()
        self.assertEqual(next(w.value for w in app.number_input if w.label == 'Franchissements initiaux par répétition'), 150)

    def test_proximity_live_title_uses_frozen_mode_and_disables_mode_choice(self):
        job = self.running_job()
        job.request = {'kinetics': {'mode': 'proximity', 'replicas': 4}}
        (job.folder / 'progress.json').write_text(json.dumps({'kinetics': {'mode': 'proximity', 'phase': 'initialization'}}))
        app = self.app()
        self.assertTrue(app.radio[0].disabled)
        self.assertTrue(any(item.value == '**Départ proche du toehold**' for item in app.markdown))
        self.assertFalse(any(item.value == '**Simulation cinétique complète**' for item in app.markdown))
        job.cancel.assert_not_called()

    def test_pause_full_launch_proximity_then_resume_each_saved_run(self):
        full = self.running_job()
        full.request = {'kinetics': {'mode': 'opening'}}
        close_folder = self.folder / 'nearby'
        close_folder.mkdir()
        nearby = SimpleNamespace(saved=self.saved, folder=close_folder, run_id='nearby',
            request={'kinetics': {'mode': 'proximity'}}, cancel=Mock(),
            snapshot=lambda: JobSnapshot(state='running', current=0, total=4, message='FFS',
                                         elapsed_seconds=15, log=()))
        records = {}

        def pause(job, mode):
            self.active.return_value = None
            self.busy.return_value = False
            records[job.run_id] = {'status': 'cancelled', 'run_id': job.run_id, 'run_dir': str(job.folder),
                                  'simulation_mode': mode, 'resumable': True,
                                  'requested_parameters': {'kinetics': {'mode': mode}}}
            self.results.return_value = list(records.values())

        def activate(job):
            self.active.return_value = job
            self.busy.return_value = True
            records.pop(job.run_id, None)
            self.results.return_value = list(records.values())
            return job, True

        full.cancel.side_effect = lambda: pause(full, 'opening')
        nearby.cancel.side_effect = lambda: pause(nearby, 'proximity')
        app = self.app()
        self.button(app, 'Mettre en pause').click().run()
        app.run()
        self.assertFalse(app.radio[0].disabled)
        app.radio[0].set_value('proximity').run()
        with patch.object(ui, 'start_ffs', side_effect=lambda *args: activate(nearby)) as start:
            self.button(app, 'Lancer depuis le toehold').click().run()
        self.assertEqual(start.call_args.args[2]['kinetics_mode'], 'proximity')
        self.assertTrue(app.radio[0].disabled)
        self.assertEqual(records['running']['simulation_mode'], 'opening')
        self.button(app, 'Mettre en pause').click().run()
        app.run()
        self.assertEqual(set(records), {'running', 'nearby'})

        def resume(saved, folder):
            return activate(full if folder == full.folder else nearby)

        with patch.object(ui, 'resume_ffs', side_effect=resume) as resumed:
            next(w for w in app.selectbox if w.label == 'Simulation enregistrée').set_value(0).run()
            self.button(app, 'Reprendre la simulation').click().run()
            self.assertEqual(resumed.call_args.args, (self.saved, full.folder))
            self.assertTrue(any(m.value == '**Simulation cinétique complète**' for m in app.markdown))
            self.button(app, 'Mettre en pause').click().run()
            app.run()
            next(w for w in app.selectbox if w.label == 'Simulation enregistrée').set_value(0).run()
            self.button(app, 'Reprendre la simulation').click().run()
            self.assertEqual(resumed.call_args.args, (self.saved, nearby.folder))
            self.assertTrue(any(m.value == '**Départ proche du toehold**' for m in app.markdown))
        self.assertFalse(app.exception)
        self.assertEqual(records['running']['simulation_mode'], 'opening')

    def test_nearby_movie_curve_follows_selected_replica_metadata(self):
        run = self.interrupted_result()
        self.results.return_value[0].update(status='completed', resumable=False,
            requested_parameters={'kinetics': {'mode': 'proximity', 'replicas': 2}, 'sampling': {'dt': .003}})
        for number in (1, 2):
            folder = run / f'replica_{number:03d}'
            folder.mkdir()
            (folder / 'progress.json').write_text(json.dumps({'index': number - 1, 'completed': True}))
            (folder / 'movie.dat').write_text('movie')
            (folder / 'movie_metadata.json').write_text(json.dumps({'final_interface_reached': True,
                'trajectory_observables': [{'step': 0, 'binding': 0, 'stem': 15, 'opening': 4},
                                          {'step': 100 * number, 'binding': number, 'stem': 14, 'opening': 3}]}))
        with patch.object(ui, 'render_oxdna_movie'), \
             patch.object(ui, '_render_trajectory_observables', wraps=ui._render_trajectory_observables) as curve:
            app = self.app()
            self.assertEqual(curve.call_args.args[0]['trajectory_observables'][-1]['binding'], 1)
            self.assertEqual(len(app.get('vega_lite_chart')), 1)
            next(w for w in app.selectbox if w.label == 'Réplique à visualiser').set_value(1).run()
            self.assertEqual(curve.call_args.args[0]['trajectory_observables'][-1]['binding'], 2)
        self.assertFalse(app.exception)
        self.assertTrue(any('Trajectoire FFS réussie assemblée' in c.value for c in app.caption))

    def test_completed_replica_movie_and_curve_are_visible_while_global_run_continues(self):
        job = self.running_job()
        job.request = {'kinetics': {'mode': 'proximity', 'replicas': 4}, 'sampling': {'dt': .003}}
        folder = job.folder / 'replica_001'
        folder.mkdir()
        (folder / 'progress.json').write_text(json.dumps({'index': 0, 'completed': True}))
        (folder / 'live.dat').write_text('live')
        (folder / 'movie.dat').write_text('successful movie')
        (folder / 'movie_metadata.json').write_text(json.dumps({'final_interface_reached': True,
            'trajectory_observables': [{'step': 0, 'binding': 0, 'stem': 15, 'opening': 4},
                                      {'step': 100, 'binding': 30, 'stem': 0, 'opening': 0}]}))
        with patch.object(ui, 'render_oxdna_movie') as movie:
            app = self.app()
            self.assertEqual([call.args[0]['movie_file'] for call in movie.call_args_list],
                             ['replica_001/live.dat', 'replica_001/movie.dat'])
        self.assertEqual(len(app.get('vega_lite_chart')), 1)
        self.assertTrue(any('Trajectoire FFS réussie assemblée' in c.value for c in app.caption))
        self.assertFalse(self.button(app, 'Mettre en pause').disabled)
        self.assertFalse(any('Aucun film de trajectoire' in c.value for c in app.caption))
        job.cancel.assert_not_called()

    def test_partial_or_complete_mode_film_never_claims_a_successful_nearby_curve(self):
        run = self.interrupted_result()
        result = self.results.return_value[0]
        result.update(movie_file='movie.dat', final_interface_reached=False,
                      requested_parameters={'kinetics': {'mode': 'proximity'}, 'sampling': {'dt': .003}},
                      trajectory_observables=[{'step': 0, 'binding': 0, 'stem': 15, 'opening': 4},
                                              {'step': 100, 'binding': 1, 'stem': 14, 'opening': 3}])
        with patch.object(ui, 'render_oxdna_movie'):
            app = self.app()
            self.assertFalse(app.get('vega_lite_chart'))
            result['final_interface_reached'] = True
            app.run()
            self.assertEqual(len(app.get('vega_lite_chart')), 1)
            result['requested_parameters']['kinetics']['mode'] = 'opening'
            app.run()
            self.assertFalse(app.get('vega_lite_chart'))


if __name__ == '__main__':
    unittest.main()
