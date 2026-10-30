from pathlib import Path
import json
import shutil
import subprocess
import tempfile
import unittest
from unittest.mock import patch

from aptaswitch_core.localization import language_context
from aptaswitch_studio.oxdna_movie import _contained, _legend_groups, _legend_html, read_movie, movie_html, apply_domain_colors


def frame(step, count=4):
    return f't = {step}\nb = 40 40 40\nE = 0 0 0\n' + ''.join(
        f'{i} 0 0 1 0 0 0 0 1 0 0 0 0 0 0\n' for i in range(count))


class OxDNAMovieTests(unittest.TestCase):
    def setUp(self):
        self.enterContext(language_context('fr'))
        self.root = Path(self.enterContext(tempfile.TemporaryDirectory()))
        self.top = self.root / 'system.top'
        self.traj = self.root / 'movie.dat'
        self.top.write_text('4 2 5->3\nAU type=RNA\nAT type=DNA\n')
        self.traj.write_text(''.join(frame(i) for i in range(10)))
        self.exports = self.enterContext(patch('aptaswitch_studio.oxdna_movie.render_export_menu'))

    def test_downsampling_preserves_actual_first_and_final_frames(self):
        movie = read_movie(self.top, self.traj, max_frames=4)
        self.assertEqual(movie['steps'], [0, 3, 6, 9])
        self.assertEqual(movie['source_frames'], 10)
        self.assertEqual(movie['strands'], [0, 0, 1, 1])
        self.assertEqual(movie['sequence'], list('AUAT'))
        self.assertEqual(movie['frames'][0][3], [3, 0, 0, 1, 0, 0, 0, 0, 1])
        self.assertEqual(movie['chemistries'], ['RNA', 'RNA', 'DNA', 'DNA'])

    def test_truncated_or_nonfinite_trajectory_is_not_rendered(self):
        for text in (frame(0)[:-30], frame(0).replace('1 0 0', 'nan 0 0', 1), 'not an oxDNA frame'):
            self.traj.write_text(text)
            with self.subTest(text=text[:30]), self.assertRaises(ValueError):
                read_movie(self.top, self.traj)

    def test_movie_is_standalone_and_exposes_actual_frame_controls_and_video_export(self):
        html = movie_html(read_movie(self.top, self.traj), reached=False)
        self.assertNotIn('<script src=', html)
        self.assertIn('MediaRecorder', html)
        self.assertIn('film.frames[index]', html)
        self.assertIn('type="range"', html)
        self.assertIn('3Dmol', html)

    def test_movie_files_cannot_escape_run_or_traverse_symlink(self):
        self.assertEqual(_contained(self.root, 'movie.dat'), self.traj)
        link = self.root / 'linked.dat'
        link.symlink_to(self.traj)
        for value in ('../outside.dat', 'missing.dat', 'linked.dat'):
            with self.subTest(value=value), self.assertRaises(ValueError):
                _contained(self.root, value)

    def test_opening_movie_is_recognized_while_ensemble_is_still_running_or_paused(self):
        from aptaswitch_studio import oxdna_movie as viewer
        for status in ('running', 'cancelled', 'incomplete'):
            with self.subTest(status=status), patch.object(viewer, 'st'), patch.object(viewer, 'components'), \
                    patch.object(viewer, 'movie_html', wraps=viewer.movie_html) as render:
                viewer.render_oxdna_movie({'status': status, 'final_interface_reached': True,
                    'movie_file': 'movie.dat', 'topology_file': 'system.top'}, self.root, 'observed_on')
                self.assertTrue(render.call_args.kwargs['reached'])

    def test_movie_endpoint_distinguishes_first_contact_from_opening_in_both_languages(self):
        movie = read_movie(self.top, self.traj)
        for language, contact, opening in (
            ('fr', 'premier appariement cible atteint', 'critère d’ouverture'),
            ('en', 'first intended base pair reached', 'opening criterion'),
        ):
            with self.subTest(language=language), language_context(language):
                html = movie_html(movie, endpoint='first_intended_base_pair', reached=True)
                self.assertIn(contact, html)
                self.assertNotIn(opening, html)
                html = movie_html(movie, endpoint='opening', reached=True)
                self.assertIn(opening, html)
                self.assertNotIn(contact, html)

    def test_completed_sampling_does_not_claim_endpoint_was_reached(self):
        from aptaswitch_studio import oxdna_movie as viewer
        for flag in (None, False, 'true'):
            with self.subTest(flag=flag), patch.object(viewer, 'st'), patch.object(viewer, 'components') as components:
                viewer.render_oxdna_movie({'status': 'completed', 'final_interface_reached': flag,
                    'trajectory_endpoint': 'opening', 'movie_file': 'movie.dat',
                    'topology_file': 'system.top'}, self.root, 'not_opened')
                html = components.html.call_args.args[0]
                self.assertIn('ouverture non atteinte', html)
                self.assertNotIn('ayant atteint le critère', html)

    def test_result_endpoint_metadata_controls_rendered_title(self):
        from aptaswitch_studio import oxdna_movie as viewer
        for metadata in ({'trajectory_endpoint': 'first_intended_base_pair'},
                         {'kinetics': {'endpoint': 'first_intended_base_pair'}},
                         {'trajectory_endpoint': 'first_intended_base_pair', 'kinetics': {'endpoint': 'opening'}}):
            with self.subTest(metadata=metadata), patch.object(viewer, 'st'), patch.object(viewer, 'components') as components:
                viewer.render_oxdna_movie({'status': 'completed', 'final_interface_reached': True,
                    'movie_file': 'movie.dat', 'topology_file': 'system.top', **metadata}, self.root, 'contact')
                html = components.html.call_args.args[0]
                self.assertIn('premier appariement cible atteint', html)
                self.assertNotIn('critère d’ouverture', html)

    def test_live_and_preparation_titles_take_priority_over_endpoint(self):
        movie = read_movie(self.top, self.traj)
        for endpoint in ('opening', 'first_intended_base_pair'):
            for options, title in (({'live': True}, 'Vue 3D en direct'),
                                   ({'preparation': True}, 'Switch et trigger · configuration gros grain préparée')):
                with self.subTest(endpoint=endpoint, options=options):
                    html = movie_html(movie, endpoint=endpoint, reached=True, **options)
                    self.assertIn(title, html)
                    self.assertNotIn('premier appariement cible atteint', html)
                    self.assertNotIn('critère d’ouverture', html)

    def test_prepared_configuration_has_a_static_view_without_film_controls(self):
        self.traj.write_text(frame(0))
        html = movie_html(read_movie(self.top, self.traj), preparation=True)
        self.assertIn('id="playback" hidden', html)
        self.assertIn('id="info" hidden', html)
        self.assertIn('id="fit">Recentrer</button>', html)
        self.assertNotIn('<script src=', html)
        self.assertIn('viewer', html)

    def test_live_and_movie_share_sequence_checked_colors_and_legend(self):
        visualization = {'sequence': list('AUAT'), 'colors': ['#33b5e5', '#00a676', '#ffd54f', '#8b5cf6'],
                         'domains': ['Toehold', 'AUG', 'Trigger', 'Extension 3′']}
        movie = apply_domain_colors(read_movie(self.top, self.traj), visualization)
        for html in (movie_html(movie, live=True, view_key='live-run'), movie_html(movie)):
            self.assertIn('film.colors', html)
            self.assertIn('class="oxdna-legend"', html)
            for color in visualization['colors']:
                self.assertIn(color, html)
        live = movie_html(movie, live=True, view_key='live-run')
        self.assertIn('id="playback" hidden', live)
        self.assertIn('oxdna-camera-', live)
        self.assertIn('sessionStorage', live)
        with self.assertRaises(ValueError):
            apply_domain_colors(movie, {**visualization, 'sequence': list('CCAT')})

    def test_legend_groups_follow_strands_even_when_colors_match(self):
        visualization = {'sequence': list('AUAT'), 'colors': ['#33b5e5', '#00a676', '#00a676', '#8b5cf6'],
                         'domains': ['Toehold', 'AUG', 'Extension 5′', 'Extension 3′']}
        movie = apply_domain_colors(read_movie(self.top, self.traj), visualization)
        self.assertEqual(_legend_groups(movie), [
            ('Toehold switch', [('Toehold', '#33b5e5'), ('AUG', '#00a676')]),
            ('Trigger', [('Extension 5′', '#00a676'), ('Extension 3′', '#8b5cf6')]),
        ])
        movie['domains'][0] = '<img src=x onerror=alert(1)>'
        legend = _legend_html(movie)
        self.assertNotIn('<img', legend)
        self.assertIn('&lt;img', legend)

    def test_embedded_view_excludes_legend_and_header_while_export_preserves_them(self):
        movie = read_movie(self.top, self.traj)
        standalone = movie_html(movie, live=True, replica_number=3)
        embedded = movie_html(movie, live=True, replica_number=3, embedded=True)
        self.assertIn('<h3>Vue 3D en direct · Réplique 3</h3>', standalone)
        self.assertIn('<section class="oxdna-legend"', standalone)
        self.assertNotIn('<h3>', embedded)
        self.assertNotIn('<section class="oxdna-legend"', embedded)
        for html in (embedded, standalone):
            self.assertNotIn('<script src=', html)
            payload = html.split('const film=', 1)[1].split(', L=', 1)[0]
            actual = json.loads(payload)
            for key in ('sequence', 'strands', 'chemistries', 'frames', 'boxes', 'steps'):
                self.assertEqual(actual[key], movie[key])

    def test_renderer_places_legend_outside_non_scrolling_view_and_exports_standalone(self):
        from aptaswitch_studio import oxdna_movie as viewer
        for kind, height in (('live', 540), ('prepared', 500), ('film', 620)):
            with self.subTest(kind=kind), patch.object(viewer, 'st') as st, patch.object(viewer, 'components') as components:
                self.exports.reset_mock()
                viewer.render_oxdna_movie({'movie_kind': kind, 'movie_file': 'movie.dat',
                    'topology_file': 'system.top', 'replica_number': 4}, self.root, 'run_replica_4')
                self.assertEqual(components.html.call_args.kwargs, {'height': height, 'scrolling': False})
                self.assertNotIn('<section class="oxdna-legend"', components.html.call_args.args[0])
                self.assertIn('<section class="oxdna-legend"', st.html.call_args.args[0])
                self.assertIn('Réplique 4', st.markdown.call_args.args[0])
                self.assertTrue(st.markdown.call_args.kwargs['help'])
                if kind == 'live':
                    self.exports.assert_not_called()
                else:
                    options = self.exports.call_args.args[0]
                    self.assertEqual(len(options), 3)
                    exported = options[0].data.decode()
                    self.assertIn('<h3>', exported)
                    self.assertIn('<section class="oxdna-legend"', exported)
                    self.assertNotIn('class="embedded"', exported)
                    self.assertEqual(options[1].data, self.traj.read_bytes())
                    self.assertEqual(options[2].data, self.top.read_bytes())

    @unittest.skipUnless(shutil.which('node'), 'Node.js is required for the viewer controller test')
    def test_live_camera_is_saved_initially_restored_per_view_and_reset_on_demand(self):
        movie = read_movie(self.top, self.traj)
        def controller(key):
            return movie_html(movie, live=True, view_key=key).rsplit('<script>', 1)[1].split('</script>', 1)[0]
        harness = r'''
const fs=require('fs'),vm=require('vm'),input=JSON.parse(fs.readFileSync(0,'utf8')),storage=new Map();
const initial=[1,2,3,50,0,0,0,1],custom=[9,8,7,40,.1,.2,.3,.9];
function render(code){
 let camera=initial.slice();const nodes={};
 const element=id=>nodes[id]||(nodes[id]={value:id==='speed'?12:0,handlers:{},addEventListener(event,fn){this.handlers[event]=fn;}});
 const viewer={addModel:()=>({addAtoms(){}}),removeAllModels(){},removeModel(){},setStyle(){},render(){},resize(){},
  zoomTo(){camera=initial.slice();},getView:()=>camera.slice(),setView(value){camera=value.slice();}};
 const context={document:{getElementById:element},$3Dmol:{createViewer:()=>viewer},
  sessionStorage:{getItem:key=>storage.get(key),setItem:(key,value)=>storage.set(key,value)},
  ResizeObserver:class{observe(){}},setTimeout:fn=>fn(),setInterval:()=>1,clearInterval(){}};
 vm.runInNewContext(code,context);return {viewer,nodes};
}
const first=render(input.script),storedInitial=JSON.parse(storage.get('oxdna-camera-run1-replica1'));
first.viewer.setView(custom);first.nodes.view.handlers.wheel();
const restored=render(input.script),restoredView=restored.viewer.getView();
const other=render(input.other),otherView=other.viewer.getView();
restored.nodes.fit.onclick();
process.stdout.write(JSON.stringify({initial:storedInitial,restored:restoredView,other:otherView,
 reset:JSON.parse(storage.get('oxdna-camera-run1-replica1')),info:restored.nodes.info.textContent}));
'''
        result = subprocess.run([shutil.which('node'), '-e', harness],
                                input=json.dumps({'script': controller('run1-replica1'), 'other': controller('run1-replica2')}),
                                text=True, capture_output=True, check=True, timeout=10)
        observed = json.loads(result.stdout)
        self.assertEqual(observed['initial'], [1, 2, 3, 50, 0, 0, 0, 1])
        self.assertEqual(observed['restored'], [9, 8, 7, 40, .1, .2, .3, .9])
        self.assertEqual(observed['other'], observed['initial'])
        self.assertEqual(observed['reset'], observed['initial'])
        self.assertEqual(observed['info'], 'Pas MD 0')

    def test_free_md_movie_has_its_own_title_and_export_names(self):
        html = movie_html(read_movie(self.top, self.traj), endpoint='free_md')
        self.assertIn('MD libre ADN–ADN', html)
        self.assertIn('oxDNA_free_MD.', html)
        self.assertNotIn('Trajectoire FFS', html)
        self.assertNotIn('oxDNA_FFS.', html)


if __name__ == '__main__':
    unittest.main()
