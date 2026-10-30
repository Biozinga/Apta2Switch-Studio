"""Per-replica summaries distinguish source collection, trials and lifecycle."""
from copy import deepcopy
import math
import unittest

from aptaswitch_core.localization import language_context
from aptaswitch_core.oxdna_kinetics import candidate_interfaces
from aptaswitch_studio.ffs_progress_charts import (
    distance_description, distance_label, replica_rows, selected_replica_chart, trajectory_observables_chart,
)


def request(mode='opening'):
    return {'kinetics': {'replicas': 4, 'mode': mode, 'source_samples': 200,
                        'trials': 1000, 'equilibration_steps': 200000},
            'plan': {'stem_pairs': [[i, 99-i] for i in range(15)],
                     'interfaces': [{'name': 'opening', 'min_binding': 24,
                                     'max_stem': 3, 'max_opening': 0}]}}


def stage(name='contact_7', *, successes=12, failures=30, target=1000, completed=False, pending=1):
    return {'name': name, 'target_trials': target, 'attempted': successes+failures+pending,
            'successes': successes, 'failures': failures, 'pending': pending, 'completed': completed}


def state(index=0, phase='production', **changes):
    value = {'replica_id': f'replica_{index+1}', 'index': index, 'phase': phase,
             'completed': False, 'burn_steps': 150000, 'production_index': 0,
             'flux': {'history': 'forward', 'crossings': 9, 'complete': False},
             'stages': [stage()] if phase == 'production' else [],
             'observation': {'binding': 0, 'distance': 7., 'stem': 15, 'opening': 2}}
    return {**value, **changes}


def chart_rows(chart):
    spec = chart.to_dict()
    data = spec['data']
    return data.get('values', spec.get('datasets', {}).get(data.get('name')))


class ReplicaSummaryTests(unittest.TestCase):
    def setUp(self):
        self.enterContext(language_context('fr'))

    def test_previous_distance_scope_remains_labelled_as_toehold_contact(self):
        current = request()
        current['plan']['distance_scope'] = 'trigger_to_toehold'
        row = replica_rows([state()], current)[0]
        self.assertEqual(row['Distance toehold (nm)'], '5.96')
        self.assertNotIn('Distance (nm)', row)

    def test_terminal_distance_labels_and_explanation_match_frozen_scope_in_both_languages(self):
        current = request()
        current['plan']['distance_scope'] = 'trigger_3prime_to_toehold_start'
        before = deepcopy(current)
        for language, heading, step, description in (
            ('fr', 'Distance 3′–toehold (nm)', '3′–toehold : ≤ 5.96 nm',
             'la première base du toehold après le leader'),
            ('en', '3′–toehold distance (nm)', '3′–toehold: ≤ 5.96 nm',
             'the first toehold base after the switch leader'),
        ):
            with self.subTest(language=language), language_context(language):
                row = replica_rows([state()], current)[0]
                self.assertEqual(row[heading], '5.96')
                self.assertIn(step, row.values())
                self.assertEqual(chart_rows(selected_replica_chart(state(), current))[0]['label'], step)
                self.assertIn(description, distance_description(current['plan']))
                self.assertIn('3′', distance_label(current['plan']))
        self.assertEqual(current, before)

    def test_legacy_or_unknown_scope_is_not_presented_as_terminal_distance(self):
        for plan in ({}, {'distance_scope': 'intended_partners'}, {'distance_scope': 'future'}):
            with self.subTest(plan=plan):
                self.assertEqual(distance_label(plan), 'Distance minimale entre partenaires cibles')
                self.assertEqual(distance_label(plan, column=True), 'Distance (nm)')
                self.assertEqual(distance_description(plan), '')

    def test_each_replica_has_its_own_phase_counts_and_queue_row(self):
        states = {0: state(0, 'equilibration'), 1: state(1, 'flux'), 2: state(2)}
        before = deepcopy(states)
        rows = replica_rows(states, request())
        self.assertEqual(len(rows), 4)
        self.assertEqual(len(rows[0]), 7)
        self.assertEqual([row['Réplique'] for row in rows], [1, 2, 3, 4])
        self.assertEqual([row['État'] for row in rows], ['Équilibrage', 'Collecte initiale', 'Essais FFS', 'En attente'])
        self.assertEqual([row['Avancement étape'] for row in rows],
                         ['150000/200000 pas', '9/200 franchissements', '42/1000 essais', '—'])
        self.assertEqual([row['Succès interface'] for row in rows], ['—', '—', '12', '—'])
        self.assertEqual(rows[2]['Étape'], 'Distance cible : ≤ 5.96 nm')
        self.assertEqual(rows[2]['Distance (nm)'], '5.96')
        self.assertEqual(rows[3]['Paires cibles'], '—')
        self.assertEqual(states, before)
        self.assertNotIn('%', str(rows))

    def test_flux_counts_geometry_crossings_without_calling_them_base_pairs(self):
        for history, expected in [('forward', 'Collecte initiale'), ('return', 'Retour vers l’état séparé'),
                                  ('B_history', 'Retour naturel depuis ON'), ('initial_A', 'Attente de l’état séparé')]:
            source = state(phase='flux', flux={'crossings': 199, 'history': history})
            row = replica_rows([source], request())[0]
            self.assertEqual(row['Étape'], expected)
            self.assertEqual(row['Avancement étape'], '199/200 franchissements')
            self.assertEqual(row['Paires cibles'], '0')
            self.assertEqual(row['Succès interface'], '—')
            self.assertIsNone(selected_replica_chart(source, request()))

    def test_termination_pause_failure_and_zero_success_remain_distinct(self):
        good = state(0, phase='completed', completed=True,
                     stages=[stage('opening', successes=10, failures=990, completed=True, pending=0)])
        zero = state(1, phase='completed', completed=True,
                     stages=[stage(successes=0, failures=1000, completed=True, pending=0)])
        active = state(2)
        rows = replica_rows([good, zero, active], request(), status='cancelled')
        self.assertEqual([row['État'] for row in rows], ['Terminée', 'Terminée sans succès', 'En pause', 'Non démarrée'])
        self.assertEqual(rows[1]['Avancement étape'], '1000/1000 essais')
        self.assertEqual(rows[1]['Succès interface'], '0')
        self.assertNotEqual(rows[1]['Étape'], 'Ouverture RBS')
        active['status'] = 'failed'
        self.assertEqual(replica_rows([active], request(), status='failed')[2]['État'], 'Échec')
        active.pop('status')
        self.assertEqual(replica_rows([active], request(), status='failed')[2]['État'], 'Interrompue')

    def test_sparse_indices_and_invalid_observations_do_not_invent_data(self):
        observation = {'binding': True, 'distance': math.nan}
        row = replica_rows([state(3, observation=observation)], request())[3]
        self.assertEqual(row['Paires cibles'], '—')
        self.assertEqual(row['Distance (nm)'], '—')
        rows = replica_rows({'replica_004': state(3)}, request())
        self.assertEqual(rows[0]['État'], 'En attente')
        self.assertEqual(rows[3]['État'], 'Essais FFS')
        self.assertEqual(replica_rows({}, {}), [])

    def test_new_distance_scope_is_labelled_as_toehold_contact(self):
        current = request()
        current['plan']['distance_scope'] = 'trigger_to_toehold'
        row = replica_rows([state()], current)[0]
        self.assertEqual(row['Distance toehold (nm)'], '5.96')
        self.assertNotIn('Distance (nm)', row)

    def test_proximity_initialization_has_no_flux_counter(self):
        current = request(mode='proximity')
        current['kinetics'].pop('source_samples')
        row = replica_rows([state(phase='initialization', flux={})], current)[0]
        self.assertEqual(row['État'], 'Préparation du départ proche')
        self.assertEqual(row['Étape'], 'Départ proche du toehold')
        self.assertEqual(row['Avancement étape'], '—')
        self.assertNotIn('franchissement', str(row))


class SelectedReplicaChartTests(unittest.TestCase):
    def setUp(self):
        self.enterContext(language_context('fr'))

    def test_one_success_visible_exact_counts_and_fixed_count_axis(self):
        current = state(stages=[stage(successes=1, failures=0)])
        chart = selected_replica_chart(current, request())
        rows, spec = chart_rows(chart), chart.to_dict()
        self.assertEqual([row['interface'] for row in rows], ['contact_7', 'contact_6'])
        self.assertEqual(rows[0]['caption'], '1 succès · 1/1000 essais')
        self.assertEqual(rows[0]['pending'], 1)
        self.assertEqual(rows[0]['label'], 'Distance cible : ≤ 5.96 nm')
        self.assertEqual(spec['encoding']['x']['scale']['domain'], [0, 1000])
        self.assertTrue({'point', 'bar', 'text'} <= {layer['mark']['type'] for layer in spec['layer']})
        self.assertEqual(spec['layer'][1]['encoding']['color']['scale']['domain'], ['En cours', 'À venir'])
        self.assertNotIn('title', spec)
        self.assertNotIn('contact_8', [row['interface'] for row in rows])

    def test_master_published_interfaces_are_used_without_reconstructing_a_new_grid(self):
        frozen = request()
        frozen['plan']['interfaces'] = candidate_interfaces(frozen['plan'])
        frozen['plan']['interfaces'][1]['distance_max'] = 7.25
        rows = chart_rows(selected_replica_chart(state(), frozen))
        self.assertEqual(rows[0]['label'], 'Distance cible : ≤ 6.18 nm')

    def test_production_chart_never_pools_another_replica_or_flux_events(self):
        first, second = state(), state(1, stages=[stage(successes=600, failures=3)])
        first['flux']['crossings'] = 200
        self.assertEqual(chart_rows(selected_replica_chart(first, request()))[0]['successes'], 12)
        self.assertEqual(chart_rows(selected_replica_chart(second, request()))[0]['successes'], 600)
        self.assertEqual(chart_rows(selected_replica_chart(first, request()))[0]['finished'], 42)

    def test_recent_window_is_bounded_and_preserves_actual_predicate_labels(self):
        interfaces = candidate_interfaces(request()['plan'])
        started = [stage(item['name'], successes=30, failures=970, pending=0, completed=True)
                   for item in interfaces[1:]]
        chart = selected_replica_chart(state(completed=True, stages=started), request())
        rows, spec = chart_rows(chart), chart.to_dict()
        self.assertEqual(len(rows), 10)
        self.assertLessEqual(spec['height'], 320)
        self.assertEqual(rows[-1]['label'], 'Ouverture RBS')
        self.assertIn('Tige : ≤ 4 pb', [row['label'] for row in rows])
        self.assertIn('RBS : ≤ 1 liaison', [row['label'] for row in rows])
        at_duplex = next(i for i, item in enumerate(interfaces) if item['name'] == 'duplex_3')
        chart = selected_replica_chart(state(stages=started[:at_duplex]), request())
        self.assertIn('Duplex : 3 pb', [row['label'] for row in chart_rows(chart)])

    def test_zero_success_stops_are_not_rendered_as_opening_or_pending(self):
        current = state(completed=True, phase='completed',
                        stages=[stage(successes=0, failures=1000, pending=0, completed=True)])
        rows = chart_rows(selected_replica_chart(current, request()))
        self.assertEqual(rows[0]['state'], 'Sans succès')
        self.assertEqual(rows[0]['caption'], '0 succès · 1000/1000 essais')
        self.assertEqual(rows[1]['state'], 'Non atteinte')
        self.assertEqual(rows[1]['finished'], 0)
        self.assertNotIn('opening', [row['interface'] for row in rows])

    def test_pause_is_not_displayed_as_active_and_incomplete_not_green(self):
        rows = chart_rows(selected_replica_chart(state(), request(), status='cancelled'))
        self.assertEqual(rows[0]['state'], 'En pause')
        current = state(stages=[stage(successes=30, failures=970, pending=0, completed=False)])
        self.assertEqual(chart_rows(selected_replica_chart(current, request()))[0]['state'], 'En cours')

    def test_proximity_graph_preserves_near_source_and_starts_at_contact_two(self):
        frozen = request(mode='proximity')
        frozen['plan']['interfaces'] = [
            {'name': 'near_toehold', 'distance_max': 2.5, 'min_binding': 0},
            {'name': 'contact_2', 'distance_max': 2., 'min_binding': 0},
            {'name': 'duplex_1', 'min_binding': 1},
            frozen['plan']['interfaces'][-1]]
        current = state(stages=[stage('contact_2')])
        rows = chart_rows(selected_replica_chart(current, frozen))
        self.assertEqual([row['interface'] for row in rows], ['contact_2', 'duplex_1'])
        self.assertEqual(rows[0]['successes'], 12)
        self.assertEqual(replica_rows([current], frozen)[0]['Étape'], 'Distance cible : ≤ 1.70 nm')


if __name__ == '__main__':
    unittest.main()


class TrajectoryObservablesChartTests(unittest.TestCase):
    def setUp(self):
        self.enterContext(language_context('fr'))
        self.points = [{'step': 100, 'binding': 0, 'stem': 15, 'opening': 4, 'distance': 2.3},
                       {'step': 200, 'binding': 1, 'stem': 14, 'opening': 3, 'distance': 1.},
                       {'step': 400, 'binding': 30, 'stem': 0, 'opening': 0, 'distance': .8}]

    def test_counts_follow_retained_frames_in_relative_model_time(self):
        before = deepcopy(self.points)
        chart = trajectory_observables_chart(self.points, {'dt': .003})
        spec, rows = chart.to_dict(), chart_rows(chart)
        self.assertEqual(spec['encoding']['x']['title'], 'Temps oxDNA le long du trajet FFS')
        self.assertEqual(len(rows), 9)
        binding = [row for row in rows if row['series'] == 'Paires trigger–switch']
        self.assertEqual([row['bonds'] for row in binding], [0, 1, 30])
        self.assertEqual([round(row['time'], 6) for row in binding], [0., .3, .9])
        self.assertEqual({row['series'] for row in rows},
                         {'Paires trigger–switch', 'Paires de la tige OFF', 'Liaisons RBS–AUG'})
        self.assertEqual(self.points, before)

    def test_invalid_or_missing_trajectory_never_invents_a_curve(self):
        for points, sampling in [(None, {'dt': .003}), ([], {'dt': .003}),
                                 (self.points[:1], {'dt': .003}),
                                 (self.points[::-1], {'dt': .003}),
                                 (self.points, {}), (self.points, {'dt': float('nan')}),
                                 (self.points, {'dt': 0.}),
                                 ([self.points[0], {**self.points[1], 'binding': True}], {'dt': .003})]:
            self.assertIsNone(trajectory_observables_chart(points, sampling))

    def test_trajectory_labels_translate(self):
        with language_context('en'):
            chart = trajectory_observables_chart(self.points, {'dt': .003})
            self.assertEqual(chart.to_dict()['encoding']['x']['title'], 'oxDNA time along the FFS path')
            self.assertIn('OFF stem base pairs', {row['series'] for row in chart_rows(chart)})
