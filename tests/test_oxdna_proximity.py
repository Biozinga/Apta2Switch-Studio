"""Conditional proximity sampling: real bookkeeping, resume and movie provenance."""
from copy import deepcopy
from pathlib import Path
import json
import tempfile
import unittest
from unittest.mock import Mock, patch

from aptaswitch_core.oxdna_ffs import _finish_result, load_ffs_resume_request
from aptaswitch_core.oxdna_ffs_worker import Observation, SamplingStopped, request_identity
from aptaswitch_core.oxdna_kinetics import (candidate_interfaces, proximity_ffs_report,
                                         source_fingerprint, validate_kinetic_protocol)
from aptaswitch_core.oxdna_kinetics_worker import FFSReplica, KineticSampler, read_json
import test_oxdna_kinetics as _kinetics
import test_oxdna_kinetics_integration as _integration
from test_oxdna_kinetics import stage
from test_oxdna_kinetics_worker import NativeScript, frame, request


def proximity_options(options):
    return {**options, 'mode': 'proximity', 'source_samples': 0, 'equilibration_steps': 0,
            'initial_distance_min': 2., 'initial_distance_max': 2.5}


class ProximityRequestTests(unittest.TestCase):
    """Retain the full integration contract while testing the additional mode."""
    setUp = _integration.KineticsIntegrationTests.setUp
    request = _integration.KineticsIntegrationTests.request

    def test_close_mode_has_no_flux_and_preserves_terminal_criterion(self):
        full = self.request()
        close = self.request({'kinetics_mode': 'proximity', 'kinetics_source_samples': 0})
        validate_kinetic_protocol(close['kinetics'])
        self.assertEqual(close['kinetics']['source_samples'], 0)
        self.assertEqual(close['kinetics']['equilibration_steps'], 200000)
        self.assertEqual(close['kinetics']['version'], 4)
        self.assertFalse(close['plan']['rate_estimation'])
        self.assertEqual(close['plan']['interfaces'], full['plan']['interfaces'])
        self.assertEqual(close['box_policy'], full['box_policy'])
        interfaces = candidate_interfaces(close['plan'], mode='proximity')
        self.assertEqual([row['name'] for row in interfaces[:3]], ['near_toehold', 'contact_2', 'duplex_1'])
        self.assertEqual(interfaces[-1], full['plan']['interfaces'][-1])
        self.assertEqual(close['kinetics']['trials'], full['kinetics']['trials'])
        with self.assertRaises(ValueError):
            self.request({'kinetics_mode': 'encounter'})

    def test_proximity_protocol_cannot_silently_change_distance_or_add_flux(self):
        options = self.request({'kinetics_mode': 'proximity'})['kinetics']
        for changes in ({'initial_distance_min': 1.}, {'initial_distance_max': 8.},
                        {'source_samples': 75}, {'equilibration_steps': 0},
                        {'source_samples': False}):
            with self.subTest(changes=changes), self.assertRaises(ValueError):
                validate_kinetic_protocol(options | changes)

    def test_finalization_keeps_conditional_result_and_mode(self):
        req = self.request({'kinetics_mode': 'proximity'})
        report = {'status': 'conditional_estimate', 'value': .2, 'unit': 'probability'}
        result = {'kinetics': {'mode': 'proximity', 'conditional_probability': report},
                  'conditional_probability': report, 'rate_estimate': None}
        finished = _finish_result(result, req, self.root, 0., None)
        self.assertEqual(finished['simulation_mode'], 'proximity')
        self.assertEqual(finished['conditional_probability'], report)
        self.assertIsNone(finished['rate_estimate'])

    def test_proximity_resume_keeps_mode_and_frozen_geometry(self):
        req = self.request({'kinetics_mode': 'proximity'})
        preparation = self.root/'preparation'
        preparation.mkdir()
        for filename in ('validated.dat', 'system.top'):
            (preparation/filename).write_text('fixture')
        req['prepared'] = {'configuration_file': str(preparation/'validated.dat'),
                           'topology_file': str(preparation/'system.top')}
        (self.root/'worker_request.json').write_text(json.dumps(req))
        resumed = load_ffs_resume_request(self.root)
        self.assertEqual(resumed['kinetics'], req['kinetics'])


class ProximityReportTests(unittest.TestCase):
    def rows(self, successes=(20, 40, 60, 80)):
        return [{'replica_id': str(i), 'completed': True, 'initial_source': 'seed',
                 'stages': [stage(count), stage(50)]} for i, count in enumerate(successes)]

    def test_averages_whole_source_probabilities_without_a_flux_or_rate(self):
        report = proximity_ffs_report(self.rows())
        self.assertEqual(report['status'], 'conditional_estimate')
        self.assertAlmostEqual(report['value'], .25)
        self.assertEqual(report['unit'], 'probability')
        self.assertFalse(report['initial_condition']['equilibrium_sample'])
        self.assertLess(report['interval95'][0], report['value'])
        self.assertGreater(report['interval95'][1], report['value'])
        for field in ('flux', 'value_model', 'inverse_model_time', 'experimental_rate', 'value_m_inv_model_time'):
            self.assertNotIn(field, report)
        scaled = self.rows()
        for row in scaled:
            for s in row['stages']:
                for name in ('successes', 'failures', 'attempted', 'target_trials'):
                    s[name] *= 100
        self.assertEqual(proximity_ffs_report(scaled)['interval95'], report['interval95'])

    def test_pending_unprepared_missing_and_all_zero_sources_are_not_success(self):
        for modification in ({'completed': False}, {'initial_source': None}, {'stages': [stage(50, pending=1)]}):
            rows = self.rows()
            rows[-1].update(modification)
            result = proximity_ffs_report(rows)
            self.assertEqual(result['status'], 'incomplete')
            self.assertIsNone(result['value'])
        rows = self.rows((0, 0, 0, 0))
        for row in rows:
            row['stages'] = row['stages'][:1]
        report = proximity_ffs_report(rows)
        self.assertEqual(report['status'], 'no_success')
        self.assertEqual(report['value'], 0.)
        self.assertIsNone(report['interval95'])
        self.assertEqual(proximity_ffs_report([])['status'], 'incomplete')

    def test_protocol_amendments_and_duplicate_sources_are_rejected(self):
        report = proximity_ffs_report(self.rows(), protocol_amendments=[{'reason': 'changed source'}])
        self.assertEqual(report['status'], 'exploratory_modified_protocol')
        self.assertIsNone(report['value'])
        self.assertTrue(all(row['probability'] is None for row in report['replicas']))
        rows = self.rows()
        rows[1]['replica_id'] = rows[0]['replica_id']
        with self.assertRaises(ValueError):
            proximity_ffs_report(rows)

    def test_source_identity_distinguishes_the_new_preparation(self):
        full = _kinetics.SourceProvenanceTests().request()
        close = deepcopy(full)
        close['kinetics'] = proximity_options(close['kinetics'])
        self.assertNotEqual(source_fingerprint(full), source_fingerprint(close))
        changed = deepcopy(close)
        changed['kinetics']['initial_distance_max'] = 3.
        self.assertNotEqual(source_fingerprint(close), source_fingerprint(changed))


class ProximityWorkerTests(unittest.TestCase):
    def setUp(self):
        self.root = Path(self.enterContext(tempfile.TemporaryDirectory()))
        self.req = request(self.root)
        self.req['kinetics'] = proximity_options(self.req['kinetics'])
        self.req['plan'].update(distance_scope='trigger_to_toehold', distance_pairs=[[0, 2], [0, 3]],
                                separation_pairs=[[0, 2], [0, 3], [1, 2], [1, 3]], rate_estimation=False)

    def replica(self, events=(), *, observation=None):
        self.script = NativeScript(events)
        def adapter(*args, **kwargs):
            value = self.script(*args, **kwargs)
            value.observe = lambda: observation or Observation(0, 1, 1, 2.25, value.step, 0, 2.25)
            return value
        return FFSReplica(self.req, self.root/'replica_001', 0, adapter_factory=adapter)

    def test_native_seed_is_unbound_and_flux_is_never_simulated(self):
        r = self.replica()
        r.prepare_proximity()
        self.assertEqual(r.state['steps'], 0)
        self.assertEqual(r.state['burn_steps'], 0)
        self.assertEqual(r.state['flux']['bank'], [])
        self.assertEqual(r.state['flux']['exposure_steps'], 0)
        self.assertEqual(r.state['production_bank'], [r.state['initial_source']])
        self.assertEqual(r.state['phase'], 'production')
        self.assertEqual(r.state['observation']['associations'], 0)
        self.assertTrue(r.state['initialization']['native_interstrand_bonds_checked'])
        self.assertTrue(self.script.starts[0]['initial'])
        with self.assertRaises(ValueError):
            r.flux()

    def test_prebound_seed_is_rejected_before_any_production(self):
        r = self.replica(observation=Observation(1, 1, 1, 2.25, 0, 1, 2.25))
        with self.assertRaisesRegex(ValueError, 'zero interstrand'):
            r.prepare_proximity()
        self.assertIsNone(r.state['initial_source'])
        self.assertEqual(r.state['production_bank'], [])

    def test_reversal_and_dissociation_count_as_failures_without_retry_until_success(self):
        r = self.replica([(3, 'auto', 0, 9, 9), (5, 'auto', 0, 9, 9)])
        result = r.run()
        self.assertEqual(result['status'], 'completed')
        self.assertEqual(result['stages'][0]['failures'], 2)
        self.assertEqual(result['stages'][0]['successes'], 0)
        self.assertEqual(result['steps'], 8)
        self.assertIsNone(result['movie'])
        self.assertEqual(proximity_ffs_report([result])['status'], 'no_success')

    def test_pause_resume_keeps_source_parent_observation_and_pending_trial(self):
        r = self.replica([(5, None, 0, 3., 3.)])
        r.prepare_proximity()
        self.script.on_advance = lambda: (self.root/'CANCEL').touch()
        with self.assertRaises(SamplingStopped):
            r.production()
        r.save(force=True)
        before = read_json(r.output/'state.json')
        (self.root/'CANCEL').unlink()
        with patch('aptaswitch_core.oxdna_kinetics_worker.random_proximity_configuration',
                   side_effect=AssertionError('Do not replace a paused source')):
            resumed = self.replica([(3, 'auto', 0, 9, 9), (4, 'auto', 0, 9, 9)])
        resumed.prepare_proximity()
        self.assertEqual(self.script.starts, [])
        self.assertEqual(resumed.state['active'], before['active'])
        self.assertEqual(resumed.state['frame'], before['frame'])
        self.assertEqual(resumed.state['initial_source'], before['initial_source'])
        self.assertEqual(resumed.rng.getstate(), r.rng.getstate())
        resumed.production()
        self.assertFalse(self.script.starts[0]['initial'])
        self.assertEqual(resumed.state['stages'][0]['attempted'], 2)
        self.assertEqual(resumed.state['stages'][0]['failures'], 2)

    def test_movie_observations_follow_successful_lineage_and_same_thinning(self):
        self.req['sampling'].update(frame_interval=1, max_movie_frames=5)
        r = self.replica()
        count = len(r.interfaces)-1
        self.script.events = ([(1, None, 0, 2.25, 2.25)] * 100
                              + [(2, 'auto', 2, 1., 1.)] * (2*count))
        result = r.run()
        self.assertEqual(result['status'], 'completed')
        metadata = read_json(r.output/'movie_metadata.json')
        frames = (r.output/'movie.dat').read_text().split('t =')[1:]
        steps = [int(segment.splitlines()[0].strip()) for segment in frames]
        observations = metadata['trajectory_observables']
        self.assertEqual(steps, [row['step'] for row in observations])
        self.assertEqual(len(steps), 5)
        self.assertEqual(observations[0]['binding'], 0)
        self.assertEqual(observations[-1]['binding'], 2)
        self.assertEqual(steps, sorted(set(steps)))
        self.assertLess(steps[-1], result['steps'])  # Discarded branches are not movie time.
        for bank in result['banks'].values():
            segment = (r.output/bank['trajectory_file']).read_text().split('t =')[1:]
            native_steps = [int(part.splitlines()[0].strip()) for part in segment]
            measured = read_json(r.output/bank['observables_file'])
            self.assertEqual(native_steps, [point['step'] for point in measured])
            self.assertLessEqual(len(measured), 32)

    def test_mode_switch_cannot_reuse_the_other_checkpoint(self):
        r = self.replica()
        r.prepare_proximity()
        full = deepcopy(self.req)
        full['kinetics'].update(mode='opening', source_samples=75)
        self.assertNotEqual(request_identity(full), request_identity(self.req))
        with self.assertRaisesRegex(ValueError, 'different simulation'):
            FFSReplica(full, r.output, 0, adapter_factory=NativeScript([]))

    def test_master_progress_and_finalization_never_publish_rate_or_source_crossings(self):
        states = {}
        factory = self.replica().adapter_factory
        for index in range(4):
            req = deepcopy(self.req)
            r = FFSReplica(req, self.root/f'replica_{index+1:03d}', index,
                           adapter_factory=factory)
            r.prepare_proximity()
            r.state.update(completed=True, phase='completed', status='completed', stages=[stage(50)])
            r.save(force=True)
            (r.output/'result.json').write_text(json.dumps({'status': 'completed'}))
            states[index] = r.state
        sampler = KineticSampler(self.req, self.root)
        progress = sampler.aggregate(states, 0)
        self.assertNotIn('flux_target', progress['kinetics'])
        self.assertNotIn('flux_crossings', progress['kinetics'])
        self.assertEqual(progress['kinetics']['initial_sources'], 4)
        self.assertEqual(sampler.report(states, 'cancelled')['status'], 'incomplete')
        self.assertIsNone(sampler.report({0: states[0]}, 'completed')['value'])
        with patch('aptaswitch_core.oxdna_kinetics_worker.subprocess.Popen',
                   return_value=Mock(poll=Mock(return_value=0), returncode=0)), \
                patch('aptaswitch_core.oxdna_kinetics_worker.time.sleep'):
            result = sampler.run()
        self.assertEqual(result['simulation_mode'], 'proximity')
        self.assertIsNone(result['rate_estimate'])
        self.assertNotIn('rate', result['kinetics'])
        self.assertEqual(result['conditional_probability']['value'], .5)
        self.assertEqual(read_json(self.root/'result.json')['conditional_probability']['value'], .5)


if __name__ == '__main__':
    unittest.main()
