"""Candidate-scoped kinetic requests and dispatch preserve production science."""
from contextlib import chdir
from copy import deepcopy
import json
import os
from pathlib import Path
import sys
import tempfile
import time
from types import SimpleNamespace
import unittest
from unittest.mock import Mock, patch

from aptaswitch_core.architecture import default_architecture_for_trigger
from aptaswitch_core.models import RunRequest, SequenceInput
from aptaswitch_core.oxdna_ffs import (
    _finish_result, load_ffs_resume_request, prepare_ffs_request,
)
from aptaswitch_core import oxdna_ffs_worker
from aptaswitch_core.oxdna_kinetics import candidate_interfaces
from aptaswitch_studio.saved_switch_library import candidate_snapshot
from tests.design_fixtures import design_result_fixture


class KineticsIntegrationTests(unittest.TestCase):
    def setUp(self):
        self.root = Path(self.enterContext(tempfile.TemporaryDirectory()))
        request = RunRequest(sequences=SequenceInput('ACGT'*3, 'ACGT'*7+'AC'),
                             architecture=default_architecture_for_trigger(30), trials=1)
        result = design_result_fixture(request)
        self.snapshot = candidate_snapshot(result, result.candidates[0], 'Kinetic candidate')

    def request(self, options=None):
        return prepare_ffs_request(self.snapshot, {}, options or {})

    def test_kinetic_defaults_are_frozen_without_adaptive_diversity(self):
        before = deepcopy(self.snapshot)
        request = self.request()
        self.assertEqual(request['kinetics'], {
            'version':3, 'mode':'opening', 'source_samples':75, 'trials':250, 'replicas':10,
            'workers':min(10, os.cpu_count() or 1), 'equilibration_steps':200000, 'chunk_steps':1000})
        self.assertNotIn('diversity',request['sampling'])
        self.assertNotIn('budget_seconds',request['sampling'])
        self.assertNotIn('max_steps_per_trial',request['sampling'])
        self.assertTrue(request['plan']['rate_estimation'])
        self.assertEqual(request['purpose'],'native_direct_ffs_flux_and_shooting')
        self.assertEqual(self.snapshot,before)

    def test_prefixed_kinetic_controls_do_not_reuse_exploratory_quotas(self):
        request = self.request({'kinetics_source_samples':150, 'kinetics_trials':500,
                                'kinetics_replicas':6, 'kinetics_workers':2})
        self.assertEqual(request['kinetics']['source_samples'],150)
        self.assertEqual(request['kinetics']['trials'],500)
        self.assertNotIn('trials_per_interface',request['sampling'])
        self.assertEqual(request['kinetics']['replicas'],6)
        self.assertEqual(request['kinetics']['workers'],2)

    def test_single_pipeline_includes_approach_full_duplex_and_opening(self):
        request = self.request({'seed': 123})
        interfaces = candidate_interfaces(request['plan'])
        self.assertEqual(interfaces[0]['name'], 'contact_8')
        self.assertTrue(any(item['name'] == 'duplex_30' for item in interfaces))
        self.assertEqual(interfaces[-1], {'name': 'opening', 'min_binding': 30,
                                         'max_stem': 0, 'max_opening': 0})
        self.assertNotIn('source_run_dir', request)
        for option in ({'simulation_mode': 'encounter'}, {'simulation_mode': 'opening'},
                       {'source_run_dir': str(self.root/'prior-contact-run')}):
            with self.subTest(option=option), self.assertRaises(ValueError):
                self.request(option)

    def test_ineligible_strand_chemistry_and_length_rejected(self):
        for changes in ({'switch_material':'dna'}, {'trigger_material':'rna'}):
            snapshot = deepcopy(self.snapshot)
            snapshot['request']['thermo'].update(changes)
            with self.subTest(changes=changes),self.assertRaises(ValueError):
                prepare_ffs_request(snapshot,{}, {})
        snapshot = deepcopy(self.snapshot)
        snapshot['candidate']['trigger_seq'] = snapshot['candidate']['trigger_seq'][:-1]
        with self.assertRaisesRegex(ValueError,'30'):
            prepare_ffs_request(snapshot,{}, {})

    def test_kinetics_rejects_time_censoring_and_bad_sampling(self):
        for options in ({'budget_minutes':1}, {'max_steps':10000}, {'budget_seconds':5},
                        {'max_steps_per_trial':1000}, {'kinetics':'true'},
                        {'kinetics_source_samples':74}, {'kinetics_source_samples':201},
                        {'kinetics_trials':249}, {'kinetics_trials':1001},
                        {'kinetics_replicas':3}, {'kinetics_replicas':101},
                        {'kinetics_workers':0}, {'kinetics_workers':65}, {'kinetics_pilot':100},
                        {'kinetics_replicas':True}, {'kinetics_trials':500.5}):
            with self.subTest(options=options),self.assertRaises(ValueError):
                self.request(options)

    def test_production_sampling_accepts_both_inclusive_boundaries(self):
        for source, trials in ((75, 250), (200, 1000)):
            with self.subTest(source=source, trials=trials):
                request = self.request({'kinetics_source_samples':source, 'kinetics_trials':trials})
                self.assertEqual(request['kinetics']['source_samples'], source)
                self.assertEqual(request['kinetics']['trials'], trials)


    def test_finalization_preserves_native_rates_interfaces_and_geometry(self):
        request = self.request()
        result = {'status':'completed', 'method':'NATIVE_FFS',
                  'kinetics':{'phase':'report','status':'model_estimate','rate':{'value_model':.123}},
                  'plan':{'interfaces':[{'name':'contact_7','distance_max':7.}]},
                  'box_size':200., 'box_policy':{'sizing_method':'finite_periodic_box'},
                  'concentration_each_m':None, 'sampling':{'production_trials':1000},
                  'rate_estimate':{'status':'model_estimate','value_model':.123,'interval95_model':[.1,.2]},
                  'trajectory_file':'movie.dat'}
        expected = deepcopy(result)
        finished = _finish_result(result,request,self.root,time.monotonic(),None)
        for key in ('rate_estimate','method','kinetics','plan','box_size','box_policy','concentration_each_m','sampling'):
            self.assertEqual(finished[key],expected[key],key)
        self.assertEqual(finished['movie_file'],'movie.dat')
        self.assertEqual(json.loads((self.root/'result.json').read_text())['rate_estimate'],expected['rate_estimate'])


    def test_resume_retains_frozen_kinetic_allocation(self):
        request = self.request()
        # A recorded run is independent of the current new-run control ranges.
        request['kinetics'].update(source_samples=300, trials=2000)
        preparation = self.root/'preparation'
        preparation.mkdir()
        for name in ('validated.dat','system.top'):
            (preparation/name).write_text('test fixture')
        request['prepared'] = {'configuration_file':str(preparation/'validated.dat'),
                               'topology_file':str(preparation/'system.top')}
        (self.root/'worker_request.json').write_text(json.dumps(request))
        resumed = load_ffs_resume_request(self.root)
        self.assertTrue(resumed['resume'])
        self.assertNotIn('diversity',resumed['sampling'])
        self.assertNotIn('upgrade_diversity',resumed)
        self.assertEqual(resumed['kinetics'],request['kinetics'])

    def test_worker_dispatches_kinetics_without_starting_exploratory_sampler(self):
        request = self.request()
        path = self.root/'request.json'
        path.write_text(json.dumps(request))
        kinetic_class = Mock(return_value=Mock(run=Mock(return_value={'status':'completed'})))
        with chdir(self.root), patch.object(sys,'argv',['worker',str(path)]), \
                patch.dict(sys.modules,{'aptaswitch_core.oxdna_kinetics_worker':SimpleNamespace(KineticSampler=kinetic_class)}):
            self.assertEqual(oxdna_ffs_worker.main(),0)
        self.assertFalse(hasattr(oxdna_ffs_worker,'Sampler'))
        self.assertEqual(kinetic_class.call_args.args[0]['kinetics'],request['kinetics'])
        self.assertEqual(kinetic_class.call_args.args[1],self.root.resolve())


if __name__ == '__main__':
    unittest.main()
