"""Shared functional families from one native partition per OFF/ON state."""
from dataclasses import replace
import json
import math
import os
from pathlib import Path
import random
import subprocess
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import numpy as np

from aptaswitch_core.structure_ensemble import EnsembleCancelled, cluster_structures, load_ensemble, save_ensemble, ensemble_key
from aptaswitch_core.thermodynamics.design_constraints import build_design_problem
from aptaswitch_core.thermodynamics.ensemble_analysis import analyze_structure_ensemble, native_ensemble_inputs
from aptaswitch_core.thermodynamics.native_backend import NativeContext, NativeEnsemble, NativeSamples, NativeHybridBackend
from aptaswitch_core.thermodynamics.types import Conditions
from test_thermodynamics_design_engine import AnalyticalBackend, request
from test_structure_ensemble import matrix_for


def candidate_for(req):
    assignment = build_design_problem(req).generate(random.Random(17))
    strands = {s.name: s for s in assignment.strands}
    return SimpleNamespace(switch_seq=strands['Switch'].sequence, trigger_seq=strands['Trigger'].sequence)


class EnsembleFixture(AnalyticalBackend):
    def ensemble(self, spec, conditions, count, *, seed):
        result = self.analyze(spec, conditions, need_pairs=True)
        samples = NativeSamples(NativeContext(spec, conditions, self.provenance, None, 1.),
                                (self.structure(spec),) * count, seed)
        return NativeEnsemble(result, samples)


class EnsembleAdapterTests(unittest.TestCase):
    def setUp(self):
        self.req = request()
        self.candidate = candidate_for(self.req)
        self.backend = EnsembleFixture(self.req)

    def test_model_seed_and_actual_chemistry_define_persistent_cache(self):
        data = analyze_structure_ensemble(self.backend, self.req, self.candidate, num_sample=100, seed=7)
        self.assertEqual(data['inputs']['sampling']['seeds'], {'off': 7, 'on': 8})
        self.assertEqual(data['inputs']['thermo']['material'], self.backend.provenance.model)
        self.assertEqual(data['inputs']['thermo']['trigger_material'], 'dna')
        self.assertEqual(data['model'], data['inputs']['thermodynamics'])
        for state in data['states'].values():
            self.assertEqual(cluster_structures(state)['families'][0]['count'], 100)
        for other in (native_ensemble_inputs(self.req, self.candidate, self.backend, num_sample=100, seed=8),
                      native_ensemble_inputs(replace(self.req, thermo=replace(self.req.thermo, sodium_m=.25)),
                                             self.candidate, self.backend, num_sample=100, seed=7)):
            self.assertNotEqual(ensemble_key(other), data['key'])
        with tempfile.TemporaryDirectory() as directory:
            path = save_ensemble(Path(directory), data)
            self.assertEqual(load_ensemble(Path(directory), data['inputs']), data)
            data['model']['version'] = 'wrong'
            # Detach inputs from model: mutation must not silently retag them.
            original_inputs = native_ensemble_inputs(self.req, self.candidate, self.backend, num_sample=100, seed=7)
            path.write_text(json.dumps({**data, 'inputs': original_inputs}))
            with self.assertRaisesRegex(ValueError, 'different thermodynamic'):
                load_ensemble(Path(directory), original_inputs)

    def test_cancel_before_and_between_states_never_returns_partial_success(self):
        with self.assertRaises(EnsembleCancelled):
            analyze_structure_ensemble(self.backend, self.req, self.candidate, cancel_check=lambda: True)
        self.assertEqual(self.backend.calls, [])
        with self.assertRaises(EnsembleCancelled):
            analyze_structure_ensemble(self.backend, self.req, self.candidate,
                                       cancel_check=lambda: bool(self.backend.calls))
        self.assertEqual(len(self.backend.calls), 1)

    def test_non_dna_trigger_and_dna_switch_are_rejected_before_sampling(self):
        for switch, trigger in (("dna", "rna"), ("rna", "rna"), ("dna", "dna")):
            req = replace(self.req, thermo=replace(self.req.thermo,
                          switch_material=switch, trigger_material=trigger))
            with self.subTest(switch=switch, trigger=trigger), self.assertRaises(ValueError):
                analyze_structure_ensemble(self.backend, req, self.candidate, num_sample=100)
        for candidate in (
            SimpleNamespace(switch_seq=self.candidate.switch_seq.replace("U", "T"), trigger_seq=self.candidate.trigger_seq),
            SimpleNamespace(switch_seq=self.candidate.switch_seq, trigger_seq=self.candidate.trigger_seq.replace("T", "U")),
        ):
            with self.subTest(candidate=candidate), self.assertRaises(ValueError):
                analyze_structure_ensemble(self.backend, self.req, candidate, num_sample=100)
        self.assertEqual(self.backend.calls, [])

    def test_mismatched_probability_sample_context_or_disconnected_samples_rejected(self):
        original = self.backend.ensemble
        def broken(spec, conditions, count, *, seed):
            record = original(spec, conditions, count, seed=seed)
            if len(spec.strands) == 2:
                samples = replace(record.samples, structures=('+'.join('.'*len(s.sequence) for s in spec.strands),)*count)
                return replace(record, samples=samples)
            return record
        with patch.object(self.backend, 'ensemble', side_effect=broken), self.assertRaisesRegex(ValueError, 'disconnected'):
            analyze_structure_ensemble(self.backend, self.req, self.candidate, num_sample=100)
        def wrong_context(spec, conditions, count, *, seed):
            record = original(spec, conditions, count, seed=seed)
            return replace(record, samples=replace(record.samples, context=replace(record.samples.context,
                                 conditions=Conditions(37., .25))))
        with patch.object(self.backend, 'ensemble', side_effect=wrong_context), self.assertRaises(ValueError):
            analyze_structure_ensemble(self.backend, self.req, self.candidate, num_sample=100)

    def test_invalid_seeds_and_sample_counts_fail_before_native_work(self):
        for seed in (True, -1, 2**64-1, .5):
            with self.assertRaises(ValueError):
                analyze_structure_ensemble(self.backend, self.req, self.candidate, seed=seed)
        for count in (0, 99, 5001, True):
            with self.assertRaises(ValueError):
                analyze_structure_ensemble(self.backend, self.req, self.candidate, num_sample=count)
        self.assertEqual(self.backend.calls, [])


@unittest.skipUnless(os.environ.get('APTASWITCH_HYBRID_NATIVE_LIBRARY'), 'Explicit native hybrid library required')
class RealEnsembleIntegrationTests(unittest.TestCase):
    def test_native_workflow_uses_no_python_reference_modules(self):
        program = '''
import importlib.abc, sys
sys.path.insert(0, sys.argv[1])
class RejectLegacy(importlib.abc.MetaPathFinder):
    def find_spec(self, fullname, path=None, target=None):
        if fullname.split('.')[0] in ('RNA', 'tests') or fullname in (
            'aptaswitch_core.tube_analysis',):
            raise AssertionError('Forbidden runtime dependency: ' + fullname)
sys.meta_path.insert(0, RejectLegacy())
from aptaswitch_core.saved_switch_analysis import analyze_saved_switch
from aptaswitch_core.saved_switch_ensemble import context_ensemble_source
from aptaswitch_core.structure_ensemble import analyze_structure_ensemble
from aptaswitch_core.thermodynamics.native_backend import NativeHybridBackend
from aptaswitch_core.thermodynamics.design_constraints import build_design_problem
from aptaswitch_core.models import RunRequest, SequenceInput
from aptaswitch_core.architecture import default_architecture_for_trigger
from types import SimpleNamespace
import random
from dataclasses import asdict
req = RunRequest(SequenceInput('AACCGGTTAACC', 'ACGT'*7+'AC'), default_architecture_for_trigger(30), trials=1)
strands = {s.name:s for s in build_design_problem(req).generate(random.Random(17)).strands}
candidate = SimpleNamespace(switch_seq=strands['Switch'].sequence, trigger_seq=strands['Trigger'].sequence)
backend = NativeHybridBackend(sys.argv[2])
result = analyze_structure_ensemble(req, candidate, backend=backend, num_sample=100)
assert result['inputs']['thermo']['trigger_material'] == 'dna'
assert sum(s['count'] for s in result['states']['on']['structures']) == 100
context = analyze_saved_switch({'request':asdict(req), 'candidate':vars(candidate)}, None, {},
                               backend=backend, cancel_check=lambda: True)
assert context['status'] == 'cancelled'
assert context['thermodynamics']['backend'] == backend.provenance.backend
print('Independent native ensemble workflow OK')
'''
        root = Path(__file__).resolve().parents[1]
        completed = subprocess.run([sys.executable, '-I', '-c', program, str(root/'src'),
                                    os.environ['APTASWITCH_HYBRID_NATIVE_LIBRARY']],
                                   capture_output=True, text=True, timeout=45)
        self.assertEqual(completed.returncode, 0, completed.stderr)
        self.assertIn('workflow OK', completed.stdout)

    def test_one_native_handle_per_state_joint_probabilities_and_reproducibility(self):
        req = request(); candidate = candidate_for(req)
        backend = NativeHybridBackend(os.environ['APTASWITCH_HYBRID_NATIVE_LIBRARY'])
        with patch.object(backend, '_request', wraps=backend._request) as native:
            result = analyze_structure_ensemble(backend, req, candidate, num_sample=1000, seed=29)
            self.assertEqual(native.call_count, 2)
        repeated = analyze_structure_ensemble(backend, req, candidate, num_sample=1000, seed=29)
        self.assertEqual(result['states'], repeated['states'])
        for name, state in result['states'].items():
            families = cluster_structures(state)
            self.assertEqual(sum(f['count'] for f in families['families']), 1000)
            self.assertAlmostEqual(sum(f['fraction'] for f in families['families']), 1.)
            empirical = sum(row['count'] * matrix_for(row['structure']) for row in state['structures']) / 1000
            n = sum(state['strand_lengths'])
            # Union bound for every pair/unpaired Bernoulli marginal, alpha .01.
            bound = math.sqrt(math.log(2*n*n/.01) / (2*1000))
            self.assertLess(np.max(np.abs(empirical - np.asarray(state['pairs']))), bound)
            self.assertEqual(state['switch_offset'], 30 if name == 'on' else 0)


if __name__ == '__main__':
    unittest.main()
