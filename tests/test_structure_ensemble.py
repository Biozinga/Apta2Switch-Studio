"""Ensemble families must preserve co-occurrence, counts and strand order."""
from __future__ import annotations

from contextlib import contextmanager
import copy
from dataclasses import asdict, replace
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import numpy as np

from aptaswitch_core.architecture import default_architecture_for_trigger, target_switch_structures
from aptaswitch_core.defect_metrics import _pair_targets
from aptaswitch_core.models import RunRequest, SequenceInput, ThermoConditions
from aptaswitch_core.structure_ensemble import (
    EnsembleCancelled, analyze_structure_ensemble, assemble_ensemble, cluster_structures, ensemble_inputs,
    ensemble_key, load_ensemble, probability_helices, save_ensemble,
)
from tests.design_fixtures import design_result_fixture
from aptaswitch_core.thermodynamics.types import ModelProvenance


def fixture(output_dir=None):
    request = RunRequest(sequences=SequenceInput('ACGTACGT', 'ACGT' * 7 + 'AC'),
                         architecture=default_architecture_for_trigger(30), trials=1,
                         thermo=ThermoConditions(temperature_c=31, sodium_m=.15))
    if output_dir is not None:
        request = replace(request, output_dir=Path(output_dir))
    return design_result_fixture(request)


def matrix_for(structure):
    lengths = list(map(len, structure.split('+')))
    partners = _pair_targets(structure, lengths)
    matrix = np.zeros((len(partners), len(partners)))
    matrix[np.arange(len(partners)), partners] = 1
    return matrix


def fake_analysis(result, num_sample=1000, *, backend=None, seed=0):
    candidate = result.candidates[0]
    off, on = target_switch_structures(result.request.architecture, len(result.request.reporter.sequence_after_start_rna))
    return assemble_ensemble(
        ensemble_inputs(result.request, candidate, num_sample, backend=backend, seed=seed),
        {"off": matrix_for(off), "on": matrix_for(on)},
        {"off": {off: num_sample}, "on": {on: num_sample}},
        asdict(backend.provenance) if backend is not None else {"engine": "saved-fixture", "version": "test"},
    )


def simple_state(structures, *, regions=None):
    lengths = list(map(len, structures[0][0].split('+')))
    return dict(strand_lengths=lengths, switch_offset=lengths[0] if len(lengths) == 2 else 0,
                target=structures[0][0], regions=regions or {'off_stem': [0]},
                structures=[{'structure': s, 'count': count} for s, count in structures])


class FamilyTests(unittest.TestCase):
    def test_differences_outside_functional_regions_share_a_family(self):
        state = simple_state([('(....)..', 60), ('(....)()', 30), ('.(....).', 10)])
        result = cluster_structures(state, tolerance=0)
        self.assertEqual([f['count'] for f in result['families']], [90, 10])
        self.assertEqual(result['families'][0]['representative'], '(....)..')
        self.assertEqual(result['families'][0]['distinct_structures'], 2)
        self.assertEqual(sum(f['fraction'] for f in result['families']), 1)
        self.assertEqual(result['families'][0]['regional_defects']['off_stem'], 0)
        self.assertEqual(result['families'][1]['regional_defects']['off_stem'], 1)

    def test_partner_identity_matters_not_only_paired_or_unpaired(self):
        state = simple_state([('(....)..', 50), ('(.....).', 50)])
        self.assertEqual(len(cluster_structures(state, tolerance=0)['families']), 2)

    def test_on_indices_refer_to_switch_and_actual_trigger_partner(self):
        state = simple_state([('(.+.)..', 65), ('.(+.)..', 35)], regions={'on_binding': [1]})
        result = cluster_structures(state, tolerance=0)
        self.assertEqual(result['selected_count'], 1)
        self.assertEqual(len(result['families']), 2)
        self.assertEqual([f['regional_defects']['on_binding'] for f in result['families']], [0, 1])

    def test_distance_is_bounded_to_seed_and_counts_are_not_lost(self):
        state = simple_state([('.'*20, 60), ('(....)' + '.'*14, 40)], regions={'off_stem': list(range(20))})
        self.assertEqual(len(cluster_structures(state, tolerance=.05)['families']), 2)
        grouped = cluster_structures(state, tolerance=.1)
        self.assertEqual(grouped['max_changes'], 2)
        self.assertEqual(len(grouped['families']), 1)
        self.assertEqual(grouped['families'][0]['count'], 100)
        self.assertAlmostEqual(grouped['families'][0]['regional_defects']['off_stem'], .04)
        state['structures'].reverse()
        self.assertEqual(grouped, cluster_structures(state, tolerance=.1))

    def test_never_constructs_consensus_by_combining_incompatible_pairs(self):
        # Same base has two alternative partners. Every displayed structure
        # must still be one of the observed complete samples.
        state = simple_state([('(....)..', 60), ('(.....).', 40)], regions={'off_stem': list(range(8))})
        for family in cluster_structures(state, tolerance=.2)['families']:
            self.assertIn(family['representative'], {row['structure'] for row in state['structures']})


class MotifTests(unittest.TestCase):
    def test_reports_competitors_and_marginal_not_joint_probabilities(self):
        matrix = .6 * matrix_for('((....))..') + .4 * matrix_for('.((....)).')
        helices = probability_helices({'pairs': matrix.tolist(), 'strand_lengths': [10]}, threshold=.3)
        self.assertEqual([h['length'] for h in helices], [2, 2])
        self.assertEqual([h['competing_motifs'] for h in helices], [[2], [1]])
        self.assertAlmostEqual(helices[0]['mean_pair_probability'], .6)
        self.assertEqual(len(probability_helices({'pairs': matrix, 'strand_lengths': [10]}, threshold=.7)), 0)

    def test_stems_do_not_cross_strand_breaks(self):
        state = {'pairs': matrix_for('((+....))'), 'strand_lengths': [2, 6]}
        self.assertEqual([h['length'] for h in probability_helices(state)], [2])
        state = {'pairs': matrix_for('(+(....))'), 'strand_lengths': [1, 7]}
        self.assertEqual([h['length'] for h in probability_helices(state)], [1, 1])


class AnalysisTests(unittest.TestCase):
    @staticmethod
    def native_analysis(*, chemistry_policy=False):
        conventions = "connected test states"
        if chemistry_policy:
            conventions += ";RNA-WC-GU;DNA-WC;hybrid-WC"
        backend = SimpleNamespace(provenance=ModelProvenance(
            "finite-test", "1", "fixture", "fixture:ensemble", conventions, "fixture-1M"))
        analysis = fake_analysis(fixture(), 100, backend=backend, seed=7)
        if chemistry_policy:
            # Targets are conceptual diagrams, not a promise about the
            # analytical fixture's actual sequence. Supply valid samples.
            AnalysisTests.set_sample_pairs(analysis, 'off', ())
            trigger, switch = analysis['inputs']['trigger'], analysis['inputs']['switch']
            AnalysisTests.set_sample_pairs(analysis, 'on', ((trigger.index('A'), len(trigger) + switch.index('U')),))
        return analysis

    @staticmethod
    def set_sample_pairs(analysis, state_name, pairs):
        state = analysis['states'][state_name]
        lengths = state['strand_lengths']
        flat = ['.'] * sum(lengths)
        for i, j in pairs:
            flat[i], flat[j] = '(', ')'
        chunks, offset = [], 0
        for length in lengths:
            chunks.append(''.join(flat[offset:offset + length]))
            offset += length
        state['structures'] = [{'structure': '+'.join(chunks), 'count': analysis['inputs']['num_sample']}]

    def test_saved_conditions_chemistry_and_full_sample_counts(self):
        result = fixture()
        data = fake_analysis(result, 100)
        self.assertEqual(data['inputs']['thermo']['temperature_c'], 31)
        self.assertEqual(data['states']['on']['switch_offset'], 30)
        self.assertEqual(cluster_structures(data['states']['on'])['families'][0]['count'], 100)
        self.assertEqual(data['model']['version'], 'test')

    def test_save_reload_portable_and_cannot_reuse_different_candidate(self):
        result = fixture()
        analysis = fake_analysis(result)
        with tempfile.TemporaryDirectory() as directory:
            path = save_ensemble(Path(directory), analysis)
            self.assertEqual(load_ensemble(Path(directory), analysis['inputs']), analysis)
            other = ensemble_inputs(result.request, replace(result.candidates[0], switch_seq='C'*len(result.candidates[0].switch_seq)))
            self.assertIsNone(load_ensemble(Path(directory), other))
            corrupt = copy.deepcopy(analysis)
            corrupt['states']['on']['regions']['on_binding'] = [0]
            path.write_text(json.dumps(corrupt))
            with self.assertRaises(ValueError):
                load_ensemble(Path(directory), analysis['inputs'])
        self.assertNotEqual(ensemble_key(analysis['inputs']), ensemble_key(ensemble_inputs(result.request, result.candidates[0], 3000)))

    def test_native_saved_samples_reject_bad_syntax_counts_and_disconnected_complexes(self):
        analysis = self.native_analysis()
        lengths = analysis['states']['on']['strand_lengths']
        disconnected = '+'.join('.' * length for length in lengths)
        # Intramolecular pairs in both strands still leave the complex graph
        # disconnected: the validation must require an interstrand connection.
        separately_folded = '+'.join('(' + '.' * (length - 2) + ')' for length in lengths)
        malformed = analysis['states']['on']['structures'][0]['structure'].replace('(', ')', 1)
        for label, structures, error in (
            ('syntax', [{'structure': malformed, 'count': 100}], 'unbalanced'),
            ('sample total', [{'structure': disconnected, 'count': 99}], 'Incomplete saved sample'),
            ('unpaired strands', [{'structure': disconnected, 'count': 100}], 'disconnected complex'),
            ('separately folded strands', [{'structure': separately_folded, 'count': 100}], 'disconnected complex'),
        ):
            with self.subTest(label=label), tempfile.TemporaryDirectory() as directory:
                corrupt = copy.deepcopy(analysis)
                corrupt['states']['on']['structures'] = structures
                save_ensemble(Path(directory), corrupt)
                with self.assertRaisesRegex(ValueError, error):
                    load_ensemble(Path(directory), analysis['inputs'])

    def test_native_saved_sample_counts_must_be_positive_integers(self):
        analysis = self.native_analysis()
        structure = analysis['states']['on']['structures'][0]['structure']
        for counts in ((True, 99), (-1, 101), (1.0, 99)):
            with self.subTest(counts=counts), tempfile.TemporaryDirectory() as directory:
                corrupt = copy.deepcopy(analysis)
                corrupt['states']['on']['structures'] = [
                    {'structure': structure, 'count': count} for count in counts]
                save_ensemble(Path(directory), corrupt)
                with self.assertRaisesRegex(ValueError, 'Invalid saved sample count'):
                    load_ensemble(Path(directory), analysis['inputs'])

    def test_native_cache_roundtrip_is_portable_without_runtime(self):
        analysis = self.native_analysis(chemistry_policy=True)
        with tempfile.TemporaryDirectory() as directory:
            original = Path(directory) / 'original'
            path = save_ensemble(original, analysis)
            self.assertEqual(load_ensemble(original, analysis['inputs']), analysis)
            moved = original.rename(Path(directory) / 'moved')
            self.assertFalse(original.exists())
            program = '''
import importlib.abc, json, pathlib, sys
sys.path.insert(0, sys.argv[1])
class RejectRuntime(importlib.abc.MetaPathFinder):
    def find_spec(self, fullname, path=None, target=None):
        if fullname.split('.')[0] in ('RNA', 'tests') or fullname in (
            'aptaswitch_core.thermodynamics.setup',
            'aptaswitch_core.thermodynamics.native_backend',
            'aptaswitch_core.thermodynamics.native_runtime'):
            raise AssertionError('Unexpected runtime dependency: ' + fullname)
sys.meta_path.insert(0, RejectRuntime())
from aptaswitch_core.structure_ensemble import load_ensemble
root = pathlib.Path(sys.argv[2])
analysis = json.loads((root / sys.argv[3]).read_text())
assert load_ensemble(root, analysis['inputs']) == analysis
'''
            completed = subprocess.run(
                [sys.executable, '-I', '-c', program, str(Path(__file__).resolve().parents[1] / 'src'),
                 str(moved), str(path.relative_to(original))], capture_output=True, text=True, timeout=20)
            self.assertEqual(completed.returncode, 0, completed.stderr)

    def test_legacy_cache_without_connected_sampling_marker_keeps_previous_reader(self):
        analysis = fake_analysis(fixture(), 100)
        self.assertNotIn('sampling', analysis['inputs'])
        trigger, switch = analysis['inputs']['trigger'], analysis['inputs']['switch']
        for pairs in ((), ((trigger.index('A'), len(trigger) + switch.index('A')),)):
            with self.subTest(pairs=pairs), tempfile.TemporaryDirectory() as directory:
                self.set_sample_pairs(analysis, 'on', pairs)
                save_ensemble(Path(directory), analysis)
                self.assertEqual(load_ensemble(Path(directory), analysis['inputs']), analysis)

    def test_declared_native_chemistry_accepts_hybrid_wc_dna_wc_and_rna_gu(self):
        analysis = self.native_analysis(chemistry_policy=True)
        trigger, switch = analysis['inputs']['trigger'], analysis['inputs']['switch']
        self.set_sample_pairs(analysis, 'off', ((switch.index('G'), switch.index('U')),))
        self.set_sample_pairs(analysis, 'on', (
            (trigger.index('A'), trigger.index('T', 4)),
            (trigger.rindex('A'), len(trigger) + switch.index('U')),
        ))
        with tempfile.TemporaryDirectory() as directory:
            save_ensemble(Path(directory), analysis)
            self.assertEqual(load_ensemble(Path(directory), analysis['inputs']), analysis)

    def test_declared_native_chemistry_rejects_aa_dna_gt_and_hybrid_gu(self):
        analysis = self.native_analysis(chemistry_policy=True)
        trigger, switch = analysis['inputs']['trigger'], analysis['inputs']['switch']
        connection = (trigger.rindex('A'), len(trigger) + switch.index('U'))
        for label, pairs in (
            ('A-A', ((trigger.index('A'), len(trigger) + switch.index('A')),)),
            ('DNA G-T', ((trigger.index('G'), trigger.index('T', 4)), connection)),
            ('hybrid G-U', ((trigger.index('G'), len(trigger) + switch.index('U')),)),
        ):
            with self.subTest(label=label), tempfile.TemporaryDirectory() as directory:
                self.set_sample_pairs(analysis, 'on', pairs)
                save_ensemble(Path(directory), analysis)
                with self.assertRaisesRegex(ValueError, 'chemically forbidden base pair'):
                    load_ensemble(Path(directory), analysis['inputs'])

    def test_native_model_without_declared_pair_policy_is_not_restricted(self):
        analysis = self.native_analysis()
        trigger, switch = analysis['inputs']['trigger'], analysis['inputs']['switch']
        self.set_sample_pairs(analysis, 'on', ((trigger.index('A'), len(trigger) + switch.index('A')),))
        with tempfile.TemporaryDirectory() as directory:
            save_ensemble(Path(directory), analysis)
            self.assertEqual(load_ensemble(Path(directory), analysis['inputs']), analysis)

    def test_cancellation_before_native_calls(self):
        result = fixture()
        from test_thermodynamics_ensemble_analysis import EnsembleFixture
        backend = EnsembleFixture(result.request)
        with patch.object(backend, "ensemble", wraps=backend.ensemble) as native:
            with self.assertRaises(EnsembleCancelled):
                analyze_structure_ensemble(result.request, result.candidates[0], backend=backend,
                                           cancel_check=lambda: True)
        native.assert_not_called()

    @unittest.skipUnless(os.environ.get('APTASWITCH_HYBRID_NATIVE_LIBRARY'), 'Native hybrid path not configured')
    def test_real_native_dna_trigger_rna_switch_order(self):
        result = fixture()
        data = analyze_structure_ensemble(result.request, result.candidates[0], None,
            thermodynamics_library_path=os.environ['APTASWITCH_HYBRID_NATIVE_LIBRARY'], num_sample=1000)
        self.assertIn('parameter_digest', data['model'])
        self.assertEqual(data['inputs']['thermo']['trigger_material'], 'dna')
        self.assertEqual(data['inputs']['thermo']['switch_material'], 'rna')
        for name, state in data['states'].items():
            with self.subTest(state=name):
                grouped = cluster_structures(state)
                self.assertEqual(sum(f['count'] for f in grouped['families']), 1000)
                self.assertAlmostEqual(sum(f['fraction'] for f in grouped['families']), 1)
                self.assertTrue(np.allclose(np.array(state['pairs']).sum(axis=1), 1))
                self.assertEqual(list(map(len, state['structures'][0]['structure'].split('+'))), state['strand_lengths'])
                empirical = sum(row['count'] * matrix_for(row['structure']) for row in state['structures']) / 1000
                # Matrices and samples must use the requested DNA trigger/RNA switch order.
                self.assertLess(float(np.max(np.abs(empirical - np.asarray(state['pairs'])))), .1)
