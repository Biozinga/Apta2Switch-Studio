"""Fixed full-reporter orchestration and origin-label conservation."""
from copy import deepcopy
from dataclasses import asdict
import math
import os
from pathlib import Path
import unittest

from aptaswitch_core.saved_switch_ensemble import context_ensemble_source
from aptaswitch_core.thermodynamics.context_analysis import analyze_saved_switch
from aptaswitch_core.thermodynamics.native_backend import NativeContext, NativeMinimumStructures, NativeHybridBackend
from aptaswitch_core.thermodynamics.tubes import StandardStateConvention
from aptaswitch_core.thermodynamics.types import ModelProvenance, PartitionResult
from test_saved_switch_ensemble import context_fixture


class FixedWeightBackend:
    """Analytical integration fixture: monomers/dimers only, not a physical model."""
    provenance = ModelProvenance('fixture', '1', 'fixture', 'fixture:1', 'connected', '1M')
    standard_state = StandardStateConvention('1M', 1., True, True, True)
    condition_limits = dict(temperature_min_c=25., temperature_max_c=55., sodium_min_m=.12,
                            sodium_max_m=1., magnesium_min_m=0., magnesium_max_m=0.)

    def __init__(self):
        self.calls = []
        self.fail_temperature = None
        self.empty_on = False

    def partition(self, spec, conditions):
        self.calls.append((spec.cyclic_key, conditions))
        if conditions.celsius == self.fail_temperature:
            raise ArithmeticError('Test failure at a later point')
        size = len(spec.strands)
        log_z = 0. if size == 1 else math.log(1e6) if size == 2 else -math.inf
        if self.empty_on and size == 2 and {s.name for s in spec.strands} == {'Switch', 'Trigger'}:
            log_z = -math.inf
        return PartitionResult(spec, conditions, self.provenance, log_z, -log_z, 1.)

    def minimum_structures(self, spec, conditions, *, maximum_count):
        parts = ['.' * len(s.sequence) for s in spec.strands]
        if len(parts) == 2:
            parts[0] = '(' + parts[0][1:]
            parts[1] = parts[1][:-1] + ')'
        structure = '+'.join(parts)
        return NativeMinimumStructures(NativeContext(spec, conditions, self.provenance, None, 1.),
                                       0., (structure,), (0.,), False, maximum_count, 1e-8)


def snapshot_fixture():
    _, snapshot, _ = context_fixture()
    return snapshot


class NativeContextOrchestrationTests(unittest.TestCase):
    def test_complete_reporter_three_strand_catalog_scans_and_mass_balance(self):
        from aptaswitch_core.saved_switch_analysis import analyze_saved_switch as analyze_public
        snapshot, backend = snapshot_fixture(), FixedWeightBackend()
        before = deepcopy(snapshot)
        result = analyze_public(snapshot, None, dict(temperature_c=37., sodium_m=.12,
            include_aptamer=True, temperature_values_c=[37., 29.], sodium_values_m=[.12, .25]), backend=backend)
        self.assertEqual(snapshot, before)
        self.assertEqual(result['status'], 'completed', result.get('error'))
        self.assertEqual(result['context']['switch_length'], 797)
        self.assertTrue(result['context']['switch_sequence'].startswith(snapshot['candidate']['switch_seq']))
        self.assertEqual(result['thermodynamics'], asdict(backend.provenance))
        baseline = result['baseline']

        self.assertEqual(result['config']['max_complex_size'], 3)
        self.assertEqual(baseline['on']['population']['equilibrium']['ordered_physical_complex_count'], 20)
        self.assertTrue(all(len(row['strand_order']) <= 3
                            for row in baseline['on']['population']['complexes']))
        self.assertEqual(len(result['sweeps']['temperature_c']), 2)
        self.assertEqual(len(result['sweeps']['sodium_m']), 2)
        self.assertEqual(len(backend.calls), 20 * 3)
        self.assertEqual(len(set(backend.calls)), len(backend.calls))
        for point in [baseline, *result['sweeps']['temperature_c'], *result['sweeps']['sodium_m']]:
            for state in ('off', 'on'):
                population = point[state]['population']
                for name, total in population['total_strand_concentrations_m'].items():
                    actual = math.fsum(row['strand_order'].count(name) * row['concentration_m']
                                       for row in population['complexes'])
                    self.assertAlmostEqual(actual / total, 1., places=8)
        empty = [row for row in baseline['on']['population']['complexes'] if len(row['strand_order']) > 2]
        self.assertTrue(all(row['mfe'] is None and row['structure_status'] == 'empty_ensemble' for row in empty))
        source, candidate = context_ensemble_source(snapshot, result, Path('/trusted/context'), backend=backend)
        self.assertEqual(len(candidate.switch_seq), 797)
        self.assertEqual(source.request.thermo.material, backend.provenance.model)
        self.assertIsNone(source.request.thermodynamics_library_path)

    def test_identical_aptamer_trigger_merge_without_losing_origin_totals(self):
        snapshot, backend = snapshot_fixture(), FixedWeightBackend()
        snapshot['request']['sequences']['aptamer_dna'] = snapshot['candidate']['trigger_seq']
        result = analyze_saved_switch(backend, snapshot, dict(include_aptamer=True,
            switch_concentration_m=1e-6, trigger_concentration_m=2e-6, aptamer_concentration_m=6e-6))
        self.assertEqual(result['status'], 'completed', result.get('error'))
        population = result['baseline']['on']['population']
        self.assertTrue(population['equilibrium']['chemical_aliases_merged'])
        self.assertEqual(population['equilibrium']['ordered_physical_complex_count'], 9)
        rows = {tuple(sorted(row['strand_order'])): row['concentration_m'] for row in population['complexes']}
        self.assertAlmostEqual(rows[('Aptamer', 'Switch')] / rows[('Switch', 'Trigger')], 3., places=8)
        self.assertAlmostEqual(rows[('Aptamer', 'Trigger')] / rows[('Trigger', 'Trigger')], 6., places=8)
        for name, total in population['total_strand_concentrations_m'].items():
            actual = math.fsum(row['strand_order'].count(name) * row['concentration_m'] for row in population['complexes'])
            self.assertAlmostEqual(actual / total, 1., places=8)

    def test_domain_checked_before_call_and_completed_reference_survives_failure(self):
        snapshot, backend = snapshot_fixture(), FixedWeightBackend()
        for bad in ({'temperature_values_c': [60]}, {'sodium_values_m': [.1]}, {'magnesium_m': .001}):
            with self.assertRaises(ValueError):
                analyze_saved_switch(backend, snapshot, bad)
        self.assertEqual(backend.calls, [])
        backend.fail_temperature = 40.
        result = analyze_saved_switch(backend, snapshot, dict(temperature_c=37., temperature_values_c=[40.]))
        self.assertEqual(result['status'], 'failed')
        self.assertIsNotNone(result['baseline'])
        self.assertEqual(result['sweeps']['temperature_c'], [])
        self.assertIn('later point', result['error'])
        context_ensemble_source(snapshot, result, Path('/tmp/context'), backend=backend)

    def test_cancel_does_not_publish_partial_point(self):
        backend = FixedWeightBackend()
        result = analyze_saved_switch(backend, snapshot_fixture(), {}, cancel_check=lambda: len(backend.calls) >= 3)
        self.assertEqual(result['status'], 'cancelled')
        self.assertIsNone(result['baseline'])
        self.assertEqual(len(backend.calls), 3)

    def test_empty_connected_on_is_reported_without_fabricating_a_structure(self):
        backend = FixedWeightBackend(); backend.empty_on = True
        result = analyze_saved_switch(backend, snapshot_fixture(), {})
        self.assertEqual(result['status'], 'completed')
        self.assertIsNone(result['baseline']['on']['mfe'])
        self.assertEqual(result['baseline']['metrics']['on_yield_pct'], 0.)
        with self.assertRaises(ValueError):
            context_ensemble_source(snapshot_fixture(), result, Path('/tmp/context'), backend=backend)

    def test_subfamilies_reject_different_backend_model(self):
        backend = FixedWeightBackend()
        result = analyze_saved_switch(backend, snapshot_fixture(), {})
        result['baseline']['thermodynamics'] = {**result['thermodynamics'], 'version': 'another'}
        with self.assertRaisesRegex(ValueError, 'modèle exact'):
            context_ensemble_source(snapshot_fixture(), result, Path('/tmp/context'), backend=backend)


@unittest.skipUnless(os.environ.get('APTASWITCH_HYBRID_NATIVE_LIBRARY'), 'Explicit native hybrid library required')
class RealContextIntegrationTests(unittest.TestCase):
    def test_full_short_reporter_and_sodium_scan(self):
        snapshot = snapshot_fixture()
        reporter = snapshot['request']['reporter']
        reporter['full_sequence_after_start_rna'] = reporter['sequence_after_start_rna'] + 'GCU' * 3
        backend = NativeHybridBackend(os.environ['APTASWITCH_HYBRID_NATIVE_LIBRARY'])
        result = analyze_saved_switch(backend, snapshot, dict(max_complex_size=2, temperature_c=37.,
                                                            sodium_m=.12, sodium_values_m=[.25]))
        self.assertEqual(result['status'], 'completed', result.get('error'))
        self.assertEqual(result['context']['switch_length'], 122)
        for point in (result['baseline'], result['sweeps']['sodium_m'][0]):
            self.assertTrue(point['on']['population']['mass_balance_checked'])
            self.assertLess(point['on']['population']['equilibrium']['max_relative_mass_balance_error'], 1e-8)
            self.assertTrue(0 <= point['metrics']['on_yield_pct'] <= 100)
            self.assertEqual(list(map(len, point['on']['mfe']['structure'].split('+'))), [30, 122])
        self.assertNotEqual(result['baseline']['off']['mfe']['energy_kcal_mol'],
                            result['sweeps']['sodium_m'][0]['off']['mfe']['energy_kcal_mol'])


if __name__ == '__main__':
    unittest.main()
