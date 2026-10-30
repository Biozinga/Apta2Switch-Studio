"""Native direct FFS interfaces, source provenance and whole-replica rates."""
from copy import deepcopy
import math
import statistics
import unittest

from aptaswitch_core.oxdna_kinetics import (
    OPENING_PROTOCOL_VERSION, candidate_interfaces, condition, direct_ffs_report,
    fixed_trial_probability, progress_index, source_fingerprint, validate_kinetic_protocol,
)


def stage(successes, target=100, *, failures=None, censored=0, pending=0):
    failures = target-successes-censored-pending if failures is None else failures
    return {"name": "contact", "target_trials": target,
            "attempted": successes+failures+censored+pending,
            "successes": successes, "failures": failures,
            "censored": censored, "pending": pending}


def replica(identity, successes, *, steps=1000, crossings=10, complete=True):
    return {'replica_id': identity, 'completed': complete,
            'flux': {'bank': [f'state_{i}' for i in range(crossings)],
                     'exposure_steps': steps, 'complete': complete},
            'stages': [stage(successes)]}


def plan():
    return {'stem_pairs': [[i, 99-i] for i in range(15)],
            'interfaces': [{'name': 'opening', 'min_binding': 30, 'max_stem': 0, 'max_opening': 0}],
            'basin_distance': 9.}


class FixedTrialStatisticsTests(unittest.TestCase):
    def test_product_uses_fixed_denominators_without_pooling_interfaces(self):
        result = fixed_trial_probability([stage(20), stage(5, target=10)])
        self.assertEqual(result["status"], "complete")
        self.assertAlmostEqual(result["beta"], .1)
        self.assertAlmostEqual(result["log_beta"], math.log(.1))

    def test_censored_and_pending_never_count_as_failures(self):
        for change in ({"censored":1}, {"pending":1}):
            result = fixed_trial_probability([stage(50, **change)])
            self.assertEqual(result["status"], "incomplete")
            self.assertIsNone(result["beta"])
            self.assertIsNone(result["stages"][0]["probability"])
        short = fixed_trial_probability([stage(50, failures=0)])
        self.assertIsNone(short["beta"])

    def test_zero_prefix_does_not_invent_downstream_trials(self):
        result = fixed_trial_probability([stage(50), stage(0), stage(0, failures=0)])
        self.assertEqual(result["status"], "no_success")
        self.assertEqual(result["beta"], 0.)
        self.assertEqual(result["zero_success_stage"], 1)
        self.assertEqual(result["stages"][2]["attempted"], 0)
        self.assertIsNone(result["stages"][2]["probability"])
        with self.assertRaises(ValueError):
            fixed_trial_probability([stage(0), stage(1)])

    def test_invalid_or_unlabelled_trials_are_rejected(self):
        bad = [[], [stage(5)|{"target_trials":False}], [stage(5)|{"successes":5.0}],
               [stage(5)|{"attempted":101}], [stage(5)|{"failures":94}],
               [stage(5)|{"censored":-1}], [{"successes":5}]]
        for stages in bad:
            with self.subTest(stages=stages), self.assertRaises(ValueError):
                fixed_trial_probability(stages)

    def test_underflow_is_not_misreported_as_observed_zero_success(self):
        result = fixed_trial_probability([stage(1, target=1000000)]*100)
        self.assertEqual(result["status"], "underflow")
        self.assertIsNone(result["beta"])
        self.assertTrue(math.isfinite(result["log_beta"]))


class NativeProtocolTests(unittest.TestCase):
    def test_only_native_version_and_supported_modes_are_accepted(self):
        self.assertEqual(OPENING_PROTOCOL_VERSION, 3)
        validate_kinetic_protocol({'mode':'opening','version':3})
        for options in ({}, {'mode': 'encounter', 'version': 3}, {'mode': 'association', 'version': 3}, {'mode': 'opening', 'version': 2}):
            with self.assertRaises(ValueError):
                validate_kinetic_protocol(options)

    def test_single_route_covers_approach_all_thirty_pairs_and_opening(self):
        interfaces = candidate_interfaces(plan())
        self.assertEqual([row['distance_max'] for row in interfaces if 'distance_max' in row],
                         [8.,7.,6.,5.,4.,3.,2.])
        self.assertEqual([row['min_binding'] for row in interfaces if row['name'].startswith('duplex_')],
                         list(range(1,31)))
        self.assertEqual(interfaces[-1], plan()['interfaces'][-1])
        self.assertEqual(interfaces[-1]['max_stem'], 0)

    def test_progress_is_nested_and_final_requires_all_on_criteria(self):
        interfaces = candidate_interfaces(plan())
        self.assertEqual(progress_index({'distance': 10., 'binding': 0}, interfaces), -1)
        self.assertEqual(progress_index({'distance': 7.5, 'binding': 0}, interfaces), 0)
        self.assertEqual(progress_index({'distance': 12., 'binding': 1}, interfaces), 7)
        final = {'distance': 1., 'binding': 30, 'stem': 0, 'opening': 0}
        self.assertEqual(progress_index(final, interfaces), len(interfaces)-1)
        for change in ({'binding': 29}, {'stem': 1}, {'opening': 1}):
            self.assertLess(progress_index({**final, **change}, interfaces), len(interfaces)-1)

    def test_native_condition_preserves_all_conjuncts(self):
        self.assertEqual(condition(plan()['interfaces'][-1]), ['binding >= 30', 'stem <= 0', 'opening <= 0'])
        self.assertEqual(condition({'min_binding': 0, 'distance_max': 8.}), ['distance <= 8.0'])


class DirectRateTests(unittest.TestCase):
    def test_modified_test_run_never_reports_kinetic_estimates(self):
        rows = [replica(str(i), value) for i, value in enumerate((20, 40, 60, 80))]
        original = deepcopy(rows)
        amendments = [{'reason': 'Changed approach coordinate during test run'}]
        for states in (rows, rows[:1], []):
            result = direct_ffs_report(states, dt=.01, concentration_m=1e-6,
                                       protocol_amendments=amendments)
            self.assertEqual(result['status'], 'exploratory_modified_protocol')
            for key in ('value_model', 'inverse_model_time', 'interval95_model',
                        'value_m_inv_model_time', 'uncertainty_method'):
                self.assertIsNone(result[key])
            for row in result['replicas']:
                self.assertIsNone(row['flux'])
                self.assertIsNone(row['probability'])
                self.assertIsNone(row['value_model'])
                self.assertTrue(all(stage['probability'] is None for stage in row['stages']))
            self.assertEqual(result['protocol_amendments'], amendments)
        self.assertEqual(rows, original)

    def report(self, rows, *, dt=.01, concentration=1e-6):
        return direct_ffs_report(rows, dt=dt, concentration_m=concentration)

    def test_report_rejects_invalid_physical_scales_or_duplicate_replicas(self):
        for kwargs in ({'dt':0}, {'dt':float('nan')}, {'concentration':0}, {'concentration':float('inf')}):
            with self.subTest(kwargs=kwargs), self.assertRaises(ValueError):
                self.report([replica('one',50)], **kwargs)
        with self.assertRaises(ValueError):
            self.report([replica('same',50),replica('same',60)])
        row = replica('one',50)
        row['flux']['complete'] = False
        self.assertEqual(self.report([row])['status'],'incomplete')
        self.assertIsNone(self.report([row])['value_model'])

    def test_rate_is_mean_of_matched_flux_probability_products(self):
        rows = [replica('fast', 20, steps=100), replica('slow', 80, steps=1000)]
        result = self.report(rows)
        expected = statistics.fmean([10*.2, 1*.8])
        self.assertAlmostEqual(result['value_model'], expected)
        wrong = statistics.fmean([10., 1.])*statistics.fmean([.2, .8])
        self.assertNotAlmostEqual(result['value_model'], wrong)
        self.assertEqual(result['formula'], 'mean_r(flux_r * product_i(P_r_i))')

    def test_all_conditional_stages_belong_to_the_same_replica_product(self):
        row = replica('first', 25)
        row['stages'].append(stage(20, target=40))
        result = self.report([row])
        self.assertAlmostEqual(result['value_model'], .125)
        self.assertAlmostEqual(result['replicas'][0]['probability'], .125)

    def test_time_normalization_and_finite_box_concentration_have_no_seconds_conversion(self):
        rows = [replica('one', 50)]
        result = self.report(rows)
        self.assertAlmostEqual(result['replicas'][0]['a_history_time'], 10.)
        self.assertAlmostEqual(result['value_model'], .5)
        self.assertAlmostEqual(result['value_m_inv_model_time'], 500000.)
        self.assertEqual(result['unit'], '1/oxDNA_time')
        self.assertIsNone(result['experimental_rate'])
        self.assertTrue(result['finite_box_estimate'])
        self.assertAlmostEqual(self.report(rows, dt=.02)['value_model'], result['value_model']/2)
        self.assertAlmostEqual(self.report(rows, concentration=2e-6)['value_m_inv_model_time'], 250000.)

    def test_zero_success_replicas_are_retained_and_all_zero_is_not_certainty(self):
        result = self.report([replica('zero', 0), replica('half', 50)])
        self.assertAlmostEqual(result['value_model'], .25)
        result = self.report([replica(str(i), 0) for i in range(4)])
        self.assertEqual(result['status'], 'no_success')
        self.assertEqual(result['value_model'], 0.)
        self.assertIsNone(result['interval95_model'])

    def test_pending_or_incomplete_replicas_prevent_cherry_picked_rate(self):
        for pending in (True, False):
            rows = [replica('finished', 50), replica('unfinished', 50)]
            if pending:
                rows[1]['stages'] = [stage(50, pending=1)]
            else:
                rows[1]['completed'] = False
            result = self.report(rows)
            self.assertEqual(result['status'], 'incomplete')
            self.assertIsNone(result['value_model'])
            self.assertIsNone(result['interval95_model'])
        self.assertEqual(self.report([])['status'], 'incomplete')

    def test_uncertainty_uses_whole_replicas_not_number_of_descendants(self):
        rows = [replica(str(i), successes) for i, successes in enumerate((20, 40, 60, 80))]
        result = self.report(rows)
        values = [row['value_model'] for row in result['replicas']]
        self.assertAlmostEqual(result['interval95_model'][1], statistics.fmean(values) + 3.182446*statistics.stdev(values)/2)
        expanded = deepcopy(rows)
        for row in expanded:
            for field in ('successes', 'failures', 'attempted', 'target_trials'):
                row['stages'][0][field] *= 1000
        self.assertEqual(self.report(expanded)['interval95_model'], result['interval95_model'])

    def test_short_or_identical_replicas_do_not_manufacture_zero_width_uncertainty(self):
        for rows in ([replica('single', 50)], [replica(str(i), 50) for i in range(4)],
                     [replica('a', 20), replica('b', 80)]):
            self.assertIsNone(self.report(rows)['interval95_model'])


class SourceProvenanceTests(unittest.TestCase):
    def test_modified_run_cannot_supply_source_for_fresh_kinetic_run(self):
        original = self.request()
        changed = deepcopy(original)
        changed['protocol_amendments'] = [{'reason': 'Exploratory continuation'}]
        self.assertNotEqual(source_fingerprint(original), source_fingerprint(changed))

    def request(self):
        return {'saved_switch_identity': 'candidate', 'switch_sequence': 'ACGUACGU',
            'trigger_sequence': 'ACGT', 'switch_material': 'rna', 'trigger_material': 'dna',
            'model': 'NA', 'box_size': 100., 'switch_guide': '........', 'trigger_guide': '....',
            'plan': plan(), 'sampling': {'temperature_c': 37., 'salt_m': .12, 'dt': .003, 'seed': 2},
            'runtime': {'parameter_sha256': {'DNA': 'dna_hash', 'RNA': 'rna_hash', 'NA': 'hybrid_hash'}},
            'kinetics': {'mode': 'opening', 'version': 3, 'trials': 1000, 'replicas': 4}}

    def test_source_fingerprint_excludes_seed_and_allocation(self):
        request = self.request()
        changed = deepcopy(request)
        changed['kinetics'].update(trials=2000, workers=2)
        changed['sampling']['seed'] = 456
        changed['runtime']['home'] = '/moved/oxdna'
        self.assertEqual(source_fingerprint(request), source_fingerprint(changed))
        explicit_defaults = deepcopy(request)
        explicit_defaults['sampling'].update(newtonian_steps=103, diff_coeff=2.5)
        self.assertEqual(source_fingerprint(request), source_fingerprint(explicit_defaults))

    def test_source_fingerprint_rejects_physical_condition_changes(self):
        request = self.request()
        changes = [('box_size', 110.), ('switch_sequence', 'UGCAUGCA'), ('switch_material', 'dna'),
                   ('switch_guide', '(......)'), ('saved_switch_identity', 'another_candidate')]
        for field, value in changes:
            self.assertNotEqual(source_fingerprint(request), source_fingerprint({**request, field: value}))
        for field, value in [('temperature_c', 40.), ('salt_m', .2), ('dt', .004),
                             ('newtonian_steps', 104), ('diff_coeff', 3.)]:
            changed = deepcopy(request)
            changed['sampling'][field] = value
            self.assertNotEqual(source_fingerprint(request), source_fingerprint(changed))
        for kind in ('basin', 'parameters'):
            changed = deepcopy(request)
            if kind == 'basin':
                changed['plan']['basin_distance'] = 6.
            else:
                changed['runtime']['parameter_sha256']['NA'] = 'different'
            self.assertNotEqual(source_fingerprint(request), source_fingerprint(changed))


if __name__ == '__main__':
    unittest.main()
