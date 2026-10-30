"""Short native checks, exact guides, and bounded preparation continuation."""
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from aptaswitch_core import oxdna_guided_worker as worker


def request():
    return {
        'seed': 37, 'configuration_file': 'initial.dat', 'switch_sequence': 'GAAAC',
        'switch_material': 'rna', 'trigger_sequence': 'AT', 'trigger_material': 'dna',
        'switch_guide': '(...)', 'trigger_guide': '..',
        'guides': {'switch': {'pairs': [[0, 4]]}, 'trigger': {'pairs': []}},
        'plan': {'stem_pairs': [[0, 4]], 'basin_min_stem': 1, 'basin_distance': 4.},
        'temperature_c': 37., 'salt_m': .12, 'initial_separation': 6.,
        'preparation_attempt': {'index': 2, 'seed': 37, 'restraint_stiffness': 3.},
    }


class FakeNative:
    """Reject simultaneous native contexts and distinguish guided from probes."""
    def __init__(self, probe_bonds, initial_exact=True, guide_becomes_exact=250):
        self.active = None
        self.instances = []
        self.probe_bonds = list(probe_bonds)
        self.initial_exact = initial_exact
        self.guide_becomes_exact = guide_becomes_exact
        self.probes_started = 0
        self.placements = []
        self.geometry = lambda instance: {'passed': True}
        self.loops = lambda instance: {'passed': True}

    def separate(self, source, target, switch_length, separation, **kwargs):
        if self.active is not None:
            raise AssertionError('Preparation context must close before translation/probe.')
        self.placements.append((str(source), str(target)))
        target.write_text(source.read_text())
        return {'box_size': 80., 'intrachain_geometry_changed': False}

    def stage(self, req, configuration, **options):
        if self.active is not None:
            raise AssertionError('Two oxpy contexts are alive simultaneously.')
        native = self
        class Instance:
            def __init__(self):
                self.source = Path(configuration)
                self.options = options
                self.steps = 0
                self.is_probe = self.source.name.startswith('probe_')
                self.probe_index = None
                if self.is_probe:
                    self.probe_index = native.probes_started
                    native.probes_started += 1
            def run(self, steps, deadline):
                self.steps += steps
            def observe(self, request):
                if self.is_probe:
                    assert not self.options.get('relaxed', False)
                    assert self.options.get('forces') is None
                    bonds = native.probe_bonds[min(self.probe_index, len(native.probe_bonds)-1)]
                else:
                    bonds = {(0, 4)} if native.initial_exact or self.steps >= native.guide_becomes_exact else set()
                # Ignore preparation distance; translation precedes actual OFF validation.
                distance = 6. if self.is_probe else .5
                observation = worker.validation_from_bonds(request, bonds, distance)
                observation['geometry'] = native.geometry(self)
                observation['loop_threading'] = native.loops(self)
                observation['ffs_eligible'] &= observation['geometry']['passed'] and observation['loop_threading']['passed']
                return observation
            def frame(self):
                return f'frame {self.source.name} steps {self.steps}\n'
            def close(self):
                if native.active is not self:
                    raise AssertionError('Context close is not balanced.')
                native.active = None
        instance = Instance()
        self.active = instance
        self.instances.append(instance)
        return instance


class GuidedWorkerTests(unittest.TestCase):
    def execute(self, native, req=None):
        req = req or request()
        with tempfile.TemporaryDirectory() as folder:
            output = Path(folder)
            (output/'initial.dat').write_text('initial seed\n')
            req = dict(req, configuration_file=str(output/'initial.dat'))
            with patch.object(worker, 'Stage', side_effect=native.stage), \
                 patch.object(worker, 'separate_monomers', side_effect=native.separate):
                result = worker.run(req, output)
            self.assertIsNone(native.active)
            self.assertEqual(json.loads((output/'guided_stages.json').read_text()), result['relaxation']['stages'])
            self.assertEqual(json.loads((output/'guided_probes.json').read_text()), result['relaxation']['probes'])
            result['force_text'] = (output/'guide_forces.txt').read_text()
            result['final_frame'] = Path(result['configuration_file']).read_text()
            return result

    def test_exact_seed_uses_real_short_unforced_check_and_stops_early(self):
        native = FakeNative([{(0, 4)}])
        result = self.execute(native)
        self.assertTrue(result['ffs_eligible'])
        self.assertEqual(len(native.instances), 2)
        self.assertEqual(native.instances[0].steps, 0)
        self.assertEqual(native.instances[1].steps, 1000)
        self.assertEqual(result['selected_checkpoint'], {'stage': 'minimize_guided', 'step': 0, 'probe_index': 1})
        self.assertFalse(result['relaxation']['equilibrium_validated'])
        self.assertFalse(result['production_external_forces'])
        self.assertEqual(result['preparation_attempt']['index'], 2)
        self.assertIn('stiff = 3', result['force_text'])
        self.assertIn('steps 1000', result['final_frame'])

    def test_probe_failure_resumes_original_checkpoint_with_new_seed(self):
        native = FakeNative([set(), {(0, 4)}])
        result = self.execute(native)
        self.assertTrue(result['ffs_eligible'])
        stages = result['relaxation']['stages']
        self.assertEqual(stages[0]['steps'], 250)
        self.assertEqual([probe['accepted'] for probe in result['relaxation']['probes']], [False, True])
        self.assertEqual(native.instances[2].source.name, 'minimize_guided.dat')
        self.assertFalse(native.instances[2].is_probe)
        self.assertEqual(len({instance.options['seed'] for instance in native.instances}), len(native.instances))
        repeat = self.execute(FakeNative([set(), {(0, 4)}]))
        self.assertEqual([p['seed'] for p in result['relaxation']['probes']], [p['seed'] for p in repeat['relaxation']['probes']])

    def test_minimization_samples_early_instead_of_waiting_full_stage(self):
        native = FakeNative([{(0, 4)}], initial_exact=False)
        result = self.execute(native)
        self.assertEqual(result['selected_checkpoint']['step'], 250)
        self.assertEqual([row['step'] for row in result['relaxation']['stages'][0]['checkpoints']], [0, 250])
        self.assertEqual(len(result['relaxation']['stages']), 1)

    def test_failed_strict_initialization_does_not_accept_restrained_candidate(self):
        native = FakeNative([{(0, 4)}])
        real_factory = native.stage
        probes = 0
        def fail_first(req, configuration, **options):
            nonlocal probes
            if Path(configuration).name.startswith('probe_'):
                probes += 1
                if probes == 1:
                    self.assertIsNone(native.active)
                    raise RuntimeError('FENE out of range')
            return real_factory(req, configuration, **options)
        native.stage = fail_first
        result = self.execute(native)
        self.assertTrue(result['ffs_eligible'])
        self.assertEqual(result['relaxation']['probes'][0]['error'], 'FENE out of range')
        self.assertFalse(result['relaxation']['probes'][0]['accepted'])
        self.assertEqual(result['selected_checkpoint']['step'], 250)

    def test_extra_pairs_and_linear_guides_require_exact_native_pair_set(self):
        req = request()
        observation = worker.validation_from_bonds(req, {(0, 4), (1, 3)}, 6.)
        self.assertFalse(observation['guide_accepted'])
        self.assertFalse(observation['ffs_eligible'])
        req['guides']['switch']['pairs'] = []
        req['plan']['stem_pairs'] = []
        req['plan']['basin_min_stem'] = 0
        self.assertTrue(worker.validation_from_bonds(req, set(), 6.)['ffs_eligible'])
        self.assertFalse(worker.validation_from_bonds(req, {(0, 4)}, 6.)['ffs_eligible'])

    def test_all_invalid_unforced_probes_keep_diagnostic_and_do_not_accept(self):
        native = FakeNative([set()])
        result = self.execute(native)
        self.assertFalse(result['guide_accepted'])
        self.assertFalse(result['ffs_eligible'])
        self.assertTrue(result['strict_force_field_initialization_passed'])
        self.assertEqual(result['status'], 'completed')
        self.assertEqual([stage['steps'] for stage in result['relaxation']['stages']], [2000, 20000, 30000, 10000])
        self.assertLess(len(native.instances), 160)
        self.assertTrue(all(probe['steps'] == 1000 for probe in result['relaxation']['probes']))

    def test_nonfinite_native_distance_is_rejected_before_report_serialization(self):
        for distance in (float('nan'), float('inf')):
            with self.subTest(distance=distance), self.assertRaisesRegex(ValueError, 'non-finite'):
                worker.validation_from_bonds(request(), {(0, 4)}, distance)

    def test_clash_at_intermediate_probe_step_cannot_be_hidden_by_good_final_frame(self):
        native = FakeNative([{(0, 4)}])
        native.geometry = lambda stage: {'passed': not (stage.is_probe and stage.probe_index == 0 and stage.steps == 500)}
        result = self.execute(native)
        probes = result['relaxation']['probes']
        self.assertFalse(probes[0]['accepted'])
        self.assertEqual(probes[0]['steps'], 500)
        self.assertIn('clash', probes[0]['error'])
        self.assertEqual([check['step'] for check in probes[0]['geometry_checks']], [100, 200, 300, 400, 500])
        self.assertTrue(probes[1]['accepted'])
        self.assertEqual(len(probes[1]['geometry_checks']), 10)
        self.assertTrue(result['ffs_eligible'])

    def test_threading_without_steric_clash_rejects_probe_and_retries(self):
        native = FakeNative([{(0, 4)}])
        native.loops = lambda stage: {'passed': not (stage.is_probe and stage.probe_index == 0 and stage.steps == 500)}
        result = self.execute(native)
        first, second = result['relaxation']['probes']
        self.assertFalse(first['accepted'])
        self.assertEqual(first['steps'], 500)
        self.assertTrue(first['geometry_checks'][-1]['geometry']['passed'])
        self.assertFalse(first['geometry_checks'][-1]['loop_threading']['passed'])
        self.assertIn('threads a guide loop', first['error'])
        self.assertTrue(second['accepted'])

    def test_uncheckable_loop_surface_cannot_start_probe(self):
        native = FakeNative([{(0, 4)}])
        native.loops = lambda stage: {'passed': not (stage.is_probe and stage.probe_index == 0), 'unresolved': ['surface']}
        result = self.execute(native)
        self.assertFalse(result['relaxation']['probes'][0]['accepted'])
        self.assertEqual(native.instances[1].steps, 0)
        self.assertTrue(result['ffs_eligible'])

    def test_initial_steric_clash_never_starts_unforced_probe_md(self):
        native = FakeNative([{(0, 4)}])
        native.geometry = lambda stage: {'passed': not (stage.is_probe and stage.probe_index == 0)}
        result = self.execute(native)
        self.assertFalse(result['relaxation']['probes'][0]['accepted'])
        self.assertEqual(native.instances[1].steps, 0)
        self.assertTrue(result['ffs_eligible'])


if __name__ == '__main__':
    unittest.main()
