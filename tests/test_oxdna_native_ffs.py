"""Native FFS serialization and opt-in actual CPU oxNA regression checks."""
from copy import deepcopy
import json
import os
from pathlib import Path
import tempfile
import unittest

from aptaswitch_core.oxdna_ffs_worker import Observation, OxpyAdapter, _same_checkpoint, validate_request
from aptaswitch_core.oxdna_native_ffs import (
    NativeFFSAdapter, conditions_text, first_condition, order_parameters_text,
)


class NativeFFSInputsTests(unittest.TestCase):
    def test_runtime_source_fix_is_idempotent_and_refuses_unknown_backend(self):
        import importlib.util
        script = Path(__file__).resolve().parents[1] / 'scripts/prepare_oxdna_ffs_runtime.py'
        spec = importlib.util.spec_from_file_location('native_ffs_runtime_fix', script)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            backends = root / 'src/Backends'
            backends.mkdir(parents=True)
            (backends / 'MD_CPUBackend.cpp').write_text(
                '_interaction->begin_energy_and_force_computation();\n'
                'p->set_initial_forces(current_step(), _box.get());\n')
            native = backends / 'FFS_MD_CPUBackend.cpp'
            original = 'void FFS_MD_CPUBackend::_ffs_compute_forces(void) {\n\t_interaction->begin_energy_computation();\n}'
            native.write_text(original)
            self.assertTrue(module.repair_source(root))
            self.assertFalse(module.repair_source(root))
            self.assertEqual(native.with_suffix('.cpp.before_force_reset_fix').read_text(), original)
            self.assertTrue((root / 'aptaswitch_ffs_force_reset.json').is_file())
            native.write_text('unknown backend')
            with self.assertRaises(ValueError):
                module.repair_source(root)
            self.assertEqual(native.read_text(), 'unknown backend')

    def test_and_blocks_are_or_conditions_and_inputs_are_not_executable(self):
        conditions = [["binding >= 2", "opening == 0"], ["associations == 0", "distance >= 9"]]
        text = conditions_text(conditions)
        self.assertIn("condition1 = {\n    binding >= 2\n    opening == 0\n}", text)
        self.assertIn("condition2 = {\n    associations == 0\n    distance >= 9\n}", text)
        self.assertEqual(first_condition(Observation(3, 1, 0, 1, 1, 3), conditions), 0)
        self.assertEqual(first_condition(Observation(0, 15, 6, 10, 1, 0), conditions), 1)
        self.assertIsNone(first_condition(Observation(0, 15, 6, 8, 1, 0), conditions))
        for bad in ("binding >= 1; distance < 3", "other == 1", "binding >= .5",
                    "opening == -1", "distance >= 1e309"):
            with self.subTest(bad=bad), self.assertRaises(ValueError):
                conditions_text([[bad]])
        with self.assertRaises(ValueError):
            conditions_text([[]])
        self.assertEqual(conditions_text([]), "")  # Uncensored equilibration.

    def test_order_parameters_include_all_interstrand_bonds_and_explicit_com(self):
        plan = {"switch_length": 2, "nucleotide_count": 4,
                "binding_pairs": [[0, 3], [1, 2]], "stem_pairs": [], "opening_pairs": [[0, 1]]}
        text = order_parameters_text(plan)
        self.assertIn("name = associations", text)
        self.assertIn("pair4 = 1, 3", text)
        self.assertIn("use_COM = true", text)
        self.assertIn("cutoff = -0.1", text)
        with self.assertRaises(ValueError):
            order_parameters_text({**plan, "hb_threshold": -.2})
        with self.assertRaises(ValueError):
            order_parameters_text({**plan, "binding_pairs": [[0, 4]]})

    def test_native_approach_and_separation_use_distinct_exact_pair_sets(self):
        from test_oxdna_ffs_worker import toehold_request
        plan = toehold_request()['plan']
        text = order_parameters_text(plan)
        blocks = text.split('\n\n')
        approach = next(block for block in blocks if 'name = distance\n' in block)
        separation = next(block for block in blocks if 'name = separation\n' in block)
        self.assertIn('pair1 = 0, 2\n    pair2 = 0, 3', approach)
        self.assertNotIn('pair3', approach)
        self.assertIn('pair4 = 1, 3', separation)
        self.assertEqual(text.count('use_COM = true'), 2)
        conditions = [['associations == 0', 'separation >= 9.5']]
        self.assertIn('separation >= 9.5', conditions_text(conditions))
        self.assertIsNone(first_condition(Observation(0, 1, 0, 10., 100, 0, 2.), conditions))
        self.assertEqual(first_condition(Observation(0, 1, 0, 10., 100, 0, 10.), conditions), 0)
        with self.assertRaisesRegex(ValueError, 'missing separation'):
            first_condition(Observation(0, 1, 0, 10., 100, 0), conditions)

    def test_terminal_approach_serializes_only_the_antiparallel_entry_pair(self):
        from test_oxdna_ffs_worker import toehold_request
        plan = toehold_request()['plan']
        plan.update(distance_scope='trigger_3prime_to_toehold_start', distance_pairs=[[0, 3]])
        text = order_parameters_text(plan)
        blocks = text.split('\n\n')
        approach = next(block for block in blocks if 'name = distance\n' in block)
        separation = next(block for block in blocks if 'name = separation\n' in block)
        self.assertIn('pair1 = 0, 3', approach)
        self.assertNotIn('pair2', approach)
        self.assertIn('use_COM = true', approach)
        self.assertIn('pair1 = 0, 2', separation)
        self.assertIn('pair4 = 1, 3', separation)
        self.assertEqual(text.count('use_COM = true'), 2)


@unittest.skipUnless(os.environ.get("APTA_NATIVE_FFS_REQUEST"), "requires external oxpy and a prepared request")
class ActualNativeFFSTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.directory = Path(self.temporary.name)
        self.request = validate_request(json.loads(Path(os.environ["APTA_NATIVE_FFS_REQUEST"]).read_text()))
        self.configuration = Path(self.request["prepared"]["configuration_file"])

    def test_ten_random_separated_replicas_are_native_valid_without_restraints(self):
        from aptaswitch_core.oxdna_kinetics_worker import FFSReplica
        request = deepcopy(self.request)
        request['kinetics'] = {'version': 3, 'mode': 'opening', 'source_samples': 200,
            'trials': 1000, 'replicas': 10, 'workers': 10, 'equilibration_steps': 1000, 'chunk_steps': 5}
        initial_frames = set()
        for index in range(10):
            replica = FFSReplica(request, self.directory/f'replica_{index+1:03d}', index)
            initial_frames.add(replica.state['frame'])
            self.assertGreaterEqual(replica.state['initialization']['minimum_interstrand_distance'],
                                    request.get('initial_separation', 10.))
            path = replica.output/'random_initial.dat'
            path.write_text(replica.state['frame'])
            adapter = NativeFFSAdapter(replica.request, path, 31+index, True,
                                       directory=replica.output/'native', conditions=[])
            try:
                initial = adapter.observe()
                self.assertEqual(initial.associations, 0)
                self.assertGreaterEqual(initial.separation or initial.distance,
                                        request.get('initial_separation', 10.))
                result = adapter.advance(5)
                self.assertEqual(result.steps, 5)
                self.assertEqual(result.observation.associations, 0)
            finally:
                adapter.close()
        self.assertEqual(len(initial_frames), 10)

    def test_both_native_distances_match_observations_at_first_crossing(self):
        # Keep these short comparisons between thermostat refreshes: oxpy's
        # global Gaussian cache is not reset when another manager is opened.
        configuration = self.directory / 'distance_comparison.dat'
        rows = self.configuration.read_text().splitlines()
        rows[0] = 't = 1'
        configuration.write_text('\n'.join(rows) + '\n')
        for scope in ('trigger_to_toehold', 'trigger_3prime_to_toehold_start'):
            request = deepcopy(self.request)
            plan = request['plan']
            n, total = plan['switch_length'], plan['nucleotide_count']
            toehold = plan.get('toehold_pairs') or plan['binding_pairs'][:1]
            exposed = sorted({pair[0] for pair in toehold})
            approach = ([[exposed[0], total - 1]] if scope == 'trigger_3prime_to_toehold_start'
                        else [[i, j] for i in exposed for j in range(n, total)])
            plan.update(distance_scope=scope, toehold_pairs=toehold, distance_pairs=approach,
                separation_pairs=[[i, j] for i in range(n) for j in range(n, total)])
            request = validate_request(request)
            plain = OxpyAdapter(request, configuration, 31, False)
            try:
                observations = [plain.observe()]
                for _ in range(25):
                    observations.append(plain.step())
            finally:
                plain.close()
            for name in ('distance', 'separation'):
                with self.subTest(scope=scope, name=name):
                    values = [getattr(item, name) for item in observations]
                    decreasing = min(values[1:]) < values[0]
                    threshold = (values[0] + (min(values[1:]) if decreasing else max(values[1:]))) / 2
                    operator = '<=' if decreasing else '>='
                    expected = next(i for i, value in enumerate(values)
                                    if (value <= threshold if decreasing else value >= threshold))
                    self.assertGreater(expected, 0)
                    native = NativeFFSAdapter(request, configuration, 31, False,
                        directory=self.directory / scope / name,
                        conditions=[[f'{name} {operator} {threshold:.17g}']])
                    try:
                        result = native.advance(100)
                        self.assertEqual((result.steps, result.condition), (expected, 0))
                        self.assertAlmostEqual(getattr(result.observation, name), values[expected], places=8)
                        self.assertIn('Reached condition1', native.log_file.read_text())
                    finally:
                        native.close()

    def test_terminal_native_threshold_ignores_near_5prime_end_including_periodic_images(self):
        from test_oxdna_ffs_worker import toehold_request
        request = toehold_request()
        request['runtime'] = self.request['runtime']
        request['plan'].update(distance_scope='trigger_3prime_to_toehold_start', distance_pairs=[[0, 3]])
        topology = self.directory / 'terminal.top'
        topology.write_text('4 2 5->3\nAC type=RNA\nGT type=DNA\n')
        request['prepared']['topology_file'] = str(topology)
        request = validate_request(request)
        # Two intact short strands. Moving the whole trigger by one box tests
        # minimum-image centres without changing its internal geometry.
        for near_end, positions, expected_distance in (
            ('5prime', (3., 3.75), 3.75), ('3prime', (3.75, 3.), 3.),
        ):
            for offset in (0., 40.):
                with self.subTest(near_end=near_end, periodic_offset=offset):
                    configuration = self.directory / f'{near_end}_{offset}.dat'
                    coordinates = (0., .75, *(x + offset for x in positions))
                    configuration.write_text('t = 1\nb = 40 40 40\nE = 0 0 0\n' + ''.join(
                        f'{x} 0 0 1 0 0 0 0 1 0 0 0 0.01 0.02 0.03\n' for x in coordinates))
                    native = NativeFFSAdapter(request, configuration, 31, False,
                        directory=self.directory / f'native_{near_end}_{offset}',
                        conditions=[['distance <= 3.4']])
                    try:
                        observation = native.observe()
                        self.assertAlmostEqual(observation.distance, expected_distance, places=9)
                        self.assertLess(observation.separation, 3.4)
                        result = native.advance(1)
                        self.assertEqual(result.condition, 0 if near_end == '3prime' else None)
                        self.assertEqual(result.steps, 0 if near_end == '3prime' else 1)
                    finally:
                        native.close()

    def test_NA_FFS_MD_matches_standard_MD_including_momenta_and_salt_forces(self):
        # An inherited DNA2Interaction supplies the NA Debye-Hückel term in
        # FFS_MD. Compare real integrator trajectories, not just input labels.
        # Avoid thermostat refresh during this deterministic force comparison:
        # oxDNA's global Gaussian cache is not reset between manager instances.
        configuration = self.directory / 'force_comparison.dat'
        rows = self.configuration.read_text().splitlines()
        rows[0] = 't = 1'
        configuration.write_text('\n'.join(rows) + '\n')
        for salt in (.12, .5):
            request = deepcopy(self.request)
            request["sampling"]["salt_m"] = salt
            plain = OxpyAdapter(request, configuration, 31, False)
            try:
                plain.manager.run(25, print_output=False)
                expected, energy = plain.frame(), plain.manager.system_energy()
            finally:
                plain.close()
            native = NativeFFSAdapter(request, configuration, 31, False,
                                      directory=self.directory / str(salt), conditions=[])
            try:
                chunk = native.advance(25)
                self.assertEqual(chunk.steps, 25)
                self.assertIsNone(chunk.condition)
                self.assertTrue(_same_checkpoint(expected, native.frame()))
                self.assertAlmostEqual(energy, native.manager.system_energy(), places=9)
            finally:
                native.close()

    def test_native_stops_inside_chunk_and_preserves_exact_phase_point_on_restart(self):
        # First find a real target distance decrease without changing geometry.
        plain = OxpyAdapter(self.request, self.configuration, 31, False)
        try:
            observations = [plain.observe()]
            for _ in range(25):
                observations.append(plain.step())
        finally:
            plain.close()
        initial = observations[0].distance
        threshold = (initial + min(item.distance for item in observations[1:])) / 2
        self.assertLess(threshold, initial, "Fixture must have a short distance decrease.")
        expected = next(index for index, item in enumerate(observations) if item.distance <= threshold)
        native = NativeFFSAdapter(self.request, self.configuration, 31, False,
                                  directory=self.directory / "stop",
                                  conditions=[["binding >= 30"], [f"distance <= {threshold:.17g}"]])
        try:
            chunk = native.advance(100)
            self.assertEqual(chunk.condition, 1)
            self.assertEqual(chunk.steps, expected)
            self.assertEqual(chunk.observation.step, observations[expected].step)
            self.assertIn("Reached condition2", native.log_file.read_text())
            frame = native.frame()
            self.assertEqual(native.advance(100).steps, 0)
        finally:
            native.close()
        resumed = self.directory / "checkpoint.dat"
        resumed.write_text(frame)
        native = NativeFFSAdapter(self.request, resumed, 37, False,
                                  directory=self.directory / "resume", conditions=[])
        try:
            self.assertTrue(_same_checkpoint(frame, native.frame()))
            continued = native.advance(5)
            self.assertEqual(continued.steps, 5)
            self.assertEqual(continued.observation.step, observations[expected].step + 5)
        finally:
            native.close()

    def test_initial_boundary_is_returned_without_one_step_overshoot(self):
        native = NativeFFSAdapter(self.request, self.configuration, 31, False,
                                  directory=self.directory / "initial", conditions=[["distance >= 0"]])
        try:
            initial = native.frame()
            result = native.advance(100)
            self.assertEqual((result.steps, result.condition), (0, 0))
            self.assertTrue(_same_checkpoint(initial, native.frame()))
        finally:
            native.close()

    def test_replica_pause_restores_native_phase_point_without_counting_failure(self):
        from aptaswitch_core.oxdna_kinetics_worker import FFSReplica
        request = deepcopy(self.request)
        request['kinetics'] = {'version': 3, 'mode': 'opening', 'source_samples': 200,
            'trials': 1000, 'replicas': 4, 'workers': 2, 'equilibration_steps': 1000, 'chunk_steps': 5}
        output = self.directory / 'replica_001'
        opened = []

        def stop_after_one_chunk(request, configuration, seed, initial, **kwargs):
            adapter = NativeFFSAdapter(request, configuration, seed, initial, **kwargs)
            opened.append(adapter.frame())
            advance = adapter.advance

            def one_chunk(steps):
                result = advance(steps)
                (self.directory / 'CANCEL').touch()
                return result

            adapter.advance = one_chunk
            return adapter

        first = FFSReplica(request, output, 0, adapter_factory=stop_after_one_chunk).run()
        self.assertEqual(first['status'], 'cancelled', first.get('error'))
        self.assertEqual(first['steps'], 5)
        self.assertEqual(first['stages'], [])
        self.assertFalse(first['flux']['complete'])
        checkpoint = first['frame']
        (self.directory / 'CANCEL').unlink()
        second = FFSReplica(request, output, 0, adapter_factory=stop_after_one_chunk).run()
        self.assertEqual(second['status'], 'cancelled', second.get('error'))
        self.assertEqual(second['steps'], 10)
        self.assertEqual(second['burn_steps'], 10)
        self.assertTrue(_same_checkpoint(checkpoint, opened[1]))
        self.assertEqual(second['flux']['exposure_steps'], 0)

    def test_two_actual_native_workers_pause_and_resume_independently(self):
        import threading
        import time
        from aptaswitch_core.oxdna_kinetics_worker import KineticSampler
        request = deepcopy(self.request)
        request['kinetics'] = {'version': 3, 'mode': 'opening', 'source_samples': 200,
            'trials': 1000, 'replicas': 4, 'workers': 2, 'equilibration_steps': 200000, 'chunk_steps': 32}
        output = self.directory / 'parallel'
        output.mkdir()
        (output / 'worker_request.json').write_text(json.dumps(request))

        def run_until_both_advanced(previous):
            deadline = time.monotonic() + 15
            advanced = threading.Event()

            def watch():
                while time.monotonic() < deadline:
                    rows = []
                    for index in range(2):
                        try:
                            rows.append(json.loads((output / f'replica_{index+1:03d}/state.json').read_text()))
                        except (OSError, ValueError):
                            rows.append({})
                    if all(row.get('steps', 0) > previous[index] for index, row in enumerate(rows)):
                        advanced.set()
                        break
                    time.sleep(.02)
                (output / 'CANCEL').touch()

            monitor = threading.Thread(target=watch, daemon=True)
            monitor.start()
            result = KineticSampler(request, output).run()
            monitor.join(timeout=16)
            self.assertTrue(advanced.is_set(), result.get('error'))
            self.assertEqual(result['status'], 'cancelled', result.get('error'))
            states = [json.loads((output / f'replica_{i+1:03d}/state.json').read_text()) for i in range(2)]
            self.assertTrue(all(state['steps'] > previous[index] for index, state in enumerate(states)))
            return [row['steps'] for row in states]

        first = run_until_both_advanced([0, 0])
        for index in range(2):
            settings = (output / f'replica_{index+1:03d}/native/native_input.txt').read_text()
            self.assertIn('sim_type = FFS_MD', settings)
        (output / 'CANCEL').unlink()
        second = run_until_both_advanced(first)
        self.assertTrue(all(new > old for old, new in zip(first, second)))


if __name__ == "__main__":
    unittest.main()
