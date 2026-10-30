"""Native FFS source history, lineage, fixed allocation and resumability."""
from dataclasses import replace
from copy import deepcopy
from pathlib import Path
import tempfile
import unittest

from aptaswitch_core.oxdna_ffs_worker import Observation, SamplingStopped
from aptaswitch_core.oxdna_kinetics_worker import FFSReplica, KineticSampler, read_json
from aptaswitch_core.oxdna_native_ffs import NativeChunk, first_condition


def frame(step=0):
    return f't = {step}\nb = 100 100 100\nE = 0 0 0\n' + ''.join(
        f'{x} {y} 3 1 0 0 0 0 1 0 0 0 0 0 0\n'
        for x, y in ((1, 2), (1, 2.761070781051), (20, 2), (20, 2.7564)))


def request(root):
    (root/'initial.dat').write_text(frame())
    (root/'system.top').write_text('4 2 5->3\nAU type=RNA\nAT type=DNA\n')
    return {'kinetics': {'version': 3, 'mode': 'opening', 'replicas': 4, 'workers': 2,
                        'source_samples': 2, 'trials': 2, 'equilibration_steps': 0, 'chunk_steps': 10},
        'prepared': {'configuration_file': str(root/'initial.dat'), 'topology_file': str(root/'system.top')},
        'switch_guide': '..', 'trigger_guide': '..',
        'runtime': {}, 'sampling': {'seed': 43, 'temperature_c': 37, 'salt_m': .12, 'dt': .003},
        'plan': {'binding_pairs': [[0, 3], [1, 2]], 'toehold_pairs': [[0, 3]], 'stem_pairs': [[0, 1]],
                 'opening_pairs': [[0, 1]], 'switch_length': 2, 'nucleotide_count': 4, 'basin_distance': 9.,
                 'interfaces': [{'name': 'opening', 'min_binding': 2, 'max_stem': 0, 'max_opening': 0}]},
        'box_size': 100., 'concentration_each_m': 1e-6}


class NativeScript:
    """Only native chunk outcomes are mocked; production bookkeeping is real."""
    def __init__(self, events):
        self.events = list(events)
        self.starts = []
        self.on_start = None
        self.on_advance = None

    def __call__(self, request, path, seed, initial, *, directory, conditions):
        self.starts.append({'frame': path.read_text(), 'seed': seed, 'initial': initial, 'conditions': conditions})
        if self.on_start:
            self.on_start()
        return ScriptManager(self, path.read_text(), conditions)


class ScriptManager:
    def __init__(self, script, config, conditions):
        self.script = script
        self.conditions = conditions
        self.step = int(config.splitlines()[0].split('=')[1])

    def advance(self, count):
        if not self.script.events:
            raise AssertionError('Unexpected native simulation beyond scripted event history')
        steps, hit, binding, distance, *extra = self.script.events.pop(0)
        if steps > count:
            raise AssertionError('Script exceeds native chunk size')
        self.step += steps
        if self.script.on_advance:
            self.script.on_advance()
        observation = Observation(binding, 0, 0, distance, self.step, binding,
                                  separation=extra[0] if extra else None)
        if hit == 'auto':
            hit = first_condition(observation, self.conditions)
        return NativeChunk(steps, hit, observation)

    def frame(self):
        return frame(self.step)

    def close(self):
        pass


class NativeProtocolTests(unittest.TestCase):
    def setUp(self):
        self.root = Path(self.enterContext(tempfile.TemporaryDirectory()))
        self.req = request(self.root)

    def replica(self, events):
        self.script = NativeScript(events)
        return FFSReplica(self.req, self.root/'replica_001', 0, adapter_factory=self.script)

    def test_source_rearms_only_after_return_and_counts_excursions(self):
        r = self.replica([(0,0,0,9), (10,0,0,8), (7,0,0,9), (8,0,0,8)])
        r.flux()
        self.assertEqual(r.state['flux']['exposure_steps'], 25)
        self.assertEqual(len(r.state['flux']['bank']), 2)
        self.assertEqual(self.script.starts[2]['conditions'][0], ['associations == 0', 'distance >= 9.0'])
        self.assertTrue(self.script.starts[0]['initial'])
        self.assertTrue(all(not x['initial'] for x in self.script.starts[1:]))
        # Continuation retains exactly the real phase point reached previously.
        self.assertEqual(self.script.starts[2]['frame'], frame(10))
        self.assertEqual(self.script.starts[3]['frame'], frame(17))

    def use_toehold_distance(self):
        self.req['plan'].update(distance_scope='trigger_to_toehold',
                                distance_pairs=[[0, 2], [0, 3]],
                                separation_pairs=[[0, 2], [0, 3], [1, 2], [1, 3]])

    def test_terminal_source_ignores_wrong_end_contact_but_keeps_physical_return_basin(self):
        self.use_toehold_distance()
        self.req['plan'].update(distance_scope='trigger_3prime_to_toehold_start', distance_pairs=[[0, 3]])
        # 5′ is near at step 10, 3′ only arrives at step 20. Recrossing without
        # whole-strand separation cannot add another source configuration.
        r = self.replica([(0, 'auto', 0, 10, 10), (10, 'auto', 0, 10, 1),
                          (10, 'auto', 0, 8, 1), (5, 'auto', 0, 10, 1),
                          (5, 'auto', 0, 8, 1), (3, 'auto', 0, 10, 10),
                          (6, 'auto', 0, 8, 7)])
        r.flux()
        self.assertEqual(r.conditions_a(), ['associations == 0', 'separation >= 9.0'])
        bank = [r.state['banks'][key]['observation'] for key in r.state['flux']['bank']]
        self.assertEqual([row['step'] for row in bank], [20, 39])
        self.assertEqual(r.state['flux']['exposure_steps'], 39)
        self.assertFalse(self.script.events)

    def test_nearby_nontarget_contact_does_not_rearm_source_as_separated(self):
        self.use_toehold_distance()
        # At step 15 the toehold is far away but another switch nucleotide is
        # near the trigger. This is still the same excursion, not a return A.
        r = self.replica([(0, 'auto', 0, 10, 10), (10, 'auto', 0, 8, 7),
                          (5, 'auto', 0, 10, 3), (5, 'auto', 0, 8, 3),
                          (3, 'auto', 0, 10, 10), (6, 'auto', 0, 8, 7)])
        r.flux()
        self.assertEqual(r.conditions_a(), ['associations == 0', 'separation >= 9.0'])
        self.assertEqual(r.state['flux']['exposure_steps'], 29)
        bank = [r.state['banks'][key]['observation'] for key in r.state['flux']['bank']]
        self.assertEqual([row['step'] for row in bank], [10, 29])
        self.assertEqual([row['separation'] for row in bank], [7, 7])
        self.assertFalse(self.script.events)

    def test_old_plan_uses_original_distance_basin_and_new_plan_preserves_binding_union(self):
        old = self.replica([])
        self.assertEqual(old.conditions_a(), ['associations == 0', 'distance >= 9.0'])
        self.assertEqual(first_condition(Observation(0, 0, 0, 10., 0, 0), [old.conditions_a()]), 0)
        self.use_toehold_distance()
        new = self.replica([])
        nearby = Observation(0, 0, 0, 10., 0, 0, separation=3.)
        self.assertIsNone(first_condition(nearby, [new.conditions_a()]))
        self.assertEqual(first_condition(replace(nearby, separation=10.), [new.conditions_a()]), 0)
        # A real intended bond must cross the source even if a distance
        # observable has not crossed it, preserving alternative approach paths.
        bonded = Observation(1, 0, 0, 20., 0, 1, separation=.5)
        self.assertIsNotNone(first_condition(bonded, new.targets(0)))
        self.assertIsNone(first_condition(bonded, [new.conditions_a()]))

    def test_b_history_is_excluded_but_return_excursion_is_included(self):
        r = self.replica([(0,0,0,9),(10,0,0,8),(5,1,2,1),(10,None,2,1),(10,0,0,9),(7,0,0,8)])
        r.flux()
        self.assertEqual(r.state['flux']['exposure_steps'], 22)
        self.assertEqual(r.state['flux']['excluded_b_steps'], 20)
        self.assertEqual(r.state['steps'], 42)
        self.assertEqual(len(r.state['flux']['bank']), 2)

    def test_native_chunk_limit_is_continuation_not_failure(self):
        r = self.replica([(0,0,0,9),(10,None,0,8.5),(4,0,0,8),(3,0,0,9),(6,0,0,8)])
        r.flux()
        self.assertEqual(r.state['flux']['exposure_steps'], 23)
        self.assertEqual(len(r.state['flux']['bank']), 2)
        self.assertEqual(len(self.script.starts), 4)

    def test_burn_in_excluded_then_physical_return_to_a(self):
        self.req['kinetics']['equilibration_steps'] = 15
        r = self.replica([(10,None,0,10),(5,None,0,10),(0,0,0,10),
                          (4,0,0,8),(3,0,0,9),(6,0,0,8)])
        r.flux()
        self.assertEqual(r.state['burn_steps'], 15)
        self.assertEqual(r.state['flux']['exposure_steps'], 13)

    def test_pause_resume_during_return_preserves_history_and_count(self):
        r = self.replica([(0,0,0,9),(10,0,0,8)])
        def pause():
            if not self.script.events:
                (self.root/'CANCEL').touch()
        self.script.on_advance = pause
        with self.assertRaises(SamplingStopped):
            r.flux()
        r.save(force=True)
        saved = read_json(r.output/'state.json')
        self.assertEqual(saved['flux']['history'], 'return')
        self.assertEqual(saved['flux']['exposure_steps'], 10)
        (self.root/'CANCEL').unlink()
        other = NativeScript([(4,0,0,9),(6,0,0,8)])
        resumed = FFSReplica(self.req, r.output, 0, adapter_factory=other)
        resumed.flux()
        self.assertEqual(resumed.state['flux']['exposure_steps'], 20)
        self.assertEqual(len(resumed.state['flux']['bank']), 2)
        self.assertEqual(other.starts[0]['frame'], frame(10))
        self.assertFalse(other.starts[0]['initial'])

    def production_ready(self, events):
        r = self.replica(events)
        r.state['initial_velocities'] = False
        r.state['observation'] = {'binding': 1, 'stem': 0, 'opening': 0, 'distance': 1., 'step': 0, 'associations': 1}
        first = r.retain(interface=0)
        r.state['flux'].update(complete=True, bank=[first], exposure_steps=10)
        r.state['production_bank'] = [first]
        r.state['phase'] = 'production'
        # Small fixed protocol with two genuine successive interfaces.
        r.interfaces = [{'name': 'source', 'min_binding': 1}, {'name': 'duplex_2', 'min_binding': 2},
                        {'name': 'opening', 'min_binding': 2, 'max_opening': 0, 'max_stem': 0}]
        return r

    def test_fixed_shots_count_failure_and_keep_exact_parent_for_next_stage(self):
        # first shot success, second physical A return; then two successes.
        r = self.production_ready([(3,0,2,1),(4,2,0,9),(2,0,2,1),(0,0,2,1)])
        r.production()
        a,b = r.state['stages']
        self.assertEqual((a['attempted'],a['successes'],a['failures']), (2,1,1))
        self.assertEqual((b['attempted'],b['successes'],b['failures']), (2,2,0))
        self.assertEqual({t['parent'] for t in b['trials']}, set(a['bank']))
        self.assertEqual(self.script.starts[2]['frame'], frame(3))
        self.assertTrue((r.output/'movie.dat').is_file())
        self.assertTrue(r.state['completed'])

    def test_each_new_shot_publishes_selected_parent_observation_before_native_work(self):
        r = self.production_ready([(3, 0, 2, 1), (4, 2, 0, 9), (2, 0, 2, 1), (0, 0, 2, 1)])
        r.state['observation'] = {'binding': 0, 'stem': 1, 'opening': 1,
                                  'distance': 99., 'step': 999, 'associations': 0}
        observed = []

        def starting():
            entry = r.state['banks'][r.state['active']['parent']]
            self.assertIsNot(r.state['observation'], entry['observation'])
            observed.append((deepcopy(entry['observation']), deepcopy(r.state['observation']),
                             read_json(r.output/'progress.json')['observation']))

        self.script.on_start = starting
        r.production()
        self.assertEqual(len(observed), 4)
        for parent, current, published in observed:
            self.assertEqual(current, parent)
            self.assertEqual(published, parent)
        self.assertEqual(observed[0][1]['step'], 0)
        self.assertEqual(observed[2][1]['step'], 3)

    def test_zero_success_stops_replica_without_fabricating_next_stage(self):
        r = self.production_ready([(3,2,0,9),(4,2,0,9)])
        r.production()
        self.assertEqual(len(r.state['stages']), 1)
        self.assertEqual(r.state['stages'][0]['failures'], 2)
        self.assertTrue(r.state['completed'])
        self.assertEqual(r.state['zero_at'], 0)
        self.assertIsNone(r.state['movie'])

    def test_pause_in_pending_shot_does_not_turn_it_into_failure(self):
        r = self.production_ready([(5,None,1,1)])
        self.script.on_advance = lambda: (self.root/'CANCEL').touch()
        with self.assertRaises(SamplingStopped):
            r.production()
        r.save(force=True)
        stage = r.state['stages'][0]
        self.assertEqual((stage['attempted'],stage['pending'],stage['failures']), (1,1,0))
        self.assertEqual(r.state['active']['steps'],5)
        self.assertEqual(r.state['frame'], frame(5))

    def test_resuming_active_shot_keeps_its_latest_observation_not_its_parent(self):
        r = self.production_ready([(5, None, 1, .7)])
        self.script.on_advance = lambda: (self.root/'CANCEL').touch()
        with self.assertRaises(SamplingStopped):
            r.production()
        r.save(force=True)
        latest = deepcopy(r.state['observation'])
        self.assertEqual(latest['step'], 5)
        (self.root/'CANCEL').unlink()
        script = NativeScript([])
        resumed = FFSReplica(self.req, r.output, 0, adapter_factory=script)
        resumed.interfaces = r.interfaces

        def starting():
            self.assertEqual(resumed.state['observation'], latest)
            self.assertEqual(read_json(resumed.output/'progress.json')['observation'], latest)
            raise SamplingStopped('stop at native boundary')

        script.on_start = starting
        with self.assertRaises(SamplingStopped):
            resumed.production()
        self.assertEqual(script.starts[0]['frame'], frame(5))
        self.assertEqual(resumed.state['stages'][0]['attempted'], 1)
        self.assertEqual(resumed.state['stages'][0]['pending'], 1)

    def test_request_mismatch_rejects_resume(self):
        r = self.replica([])
        r.save(force=True)
        self.req['sampling']['salt_m'] = .2
        with self.assertRaisesRegex(ValueError,'different simulation'):
            self.replica([])

    def test_master_progress_aggregates_parallel_success_counts(self):
        sampler = KineticSampler(self.req, self.root)
        states = {i: {'replica_id':f'replica_{i+1}', 'index':i, 'phase':'production', 'completed':False,
                      'steps':10, 'flux':{'bank':['a','b'], 'complete':True},
                      'stages':[{'name':'duplex_2','successes':i+1,'failures':0,'pending':1}],
                      'observation':{'binding':1, 'stem':0,'opening':0,'distance':1.,'step':10}}
                  for i in range(2)}
        progress = sampler.aggregate(states,2)
        self.assertEqual(progress['kinetics']['workers'],2)
        self.assertEqual(progress['kinetics']['flux_crossings'],4)
        self.assertEqual(progress['kinetics']['flux_target'],8)
        self.assertEqual(progress['hybridization_progress']['production']['counts'][0],3)



if __name__ == '__main__':
    unittest.main()
