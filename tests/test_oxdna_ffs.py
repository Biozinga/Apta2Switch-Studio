"""PDB-free sequence/guide provenance and preparation/FFS lifecycle guards."""
from copy import deepcopy
from dataclasses import replace
from contextlib import contextmanager
from pathlib import Path
from types import SimpleNamespace
import tempfile
import json
import sys
import unittest
from unittest.mock import Mock, patch

from aptaswitch_core.architecture import architecture_layout, default_architecture_for_trigger
from aptaswitch_core.models import RunRequest, SequenceInput
from aptaswitch_core.oxdna_ffs import (build_ffs_plan, prepare_ffs_request, resolve_ffs_guides, ffs_visualization,
                                     run_oxdna_ffs)
from aptaswitch_studio.saved_switch_library import candidate_snapshot
from tests.design_fixtures import design_result_fixture
from tests.guide_fixtures import GuideBackend


class OxDNAFFSTests(unittest.TestCase):
    def setUp(self):
        self.root = Path(self.enterContext(tempfile.TemporaryDirectory()))
        self.backend = GuideBackend([])
        self.enterContext(patch('aptaswitch_core.thermodynamics.setup.get_backend', return_value=self.backend))
        request = RunRequest(sequences=SequenceInput('ACGT' * 3, 'ACGT'*7+'AC'),
                             architecture=default_architecture_for_trigger(30), trials=1)
        result = design_result_fixture(request)
        self.snapshot = candidate_snapshot(result, result.candidates[0], 'FFS candidate')
        self.guides = {'switch': self.snapshot['candidate']['structure_off'],
                       'trigger': '.' * len(self.snapshot['candidate']['trigger_seq'])}

    def request(self, options=None, *, guides=True):
        return prepare_ffs_request(self.snapshot, {'python':sys.executable}, options or {}, guides=self.guides if guides else None)

    def test_binding_register_and_opening_regions_follow_saved_architecture(self):
        plan = build_ffs_plan(self.snapshot)
        n = len(self.snapshot['candidate']['switch_seq'])
        arch = default_architecture_for_trigger(30)
        layout = architecture_layout(arch)
        self.assertEqual(sorted(plan['binding_pairs']), sorted([[layout['len_t7'] + 29 - i, n + i] for i in range(30)]))
        self.assertEqual(len(plan['stem_pairs']), arch.stem_length + arch.upper_stem_length + arch.upper_stem2_length)
        self.assertEqual(plan['opening_region'], list(range(layout['rbs_prefix_start'], layout['aug_start'] + 3)))
        self.assertEqual(plan['interfaces'][-1], {'name':'opening','min_binding':30,'max_stem':0,'max_opening':0})
        self.assertTrue(plan['rate_estimation'])
        self.assertEqual(len(plan['interfaces']),1)
        self.assertEqual(plan['distance_scope'], 'trigger_3prime_to_toehold_start')
        self.assertEqual(plan['distance_pairs'], [[layout['len_t7'], n + 29]])
        self.assertEqual(plan['separation_pairs'], [[i, j] for i in range(n) for j in range(n, n + 30)])

    def test_approach_skips_actual_leader_length_and_targets_antiparallel_terminal_base(self):
        from aptaswitch_core.oxdna_ffs_worker import validate_request
        from aptaswitch_core.oxdna_kinetics import source_fingerprint
        for leader in ('', 'GGG', 'GGGAGA'):
            with self.subTest(leader=leader):
                arch = replace(default_architecture_for_trigger(30), t7_leader=leader)
                request = RunRequest(sequences=SequenceInput('ACGT' * 3, 'ACGT'*7+'AC'),
                                     architecture=arch, trials=1)
                result = design_result_fixture(request)
                snapshot = candidate_snapshot(result, result.candidates[0], 'Terminal contact')
                prepared = prepare_ffs_request(snapshot, {'python':sys.executable}, {'seed':123})
                plan = validate_request(prepared)['plan']
                pair = (len(leader), plan['nucleotide_count'] - 1)
                self.assertEqual(plan['distance_pairs'], (pair,))
                self.assertIn(pair, plan['toehold_pairs'])
                self.assertIn(pair, plan['binding_pairs'])
                original = deepcopy(prepared)
                original['plan'].update(distance_scope='trigger_to_toehold',
                    distance_pairs=[[i, j] for i in sorted({p[0] for p in plan['toehold_pairs']})
                                    for j in range(plan['switch_length'], plan['nucleotide_count'])])
                self.assertNotEqual(source_fingerprint(prepared), source_fingerprint(original))

    def test_one_full_run_defaults_to_ten_replicas_without_mutating_candidate(self):
        original = deepcopy(self.snapshot)
        request = self.request({'seed':123})
        self.assertEqual(request['kinetics']['mode'], 'opening')
        self.assertEqual(request['kinetics']['replicas'], 10)
        self.assertEqual(request['kinetics']['source_samples'], 75)
        self.assertEqual(request['kinetics']['trials'], 250)
        self.assertLessEqual(request['kinetics']['workers'], 10)
        self.assertEqual(request['purpose'], 'native_direct_ffs_flux_and_shooting')
        self.assertEqual(self.snapshot, original)
        layout = architecture_layout(default_architecture_for_trigger(30))
        expected = [pair for pair in request['plan']['binding_pairs']
                    if layout['len_t7'] <= pair[0] < layout['len_t7'] + layout['len_toehold']]
        self.assertEqual(request['plan']['toehold_pairs'], expected)
        # Removed two-run choices cannot silently create a different experiment.
        for mode in ('encounter', 'opening', None):
            with self.subTest(mode=mode), self.assertRaises(ValueError):
                self.request({'simulation_mode':mode})

    def test_full_endpoint_matches_design_and_has_no_conflicting_pair_targets(self):
        from aptaswitch_core.preview import parse_structure
        from aptaswitch_core.architecture import target_switch_structures
        arch = default_architecture_for_trigger(30)
        layout = architecture_layout(arch)
        n = len(self.snapshot['candidate']['switch_seq'])
        _, target = target_switch_structures(arch, n-layout['reporter_start'])
        plan = self.request()['plan']
        # Target ON is trigger+switch; the MD topology is switch+trigger.
        remap = lambda i: n+i if i<30 else i-30
        on_pairs = {tuple(sorted((remap(i),remap(j)))) for i,j in parse_structure(target)[0]}
        self.assertEqual(on_pairs, {tuple(pair) for pair in plan['binding_pairs']})
        self.assertFalse(on_pairs & {tuple(pair) for pair in plan['stem_pairs']})
        self.assertFalse(on_pairs & {tuple(pair) for pair in plan['opening_pairs']})

    def test_final_artifacts_identify_the_full_simulation(self):
        with patch('aptaswitch_core.oxdna_guided.prepare_guided_system', side_effect=self.prepared), \
                patch('aptaswitch_core.oxdna_ffs.subprocess.Popen', side_effect=self.worker):
            result = run_oxdna_ffs(self.request(), self.root / 'kinetics')
        self.assertEqual(result['status'], 'completed')
        self.assertEqual(result['method'], 'NATIVE_FFS')
        self.assertEqual(result['kinetics']['mode'], 'opening')

    def test_proximity_keeps_the_same_guided_3d_preparation_and_movie_as_full_ffs(self):
        prepared_inputs = {}
        for mode in ('opening', 'proximity'):
            with self.subTest(mode=mode):
                request = self.request({'kinetics_mode': mode, 'seed': 41})
                folder = self.root / mode
                with patch('aptaswitch_core.oxdna_guided.prepare_guided_system', side_effect=self.prepared) as prepare, \
                        patch('aptaswitch_core.oxdna_ffs.subprocess.Popen', side_effect=self.worker) as worker:
                    result = run_oxdna_ffs(request, folder)
                prepare.assert_called_once()
                prepared_request = prepare.call_args.args[0]
                prepared_inputs[mode] = {key: prepared_request[key] for key in (
                    'switch_sequence', 'trigger_sequence', 'switch_guide', 'trigger_guide',
                    'guide_conditions', 'temperature_c', 'salt_m', 'initial_separation',
                )}
                self.assertEqual(prepare.call_args.args[1], folder.resolve() / 'preparation')
                self.assertEqual(result['status'], 'completed')
                self.assertTrue(result['preparation']['ffs_eligible'])
                self.assertEqual(result['movie_file'], 'preparation/validated.dat')
                self.assertEqual(result['movie_kind'], 'prepared')
                worker.assert_called_once()
                frozen = json.loads((folder / 'worker_request.json').read_text())
                self.assertEqual(frozen['prepared'], result['preparation'])
                if mode == 'proximity':
                    self.assertFalse(frozen['plan']['rate_estimation'])
                    self.assertEqual(frozen['kinetics']['source_samples'], 0)
        self.assertEqual(prepared_inputs['opening'], prepared_inputs['proximity'])

    def test_proximity_cannot_start_native_ffs_when_guided_3d_is_rejected(self):
        request = self.request({'kinetics_mode': 'proximity'})
        def reject(*args, **kwargs):
            return {**self.prepared(*args, **kwargs), 'guide_accepted': False,
                    'ffs_eligible': False, 'attempts': [1, 2, 3]}
        with patch('aptaswitch_core.oxdna_guided.prepare_guided_system', side_effect=reject), \
                patch('aptaswitch_core.oxdna_ffs.subprocess.Popen') as worker:
            result = run_oxdna_ffs(request, self.root / 'rejected_proximity')
        self.assertEqual(result['status'], 'failed')
        self.assertIn('3', result['error'])
        self.assertEqual(result['simulation_mode'], 'proximity')
        self.assertEqual(result['movie_file'], 'preparation/validated.dat')
        worker.assert_not_called()

    def test_new_proximity_request_equilibrates_monomers_without_source_flux(self):
        from aptaswitch_core.oxdna_kinetics import (
            OPENING_PROTOCOL_VERSION, PROXIMITY_PROTOCOL_VERSION, validate_kinetic_protocol,
        )
        request = self.request({'kinetics_mode': 'proximity'})
        self.assertEqual(request['kinetics']['version'], PROXIMITY_PROTOCOL_VERSION)
        self.assertEqual(request['kinetics']['equilibration_steps'], 200000)
        self.assertEqual(request['kinetics']['equilibration_method'],
                         'isolated_monomers_before_close_placement_v1')
        self.assertEqual(request['kinetics']['source_samples'], 0)
        self.assertFalse(request['plan']['rate_estimation'])
        validate_kinetic_protocol(request['kinetics'])
        opening = self.request()['kinetics']
        self.assertEqual(opening['version'], OPENING_PROTOCOL_VERSION)
        self.assertEqual(opening['equilibration_steps'], 200000)
        self.assertNotIn('equilibration_method', opening)

    def test_freeze_candidate_uses_native_coordinates_without_pdb_converter(self):
        data = self.request()
        self.assertEqual(data['version'], 2)
        self.assertEqual((data['switch_material'], data['trigger_material']), ('rna', 'dna'))
        self.assertEqual(data['trigger_guide'], self.guides['trigger'])
        self.assertNotIn('switch_pdb', data)
        self.assertNotIn('sources', data)
        self.assertFalse(data['use_average_seq'])
        self.assertFalse(data['magnesium_modeled'])
        self.assertNotIn('budget_seconds',data['sampling'])
        self.assertNotIn('max_steps_per_trial',data['sampling'])
        self.assertGreater(data['concentration_each_m'], 0)
        self.assertEqual(data['guide_conditions']['magnesium_m'], 0)


    def test_new_runs_use_fresh_seeds_while_explicit_seed_is_reproducible(self):
        with patch('aptaswitch_core.oxdna_ffs.secrets.randbelow', side_effect=[12, 91]):
            self.assertEqual(self.request()['sampling']['seed'], 13)
            self.assertEqual(self.request()['sampling']['seed'], 92)
        self.assertEqual(self.request({'seed': 123})['sampling']['seed'], 123)


    def test_lightweight_request_defers_folding_and_retains_candidate_identity(self):
        original = deepcopy(self.snapshot)
        with patch.object(self.backend, 'minimum') as folding:
            data = self.request(guides=False)
        folding.assert_not_called()
        self.assertNotIn('switch_guide', data)
        self.assertEqual(self.snapshot, original)
        self.snapshot['snapshot_id'] = ''
        with self.assertRaises(ValueError):
            self.request()

    def test_box_is_computed_from_design_strands_and_conditions_without_manual_override(self):
        from aptaswitch_core.oxdna_box import estimate_ffs_box
        data = self.request()
        expected = estimate_ffs_box(len(data['switch_sequence']), len(data['trigger_sequence']),
            data['switch_material'], data['trigger_material'], data['temperature_c'], data['salt_m'])
        self.assertEqual(data['box_size'], expected['box_size'])
        self.assertEqual(data['concentration_each_m'], expected['concentration_each_m'])
        self.assertEqual(data['box_policy']['estimated_box_size'], expected['box_size'])
        for value in (20, 100, None):
            with self.subTest(value=value), self.assertRaises(ValueError):
                self.request({'box_size': value})

    def test_native_guides_use_exact_strands_and_chosen_conditions(self):
        data = self.request({'temperature_c': 25., 'salt_m': .2}, guides=False)
        self.backend.records = iter((self.guides[role], -2.) for role in ('switch', 'trigger'))
        resolved = resolve_ffs_guides(data, backend=self.backend)
        self.assertEqual([(s.strands[0].sequence, s.strands[0].material) for s, _ in self.backend.calls],
                         [(data['switch_sequence'], 'rna'), (data['trigger_sequence'], 'dna')])
        for _, c in self.backend.calls:
            self.assertEqual((c.celsius, c.sodium, c.magnesium), (25., .2, 0.))
        self.assertEqual(resolved['switch_guide'], self.guides['switch'])
        self.assertEqual(resolved['guide_source'], 'native_hybrid_mfe_isolated_strands')
        self.assertEqual(resolved['guide_thermodynamics']['model'], self.backend.provenance.model)
        self.assertNotIn('switch_guide', data)

    def test_changed_model_malformed_guide_and_noncomplementary_trigger_rejected(self):
        data = self.request(guides=False)
        data['guide_model'] = 'another-model'
        with self.assertRaises(ValueError):
            resolve_ffs_guides(data, backend=self.backend)
        self.assertEqual(self.backend.calls, [])
        for guide in ('..', ')' * len(self.guides['trigger']), '[' * len(self.guides['trigger'])):
            with self.subTest(guide=guide), self.assertRaises(ValueError):
                prepare_ffs_request(self.snapshot, {}, {}, guides={**self.guides, 'trigger': guide})
        self.snapshot['candidate']['trigger_seq'] = 'A' * 30
        with self.assertRaises(ValueError):
            build_ffs_plan(self.snapshot)

    def test_saved_legacy_guides_do_not_load_or_recompute_thermodynamics(self):
        data = self.request()
        data.update(guide_model='external-archive', guide_source='archived_mfe_isolated_strands',
                    guide_thermodynamics={}, thermodynamics_library_path='/missing')
        original = deepcopy(data)
        with patch('aptaswitch_core.thermodynamics.setup.get_backend', side_effect=AssertionError('must not load')):
            resolved = resolve_ffs_guides(data)
        self.assertEqual(resolved, original)
        self.assertEqual(data, original)

    def test_parameters_are_finite_bounded_and_integer_where_needed(self):
        for options in ({'budget_minutes': float('nan')}, {'trials': 2.5}, {'flux_crossings': True},
                        {'max_steps': 10}, {'box_size': -1}, {'salt_m': 0}, {'salt_m': .02},
                        {'seed': -1}, {'prepare_only': 'yes'}, {'switch_target': 'other'}):
            with self.subTest(options=options), self.assertRaises(ValueError):
                self.request(options)

    def test_md_uses_only_design_fragment_and_rejects_full_reporter_options(self):
        original = deepcopy(self.snapshot)
        prefix = self.snapshot['request']['reporter']['sequence_after_start_rna']
        data = self.request(guides=False)
        self.assertEqual(data['switch_sequence'], original['candidate']['switch_seq'])
        self.assertEqual(data['plan']['switch_length'], len(original['candidate']['switch_seq']))
        self.assertEqual(min(j for i, j in data['plan']['binding_pairs']), data['plan']['switch_length'])
        self.assertEqual(self.snapshot, original)
        for options in ({'switch_target': 'full_transcript'},
                        {'full_reporter_sequence_after_start': prefix + 'ACG' * 100}):
            with self.subTest(options=options), self.assertRaises(ValueError):
                self.request(options, guides=False)
        extended = deepcopy(self.snapshot)
        extended['candidate']['switch_seq'] += 'ACG' * 100
        with self.assertRaises(ValueError):
            prepare_ffs_request(extended, {}, {})

    def prepared(self, request, output, **kwargs):
        output.mkdir(parents=True)
        configuration = output / 'validated.dat'
        configuration.write_text('validated by the preparation fixture')
        return {'status': 'completed', 'configuration_file': str(configuration),
                'topology_file': str(output / 'system.top'), 'ffs_eligible': True}


    def worker(self, command, **kwargs):
        folder=Path(kwargs['cwd'])
        (folder/'result.json').write_text(json.dumps({'status':'completed','method':'NATIVE_FFS'}))
        return SimpleNamespace(poll=lambda:0,returncode=0)

    def test_preparation_box_adjustment_updates_reported_volume_and_provenance(self):
        import json
        data = self.request()
        initial = data['box_size']
        def expanded(request, output, **kwargs):
            return {**self.prepared(request, output, **kwargs), 'box_size': initial+5}
        with patch('aptaswitch_core.oxdna_guided.prepare_guided_system', side_effect=expanded), \
             patch('aptaswitch_core.oxdna_ffs.subprocess.Popen',side_effect=self.worker):
            result = run_oxdna_ffs(data, self.root / 'expanded')
        self.assertEqual(result['status'], 'completed')
        self.assertEqual(result['box_size'], initial+5)
        self.assertTrue(result['box_policy']['placement_adjusted'])
        self.assertEqual(result['box_policy']['estimated_box_size'], initial)
        self.assertAlmostEqual(result['concentration_each_m'], data['concentration_each_m'] * (initial/(initial+5))**3)
        resolved = json.loads((self.root/'expanded/resolved_request.json').read_text())
        self.assertEqual(resolved['box_size'], result['box_size'])
        self.assertEqual(resolved['concentration_each_m'], result['concentration_each_m'])
        self.assertEqual(data['box_size'], initial)  # Original frozen request stays intact.

    def test_preparation_progress_uses_replica_target_not_number_of_interfaces(self):
        seen = []
        data = self.request({'kinetics_replicas':6})
        def prepare(*args,**kwargs):
            kwargs['progress_callback']('Preparation fixture')
            return self.prepared(*args,**kwargs)
        with patch('aptaswitch_core.oxdna_guided.prepare_guided_system',side_effect=prepare), \
             patch('aptaswitch_core.oxdna_ffs.subprocess.Popen',side_effect=self.worker):
            result=run_oxdna_ffs(data,self.root/'progress',progress_callback=lambda *values:seen.append(values))
        self.assertEqual(result['status'],'completed')
        self.assertEqual(seen,[(0,6,'Preparation fixture')])

    def test_unlimited_preparation_has_no_hidden_deadline_after_hours(self):
        now = [0.]
        def prepare(*args, **kwargs):
            self.assertIsNone(kwargs['deadline'])
            now[0] = 36000.
            return self.prepared(*args, **kwargs)
        with patch('aptaswitch_core.oxdna_ffs.time.monotonic', side_effect=lambda: now[0]), \
             patch('aptaswitch_core.oxdna_guided.prepare_guided_system', side_effect=prepare), \
             patch('aptaswitch_core.oxdna_ffs.subprocess.Popen',side_effect=self.worker):
            result = run_oxdna_ffs(self.request(), self.root / 'long')
        self.assertEqual(result['status'], 'completed', result)
        self.assertNotIn('budget_seconds',result['sampling'])

    def test_colors_follow_2d_domains_and_trigger_extensions(self):
        from aptaswitch_core.preview import build_switch_model, TRIGGER_COLOR, EXTENSION_COLOR, EXTENSION_3PRIME_COLOR
        from aptaswitch_core.models import SwitchArchitecture
        data = ffs_visualization(self.snapshot)
        request = self.snapshot['request']
        model = build_switch_model(SwitchArchitecture(**request['architecture']), self.snapshot['candidate']['trigger_seq'],
                                   request['reporter']['sequence_after_start_rna'], request['reporter']['label'])
        n = len(self.snapshot['candidate']['switch_seq'])
        self.assertEqual(data['colors'][:n], model.switch_colors)
        self.assertEqual(data['colors'][n:], [TRIGGER_COLOR] * 30)
        extended = deepcopy(self.snapshot)
        extended['request']['sequences']['trigger_dna'] = extended['candidate']['trigger_seq'][3:-3]
        extended['request']['inverse_extension'] = {'extension_side': 'both'}
        colors = ffs_visualization(extended)['colors'][n:]
        self.assertEqual(colors, [EXTENSION_COLOR] * 3 + [TRIGGER_COLOR] * 24 + [EXTENSION_3PRIME_COLOR] * 3)

    def test_invalid_off_does_not_start_ffs_but_preserves_prepared_view(self):
        def prepare(*args, **kwargs):
            return {**self.prepared(*args, **kwargs), 'ffs_eligible': False}
        with patch('aptaswitch_core.oxdna_guided.prepare_guided_system', side_effect=prepare), \
             patch('aptaswitch_core.oxdna_ffs.subprocess.Popen') as worker:
            result = run_oxdna_ffs(self.request(), self.root / 'run')
        self.assertEqual(result['status'], 'failed')
        self.assertEqual(result['movie_kind'], 'prepared')
        self.assertTrue(result['prepared_movie_file'])
        worker.assert_not_called()

    def test_preparation_timeout_and_cancellation_are_terminal(self):
        for error, cancelled, expected in [(RuntimeError('preparation failed'), False, 'failed'),
                                           (RuntimeError('cancelled'), True, 'cancelled')]:
            with self.subTest(expected=expected), patch('aptaswitch_core.oxdna_guided.prepare_guided_system', side_effect=error):
                result = run_oxdna_ffs(self.request(), self.root / expected, cancel_check=lambda: cancelled)
            self.assertEqual(result['status'], expected)
            self.assertFalse(result['movie_file'])
            self.assertIsNone(result['rate_estimate'])

    def test_resume_uses_frozen_box_guides_and_preparation_without_a_current_runtime(self):
        import sys
        data = self.request()
        folder = self.root / 'resume'
        folder.mkdir()
        prepared = self.prepared(data, folder / 'preparation')
        Path(prepared['topology_file']).write_text('topology')
        data.update(runtime={'python': sys.executable}, prepared=prepared)
        # This is the saved experiment, even if automatic sizing changed since.
        data['box_size'] = 100.
        data.pop('box_policy')
        (folder / 'worker_request.json').write_text(json.dumps(data))
        (folder / 'CANCEL').touch()
        (folder / 'bank_00000.dat').write_text('preserved bank')
        (folder / 'result.json').write_text(json.dumps({'status': 'cancelled', 'elapsed_seconds': 123.}))

        def launch(command, **kwargs):
            self.assertFalse((folder / 'CANCEL').exists())
            sent = json.loads((folder / 'worker_request.json').read_text())
            self.assertTrue(sent['resume'])
            self.assertEqual(sent['box_size'], 100.)
            self.assertEqual(sent['switch_guide'], data['switch_guide'])
            (folder / 'result.json').write_text(json.dumps({'status': 'completed', 'steps': 1234}))
            return SimpleNamespace(poll=lambda: 0, returncode=0)

        with patch('aptaswitch_core.oxdna_ffs.resolve_ffs_guides') as guides, \
             patch('aptaswitch_core.oxdna_guided.prepare_guided_system') as prepare, \
             patch('aptaswitch_core.oxdna_ffs.subprocess.Popen', side_effect=launch):
            result = run_oxdna_ffs({**data, 'resume': True}, folder)
        guides.assert_not_called()
        prepare.assert_not_called()
        self.assertEqual(result['status'], 'completed', result)
        self.assertGreaterEqual(result['elapsed_seconds'], 123.)
        self.assertEqual((folder / 'bank_00000.dat').read_text(), 'preserved bank')

    def test_resume_rebases_only_preparation_files_after_folder_move(self):
        from aptaswitch_core.oxdna_ffs import load_ffs_resume_request
        folder = self.root / 'before'
        folder.mkdir()
        data = self.request()
        prepared = self.prepared(data, folder / 'preparation')
        Path(prepared['topology_file']).write_text('topology')
        (folder / 'worker_request.json').write_text(json.dumps({**data, 'prepared': prepared}))
        moved = self.root / 'after'
        folder.rename(moved)
        request = load_ffs_resume_request(moved)
        self.assertTrue(Path(request['prepared']['configuration_file']).is_relative_to(moved.resolve()))
        (moved / 'preparation' / 'system.top').unlink()
        (moved / 'preparation' / 'system.top').symlink_to(moved / 'preparation' / 'validated.dat')
        with self.assertRaises(ValueError):
            load_ffs_resume_request(moved)

    def test_resume_rebases_selected_attempt_without_losing_root_topology(self):
        from aptaswitch_core.oxdna_ffs import load_ffs_resume_request
        folder = self.root / 'before'
        folder.mkdir()
        data = self.request()
        prepared = self.prepared(data, folder / 'preparation')
        Path(prepared['topology_file']).write_text('topology')
        attempt = folder / 'preparation' / 'attempt_02'
        attempt.mkdir()
        final = attempt / 'validated.dat'
        Path(prepared['configuration_file']).rename(final)
        prepared['configuration_file'] = str(final)
        (folder / 'worker_request.json').write_text(json.dumps({**data, 'prepared': prepared}))
        moved = self.root / 'after'
        folder.rename(moved)
        request = load_ffs_resume_request(moved)
        self.assertEqual(Path(request['prepared']['configuration_file']), moved.resolve() / 'preparation' / 'attempt_02' / 'validated.dat')
        self.assertEqual(Path(request['prepared']['topology_file']), moved.resolve() / 'preparation' / 'system.top')


if __name__ == '__main__':
    unittest.main()
