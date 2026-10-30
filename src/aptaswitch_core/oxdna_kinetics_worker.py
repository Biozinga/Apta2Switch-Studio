"""Native direct FFS: continuous source flux, real crossing banks and shooting.

Published protocol: https://lorenzo-rovigatti.github.io/oxDNA/ffs.html .
Native FFS_MD owns every stopping condition. Replicas run in isolated processes;
finite integration chunks are checkpoints, never physical failure boundaries.
"""
from __future__ import annotations

from dataclasses import asdict
import copy
import gc
import json
import math
import os
from pathlib import Path
import random
import shutil
import subprocess
import sys
import time

try:
    from .oxdna_ffs_worker import (FrameBuffer, SamplingStopped, atomic_json, validate_request,
                                    request_identity, _tuples, _same_checkpoint)
    from .oxdna_initialization import (random_separated_configuration, random_proximity_configuration,
                                       split_monomer_frames, combine_monomer_frames)
    from .oxdna_native_ffs import NativeFFSAdapter
    from .oxdna_box import check_periodic_self_clearance
    from .oxdna_kinetics import (candidate_interfaces, progress_index, condition, source_fingerprint,
                                  direct_ffs_report, proximity_ffs_report, validate_kinetic_protocol,
                                  PROXIMITY_PROTOCOL_VERSION)
except ImportError:
    from oxdna_ffs_worker import (FrameBuffer, SamplingStopped, atomic_json, validate_request,
                                  request_identity, _tuples, _same_checkpoint)
    from oxdna_initialization import (random_separated_configuration, random_proximity_configuration,
                                      split_monomer_frames, combine_monomer_frames)
    from oxdna_native_ffs import NativeFFSAdapter
    from oxdna_box import check_periodic_self_clearance
    from oxdna_kinetics import (candidate_interfaces, progress_index, condition, source_fingerprint,
                                direct_ffs_report, proximity_ffs_report, validate_kinetic_protocol,
                                PROXIMITY_PROTOCOL_VERSION)


def read_json(path):
    try:
        return json.loads(Path(path).read_text())
    except (OSError, ValueError):
        return {}


class IsolatedMonomerAdapter:
    """Strict, unrestrained NA MD of one monomer, with no partner present."""

    def __init__(self, request, configuration, seed, initial, *, topology_file, directory):
        import oxpy
        self.context = None
        self.manager = self.info = self.observable = None
        self.box_policy = request.get('box_policy')
        self.periodic_clearance = None
        directory = Path(directory)
        directory.mkdir(parents=True, exist_ok=True)
        sampling = request['sampling']
        values = {
            'backend': 'CPU', 'sim_type': 'MD', 'interaction_type': 'NA',
            'topology': str(Path(topology_file).resolve()), 'conf_file': str(Path(configuration).resolve()),
            'seed': str(seed), 'T': f"{sampling['temperature_c']}C", 'salt_concentration': str(sampling['salt_m']),
            'dt': str(sampling['dt']), 'thermostat': 'brownian',
            'newtonian_steps': str(sampling['newtonian_steps']), 'diff_coeff': str(sampling['diff_coeff']),
            'refresh_vel': 'true' if initial else 'false', 'restart_step_counter': 'false',
            'fix_diffusion': 'false', 'reset_initial_com_momentum': 'false',
            'verlet_skin': '0.2', 'steps': '0', 'time_scale': 'linear',
            'print_conf_interval': '1000000000', 'print_energy_every': '1000000000',
            'external_forces': 'false', 'use_average_seq': 'false',
            'trajectory_file': os.devnull, 'lastconf_file': os.devnull, 'energy_file': os.devnull,
            'log_file': str(directory / 'native_md.log'),
            'data_output_1': '{\nname = ' + os.devnull + '\nprint_every = 1000000000\nonly_last = true\ncol_1 = {\ntype = Configuration\nid = monomer_checkpoint\n}\n}',
        }
        for field in ('seq_dep_file_DNA', 'seq_dep_file_RNA', 'seq_dep_file_NA'):
            values[field] = str(Path(request['runtime'][field]).resolve(strict=True))
        (directory / 'native_input.txt').write_text('\n'.join(f'{key} = {value}' for key, value in values.items()) + '\n')
        self.context = oxpy.Context()
        self.context.__enter__()
        try:
            settings = oxpy.InputFile()
            for key, value in values.items():
                settings[key] = value
            self.manager = oxpy.OxpyManager(settings)
            self.info = self.manager.config_info()
            self.observable = self.info.get_observable_by_id('monomer_checkpoint')
            self.nucleotide_count = len(self.info.particles())
        except BaseException:
            self.close()
            raise

    def advance(self, steps):
        before = self.manager.current_step
        self.manager.run(steps, print_output=False)
        actual = self.manager.current_step - before
        if actual != steps:
            raise RuntimeError('Native isolated-monomer MD stopped before its requested burn-in chunk.')
        return actual

    def frame(self):
        self.manager.update_observable_data()
        frame = self.observable.get_output_string(self.manager.current_step).rstrip() + '\n'
        lines = frame.splitlines()
        rows = [[float(value) for value in line.split()] for line in lines[3:] if line.strip()]
        if (len(rows) != self.nucleotide_count or any(len(row) != 15 for row in rows)
                or any(not math.isfinite(value) for row in rows for value in row)):
            raise ValueError('Native isolated-monomer equilibration produced a nonfinite phase point.')
        if self.box_policy:
            box = [float(value) for value in lines[1].split('=', 1)[1].split()]
            self.periodic_clearance = check_periodic_self_clearance(
                [row[:3] for row in rows], [self.nucleotide_count], box, self.box_policy['interaction_cutoff'])
        return frame

    def close(self):
        self.observable = self.info = self.manager = None
        gc.collect()
        context, self.context = self.context, None
        if context is not None:
            context.__exit__(None, None, None)


class FFSReplica:
    """One independently randomized source trajectory and its full FFS lineage."""

    def __init__(self, request, output, replica_index, *, adapter_factory=NativeFFSAdapter,
                 monomer_adapter_factory=IsolatedMonomerAdapter, clock=time.monotonic):
        self.request = validate_request(request)
        self.options = self.request['kinetics']
        validate_kinetic_protocol(self.options)
        self.output = Path(output)
        self.output.mkdir(parents=True, exist_ok=True)
        (self.output/'banks').mkdir(exist_ok=True)
        self.root = self.output.parent
        self.index = replica_index
        self.adapter_factory = adapter_factory
        self.monomer_adapter_factory = monomer_adapter_factory
        self.clock = clock
        self.last_save = -float('inf')
        self.identity = request_identity(self.request)
        self.proximity = self.options['mode'] == 'proximity'
        self.equilibrated_proximity = self.proximity and self.options['version'] == PROXIMITY_PROTOCOL_VERSION
        self.interfaces = candidate_interfaces(self.request['plan'], mode=self.options['mode'])
        self.rng = random.Random(self.request['sampling']['seed'] + 104729 * (replica_index+1))
        self.state = read_json(self.output/'state.json')
        if self.state:
            if self.state.get('identity') != self.identity:
                raise ValueError('The saved FFS replica belongs to different simulation conditions.')
            self.rng.setstate(_tuples(self.state['random_state']))
        else:
            place = random_proximity_configuration if self.proximity and not self.equilibrated_proximity else random_separated_configuration
            initial_frame, initialization = place(self.request, self.rng)
            initialization.update(replica_index=replica_index,
                                  seed=self.request['sampling']['seed'] + 104729 * (replica_index+1))
            self.state = {'schema_version': 3, 'identity': self.identity,
                'replica_id': f'replica_{replica_index+1}', 'index': replica_index,
                'phase': 'equilibration', 'frame': initial_frame, 'initialization': initialization,
                'initial_velocities': True, 'steps': 0, 'burn_steps': 0, 'active': None,
                'observation': {}, 'banks': {}, 'stages': [], 'production_index': 0,
                'flux': {'history': 'initial_A', 'exposure_steps': 0, 'excluded_b_steps': 0,
                         'bank': [], 'complete': False}, 'production_bank': [], 'completed': False,
                'movie': None, 'source_fingerprint': source_fingerprint(self.request)}
            if self.proximity:
                self.state.update(phase='equilibration_switch' if self.equilibrated_proximity else 'initialization', initial_source=None)
            if self.equilibrated_proximity:
                frames = split_monomer_frames(initial_frame, self.request['plan']['switch_length'])
                self.state['monomer_equilibration'] = {
                    'method': self.options['equilibration_method'],
                    'target_steps_per_monomer': self.options['equilibration_steps'],
                    'completed': False, 'monomers': [
                        {'role': role, 'steps': 0, 'frame': frame, 'initial_velocities': True}
                        for role, frame in zip(('switch', 'trigger'), frames)]}

    def check(self):
        if (self.root/'CANCEL').exists():
            raise SamplingStopped('cancelled')

    def save(self, force=False):
        now = self.clock()
        if not force and now-self.last_save < 2:
            return
        self.state['random_state'] = self.rng.getstate()
        self.state['updated_at'] = time.time()
        atomic_json(self.output/'state.json', self.state)
        atomic_json(self.output/'progress.json', self.summary())
        (self.output/'live.dat.tmp').write_text(self.state['frame'])
        (self.output/'live.dat.tmp').replace(self.output/'live.dat')
        self.last_save = now

    def summary(self):
        value = {key: self.state[key] for key in ('replica_id', 'index', 'phase', 'steps', 'burn_steps',
                  'production_index', 'completed', 'observation', 'movie')}
        value['initialization'] = self.state.get('initialization')
        if self.equilibrated_proximity:
            equilibration = self.state['monomer_equilibration']
            steps = {row['role']: row['steps'] for row in equilibration['monomers']}
            value['monomer_equilibration'] = {
                'method': equilibration['method'], 'target_steps_each': equilibration['target_steps_per_monomer'],
                'total_steps': 2 * equilibration['target_steps_per_monomer'],
                'completed_steps': sum(steps.values()), 'switch_steps': steps['switch'], 'trigger_steps': steps['trigger'],
                'completed': equilibration['completed'], 'equilibrium_validated': False,
            }
        if self.proximity:
            value['initial_source'] = self.state.get('initial_source')
        value['flux'] = {key: item for key, item in self.state['flux'].items() if key != 'bank'}
        value['flux']['crossings'] = len(self.state['flux']['bank'])
        value['stages'] = [{key: item for key, item in stage.items() if key not in ('bank', 'trials')}
                           for stage in self.state['stages']]
        return value

    def conditions_a(self):
        distance = 'separation' if self.request['plan'].get('separation_pairs') else 'distance'
        return ['associations == 0', f"{distance} >= {self.request['plan']['basin_distance']}"]

    def conditions_b(self):
        return condition(self.request['plan']['interfaces'][-1])

    def targets(self, index):
        # Union with all later predicates makes interfaces nested even if
        # geometry and HB-count coordinates disagree or an MD step overshoots.
        return [condition(item) for item in self.interfaces[index:]]

    def drive(self, conditions, *, exposure=None, max_steps=None, record=False):
        """Continue the same physical phase point until a native condition fires."""
        self.check()
        path = self.output/'active.dat'
        path.write_text(self.state['frame'])
        initial = self.state['initial_velocities']
        adapter = self.adapter_factory(self.request, path, self.rng.randrange(1, 2**31), initial,
                                       directory=self.output/'native', conditions=conditions)
        self.state['initial_velocities'] = False
        advanced = 0
        buffer = None
        if record:
            active = self.state['active']
            buffer = FrameBuffer(active['frame_interval'], 32)
            buffer.frames = [tuple(row) for row in active['frames']]
            if buffer.frames:
                buffer.next_step = buffer.frames[-1][0] + buffer.interval
        try:
            while True:
                self.check()
                chunk = min(self.options['chunk_steps'], max_steps-advanced) if max_steps is not None else self.options['chunk_steps']
                if chunk <= 0:
                    return None
                result = adapter.advance(chunk)
                advanced += result.steps
                self.state['steps'] += result.steps
                self.state['frame'] = adapter.frame()
                self.state['observation'] = asdict(result.observation)
                if exposure == 'A':
                    self.state['flux']['exposure_steps'] += result.steps
                elif exposure == 'B':
                    self.state['flux']['excluded_b_steps'] += result.steps
                elif exposure == 'burn':
                    self.state['burn_steps'] += result.steps
                if record:
                    active['steps'] += result.steps
                    buffer.add(active['steps'], lambda: self.state['frame'], force=result.condition is not None)
                    active['frames'] = buffer.frames
                    active['frame_interval'] = buffer.interval
                    if self.proximity:
                        observations = active.setdefault('observations', {})
                        observations[str(active['steps'])] = asdict(result.observation)
                        active['observations'] = {str(step): observations[str(step)] for step, _ in buffer.frames
                                                  if str(step) in observations}
                self.save()
                if result.condition is not None:
                    return result.condition
                if max_steps is not None and advanced >= max_steps:
                    return None
                if result.steps == 0:
                    raise RuntimeError('Native FFS made no progress and returned no stopping condition.')
        finally:
            adapter.close()

    def retain(self, parent=None, frames=None, interface=0, observations=None):
        key = f"state_{len(self.state['banks']):07d}"
        config = f'banks/{key}.dat'
        trajectory = f'banks/{key}.trajectory.dat'
        frame = self.state['frame']
        (self.output/config).write_text(frame)
        (self.output/trajectory).write_text(''.join(row[1] for row in frames) if frames else frame)
        self.state['banks'][key] = {'id': key, 'parent': parent, 'configuration_file': config,
            'trajectory_file': trajectory, 'interface_index': interface,
            'observation': copy.deepcopy(self.state['observation'])}
        if self.proximity:
            # Keep the controller checkpoint small as the branching bank grows.
            observables_file = f'banks/{key}.observables.json'
            values = (
                [copy.deepcopy(observations[str(step)]) for step, _ in frames]
                if frames and observations else [copy.deepcopy(self.state['observation'])])
            atomic_json(self.output/observables_file, values)
            self.state['banks'][key]['observables_file'] = observables_file
        return key

    def flux(self):
        if self.proximity:
            raise ValueError('The proximity protocol has no source-flux phase.')
        flux = self.state['flux']
        if flux['complete']:
            return
        if self.state['burn_steps'] < self.options['equilibration_steps']:
            self.state['phase'] = 'equilibration'
            # Unattainable native stopping threshold: integrate a fixed burn-in.
            self.drive([[f"binding > {len(self.request['plan']['binding_pairs'])}"]], exposure='burn',
                       max_steps=self.options['equilibration_steps']-self.state['burn_steps'])
        self.state['phase'] = 'flux'
        while len(flux['bank']) < self.options['source_samples']:
            history = flux['history']
            if history == 'initial_A':
                self.drive([self.conditions_a()])
                flux['history'] = 'forward'
            elif history == 'forward':
                self.drive(self.targets(0), exposure='A')
                flux['bank'].append(self.retain(interface=0))
                flux['history'] = 'return'
            elif history == 'return':
                hit = self.drive([self.conditions_a(), self.conditions_b()], exposure='A')
                flux['history'] = 'forward' if hit == 0 else 'B_history'
            elif history == 'B_history':
                # No physical reset/reinjection. Time whose most recent basin
                # is B must be excluded from the effective A-source flux.
                self.drive([self.conditions_a()], exposure='B')
                flux['history'] = 'forward'
            else:
                raise ValueError('Unknown source-flux history in checkpoint.')
            self.save(force=True)
        if flux['exposure_steps'] <= 0:
            raise ValueError('A source flux requires positive physical exposure time.')
        flux['complete'] = True
        self.state['production_bank'] = flux['bank'][:]
        self.state['phase'] = 'production'
        self.save(force=True)

    def equilibrate_proximity_monomers(self):
        """Burn in each isolated monomer; partner interactions are absent."""
        if not self.equilibrated_proximity:
            return
        equilibration = self.state['monomer_equilibration']
        target = equilibration['target_steps_per_monomer']
        directory = self.output / 'monomer_equilibration'
        directory.mkdir(exist_ok=True)
        for monomer in equilibration['monomers']:
            if monomer['steps'] >= target:
                continue
            self.check()
            role = monomer['role']
            self.state['phase'] = 'equilibration_' + role
            self.save(force=True)
            topology, configuration = directory / (role + '.top'), directory / (role + '.dat')
            sequence, material = self.request[role + '_sequence'], self.request[role + '_material']
            topology.write_text(f"{len(sequence)} 1 5->3\n{sequence} type={material.upper()}\n")
            configuration.write_text(monomer['frame'])
            adapter = self.monomer_adapter_factory(
                self.request, configuration, self.rng.randrange(1, 2**31), monomer['initial_velocities'],
                topology_file=topology, directory=directory / role)
            try:
                # Persist thermally drawn momenta even if interrupted before
                # the first integration chunk; a restart never redraws them.
                monomer['frame'] = adapter.frame()
                monomer['initial_velocities'] = False
                self.state['frame'] = combine_monomer_frames([part['frame'] for part in equilibration['monomers']])
                self.save(force=True)
                while monomer['steps'] < target:
                    self.check()
                    chunk = min(self.options['chunk_steps'], target - monomer['steps'])
                    actual = adapter.advance(chunk)
                    if type(actual) is not int or actual != chunk:
                        raise RuntimeError('Isolated monomer burn-in did not complete its requested native steps.')
                    monomer['frame'] = adapter.frame()
                    monomer['steps'] += actual
                    self.state['burn_steps'] += actual
                    self.state['steps'] += actual
                    self.state['frame'] = combine_monomer_frames([part['frame'] for part in equilibration['monomers']])
                    self.save(force=True)
            finally:
                adapter.close()
        equilibration['completed'] = True
        self.state['phase'] = 'initialization'
        self.save(force=True)

    def prepare_proximity(self):
        """Accept a native-unbound near source after the frozen preparation.

        New sources equilibrate isolated monomers first, then rotate/translate
        their entire phase points close to the toehold. Legacy protocol 3 keeps
        its original no-burn-in source and draws thermal velocities once.
        """
        if not self.proximity:
            raise ValueError('Only the proximity protocol accepts an artificial source.')
        if self.state.get('initial_source'):
            return
        if self.equilibrated_proximity:
            self.equilibrate_proximity_monomers()
            if not self.state.get('proximity_placement'):
                self.check()
                source = self.output / 'equilibrated_monomers.dat'
                source.write_text(combine_monomer_frames(
                    [row['frame'] for row in self.state['monomer_equilibration']['monomers']]))
                proposal = copy.deepcopy(self.request)
                proposal['prepared']['configuration_file'] = str(source)
                frame, initialization = random_proximity_configuration(proposal, self.rng)
                initialization.update(
                    monomer_equilibration=self.summary()['monomer_equilibration'],
                    separated_initialization=self.state['initialization'],
                    replica_index=self.index,
                    seed=self.request['sampling']['seed'] + 104729 * (self.index+1))
                self.state.update(frame=frame, initialization=initialization, proximity_placement=True,
                                  initial_velocities=False)
                self.save(force=True)
        self.check()
        path = self.output/'active.dat'
        path.write_text(self.state['frame'])
        adapter = self.adapter_factory(self.request, path, self.rng.randrange(1, 2**31), not self.equilibrated_proximity,
                                       directory=self.output/'native', conditions=[])
        try:
            observation = adapter.observe()
            minimum, maximum = self.options['initial_distance_min'], self.options['initial_distance_max']
            if (observation.associations != 0 or observation.binding != 0
                    or observation.separation is None or observation.separation < minimum-1e-8
                    or not minimum-1e-8 <= observation.distance <= maximum+1e-8):
                raise ValueError('The native proximity seed must be close to the toehold with zero interstrand bonds.')
            self.state['frame'] = adapter.frame()
            self.state['observation'] = asdict(observation)
        finally:
            adapter.close()
        self.state['initial_velocities'] = False
        self.state['initialization'].update(native_unbound_verified=True,
                                             native_interstrand_bonds_checked=True,
                                             native_interstrand_bonds=0)
        key = self.retain(interface=0)
        self.state.update(initial_source=key, production_bank=[key], phase='production')
        self.save(force=True)

    def production(self):
        # λ0 comes from the continuous source trajectory; shooting starts at λ1.
        stages = self.state['stages']
        while self.state['production_index'] < len(self.interfaces)-1:
            index = self.state['production_index']
            level = index+1
            if len(stages) == index:
                stages.append({'name': self.interfaces[level]['name'], 'target_trials': self.options['trials'],
                    'attempted': 0, 'successes': 0, 'failures': 0, 'pending': 0, 'censored': 0,
                    'completed': False, 'bank': [], 'trials': []})
            stage = stages[index]
            while stage['successes']+stage['failures'] < stage['target_trials']:
                self.check()
                if self.state['active'] is None:
                    parent = self.rng.choice(self.state['production_bank'])
                    entry = self.state['banks'][parent]
                    self.state['frame'] = (self.output/entry['configuration_file']).read_text()
                    self.state['observation'] = copy.deepcopy(entry['observation'])
                    self.state['active'] = {'parent': parent, 'steps': 0,
                        'frames': [[0, self.state['frame']]], 'frame_interval': self.request['sampling']['frame_interval']}
                    if self.proximity:
                        self.state['active']['observations'] = {'0': copy.deepcopy(entry['observation'])}
                    stage['attempted'] += 1
                    stage['pending'] = 1
                    self.save(force=True)
                targets = self.targets(level)
                hit = self.drive([*targets, self.conditions_a()], record=True)
                success = hit < len(targets)
                active = self.state['active']
                stage['successes' if success else 'failures'] += 1
                stage['pending'] = 0
                stage['trials'].append({'parent': active['parent'], 'steps': active['steps'],
                                        'outcome': 'success' if success else 'failure'})
                if success:
                    key = self.retain(active['parent'], active['frames'], level, active.get('observations'))
                    stage['bank'].append(key)
                    if level == len(self.interfaces)-1 and self.state['movie'] is None:
                        self.state['movie'] = key
                        self.movie()
                self.state['active'] = None
                self.save(force=True)
            stage['completed'] = True
            self.state['production_index'] += 1
            if not stage['successes']:
                self.state['zero_at'] = index
                break
            self.state['production_bank'] = stage['bank'][:]
            self.save(force=True)
        self.state['completed'] = True
        self.state['phase'] = 'completed'
        self.save(force=True)

    def movie(self):
        key = self.state['movie']
        lineage, frames, observations = [], [], []
        while key is not None:
            row = self.state['banks'][key]
            lineage.append(row)
            key = row['parent']
        for row in reversed(lineage):
            data = (self.output/row['trajectory_file']).read_text()
            segment = ['t ='+part for part in data.split('t =')[1:]]
            if frames and segment and not _same_checkpoint(frames[-1], segment[0]):
                raise ValueError('A saved FFS branch does not start from its physical parent.')
            if self.proximity:
                series = read_json(self.output/row['observables_file']) if row.get('observables_file') else []
                if len(series) != len(segment):
                    raise ValueError('The saved proximity observations do not match the trajectory frames.')
                observations.extend(series[1:] if frames else series)
            frames.extend(segment[1:] if frames else segment)
        cap = self.request['sampling']['max_movie_frames']
        if len(frames) > cap:
            indices = [round(i*(len(frames)-1)/(cap-1)) for i in range(cap)]
            frames = [frames[i] for i in indices]
            if self.proximity:
                observations = [observations[i] for i in indices]
        (self.output/'movie.dat').write_text(''.join(frames))
        atomic_json(self.output/'movie_metadata.json', {'trajectory_frames': len(frames),
            'lineage': [row['id'] for row in reversed(lineage)], 'trajectory_kind': 'first_successful_FFS_lineage',
            'trajectory_is_unbiased_waiting_time': False,
            'trajectory_endpoint': 'opening',
            'final_interface_reached': True})
        if self.proximity:
            metadata = read_json(self.output/'movie_metadata.json')
            metadata.update(simulation_mode='proximity', trajectory_start='prepared_unbound_proximity',
                            trajectory_observables=observations)
            if self.equilibrated_proximity:
                metadata.update(trajectory_start='equilibrated_monomers_then_unbound_proximity',
                                monomer_equilibration=self.summary()['monomer_equilibration'])
            atomic_json(self.output/'movie_metadata.json', metadata)

    def run(self):
        try:
            if not self.state['completed']:
                if self.proximity:
                    self.prepare_proximity()
                else:
                    self.flux()
                self.production()
            self.state['status'] = 'completed'
        except SamplingStopped:
            self.state['status'] = 'cancelled'
        except Exception as exc:
            self.state.update(status='failed', error=str(exc))
        finally:
            self.save(force=True)
            atomic_json(self.output/'result.json', {**self.summary(), 'status': self.state['status'],
                                                   'error': self.state.get('error')})
        return self.state


class KineticSampler:
    """Supervise isolated native replicas; navigation never owns their lifetime."""

    def __init__(self, request, output, *, started=None, clock=time.monotonic):
        self.request = validate_request(request)
        self.options = self.request['kinetics']
        validate_kinetic_protocol(self.options)
        self.output = Path(output)
        self.started = started or clock()
        self.clock = clock
        self.proximity = self.options['mode'] == 'proximity'
        self.interfaces = candidate_interfaces(self.request['plan'], mode=self.options['mode'])
        self.identity = request_identity(self.request)
        checkpoint = read_json(self.output/'resume_state.json')
        if checkpoint and checkpoint.get('identity') != self.identity:
            raise ValueError('This checkpoint has different native FFS conditions.')
        self.previous_elapsed = checkpoint.get('elapsed_seconds', 0.)
        atomic_json(self.output/'resume_state.json', {'schema_version': 3, 'identity': self.identity,
            'source_fingerprint': source_fingerprint(self.request), 'elapsed_seconds': self.previous_elapsed})

    def aggregate(self, states, active_workers, *, status='running'):
        rows = list(states.values())
        count = self.options['replicas']
        completed = sum(bool(row.get('completed')) for row in rows)
        representative = next((row for row in rows if not row.get('completed')), rows[-1] if rows else {})
        stages = []
        for level, interface in enumerate(self.interfaces):
            target = (1 if self.proximity else self.options['source_samples']) if level == 0 else self.options['trials']
            stage = {'name': interface['name'], 'target_trials': target*count, 'attempted': 0,
                     'successes': 0, 'failures': 0, 'pending': 0, 'completed': len(rows) == count}
            for row in rows:
                if level == 0:
                    if self.proximity:
                        ready = bool(row.get('initial_source'))
                        part = {'successes': int(ready), 'attempted': int(ready), 'completed': ready}
                    else:
                        flux = row.get('flux', {})
                        crossing_count = len(flux['bank']) if 'bank' in flux else flux.get('crossings', 0)
                        part = {'successes': crossing_count, 'attempted': crossing_count,
                                'completed': flux.get('complete', False)}
                else:
                    part = next((s for s in row.get('stages', []) if s['name'] == interface['name']), {})
                for key in ('attempted', 'successes', 'failures', 'pending'):
                    stage[key] += part.get(key, 0)
                stage['completed'] &= bool(part.get('completed'))
            stages.append(stage)
        thresholds = [{'name': item['name'], 'binding_pairs': item.get('min_binding') if item['name'] != 'opening' else None}
                      for item in self.interfaces if (item['name'].startswith('duplex_') and item['min_binding'] >= 2) or item['name'] == 'opening']
        counts = {stage['name']: stage['successes'] for stage in stages}
        progress = {'status': status, 'method': 'NATIVE_FFS', 'simulation_mode': self.options['mode'],
            'current': completed, 'total': count, 'steps': sum(r.get('steps', 0) for r in rows),
            'plan': {**self.request['plan'], 'interfaces': self.interfaces}, 'sampling': self.request['sampling'],
            'box_size': self.request['box_size'], 'box_policy': self.request.get('box_policy'),
            'concentration_each_m': self.request['concentration_each_m'],
            'stages': stages, 'observation': representative.get('observation', {}),
            'current_interface_index': progress_index(representative.get('observation', {}), self.interfaces),
            'kinetics': {'mode': self.options['mode'], 'phase': representative.get('phase', 'equilibration'),
                'replica': representative.get('index', 0)+1, 'replicas': count, 'workers': active_workers,
                'completed_replicas': completed, 'flux_crossings': stages[0]['successes'],
                'flux_target': stages[0]['target_trials'], 'flux_history': representative.get('flux', {}).get('history')},
            'hybridization_progress': {'thresholds': thresholds, 'production': {'counts': [counts[t['name']] for t in thresholds],
                'complete': True, 'completed_trials': sum(s['successes']+s['failures'] for s in stages[1:]),
                'active_trials': sum(s['pending'] for s in stages[1:])}},
            'elapsed_seconds': self.previous_elapsed+self.clock()-self.started}
        if self.request.get('protocol_amendments'):
            progress['protocol_amendments'] = self.request['protocol_amendments']
        phase = progress['kinetics']['phase']
        description = ('attente du retour naturel vers les brins libres' if progress['kinetics']['flux_history'] == 'B_history'
                       and phase == 'flux' else phase)
        progress['message'] = (f"FFS native · {description} · {active_workers} processus · "
                               f"{completed}/{count} répétitions · {stages[0]['successes']}/{stages[0]['target_trials']} contacts source")
        if self.proximity:
            for field in ('flux_crossings', 'flux_target', 'flux_history'):
                progress['kinetics'].pop(field, None)
            progress['kinetics'].update(initial_sources=stages[0]['successes'], initial_sources_target=count,
                                       phase=representative.get('phase', 'initialization'))
            progress['message'] = (f"FFS · départ proche du toehold · {progress['kinetics']['phase']} · "
                                   f"{active_workers} processus · {completed}/{count} répétitions")
            if self.options['version'] == PROXIMITY_PROTOCOL_VERSION:
                progress['monomer_equilibration'] = {
                    'method': self.options['equilibration_method'],
                    'target_steps_each': self.options['equilibration_steps'],
                    'total_steps': 2 * self.options['equilibration_steps'] * count,
                    'completed_steps': sum(row.get('burn_steps', 0) for row in rows),
                    'completed_monomers': sum(
                        sum(part['steps'] >= self.options['equilibration_steps'] for part in eq['monomers'])
                        if 'monomers' in eq else int(eq.get('switch_steps', 0) >= self.options['equilibration_steps'])
                        + int(eq.get('trigger_steps', 0) >= self.options['equilibration_steps'])
                        for eq in (row.get('monomer_equilibration', {}) for row in rows)),
                    'total_monomers': 2 * count, 'equilibrium_validated': False,
                }
        return progress

    def report(self, states, status):
        if self.proximity:
            report = proximity_ffs_report(list(states.values()),
                                          protocol_amendments=self.request.get('protocol_amendments'))
            if (not self.request.get('protocol_amendments')
                    and (len(states) != self.options['replicas'] or status != 'completed')):
                report.update(status='incomplete', value=None, interval95=None)
        else:
            report = direct_ffs_report(list(states.values()), dt=self.request['sampling']['dt'],
                                      concentration_m=self.request['concentration_each_m'],
                                      protocol_amendments=self.request.get('protocol_amendments'))
            if (not self.request.get('protocol_amendments')
                    and (len(states) != self.options['replicas'] or status != 'completed')):
                report.update(status='incomplete', value_model=None, inverse_model_time=None,
                              interval95_model=None, value_m_inv_model_time=None)
        report['endpoint'] = 'opening'
        return report

    def run(self):
        processes, logs = {}, {}
        finished, states = set(), {}
        final_status, error = 'completed', None
        try:
            while len(finished) < self.options['replicas']:
                cancelled = (self.output/'CANCEL').exists()
                for index, process in list(processes.items()):
                    if process.poll() is not None:
                        folder = self.output/f'replica_{index+1:03d}'
                        result = read_json(folder/'result.json')
                        if result.get('status') == 'failed' or process.returncode:
                            raise RuntimeError(result.get('error') or f'Native FFS replica {index+1} stopped unexpectedly.')
                        finished.add(index)
                        logs.pop(index).close()
                        del processes[index]
                if cancelled:
                    final_status = 'cancelled'
                    if not processes:
                        break
                else:
                    for index in range(self.options['replicas']):
                        if len(processes) >= self.options['workers']:
                            break
                        if index in finished or index in processes:
                            continue
                        folder = self.output/f'replica_{index+1:03d}'
                        folder.mkdir(exist_ok=True)
                        logs[index] = (folder/'engine.log').open('ab')
                        processes[index] = subprocess.Popen([sys.executable, str(Path(__file__).resolve()),
                            str(self.output/'worker_request.json'), str(index)], cwd=folder,
                            stdout=logs[index], stderr=logs[index], stdin=subprocess.DEVNULL)
                for index in range(self.options['replicas']):
                    folder = self.output/f'replica_{index+1:03d}'
                    state = read_json(folder/'progress.json')
                    if state:
                        states[index] = state
                progress = self.aggregate(states, len(processes), status='running')
                atomic_json(self.output/'progress.json', progress)
                representative = progress['kinetics']['replica']
                live = self.output/f'replica_{representative:03d}'/'live.dat'
                if live.is_file():
                    shutil.copyfile(live, self.output/'live.dat.tmp')
                    (self.output/'live.dat.tmp').replace(self.output/'live.dat')
                # Persist the master identity independently of browser state.
                atomic_json(self.output/'resume_state.json', {'schema_version': 3, 'identity': self.identity,
                    'source_fingerprint': source_fingerprint(self.request), 'elapsed_seconds': progress['elapsed_seconds']})
                if len(finished) < self.options['replicas']:
                    time.sleep(.5)
        except Exception as exc:
            final_status, error = 'failed', str(exc)
            (self.output/'CANCEL').touch()
        finally:
            for process in processes.values():
                try:
                    process.wait(timeout=10)
                except subprocess.TimeoutExpired:
                    process.terminate()
                    process.wait(timeout=10)
            for log in logs.values():
                log.close()
        states = {i: read_json(self.output/f'replica_{i+1:03d}'/'state.json') for i in range(self.options['replicas'])}
        states = {i: row for i, row in states.items() if row}
        result = self.aggregate(states, 0, status=final_status)
        report = self.report(states, final_status)
        result.update(error=error, sampling_complete=final_status == 'completed',
                      rate_estimate=report, trajectory_file=None, repetitions=[{
                          'replica_id': r['replica_id'], 'flux': {k:v for k,v in r['flux'].items() if k != 'bank'},
                          'stages': [{k:v for k,v in stage.items() if k not in ('bank','trials')} for stage in r['stages']],
                          'steps': r['steps'], 'completed': r['completed']} for r in states.values()])
        result['kinetics'].update(phase='report', status=report['status'], rate=report)
        if self.proximity:
            result.update(rate_estimate=None, conditional_probability=report)
            result['kinetics'].pop('rate', None)
            result['kinetics']['conditional_probability'] = report
            for row in result['repetitions']:
                row.pop('flux', None)
        for index in sorted(states):
            movie = self.output/f'replica_{index+1:03d}'/'movie.dat'
            if movie.is_file():
                shutil.copyfile(movie, self.output/'movie.dat')
                shutil.copyfile(movie.with_name('movie_metadata.json'), self.output/'movie_metadata.json')
                result.update(read_json(movie.with_name('movie_metadata.json')))
                result['trajectory_file'] = 'movie.dat'
                break
        atomic_json(self.output/'kinetics_report.json', report)
        atomic_json(self.output/'result.json', result)
        atomic_json(self.output/'progress.json', result)
        return result


if __name__ == '__main__':
    request_path, replica = Path(sys.argv[1]).resolve(), int(sys.argv[2])
    worker = FFSReplica(read_json(request_path), request_path.parent/f'replica_{replica+1:03d}', replica)
    result = worker.run()
    raise SystemExit(1 if result['status'] == 'failed' else 0)
