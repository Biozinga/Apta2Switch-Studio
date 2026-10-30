"""Construct oxNA coordinates from immutable sequences and monomer 2D guides.

Geometric seeds are preparation configurations, never predicted equilibria.
Native restraints apply only during preparation; acceptance uses the unforced
NA energy terms and retains an explicit distinction between guide and OFF basin.
"""
from __future__ import annotations
import hashlib
import json
import math
from pathlib import Path
import shutil
import time
import numpy as np
from .oxdna_geometry import (_number, _place_strands, _read_configuration,
                             _write_report, OxDNAPreparationCancelled)
from .oxdna_loop_geometry import configuration_loop_threading
from .localization import tr
from .scientific_process import CalculationCancelled, external_process_environment, run_stage


def guide_pairs(sequence: str, guide: str, material: str) -> list[tuple[int, int]]:
    if material not in ('dna', 'rna') or not sequence or set(sequence) - set('ACGT' if material == 'dna' else 'ACGU'):
        raise ValueError('The oxDNA sequence must match its DNA/RNA chemistry.')
    if len(guide) != len(sequence) or set(guide) - set('.()'):
        raise ValueError('A monomer dot-bracket guide of the same length is required.')
    stack, pairs = [], []
    allowed = {'AT', 'TA', 'GC', 'CG'} if material == 'dna' else {'AU', 'UA', 'GC', 'CG', 'GU', 'UG'}
    for i, mark in enumerate(guide):
        if mark == '(':
            stack.append(i)
        elif mark == ')':
            if not stack:
                raise ValueError('Unbalanced monomer guide.')
            j = stack.pop()
            if i-j < 4 or sequence[j] + sequence[i] not in allowed:
                raise ValueError('The monomer guide contains an unsupported base pair.')
            pairs.append((j, i))
    if stack:
        raise ValueError('Unbalanced monomer guide.')
    return sorted(pairs)


def _unit(vector):
    length = np.linalg.norm(vector)
    if length < 1e-10:
        raise ValueError('Degenerate geometric orientation.')
    return vector / length


def _rotation(theta):
    c, s = math.cos(theta), math.sin(theta)
    return np.array([[c, -s, 0], [s, c, 0], [0, 0, 1.]])


def _axis_rotation(axis, angle):
    """Proper rotation; a seed variation must never reflect a nucleotide."""
    axis = _unit(np.asarray(axis, dtype=float))
    x, y, z = axis
    cross = np.array([[0., -z, y], [z, 0., -x], [-y, x, 0.]])
    return np.eye(3)*math.cos(angle) + (1-math.cos(angle))*np.outer(axis, axis) + math.sin(angle)*cross


def _perpendicular(vector):
    vector = _unit(vector)
    axis = np.eye(3)[int(np.argmin(np.abs(vector)))]
    return _unit(np.cross(vector, axis))


def _helix_pair(depth, material):
    # Ideal A-form RNA and B-form DNA seed; the native NA potential relaxes
    # these approximations. Sequence ordering is explicitly 5-prime to 3-prime.
    if material == 'rna':
        inclination = math.radians(15.5)
        c, s = math.cos(inclination), math.sin(inclination)
        radial = math.sqrt(1.175**2-c*c)
        a1, a3 = np.array([0., -c, s]), np.array([0., s, c])
        left = np.array([radial, c, -s]) + .4*a1
        right = np.array([radial, -c, s]) - .4*a1
        rise, twist = .3287, math.radians(32.7)
    else:
        a1, a3 = np.array([1., 0, 0]), np.array([0., 0, 1.])
        left, right = -.6*a1, .6*a1
        rise, twist = .3897628551303122, math.radians(35.9)
    rotation = _rotation(-twist*depth)
    translation = np.array([0., 0, -rise*depth])
    return (np.r_[rotation@left+translation, rotation@a1, rotation@a3, np.zeros(6)],
            np.r_[rotation@right+translation, -rotation@a1, -rotation@a3, np.zeros(6)])


def _backbone(rows, material):
    a1, a3 = rows[..., 3:6], rows[..., 6:9]
    if material == 'rna':
        return rows[..., :3] - .4*a1 + .2*a3
    return rows[..., :3] - .34*a1 + .3408*np.cross(a3, a1)


def seed_strand(sequence: str, guide: str, material: str, *, seed: int | None = None) -> np.ndarray:
    """Seed helices and contour-aware loops, with reproducible geometric variants.

    These are construction guesses, not thermodynamic samples. The native
    relaxation and strict force-field checks remain authoritative.
    """
    if seed is not None and (type(seed) is not int or not 0 <= seed < 2**31):
        raise ValueError('A geometric seed must be an integer in [0, 2**31).')
    rng = np.random.default_rng(seed) if seed is not None else None
    pairs = guide_pairs(sequence, guide, material)
    partners = {i: j for a, b in pairs for i, j in ((a, b), (b, a))}
    n = len(sequence)
    rows = np.full((n, 15), np.nan)
    children = {pair: [] for pair in pairs}
    roots = []
    for pair in pairs:
        enclosed = [p for p in pairs if p[0] < pair[0] < pair[1] < p[1]]
        if enclosed:
            parent = min(enclosed, key=lambda p: p[1]-p[0])
            children[parent].append(pair)
        else:
            roots.append(pair)
    def place(pair, origin, rotation, depth=0.):
        first, second = _helix_pair(depth, material)
        for i, row in zip(pair, (first, second)):
            row[:3] = rotation@row[:3]+origin
            row[3:6] = rotation@row[3:6]
            row[6:9] = rotation@row[6:9]
            rows[i] = row
        subs = children[pair]
        if len(subs) == 1:
            sub = subs[0]
            # Internal bulges create a short gap along a coaxial stem.
            # The short arm limits separation: a long bulge on the other arm
            # must not stretch a directly bonded pair of stem nucleotides.
            gap = 1 + .5*min(sub[0]-pair[0]-1, pair[1]-sub[1]-1)
            place(sub, origin, rotation, depth+gap)
        else:
            if not subs:
                return
            # A junction is a polygon whose pair edges are rigid helix ends
            # and whose intervening edges fit the available loop contour.
            # A fixed fan radius overstretched short links between branches.
            parent_left, parent_right = _backbone(rows[list(pair)], material)
            first, second = _helix_pair(0, material)
            first_backbone, second_backbone = _backbone(np.vstack([first, second]), material)
            pair_length = np.linalg.norm(second_backbone-first_backbone)
            lengths = []
            previous = pair[0]
            for sub in subs:
                contour = (sub[0]-previous)*.76
                fraction = rng.uniform(.72, .9) if rng is not None else .82
                lengths.extend([min(3., max(.76, fraction*contour)), pair_length])
                previous = sub[1]
            contour = (pair[1]-previous)*.76
            fraction = rng.uniform(.72, .9) if rng is not None else .82
            lengths.extend([min(3., max(.76, fraction*contour)), pair_length])
            low, high = max(lengths)/2, sum(lengths)
            for _ in range(48):
                radius = (low+high)/2
                total = sum(2*math.asin(min(1., length/(2*radius))) for length in lengths)
                if total > 2*math.pi:
                    low = radius
                else:
                    high = radius
            radius = (low+high)/2
            x = _unit(parent_right-parent_left)
            z = -rows[pair[0], 6:9]
            z = _unit(z-x*np.dot(z, x))
            height = math.sqrt(max(0., radius*radius-pair_length*pair_length/4))
            center = (parent_left+parent_right)/2 + height*z
            angle = math.atan2(-height, -pair_length/2)
            points = [parent_left]
            for length in lengths[:-1]:
                angle -= 2*math.asin(min(1., length/(2*radius)))
                points.append(center+radius*(math.cos(angle)*x+math.sin(angle)*z))
            source_x = _unit(second_backbone-first_backbone)
            source_z = -first[6:9]
            source_z = _unit(source_z-source_x*np.dot(source_z, source_x))
            source_frame = np.column_stack([source_x, np.cross(source_z, source_x), source_z])
            for k, sub in enumerate(subs):
                start, end = points[2*k+1:2*k+3]
                target_x = _unit(end-start)
                target_z = (start+end)/2-center
                target_z = _unit(target_z-target_x*np.dot(target_z, target_x))
                target_frame = np.column_stack([target_x, np.cross(target_z, target_x), target_z])
                branch_rotation = target_frame@source_frame.T
                if rng is not None:
                    branch_rotation = _axis_rotation(target_x, rng.uniform(-.2, .2))@branch_rotation
                place(sub, start-branch_rotation@first_backbone, branch_rotation)
    for k, pair in enumerate(roots):
        if not k:
            place(pair, np.zeros(3), np.eye(3))
            continue
        previous = roots[k-1][1]
        first, second = _helix_pair(0, material)
        # Progress along one exterior axis in sequence order. Aligning each
        # new helix with the preceding return strand reversed that progress:
        # every third identical hairpin occupied the first one's coordinates.
        exterior = _unit(_backbone(second, material)-_backbone(first, material))
        rotation = np.eye(3)
        if rng is not None:
            rotation = _axis_rotation(exterior, rng.uniform(-.6, .6))
        # Use the available covalent contour, including a single link for
        # adjacent hairpins. A fixed cap folded long exterior loops back over
        # existing stems and tails. No backbone link is stretched to separate
        # the helices, and each helix keeps its exact guide-pair geometry.
        gap = (pair[0]-previous)*.76
        first_backbone = _backbone(rows[previous], material) + exterior*gap
        origin = first_backbone - rotation@_backbone(first, material)
        place(pair, origin, rotation)
    assigned = [i for i in range(n) if np.isfinite(rows[i, 0])]
    if not assigned:
        # Open strands begin as extended, non-overlapping backbones, not duplexes.
        rows[:] = 0
        rows[:, 3] = 1
        rows[:, 8] = 1
        rows[:, 2] = -np.arange(n)*.76
        if rng is not None and n > 1:
            axis = _perpendicular(rows[0, 6:9])
            bend = rng.uniform(-.8, .8)
            backbone = _backbone(rows[0], material)
            for k in range(1, n):
                rotation = _axis_rotation(axis, bend*k/(n-1))
                rows[k, 3:6] = rotation@rows[0, 3:6]
                rows[k, 6:9] = rotation@rows[0, 6:9]
                backbone -= rows[k, 6:9]*.76
                rows[k, :3] = backbone - _backbone(np.r_[np.zeros(3), rows[k, 3:]], material)
        return rows
    boundaries = [-1, *assigned, n]
    for left, right in zip(boundaries, boundaries[1:]):
        count = right-left-1
        if not count:
            continue
        if left < 0 or right == n:
            anchor = right if left < 0 else left
            direction = 1 if left < 0 else -1
            tangent = rows[anchor, 6:9]*direction
            # Tail variations bend a backbone of unchanged contour length.
            axis = _perpendicular(rows[anchor, 6:9])
            bend = rng.uniform(-.8, .8) if rng is not None else 0.
            backbone = _backbone(rows[anchor], material)
            for k in range(1, count+1):
                index = anchor-k if left < 0 else anchor+k
                row = rows[anchor].copy()
                rotation = _axis_rotation(axis, bend*k/count)
                row[3:6] = rotation@row[3:6]
                row[6:9] = rotation@row[6:9]
                backbone += rotation@tangent*.76
                row[:3] = backbone - _backbone(np.r_[np.zeros(3), row[3:]], material)
                rows[index] = row
            continue
        start, end = _backbone(rows[[left, right]], material)
        chord = end-start
        midpoint = (start+end)/2
        # Hairpins continue beyond the stem tip. Internal loops instead
        # protrude radially away from the partners of their two anchors;
        # following -a3 there can wrap a bulge around a neighbouring turn
        # of the helix, even when no nucleotide sites collide.
        outward = -rows[left, 6:9]
        if partners[left] != right:
            opposite = _backbone(rows[[partners[left], partners[right]]], material).mean(axis=0)
            outward = midpoint-opposite
        outward -= _unit(chord)*np.dot(outward, _unit(chord))
        if np.linalg.norm(outward) < .01:
            outward = _perpendicular(chord)
        outward = _unit(outward)
        if rng is not None:
            outward = _axis_rotation(chord, rng.uniform(-.65, .65))@outward
        # A discrete circular arc keeps every backbone link at its nominal
        # length. Sampling a smooth sinusoid by arclength compressed links
        # spanning a tight turn, especially in long loops with nearby ends.
        segments = count+1
        chord_length = np.linalg.norm(chord)
        if chord_length >= segments*.76-1e-10:
            points = np.array([start+chord*t for t in np.linspace(0, 1, count+2)[1:-1]])
            loop_center = midpoint
        else:
            low, high = 0., 2*math.pi/segments
            for _ in range(48):
                angle = (low+high)/2
                span = .76*math.sin(segments*angle/2)/math.sin(angle/2)
                if span > chord_length:
                    low = angle
                else:
                    high = angle
            angle = (low+high)/2
            radius = .76/(2*math.sin(angle/2))
            half_angle = segments*angle/2
            loop_center = midpoint-radius*math.cos(half_angle)*outward
            points = np.array([loop_center+radius*(math.sin(-half_angle+k*angle)*_unit(chord)
                        + math.cos(-half_angle+k*angle)*outward) for k in range(1, segments)])
        all_points = np.vstack([start, points, end])
        for k, point in enumerate(points, 1):
            tangent = _unit(all_points[k-1]-all_points[k+1])
            a1 = point-loop_center
            a1 -= tangent*np.dot(a1, tangent)
            if np.linalg.norm(a1) < .01:
                a1 = _perpendicular(tangent)
            a1 = _unit(a1)
            row = np.r_[point, a1, tangent, np.zeros(6)]
            row[:3] += point-_backbone(row, material)
            rows[left+k] = row
    return rows


def _configuration(path, rows, box):
    path.write_text(f't = 0\nb = {box:.12g} {box:.12g} {box:.12g}\nE = 0 0 0\n' +
                    '\n'.join(' '.join(f'{v:.15g}' for v in row) for row in rows)+'\n', encoding='utf-8')


PREPARATION_ATTEMPTS = 10
PREPARATION_ATTEMPT_SECONDS = 180.


def prepare_guided_system(request: dict, directory: Path, *, progress_callback=None,
                          cancel_check=None, deadline=None) -> dict:
    """Try independent geometric/thermal seeds; accept only a strict NA check.

    The attempt budget covers preparation only, never the production FFS.
    Every attempt keeps its inputs and diagnostics. Sequences, guides and
    physical conditions remain immutable across retries.
    """
    directory = Path(directory).resolve()
    directory.mkdir(parents=True, exist_ok=True)
    if any(directory.iterdir()):
        raise ValueError('oxDNA preparation requires an empty output directory')
    started = time.monotonic()
    report = {'status': 'running', 'method': 'sequence_and_recorded_2d', 'model': 'NA',
              'external_forces': False, 'equilibrium_validated': False,
              'secondary_structure_constrained': True, 'production_external_forces': False,
              'report_path': str(directory/'preparation.json'), 'sequence_parameter_files': {},
              'guide_model': request.get('guide_model'),
              'guide_source': request.get('guide_source'),
              'guide_thermodynamics': request.get('guide_thermodynamics', {}),
              'attempts': [], 'attempt_limit': PREPARATION_ATTEMPTS,
              'attempt_budget_seconds': PREPARATION_ATTEMPT_SECONDS}

    def save():
        report['elapsed_seconds'] = time.monotonic()-started
        _write_report(directory/'preparation.json', report)

    def check():
        if cancel_check and cancel_check():
            raise OxDNAPreparationCancelled('oxDNA preparation cancelled')
        if deadline is not None and time.monotonic() >= deadline:
            raise TimeoutError('oxDNA preparation exceeded its time limit')

    try:
        check()
        for name in ('temperature_c', 'salt_m', 'box_size', 'initial_separation'):
            value = _number(request.get(name, {'box_size':100, 'initial_separation':6}.get(name)), name,
                            positive=name != 'temperature_c')
            if name == 'temperature_c' and value <= -273.15:
                raise ValueError('Temperature must exceed absolute zero.')
        if type(request.get('seed')) is not int or not 0 <= request['seed'] < 2**31:
            raise ValueError('oxDNA seed must be an integer in [0, 2**31).')
        if request.get('use_average_seq', False) is not False:
            raise ValueError('Guided NA preparation requires sequence-dependent parameters.')
        guides, offset = {}, 0
        for role in ('switch', 'trigger'):
            sequence, material = request[role+'_sequence'], request[role+'_material']
            guide = request[role+'_guide']
            pairs = guide_pairs(sequence, guide, material)
            guides[role] = {'sequence': sequence, 'material': material, 'dot_bracket': guide,
                'pairs': [[a+offset,b+offset] for a,b in pairs],
                'sha256': hashlib.sha256((material+'\n'+sequence+'\n'+guide).encode()).hexdigest()}
            offset += len(sequence)
        report['guides'] = guides
        # No amount of restarting can reconcile contradictory input targets.
        stem = {tuple(p) for p in request.get('plan', {}).get('stem_pairs', [])}
        expected = {tuple(p) for g in guides.values() for p in g['pairs']}
        if len(stem & expected) < request.get('plan', {}).get('basin_min_stem', 0):
            raise ValueError(tr('Le guide 2D ne contient pas assez de paires de la tige OFF requise. Les réessais ne peuvent pas résoudre cette incompatibilité.'))
        topology = directory/'system.top'
        topology.write_text(f'{offset} 2 5->3\n'+ '\n'.join(
            f'{g["sequence"]} type={g["material"].upper()}' for g in guides.values())+'\n')
        report['topology_file'] = str(topology)
        local_request = json.loads(json.dumps(request))
        for chemistry in ('DNA','RNA','NA'):
            key = 'seq_dep_file_'+chemistry
            source = Path(request['runtime'][key]).resolve()
            target = directory/(key+'.txt')
            shutil.copyfile(source, target)
            report['sequence_parameter_files'][chemistry] = {'source':str(source), 'copy':str(target),
                'sha256':hashlib.sha256(target.read_bytes()).hexdigest()}
            local_request['runtime'][key] = str(target)
        local_request.update(guides=guides, topology_file=str(topology))
        policy = request.get('box_policy')
        placement_options = ({'periodic_margin': policy['minimum_periodic_clearance'],
                              'auto_size': True, 'reorient': True} if policy else {})
        worker = Path(__file__).with_name('oxdna_guided_worker.py')
        env = external_process_environment()
        env['PYTHONDONTWRITEBYTECODE'] = '1'
        phases = {'minimize_guided': 'Relaxation initiale', 'relax_guided': 'Ajustement des boucles',
                  'guide_physical': 'Stabilisation des appariements', 'release': 'Retrait des contraintes',
                  'probe': 'Contrôle sans contraintes'}
        best, best_rank = None, None
        save()
        for index in range(1, PREPARATION_ATTEMPTS+1):
            check()
            attempt_started = time.monotonic()
            attempt_end = attempt_started + PREPARATION_ATTEMPT_SECONDS
            if deadline is not None:
                attempt_end = min(attempt_end, deadline)
            folder = directory/f'attempt_{index:02d}'
            folder.mkdir()
            seed = (request['seed'] + 104729*(index-1)) % (2**31-1)
            attempt = {'index': index, 'seed': seed,
                       'geometry_variant': 'baseline' if index == 1 else 'varied',
                       'restraint_stiffness': (.9, 3., 6.)[(index-1) % 3],
                       'directory': str(folder), 'status': 'running'}
            report['attempts'].append(attempt)
            seeds = [seed_strand(g['sequence'], g['dot_bracket'], g['material'],
                        seed=None if index == 1 else (seed+role*65537) % (2**31-1))
                     for role, g in enumerate(guides.values())]
            rows, box = _place_strands(*seeds, float(request.get('initial_separation', 6)),
                                      float(request.get('box_size', 100)), **placement_options)
            initial = folder/'system.dat'
            _configuration(initial, rows, box)
            _read_configuration(initial, offset, box_size=box)
            attempt_request = {**local_request, 'seed': seed, 'configuration_file': str(initial),
                'box_size': box, 'preparation_attempt': {key: attempt[key] for key in
                    ('index', 'seed', 'geometry_variant', 'restraint_stiffness')},
                'budget_seconds': max(.01, attempt_end-time.monotonic())}
            _write_report(folder/'worker_request.json', attempt_request)
            attempt['loop_threading'] = configuration_loop_threading(initial.read_text(),
                [len(g['sequence']) for g in guides.values()],
                [g['material'] for g in guides.values()],
                [g['dot_bracket'] for g in guides.values()])
            if not attempt['loop_threading']['passed']:
                attempt.update(status='rejected', failure_kind='loop_threading',
                               elapsed_seconds=time.monotonic()-attempt_started)
                save()
                continue
            save()

            def notify(name='', elapsed=0, detail=''):
                if progress_callback:
                    try:
                        phase = json.loads((folder/'guided_progress.json').read_text()).get('phase')
                    except (OSError, ValueError):
                        phase = None
                    label = tr(phases.get(phase, 'Construction de la conformation'))
                    progress_callback(tr('Préparation 3D guidée · essai {index}/{total} · {phase} · {seconds:.0f} s',
                        index=index, total=PREPARATION_ATTEMPTS, phase=label, seconds=elapsed))

            notify()
            result = None
            try:
                run_stage([request['runtime']['python'], str(worker), str(folder/'worker_request.json')],
                          cwd=folder, env=env, stage='oxdna_guided', log_path=folder/'guided.log',
                          deadline=attempt_end, cancel_check=cancel_check, progress_callback=notify)
                check()
                result = json.loads((folder/'guided_result.json').read_text())
                if result.get('status') != 'completed':
                    raise RuntimeError(result.get('error', 'Native preparation failed.'))
                final_box = float(result.get('box_size', box))
                if final_box < box or not math.isfinite(final_box):
                    raise ValueError('The prepared box changed below its initial contour and placement requirement.')
                configuration = Path(result['configuration_file']).resolve()
                if not configuration.is_relative_to(folder):
                    raise ValueError('The prepared configuration must belong to its attempt.')
                _read_configuration(configuration, offset, box_size=final_box)
                validation = result['validation']
                # Never infer acceptance from file existence or a relaxed stage.
                eligible = (result.get('ffs_eligible') is True
                    and result.get('strict_force_field_initialization_passed') is True
                    and validation.get('guide_exact_match') is True
                    and validation.get('geometry', {}).get('passed') is True
                    and validation.get('loop_threading', {}).get('passed') is True
                    and validation.get('off_basin_eligible') is True)
                result['ffs_eligible'] = eligible
                attempt.update(status='accepted' if eligible else 'rejected', validation=validation,
                               selected_checkpoint=result.get('selected_checkpoint'))
                missing = sum(len(g['missing_pairs']) for g in validation['guides'].values())
                rank = (not eligible, len(validation.get('interstrand_bonds', [])), missing,
                        len(validation.get('unexpected_pairs', [])), not validation.get('off_basin_eligible', False))
                if result.get('strict_force_field_initialization_passed') and (best_rank is None or rank < best_rank):
                    best, best_rank = {**result, 'selected_attempt': index,
                                      'initial_configuration_file': str(initial)}, rank
                if eligible:
                    report.update(best, retry_exhausted=False)
                    return report
            except CalculationCancelled:
                attempt['status'] = 'cancelled'
                raise
            except OxDNAPreparationCancelled:
                attempt['status'] = 'cancelled'
                raise
            except (FileNotFoundError, PermissionError):
                attempt.update(status='failed', failure_kind='environment')
                raise
            except (RuntimeError, TimeoutError, ValueError, KeyError) as exc:
                failure = result or {}
                if not failure:
                    try:
                        failure = json.loads((folder/'guided_result.json').read_text())
                    except (OSError, ValueError):
                        pass
                kind = 'timeout' if isinstance(exc, TimeoutError) else failure.get('failure_kind', 'validation')
                attempt.update(status='failed', error=failure.get('error', str(exc)), failure_kind=kind)
                try:
                    check()
                except OxDNAPreparationCancelled:
                    attempt['status'] = 'cancelled'
                    raise
                if kind == 'environment':
                    raise RuntimeError(attempt['error']) from exc
            finally:
                attempt['elapsed_seconds'] = time.monotonic()-attempt_started
                save()
        check()
        report['retry_exhausted'] = True
        if best is not None:
            report.update(best)
            return report
        raise RuntimeError(tr('Aucune conformation valide après {count} essais de préparation. Consultez les diagnostics de préparation.', count=len(report['attempts'])))
    except CalculationCancelled as exc:
        report.update(status='cancelled',error=str(exc))
        raise OxDNAPreparationCancelled(str(exc)) from exc
    except Exception as exc:
        report.update(status='cancelled' if isinstance(exc,OxDNAPreparationCancelled) else 'failed',error=str(exc))
        raise
    finally:
        if report['status'] != 'completed':
            report.pop('configuration_file', None)
            report['ffs_eligible'] = False
        save()
