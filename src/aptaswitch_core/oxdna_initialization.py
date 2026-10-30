"""Independent, separated rigid-body placements before native FFS burn-in."""
from __future__ import annotations

import math
from pathlib import Path

import numpy as np

try:
    from .oxdna_box import check_periodic_self_clearance, estimate_ffs_box
    from .oxdna_geometry import _read_configuration, configuration_geometry
    from .oxdna_loop_geometry import configuration_loop_threading
    from .oxdna_kinetics import PROXIMITY_PROTOCOL_VERSION
except ImportError:  # Direct worker invocation in the external oxpy interpreter.
    from oxdna_box import check_periodic_self_clearance, estimate_ffs_box
    from oxdna_geometry import _read_configuration, configuration_geometry
    from oxdna_loop_geometry import configuration_loop_threading
    from oxdna_kinetics import PROXIMITY_PROTOCOL_VERSION


def _rotation(rng):
    """Uniform proper SO(3) rotation from a uniform unit quaternion."""
    u, v, w = rng.random(), 2 * math.pi * rng.random(), 2 * math.pi * rng.random()
    x, y = math.sqrt(1 - u) * math.sin(v), math.sqrt(1 - u) * math.cos(v)
    z, q = math.sqrt(u) * math.sin(w), math.sqrt(u) * math.cos(w)
    return np.array([[1-2*(y*y+z*z), 2*(x*y-z*q), 2*(x*z+y*q)],
                     [2*(x*y+z*q), 1-2*(x*x+z*z), 2*(y*z-x*q)],
                     [2*(x*z-y*q), 2*(y*z+x*q), 1-2*(x*x+y*y)]])


def split_monomer_frames(frame, switch_length):
    """Split coordinates and momenta for physically independent MD managers."""
    lines = frame.splitlines()
    rows = [line for line in lines[3:] if line.strip()]
    if type(switch_length) is not int or not 0 < switch_length < len(rows):
        raise ValueError('Isolated equilibration requires two nonempty monomers.')
    for row in rows:
        values = [float(value) for value in row.split()]
        if len(values) != 15 or not all(math.isfinite(value) for value in values):
            raise ValueError('Isolated equilibration requires finite native phase points.')
    header = '\n'.join(lines[:3]) + '\n'
    return [header + '\n'.join(part) + '\n' for part in (rows[:switch_length], rows[switch_length:])]


def combine_monomer_frames(frames):
    """Assemble independent monomers for a new rigid-body proposal, not a path."""
    if len(frames) != 2:
        raise ValueError('A proximity proposal requires exactly two monomer frames.')
    parts = [frame.splitlines() for frame in frames]
    boxes = [[float(value) for value in part[1].split('=', 1)[1].split()] for part in parts]
    if boxes[0] != boxes[1] or len(boxes[0]) != 3:
        raise ValueError('The isolated monomers must retain the same preparation box.')
    # The artificial assembly starts its own clock; independent burn-ins are
    # recorded separately, never concatenated into an apparent physical path.
    return '\n'.join(['t = 0', parts[0][1], 'E = 0 0 0',
                      *parts[0][3:], *parts[1][3:]]) + '\n'


def random_separated_configuration(request, rng, *, max_attempts=10000):
    """Draw positions/orientations uniformly, conditional on a safe free basin.

    The monomers keep their prepared internal conformations. This draw is
    preparation, not an equilibrium sample or a segment of a kinetic path.
    Native equilibration and continuous source-flux sampling must follow it.
    """
    n, total = request['plan']['switch_length'], request['plan']['nucleotide_count']
    if type(n) is not int or type(total) is not int or not 0 < n < total:
        raise ValueError('Random FFS initialization requires exactly two nonempty strands.')
    if type(max_attempts) is not int or max_attempts < 1:
        raise ValueError('Random FFS initialization requires a positive placement attempt count.')
    source = Path(request['prepared']['configuration_file'])
    rows, box = _read_configuration(source, total, box_size=request['box_size'])
    policy = request.get('box_policy') or estimate_ffs_box(
        n, total-n, 'rna', 'dna', request['sampling']['temperature_c'], request['sampling']['salt_m'])
    cutoff = float(policy['interaction_cutoff'])
    if not math.isfinite(cutoff) or cutoff <= 0:
        raise ValueError('The periodic interaction cutoff must be finite and positive.')
    initial_separation = request.get('initial_separation', 10.)
    separation = float(initial_separation)
    if isinstance(initial_separation, bool) or not math.isfinite(separation) or separation <= 0:
        raise ValueError('The initial interstrand separation must be finite and positive.')
    basin_distance = float(request['plan']['basin_distance'])
    if not math.isfinite(basin_distance) or basin_distance <= 0:
        raise ValueError('The separated basin distance must be finite and positive.')
    separation = max(separation, basin_distance, cutoff + 1.)
    if separation >= float(np.linalg.norm(box / 2)):
        raise ValueError('The periodic box is too small for the requested initial interstrand separation.')

    header = source.read_text().splitlines()[:3]

    def serialize(coordinates):
        return '\n'.join(header) + '\n' + '\n'.join(
            ' '.join(f'{value:.17g}' for value in row) for row in coordinates) + '\n'

    strands = []
    guides = [request[role+'_guide'] for role in ('switch', 'trigger')]
    for part, chemistry, guide in ((rows[:n].copy(), 'rna', guides[0]), (rows[n:].copy(), 'dna', guides[1])):
        # Native checkpoints may wrap each particle independently. Unwrap the
        # covalent chain before a rotation so bonds never become box-length.
        delta = np.diff(part[:, :3], axis=0)
        delta -= box * np.floor(delta / box + .5)
        part[1:, :3] = part[0, :3] + np.cumsum(delta, axis=0)
        geometry = configuration_geometry(serialize(part), [len(part)], [chemistry])
        if not geometry['passed']:
            raise ValueError('The prepared monomer has severe steric contacts or invalid backbone bonds; '
                             'prepare a new conformation before starting FFS.')
        loops = configuration_loop_threading(serialize(part), [len(part)], [chemistry], [guide])
        if not loops['passed']:
            raise ValueError('The prepared monomer threads a guide loop or its loop surface cannot be checked; '
                             'prepare a new conformation before starting FFS.')
        part[:, :3] -= part[:, :3].mean(axis=0)
        strands.append(part)

    for attempt in range(1, max_attempts+1):
        placed, rotations, centres = [], [], []
        for part in strands:
            rotation = _rotation(rng)
            centre = box * np.array([rng.random() for _ in range(3)])
            transformed = np.concatenate([part[:, i:i+3] @ rotation.T for i in range(0, 15, 3)], axis=1)
            transformed[:, :3] += centre
            placed.append(transformed)
            rotations.append(rotation.tolist())
            centres.append(centre.tolist())
        coordinates = np.vstack(placed)
        try:
            guard = check_periodic_self_clearance(coordinates[:, :3].tolist(),
                                                 [n, total-n], box.tolist(), cutoff)
        except ValueError:
            continue
        delta = placed[0][:, None, :3] - placed[1][None, :, :3]
        delta -= box * np.floor(delta / box + .5)
        minimum = float(np.sqrt(np.min(np.sum(delta * delta, axis=2))))
        if minimum < separation:
            continue
        frame = serialize(coordinates)
        geometry = configuration_geometry(frame, [n, total-n], ['rna', 'dna'])
        if not geometry['passed']:
            continue
        loops = configuration_loop_threading(frame, [n, total-n], ['rna', 'dna'], guides)
        if not loops['passed']:
            continue
        return frame, {
            'method': 'uniform_periodic_positions_and_uniform_SO3_rotations_with_separation_rejection',
            'attempts': attempt, 'minimum_interstrand_distance': minimum,
            'required_interstrand_distance': separation, 'box_sides': box.tolist(),
            'strand_centres': centres, 'strand_rotations': rotations,
            'periodic_self_check': guard, 'internal_conformations_changed': False,
            'geometry_check': geometry,
            'loop_threading_check': loops,
            'equilibrium_sample': False,
        }
    raise ValueError(
        f'Unable to place two separated strands safely in the fixed periodic box after {max_attempts} attempts. '
        'Increase the preparation box or inspect the prepared monomers; FFS has not started.')


def random_proximity_configuration(request, rng, *, max_attempts=10000):
    """Prepare an artificial close encounter without deforming either monomer.

    Independent proper rotations precede a random anchor-to-anchor displacement
    sampled uniformly in a spherical shell. Rejection enforces the distance to
    *every* interstrand nucleotide pair, as well as proximity to the toehold.
    This proposal is not an equilibrium or first-crossing ensemble. Native
    hydrogen-bond validation must still confirm zero interstrand associations
    before the frame is accepted for conditional sampling.
    """
    n, total = request['plan']['switch_length'], request['plan']['nucleotide_count']
    if type(n) is not int or type(total) is not int or not 0 < n < total:
        raise ValueError('Proximity initialization requires exactly two nonempty strands.')
    if type(max_attempts) is not int or max_attempts < 1:
        raise ValueError('Proximity initialization requires a positive placement attempt count.')
    options = request.get('kinetics', {})
    if options.get('mode') != 'proximity':
        raise ValueError('An artificial close encounter requires the proximity protocol.')
    distances = []
    for key, default in (('initial_distance_min', 2.), ('initial_distance_max', 2.5)):
        value = options.get(key, default)
        if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or value <= 0:
            raise ValueError('Proximity distances must be finite and positive.')
        distances.append(float(value))
    lower, upper = distances
    if lower > upper:
        raise ValueError('The minimum proximity distance cannot exceed its maximum.')
    pairs = request['plan'].get('toehold_pairs', [])
    if not pairs or any(not isinstance(pair, (list, tuple)) or len(pair) != 2
                        or any(type(index) is not int for index in pair)
                        or not 0 <= pair[0] < n <= pair[1] < total for pair in pairs):
        raise ValueError('Proximity initialization requires explicit interstrand toehold pairs.')
    exposed = sorted({pair[0] for pair in pairs})
    source = Path(request['prepared']['configuration_file'])
    rows, box = _read_configuration(source, total, box_size=request['box_size'])
    policy = request.get('box_policy') or estimate_ffs_box(
        n, total-n, 'rna', 'dna', request['sampling']['temperature_c'], request['sampling']['salt_m'])
    cutoff = float(policy['interaction_cutoff'])
    if not math.isfinite(cutoff) or cutoff <= 0:
        raise ValueError('The periodic interaction cutoff must be finite and positive.')
    if lower >= float(np.linalg.norm(box / 2)):
        raise ValueError('The periodic box is too small for the minimum proximity distance.')

    header = source.read_text().splitlines()[:3]

    def serialize(coordinates):
        return '\n'.join(header) + '\n' + '\n'.join(
            ' '.join(f'{value:.17g}' for value in row) for row in coordinates) + '\n'

    strands = []
    guides = [request[role+'_guide'] for role in ('switch', 'trigger')]
    for part, chemistry, guide in ((rows[:n].copy(), 'rna', guides[0]), (rows[n:].copy(), 'dna', guides[1])):
        # Unwrap before rotating, so a periodic boundary never stretches bonds.
        delta = np.diff(part[:, :3], axis=0)
        delta -= box * np.floor(delta / box + .5)
        part[1:, :3] = part[0, :3] + np.cumsum(delta, axis=0)
        geometry = configuration_geometry(serialize(part), [len(part)], [chemistry])
        if not geometry['passed']:
            raise ValueError('The prepared monomer has severe steric contacts or invalid backbone bonds; prepare a new conformation.')
        loops = configuration_loop_threading(serialize(part), [len(part)], [chemistry], [guide])
        if not loops['passed']:
            raise ValueError('The prepared monomer threads a guide loop or its loop surface cannot be checked.')
        part[:, :3] -= part[:, :3].mean(axis=0)
        strands.append(part)

    for attempt in range(1, max_attempts+1):
        rotations = [_rotation(rng), _rotation(rng)]
        placed = [np.concatenate([part[:, i:i+3] @ rotation.T for i in range(0, 15, 3)], axis=1)
                  for part, rotation in zip(strands, rotations)]
        if options.get('version') == PROXIMITY_PROTOCOL_VERSION:
            # Native oxDNA stores L in the particle's body frame. A global
            # placement rotates r, axes and v, while its L_body components stay
            # unchanged; the rotated orientation already rotates L_lab.
            # Legacy sources redraw momenta after placement and retain their
            # original transformation for exact checkpoint reproducibility.
            for coordinates, part in zip(placed, strands):
                coordinates[:, 12:15] = part[:, 12:15]
        switch_anchor, trigger_anchor = rng.choice(exposed), rng.randrange(total-n)
        z, angle = 2*rng.random()-1, 2*math.pi*rng.random()
        radial = math.sqrt(max(0., 1-z*z))
        direction = np.array([radial*math.cos(angle), radial*math.sin(angle), z])
        radius = (lower**3 + rng.random()*(upper**3-lower**3))**(1/3)
        centre = box * np.array([rng.random() for _ in range(3)])
        placed[0][:, :3] += centre
        displacement = placed[0][switch_anchor, :3] + radius*direction - placed[1][trigger_anchor, :3]
        placed[1][:, :3] += displacement
        coordinates = np.vstack(placed)
        try:
            guard = check_periodic_self_clearance(coordinates[:, :3].tolist(),
                                                 [n, total-n], box.tolist(), cutoff)
        except ValueError:
            continue
        delta = placed[0][:, None, :3] - placed[1][None, :, :3]
        delta -= box * np.floor(delta / box + .5)
        squared_distances = np.sum(delta*delta, axis=2)
        minimum = float(np.sqrt(squared_distances.min()))
        toehold_distance = float(np.sqrt(squared_distances[exposed].min()))
        if minimum < lower or toehold_distance > upper:
            continue
        frame = serialize(coordinates)
        geometry = configuration_geometry(frame, [n, total-n], ['rna', 'dna'])
        if not geometry['passed']:
            continue
        loops = configuration_loop_threading(frame, [n, total-n], ['rna', 'dna'], guides)
        if not loops['passed']:
            continue
        return frame, {
            'method': 'artificial_toehold_encounter_with_uniform_SO3_rotations_and_shell_anchor_rejection',
            'artificial': True, 'equilibrium_sample': False, 'first_crossing_sample': False,
            'native_interstrand_bonds_checked': False,
            'attempts': attempt, 'minimum_interstrand_distance': minimum,
            'minimum_toehold_distance': toehold_distance,
            'required_interstrand_distance': lower, 'maximum_toehold_distance': upper,
            'anchor_pair': [switch_anchor, n+trigger_anchor], 'anchor_distance': radius,
            'box_sides': box.tolist(),
            'strand_centres': [part[:, :3].mean(axis=0).tolist() for part in placed],
            'strand_rotations': [rotation.tolist() for rotation in rotations],
            'periodic_self_check': guard, 'internal_conformations_changed': False,
            'geometry_check': geometry, 'loop_threading_check': loops,
        }
    raise ValueError(
        f'Unable to place an artificial toehold encounter safely in the fixed periodic box after {max_attempts} attempts. '
        'Inspect the prepared monomers and proximity distances; conditional sampling has not started.')
