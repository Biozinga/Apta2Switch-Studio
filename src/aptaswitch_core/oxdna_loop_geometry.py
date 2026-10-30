"""Conservative local loop-surface diagnostic for oxNA preparation only.

Guide pairs close backbone polygons with virtual straight chords. A simple
projection supplies a triangulated spanning surface, lifted to the actual 3D
vertices. This convention detects a foreign backbone passing through a loop
aperture; it is neither a knot invariant nor an equilibrium/trajectory filter.
"""
from itertools import product

import numpy as np

try:
    from .oxdna_geometry import _parse_configuration
except ImportError:  # Standalone bundled worker.
    from oxdna_geometry import _parse_configuration


def _guide_faces(guide, offset=0):
    if not isinstance(guide, str) or set(guide)-set('.()'):
        raise ValueError('Loop geometry requires a pseudoknot-free dot-bracket guide.')
    stack, partners = [], {}
    for i, symbol in enumerate(guide):
        if symbol == '(':
            stack.append(i)
        elif symbol == ')':
            if not stack:
                raise ValueError('Unbalanced loop geometry guide.')
            j = stack.pop()
            partners[i], partners[j] = j, i
    if stack:
        raise ValueError('Unbalanced loop geometry guide.')
    faces = []
    for i, j in sorted(partners.items()):
        if i >= j:
            continue
        cycle, unpaired, k = [i], False, i+1
        while k <= j:
            cycle.append(k)
            if k in partners and k < partners[k] < j:
                cycle.append(partners[k])
                k = partners[k]+1
            else:
                unpaired |= k not in partners
                k += 1
        if unpaired:  # Exclude stacked duplex quadrilaterals, retain bulges.
            faces.append(((i+offset, j+offset), [k+offset for k in cycle]))
    return faces


def _cross(a, b):
    return a[0]*b[1]-a[1]*b[0]


def _simple_polygon(points, epsilon):
    if np.any(np.linalg.norm(np.roll(points, -1, axis=0)-points, axis=1) <= epsilon):
        return False
    for i, a in enumerate(points):
        b = points[(i+1) % len(points)]
        for j in range(i+1, len(points)):
            if j == i+1 or (i == 0 and j == len(points)-1):
                continue
            c, d = points[j], points[(j+1) % len(points)]
            if (np.any(np.maximum(a,b) < np.minimum(c,d)-epsilon)
                    or np.any(np.maximum(c,d) < np.minimum(a,b)-epsilon)):
                continue
            ab = (_cross(b-a,c-a), _cross(b-a,d-a))
            cd = (_cross(d-c,a-c), _cross(d-c,b-c))
            if min(ab) <= epsilon and max(ab) >= -epsilon and min(cd) <= epsilon and max(cd) >= -epsilon:
                return False
    return True


def _projection(points):
    centered = points-points.mean(axis=0)
    _, singular, axes = np.linalg.svd(centered, full_matrices=False)
    scale = max(1., float(np.ptp(centered, axis=0).max()))
    epsilon = 1e-10*scale*scale
    if len(singular) < 3 or singular[1] <= 1e-9*scale:
        raise ValueError('The loop boundary is degenerate or nearly collinear.')
    # A warped loop can overlap itself in its best-fit projection while having
    # a simple projection in another direction. Search a bounded PCA-frame grid.
    candidates = [(0.,0.), *[p for p in product((-2.,-1.,-.5,0.,.5,1.,2.), repeat=2) if p != (0.,0.)]]
    best = None
    for a, b in candidates:
        normal = axes[2]+a*axes[0]+b*axes[1]
        normal /= np.linalg.norm(normal)
        x = axes[0]-normal*np.dot(normal, axes[0])
        x /= np.linalg.norm(x)
        projected = np.column_stack((centered@x, centered@np.cross(normal,x)))
        area = abs(sum(_cross(p,projected[(i+1) % len(points)]) for i,p in enumerate(projected)))/2
        if area <= epsilon or not _simple_polygon(projected, epsilon):
            continue
        if a == b == 0:
            return projected, epsilon
        if best is None or area > best[0]:
            best = area, projected
    if best is None:
        raise ValueError('No simple nondegenerate projection of the loop boundary was found.')
    return best[1], epsilon


def _triangles(points, epsilon, reverse=False):
    sign = np.sign(sum(_cross(p,points[(i+1) % len(points)]) for i,p in enumerate(points)))
    pending, triangles = list(range(len(points))), []
    while len(pending) > 3:
        for k in (range(len(pending)-1,-1,-1) if reverse else range(len(pending))):
            a, b, c = pending[k-1], pending[k], pending[(k+1) % len(pending)]
            if sign*_cross(points[b]-points[a],points[c]-points[b]) <= epsilon:
                continue
            if any(min(sign*_cross(points[y]-points[x],points[p]-points[x])
                       for x,y in ((a,b),(b,c),(c,a))) >= -epsilon
                   for p in pending if p not in (a,b,c)):
                continue
            triangles.append((a,b,c))
            pending.pop(k)
            break
        else:
            raise ValueError('The loop boundary could not be triangulated reliably.')
    if abs(_cross(points[pending[1]]-points[pending[0]],points[pending[2]]-points[pending[0]])) <= epsilon:
        raise ValueError('The loop triangulation ends in a degenerate face.')
    return triangles+[tuple(pending)]


def _intersections(points, triangles, starts, ends, segment_ids):
    """Vectorized segment tests; include internal diagonals, exclude true grazes."""
    direction = ends-starts
    boundary_delta = np.roll(points,-1,axis=0)-points
    boundary_length2 = np.sum(boundary_delta**2,axis=1)
    hits, endpoints = {}, {}
    for indices in triangles:
        a,b,c = points[list(indices)]
        e1,e2 = b-a,c-a
        h = np.cross(direction,e2)
        determinant = h@e1
        inverse = np.divide(1.,determinant,out=np.zeros_like(determinant),where=abs(determinant)>1e-10)
        s = starts-a
        u = inverse*np.einsum('ij,ij->i',s,h)
        q = np.cross(s,e1)
        v = inverse*np.einsum('ij,ij->i',direction,q)
        t = inverse*(q@e2)
        valid = (abs(determinant)>1e-10)&(u>=-1e-9)&(v>=-1e-9)&(u+v<=1+1e-9)&(t>=-1e-8)&(t<=1+1e-8)
        for k in np.flatnonzero(valid):
            point = starts[k]+t[k]*direction[k]
            fractions = np.clip(np.sum((point-points)*boundary_delta,axis=1)/boundary_length2,0,1)
            if np.min(np.linalg.norm(point-points-fractions[:,None]*boundary_delta,axis=1)) <= 1e-7:
                continue  # A true polygon-boundary touch is not an aperture crossing.
            if t[k] <= 1e-8 or t[k] >= 1-1e-8:
                end = int(t[k] > .5)
                vertex = int(segment_ids[k,end])
                outward_sign = int(np.sign(determinant[k]))*(1 if end == 0 else -1)
                endpoints.setdefault(vertex,set()).add((int(k),outward_sign))
                continue
            hits.setdefault(int(k), []).append((float(t[k]), int(np.sign(determinant[k]))))
    crossings = set()
    for k, contacts in hits.items():
        groups = []
        for t,sign in sorted(contacts):
            if not groups or abs(t-groups[-1][0]) > 1e-7:
                groups.append((t,{sign}))
            else:
                groups[-1][1].add(sign)
        if any(len(signs)==1 for _,signs in groups):
            crossings.add(k)  # Opposite signs at one folded edge describe a graze.
    for incident in endpoints.values():
        # A free endpoint or V-shaped touch is not a piercing. Two covalent
        # links on opposite sides can cross exactly at their shared vertex.
        for k,sign in incident:
            if any(other != k and other_sign != sign for other,other_sign in incident):
                crossings.add(k)
    return crossings


def configuration_loop_threading(frame, strand_lengths, chemistries, guides):
    """Reject pierced or unresolved guide faces during preparation only.

    Covalent input must already be unwrapped. Pair closures are guide-based
    virtual chords: callers must independently validate the guide's native
    base pairs before accepting preparation. An empty set of guide loops is
    vacuously clear. Failure to construct a surface is never a passing result.
    """
    lengths, acids, guides = tuple(strand_lengths), tuple(chemistries), tuple(guides)
    if (not lengths or len(lengths)!=len(acids) or len(lengths)!=len(guides)
            or any(type(n) is not int or n<1 for n in lengths)
            or any(acid not in ('rna','dna') for acid in acids)
            or any(not isinstance(g,str) or len(g)!=n for g,n in zip(guides,lengths))):
        raise ValueError('Loop geometry requires matching strands, RNA/DNA chemistries and guides.')
    rows,box = _parse_configuration(frame,sum(lengths))
    backbone = np.empty((len(rows),3))
    faces, bonds, offset = [], [], 0
    for n,acid,guide in zip(lengths,acids,guides):
        part = rows[offset:offset+n]
        a1,a3 = part[:,3:6],part[:,6:9]
        backbone[offset:offset+n] = part[:,:3]+(-.4*a1+.2*a3 if acid=='rna' else -.34*a1+.3408*np.cross(a3,a1))
        faces.extend(_guide_faces(guide,offset))
        bonds.extend((i,i+1) for i in range(offset,offset+n-1))
        offset += n
    bonds = np.asarray(bonds,dtype=int).reshape(-1,2)
    result = {'passed':True,'faces_checked':0,'piercings':[],'unresolved':[],
              'method':'guide_backbone_chords; simple_PCA_projection_search; two_lifted_ear_clippings; periodic_segment_intersections',
              'knot_invariant':False}
    for closing,cycle in faces:
        identity = {'loop_closing_pair':[i+1 for i in closing], 'loop_boundary':[i+1 for i in cycle]}
        points = backbone[cycle]
        try:
            if np.any(np.ptp(points,axis=0)>=box):
                raise ValueError('The loop boundary spans an entire periodic box.')
            projected,epsilon = _projection(points)
            triangulations = [_triangles(projected,epsilon,reverse) for reverse in (False,True)]
        except ValueError as exc:
            result['unresolved'].append({**identity,'reason':str(exc)})
            continue
        result['faces_checked'] += 1
        foreign = bonds[~np.isin(bonds,cycle).any(axis=1)]
        starts,ends = backbone[foreign[:,0]],backbone[foreign[:,1]]
        shift = box*np.floor((points.mean(axis=0)-(starts+ends)/2)/box+.5)
        starts,ends = starts+shift,ends+shift
        detected = [set(),set()]
        for image in product((-1,0,1),repeat=3):
            p,q = starts+box*np.asarray(image),ends+box*np.asarray(image)
            near = np.all(np.maximum(p,q)>=points.min(axis=0)-1e-8,axis=1)&np.all(np.minimum(p,q)<=points.max(axis=0)+1e-8,axis=1)
            indices = np.flatnonzero(near)
            if len(indices):
                for k,triangles in enumerate(triangulations):
                    detected[k].update(indices[list(_intersections(points,triangles,p[near],q[near],foreign[near]))].tolist())
        for index in sorted(detected[0]|detected[1]):
            result['piercings'].append({**identity,'segment':(foreign[index]+1).tolist(),
                                       'triangulations_detecting':sum(index in found for found in detected)})
        if bool(detected[0]) != bool(detected[1]):
            result['unresolved'].append({**identity,'reason':'The two spanning surfaces disagree about a piercing.'})
    result['passed'] = not result['piercings'] and not result['unresolved']
    return result
