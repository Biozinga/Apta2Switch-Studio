"""Planar loop decomposition with strand cuts and explicit local chemistry.

This is the structural input to a nearest-neighbor energy model, not an
energy model itself. A loop's chemistry includes only its boundary and free
bases, never the interiors of its child helices. Strand cuts turn the face
they lie on into an exterior loop; they are not counted as backbone bonds.
"""
from __future__ import annotations

from dataclasses import dataclass

from .structures import is_connected, parse_structure
from .types import ComplexSpec


@dataclass(frozen=True)
class Loop:
    kind: str
    closing_pair: tuple[int, int] | None
    child_pairs: tuple[tuple[int, int], ...]
    unpaired_runs: tuple[tuple[int, int], ...]  # half-open intervals, including empty runs
    boundary_segments: tuple[tuple[int, int], ...]  # inclusive backbone paths
    cutpoints: tuple[int, ...]  # break between k-1 and k on this face
    materials: frozenset[str]

    @property
    def unpaired_count(self) -> int:
        return sum(b - a for a, b in self.unpaired_runs)


def decompose_structure(complex_spec: ComplexSpec, structure: str, *,
                        require_connected: bool = True) -> tuple[Loop, ...]:
    """Return the root exterior face followed by pair-bounded faces in 5' order.

    This validates topology only; allowed pair types and loop-size energy
    rules belong to the chosen physical model. Disconnected structures can
    be inspected explicitly, but cannot be evaluated as a connected complex.
    """
    lengths = tuple(len(s.sequence) for s in complex_spec.strands)
    state = parse_structure(structure, lengths)
    if require_connected and not is_connected(state, lengths):
        raise ValueError("The structure does not form one connected complex")
    n = len(state.partners)
    materials = tuple(s.material for s in complex_spec.strands for _ in s.sequence)
    children: dict[tuple[int, int] | None, list[tuple[int, int]]] = {None: []}
    stack = []
    for i, j in enumerate(state.partners):
        if j > i:
            pair = (i, j)
            children[stack[-1] if stack else None].append(pair)
            children[pair] = []
            stack.append(pair)
        elif j < i:
            stack.pop()
    loops = []
    for closing, nested in children.items():
        first, stop = (0, n) if closing is None else (closing[0] + 1, closing[1])
        cursor = first
        runs = []
        for i, j in nested:
            runs.append((cursor, i))
            cursor = j + 1
        runs.append((cursor, stop))
        # Pair ends bound each segment. Root strand ends replace the absent
        # outer closing pair. Adjacent child helices leave an empty free run,
        # but still contribute a backbone bond unless that bond is a cut.
        left = 0 if closing is None else closing[0]
        segments = []
        for i, j in nested:
            segments.append((left, i))
            left = j
        segments.append((left, n - 1 if closing is None else closing[1]))
        cuts = tuple(k for k in complex_spec.cutpoints if any(a < k <= b for a, b in segments))
        if closing is None or cuts:
            kind = "exterior"
        elif not nested:
            kind = "hairpin"
        elif len(nested) == 1:
            left_size, right_size = (b-a for a, b in runs)
            kind = "stack" if left_size == right_size == 0 else (
                "bulge" if min(left_size, right_size) == 0 else "interior")
        else:
            kind = "multiloop"
        boundary = {i for a, b in runs for i in range(a, b)}
        if closing is not None:
            boundary.update(closing)
        for pair in nested:
            boundary.update(pair)
        loops.append(Loop(kind, closing, tuple(nested), tuple(runs), tuple(segments),
                          cuts, frozenset(materials[i] for i in boundary)))
    return tuple(loops)
