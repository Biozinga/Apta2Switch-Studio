"""Experimental connected PF grammar, with local affine face weights.

This Python proof is NOT an application backend or a physical energy model.
Its scalar arithmetic is deliberately unscaled and is intended for small
numerical proofs before a native, numerically scaled implementation. Forward
PF, outside marginals, weighted traceback and minimum-total-energy traceback
use the same productions; none establishes a physical energy model's validity.

The topology criterion is Dirks et al. (2007), sections 2.1 and 2.3:
https://www.its.caltech.edu/~niles/refs/dirks07.pdf
Every pair-bounded face has at most one nick; the root has no physical nick
because it already contains the implicit final-to-first strand break.

The four paired productions are disjoint: hairpin (zero nick/no children),
interior (zero nick/one child), multiloop (zero nick/two or more children),
and exterior (one nick/any children). Interior size is explicitly capped.

An affine forest state records only its own face's physical nick count,
number of child helices (saturated at two), and material mask. Child interiors
never contribute to that mask. A forest walks from a to b by taking the
backbone edge a->a+1, then either consuming a+1 unpaired or jumping over a
child pair (a+1,k). A terminal edge a->b finishes the empty walk. This handles
nicks next to closing pairs and between immediately adjacent helices exactly
once, without summing over alternative placements of the same nick.

For fixed-size masks and affine coefficient families the grammar has
O(n^3+n^2*L^2) time and O(n^2) states (L = max_interior_size). This prototype
uses Python callbacks and dictionaries and makes no performance claim.
"""
from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
from functools import lru_cache
import heapq
from itertools import count
import math
from numbers import Real
import random
from typing import Protocol

from aptaswitch_core.thermodynamics.structures import PairingConstraints, format_structure, is_connected
from aptaswitch_core.thermodynamics.types import ComplexSpec, EmptyEnsembleError


DNA, RNA, MIXED = 1, 2, 3
FaceKey = tuple[int, int, int]  # own-face nicks, child count saturated at 2, material mask
Node = tuple  # R; B,i,j; F,mode,family,a,b,c,h,material_mask


@dataclass(frozen=True)
class GrammarPairProbabilities:
    matrix: tuple[tuple[float, ...], ...]
    partition_function: float
    maximum_endpoint_adjustment: float


class LocalWeights(Protocol):
    """Deterministic nonnegative weights; no hidden global corrections.

    All coefficients must remain unchanged for the lifetime of a grammar.
    ``family`` is the material mask of the completed face, not of the whole
    complex. Branch arguments are always i<j. A closing weight includes the
    outer pair's correctly oriented terminal penalty and, for a multiloop,
    its initiation penalty. It must not charge child branches or free bases.
    """

    def pair_allowed(self, spec: ComplexSpec, i: int, j: int) -> bool: ...
    def hairpin(self, spec: ComplexSpec, i: int, j: int) -> float: ...
    def interior(self, spec: ComplexSpec, i: int, j: int, k: int, l: int) -> float: ...
    def unpaired(self, spec: ComplexSpec, mode: str, family: int, i: int) -> float: ...
    def branch(self, spec: ComplexSpec, mode: str, family: int, i: int, j: int) -> float: ...
    def closing(self, spec: ComplexSpec, mode: str, family: int, i: int, j: int) -> float: ...


def _weight(value: object) -> float:
    if isinstance(value, bool) or not isinstance(value, Real):
        raise ValueError("A local weight must be a finite nonnegative real number")
    result = float(value)
    if not math.isfinite(result) or result < 0:
        raise ValueError("A local weight must be a finite nonnegative real number")
    return result


class ConnectedGrammar:
    """One ordered complex and one fixed set of constraints/model coefficients.

    ``global_factor`` applies an explicitly supplied association/symmetry
    correction once, after the labelled-state recurrence. Material labels and
    nucleotide coordinates retain the input order. No T/U conversion occurs.
    """

    def __init__(self, spec: ComplexSpec, weights: LocalWeights, *, min_hairpin_size: int,
                 max_interior_size: int = 30, constraints: PairingConstraints | None = None,
                 global_factor: float = 1.0, max_nucleotides: int = 80):
        if not isinstance(spec, ComplexSpec):
            raise ValueError("spec must be a ComplexSpec")
        for value, name, minimum in ((min_hairpin_size, "min_hairpin_size", 0),
                                     (max_interior_size, "max_interior_size", 0),
                                     (max_nucleotides, "max_nucleotides", 1)):
            if type(value) is not int or value < minimum:
                raise ValueError(f"{name} must be an integer >= {minimum}")
        self.spec, self.weights, self.n = spec, weights, spec.n_nucleotides
        if self.n > max_nucleotides:
            raise ValueError("Request exceeds the explicit small-prototype nucleotide limit")
        self.min_hairpin_size, self.max_interior = min_hairpin_size, max_interior_size
        self.global_factor = _weight(global_factor)
        self.constraints = constraints if constraints is not None else PairingConstraints.unconstrained(self.n)
        if not isinstance(self.constraints, PairingConstraints) or len(self.constraints.partners) != self.n:
            raise ValueError("Constraints must describe every nucleotide")
        self.materials = tuple(DNA if s.material == "dna" else RNA for s in spec.strands for _ in s.sequence)
        strands = tuple(k for k, s in enumerate(spec.strands) for _ in s.sequence)
        self.cuts = frozenset(spec.cutpoints)
        cut_prefix, blocked_prefix = [0], [0]
        required = self.constraints.partners
        for i in range(self.n):
            cut_prefix.append(cut_prefix[-1] + (i + 1 in self.cuts))
            blocked_prefix.append(blocked_prefix[-1] + (required[i] not in (-1, i)))
        self.cut_prefix, self.blocked_prefix = tuple(cut_prefix), tuple(blocked_prefix)
        self.allowed = set()
        for i in range(self.n):
            for j in range(i + 1, self.n):
                if (required[i] not in (-1, j) or required[j] not in (-1, i)
                        or (strands[i] == strands[j] and j - i - 1 < self.min_hairpin_size)):
                    continue
                accepts = weights.pair_allowed(spec, i, j)
                if type(accepts) is not bool:
                    raise ValueError("pair_allowed must return bool")
                if accepts:
                    self.allowed.add((i, j))
        # Per-instance memoization: no global cache retains previous models.
        self._B = lru_cache(maxsize=None)(self._paired)
        self._F = lru_cache(maxsize=None)(self._forest)

    def _free(self, start: int, stop: int) -> bool:
        return self.blocked_prefix[start] == self.blocked_prefix[stop]

    def _nicks(self, left: int, right: int) -> int:
        """Number of actual backbone cuts in (left, right]."""
        return self.cut_prefix[right] - self.cut_prefix[left]

    def _forest(self, mode: str, family: int, a: int, b: int) -> tuple[tuple[FaceKey, float], ...]:
        terms: dict[FaceKey, list[float]] = defaultdict(list)
        for key, coefficient, children in self._forest_terms(mode, family, a, b):
            terms[key].append(coefficient * math.prod(self._inside(child) for child in children))
        return tuple((key, math.fsum(values)) for key, values in sorted(terms.items()))

    def _forest_terms(self, mode: str, family: int, a: int, b: int):
        if b == a + 1:
            yield (int(b in self.cuts), 0, 0), 1.0, ()
            return
        t, edge_nick = a + 1, int(a + 1 in self.cuts)
        if self._free(t, t + 1):
            factor = _weight(self.weights.unpaired(self.spec, mode, family, t))
            if factor:
                for (c, h, mask), value in self._F(mode, family, t, b):
                    if c + edge_nick <= 1 and value:
                        yield ((c + edge_nick, h, mask | self.materials[t]), factor,
                               (("F", mode, family, t, b, c, h, mask),))
        for k in range(t + 1, b):
            inside = self._B(t, k)
            if not inside:
                continue
            factor = _weight(self.weights.branch(self.spec, mode, family, t, k))
            if factor:
                for (c, h, mask), value in self._F(mode, family, k, b):
                    if c + edge_nick <= 1 and value:
                        yield ((c + edge_nick, min(2, h + 1), mask | self.materials[t] | self.materials[k]),
                               factor, (("B", t, k), ("F", mode, family, k, b, c, h, mask)))

    def _paired(self, i: int, j: int) -> float:
        return math.fsum(coefficient * math.prod(self._inside(child) for child in children)
                         for coefficient, children in self._paired_terms(i, j))

    def _paired_terms(self, i: int, j: int):
        if (i, j) not in self.allowed:
            return
        if self._nicks(i, j) == 0 and self._free(i + 1, j) and j - i - 1 >= self.min_hairpin_size:
            yield _weight(self.weights.hairpin(self.spec, i, j)), ()
        # One child. Only the two exposed backbone runs belong to this face;
        # cuts inside the child's interval are handled by that B state.
        for left_size in range(min(self.max_interior, j - i - 3) + 1):
            k = i + 1 + left_size
            for right_size in range(min(self.max_interior - left_size, j - i - 3 - left_size) + 1):
                l = j - 1 - right_size
                if (self._nicks(i, k) + self._nicks(l, j) != 0
                        or not self._free(i + 1, k) or not self._free(l + 1, j)):
                    continue
                inside = self._B(k, l)
                if inside:
                    yield _weight(self.weights.interior(self.spec, i, j, k, l)), (("B", k, l),)
        boundary = self.materials[i] | self.materials[j]
        for family in (DNA, RNA, MIXED):
            for (c, h, mask), value in self._F("multi", family, i, j):
                if value and c == 0 and h == 2 and (mask | boundary) == family:
                    yield (_weight(self.weights.closing(self.spec, "multi", family, i, j)),
                           (("F", "multi", family, i, j, c, h, mask),))
            for (c, h, mask), value in self._F("exterior", family, i, j):
                if value and c == 1 and (mask | boundary) == family:
                    yield (_weight(self.weights.closing(self.spec, "exterior", family, i, j)),
                           (("F", "exterior", family, i, j, c, h, mask),))

    def _root_terms(self):
        for family in (DNA, RNA, MIXED):
            for (c, h, mask), value in self._F("exterior", family, -1, self.n):
                if value and c == 0 and mask == family:
                    yield self.global_factor, (("F", "exterior", family, -1, self.n, c, h, mask),)

    def _inside(self, node: Node) -> float:
        if node[0] == "R":
            return self.partition_function()
        if node[0] == "B":
            return self._B(node[1], node[2])
        _, mode, family, a, b, c, h, mask = node
        return next((value for key, value in self._F(mode, family, a, b) if key == (c, h, mask)), 0.0)

    def _productions(self, node: Node):
        if node[0] == "R":
            yield from self._root_terms()
        elif node[0] == "B":
            yield from self._paired_terms(node[1], node[2])
        else:
            _, mode, family, a, b, c, h, mask = node
            for key, coefficient, children in self._forest_terms(mode, family, a, b):
                if key == (c, h, mask):
                    yield coefficient, children

    def partition_function(self) -> float:
        """Return Z with no truncation; reject numerical overflow explicitly.

        Relative/absolute underflow remains a limitation of this deliberately
        small unscaled prototype. Production code requires stable arithmetic.
        """
        value = math.fsum(coefficient * math.prod(self._inside(child) for child in children)
                          for coefficient, children in self._root_terms())
        if not math.isfinite(value):
            raise ArithmeticError("Unscaled prototype partition overflowed")
        return value

    def _nonempty_partition(self) -> float:
        value = self.partition_function()
        if value == 0:
            raise EmptyEnsembleError("Connected prototype ensemble is empty or its unscaled PF underflowed")
        return value

    def pairs(self) -> GrammarPairProbabilities:
        """Reverse-mode outside pass; P(i,j) = B(i,j)*adjoint(B(i,j))/Z.

        No cubic production tape is stored: productions are regenerated in
        reverse topological order. B precedes forests of the same span. Tiny
        endpoint overshoots at most 1e-12 are reported; larger errors fail.
        """
        z = self._nonempty_partition()
        root, serial = ("R",), count()
        adjoints = {root: 1.0}
        queue = []
        def enqueue(node):
            if node[0] == "R":
                rank = (-self.n - 2, 0)
            elif node[0] == "B":
                rank = (node[1] - node[2], 0)
            else:
                rank = (node[3] - node[4], 1)
            heapq.heappush(queue, (*rank, next(serial), node))
        enqueue(root)
        while queue:
            *_, node = heapq.heappop(queue)
            adjoint = adjoints[node]
            for coefficient, children in self._productions(node):
                values = tuple(self._inside(child) for child in children)
                if not coefficient or any(value == 0 for value in values):
                    continue
                for k, child in enumerate(children):
                    contribution = adjoint * coefficient * math.prod(values[:k] + values[k + 1:])
                    if contribution:
                        if child not in adjoints:
                            adjoints[child] = 0.0
                            enqueue(child)
                        adjoints[child] += contribution
        matrix = [[0.0] * self.n for _ in range(self.n)]
        adjustment, tolerance = 0.0, 1e-12
        for i, j in sorted(self.allowed):
            probability = self._B(i, j) * adjoints.get(("B", i, j), 0.0) / z
            if not math.isfinite(probability) or probability < -tolerance or probability > 1 + tolerance:
                raise ArithmeticError("Outside pair probability is outside [0,1]")
            corrected = min(1.0, max(0.0, probability))
            adjustment = max(adjustment, abs(corrected - probability))
            matrix[i][j] = matrix[j][i] = corrected
        for i, row in enumerate(matrix):
            unpaired = 1 - math.fsum(row)
            if not math.isfinite(unpaired) or unpaired < -tolerance:
                raise ArithmeticError("Outside pair probabilities sum to more than one")
            row[i] = max(0.0, unpaired)
            adjustment = max(adjustment, abs(row[i] - unpaired))
        return GrammarPairProbabilities(tuple(tuple(row) for row in matrix), z, adjustment)

    def _structure_from_traceback(self, choose) -> str:
        partners = list(range(self.n))
        pending = [("R",)]
        while pending:
            node = pending.pop()
            if node[0] == "B":
                i, j = node[1:]
                partners[i], partners[j] = j, i
            _, children = choose(node)
            pending.extend(reversed(children))
        state = PairingConstraints(tuple(partners))
        lengths = tuple(len(s.sequence) for s in self.spec.strands)
        if not is_connected(state, lengths) or any(
                j >= 0 and state.partners[i] != j for i, j in enumerate(self.constraints.partners)):
            raise ArithmeticError("Grammar traceback violated connectedness or constraints")
        return format_structure(state, lengths)

    def sample(self, sample_count: int, *, seed: int) -> list[str]:
        """Independent seeded draws by complete weighted grammar traceback."""
        if type(sample_count) is not int or sample_count < 0 or type(seed) is not int:
            raise ValueError("sample_count must be nonnegative and seed must be an integer")
        if sample_count == 0:
            return []
        self._nonempty_partition()
        rng = random.Random(seed)
        def choose(node):
            choices, weights = [], []
            for coefficient, children in self._productions(node):
                weight = coefficient * math.prod(self._inside(child) for child in children)
                if weight > 0:
                    choices.append((coefficient, children))
                    weights.append(weight)
            if not choices:
                raise ArithmeticError("No positive traceback production")
            return rng.choices(choices, weights=weights, k=1)[0]
        return [self._structure_from_traceback(choose) for _ in range(sample_count)]

    def minimum(self, *, thermal_energy_kcal_mol: float = 1.0) -> tuple[float, str]:
        """One deterministic minimum of the supplied total labelled-state energy.

        Converts each positive local factor w to -RT ln(w), then uses min-plus
        over the same productions. This includes any global_factor correction,
        and is not a structure-specific physical symmetry/MFE correction.
        """
        thermal = _weight(thermal_energy_kcal_mol)
        if thermal == 0:
            raise ValueError("thermal_energy_kcal_mol must be positive")
        self._nonempty_partition()
        choices = {}
        @lru_cache(maxsize=None)
        def minimum(node):
            candidates = []
            for coefficient, children in self._productions(node):
                if coefficient > 0 and all(self._inside(child) > 0 for child in children):
                    candidates.append((-thermal * math.log(coefficient)
                                       + math.fsum(minimum(child) for child in children), coefficient, children))
            if not candidates:
                raise ArithmeticError("No finite minimum-energy production")
            energy, coefficient, children = min(candidates, key=lambda entry: entry[0])
            choices[node] = coefficient, children
            return energy
        energy = minimum(("R",))
        return energy, self._structure_from_traceback(lambda node: choices[node])

    def cache_sizes(self) -> tuple[int, int]:
        """Numbers of memoized paired and affine-forest states, for diagnostics."""
        return self._B.cache_info().currsize, self._F.cache_info().currsize
