"""Topological connectivity criterion for a future connected-ensemble DP.

This deliberately allows every pair type and arbitrarily short loops: the
invariant concerns planar topology, independent of any energetic alphabet.
Dirks et al., SIAM Review 49 (2007), section 2.1, pp. 70–71, states the
one-nick-per-loop criterion for a connected complex:
https://www.its.caltech.edu/~niles/refs/dirks07.pdf

Our root face already contains the implicit nick from the final 3' end back
to the initial 5' end. Thus it may contain no additional physical cutpoint.
"""
from itertools import combinations, product
import unittest

from aptaswitch_core.thermodynamics.loops import decompose_structure
from aptaswitch_core.thermodynamics.types import ComplexSpec, Strand


def planar_words(length):
    """Independent exhaustive word search and stack parser; no app enumerator."""
    for word in product(".()", repeat=length):
        stack, pairs = [], []
        for position, symbol in enumerate(word):
            if symbol == "(":
                stack.append(position)
            elif symbol == ")":
                if not stack:
                    break
                pairs.append((stack.pop(), position))
        else:
            if not stack:
                yield "".join(word), pairs


def component_count(length, cuts, pairs):
    """Nucleotide graph with real backbone bonds and pairs; independent union-find."""
    parents = list(range(length))
    def root(i):
        while parents[i] != i:
            i = parents[i]
        return i
    for i, j in [(k - 1, k) for k in range(1, length) if k not in cuts] + pairs:
        parents[root(i)] = root(j)
    return len({root(i) for i in range(length)})


def composition(length, cuts):
    boundaries = (0,) + cuts + (length,)
    return ComplexSpec(tuple(Strand("A" * (b - a), "dna" if i % 2 else "rna", f"s{i}")
                             for i, (a, b) in enumerate(zip(boundaries, boundaries[1:]))))


class LoopConnectivityTests(unittest.TestCase):
    def test_every_planar_matching_and_cut_placement_through_eight_nucleotides(self):
        checked = 0
        for length in range(1, 9):
            states = tuple(planar_words(length))
            for n_cuts in range(min(length - 1, 4) + 1):
                for cuts in combinations(range(1, length), n_cuts):
                    current = composition(length, cuts)
                    for word, pairs in states:
                        loops = decompose_structure(current, word, require_connected=False)
                        count = component_count(length, cuts, pairs)
                        face_connected = not loops[0].cutpoints and all(
                            len(face.cutpoints) <= 1 for face in loops[1:])
                        # Stronger Euler/dual-graph identity: each additional
                        # nick on the same face creates another component.
                        predicted_components = 1 + len(loops[0].cutpoints) + sum(
                            max(0, len(face.cutpoints) - 1) for face in loops[1:])
                        context = (length, cuts, word)
                        self.assertEqual(face_connected, count == 1, context)
                        self.assertEqual(predicted_components, count, context)
                        self.assertEqual(len(loops), len(pairs) + 1, context)
                        self.assertEqual(sorted(k for face in loops for k in face.cutpoints),
                                         list(cuts), context)
                        checked += 1
        self.assertEqual(checked, 41226)

    def test_omitting_the_implicit_root_nick_misclassifies_disconnected_strands(self):
        loops = decompose_structure(composition(2, (1,)), "..", require_connected=False)
        self.assertTrue(all(len(face.cutpoints) <= 1 for face in loops))
        self.assertEqual(component_count(2, (1,), []), 2)
        self.assertEqual(loops[0].cutpoints, (1,))

    def test_internal_two_nicks_disconnect_an_enclosed_strand(self):
        loops = decompose_structure(composition(3, (1, 2)), "(.)", require_connected=False)
        self.assertEqual(loops[0].cutpoints, ())
        self.assertEqual(loops[1].cutpoints, (1, 2))
        self.assertEqual(component_count(3, (1, 2), [(0, 2)]), 2)


if __name__ == "__main__":
    unittest.main()
