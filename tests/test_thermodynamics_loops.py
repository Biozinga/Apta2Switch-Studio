import unittest

from aptaswitch_core.thermodynamics.loops import decompose_structure
from aptaswitch_core.thermodynamics.structures import PairingConstraints, format_structure
from aptaswitch_core.thermodynamics.types import ComplexSpec, Strand


class LoopDecompositionTests(unittest.TestCase):
    def test_hybrid_duplex_has_exterior_faces_and_stacks(self):
        spec = ComplexSpec((Strand("GGG", "rna", "r"), Strand("CCC", "dna", "d")))
        loops = decompose_structure(spec, "(((+)))")
        self.assertEqual([x.kind for x in loops], ["exterior", "stack", "stack", "exterior"])
        self.assertEqual([x.cutpoints for x in loops], [(), (), (), (3,)])
        self.assertTrue(all(x.materials == {"rna", "dna"} for x in loops))

    def test_monomer_hairpin_and_free_bases(self):
        spec = ComplexSpec((Strand("AGGAAACCU", "rna", "r"),))
        loops = decompose_structure(spec, ".((...)).")
        self.assertEqual([x.kind for x in loops], ["exterior", "stack", "hairpin"])
        self.assertEqual([x.unpaired_count for x in loops], [2, 0, 3])

    def test_bulge_and_interior_are_distinct(self):
        for structure, kind in (("(.(...))", "bulge"), ("(.(...).)", "interior")):
            spec = ComplexSpec((Strand("A" * len(structure), "rna", "r"),))
            self.assertEqual(decompose_structure(spec, structure)[1].kind, kind)

    def test_child_interior_chemistry_does_not_leak_into_parent_loop(self):
        # DNA 5..8 is inside child helix 2..16, not on the 0..24 multiloop.
        letters = list("A" * 25)
        pairs = ((0, 24), (2, 16), (4, 8), (18, 22))
        partners = list(range(25))
        for i, j in pairs:
            letters[i], letters[j] = "G", "C"
            partners[i], partners[j] = j, i
        spec = ComplexSpec((Strand("".join(letters[:5]), "rna", "r1"),
                            Strand("".join(letters[5:9]), "dna", "d"),
                            Strand("".join(letters[9:]), "rna", "r2")))
        structure = format_structure(PairingConstraints(tuple(partners)), (5, 4, 16))
        loops = decompose_structure(spec, structure)
        parent = next(x for x in loops if x.closing_pair == (0, 24))
        self.assertEqual(parent.kind, "multiloop")
        self.assertEqual(parent.materials, {"rna"})
        self.assertEqual(parent.cutpoints, ())
        # Each cut belongs to one face, despite being inside several intervals.
        self.assertEqual(sorted(k for x in loops for k in x.cutpoints), [5, 9])
        free = [i for x in loops for a, b in x.unpaired_runs for i in range(a, b)]
        self.assertEqual(sorted(free), [i for i, j in enumerate(partners) if i == j])

    def test_disconnected_structures_require_explicit_diagnostic_mode(self):
        spec = ComplexSpec((Strand("A", "rna", "r"), Strand("A", "dna", "d")))
        with self.assertRaisesRegex(ValueError, "connected"):
            decompose_structure(spec, ".+.")
        loops = decompose_structure(spec, ".+.", require_connected=False)
        self.assertEqual(loops[0].cutpoints, (1,))


if __name__ == "__main__":
    unittest.main()
