"""Origin-label allocation leaves chemical equilibria and cyclic orders intact."""
from collections import defaultdict
import math
import unittest

from aptaswitch_core.thermodynamics.tube_labels import labelled_populations, merge_species
from aptaswitch_core.thermodynamics.tubes import solve_tube
from aptaswitch_core.thermodynamics.types import Strand, enumerate_complexes
from test_thermodynamics_tubes import convention, partition


def solved(species, totals, max_size):
    return solve_tube(species, totals,
        [partition(spec.strands, (len(spec.strands)-1)*10.) for spec in enumerate_complexes(species, max_size)],
        standard_state=convention())


class OriginLabelsTests(unittest.TestCase):
    def test_all_labelled_orders_through_five_and_chemical_totals_are_preserved(self):
        origins = (Strand('ACG', 'dna', 'A'), Strand('ACG', 'dna', 'T'), Strand('ACG', 'rna', 'S'))
        totals = (2e-6, 6e-6, 3e-6)
        species, amounts = merge_species(origins, totals)
        self.assertEqual(len(species), 2)  # Equal alphabet alone does not merge chemistries.
        tube = solved(species, amounts, 5)
        populations = labelled_populations(tube, origins, totals)
        self.assertEqual(len(populations), 95)
        chemical = defaultdict(list)
        for spec, value in populations.items():
            chemical[spec.cyclic_key].append(value)
        for record, expected in zip(tube.partitions, tube.concentrations_molar):
            self.assertAlmostEqual(math.fsum(chemical[record.complex_spec.cyclic_key])/expected, 1., places=12)
        for origin, amount in zip(origins, totals):
            recovered = math.fsum(sum(s.name == origin.name for s in spec.strands)*c for spec, c in populations.items())
            self.assertAlmostEqual(recovered/amount, 1., places=8)

    def test_distinct_mirror_orders_are_not_merged(self):
        species = (Strand('A', 'dna', 'A'), Strand('C', 'dna', 'B'), Strand('G', 'rna', 'C'))
        totals = (1e-6, 2e-6, 3e-6)
        populations = labelled_populations(solved(species, totals, 3), species, totals)
        keys = {tuple(s.name for s in spec.strands) for spec in populations}
        self.assertIn(('A', 'B', 'C'), keys)
        self.assertIn(('A', 'C', 'B'), keys)

    def test_zero_origin_and_completely_absent_pool_remain_zero(self):
        a, b = Strand('AC', 'dna', 'A'), Strand('AC', 'dna', 'B')
        for totals in ((0., 1e-6), (0., 0.)):
            species, amounts = merge_species((a, b), totals)
            result = labelled_populations(solved(species, amounts, 3), (a, b), totals)
            for spec, value in result.items():
                if any(s.name == 'A' for s in spec.strands):
                    self.assertEqual(value, 0.)

    def test_invalid_labels_totals_or_wrong_solved_tube_fail(self):
        a, b = Strand('AC', 'dna', 'A'), Strand('AC', 'dna', 'B')
        for species, totals in (((a, a), (1., 1.)), ((a,), ()), ((a,), (-1.,)), ((a,), (True,)), ((a,), (math.nan,))):
            with self.assertRaises(ValueError):
                merge_species(species, totals)
        tube = solved((a,), (1e-6,), 2)
        with self.assertRaises(ValueError):
            labelled_populations(tube, (a, b), (1e-6, 1e-6))


if __name__ == '__main__':
    unittest.main()
