import math
import unittest

from aptaswitch_core.oxdna_box import check_periodic_self_clearance, estimate_ffs_box


class BoxSizingTests(unittest.TestCase):
    def estimate(self, **options):
        values = dict(switch_length=104, trigger_length=30, switch_material="rna",
                      trigger_material="dna", temperature_c=37., salt_m=.12)
        values.update(options)
        return estimate_ffs_box(**values)

    def test_reference_candidate_uses_extended_switch_not_fixed_hundred(self):
        box = self.estimate()
        self.assertEqual(box["box_size"], 88.)
        self.assertAlmostEqual(box["contour_lengths"]["switch"], 103 * .761070781051)
        self.assertGreaterEqual(box["box_size"], box["nominal_diameters"]["switch"]
                                + box["interaction_cutoff"] + box["periodic_margin"])
        self.assertLess(box["box_size"], 100.)
        self.assertTrue(box["geometry_monitor_required"])
        self.assertAlmostEqual(box["box_size_nm"], 88 * .8518)
        self.assertAlmostEqual(box["concentration_each_m"],
                               1 / (6.02214076e23 * (88 * .8518e-9) ** 3 * 1000))

    def test_longest_strand_can_be_trigger_and_has_no_box_size_cap(self):
        forward = self.estimate(switch_length=1200, trigger_length=30)
        reversed_roles = self.estimate(switch_length=30, trigger_length=1200,
                                      switch_material="dna", trigger_material="rna")
        self.assertEqual(forward["box_size"], reversed_roles["box_size"])
        self.assertGreater(forward["box_size"], 500.)
        self.assertLess(forward["box_size"], 2 * forward["contour_lengths"]["switch"])

    def test_salt_and_temperature_change_electrostatic_clearance(self):
        normal = self.estimate()
        low_salt = self.estimate(salt_m=.01)
        warm = self.estimate(temperature_c=90.)
        self.assertGreater(low_salt["interaction_cutoff"], normal["interaction_cutoff"])
        self.assertGreater(low_salt["box_size"], normal["box_size"])
        self.assertGreater(warm["interaction_cutoff"], normal["interaction_cutoff"])

    def test_material_changes_contour_and_cutoff_covers_all_pair_types(self):
        mixed = self.estimate()
        dna = self.estimate(switch_material="dna")
        rna = self.estimate(trigger_material="rna")
        self.assertLess(dna["contour_lengths"]["switch"], rna["contour_lengths"]["switch"])
        self.assertGreaterEqual(mixed["interaction_cutoff"], dna["interaction_cutoff"])
        self.assertGreaterEqual(mixed["interaction_cutoff"], rna["interaction_cutoff"])
        hybrid = math.hypot(.34, .3408) + math.hypot(.4, .2)
        hybrid += 4.5 * .3667258 * math.sqrt(310.15 / 300) / math.sqrt(.12)
        self.assertGreaterEqual(mixed["interaction_cutoff"], hybrid)

    def test_small_strands_keep_minimum_image_and_placement_room(self):
        box = self.estimate(switch_length=1, trigger_length=1, salt_m=.01)
        self.assertGreaterEqual(box["box_size"], 2 * box["minimum_periodic_clearance"])
        self.assertGreaterEqual(box["box_size"], 2 * box["initial_separation"] + 2)

    def test_invalid_parameters_are_rejected(self):
        for invalid in ({"switch_length": 0}, {"trigger_length": True},
                        {"switch_material": "experimental-Sugimoto95-Vienna04-loopmean"}, {"salt_m": 0},
                        {"temperature_c": float("nan")}, {"temperature_c": -273.15},
                        {"initial_separation": -1}):
            with self.subTest(invalid=invalid), self.assertRaises(ValueError):
                self.estimate(**invalid)


class PeriodicClearanceTests(unittest.TestCase):
    def test_boundary_crossing_is_unwrapped_before_check(self):
        wrapped = [(x % 20, 4., 3.) for x in (18., 19., 20., 21., 22.)]
        result = check_periodic_self_clearance(wrapped, [5], [20.] * 3, 5.)
        self.assertEqual(result["strand_spans"], [[4., 0., 0.]])
        self.assertEqual(result["minimum_clearance"], 16.)

    def test_extended_chain_is_rejected_before_reaching_image_cutoff(self):
        for axis in range(3):
            positions = [[float(i) if j == axis else 0. for j in range(3)] for i in range(16)]
            with self.subTest(axis=axis), self.assertRaisesRegex(ValueError, "periodic image"):
                check_periodic_self_clearance(positions, [16], [20.] * 3, 5.)

    def test_wrapped_extended_chain_is_still_rejected(self):
        positions = [(float((i + 17) % 20), 0., 0.) for i in range(16)]
        with self.assertRaisesRegex(ValueError, "must not be retained"):
            check_periodic_self_clearance(positions, [16], [20.] * 3, 5.)

    def test_each_chain_is_independent_of_other_chain_translation(self):
        positions = [(0., 0., 0.), (1., 0., 0.), (18., 0., 0.), (19., 0., 0.)]
        result = check_periodic_self_clearance(positions, [2, 2], [20.] * 3, 5.)
        self.assertEqual(result["strand_spans"], [[1., 0., 0.], [1., 0., 0.]])

    def test_geometry_at_safe_boundary_is_allowed(self):
        positions = [(float(i), 0., 0.) for i in range(15)]
        result = check_periodic_self_clearance(positions, [15], [20.] * 3, 5.)
        self.assertEqual(result["minimum_clearance"], result["required_clearance"])

    def test_malformed_geometry_is_rejected(self):
        cases = [([], [], [20.] * 3, 5.), ([(0., 0., 0.)], [2], [20.] * 3, 5.),
                 ([(0., 0., 0.)], [1], [0., 20., 20.], 5.),
                 ([(float("inf"), 0., 0.)], [1], [20.] * 3, 5.),
                 ([(0., 0.)], [1], [20.] * 3, 5.),
                 ([(0., 0., 0.)], [1], [20.] * 3, float("nan"))]
        for args in cases:
            with self.subTest(args=args), self.assertRaises(ValueError):
                check_periodic_self_clearance(*args)


if __name__ == "__main__":
    unittest.main()
