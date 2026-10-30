"""Joint window probabilities from exact enumeration and the native complex."""
from copy import deepcopy
from dataclasses import asdict, FrozenInstanceError, replace
import json
import math
import os
from pathlib import Path
from types import SimpleNamespace
import unittest
from unittest.mock import patch

from aptaswitch_core.thermodynamics.accessibility import UnpairedWindowAccessibility, unpaired_window_accessibility
from aptaswitch_core.thermodynamics.enumeration import ExhaustiveBackend
from aptaswitch_core.thermodynamics.native_backend import NativeHybridBackend
from aptaswitch_core.thermodynamics.structures import PairingConstraints
from aptaswitch_core.thermodynamics.types import ComplexSpec, Conditions, EmptyEnsembleError, ModelProvenance, Strand


def complex_spec(rna="AAA", dna="T"):
    return ComplexSpec((Strand(rna, "rna", "switch"), Strand(dna, "dna", "trigger")))


def equal_weight_backend():
    # Each admissible interstrand A-T pair has zero energy. Short strands
    # cannot make hairpins; connected structures are enumerated independently.
    return ExhaustiveBackend(
        pair_allowed=lambda spec, i, j: spec.sequence[i] == "A" and spec.sequence[j] == "T",
        energy=lambda spec, conditions, state: 0., thermal_energy=lambda conditions: .6,
        provenance=ModelProvenance("test-enumeration", "1", "equal-weight-AT", "sha256:fixture",
                                   "connected-pseudoknot-free", "1 M"), min_hairpin_size=3)


class JointWindowAccessibilityTests(unittest.TestCase):
    def setUp(self):
        self.backend = equal_weight_backend()
        self.spec, self.conditions = complex_spec(), Conditions()
        self.partition = self.backend.partition(self.spec, self.conditions)

    def test_analytic_joint_ratio_differs_from_product_of_unpaired_marginals(self):
        result = unpaired_window_accessibility(self.backend, self.spec, self.conditions,
                                               start=0, end=2, partition=self.partition)
        self.assertAlmostEqual(result.probability, 1 / 3, places=15)
        self.assertAlmostEqual(result.log_probability, -math.log(3), places=15)
        self.assertAlmostEqual(result.opening_free_energy_kcal_mol, .6 * math.log(3), places=15)
        marginal = self.backend.pairs(self.spec, self.conditions).matrix
        self.assertAlmostEqual(marginal[0][0] * marginal[1][1], 4 / 9, places=15)
        self.assertNotAlmostEqual(result.probability, marginal[0][0] * marginal[1][1])
        self.assertEqual(result.complex_spec, self.spec)
        self.assertEqual(result.conditions, self.conditions)
        self.assertEqual(result.provenance, self.backend.provenance)
        self.assertFalse(result.constrained_ensemble_empty)
        self.assertEqual(result.to_dict(), json.loads(json.dumps(asdict(result), allow_nan=False)))
        json.dumps(result.to_dict(), allow_nan=False)
        restored = UnpairedWindowAccessibility.from_dict(json.loads(json.dumps(result.to_dict(), allow_nan=False)))
        self.assertEqual(restored, result)
        with self.assertRaises(FrozenInstanceError):
            result.probability = 1.

    def test_reuses_denominator_and_preserves_the_whole_complex_and_partial_constraints(self):
        with patch.object(self.backend, "partition", wraps=self.backend.partition) as calculate:
            unpaired_window_accessibility(self.backend, self.spec, self.conditions,
                                          start=0, end=2, partition=self.partition)
        calculate.assert_called_once_with(self.spec, self.conditions, PairingConstraints((0, 1, -1, -1)))
        self.assertIs(calculate.call_args.args[0], self.spec)
        self.assertEqual(calculate.call_args.args[0].ordered_key, (("rna", "AAA"), ("dna", "T")))

    def test_window_can_cross_a_strand_break_without_isolating_any_bases(self):
        spec = complex_spec("AA", "TT")
        result = unpaired_window_accessibility(self.backend, spec, self.conditions, start=1, end=3)
        self.assertAlmostEqual(result.probability, 1 / 5, places=15)
        self.assertEqual(result.complex_spec.cutpoints, (2,))

    def test_impossible_connected_event_is_zero_and_strict_json_safe(self):
        result = unpaired_window_accessibility(self.backend, self.spec, self.conditions,
                                               start=0, end=3, partition=self.partition)
        self.assertEqual(result.probability, 0.)
        self.assertTrue(result.constrained_ensemble_empty)
        self.assertIsNone(result.log_probability)
        self.assertIsNone(result.opening_free_energy_kcal_mol)
        self.assertIsNone(result.constrained_log_partition)
        self.assertIsNone(json.loads(json.dumps(result.to_dict(), allow_nan=False))["opening_free_energy_kcal_mol"])
        self.assertEqual(UnpairedWindowAccessibility.from_dict(result.to_dict()), result)

    def test_empty_unrestricted_ensemble_is_undefined_before_constrained_call(self):
        spec = complex_spec("A", "A")
        empty = self.backend.partition(spec, self.conditions)
        with patch.object(self.backend, "partition") as calculate:
            with self.assertRaises(EmptyEnsembleError):
                unpaired_window_accessibility(self.backend, spec, self.conditions, start=0, end=1, partition=empty)
            calculate.assert_not_called()

    def test_invalid_window_or_unrestricted_context_is_rejected_before_calculation(self):
        bad_partitions = [True,
            replace(self.partition, complex_spec=complex_spec("AAU")),
            replace(self.partition, complex_spec=ComplexSpec((Strand("AAA", "dna", "switch"), self.spec.strands[1]))),
            replace(self.partition, complex_spec=ComplexSpec(tuple(reversed(self.spec.strands)))),
            replace(self.partition, complex_spec=ComplexSpec((replace(self.spec.strands[0], name="alias"), self.spec.strands[1]))),
            replace(self.partition, conditions=Conditions(celsius=38)),
            replace(self.partition, provenance=replace(self.backend.provenance, version="different")),
            replace(self.partition, constraints=PairingConstraints((0, -1, -1, -1)))]
        with patch.object(self.backend, "partition") as calculate:
            for start, end in ((True, 2), (0, False), (-1, 2), (2, 2), (3, 2), (0, 5), (0., 2)):
                with self.subTest(window=(start, end)), self.assertRaises(ValueError):
                    unpaired_window_accessibility(self.backend, self.spec, self.conditions,
                                                  start=start, end=end, partition=self.partition)
            for partition in bad_partitions:
                with self.subTest(partition=partition), self.assertRaises(ValueError):
                    unpaired_window_accessibility(self.backend, self.spec, self.conditions,
                                                  start=0, end=2, partition=partition)
            calculate.assert_not_called()
        for spec, conditions in ((None, self.conditions), (self.spec, {}),
                                 (self.spec, Conditions(magnesium=.001))):
            with self.assertRaises(ValueError):
                unpaired_window_accessibility(self.backend, spec, conditions, start=0, end=2)
        with self.assertRaises(ValueError):
            unpaired_window_accessibility(SimpleNamespace(provenance=None), self.spec,
                                          self.conditions, start=0, end=2)

    def test_constrained_result_must_preserve_context_constraints_and_thermal_energy(self):
        target = self.backend.partition(self.spec, self.conditions, PairingConstraints((0, 1, -1, -1)))
        invalid = [None, replace(target, constraints=None),
            replace(target, complex_spec=complex_spec("AAU")),
            replace(target, conditions=Conditions(celsius=38)),
            replace(target, provenance=replace(target.provenance, version="different")),
            replace(target, thermal_energy_kcal_mol=.7, free_energy_kcal_mol=-.7 * target.log_partition)]
        for target in invalid:
            with self.subTest(target=target), patch.object(self.backend, "partition", return_value=target):
                with self.assertRaises(ValueError):
                    unpaired_window_accessibility(self.backend, self.spec, self.conditions,
                                                  start=0, end=2, partition=self.partition)

    def test_ratio_over_one_is_rejected_but_endpoint_roundoff_is_absorbed(self):
        target = self.backend.partition(self.spec, self.conditions, PairingConstraints((0, 1, -1, -1)))
        for excess in (1., 1e-12):
            log_z = self.partition.log_partition + excess
            altered = replace(target, log_partition=log_z, free_energy_kcal_mol=-.6 * log_z)
            with patch.object(self.backend, "partition", return_value=altered):
                if excess == 1.:
                    with self.assertRaisesRegex(ArithmeticError, "exceeds"):
                        unpaired_window_accessibility(self.backend, self.spec, self.conditions,
                                                      start=0, end=2, partition=self.partition)
                else:
                    result = unpaired_window_accessibility(self.backend, self.spec, self.conditions,
                                                           start=0, end=2, partition=self.partition)
                    self.assertEqual((result.probability, result.log_probability, result.opening_free_energy_kcal_mol),
                                     (1., 0., 0.))
                    self.assertEqual(UnpairedWindowAccessibility.from_dict(result.to_dict()), result)

    def test_underflow_keeps_finite_log_probability_and_opening_energy(self):
        target = self.backend.partition(self.spec, self.conditions, PairingConstraints((0, 1, -1, -1)))
        log_z = self.partition.log_partition - 1000.
        target = replace(target, log_partition=log_z, free_energy_kcal_mol=-.6 * log_z)
        with patch.object(self.backend, "partition", return_value=target):
            result = unpaired_window_accessibility(self.backend, self.spec, self.conditions,
                                                   start=0, end=2, partition=self.partition)
        self.assertEqual(result.probability, 0.)
        self.assertFalse(result.constrained_ensemble_empty)
        self.assertAlmostEqual(result.log_probability, -1000.)
        self.assertAlmostEqual(result.opening_free_energy_kcal_mol, 600.)
        json.dumps(result.to_dict(), allow_nan=False)
        self.assertEqual(UnpairedWindowAccessibility.from_dict(result.to_dict()), result)

    def test_recorded_diagnostics_reject_inconsistent_or_nonfinite_numbers(self):
        result = unpaired_window_accessibility(self.backend, self.spec, self.conditions, start=0, end=2)
        invalid = (("probability", True), ("probability", "0.333"), ("probability", -1.),
                   ("probability", 1.1), ("probability", 0.), ("probability", math.nan),
                   ("log_probability", math.inf), ("log_probability", .1),
                   ("log_probability", None), ("log_probability", -10.),
                   ("thermal_energy_kcal_mol", 0.), ("thermal_energy_kcal_mol", -.6),
                   ("opening_free_energy_kcal_mol", -1.), ("opening_free_energy_kcal_mol", 5.),
                   ("unrestricted_log_partition", math.inf), ("constrained_log_partition", 10.),
                   ("constrained_log_partition", -10.), ("constrained_ensemble_empty", 1),
                   ("constrained_ensemble_empty", True), ("start", True), ("end", 5))
        for name, value in invalid:
            with self.subTest(name=name, value=value), self.assertRaises(ValueError):
                UnpairedWindowAccessibility.from_dict({**result.to_dict(), name: value})
        impossible = unpaired_window_accessibility(self.backend, self.spec, self.conditions, start=0, end=3)
        for name, value in (("constrained_ensemble_empty", False), ("probability", .1),
                            ("log_probability", -math.inf), ("opening_free_energy_kcal_mol", math.inf)):
            with self.subTest(impossible=name), self.assertRaises(ValueError):
                UnpairedWindowAccessibility.from_dict({**impossible.to_dict(), name: value})

    def test_recorded_diagnostics_require_complete_typed_context_and_fields(self):
        result = unpaired_window_accessibility(self.backend, self.spec, self.conditions, start=0, end=2)
        invalid = [None, {**result.to_dict(), "extra": "field"},
                   {name: value for name, value in result.to_dict().items() if name != "probability"}]
        for field, replacement in (("complex_spec", None), ("complex_spec", {"strands": "AAA+T"}),
                ("complex_spec", {"strands": [{"sequence": "T", "material": "rna", "name": "wrong"}]}),
                ("conditions", {"celsius": 37}), ("conditions", {"celsius": 37, "sodium": .12, "magnesium": .001}),
                ("provenance", {})):
            invalid.append({**result.to_dict(), field: replacement})
        for payload in invalid:
            with self.subTest(payload=payload), self.assertRaises(ValueError):
                UnpairedWindowAccessibility.from_dict(payload)
        # The record validates its own context; callers additionally compare it
        # to their candidate. Valid but different labels remain visible.
        altered = deepcopy(result.to_dict())
        altered["complex_spec"]["strands"][0]["name"] = "another label"
        restored = UnpairedWindowAccessibility.from_dict(altered)
        self.assertNotEqual(restored.complex_spec, self.spec)


ROOT = Path(__file__).resolve().parents[1]
LOCAL_LIBRARY = ROOT / ".local_deps/thermodynamics-vienna-off-20261008/libaptaswitch_hybrid.dylib"
LIBRARY = os.environ.get("APTASWITCH_HYBRID_NATIVE_LIBRARY") or (str(LOCAL_LIBRARY) if LOCAL_LIBRARY.is_file() else None)


@unittest.skipUnless(LIBRARY, "explicit native hybrid library required")
class NativeJointWindowAccessibilityTests(unittest.TestCase):
    def setUp(self):
        self.backend = NativeHybridBackend(LIBRARY)

    def test_real_connected_hybrid_matches_analytic_joint_ratio_at_two_conditions(self):
        spec = complex_spec()
        for conditions in (Conditions(celsius=25), Conditions(celsius=55, sodium=1)):
            with self.subTest(conditions=conditions):
                partition = self.backend.partition(spec, conditions)
                result = unpaired_window_accessibility(self.backend, spec, conditions,
                                                       start=0, end=2, partition=partition)
                self.assertAlmostEqual(result.probability, 1 / 3, places=12)
                self.assertAlmostEqual(result.opening_free_energy_kcal_mol,
                                       partition.thermal_energy_kcal_mol * math.log(3), places=12)
                pairs = self.backend.pairs(spec, conditions).matrix
                self.assertNotAlmostEqual(result.probability, pairs[0][0] * pairs[1][1])
                self.assertEqual(result.provenance, partition.provenance)

    def test_real_joint_impossibility_keeps_connected_complex_and_unit_event(self):
        spec, conditions = complex_spec("AA"), Conditions()
        result = unpaired_window_accessibility(self.backend, spec, conditions, start=0, end=2)
        self.assertTrue(result.constrained_ensemble_empty)
        self.assertEqual(result.probability, 0.)
        rna = ComplexSpec((Strand("AAA", "rna", "isolated"),))
        unit = unpaired_window_accessibility(self.backend, rna, conditions, start=0, end=3)
        self.assertEqual((unit.probability, unit.log_probability, unit.opening_free_energy_kcal_mol), (1., 0., 0.))
        self.assertFalse(unit.constrained_ensemble_empty)


if __name__ == "__main__":
    unittest.main()
