"""Reference OSTIR calculations and explicit ideal-occupancy projection checks."""
from importlib.metadata import PackageNotFoundError, version
from copy import deepcopy
import json
import math
from types import SimpleNamespace
import unittest
from unittest.mock import patch

from aptaswitch_core import rbs_calculator as rbs


PAPER_SEQUENCE = "ACUUCUAAUUUAUUCUAUUUAUUCGCGGAUAUGCAUAGGAGUGCUUCGAUGUCAU"
CLOSED_SEQUENCE = "GGGCUCCUUAAAAAAGGAGAAAAAAAUG" + "GCG" * 15


def runtime_available():
    try:
        return version("ostir") == "1.1.3" and version("ViennaRNA") == "2.7.2"
    except PackageNotFoundError:
        return False


class RBSMeritTests(unittest.TestCase):
    def test_bounded_finite_and_no_reward_for_reversed_or_equal_contrast(self):
        for off in (-100., 0., 100.):
            for on in (-100., 0., 100.):
                result = rbs.compute_rbs_merit(off, on)
                self.assertTrue(0 <= result["merit"] <= 1)
                self.assertTrue(all(math.isfinite(result[key]) for key in (
                    "off_rate", "on_rate", "merit", "log10_on_off")))
                if on >= off:
                    self.assertEqual(result["merit"], 0.)
                json.dumps(result, allow_nan=False)

    def test_ratio_alone_does_not_reward_two_inactive_states(self):
        weak = rbs.compute_rbs_merit(100., 80.)
        strong = rbs.compute_rbs_merit(0., -20.)
        self.assertAlmostEqual(weak["log10_on_off"], strong["log10_on_off"])
        self.assertLess(weak["merit"], 1e-10)
        self.assertGreater(strong["merit"], .99)

    def test_energy_zero_is_only_a_normalization_pivot(self):
        result = rbs.compute_rbs_merit(10., 0.)
        self.assertEqual(result["on_fraction"], .5)
        self.assertAlmostEqual(result["on_rate"], result["pivot_rate"])
        self.assertAlmostEqual(result["contrast"],
            (result["on_rate"] - result["off_rate"]) / (result["on_rate"] + result["off_rate"]))
        self.assertAlmostEqual(result["merit"], result["on_fraction"] * result["contrast"])

    def test_invalid_or_unrepresentable_rates_raise(self):
        for energy in (True, "3", math.nan, math.inf, -math.inf, -1e5, 1e5):
            with self.subTest(energy=energy), self.assertRaises(ValueError):
                rbs.compute_rbs_merit(energy, 0.)
            with self.subTest(on_energy=energy), self.assertRaises(ValueError):
                rbs.compute_rbs_merit(0., energy)


class RBSInputTests(unittest.TestCase):
    def test_invalid_sequence_site_or_mask_never_invokes_model(self):
        sequence = CLOSED_SEQUENCE
        aug = sequence.index("AUG")
        cases = [(sequence.replace("U", "T"), aug, ()), ("", 1, ()),
                 (sequence, True, ()), (sequence, aug + 1, ()),
                 (sequence, aug, [3]), (sequence, aug, (True,)),
                 (sequence, aug, (-1,)), (sequence, aug, (aug,)),
                 (sequence, aug, (3, 3)), (sequence, aug, (len(sequence),))]
        with patch.object(rbs, "_prediction") as calculate:
            for args in cases:
                with self.subTest(args=args), self.assertRaises(ValueError):
                    rbs.predict_rbs_switch(*args)
            calculate.assert_not_called()

    def test_no_ostir_result_is_not_a_zero_score(self):
        class EmptyFactory:
            def __init__(self, *args, **kwargs):
                self.results = []

            def calc_dG(self):
                pass

        with self.assertRaises(ValueError):
            rbs._model_result(EmptyFactory, PAPER_SEQUENCE, PAPER_SEQUENCE.rindex("AUG"), None, "none")

    def test_wrong_aug_result_is_rejected(self):
        class WrongFactory:
            def __init__(self, *args, **kwargs):
                self.results = [SimpleNamespace(results=lambda: {"start_position": 1, "start_codon": "AUG"})]

            def calc_dG(self):
                pass

        with self.assertRaisesRegex(ValueError, "different initiation site"):
            rbs._model_result(WrongFactory, PAPER_SEQUENCE, PAPER_SEQUENCE.rindex("AUG"), None, "none")


@unittest.skipUnless(runtime_available(), "Pinned OSTIR/ViennaRNA runtime unavailable")
class RBSNativeTests(unittest.TestCase):
    def setUp(self):
        rbs._prediction.cache_clear()

    def test_off_matches_ordinary_ostir_on_complete_sequence(self):
        from ostir import run_ostir
        for sequence in (PAPER_SEQUENCE, CLOSED_SEQUENCE):
            aug = sequence.rindex("AUG")
            predicted = rbs.predict_rbs_switch(sequence, aug, ())
            reference = run_ostir(sequence, start=aug + 1, end=aug + 1,
                                  threads=1, decimal_places=8)[0]
            self.assertEqual(predicted["off"]["dG_total"], reference["dG_total"])
            self.assertAlmostEqual(predicted["off_rate"], reference["expression"], places=7)
            self.assertEqual(predicted["reference"]["dangles"], "none" if aug > 35 else "all")

    def test_exact_window_coordinates_and_only_requested_bases_are_masked(self):
        aug = PAPER_SEQUENCE.rindex("AUG")
        occupied = (2, 15, 16, 20)
        result = rbs.predict_rbs_switch(PAPER_SEQUENCE, aug, occupied)
        start, end = max(0, aug - 35), min(len(PAPER_SEQUENCE), aug + 35)
        self.assertEqual((result["window_start"], result["window_end"]), (start, end))
        self.assertEqual(result["window_sequence"], PAPER_SEQUENCE[start:end])
        self.assertEqual(result["sequence"], PAPER_SEQUENCE)
        self.assertEqual(result["occupied_indices"], list(occupied))
        self.assertEqual(result["window_aug_start"], aug - start)
        self.assertEqual(result["window_occupied_indices"], [i - start for i in occupied if start <= i < end])
        for local, symbol in enumerate(result["on_constraints"]):
            self.assertEqual(symbol, "x" if local + start in occupied else ".")
        self.assertEqual(result["on_constraints"][aug - start:aug - start + 3], "...")
        self.assertEqual(result["reference"]["sodium_m"], 1.021)
        self.assertEqual(result["reference"]["celsius"], 37.)
        json.dumps(result, allow_nan=False)

    def test_no_occupancy_produces_identical_on_off(self):
        aug = PAPER_SEQUENCE.rindex("AUG")
        for occupied in ((), (0, 1, 2)):
            # Sites entirely outside the model window cannot modify its rate.
            result = rbs.predict_rbs_switch(PAPER_SEQUENCE, aug, occupied)
            self.assertEqual(result["off"], result["on"])
            self.assertEqual(result["off_rate"], result["on_rate"])
            self.assertEqual(result["merit"], 0.)

    def test_occupying_the_competitor_releases_the_rbs_in_synthetic_control(self):
        aug = CLOSED_SEQUENCE.index("AUG")
        result = rbs.predict_rbs_switch(CLOSED_SEQUENCE, aug, tuple(range(3, 9)))
        self.assertGreater(result["on_rate"], result["off_rate"])
        self.assertGreater(result["merit"], 0.)
        self.assertEqual(result["sequence"], CLOSED_SEQUENCE)

    def test_dna_table_interleaving_cannot_contaminate_cached_ostir_details(self):
        import RNA
        aug = CLOSED_SEQUENCE.index("AUG")
        args = (CLOSED_SEQUENCE, aug, tuple(range(3, 9)))
        reference = rbs.predict_rbs_switch(*args)
        RNA.params_load_DNA_Mathews2004()
        rbs._prediction.cache_clear()
        self.assertEqual(rbs.predict_rbs_switch(*args), reference)

    def test_cache_is_bounded_and_return_values_are_independent(self):
        aug = CLOSED_SEQUENCE.index("AUG")
        with patch.object(rbs, "_model_result", wraps=rbs._model_result) as calculate:
            first = rbs.predict_rbs_switch(CLOSED_SEQUENCE, aug, (3, 4))
            first["merit"] = -100
            first["off"]["dG_total"] = 123456
            second = rbs.predict_rbs_switch(CLOSED_SEQUENCE, aug, (4, 3))
            self.assertEqual(calculate.call_count, 2)
            self.assertGreaterEqual(second["merit"], 0)
            self.assertNotEqual(second["off"]["dG_total"], first["off"]["dG_total"])
        self.assertEqual(rbs._prediction.cache_info().maxsize, 256)

    def test_runtime_version_mismatch_is_explicit(self):
        rbs._runtime.cache_clear()
        try:
            with patch.object(rbs, "version", return_value="9.9.9"), self.assertRaisesRegex(RuntimeError, "requires OSTIR"):
                rbs._runtime()
        finally:
            rbs._runtime.cache_clear()

    def test_saved_validation_uses_no_native_runtime(self):
        aug = CLOSED_SEQUENCE.index("AUG")
        occupied = tuple(range(3, 9))
        result = rbs.predict_rbs_switch(CLOSED_SEQUENCE, aug, occupied)
        with patch.object(rbs, "_runtime", side_effect=AssertionError("must not load OSTIR")):
            self.assertIsNone(rbs.validate_rbs_prediction(result, CLOSED_SEQUENCE, aug, occupied))

    def test_modified_recorded_context_rates_and_merit_are_rejected(self):
        aug = CLOSED_SEQUENCE.index("AUG")
        occupied = tuple(range(3, 9))
        result = rbs.predict_rbs_switch(CLOSED_SEQUENCE, aug, occupied)
        changes = [(("window_start",), 1), (("aug_start",), aug + 1),
                   (("on_constraints",), "." * len(result["on_constraints"])),
                   (("off_rate",), 0.), (("on_rate",), float("nan")),
                   (("merit",), 1.), (("formula",), "on/off"),
                   (("reference", "sodium_m"), .12),
                   (("reference", "magnesium_m"), False),
                   (("reference", "ostir_version"), "1.0"),
                   (("reference", "parameter_sha256"), "a" * 64),
                   (("off", "dG_mRNA"), 0.)]
        with patch.object(rbs, "_runtime", side_effect=AssertionError("must not load OSTIR")):
            for path, value in changes:
                altered = deepcopy(result)
                target = altered
                for part in path[:-1]:
                    target = target[part]
                target[path[-1]] = value
                with self.subTest(path=path), self.assertRaises(ValueError):
                    rbs.validate_rbs_prediction(altered, CLOSED_SEQUENCE, aug, occupied)
            with self.assertRaises(ValueError):
                rbs.validate_rbs_prediction({**result, "surprise": 1}, CLOSED_SEQUENCE, aug, occupied)


if __name__ == "__main__":
    unittest.main()
