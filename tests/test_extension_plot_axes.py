from __future__ import annotations

import copy
import math
import re
import sys
import unittest
import xml.etree.ElementTree as ET
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from aptaswitch_core.extension_visuals import (
    extension_selection_landscape_svg,
    reference_pair_probability_svg,
)

NS = {"s": "http://www.w3.org/2000/svg"}


def candidate(x: float, y: float, score: float, rank: int = 1) -> dict:
    return {
        "rank": rank,
        "trigger_length": 24,
        "extension_5prime": "AACC",
        "extension_3prime": "",
        "extension_pair_probability_sum": x,
        "extension_to_aptamer_probability_sum": 0.0,
        "reference_pair_probability_loss": y,
        "final_score": score,
    }


def elements(root: ET.Element, tag: str, class_name: str) -> list[ET.Element]:
    return root.findall(f'.//s:{tag}[@class="{class_name}"]', NS)


def span(nodes: list[ET.Element], attribute: str) -> float:
    coordinates = [float(node.attrib[attribute]) for node in nodes]
    return max(coordinates) - min(coordinates)


class ExtensionPlotAxesTests(unittest.TestCase):
    def setUp(self):
        from aptaswitch_core.localization import language_context
        self.enterContext(language_context("fr"))

    def _read_compact_axis(self, root: ET.Element, name: str) -> tuple[list[ET.Element], list[float]]:
        ticks = elements(root, "text", f"{name}-tick")
        labels = [node.text for node in ticks]
        self.assertGreaterEqual(len(labels), 3)
        self.assertEqual(len(labels), len(set(labels)))
        self.assertTrue(all(len(label) <= 9 for label in labels), labels)
        note = " ".join(node.text for node in elements(root, "text", f"{name}-scale"))
        scientific = re.search(r"× 10([⁻⁰¹²³⁴⁵⁶⁷⁸⁹]+)", note)
        self.assertIsNotNone(scientific, note)
        exponent = int(scientific.group(1).translate(str.maketrans("⁻⁰¹²³⁴⁵⁶⁷⁸⁹", "-0123456789")))
        offset = 0.0
        if "+" in note:
            offset = float(note.split("+", 1)[0].strip())
            self.assertRegex(note, r"grad(?:uation|\.)")
        return ticks, [offset + float(label) * 10 ** exponent for label in labels]

    def test_close_values_fill_landscape_at_different_numerical_scales(self):
        # The first case needs a zoom around a nonzero value; the others catch
        # a fixed small-range floor and the mostly empty y axis in real runs.
        cases = (
            (0.012345, 1e-10, 0.0008765, 2e-11, 18.912, 1e-7),
            (1e-10, 1e-12, 4e-12, 2e-14, 1e-10, 2e-13),
            (0.0018, 0.005, 0.0, 0.00013, 17.23, 0.423),
        )
        width, height = 1000, 640
        plot_width, plot_height = width - 150 - 115, height - 140 - 85
        for x, dx, y, dy, score, ds in cases:
            with self.subTest(x=x, y=y):
                rows = [candidate(x + i * dx, y + i * dy, score + i * ds, i + 1) for i in range(4)]
                root = ET.fromstring(extension_selection_landscape_svg(rows, width=width, height=height, fraction=1.0))
                points = elements(root, "circle", "candidate-point")
                self.assertEqual(len(points), len(rows))
                self.assertGreater(span(points, "cx"), 0.65 * plot_width)
                self.assertGreater(span(points, "cy"), 0.65 * plot_height)
                for point, row in zip(points, rows):
                    self.assertTrue(math.isclose(float(point.attrib["data-x"]), row["extension_pair_probability_sum"], rel_tol=1e-13))
                    self.assertTrue(math.isclose(float(point.attrib["data-y"]), row["reference_pair_probability_loss"], rel_tol=1e-13))

    def test_ticks_distinguish_close_axis_values_and_scores(self):
        rows = [candidate(0.012345 + i * 1e-10, 0.0008765 + i * 2e-11, 18.912 + i * 1e-7) for i in range(4)]
        root = ET.fromstring(extension_selection_landscape_svg(rows, fraction=1.0))
        for class_name in ("x-tick", "y-tick", "score-tick"):
            with self.subTest(axis=class_name):
                labels = [node.text for node in elements(root, "text", class_name)]
                self.assertGreaterEqual(len(labels), 3)
                self.assertEqual(len(labels), len(set(labels)), labels)

    def test_compact_landscape_ticks_show_units_and_reconstruct_true_values(self):
        cases = (
            (0.012345, 1e-12, 0.0008765, 2e-12, 18.912, 1e-9),
            (1e-10, 1e-12, 4e-12, 2e-14, 1e-10, 2e-13),
        )
        for x, dx, y, dy, score, ds in cases:
            root = ET.fromstring(extension_selection_landscape_svg([
                candidate(x + i * dx, y + i * dy, score + i * ds) for i in range(4)
            ], fraction=1.0))
            for axis in ("x", "y", "score"):
                with self.subTest(x=x, axis=axis):
                    ticks, decoded = self._read_compact_axis(root, axis)
                    actual = [float(tick.attrib["data-value"]) for tick in ticks]
                    # Visible graduations plus the displayed offset/unit must
                    # recover the data, within normal tick-label rounding.
                    tolerance = (max(actual) - min(actual)) * 0.005 + max(math.ulp(value) for value in actual)
                    for value, expected in zip(decoded, actual):
                        self.assertAlmostEqual(value, expected, delta=tolerance)

    def test_compact_reference_ticks_keep_probability_and_delta_units_explicit(self):
        reference, difference = 0.92, 1e-12
        rows = [
            {"aptamer_index": 1, "trigger_index": 2, "reference_probability": reference, "candidate_probability": reference + difference},
            {"aptamer_index": 2, "trigger_index": 1, "reference_probability": reference, "candidate_probability": reference - difference},
        ]
        root = ET.fromstring(reference_pair_probability_svg(rows))
        for axis, values in (
            ("probability", [row["candidate_probability"] for row in rows]),
            ("delta", [row["candidate_probability"] - reference for row in rows]),
        ):
            with self.subTest(axis=axis):
                _, decoded = self._read_compact_axis(root, axis)
                self.assertLess(min(decoded), min(values))
                self.assertGreater(max(decoded), max(values))
                self.assertLess(max(decoded) - min(decoded), 1.3 * (max(values) - min(values)))

    def test_identical_and_zero_values_have_finite_axes_and_no_jitter(self):
        for value in (0.0, 1e-12, 0.38):
            for count in (1, 3):
                with self.subTest(value=value, count=count):
                    rows = [candidate(value, value, value)] * count
                    root = ET.fromstring(extension_selection_landscape_svg(rows, fraction=1.0))
                    points = elements(root, "circle", "candidate-point")
                    self.assertEqual(len(points), count)
                    self.assertEqual(span(points, "cx"), 0)
                    self.assertEqual(span(points, "cy"), 0)
                    for point in points:
                        cx, cy = float(point.attrib["cx"]), float(point.attrib["cy"])
                        self.assertTrue(math.isfinite(cx) and math.isfinite(cy))
                        self.assertGreaterEqual(cx, 115)
                        self.assertLessEqual(cx, 850)
                        self.assertGreaterEqual(cy, 85)
                        self.assertLessEqual(cy, 500)
                    for class_name in ("x-tick", "y-tick"):
                        labels = [node.text for node in elements(root, "text", class_name)]
                        self.assertGreaterEqual(len(labels), 3)
                        self.assertEqual(len(labels), len(set(labels)))

    def test_autoscaling_keeps_every_candidate_and_highlights_actual_best(self):
        rows = [candidate(0.002 + (i % 7) * 1e-9, 0.0001 + (i % 5) * 1e-10, 20 - i * 1e-8, i + 1) for i in range(100)]
        original = copy.deepcopy(rows)
        root = ET.fromstring(extension_selection_landscape_svg(rows, fraction=1.0))
        points = elements(root, "circle", "candidate-point")
        self.assertEqual(len(points), 100)
        self.assertEqual(rows, original)
        expected = sorted(rows, key=lambda row: row["final_score"])
        for point, row in zip(points, expected):
            for attribute, field in (("data-x", "extension_pair_probability_sum"), ("data-y", "reference_pair_probability_loss"), ("data-score", "final_score")):
                self.assertTrue(math.isclose(float(point.attrib[attribute]), row[field], rel_tol=1e-13))
        best = elements(root, "path", "best-candidate")
        self.assertEqual(len(best), 1)
        diamond = re.findall(r"[-+]?(?:\d*\.)?\d+", best[0].attrib["d"])
        self.assertAlmostEqual(float(diamond[0]), float(points[0].attrib["cx"]), places=2)
        self.assertAlmostEqual(float(diamond[1]) + 11, float(points[0].attrib["cy"]), places=2)
        # Repeated measured coordinates remain coincident rather than getting
        # synthetic offsets to make the picture appear more diverse.
        locations = {}
        for point in points:
            key = point.attrib["data-x"], point.attrib["data-y"]
            position = point.attrib["cx"], point.attrib["cy"]
            self.assertEqual(locations.setdefault(key, position), position)

    def test_default_landscape_selects_best_ten_percent_and_fits_only_them(self):
        for count in (0, 1, 11, 100, 1001):
            with self.subTest(count=count):
                expected_count = math.ceil(count * 0.1)
                rows = [
                    candidate(
                        0.002 + i * 1e-10 if i < expected_count else 1000 + i,
                        0.0001 + i * 1e-10 if i < expected_count else 1000 + i,
                        float(i), i + 1,
                    )
                    for i in range(count)
                ]
                root = ET.fromstring(extension_selection_landscape_svg(list(reversed(rows))))
                points = elements(root, "circle", "candidate-point")
                self.assertEqual(len(points), expected_count)
                self.assertEqual([float(point.attrib["data-score"]) for point in points], list(range(expected_count)))
                for point, row in zip(points, rows[:expected_count]):
                    self.assertTrue(math.isclose(float(point.attrib["data-x"]), row["extension_pair_probability_sum"], rel_tol=1e-13))
                    self.assertTrue(math.isclose(float(point.attrib["data-y"]), row["reference_pair_probability_loss"], rel_tol=1e-13))
                if expected_count > 1:
                    self.assertGreater(span(points, "cx"), 0.65 * (850 - 115))
                    self.assertGreater(span(points, "cy"), 0.65 * (500 - 85))
                    self.assertLess(max(float(tick.attrib["data-value"]) for tick in elements(root, "text", "x-tick")), 0.003)
                    self.assertLess(max(float(tick.attrib["data-value"]) for tick in elements(root, "text", "y-tick")), 0.0002)

    def test_tiny_reference_probability_differences_remain_visible(self):
        for difference in (1e-8, 1e-12):
            with self.subTest(difference=difference):
                rows = [
                    {"aptamer_index": 1, "trigger_index": 2, "reference_probability": 0.92, "candidate_probability": 0.92 + difference},
                    {"aptamer_index": 2, "trigger_index": 1, "reference_probability": 0.92, "candidate_probability": 0.92 - difference},
                ]
                root = ET.fromstring(reference_pair_probability_svg(rows))
                reference = elements(root, "circle", "reference-probability")
                candidates = elements(root, "circle", "candidate-probability")
                deltas = elements(root, "circle", "probability-delta")
                self.assertGreater(span(reference + candidates, "cx"), 0.65 * (630 - 118))
                self.assertGreater(span(deltas, "cx"), 0.65 * (955 - 730))
                self.assertEqual(len(reference), len(rows))
                self.assertEqual(len(candidates), len(rows))
                for delta, row in zip(deltas, rows):
                    expected = row["candidate_probability"] - row["reference_probability"]
                    self.assertTrue(math.isclose(float(delta.attrib["data-delta"]), expected, rel_tol=1e-13))

    def test_constant_reference_probabilities_at_bounds_are_finite(self):
        for probability in (0.0, 1e-12, 0.92, 1.0):
            with self.subTest(probability=probability):
                rows = [{"aptamer_index": 1, "trigger_index": 1, "reference_probability": probability, "candidate_probability": probability}]
                root = ET.fromstring(reference_pair_probability_svg(rows))
                reference = elements(root, "circle", "reference-probability")[0]
                candidate_point = elements(root, "circle", "candidate-probability")[0]
                delta = elements(root, "circle", "probability-delta")[0]
                self.assertEqual(reference.attrib["cx"], candidate_point.attrib["cx"])
                self.assertEqual(float(delta.attrib["data-delta"]), 0.0)
                for point in (reference, candidate_point, delta):
                    self.assertTrue(math.isfinite(float(point.attrib["cx"])))
                    self.assertTrue(math.isfinite(float(point.attrib["cy"])))
                self.assertGreaterEqual(float(reference.attrib["cx"]), 118)
                self.assertLessEqual(float(reference.attrib["cx"]), 630)


if __name__ == "__main__":
    unittest.main()
