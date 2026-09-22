from __future__ import annotations

import sys
import tempfile
import unittest
import xml.etree.ElementTree as ET
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from aptaswitch_core.extension_visuals import (
    export_extension_visuals,
    extension_selection_landscape_svg,
    reference_pair_probability_svg,
    render_extension_structure_svg,
)

NS = {"s": "http://www.w3.org/2000/svg"}


class ExtensionVisualTests(unittest.TestCase):
    def setUp(self):
        from aptaswitch_core.localization import language_context
        self.enterContext(language_context("en"))

    def test_structure_preserves_actual_sequence_and_rejects_missing_data(self):
        svg = render_extension_structure_svg("AAGC", "GCTT", "((((+))))", title="Test <référence>")
        root = ET.fromstring(svg)
        bases = [node.text for node in root.findall(".//s:text", NS) if node.text in "ACGT"]
        self.assertEqual(bases, list("AAGCGCTT"))
        self.assertIn("Test &lt;référence&gt;", svg)
        with self.assertRaisesRegex(ValueError, "strand lengths"):
            render_extension_structure_svg("AAGC", "GCTT", "(((+)))")
        with self.assertRaisesRegex(ValueError, "Unbalanced"):
            render_extension_structure_svg("AAGC", "GCTT", "((((+))).")

    def test_probability_colors_require_one_true_probability_per_base(self):
        from aptaswitch_core.localization import language_context
        self.enterContext(language_context("fr"))
        svg = render_extension_structure_svg("GC", "GC", "((+))", nucleotide_probabilities=[0.2, 0.9, 0.9, 0.2])
        ET.fromstring(svg)
        self.assertIn("P(paire MFE)", svg)
        for probabilities in ([0.2], [0.2, 1.1, 0.9, 0.2], [0.2, float("nan"), 0.9, 0.2]):
            with self.assertRaises(ValueError):
                render_extension_structure_svg("GC", "GC", "((+))", nucleotide_probabilities=probabilities)

    def test_landscape_plots_sums_and_highlights_minimum_score(self):
        candidates = [
            {"rank": 2, "trigger_length": 30, "extension_5prime": "AAA", "final_score": 11, "extension_pair_probability_sum": 0.2, "extension_to_aptamer_probability_sum": 0.15, "reference_pair_probability_loss": 0.012},
            {"rank": 1, "trigger_length": 30, "extension_5prime": "AAT", "final_score": 7, "extension_pair_probability_sum": 0.1, "extension_to_aptamer_probability_sum": 0.05, "reference_pair_probability_loss": 0.008},
        ]
        root = ET.fromstring(extension_selection_landscape_svg(candidates, fraction=1.0))
        points = root.findall('.//s:circle[@class="candidate-point"]', NS)
        self.assertEqual([float(point.attrib["data-x"]) for point in points], [0.15, 0.35])
        self.assertEqual([float(point.attrib["data-y"]) for point in points], [0.008, 0.012])
        self.assertEqual([float(point.attrib["data-score"]) for point in points], [7, 11])
        self.assertEqual(len(root.findall('.//s:path[@class="best-candidate"]', NS)), 1)

    def test_degenerate_landscapes_are_valid(self):
        row = {"final_score": 0, "extension_pair_probability_sum": 0, "extension_to_aptamer_probability_sum": 0, "reference_pair_probability_loss": 0}
        for rows in ([], [row], [row, row]):
            svg = extension_selection_landscape_svg(rows)
            ET.fromstring(svg)
            self.assertNotIn("nan", svg.lower())

    def test_reference_probabilities_preserve_signed_deltas_and_pair_labels(self):
        pairs = [
            {"aptamer_index": 3, "trigger_index": 2, "reference_probability": 0.6, "candidate_probability": 0.5},
            {"reference_position_a": 1, "reference_position_b": 5, "aptamer_index": None, "trigger_index": None, "reference_probability": 0.2, "candidate_probability": 0.25},
        ]
        root = ET.fromstring(reference_pair_probability_svg(pairs))
        deltas = root.findall('.//s:circle[@class="probability-delta"]', NS)
        self.assertAlmostEqual(float(deltas[0].attrib["data-delta"]), -0.1)
        self.assertAlmostEqual(float(deltas[1].attrib["data-delta"]), 0.05)
        labels = [element.text for element in root.findall(".//s:text", NS)]
        self.assertIn("A3 – T2", labels)
        self.assertIn("1 – 5", labels)
        ET.fromstring(reference_pair_probability_svg([]))

    def test_exports_real_reference_and_best_candidate_for_each_length(self):
        pair = {"aptamer_index": 1, "trigger_index": 2, "reference_probability": 0.6, "candidate_probability": 0.55}
        result = {
            "input": {"aptamer_dna": "GC", "trigger_dna": "GC"},
            "reference": {"mfe_structure": "((+))", "mfe_energy_kcal_mol": -2.0},
            "candidates": [
                {"extended_trigger": "AAGC", "extension_5prime": "AA", "extension_length": 2, "trigger_length": 4, "mfe_structure": "((+..))", "mfe_energy_kcal_mol": -2.1, "final_score": 4, "extension_pair_probability_sum": 0.1, "extension_to_aptamer_probability_sum": 0.02, "reference_pair_probability_loss": 0.05, "reference_pairs": [pair]},
                {"extended_trigger": "AGC", "extension_5prime": "A", "extension_length": 1, "trigger_length": 3, "mfe_structure": "((+.))", "mfe_energy_kcal_mol": -2.2, "final_score": 3, "extension_pair_probability_sum": 0.08, "extension_to_aptamer_probability_sum": 0.01, "reference_pair_probability_loss": 0.04, "reference_pairs": [pair]},
            ],
        }
        with tempfile.TemporaryDirectory() as folder:
            outputs = export_extension_visuals(result, folder)
            self.assertEqual(set(outputs), {"svg_reference", "svg_landscape", "svg_landscape_3", "svg_landscape_4", "svg_comparison", "svg_candidate_3", "svg_candidate_4", "svg_probabilities_3", "svg_probabilities_4"})
            for filename in outputs.values():
                ET.parse(filename)
            comparison = ET.parse(outputs["svg_comparison"])
            self.assertEqual(len(comparison.findall('.//s:g[@class="structure-panel"]', NS)), 3)
            for length in (3, 4):
                landscape = ET.parse(outputs[f"svg_landscape_{length}"])
                self.assertEqual(len(landscape.findall('.//s:circle[@class="candidate-point"]', NS)), 1)
            result["reference"] = None
            self.assertEqual(export_extension_visuals(result, folder), {})


if __name__ == "__main__":
    unittest.main()
