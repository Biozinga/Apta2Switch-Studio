from __future__ import annotations

import unittest
import xml.etree.ElementTree as ET

from aptaswitch_core.extension_preview import extension_preview_svg, ligand_positions
from aptaswitch_core.extension_visuals import render_extension_structure_svg

NS = {"s": "http://www.w3.org/2000/svg"}


class ExtensionPreviewTests(unittest.TestCase):
    def test_each_side_preserves_core_and_ligand_positions(self):
        for side, left, right in (("5prime", 5, 0), ("3prime", 0, 5), ("both", 3, 2)):
            svg, width, height = extension_preview_svg("ACGTAC", "GCTA", 9, side, [1, 4, 6])
            root = ET.fromstring(svg)
            nodes = root.findall('.//s:g', NS)
            self.assertEqual([int(n.attrib['data-position']) for n in nodes if n.attrib['data-domain'] == 'ligand'], [1, 4, 6])
            self.assertEqual(sum(n.attrib['data-domain'] == 'extension-5prime' for n in nodes), left)
            self.assertEqual(sum(n.attrib['data-domain'] == 'extension-3prime' for n in nodes), right)
            sequence = ''.join(n.find('s:text', NS).text for n in nodes)
            self.assertEqual(sequence, "ACGTAC" + "N" * left + "GCTA" + "N" * right)
            self.assertEqual(int(root.attrib['height']), height)
            self.assertEqual(int(root.attrib['width']), width)

    def test_rna_and_large_requests_are_safe_to_preview(self):
        svg, _, _ = extension_preview_svg("a c u g", "ttag", 10**9, "both", [3], rna=True)
        root = ET.fromstring(svg)
        letters = [n.text for n in root.findall('.//s:g/s:text', NS)]
        self.assertIn('U', letters)
        self.assertNotIn('T', letters)
        self.assertLess(len(svg), 25000)
        self.assertIn('1000000000 nt', svg)

    def test_invalid_ligand_positions_rejected(self):
        for invalid in ([0], [5], [-1], [True], [1.5]):
            with self.assertRaises(ValueError):
                ligand_positions(invalid, 4)
        self.assertEqual(ligand_positions([4, 2, 4], 4), (2, 4))

    def test_result_annotations_do_not_replace_probability_colors(self):
        kwargs = dict(nucleotide_probabilities=[0.2] * 8)
        original = ET.fromstring(render_extension_structure_svg('AAGC', 'GCTT', '((((+))))', **kwargs))
        annotated = ET.fromstring(render_extension_structure_svg('AAGC', 'GCTT', '((((+))))', ligand_aptamer_positions=[1, 4], **kwargs))
        rings = annotated.findall('.//s:circle[@class="base-annotation"]', NS)
        self.assertEqual([n.attrib['data-position'] for n in rings], ['1', '4'])
        base_fills = lambda root: [n.attrib.get('fill') for n in root.findall('.//s:circle', NS) if n.attrib.get('class') != 'base-annotation']
        self.assertEqual(base_fills(original), base_fills(annotated))


if __name__ == '__main__':
    unittest.main()
