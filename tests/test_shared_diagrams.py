"""Keep the original tSwitch drawings authoritative for every shared diagram."""
from __future__ import annotations

import hashlib
import unittest
from unittest.mock import patch
import xml.etree.ElementTree as ET

from aptaswitch_core import preview
from aptaswitch_core.architecture import default_architecture_for_trigger
from aptaswitch_core.extension_preview import aptamer_trigger_linear_svg, extension_preview_svg
from aptaswitch_core.extension_visuals import render_extension_structure_svg
from aptaswitch_core.models import SwitchArchitecture

NS = {"s": "http://www.w3.org/2000/svg"}
TRIGGER = "CTCTCCTTACGCCACCCACACCCG"


class SharedDiagramTests(unittest.TestCase):
    def setUp(self):
        from aptaswitch_core.localization import language_context
        self.enterContext(language_context("fr"))

    def test_original_switch_svg_output_is_preserved_exactly(self):
        # Captured from the original tSwitch functions before extracting the
        # common renderer: geometry, ordering, typography and colors stay exact.
        # Explicit historical settings retain the renderer regression check
        # even when the standard activator architecture changes.
        architecture = SwitchArchitecture(12, 18, upper_stem_length=5, upper_stem2_length=4,
                                          rbs_loop="UAGAGGAGAAC", rbs_prefix_length=0, frame_linker="CC")
        legacy_trigger = "CCGATGCTCTCCTTACGCCACCCACACCCG"
        outputs = (
            (preview.linear_sequence_svg(architecture, legacy_trigger)[0], "4d2c373cd821ade0d624b4c7886c1e14745047c2931b144740b296d1fc696b2e"),
            (preview.switch_preview_svg(architecture, legacy_trigger), "0dc3b105649c5cc4564aa4ad86496578d2ff83339f2e7f10092ba3a846d23ff7"),
        )
        for svg, expected in outputs:
            self.assertEqual(hashlib.sha256(svg.encode()).hexdigest(), expected)

    def test_all_linear_views_use_the_same_renderer_and_box_geometry(self):
        with patch.object(preview, "render_linear_svg", wraps=preview.render_linear_svg) as shared:
            architecture = default_architecture_for_trigger(len(TRIGGER))
            outputs = (
                preview.linear_sequence_svg(architecture, TRIGGER),
                extension_preview_svg("ACGTAC", "GCTA", 9, "both", [2]),
                aptamer_trigger_linear_svg("ACGTAC", "GCTA", ligand_aptamer_positions=[2]),
            )
        self.assertEqual(shared.call_count, 3)
        for svg, width, height in outputs:
            root = ET.fromstring(svg)
            self.assertEqual(root.attrib["viewBox"], f"0 0 {width} {height}")
            boxes = [box for box in root.findall('.//s:rect', NS) if box.attrib.get('rx') == '3']
            self.assertTrue(boxes)
            self.assertTrue(all((box.attrib['width'], box.attrib['height']) == ('20.0', '26.0') for box in boxes))
            self.assertEqual(float(boxes[1].attrib['x']) - float(boxes[0].attrib['x']), 22)

    def test_actual_switch_sequence_replaces_only_placeholder_letters(self):
        architecture = default_architecture_for_trigger(len(TRIGGER))
        model = preview.build_switch_model(architecture, TRIGGER)
        resolved = model.switch_sequence.replace('N', 'A')
        svg, _, _ = preview.linear_sequence_svg(architecture, TRIGGER, switch_sequence=resolved)
        root = ET.fromstring(svg)
        bases = [node.text for node in root.findall('.//s:text', NS) if node.attrib.get('font-size') == '13']
        self.assertEqual(''.join(bases), model.trigger_sequence[::-1] + resolved)
        for invalid in (resolved[:-1], model.switch_sequence):
            with self.assertRaises(ValueError):
                preview.linear_sequence_svg(architecture, TRIGGER, switch_sequence=invalid)

    def test_measured_linear_sequences_preserve_chemistry_order_and_local_annotations(self):
        svg, _, _ = aptamer_trigger_linear_svg(
            'aagc u', 'uagc uu', extension_length=3, extension_3prime_length=2,
            ligand_aptamer_positions=[2, 5],
        )
        root = ET.fromstring(svg)
        groups = root.findall('.//s:g', NS)
        sequence = ''.join(g.find('s:text', NS).text for g in groups)
        self.assertEqual(sequence, 'AAGCUUAGCUU')
        self.assertEqual([g.attrib['data-position'] for g in groups if g.attrib['data-domain'] == 'ligand'], ['2', '5'])
        self.assertEqual([g.find('s:text', NS).text for g in groups if g.attrib['data-domain'] == 'extension-5prime'], ['U'])
        self.assertEqual([g.find('s:text', NS).text for g in groups if g.attrib['data-domain'] == 'extension-3prime'], ['U', 'U'])
        core = [g for g in groups if g.attrib['data-domain'] == 'trigger-initial']
        self.assertEqual(''.join(g.find('s:text', NS).text for g in core), 'AGC')
        self.assertEqual([g.attrib['data-position'] for g in core], ['1', '2', '3'])
        self.assertTrue(all(g.find('s:rect', NS).attrib['fill'] == preview.TRIGGER_COLOR for g in core))
        self.assertFalse(root.findall('.//s:line[@stroke-width="1"]', NS), 'No aptamer-trigger pairing may be invented')

    def test_unresolved_or_inconsistent_measured_linear_data_is_rejected(self):
        for aptamer, trigger, kwargs in (
            ('ACGN', 'ACGT', {}), ('ACGT', '', {}), ('ACGT', 'ACGT', {'extension_length': 5}),
            ('ACGT', 'ACGT', {'extension_length': 1, 'extension_3prime_length': 2}),
            ('ACGT', 'ACGT', {'ligand_aptamer_positions': [5]}),
        ):
            with self.assertRaises(ValueError):
                aptamer_trigger_linear_svg(aptamer, trigger, **kwargs)

    def test_measured_2d_uses_one_canonical_svg_with_probability_and_ligand_annotations(self):
        with patch.object(preview, 'render_structure_svg', wraps=preview.render_structure_svg) as shared:
            svg = render_extension_structure_svg(
                'AAGC', 'GCUU', '((((+))))', title='Mesure réelle',
                nucleotide_probabilities=[.5] * 8, ligand_aptamer_positions=[1, 4],
            )
        shared.assert_called_once()
        self.assertEqual(shared.call_args.args[:2], ('AAGCGCUU', '((((+))))'))
        root = ET.fromstring(svg)
        self.assertFalse(root.findall('.//s:svg', NS), 'The 2D adapter must not wrap a second SVG')
        self.assertEqual(root.find('s:title', NS).text, 'Mesure réelle')
        rings = root.findall('.//s:circle[@class="base-annotation"]', NS)
        self.assertEqual([ring.attrib['data-position'] for ring in rings], ['1', '4'])
        self.assertIn('P(paire MFE)', svg)
        self.assertIn('annotation utilisateur', svg)
        self.assertEqual(len(shared.call_args.kwargs['probability_scale']), 100)

    def test_measured_2d_never_replaces_incomplete_data_with_placeholder_bases(self):
        for aptamer, trigger, structure in (
            ('AANG', 'GCTT', '((((+))))'), ('AAGC', 'GCTT', '(((+)))'),
            ('AAGC', 'GCTT', '((((+....'),
        ):
            with patch.object(preview, 'render_structure_svg') as shared:
                with self.assertRaises(ValueError):
                    render_extension_structure_svg(aptamer, trigger, structure)
                shared.assert_not_called()

    def test_live_structure_hides_export_annotations_but_preserves_molecular_content(self):
        options = dict(
            title="Titre export", subtitle="ΔG = -4.2 kcal/mol",
            footer_notes=[("Conditions de calcul", "#64748b")],
            description="Structure MFE : ((((+))))",
            legend=[("Aptamère", preview.APTAMER_COLOR)],
            base_outlines={1: preview.LIGAND_COLOR},
        )
        exported = preview.render_structure_svg('AAGCGCUU', '((((+))))', [preview.APTAMER_COLOR] * 8, **options)
        live = preview.render_structure_svg('AAGCGCUU', '((((+))))', [preview.APTAMER_COLOR] * 8,
                                            **options, show_annotations=False)
        for annotation in ('Titre export', 'ΔG = -4.2 kcal/mol', 'Conditions de calcul', '(((( ))))'):
            self.assertIn(annotation, exported)
            self.assertNotIn(annotation, live)
        self.assertNotIn('((((+))))', live)
        for svg in (exported, live):
            root = ET.fromstring(svg)
            text = [node.text for node in root.findall('.//s:text', NS)]
            self.assertEqual(''.join(value for value in text if value in 'ACGTU'), 'AAGCGCUU')
            self.assertEqual(text.count('5′'), 2)
            self.assertEqual(text.count('3′'), 2)
            self.assertIn('Aptamère', text)
            self.assertEqual(len(root.findall('.//s:circle[@class="base-annotation"]', NS)), 1)

    def test_measured_live_structure_keeps_probability_scale_and_ligand_positions(self):
        kwargs = dict(nucleotide_probabilities=[.5] * 8, ligand_aptamer_positions=[1, 4])
        exported = render_extension_structure_svg('AAGC', 'GCUU', '((((+))))', **kwargs)
        live = render_extension_structure_svg('AAGC', 'GCUU', '((((+))))', **kwargs, show_annotations=False)
        self.assertIn('P(paire MFE)', exported)
        self.assertNotIn('P(paire MFE)', live)
        self.assertNotIn('Structure NUPACK', live)
        self.assertNotIn('((((', live)
        root = ET.fromstring(live)
        self.assertEqual([node.attrib['data-position'] for node in root.findall('.//s:circle[@class="base-annotation"]', NS)], ['1', '4'])
        self.assertEqual(len(root.findall('.//s:rect[@height="10"]', NS)), 100)

    def test_linear_live_modes_remove_captions_without_removing_domains_or_bases(self):
        kwargs = dict(title="Séquences exportées", subtitle="ΔG = -4.2 kcal/mol", ligand_aptamer_positions=[2])
        exported, _, export_height = aptamer_trigger_linear_svg('AAGC', 'GCUU', **kwargs)
        live, _, live_height = aptamer_trigger_linear_svg('AAGC', 'GCUU', **kwargs, show_annotations=False)
        self.assertIn('ΔG = -4.2 kcal/mol', exported)
        self.assertNotIn('ΔG', live)
        self.assertNotIn('Séquences exportées', live)
        self.assertLess(live_height, export_height)
        root = ET.fromstring(live)
        groups = root.findall('.//s:g', NS)
        self.assertEqual(''.join(g.find('s:text', NS).text for g in groups), 'AAGCGCUU')
        self.assertEqual([g.attrib['data-position'] for g in groups if g.attrib['data-domain'] == 'ligand'], ['2'])
        self.assertIn('Trigger initial', live)
        extension, _, _ = extension_preview_svg('AAGC', 'GCUU', 8, 'both', [2], show_annotations=False)
        self.assertNotIn('Extension · trigger final', extension)
        extension_root = ET.fromstring(extension)
        self.assertEqual(sum(node.text == 'N' for node in extension_root.findall('.//s:text', NS)), 4)

    def test_switch_live_adapters_hide_the_selected_structure_not_its_nucleotides(self):
        architecture = default_architecture_for_trigger(len(TRIGGER))
        model = preview.build_switch_model(architecture, TRIGGER)
        for state, structure in (('off', model.off_structure), ('on', model.on_structure)):
            exported = preview.switch_preview_svg(architecture, TRIGGER, state=state)
            live = preview.switch_preview_svg(architecture, TRIGGER, state=state, show_annotations=False)
            self.assertIn(structure.replace('+', ' '), exported)
            self.assertNotIn(structure.replace('+', ' '), live)
            candidate = preview.candidate_structure_svg(
                architecture, model.trigger_sequence, model.switch_sequence.replace('N', 'A'), structure,
                state=state, show_annotations=False,
            )
            for svg in (live, candidate):
                self.assertNotIn('Switch ON', svg)
                self.assertNotIn('Switch OFF', svg)
                self.assertNotIn(structure.replace('+', ' '), svg)
                self.assertIn('<circle', svg)
        self.assertEqual(preview.linear_sequence_svg(architecture, TRIGGER),
                         preview.linear_sequence_svg(architecture, TRIGGER, show_annotations=False))

    def test_live_animation_uses_clean_frames_through_the_final_state(self):
        architecture = default_architecture_for_trigger(len(TRIGGER))
        kwargs = dict(approach_steps=2, tween_steps=1, hold_start=1, hold_end=1)
        exported = preview.switch_animation_frames(architecture, TRIGGER, **kwargs)
        live = preview.switch_animation_frames(architecture, TRIGGER, **kwargs, show_annotations=False)
        self.assertEqual(len(live), len(exported))
        self.assertIn('Fixation du trigger', exported[0])
        self.assertIn('Switch ON', exported[-1])
        for frame in live:
            root = ET.fromstring(frame)
            self.assertTrue(root.findall('.//s:circle', NS))
            self.assertNotIn('Fixation du trigger', frame)
            self.assertNotIn('Switch ON', frame)
            self.assertFalse(root.findall('.//s:text[@font-size="12"]', NS))


if __name__ == '__main__':
    unittest.main()
