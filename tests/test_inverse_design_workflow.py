from __future__ import annotations

import copy
import tempfile
import unittest
from dataclasses import replace
from pathlib import Path
from unittest.mock import patch

from aptaswitch_core.architecture import default_architecture_for_trigger
from aptaswitch_core.models import InverseExtensionSettings, RunRequest, RunResult, SequenceInput
from aptaswitch_core.thermodynamics.design_constraints import BaseSite, build_design_problem
from aptaswitch_core.preview import build_switch_model
from aptaswitch_core.sequences import reverse_complement_rna_pattern
from aptaswitch_core.workflow import reanalyze_inverse_and_export
from aptaswitch_studio.run_library import load_run
from tests.test_inverse_extension import candidate
from tests.design_fixtures import synthetic_native_metadata
from tests.ui_helpers import FrenchAppTest

ROOT = Path(__file__).resolve().parents[1]


def request(side='5prime', target=7):
    return RunRequest(SequenceInput('ACGTACGT', 'AGTC'), default_architecture_for_trigger(target),
                      inverse_extension=InverseExtensionSettings(target, side, 2), trials=2,
                      thermodynamics_library_path='/fixture/native')


class PatternAndConstraintsTests(unittest.TestCase):
    def setUp(self):
        from tests.native_ui_fixtures import patch_native_runtime
        patch_native_runtime(self)

    def test_job_eta_resets_for_comparison_after_design(self):
        from aptaswitch_studio.web_runtime import DesignRunJob
        with patch('aptaswitch_studio.web_runtime.time.monotonic', side_effect=[100., 200., 210.]):
            job = DesignRunJob(request())
            job._progress(2, 2, 'Design complete')
            job._progress(0, 3, 'Comparison')
            job._progress(1, 3, 'Reference complete')
            self.assertEqual(job.snapshot().eta_seconds, 20)

    def test_patterns_preserve_entire_original_and_place_odd_extra_base_at_five_prime(self):
        for side, pattern in [('5prime', 'XXXAGTC'), ('3prime', 'AGTCXXX'), ('both', 'XXAGTCX')]:
            with self.subTest(side=side):
                req = request(side).normalized()
                self.assertEqual(req.sequences.trigger_dna, 'AGTC')
                self.assertEqual(req.design_trigger_pattern(), pattern)
                model = build_switch_model(req.architecture, pattern)
                self.assertEqual(model.trigger_sequence, pattern)
                self.assertEqual(model.switch_sequence[3:10], reverse_complement_rna_pattern(pattern))
                self.assertEqual(model.switch_sequence.count('GACU'), 1)

    def test_unknown_bases_are_only_allowed_in_internal_patterns(self):
        with self.assertRaises(ValueError):
            replace(request(), sequences=SequenceInput('ACGT', 'XAGTC')).normalized()
        with self.assertRaises(ValueError):
            replace(request(), sequences=SequenceInput('ACGT', 'AGTC', binding_start=1)).normalized()
        with self.assertRaises(ValueError):
            replace(request(), inverse_extension=InverseExtensionSettings(4)).normalized()
        with self.assertRaises(ValueError):
            replace(request(), thermo=replace(request().thermo, trigger_material='rna')).normalized()

    def test_libraries_bind_each_free_trigger_position_to_reverse_switch_position(self):
        for side in ('5prime', '3prime', 'both'):
            req = request(side, 30).normalized()
            design = build_design_problem(req)
            pattern = req.design_trigger_pattern()
            strands = {strand.name: strand for strand in design.strands}
            start = design.layout['len_t7']
            self.assertEqual(strands['Trigger'].pattern, pattern)
            self.assertEqual(strands['Switch'].pattern[start:start+30], reverse_complement_rna_pattern(pattern))
            constraints = [item for item in design.libraries if item.name.startswith('switch_on')]
            self.assertEqual(len(constraints), 30)
            free_sites = [i for i, base in enumerate(pattern) if base == 'X']
            self.assertEqual(len(free_sites), 26)
            for i in free_sites:
                constraint = next(item for item in constraints if BaseSite('Trigger', i) in item.sites)
                self.assertEqual(constraint.sites, (BaseSite('Trigger', i), BaseSite('Switch', start+29-i)))
                self.assertEqual(set(constraint.allowed_words), {'AU', 'CG', 'GC', 'TA'})

    def test_recomparison_exports_new_run_without_redesign_or_mutating_original(self):
        with tempfile.TemporaryDirectory() as tmp:
            req = request(target=6)
            # This synthetic candidate has no OSTIR measurement to score.
            req = replace(req, output_dir=Path(tmp), trials=1,
                          scoring=replace(req.scoring, rbs_calculator=0.))
            source = RunResult('original', req, [candidate(1, 'AAAGTC')], Path(tmp)/'original', 'start', 'end',
                               **synthetic_native_metadata(req))
            before = copy.deepcopy(source.to_dict())
            def analyze(result, *args, **kwargs):
                result.inverse_extension_analysis = {'status':'completed', 'model':{'material':'experimental-Sugimoto95-Vienna04-loopmean'}}
                result.candidates[0].inverse_extension_metrics = {
                    'status':'completed', 'final_score':0.5, 'extended_trigger':'AAAGTC',
                    'extension_5prime':'AA','extension_3prime':'',
                    'reference_mfe_structure':'........+....','mfe_structure':'........+......'}
            with patch('aptaswitch_core.workflow.run_design', side_effect=AssertionError('Must not redesign')), \
                 patch('aptaswitch_core.inverse_extension.analyze_inverse_extensions', side_effect=analyze):
                result, exports = reanalyze_inverse_and_export(source, req)
            self.assertEqual(source.to_dict(), before)
            self.assertNotEqual(result.run_id, source.run_id)
            self.assertTrue(Path(exports['inverse_csv']).is_file())
            restored = load_run(Path(exports['json'])).result
            self.assertEqual(restored.request.inverse_extension, req.inverse_extension)
            self.assertEqual(restored.inverse_extension_analysis['source_run_id'], 'original')
            self.assertFalse(restored.inverse_extension_analysis['switch_design_repeated'])
            self.assertEqual(restored.candidates[0].inverse_extension_metrics['final_score'], .5)


class InverseFullAppTests(unittest.TestCase):
    def setUp(self):
        from tests.native_ui_fixtures import patch_native_runtime
        patch_native_runtime(self)

    def test_short_original_can_design_standard_switch_with_free_extension(self):
        with patch('aptaswitch_core.settings.load_settings', return_value={}), \
             patch('aptaswitch_core.settings.save_settings'), \
             patch('aptaswitch_core.thermodynamics.setup.configured_library_path', return_value='/fixture/native'), \
             patch('aptaswitch_studio.extension_ui.configured_library_path', return_value='/fixture/native'):
            app = FrenchAppTest.from_file(str(ROOT/'src/aptaswitch_studio/web_app.py'), default_timeout=30)
            app.session_state['dependency_detection_done'] = True
            app.run()
            self.assertEqual(app.radio(key='design_section').value, 'system')
            self.assertFalse(any(button.label == 'Lancer la conception' for button in app.button))
            app.text_area(key='aptamer_input').set_value('CAGCTTATTCAATTATACCAGCAC')
            app.text_area(key='trigger_input').set_value('GTGCTGGTATAA')
            app.run()
            app.checkbox(key='extension_enabled').check().run()
            self.assertEqual(app.radio(key='design_section').value, 'extension')
            app.radio(key='inverse_extension_enabled').set_value(True).run()
            self.assertFalse(any(button.label == 'Étendre le trigger' for button in app.button))
            app.radio(key='design_section').set_value('switch').run()
            self.assertFalse(app.exception)
            start = next(button for button in app.button if button.label == 'Lancer la conception')
            self.assertFalse(start.disabled)
            with patch('aptaswitch_studio.web_runtime.DesignRunJob.start') as begin:
                start.click().run()
            self.assertFalse(app.exception)
            begin.assert_called_once()
            req = app.session_state['design_job'].request
            self.assertEqual(req.sequences.trigger_dna, 'GTGCTGGTATAA')
            self.assertEqual(req.sequences.effective_trigger(), 'GTGCTGGTATAA')
            self.assertEqual(req.design_trigger_pattern(), 'X'*18+'GTGCTGGTATAA')
            self.assertEqual(req.architecture.trigger_region_length, 30)
            self.assertEqual(req.inverse_extension.thermo.temperature_c, 37)
            self.assertEqual(req.thermo.temperature_c, 29)
            self.assertEqual(app.session_state['trigger_input'], 'GTGCTGGTATAA')
