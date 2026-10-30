"""The public system route uses an identified native backend, including archived inputs."""
from dataclasses import asdict, replace
import json
import unittest
from unittest.mock import patch
from aptaswitch_core import system_analysis
from aptaswitch_core.models import ThermoConditions
from tests.guide_fixtures import GuideBackend

class SystemAnalysisTests(unittest.TestCase):
    def test_native_runtime_selection_normalization_and_conditions(self):
        backend = GuideBackend([('((..+..))', -3.25)])
        conditions = ThermoConditions(temperature_c=26.75, sodium_m=.242, switch_material="dna")
        with patch('aptaswitch_core.thermodynamics.setup.get_backend', return_value=backend) as load:
            result = system_analysis.analyze_system('/old/native', 'ac gu\n', 'tg ca', conditions,
                                                   thermodynamics_library_path='/native/library')
        load.assert_called_once_with('/native/library')
        spec, c = backend.calls[0]
        self.assertEqual([(s.sequence, s.material) for s in spec.strands], [('ACGT','dna'), ('TGCA','dna')])
        self.assertEqual((c.celsius, c.sodium, c.magnesium), (26.75, .242, 0.))
        self.assertEqual(result['paired_trigger_indices'], [2, 3])
        self.assertEqual(result['thermodynamics'], asdict(backend.provenance))

        self.assertEqual(json.loads(json.dumps(result, allow_nan=False)), result)

    def test_other_system_chemistries_are_rejected_without_folding(self):
        backend = GuideBackend([('((..+..))', -3.25)])
        for aptamer, trigger in (("rna", "rna"), ("rna", "dna"), ("dna", "rna")):
            with self.subTest(aptamer=aptamer, trigger=trigger), self.assertRaises(ValueError):
                system_analysis.analyze_system(None, 'ACGT', 'TGCA',
                    ThermoConditions(switch_material=aptamer, trigger_material=trigger), backend=backend)
        self.assertEqual(backend.calls, [])

    def test_self_pairs_do_not_count_as_aptamer_contacts(self):
        backend = GuideBackend([('((....((+))....))', -3.25)])
        result = system_analysis.analyze_system(None, 'ACGTACGT', 'ACGTACGT', ThermoConditions(switch_material="dna"), backend=backend)
        self.assertEqual(result['paired_trigger_indices'], [0, 1, 6, 7])

    def test_invalid_or_disconnected_minimum_is_never_rendered_as_a_complex(self):
        for structure, energy in [('....+...', -1), ('........', -1), ('((..+....', -1),
                                  ('....+....', -1), ('((..+..))', float('nan'))]:
            with self.subTest(structure=structure), self.assertRaises((ValueError, RuntimeError)):
                system_analysis.analyze_system(None, 'ACGT', 'TGCA', ThermoConditions(switch_material="dna"),
                                                backend=GuideBackend([(structure, energy)]))

    def test_invalid_sequence_is_rejected_before_folding(self):
        for aptamer, trigger in [('', 'ACGT'), ('ACGN', 'TGCA'), ('ACGT', '')]:
            backend = GuideBackend([])
            with self.assertRaises(ValueError):
                system_analysis.analyze_system(None, aptamer, trigger, ThermoConditions(switch_material="dna"), backend=backend)
            self.assertEqual(backend.calls, [])
