"""The four revised controls and saved diagnostics describe the same regions."""
from dataclasses import replace
import copy
import tempfile
from pathlib import Path
import unittest

import pandas as pd

from aptaswitch_core.architecture import target_switch_structures
from aptaswitch_core.defect_metrics import functional_defect_metrics
from aptaswitch_core.models import FunctionalDefectWeights
from aptaswitch_core.structure_ensemble import ensemble_inputs, ensemble_key, load_ensemble, save_ensemble
from aptaswitch_studio.design_session import open_design_into_session, validate_design_state
from aptaswitch_studio.design_library import create_design
from tests.test_defect_ui import make_result, result_app
from tests.defect_fixtures import target_matrix
from tests.test_structure_ensemble import fake_analysis
from tests.ui_helpers import FrenchAppTest as AppTest


def functional_result(weights=None):
    result = make_result()
    result.request = replace(result.request, defect_metric='functional',
                             functional_defect_weights=weights or FunctionalDefectWeights())
    off, on = target_switch_structures(result.request.architecture, len(result.request.reporter.sequence_after_start_rna))
    result.candidates[0].regional_defects = functional_defect_metrics(
        result.request.architecture, len(result.request.reporter.sequence_after_start_rna),
        target_matrix(off), target_matrix(on), result.request.functional_defect_weights)
    return result


class FunctionalMetricUITests(unittest.TestCase):
    def test_four_new_regions_and_two_off_stem_details(self):
        result = functional_result()
        original = copy.deepcopy(result)
        app = AppTest.from_function(result_app)
        app.session_state.saved_result = result
        app.session_state.weight_functional_on_rbs_linker = 20.
        app.run()
        self.assertFalse(app.exception)
        zones = app.table[1].value
        self.assertEqual(list(zones['Zone']), ['Toehold libre', 'Tiges avant/après AUG',
                                              'Préfixe RBS → AUG libre', 'Préfixe RBS → fin linker libre'])
        self.assertEqual(list(zones['Part du compromis (%)']), ['26.67 %']*3 + ['20 %'])
        self.assertEqual(list(zones['Défaut d’ensemble (0–1)']), ['0']*4)
        self.assertEqual(list(app.table[2].value['Paires cibles']), ['6', '9'])
        captions = ' '.join(c.value for c in app.caption)
        help_text = ' '.join(c.proto.help for c in app.markdown)
        self.assertIn('séquestration indirecte', help_text)
        self.assertIn('se recouvrent', help_text)
        self.assertNotIn('ancienne définition', captions.lower())
        self.assertEqual(app.session_state.saved_result, original)

    def test_disabled_on_has_two_measured_local_defects_but_no_aggregate(self):
        result = functional_result(FunctionalDefectWeights(1, 1, 0, 0))
        app = AppTest.from_function(result_app)
        app.session_state.saved_result = result
        app.run()
        self.assertFalse(app.exception)
        self.assertTrue(pd.isna(app.table[0].value['Défaut ON'].iloc[0]))
        self.assertEqual(list(app.table[1].value['Part du compromis (%)']), ['50 %', '50 %', '0 %', '0 %'])
        self.assertEqual(list(app.table[1].value['Défaut d’ensemble (0–1)']), ['0']*4)

    def test_new_weights_restore_and_validate_as_project_fields(self):
        initial = {f'weight_functional_{key}': value for key, value in
                   {'off_toehold': 2, 'off_stem': 3, 'on_rbs_aug': 4, 'on_rbs_linker': .5}.items()}
        with tempfile.TemporaryDirectory() as folder:
            design = create_design(Path(folder), 'Current design', state={'form': initial})
            form = {}
            open_design_into_session(design.path, form)
        for key, value in {'off_toehold': 2, 'off_stem': 3, 'on_rbs_aug': 4, 'on_rbs_linker': .5}.items():
            self.assertEqual(form[f'weight_functional_{key}'], value)
        validate_design_state({'form': form})
        form['weight_functional_on_rbs_aug'] = -1
        with self.assertRaises(ValueError):
            validate_design_state({'form': form})

    def test_family_cache_uses_functional_regions_and_ignores_only_score_weights(self):
        result = functional_result()
        analysis = fake_analysis(result)
        self.assertEqual(analysis['inputs']['functional_scheme'], 'functional_v2')
        self.assertEqual(set(analysis['states']['on']['regions']), {'on_rbs_aug', 'on_rbs_linker'})
        with tempfile.TemporaryDirectory() as directory:
            save_ensemble(Path(directory), analysis)
            self.assertEqual(load_ensemble(Path(directory), analysis['inputs']), analysis)
        other = replace(result.request, functional_defect_weights=FunctionalDefectWeights(1, 2, 3, 4))
        self.assertEqual(analysis['key'], ensemble_key(ensemble_inputs(other, result.candidates[0])))
