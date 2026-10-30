"""Saved candidates display their original design without refolding or the source run."""
from copy import deepcopy
from dataclasses import asdict, fields, replace
from pathlib import Path
from types import SimpleNamespace
import shutil
import tempfile
import unittest
from unittest.mock import patch

from aptaswitch_core.architecture import architecture_layout
from aptaswitch_core.models import DesignCandidate, FunctionalDefectWeights, InverseExtensionSettings, ScoringWeights
from aptaswitch_core.preview import DOMAIN_COLORS, TRIGGER_COLOR, render_structure_svg
from aptaswitch_core.scoring import selection_score_100
from aptaswitch_core.structure_ensemble import load_ensemble, save_ensemble
from aptaswitch_core.thermodynamics.types import ModelProvenance
from aptaswitch_studio import saved_switch_library as library
from aptaswitch_studio.saved_switch_runtime import save_context
from tests.test_structure_ensemble import fake_analysis, fixture
from tests.ui_helpers import FrenchAppTest as AppTest


CANDIDATE_SCRIPT = """
import streamlit as st
from aptaswitch_studio.saved_switch_library import load_switch
from aptaswitch_studio.saved_switch_ui import _render_candidate
_render_candidate(load_switch(st.session_state.candidate_path))
"""

SAVE_SCRIPT = """
import streamlit as st
from aptaswitch_studio.saved_switch_ui import render_save_candidate
result = st.session_state.result
render_save_candidate(result, result.candidates[0])
"""


class SavedDesignUITests(unittest.TestCase):
    def setUp(self):
        self.folder = Path(self.enterContext(tempfile.TemporaryDirectory())).resolve()
        self.design = self.folder / 'Design'
        self.result = fixture(self.design / 'calculations')
        self.candidate = self.result.candidates[0]
        # These synthetic MFE records intentionally differ from the design
        # targets and family representatives, so replacing either is detectable.
        length = len(self.candidate.switch_seq)
        self.candidate.structure_off = '.' * (length - 12) + '(((......)))'
        self.candidate.structure_on = '.' * len(self.candidate.trigger_seq) + '+' + '.' * length
        self.native = self.enterContext(patch(
            'aptaswitch_core.thermodynamics.native_backend.NativeHybridBackend.analyze',
            side_effect=AssertionError('Viewing saved design must not invoke a calculation')))
        self.job = self.enterContext(patch(
            'aptaswitch_studio.ensemble_ui.EnsembleJob',
            side_effect=AssertionError('Saved design families are read-only')))
        self.start_context = self.enterContext(patch(
            'aptaswitch_studio.saved_switch_ui.start_analysis',
            side_effect=AssertionError('Viewing design must not start a context calculation')))

    def save(self, *, with_source=False):
        snapshot = library.candidate_snapshot(self.result, self.candidate, 'Retained candidate')
        kwargs = {'source_output_dir': self.result.output_dir} if with_source else {}
        return library.save_switch(library.saved_switch_root(self.design), snapshot['name'], snapshot, **kwargs)

    def family_analysis(self, count, *, seed=0):
        backend = SimpleNamespace(provenance=ModelProvenance(**{
            field.name: self.result.thermodynamics[field.name] for field in fields(ModelProvenance)}))
        return fake_analysis(self.result, count, backend=backend, seed=seed)

    def app(self, saved):
        app = AppTest.from_string(CANDIDATE_SCRIPT, default_timeout=30)
        app.session_state.candidate_path = str(saved.path)
        app.run()
        self.assertFalse(app.exception)
        self.assertFalse(app.error)
        return app

    def assert_no_calculation(self, app):
        self.assertFalse(any('Analyser les sous-structures' in button.label
                             or 'Lancer le test' in button.label for button in app.button))
        self.native.assert_not_called()
        self.job.assert_not_called()
        self.start_context.assert_not_called()

    def provenance(self, app):
        return next(item for item in app.expander if item.label == 'Conditions du run et réglages du score')

    def table_rows(self, container):
        return [row for table in container.get('table') for row in table.value.to_dict(orient='records')]

    def row_for(self, container, label):
        rows = [row for row in self.table_rows(container) if label in row.values()]
        self.assertEqual(len(rows), 1, f'Expected a single {label!r} row: {rows!r}')
        return rows[0]

    def with_regional_metrics(self):
        weights = FunctionalDefectWeights(off_toehold=2., off_stem=1., on_rbs_aug=0., on_rbs_linker=1.)
        self.result.request = replace(self.result.request, functional_defect_weights=weights,
                                      scoring=ScoringWeights(on_yield=4.5, leak=.2, defect=8., ddg=.75,
                                                             rbs=3.5, trigger_linearity=2.25, rbs_calculator=0.))
        from tests.design_fixtures import synthetic_native_metadata
        self.result.optimization = synthetic_native_metadata(self.result.request)['optimization']
        self.candidate.regional_defects = {
            'scheme': 'functional_v2', 'weights': asdict(weights), 'normalized_weights': weights.fractions(),
            'off_toehold': .2, 'off_stem': .4, 'on_rbs_aug': .8, 'on_rbs_linker': .1,
            'off': .8 / 3., 'on': .1, 'balanced': .225,
            'regions': {key: {} for key in asdict(weights)},
        }
        self.candidate.defect_structural = .06
        self.candidate.defect_concentration = .04
        self.candidate.score_raw = 199.
        self.candidate.score = selection_score_100(199.)

    def assert_rendered_design(self, calls, *, structure, state):
        rendered = next(call for call in reversed(calls) if call.args[1] == structure)
        offset = len(self.candidate.trigger_seq) if state == 'on' else 0
        sequence, actual_structure, colors = rendered.args[:3]
        self.assertEqual(sequence, (self.candidate.trigger_seq if offset else '') + self.candidate.switch_seq)
        self.assertEqual(actual_structure, structure)
        self.assertEqual(len(colors), len(sequence))
        layout = architecture_layout(self.result.request.architecture)
        reporter_start = offset + layout['reporter_start']
        self.assertEqual(len(sequence) - reporter_start, 30)
        self.assertEqual(sequence[reporter_start:], self.result.request.reporter.sequence_after_start_rna)
        self.assertEqual(colors[reporter_start:], [DOMAIN_COLORS['Reporter']] * 30)
        self.assertEqual(colors[offset + layout['aug_start']:offset + layout['aug_start'] + 3],
                         [DOMAIN_COLORS['AUG']] * 3)
        self.assertEqual(colors[:offset], [TRIGGER_COLOR] * offset)

    def test_exact_off_on_mfe_and_reporter_prefix_survive_changed_full_context(self):
        saved = self.save()
        original = deepcopy(saved.snapshot)
        context = save_context(saved, {'temperature_c': 45., 'sodium_m': .3})
        self.assertGreater(context['context']['switch_length'], len(self.candidate.switch_seq))
        with patch('aptaswitch_core.preview.render_structure_svg', wraps=render_structure_svg) as render:
            app = self.app(saved)
            self.assertEqual([tab.label for tab in app.tabs], ['Structures MFE', 'Sous-familles du design'])
            provenance = self.provenance(app)
            self.assertEqual(self.row_for(provenance, 'Température (°C)')['Valeur'], '31.00')
            self.assertEqual(self.row_for(provenance, 'Na⁺ (mM)')['Valeur'], '150.000')
            self.assertEqual(self.row_for(provenance, 'Rapporteur utilisé pour le design')['Valeur'], 'sfGFP · 30 nt')
            self.assertFalse(any('30 nt de rapporteur' in item.value and '31 °C' in item.value
                                 for item in app.caption))
            self.assert_rendered_design(render.call_args_list, structure=self.candidate.structure_off, state='off')
            next(radio for radio in app.radio if radio.label == 'État MFE du design').set_value('on').run()
            self.assertFalse(app.exception)
            self.assertFalse(app.error)
            self.assert_rendered_design(render.call_args_list, structure=self.candidate.structure_on, state='on')
            self.assertIn(self.candidate.structure_on, [item.value for item in app.code])
        self.assertEqual(library.load_switch(saved.path).snapshot, original)
        self.assert_no_calculation(app)

    def test_native_provenance_and_selected_trial_remain_visible_offline(self):
        import json
        from aptaswitch_core.thermodynamics.optimization import TrialOptions
        from aptaswitch_core.thermodynamics.setup import MODEL_ID
        from tests.design_fixtures import synthetic_native_metadata
        options = asdict(TrialOptions(max_evaluations=37))
        self.result.request = replace(self.result.request,
                                      thermo=replace(self.result.request.thermo, material=MODEL_ID),
                                      engine='native-hybrid', trials=2, optimization_options=options)
        metadata = synthetic_native_metadata(self.result.request)
        provenance = metadata['thermodynamics']
        self.result.thermodynamics = provenance
        self.result.optimization = metadata['optimization']
        outcome = self.result.optimization['trials'][0]
        outcome.update(evaluations=37, best_loss=.2, stop_reason='budget_exhausted', budget_kind='evaluations')
        outcome['best_objective']['local_defects'] = [.2] * 4
        saved = self.save()
        before = saved.manifest_path.read_bytes()
        with patch('aptaswitch_studio.ensemble_ui.get_backend', side_effect=AssertionError('Runtime removed')):
            app = self.app(saved)
        settings = self.provenance(app)
        self.assertEqual(self.row_for(settings, 'Modèle thermodynamique')['Valeur'], MODEL_ID)
        self.assertEqual(self.row_for(settings, 'Moteur du résultat')['Valeur'], 'synthetic-rendering-fixture')
        self.assertEqual(self.row_for(settings, 'Défaut local du candidat (%)')['Valeur'], '0.000')
        self.assertEqual([json.loads(item.value) for item in settings.json],
                         [provenance, {**{key: self.result.optimization[key] for key in ('algorithm', 'objective', 'functional_weights')},
                                       'options': options, 'trial': outcome}])
        self.assertEqual(saved.manifest_path.read_bytes(), before)
        self.assert_no_calculation(app)

    def test_save_button_copies_both_family_caches_and_reopens_after_source_deleted(self):
        analyses = {count: self.family_analysis(count) for count in (1000, 3000)}
        for analysis in analyses.values():
            save_ensemble(self.result.output_dir, analysis)
        app = AppTest.from_string(SAVE_SCRIPT, default_timeout=30)
        app.session_state.result = self.result
        app.session_state.active_design_path = str(self.design)
        app.run()
        self.assertFalse(app.exception)
        next(button for button in app.button if button.label == 'Enregistrer ce toehold switch').click().run()
        self.assertFalse(app.exception)
        self.assertFalse(app.error)
        saved = library.list_saved_switches(library.saved_switch_root(self.design))[0]
        save_context(saved, {'temperature_c': 45., 'sodium_m': .3})
        shutil.rmtree(self.result.output_dir)
        self.assertFalse(self.result.output_dir.exists())
        from aptaswitch_studio.saved_switch_design import load_saved_design
        reopened_result, reopened_candidate = load_saved_design(library.load_switch(saved.path))
        self.assertTrue(reopened_result.output_dir.is_relative_to(saved.path))
        self.assertEqual(reopened_candidate.to_dict(), self.candidate.to_dict())
        for analysis in analyses.values():
            self.assertEqual(load_ensemble(reopened_result.output_dir, analysis['inputs']), analysis)
        before = {path.relative_to(saved.path): path.read_bytes() for path in saved.path.rglob('*') if path.is_file()}
        with patch('aptaswitch_core.preview.render_structure_svg', wraps=render_structure_svg) as render:
            app = self.app(saved)
            count = next(item for item in app.selectbox if item.label == 'Échantillonnage enregistré (structures par état)')
            self.assertEqual(count.options, ['1000', '3000'])
            count.set_value(1000).run()
            self.assertFalse(app.exception)
            for state in ('off', 'on'):
                next(item for item in app.radio if item.label == 'Ensemble à explorer').set_value(state).run()
                self.assertFalse(app.exception)
                self.assertFalse(app.error)
                self.assertEqual(app.get("table")[0].value['Tirages'].tolist(), [1000])
                self.assertEqual(app.get("table")[0].value['Fréquence estimée (%)'].tolist(), [100.])
                structure = analyses[1000]['states'][state]['structures'][0]['structure']
                self.assert_rendered_design(render.call_args_list, structure=structure, state=state)
            self.assertTrue(any('Analyse enregistrée' in item.value and '31' in item.value and '0.15' in item.value
                                for item in app.caption))
            next(item for item in app.selectbox if item.label == 'Échantillonnage enregistré (structures par état)').set_value(3000).run()
            self.assertFalse(app.exception)
            self.assertEqual(app.get("table")[0].value['Tirages'].tolist(), [3000])
            next(item for item in app.slider if item.label == 'Tolérance du regroupement (%)').set_value(0).run()
            next(item for item in app.checkbox if item.label == 'Afficher la matrice de probabilités').check().run()
            self.assertFalse(app.exception)
            self.assertTrue(app.get('vega_lite_chart'))
        self.assertEqual(before, {path.relative_to(saved.path): path.read_bytes()
                                  for path in saved.path.rglob('*') if path.is_file()})
        self.assert_no_calculation(app)

    def test_missing_families_explain_absence_without_calculation_controls(self):
        app = self.app(self.save())
        self.assertTrue(any(item.value == 'Aucune sous-famille enregistrée pour ce candidat dans le calcul de design.'
                            for item in app.info))
        self.assertFalse(any(item.label == 'Échantillonnage enregistré (structures par état)' for item in app.selectbox))
        self.assertFalse(any(item.label == 'Ensemble à explorer' for item in app.radio))
        self.assert_no_calculation(app)

    def test_saved_custom_sample_count_is_discovered_and_displayed_offline(self):
        from aptaswitch_core.thermodynamics.setup import MODEL_ID
        from tests.design_fixtures import synthetic_native_metadata
        self.result.request = replace(self.result.request, engine='native-hybrid',
                                      thermo=replace(self.result.request.thermo, material=MODEL_ID))
        metadata = synthetic_native_metadata(self.result.request)
        self.result.thermodynamics = metadata['thermodynamics']
        self.result.optimization = metadata['optimization']
        analysis = self.family_analysis(137, seed=23)
        save_ensemble(self.result.output_dir, analysis)
        saved = self.save(with_source=True)
        shutil.rmtree(self.result.output_dir)
        app = self.app(saved)
        count = next(item for item in app.selectbox if item.label == 'Échantillonnage enregistré (structures par état)')
        self.assertEqual(count.options, ['137'])
        self.assertEqual(app.get('table')[0].value['Tirages'].tolist(), [137])
        self.assert_no_calculation(app)

    def test_missing_mfe_does_not_substitute_design_target(self):
        self.candidate.structure_off = self.candidate.structure_on = ''
        with patch('aptaswitch_core.preview.render_structure_svg', wraps=render_structure_svg) as render:
            app = self.app(self.save())
            self.assertTrue(any(item.value == 'Aucune structure MFE n’a été conservée pour ce candidat.' for item in app.info))
            self.assertFalse(any(item.label == 'État MFE du design' for item in app.radio))
            render.assert_not_called()
        self.assert_no_calculation(app)

    def test_summary_follows_both_structure_tabs_and_does_not_duplicate_sequence_blocks(self):
        app = self.app(self.save())
        nodes = list(app.main)
        heading = next(item for item in app.subheader if item.value.startswith('Métriques du candidat'))
        position = next(index for index, node in enumerate(nodes) if node is heading)
        for tab in app.tabs:
            for child in tab:
                self.assertLess(next(index for index, node in enumerate(nodes) if node is child), position)
        self.assertEqual([item.label for item in app.metric],
                         ['Score net (/100)', 'Association ON (%)', 'Trigger–switch avec aptamère (%)', 'Défaut OFF/ON (0–1)'])
        self.assertFalse(any(item.label in {'Séquences et résultat d’origine',
                                           'Séquence complète du rapporteur enregistrée',
                                           'Séquence du transcrit complet', 'Aptamère · structure 3D'}
                             for item in app.expander))
        self.assertFalse(app.dataframe)
        self.assert_no_calculation(app)

    def test_absent_regional_metrics_remain_missing_instead_of_using_global_defect(self):
        app = self.app(self.save())
        metrics = {item.label: item.value for item in app.metric}
        self.assertEqual(metrics['Défaut OFF/ON (0–1)'], '—')
        self.assertEqual(float(metrics['Score net (/100)']), 50.)
        local_rows = [row for row in self.table_rows(app) if row.get('Critère', '').startswith(('OFF ·', 'ON ·'))]
        self.assertEqual(len(local_rows), 4)
        self.assertTrue(all(row['Valeur'] == '—' for row in local_rows))
        self.assertTrue(any('Défaut global : 0.100 · écart structural : — · manque de complexes : —'
                            in item.value for item in app.caption))
        self.assert_no_calculation(app)

    def test_original_metrics_and_weights_are_displayed_without_rescoring_or_context_overrides(self):
        self.with_regional_metrics()
        self.result.request = replace(self.result.request, thermo=replace(
            self.result.request.thermo, switch_concentration_m=7e-6, trigger_concentration_m=2e-6,
            screening_max_complex_size=2, validation_max_complex_size=3))
        saved = self.save()
        save_context(saved, {'temperature_c': 45., 'sodium_m': .3,
                             'switch_concentration_m': 50e-6, 'trigger_concentration_m': 60e-6})
        before = saved.manifest_path.read_bytes()
        with patch('aptaswitch_core.scoring.compute_selection_score', side_effect=AssertionError('Do not rescore')):
            app = self.app(saved)
        metrics = {item.label: item.value for item in app.metric}
        self.assertEqual(metrics, {'Score net (/100)': f'{self.candidate.score:.2f}',
                                  'Association ON (%)': '60.00', 'Trigger–switch avec aptamère (%)': '10.00',
                                  'Défaut OFF/ON (0–1)': '0.225'})
        measures = {row['Critère']: row['Valeur'] for row in self.table_rows(app) if 'Valeur' in row and 'Critère' in row}
        self.assertEqual({name: measures[name] for name in ('ΔG OFF', 'ΔG ON', 'ΔΔG', 'ΔG RBS–linker')},
                         {'ΔG OFF': '-15.000', 'ΔG ON': '-25.000', 'ΔΔG': '-8.000', 'ΔG RBS–linker': '-2.000'})
        self.assertEqual({name: measures[name] for name in measures if name.startswith(('OFF ·', 'ON ·'))},
                         {'OFF · Toehold libre': '0.200', 'OFF · Tiges avant/après AUG': '0.400',
                          'ON · Préfixe RBS → AUG libre': '0.800', 'ON · Préfixe RBS → fin linker libre': '0.100'})
        self.assertTrue(any('Défaut global : 0.100 · écart structural : 0.060 · manque de complexes : 0.040'
                            in item.value for item in app.caption))
        provenance = self.provenance(app)
        expected_conditions = {'Température (°C)': '31.00', 'Na⁺ (mM)': '150.000', 'Mg²⁺ (mM)': '0.000',
                               'Switch (µM)': '7.000', 'Trigger (µM)': '2.000',
                               'Complexes · criblage / validation': '2 / 3 brins maximum'}
        for name, expected in expected_conditions.items():
            self.assertEqual(self.row_for(provenance, name)['Valeur'], expected)
        for name, expected in {'ON % · bonus': '4.500', 'Trigger–switch avec aptamère · malus': '0.200',
                               'Défaut OFF/ON · malus': '8.000', 'ΔΔG · liaison': '0.750',
                               'ΔG RBS–linker': '3.500'}.items():
            self.assertEqual(self.row_for(provenance, name)['Poids enregistré'], expected)
        for name, weight, fraction in (('OFF · Toehold libre', '2.000', '50.0'),
                                      ('OFF · Tiges avant/après AUG', '1.000', '25.0'),
                                      ('ON · Préfixe RBS → AUG libre', '0.000', '0.0'),
                                      ('ON · Préfixe RBS → fin linker libre', '1.000', '25.0')):
            row = self.row_for(provenance, name)
            self.assertEqual((row['Poids enregistré'], row['Part du compromis (%)']), (weight, fraction))
        self.assertTrue(any('Score brut enregistré : 199.0000' in item.value for item in provenance.caption))
        self.assertNotIn('Trigger seul · défaut de linéarité', measures)
        self.assertEqual(saved.manifest_path.read_bytes(), before)
        self.assert_no_calculation(app)

    def test_unrecorded_score_weight_is_not_invented_from_current_defaults(self):
        snapshot = library.candidate_snapshot(self.result, self.candidate, 'Partial settings')
        del snapshot['request']['scoring']['leak']
        saved = library.save_switch(library.saved_switch_root(self.design), snapshot['name'], snapshot)
        app = self.app(saved)
        provenance = self.provenance(app)
        self.assertEqual(self.row_for(provenance, 'Trigger–switch avec aptamère · malus')['Poids enregistré'], '—')
        self.assertEqual(self.row_for(provenance, 'ON % · bonus')['Poids enregistré'], '3.000')

    def test_inverse_trigger_linearity_keeps_its_original_value_and_weight(self):
        self.with_regional_metrics()
        self.result.request = replace(self.result.request,
            sequences=replace(self.result.request.sequences, trigger_dna=self.candidate.trigger_seq[:24]),
            inverse_extension=InverseExtensionSettings(target_length=30))
        self.candidate.trigger_linearity = {
            'normalized_defect': .2, 'unpaired_fraction': .8,
            'per_base_unpaired': [.8] * len(self.candidate.trigger_seq),
            'structure_mfe': '.' * len(self.candidate.trigger_seq),
        }
        app = self.app(self.save())
        rows = [row for row in self.table_rows(app)
                if row.get('Critère') == 'Trigger seul · défaut de linéarité']
        self.assertEqual(len(rows), 2)
        self.assertEqual(next(row['Valeur'] for row in rows if 'Valeur' in row), '0.200')
        self.assertEqual(next(row['Poids enregistré'] for row in rows if 'Poids enregistré' in row), '2.250')
        self.assertTrue(any('− 100 × wlinéarité × défaut du trigger seul' in item.value
                            for item in self.provenance(app).caption))
        self.assert_no_calculation(app)

    def test_excluded_uncomputed_candidate_displays_dashes_without_invented_metrics(self):
        excluded = DesignCandidate.excluded_for_stop(
            trial=self.candidate.trial, switch_seq=self.candidate.switch_seq,
            trigger_seq=self.candidate.trigger_seq, first_stop=2,
            reporter_label=self.candidate.reporter_label,
        )
        self.result.candidates = [excluded]
        self.candidate = excluded
        saved = self.save()
        before = saved.manifest_path.read_bytes()
        app = self.app(saved)
        self.assertEqual([item.value for item in app.metric], ['—'] * 4)
        measures = [row for row in self.table_rows(app) if 'Critère' in row and 'Valeur' in row]
        self.assertEqual(len(measures), 8)
        self.assertTrue(all(row['Valeur'] == '—' for row in measures))
        self.assertEqual(saved.manifest_path.read_bytes(), before)
        self.assert_no_calculation(app)


if __name__ == '__main__':
    unittest.main()
