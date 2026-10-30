"""User journey: retain, reopen, configure, and inspect a full-context test."""
from copy import deepcopy
from dataclasses import asdict, replace
from pathlib import Path
from types import SimpleNamespace
import tempfile
import unittest
from unittest.mock import Mock, patch

from streamlit.testing.v1 import AppTest

from aptaswitch_core.models import RunRequest, SequenceInput, SwitchArchitecture, REPORTERS
from aptaswitch_core.saved_switch_analysis import build_saved_switch_context, validate_analysis_config
from aptaswitch_studio import saved_switch_library as library
from aptaswitch_studio import saved_switch_runtime as runtime
from aptaswitch_studio.web_runtime import JobSnapshot
from design_fixtures import design_result_fixture

SCRIPT = """
import streamlit as st
from aptaswitch_core.localization import set_language
from aptaswitch_studio.saved_switch_ui import render_save_candidate, render_saved_switches_page
set_language(st.session_state.get('language', 'fr'))
if st.session_state.get('nav') == 'Toehold switches enregistrés':
    render_saved_switches_page()
else:
    result = st.session_state['result']
    render_save_candidate(result, result.candidates[0])
"""


class SavedSwitchUITests(unittest.TestCase):
    def setUp(self):
        from tests.native_ui_fixtures import patch_native_runtime
        self.native_backend = patch_native_runtime(self)
        self.folder = Path(self.enterContext(tempfile.TemporaryDirectory())).resolve()
        self.design = self.folder / "Design"
        self.design.mkdir()
        request = RunRequest(sequences=SequenceInput("AACCGGTTAACC", "AACC"),
                             architecture=SwitchArchitecture(toehold_length=2, stem_length=2),
                             trials=1, output_dir=self.design / "calculations")
        self.result = design_result_fixture(request)
        self.enterContext(patch("aptaswitch_studio.saved_switch_ui.configured_library_path", return_value="/native"))
        self.addCleanup(runtime._JOBS.clear)

    def app(self, *, page="results", state=None):
        app = AppTest.from_string(SCRIPT, default_timeout=10)
        app.session_state["result"] = self.result
        app.session_state["active_design_path"] = str(self.design)
        app.session_state["nav"] = page
        for key, value in (state or {}).items():
            app.session_state[key] = value
        app.run()
        self.assertFalse(app.exception)
        return app

    def button(self, app, label):
        return next(button for button in app.button if button.label == label)

    def section(self, app, section):
        app.radio(key="saved_switch_section_selector").set_value(section).run()
        self.assertFalse(app.exception)
        return app

    def save(self, name="Premier choix"):
        snap = library.candidate_snapshot(self.result, self.result.candidates[0], name)
        return library.save_switch(library.saved_switch_root(self.design), name, snap)

    def assert_context_sequences_hidden(self, app):
        self.assertFalse(app.text_area)
        self.assertFalse(app.code)
        self.assertFalse(any(item.label in {"Séquence complète du rapporteur enregistrée", "Séquence du transcrit complet"} for item in app.expander))

    def test_saved_extension_comparison_is_displayed_from_detached_snapshot(self):
        candidate = self.result.candidates[0]
        aptamer = self.result.request.sequences.aptamer_dna
        trigger = candidate.trigger_seq
        reference = {"aptamer_sequence": aptamer, "trigger_sequence": trigger,
                     "mfe_structure": "." * len(aptamer) + "+" + "." * len(trigger),
                     "mfe_energy_kcal_mol": 0.}
        candidate.inverse_extension_metrics = {
            "status": "completed", "original_switch_rank": 14, "final_score": .25,
            "extended_trigger": trigger, "extension_length": 0,
            "mfe_structure": reference["mfe_structure"], "mfe_energy_kcal_mol": 0.,
        }
        self.result.inverse_extension_analysis = {
            "status": "completed", "reference": reference,
            "selected_count": 20, "analyzed_count": 20,
        }
        saved = self.save()
        original = saved.manifest_path.read_bytes()
        with patch("aptaswitch_core.inverse_extension.analyze_inverse_extensions",
                   side_effect=AssertionError("Saved comparison must not be recomputed")):
            app = self.app(page="Toehold switches enregistrés")
        self.assertIn("Comparaison de l’extension", [tab.label for tab in app.tabs])
        metrics = {metric.label: metric.value for metric in app.metric}
        self.assertEqual(metrics["Rang switch initial"], "14")
        self.assertEqual(metrics["Impact structurel"], "0.25")
        self.assertFalse(any(button.label == "Lancer une nouvelle comparaison" for button in app.button))
        self.assertEqual(saved.manifest_path.read_bytes(), original)

    def test_save_open_and_reload_preserves_exact_candidate_and_full_reporter(self):
        app = self.app()
        app.text_input[0].set_value("Mon candidat").run()
        self.button(app, "Enregistrer ce toehold switch").click().run()
        saved = library.list_saved_switches(library.saved_switch_root(self.design))[0]
        self.assertEqual(saved.name, "Mon candidat")
        self.assertEqual(saved.snapshot["candidate"], self.result.candidates[0].to_dict())
        self.assertEqual(saved.snapshot["request"]["reporter"]["full_sequence_after_start_rna"], REPORTERS["sfgfp"].full_sequence_after_start_rna)
        self.button(app, "Ouvrir les toehold switches enregistrés").click().run()
        self.assertFalse(app.exception)
        self.assertEqual(app.selectbox(key="saved_switch_selector").options, ["Mon candidat"])
        self.assertFalse(any(button.label == "Lancer l’analyse thermodynamique" for button in app.button))
        self.section(app, "Thermodynamique")
        self.assert_context_sequences_hidden(app)
        self.button(app, "Enregistrer le contexte").click().run()
        self.assertFalse(app.exception)
        context = runtime.load_context(saved)
        self.assertEqual(context["context"]["switch_sequence"], self.result.candidates[0].switch_seq[:-30] + REPORTERS["sfgfp"].full_sequence_after_start_rna)
        restarted = self.app(page="Toehold switches enregistrés")
        self.section(restarted, "Thermodynamique")
        self.assert_context_sequences_hidden(restarted)

    def test_full_application_routes_results_save_button_to_the_new_tab(self):
        from ui_helpers import prepare_test_design

        for name, value in (("load_settings", {}),):
            self.enterContext(patch("aptaswitch_core.settings." + name, return_value=value))
        app = AppTest.from_file(str(Path(__file__).parents[1] / "src/aptaswitch_studio/web_app.py"), default_timeout=15)
        prepare_test_design(app)
        app.session_state["ui_language"] = "fr"
        app.session_state["nav"] = "Résultats"
        app.session_state["last_result"] = self.result
        app.run()
        self.assertFalse(app.exception)
        self.button(app, "Enregistrer ce toehold switch").click().run()
        self.assertFalse(app.exception)
        self.button(app, "Ouvrir les toehold switches enregistrés").click().run()
        self.assertFalse(app.exception)
        self.assertEqual(app.radio(key="nav").value, "Toehold switches enregistrés")
        self.assertFalse("design_save_error" in app.session_state and app.session_state["design_save_error"])
        self.section(app, "Thermodynamique")
        self.assert_context_sequences_hidden(app)

    def test_selector_preserves_candidate_context_and_navigation(self):
        first = self.save("Premier choix")
        second = self.save("Second choix")
        runtime.save_context(first, {"temperature_c": 34.})
        runtime.save_context(second, {"temperature_c": 42.})
        script = SCRIPT.replace(
            "if st.session_state.get('nav')",
            "st.sidebar.radio('Page', ['results', 'Toehold switches enregistrés'], key='nav')\nif st.session_state.get('nav')",
        )
        app = AppTest.from_string(script, default_timeout=10)
        app.session_state["result"] = self.result
        app.session_state["active_design_path"] = str(self.design)
        app.session_state["nav"] = "Toehold switches enregistrés"
        app.session_state["saved_switch_selected"] = str(first.path)
        app.session_state["saved_switch_section"] = "Thermodynamique"
        with patch("aptaswitch_studio.saved_switch_ui.active_ffs", side_effect=lambda saved: saved.path == second.path):
            app.run()
            self.assertFalse(app.exception)
            self.assertIn("Second choix · Calcul en cours", app.selectbox(key="saved_switch_selector").options)
            for saved, temperature in ((second, 42.), (first, 34.), (second, 42.)):
                app.selectbox(key="saved_switch_selector").set_value(str(saved.path)).run()
                self.assertFalse(app.exception)
                self.assertEqual(app.session_state["saved_switch_selected"], str(saved.path))
                self.assertEqual(next(item for item in app.number_input if item.label == "Température (°C)").value, temperature)
            app.radio(key="nav").set_value("results").run()
            self.assertFalse(app.exception)
            self.assertEqual(app.session_state["saved_switch_selected"], str(second.path))
            app.radio(key="nav").set_value("Toehold switches enregistrés").run()
            self.assertFalse(app.exception)
            self.assertEqual(app.selectbox(key="saved_switch_selector").value, str(second.path))
            self.assertEqual(app.radio(key="saved_switch_section_selector").value, "Thermodynamique")


    def test_duplicate_requires_explicit_overwrite_and_retains_backup(self):
        app = self.app()
        self.button(app, "Enregistrer ce toehold switch").click().run()
        saved = library.list_saved_switches(library.saved_switch_root(self.design))[0]
        old_id = saved.snapshot["snapshot_id"]
        self.button(app, "Enregistrer ce toehold switch").click().run()
        self.assertEqual(library.load_switch(saved.path).snapshot["snapshot_id"], old_id)
        self.button(app, "Confirmer l’écrasement du toehold switch").click().run()
        self.assertFalse(app.exception)
        self.assertNotEqual(library.load_switch(saved.path).snapshot["snapshot_id"], old_id)
        self.assertTrue((saved.path.parent / ".backups").is_dir())

    def test_mixed_magnesium_disabled_optional_aptamer_and_sweep_are_persisted(self):
        saved = self.save()
        app = self.app(page="Toehold switches enregistrés")
        self.section(app, "Thermodynamique")
        magnesium = next(item for item in app.number_input if item.label == "Mg²⁺")
        self.assertTrue(magnesium.disabled)
        self.assertEqual(magnesium.value, 0.)
        self.assertEqual(magnesium.proto.help, "Mg²⁺ : 0 mM pour ce modèle.")
        self.assertFalse(any(item.label == "Balayer Mg²⁺" for item in app.checkbox))
        self.assertFalse(any(item.label == "Aptamère" for item in app.number_input))
        aptamer = next(item for item in app.checkbox if item.label == "Ajouter l’aptamère aux deux tubes")
        self.assertFalse(aptamer.value)
        aptamer.check().run()
        self.assertTrue(any(item.label == "Aptamère" for item in app.number_input))
        next(item for item in app.checkbox if item.label == "Balayer la température (°C)").check().run()
        self.button(app, "Enregistrer le contexte").click().run()
        self.assertFalse(app.exception)
        config = runtime.load_context(saved)["config"]
        self.assertTrue(config["include_aptamer"])
        self.assertEqual(len(config["temperature_values_c"]), 5)
        self.assertEqual(config["magnesium_values_m"], [])


    def test_rna_candidate_keeps_full_reporter_in_context_without_sequence_blocks(self):
        saved = self.save()
        app = self.app(page='Toehold switches enregistrés')
        self.section(app, "Thermodynamique")
        self.assert_context_sequences_hidden(app)
        self.assertFalse(self.button(app, 'Lancer l’analyse thermodynamique').disabled)
        self.button(app, 'Enregistrer le contexte').click().run()
        self.assertFalse(app.exception)
        context = runtime.load_context(saved)['context']
        expected = self.result.candidates[0].switch_seq[:-30] + REPORTERS['sfgfp'].full_sequence_after_start_rna
        self.assertEqual(context['switch_sequence'], expected)
        self.assertNotIn('T', context['switch_sequence'])
        self.assertNotIn('U', context['trigger_sequence'])

    def test_context_units_and_salt_sweep_survive_navigation_save_and_restart(self):
        saved = self.save()
        key = "saved_switch_" + saved.snapshot["snapshot_id"]
        with patch("aptaswitch_studio.saved_switch_ui.start_analysis") as start:
            app = self.app(page="Toehold switches enregistrés")
            self.section(app, "Thermodynamique")
            app.selectbox(key=key + "_sodium_unit").set_value("M").run()
            app.number_input(key=key + "_sodium").set_value(.175).run()
            for field, value in (("switch_concentration_m", 250.), ("trigger_concentration_m", 125.)):
                app.selectbox(key=key + field + "_unit").set_value("nM").run()
                app.number_input(key=key + field).set_value(value).run()
            app.checkbox(key=key + "_aptamer").check().run()
            app.selectbox(key=key + "_aptamer_concentration_unit").set_value("nM").run()
            app.number_input(key=key + "_aptamer_concentration").set_value(75.).run()
            app.checkbox(key=key + "_na_scan_enabled").check().run()
            app.selectbox(key=key + "_na_scan_unit").set_value("M").run()
            app.number_input(key=key + "_na_scan_min").set_value(.15).run()
            app.number_input(key=key + "_na_scan_max").set_value(.25).run()
            app.number_input(key=key + "_na_scan_count").set_value(3).run()
            self.assertFalse(app.exception)
            self.section(app, "Candidat")
            self.section(app, "Thermodynamique")
            self.assertEqual(app.selectbox(key=key + "_sodium_unit").value, "M")
            self.assertEqual(app.number_input(key=key + "_sodium").value, .175)
            self.assertEqual(app.selectbox(key=key + "switch_concentration_m_unit").value, "nM")
            self.assertEqual(app.number_input(key=key + "switch_concentration_m").value, 250.)
            self.assertEqual(app.selectbox(key=key + "_na_scan_unit").value, "M")
            self.assertEqual(app.number_input(key=key + "_na_scan_min").value, .15)
            self.assertEqual(app.number_input(key=key + "_na_scan_max").value, .25)
            self.button(app, "Enregistrer le contexte").click().run()
            self.assertFalse(app.exception)
            config = runtime.load_context(saved)["config"]
            for field, value in (("sodium_m", .175), ("switch_concentration_m", 250e-9),
                                 ("trigger_concentration_m", 125e-9), ("aptamer_concentration_m", 75e-9)):
                self.assertEqual(config[field], value)
            self.assertEqual(config["sodium_values_m"], [.15, .2, .25])
            self.assertEqual(config["magnesium_m"], 0.)
            restarted = self.app(page="Toehold switches enregistrés", state={"saved_switch_section": "Thermodynamique"})
            restored = restarted.session_state["saved_switch_context_drafts"][saved.snapshot["snapshot_id"]]
            self.assertEqual(restored, config)
            start.assert_not_called()

    def test_concentration_unit_change_does_not_change_the_saved_calculation_key(self):
        saved = self.save()
        key = "saved_switch_" + saved.snapshot["snapshot_id"]
        app = self.app(page="Toehold switches enregistrés", state={"saved_switch_section": "Thermodynamique"})
        initial = dict(app.session_state["saved_switch_context_drafts"][saved.snapshot["snapshot_id"]])
        signature = runtime.configuration_key(saved, initial)
        for unit in ("M", "nM", "µM", "mM"):
            for field in ("_sodium", "switch_concentration_m", "trigger_concentration_m"):
                app.selectbox(key=key + field + "_unit").set_value(unit).run()
                self.assertFalse(app.exception)
            draft = app.session_state["saved_switch_context_drafts"][saved.snapshot["snapshot_id"]]
            self.assertEqual(draft, initial)
            self.assertEqual(runtime.configuration_key(saved, draft), signature)

    def test_active_ffs_prevents_overwrite(self):
        saved = self.save("toehold switch 1")
        app = self.app()
        job = SimpleNamespace(cancel=Mock())
        with patch("aptaswitch_studio.saved_switch_ui.active_ffs", return_value=job):
            self.button(app, "Enregistrer ce toehold switch").click().run()
        self.assertTrue(app.error)
        self.assertEqual(library.load_switch(saved.path).snapshot["snapshot_id"], saved.snapshot["snapshot_id"])

    def test_running_ffs_allows_full_context_analysis_without_changing_the_simulation(self):
        saved = self.save()
        folder = self.design / 'running-ffs'
        folder.mkdir()
        request_file = folder / 'worker_request.json'
        request_file.write_text('{"replicas":10,"guide":"frozen"}')
        before = request_file.read_bytes()
        frozen = {'kinetics': {'replicas': 10}, 'prepared': {'guide': 'frozen'}}
        ffs = SimpleNamespace(folder=folder, request=deepcopy(frozen), cancel=Mock(),
                              snapshot=Mock(return_value=JobSnapshot('running', 2, 10, 'FFS running')))
        analysis = SimpleNamespace(run_id='parallel_context',
                                   snapshot=lambda: JobSnapshot('done', 1, 1, ''))
        with patch('aptaswitch_studio.saved_switch_ui.active_ffs', return_value=ffs), \
             patch('aptaswitch_studio.saved_switch_ui.start_analysis', return_value=(analysis, True)) as start:
            app = self.app(page='Toehold switches enregistrés', state={
                'saved_switch_section': 'Thermodynamique', 'saved_switch_ffs_job': ffs})
            self.assertFalse(self.button(app, 'Lancer l’analyse thermodynamique').disabled)
            self.button(app, 'Lancer l’analyse thermodynamique').click().run()
            self.assertFalse(app.exception)
            start.assert_called_once()
            self.assertEqual(start.call_args.args[1]['max_complex_size'], 3)
            self.assertIs(app.session_state['saved_switch_ffs_job'], ffs)
        ffs.cancel.assert_not_called()
        self.assertEqual(ffs.request, frozen)
        self.assertEqual(request_file.read_bytes(), before)
        self.assertEqual(library.load_switch(saved.path).snapshot['snapshot_id'], saved.snapshot['snapshot_id'])

    def test_other_running_jobs_still_disable_context_analysis_even_with_an_active_ffs(self):
        self.save()
        for key in ('saved_switch_free_md_job', 'extension_job', 'design_job'):
            with self.subTest(job=key):
                ffs = SimpleNamespace(cancel=Mock(), snapshot=lambda: JobSnapshot('running', 2, 10, 'FFS running'))
                other = SimpleNamespace(cancel=Mock(), snapshot=lambda: JobSnapshot('running', 0, 1, 'Other running'))
                with patch('aptaswitch_studio.saved_switch_ui.active_ffs', return_value=ffs), \
                     patch('aptaswitch_studio.saved_switch_ui.start_analysis') as start:
                    app = self.app(page='Toehold switches enregistrés', state={
                        'saved_switch_section': 'Thermodynamique', 'saved_switch_ffs_job': ffs, key: other})
                    self.assertTrue(self.button(app, 'Lancer l’analyse thermodynamique').disabled)
                    start.assert_not_called()
                ffs.cancel.assert_not_called()
                other.cancel.assert_not_called()

    def test_running_ffs_allows_calculating_missing_full_context_structures(self):
        saved = self.save()
        config = validate_analysis_config(saved.snapshot, {})
        context = build_saved_switch_context(saved.snapshot, config)
        sequence = context['switch_sequence']
        trigger = context['trigger_sequence']
        point = {'conditions': {key: config[key] for key in ('temperature_c', 'sodium_m', 'magnesium_m')},
                 'metrics': {'on_yield_pct': 65., 'on_trigger_bound_pct': 60., 'off_switch_monomer_pct': 75.}}
        for state in ('off', 'on'):
            names = ['Switch'] if state == 'off' else ['Trigger', 'Switch']
            sequences = [sequence] if state == 'off' else [trigger, sequence]
            mfe = {'strand_order': names, 'sequences': sequences, 'energy_kcal_mol': 0.,
                   'structure': '+'.join('.' * len(value) for value in sequences)}
            rows = [{'label': ' + '.join(names), 'strand_order': names, 'complex_key': '|'.join(names),
                     'switch_count': 1, 'trigger_count': int(state == 'on'), 'concentration_m': 1e-7, 'mfe': mfe}]
            if state == 'off':
                rows.append({'label': '2 Switch', 'strand_order': ['Switch', 'Switch'],
                             'complex_key': 'Switch|Switch', 'switch_count': 2, 'trigger_count': 0,
                             'concentration_m': 1e-8})
            point[state] = {'mfe': mfe, 'population': {'complexes': rows}}
        recorded = {'status': 'completed', 'run_id': 'missing_structure', 'config': config,
                    'context': context, 'baseline': point, 'sweeps': {}}
        before = deepcopy(recorded)
        ffs = SimpleNamespace(cancel=Mock(), snapshot=lambda: JobSnapshot('running', 2, 10, 'FFS running'))
        analysis = SimpleNamespace(run_id='new_structures', snapshot=lambda: JobSnapshot('done', 1, 1, ''))
        script = """
import streamlit as st
from aptaswitch_core.localization import set_language
from aptaswitch_studio.saved_switch_ui import _render_analysis
set_language('fr')
_render_analysis(st.session_state.recorded, st.session_state.saved)
"""
        with patch('aptaswitch_studio.saved_switch_ui.start_analysis', return_value=(analysis, True)) as start:
            app = AppTest.from_string(script, default_timeout=10)
            app.session_state['recorded'] = recorded
            app.session_state['saved'] = saved
            app.session_state['saved_switch_ffs_job'] = ffs
            app.run()
            self.assertFalse(app.exception)
            next(widget for widget in app.selectbox if widget.key.endswith('_off_complex')).set_value('Switch|Switch').run()
            self.assertFalse(app.exception)
            self.assertFalse(self.button(app, 'Calculer les structures de ces complexes').disabled)
            self.button(app, 'Calculer les structures de ces complexes').click().run()
            self.assertFalse(app.exception)
            start.assert_called_once()
            self.assertEqual(start.call_args.args[1]['max_complex_size'], 3)
            self.assertIs(app.session_state['saved_switch_ffs_job'], ffs)
        ffs.cancel.assert_not_called()
        self.assertEqual(recorded, before)

    def test_sections_render_lazily_without_starting_or_cancelling_jobs(self):
        self.save()
        live = SimpleNamespace(cancel=Mock())
        with patch("aptaswitch_studio.saved_switch_ui._render_candidate") as candidate, \
                patch("aptaswitch_studio.saved_switch_ui._render_thermodynamics") as external, \
                patch("aptaswitch_studio.saved_switch_ui.render_saved_switch_ffs") as kinetics, \
                patch("aptaswitch_studio.saved_switch_ui.active_ffs", return_value=live), \
                patch("aptaswitch_studio.saved_switch_ui.start_analysis") as start:
            renderers = (candidate, external, kinetics)
            app = self.app(page="Toehold switches enregistrés")
            candidate.assert_called_once()
            for renderer in renderers[1:]:
                renderer.assert_not_called()
            for selected, renderer in zip(("Thermodynamique", "Cinétique et film", "Candidat"),
                                          (external, kinetics, candidate)):
                for item in renderers:
                    item.reset_mock()
                self.section(app, selected)
                renderer.assert_called_once()
                for item in renderers:
                    if item is not renderer:
                        item.assert_not_called()
            live.cancel.assert_not_called()
            start.assert_not_called()

    def test_candidate_summary_only_appears_once_in_sequences_and_structures(self):
        self.save()
        with patch("aptaswitch_studio.saved_switch_ui.render_saved_switch_ffs"):
            app = self.app(page="Toehold switches enregistrés")
            for section in ("Candidat", "Thermodynamique", "Cinétique et film", "Candidat"):
                self.section(app, section)
                summaries = [item for item in app.metric if item.label == "Score net (/100)"]
                self.assertEqual(len(summaries), 1 if section == "Candidat" else 0)
                self.assertFalse(any(" · score " in item.value for item in (*app.caption, *app.markdown)))
                self.assertFalse(any(item.label == "Aptamère · structure 3D" for item in app.expander))

    def test_external_draft_survives_sections_but_kinetics_uses_saved_conditions(self):
        saved = self.save()
        runtime.save_context(saved, {"temperature_c": 34., "sodium_m": .12})
        with patch("aptaswitch_studio.saved_switch_ui.render_saved_switch_ffs") as kinetics, \
                patch("aptaswitch_studio.saved_switch_ui.start_analysis") as start:
            app = self.app(page="Toehold switches enregistrés")
            self.section(app, "Thermodynamique")
            next(item for item in app.number_input if item.label == "Température (°C)").set_value(42.).run()
            self.section(app, "Cinétique et film")
            self.assertEqual(kinetics.call_args.kwargs["thermo"].temperature_c, 34.)
            self.assertEqual(kinetics.call_args.kwargs["thermo"].sodium_m, .12)
            self.assertFalse(any(item.label == "Température (°C)" for item in app.number_input))
            self.section(app, "Thermodynamique")
            self.assertEqual(next(item for item in app.number_input if item.label == "Température (°C)").value, 42.)
            self.assertEqual(runtime.load_context(saved)["config"]["temperature_c"], 34.)
            self.button(app, "Enregistrer le contexte").click().run()
            self.section(app, "Cinétique et film")
            self.assertEqual(kinetics.call_args.kwargs["thermo"].temperature_c, 42.)
            start.assert_not_called()

    def test_section_labels_translate_without_changing_persisted_selection(self):
        self.save()
        app = self.app(page="Toehold switches enregistrés", state={"language": "en", "saved_switch_section": "Thermodynamique"})
        selector = app.radio(key="saved_switch_section_selector")
        self.assertEqual(selector.options, ["Sequences and structures", "Toehold switch context with the full reporter",
                                            "Kinetic simulation and movie"])
        self.assertEqual(selector.value, "Thermodynamique")

    def test_result_curves_structures_and_stored_conditions_render_after_restart(self):
        saved = self.save()
        config = validate_analysis_config(saved.snapshot, {"temperature_c": 29.})
        job = runtime.SavedSwitchAnalysisJob(saved, config, "/native")
        sequence = self.result.candidates[0].switch_seq[:-30] + REPORTERS["sfgfp"].full_sequence_after_start_rna
        point = {"conditions": {"temperature_c": 29., "sodium_m": .12, "magnesium_m": 0.}, "thermodynamics": asdict(self.native_backend.provenance),
                 "metrics": {"on_yield_pct": 65., "on_trigger_bound_pct": 60., "off_switch_monomer_pct": 75.}}
        for state in ("off", "on"):
            seqs = [sequence] if state == "off" else ["AACC", sequence]
            point[state] = {"population": {"complexes": [{"label": "Switch", "concentration_m": 5e-6}]},
                            "mfe": {"structure": "+".join("." * len(seq) for seq in seqs), "sequences": seqs,
                                    "strand_order": ["Switch"] if state == "off" else ["Trigger", "Switch"], "energy_kcal_mol": -10.}}
        with patch.object(runtime, "analyze_saved_switch", return_value={"status": "completed", "baseline": point, "sweeps": {"temperature_c": [point]}, "config": config, "model_limits": dict(self.native_backend.condition_limits), "context": build_saved_switch_context(saved.snapshot, config, backend=self.native_backend)}):
            job._run()
        app = self.app(page="Toehold switches enregistrés")
        self.section(app, "Thermodynamique")
        self.assertFalse(app.exception)
        self.assertEqual([item.value for item in app.metric], ["65.00 %", "60.00 %", "75.00 %"])
        self.assertEqual(len(app.get("vega_lite_chart")), 1)
        self.assertFalse(app.error)
        self.assertEqual([tab.label for tab in app.tabs], ["Structure 2D", "Séquence linéaire"] * 2)
        self.assertEqual(len([item for item in app.get("popover") if item.proto.popover.label == "Exporter"]), 4)
        self.assertEqual(len([item for item in app.get("popover") if item.proto.popover.label == "Structure"]), 4)
        # Editing physical inputs cannot recolor/refold the stored result.
        with patch("aptaswitch_core.thermodynamics.native_backend.NativeHybridBackend.analyze", side_effect=AssertionError("No native call while viewing")):
            next(item for item in app.number_input if item.label == "Température (°C)").set_value(42.).run()
        self.assertFalse(app.exception)
        self.assertFalse(app.error)
        self.assertTrue(any("29 °C" in item.value for item in app.caption))


if __name__ == "__main__":
    unittest.main()
