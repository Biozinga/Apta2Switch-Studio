"""Exercise real named-design navigation in a wholly temporary library."""

from __future__ import annotations

from contextlib import ExitStack
import io
import json
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import Mock, patch
import zipfile

from streamlit.testing.v1 import AppTest

from aptaswitch_studio.design_library import create_design, design_output_dir, load_design
from aptaswitch_studio.web_runtime import JobSnapshot
from tests.export_fixtures import export_extension_result


APP_PATH = Path(__file__).resolve().parents[1] / "src" / "aptaswitch_studio" / "web_app.py"
TRIGGER = "GTCCAGGCTGGTATAATTAGATCCACGTAC"


class DesignLibraryUITests(unittest.TestCase):
    def setUp(self):
        from tests.native_ui_fixtures import patch_native_runtime
        patch_native_runtime(self)
        self.context = ExitStack()
        self.addCleanup(self.context.close)
        self.workspace = Path(self.context.enter_context(tempfile.TemporaryDirectory())).resolve()
        self.library = self.workspace / "library"
        for target, value in (
            ("aptaswitch_core.settings.load_settings", {}),
            ("aptaswitch_core.thermodynamics.setup.configured_library_path", None),
            ("aptaswitch_studio.extension_ui.configured_library_path", None),
            ("aptaswitch_studio.design_library.default_design_root", self.library),
            ("aptaswitch_studio.design_library_ui.default_design_root", self.library),
        ):
            self.context.enter_context(patch(target, return_value=value))
        for target in ("save_settings",):
            self.context.enter_context(patch(
                f"aptaswitch_core.settings.{target}",
                side_effect=AssertionError("Design tests must not write user preferences"),
            ))

    def app(self, **state):
        app = AppTest.from_file(str(APP_PATH), default_timeout=30)
        app.session_state["dependency_detection_done"] = True
        app.session_state["ui_language"] = "fr"
        for key, value in state.items():
            app.session_state[key] = value
        app.run()
        self.assertFalse(app.exception)
        return app

    def create(self, app, name):
        app.text_input(key="design_new_name").set_value(name)
        app.button(key="create_design").click().run()
        self.assertFalse(app.exception)
        return load_design(Path(app.session_state["active_design_path"]))

    def fill_sequences(self, app, aptamer="ACGTACGT", trigger=TRIGGER):
        app.text_area(key="aptamer_input").set_value(aptamer)
        app.text_area(key="trigger_input").set_value(trigger).run()
        self.assertFalse(app.exception)

    def assert_no_geometry_controls(self, app):
        retired = {"customize_tswitch", "toehold_length", "rbs_aug_distance", "binding_start", "binding_end"}
        rendered = {item.key for kind in ("checkbox", "slider", "number_input") for item in app.get(kind)}
        self.assertTrue(retired.isdisjoint(rendered))
        self.assertTrue(retired.isdisjoint(app.session_state.filtered_state))

    def test_empty_session_lands_on_named_library_and_gates_scientific_pages(self):
        for page in (None, "Conception", "Résultats"):
            with self.subTest(page=page):
                app = self.app(**({"nav": page} if page else {}))
                self.assertEqual(app.radio(key="nav").value, "Mes designs")
                self.assertEqual(app.radio(key="nav").options, ["Mes designs", "Réglages"])
                self.assertTrue(any(button.key == "create_design" for button in app.button))
                self.assertFalse(any(area.key == "aptamer_input" for area in app.text_area))
                self.assertFalse(app.session_state.filtered_state.get("active_design_path"))
                app.radio(key="nav").set_value("Réglages").run()
                self.assertFalse(app.exception)
                self.assertEqual(app.radio(key="nav").value, "Réglages")

    def test_library_has_a_native_page_title_without_an_anchor(self):
        app = self.app()
        title = next(item for item in app.title if item.value == "Mes designs")
        self.assertTrue(title.proto.hide_anchor)
        self.assertFalse(any(item.value.startswith("## Mes designs") for item in app.markdown))

    def test_active_design_banner_shows_full_escaped_name_without_a_heading(self):
        project = create_design(self.library, "Toolbar")
        name = 'Un nom de design très long — ' * 6 + '<img src=x onerror="alert(1)"> & final'
        script = """
from aptaswitch_core.localization import set_language
from aptaswitch_studio.design_library_ui import render_design_toolbar
set_language('fr')
render_design_toolbar()
"""
        app = AppTest.from_string(script)
        app.session_state["active_design_path"] = str(project.path)
        app.session_state["active_design_name"] = name
        app.session_state["design_restore_warnings"] = ["Ancien réglage restauré"]
        app.run()
        self.assertFalse(app.exception)
        self.assertFalse(app.title)
        self.assertFalse(app.header)
        self.assertFalse(app.subheader)
        body = next(item.proto.body for item in app.get("html") if "active-design-name" in item.proto.body)
        from html import escape
        self.assertIn(escape(name), body)
        self.assertIn('class="active-design-label">Design actif</div>', body)
        self.assertNotIn("<img", body)
        self.assertIn("Sauvegarde automatique dans le dossier de ce design.", body)
        self.assertEqual(app.button(key="design_home").label, "Mes designs")
        self.assertEqual(app.button(key="show_design_folder").label, "Afficher le dossier")
        self.assertFalse(app.button(key="prepare_design_export").disabled)
        self.assertEqual([warning.value for warning in app.warning], ["Ancien réglage restauré"])

    def test_creation_requires_name_then_opens_design_with_scoped_output(self):
        app = self.app()
        app.button(key="create_design").click().run()
        self.assertFalse(app.exception)
        self.assertTrue(app.error)
        self.assertFalse(app.session_state.filtered_state.get("active_design_path"))
        project = self.create(app, "Adenosine experiment")
        self.assertEqual(project.name, "Adenosine experiment")
        self.assertEqual(project.path.parent, self.library)
        self.assertEqual(app.session_state["active_design_name"], project.name)
        self.assertEqual(app.radio(key="nav").value, "Conception")
        self.assertEqual(Path(app.session_state["output_dir"]), design_output_dir(project))
        self.assertTrue(design_output_dir(project).is_dir())
        self.assertTrue((project.path / "sequences").is_dir())
        self.assertEqual(app.text_input(key="molecule_name").value, project.name)

    def test_two_designs_keep_separate_sequences_and_reopen_saved_state(self):
        app = self.app()
        first = self.create(app, "First design")
        self.fill_sequences(app)
        app.number_input(key="extension_temperature").set_value(25.5).run()
        app.button(key="design_home").click().run()
        self.assertFalse(app.exception)
        saved = load_design(first.path)
        self.assertEqual(saved.state["form"]["aptamer_input"], "ACGTACGT")
        self.assertEqual(saved.state["form"]["trigger_input"], TRIGGER)
        second = self.create(app, "Second design")
        self.assertNotEqual(second.path, first.path)
        self.assertEqual(app.text_area(key="aptamer_input").value, "")
        self.assertEqual(app.text_area(key="trigger_input").value, "")
        self.assertIsNone(app.session_state["last_result"])
        self.assertIsNone(app.session_state["extension_result"])
        self.fill_sequences(app, "GGGGCCCC", "AACCGGTT")
        app.button(key="design_home").click().run()
        app.selectbox(key="design_to_open").set_value(str(first.path))
        app.button(key="open_design").click().run()
        self.assertFalse(app.exception)
        self.assertEqual(app.text_area(key="aptamer_input").value, "ACGTACGT")
        self.assertEqual(app.text_area(key="trigger_input").value, TRIGGER)
        self.assertEqual(app.number_input(key="extension_temperature").value, 25.5)
        self.assertEqual(Path(app.session_state["output_dir"]), design_output_dir(first))
        self.assertEqual(load_design(second.path).state["form"]["aptamer_input"], "GGGGCCCC")

    def test_existing_name_requires_explicit_overwrite_and_cancel_preserves_it(self):
        original = create_design(self.library, "Existing", state={"form": {"aptamer_input": "AAAACCCC"}})
        marker = design_output_dir(original) / "existing.txt"
        marker.write_text("preserve me", encoding="utf-8")
        app = self.app()
        app.text_input(key="design_new_name").set_value("Existing")
        app.button(key="create_design").click().run()
        self.assertFalse(app.exception)
        self.assertTrue(app.warning)
        self.assertTrue(marker.is_file())
        self.assertEqual(load_design(original.path).state["form"]["aptamer_input"], "AAAACCCC")
        app.button(key="overwrite_design_cancel").click().run()
        self.assertFalse(app.exception)
        self.assertTrue(marker.is_file())
        self.assertFalse(app.session_state.filtered_state.get("design_pending_overwrite"))
        app.button(key="create_design").click().run()
        app.button(key="overwrite_design_confirm").click().run()
        self.assertFalse(app.exception)
        self.assertEqual(app.radio(key="nav").value, "Conception")
        self.assertEqual(app.text_area(key="aptamer_input").value, "")
        self.assertFalse(marker.exists())
        self.assertTrue(any(path.read_text(encoding="utf-8") == "preserve me"
                            for path in self.library.rglob("existing.txt")))

    def test_reopening_legacy_design_keeps_leader_rbs_linker_and_discards_geometry_controls(self):
        project = create_design(self.library, "Editable motifs", state={"form": {
            "aptamer_input": "ACGTACGT", "trigger_input": TRIGGER,
            "customize_tswitch": True, "binding_start": 4, "binding_end": 25,
            "toehold_length": 10, "rbs_aug_distance": 11,
        }})
        app = self.app()
        app.selectbox(key="design_to_open").set_value(str(project.path))
        app.button(key="open_design").click().run()
        app.radio(key="design_section").set_value("switch").run()
        self.assertFalse(app.exception)
        self.assert_no_geometry_controls(app)
        app.text_input(key="t7_leader").set_value("AAA")
        app.text_input(key="rbs_loop").set_value("AGGAGGA")
        app.text_input(key="frame_linker").set_value("CCX").run()
        app.button(key="design_home").click().run()
        saved = load_design(project.path)
        self.assertTrue({"customize_tswitch", "binding_start", "binding_end", "toehold_length",
                         "rbs_aug_distance"}.isdisjoint(saved.state['form']))
        app.selectbox(key="design_to_open").set_value(str(project.path))
        app.button(key="open_design").click().run()
        self.assertFalse(app.exception)
        app.radio(key="design_section").set_value("switch").run()
        self.assert_no_geometry_controls(app)
        self.assertEqual(app.text_input(key="t7_leader").value, "AAA")
        self.assertEqual(app.text_input(key="rbs_loop").value, "AGGAGGA")
        self.assertFalse(app.session_state["linker_auto"])
        self.assertEqual(app.text_input(key="frame_linker").value, "CCX")

    def test_standard_linker_is_editable_persisted_and_used_by_design(self):
        from aptaswitch_core.architecture import CLASSIC_30NT_LINKER
        from aptaswitch_studio.web_app import _build_request
        from aptaswitch_core.thermodynamics.design_constraints import build_design_problem
        app = self.app()
        project = self.create(app, "Editable linker")
        self.fill_sequences(app)
        app.radio(key="design_section").set_value("switch").run()
        self.assert_no_geometry_controls(app)
        self.assertEqual(app.text_input(key="frame_linker").value, CLASSIC_30NT_LINKER)
        self.assertFalse(app.text_input(key="frame_linker").disabled)
        app.text_input(key="frame_linker").set_value(" gct aac ").run()
        self.assertFalse(app.exception)
        state = SimpleNamespace(**app.session_state.filtered_state)
        state.get = vars(state).get
        with patch("aptaswitch_studio.web_app.st.session_state", state), \
                patch("aptaswitch_studio.web_app.configured_library_path", return_value="ui-fixture"), \
                patch("aptaswitch_studio.web_app.native_backend_identity", return_value={"model": "ui-fixture"}):
            request = _build_request()
        self.assertEqual(request.architecture.frame_linker, "GCUAAC")
        # The scientific design domain must receive this sequence, not Green's linker.
        spec = build_design_problem(request)
        switch = next(s for s in spec.strands if s.name == "Switch")
        region = next(r for r in spec.regions if r.name == "frame_linker")
        self.assertEqual(switch.pattern[region.start:region.stop], "GCUAAC")
        app.button(key="design_home").click().run()
        app.selectbox(key="design_to_open").set_value(str(project.path))
        app.button(key="open_design").click().run()
        app.radio(key="design_section").set_value("switch").run()
        self.assertFalse(app.exception)
        self.assert_no_geometry_controls(app)
        self.assertEqual(app.text_input(key="frame_linker").value, " gct aac ")
        # Empty is an explicit choice; default restoration is a separate action.
        app.text_input(key="frame_linker").set_value("").run()
        from aptaswitch_studio.web_app import _current_architecture
        state = SimpleNamespace(**app.session_state.filtered_state)
        state.get = vars(state).get
        with patch("aptaswitch_studio.web_app.st.session_state", state):
            self.assertEqual(_current_architecture(len(TRIGGER)).frame_linker, "")
        app.button(key="restore_architecture_defaults").click().run()
        self.assertFalse(app.exception)
        self.assertTrue(app.session_state["linker_auto"])
        self.assertEqual(app.text_input(key="frame_linker").value, CLASSIC_30NT_LINKER)

    def test_folder_import_copies_complete_design_and_uses_optional_new_name(self):
        source = create_design(self.workspace / "external", "Imported source",
                               state={"form": {"aptamer_input": "ACGTACGT", "trigger_input": TRIGGER}})
        marker = design_output_dir(source) / "engine.log"
        marker.write_text("local completed calculation", encoding="utf-8")
        app = self.app()
        app.text_input(key="design_import_path").set_value(str(source.path))
        app.text_input(key="design_import_name").set_value("Imported copy")
        app.button(key="import_design_path").click().run()
        self.assertFalse(app.exception)
        imported = load_design(Path(app.session_state["active_design_path"]))
        self.assertEqual(imported.name, "Imported copy")
        self.assertEqual(imported.path.parent, self.library)
        self.assertNotEqual(imported.path, source.path)
        self.assertEqual(app.text_area(key="aptamer_input").value, "ACGTACGT")
        self.assertEqual((design_output_dir(imported) / "engine.log").read_text(), marker.read_text())
        self.assertTrue(marker.is_file())

    def test_extended_sequence_and_originals_survive_reopen_and_can_be_restored(self):
        app = self.app()
        project = self.create(app, "Extended trigger")
        fixture = APP_PATH.parents[2] / "tests" / "fixtures" / "extension_native_fixture.json"
        result = json.loads(fixture.read_text())
        original = result["input"]
        self.fill_sequences(app, original["aptamer_dna"], original["trigger_dna"])
        result["output_dir"] = str(design_output_dir(project) / "extensions" / result["run_id"])
        exports = export_extension_result(result)
        app.session_state["extension_result"] = result
        app.session_state["extension_exports"] = exports
        app.session_state["results_kind"] = "Extension de triggers"
        app.radio(key="nav").set_value("Résultats").run()
        app.button(key="extension_use_1").click().run()
        self.assertFalse(app.exception)
        self.assertEqual(app.radio(key="design_section").value, "switch")
        app.radio(key="design_section").set_value("system").run()
        extended = app.text_area(key="trigger_input").value
        self.assertNotEqual(extended, original["trigger_dna"])
        app.button(key="design_home").click().run()
        app.selectbox(key="design_to_open").set_value(str(project.path))
        app.button(key="open_design").click().run()
        self.assertFalse(app.exception)
        self.assertEqual(app.text_area(key="trigger_input").value, extended)
        self.assertEqual(app.session_state["design_original_trigger"], original["trigger_dna"])
        app.button(key="restore_original_sequences").click().run()
        self.assertFalse(app.exception)
        self.assertEqual(app.text_area(key="trigger_input").value, original["trigger_dna"])
        self.assertEqual(app.text_area(key="aptamer_input").value, original["aptamer_dna"])
        self.assertFalse(app.session_state["design_using_extension"])
        self.assertTrue(Path(exports["json"]).is_file())

    def test_import_collision_replaces_active_project_only_after_confirmation(self):
        app = self.app()
        active = self.create(app, "Same name")
        self.fill_sequences(app, "AAAACCCC", TRIGGER)
        source = create_design(self.workspace / "incoming", "Same name",
                               state={"form": {"aptamer_input": "GGGGCCCC", "trigger_input": "AACCGGTT"}})
        app.button(key="design_home").click().run()
        app.text_input(key="design_import_path").set_value(str(source.path))
        app.button(key="import_design_path").click().run()
        self.assertFalse(app.exception)
        self.assertEqual(load_design(active.path).state["form"]["aptamer_input"], "AAAACCCC")
        app.button(key="overwrite_design_confirm").click().run()
        self.assertFalse(app.exception)
        self.assertEqual(app.session_state["active_design_path"], str(active.path))
        self.assertEqual(app.text_area(key="aptamer_input").value, "GGGGCCCC")
        self.assertEqual(load_design(active.path).state["form"]["aptamer_input"], "GGGGCCCC")
        self.assertEqual(load_design(source.path).state["form"]["aptamer_input"], "GGGGCCCC")

    def test_zip_export_contains_project_and_can_be_reimported_from_path(self):
        app = self.app()
        project = self.create(app, "Portable project")
        self.fill_sequences(app)
        app.button(key="prepare_design_export").click().run()
        self.assertFalse(app.exception)
        self.assertTrue(any(item.key == "export_design" for item in app.get("download_button")))
        data = app.session_state["design_export_bytes"]
        with zipfile.ZipFile(io.BytesIO(data)) as archive:
            manifest = next(name for name in archive.namelist() if name.endswith("design.json"))
            state = json.loads(archive.read(manifest))["state"]
            self.assertEqual(state["form"]["aptamer_input"], "ACGTACGT")
            self.assertEqual(state["original_sequences"]["trigger"], TRIGGER)
        path = self.workspace / "project.zip"
        path.write_bytes(data)
        app.button(key="design_home").click().run()
        app.text_input(key="design_import_path").set_value(str(path))
        app.text_input(key="design_import_name").set_value("Portable copy")
        app.button(key="import_design_path").click().run()
        self.assertFalse(app.exception)
        self.assertNotEqual(app.session_state["active_design_path"], str(project.path))
        self.assertEqual(app.session_state["active_design_name"], "Portable copy")
        self.assertEqual(app.text_area(key="trigger_input").value, TRIGGER)

    def test_running_calculation_blocks_creating_opening_and_importing_another_design(self):
        app = self.app()
        project = self.create(app, "Working project")
        other = create_design(self.library, "Other project")
        running = SimpleNamespace(request=SimpleNamespace(engine="native-hybrid"),
                                  snapshot=Mock(return_value=JobSnapshot(state="running", current=0, total=1, message="Running")))
        app.session_state["design_job"] = running
        app.button(key="design_home").click().run()
        self.assertFalse(app.exception)
        for key in ("create_design", "open_design", "import_design_path"):
            self.assertTrue(app.button(key=key).disabled)
        self.assertEqual(app.session_state["active_design_path"], str(project.path))
        self.assertNotEqual(app.session_state["active_design_path"], str(other.path))


if __name__ == "__main__":
    unittest.main()
