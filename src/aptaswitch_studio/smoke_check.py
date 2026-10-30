"""Internal release check, executed by the packaged application itself."""

from __future__ import annotations

import json
import hashlib
import math
import os
import subprocess
import sys
import tempfile
import traceback
from contextlib import contextmanager
from io import BytesIO
from pathlib import Path
from unittest.mock import patch


@contextmanager
def _isolated_ui_settings(directory: Path):
    """Keep settings reads and writes isolated, including imported aliases."""
    from aptaswitch_core import settings

    directory = directory.resolve()
    settings_dir = directory / "settings"
    with (
        patch.object(settings, "SETTINGS_DIR", settings_dir),
        patch.object(settings, "SETTINGS_PATH", settings_dir / "settings.json"),
    ):
        settings.save_settings({
            "design_library_dir": str(directory / "designs"),
            "default_output_dir": str(directory / "runs"),
        })
        yield


def _check_rbs_runtime() -> dict:
    """Exercise the shipped OSTIR tables and RNA binding with a constrained ON."""
    from aptaswitch_core import rbs_calculator

    # This synthetic competitor closes the RBS in OFF. Masking its six bases
    # must release it in ON, so a silently ignored constraint cannot pass.
    sequence = "GGGCUCCUUAAAAAAGGAGAAAAAAAUG" + "GCG" * 15
    aug_start = sequence.index("AUG")
    occupied = tuple(range(3, 9))
    prediction = rbs_calculator.predict_rbs_switch(sequence, aug_start, occupied)
    assert prediction["reference"]["ostir_version"] == "1.1.3"
    assert prediction["reference"]["vienna_version"] == "2.7.2"
    assert all(math.isfinite(prediction[state]["dG_total"]) for state in ("off", "on"))
    assert all(math.isfinite(prediction[key]) for key in ("off_rate", "on_rate", "merit"))
    assert prediction["on_rate"] > prediction["off_rate"] > 0
    assert 0 < prediction["merit"] <= 1
    # Historical validation must preserve the exact RNA/mask and use only the
    # recorded values; it must not need a second native calculation.
    with patch.object(rbs_calculator, "_runtime", side_effect=AssertionError("Unexpected RBS recalculation")):
        rbs_calculator.validate_rbs_prediction(prediction, sequence, aug_start, occupied)
    return {key: prediction[key] for key in (
        "reference", "on_constraints", "off_rate", "on_rate", "merit",
    )}


def run_smoke_check(output_path: str) -> int:
    """Check shipped dependencies and every main page without running a design."""
    from aptaswitch_studio import __version__

    report: dict = {"ok": False, "frozen": bool(getattr(sys, "frozen", False)),
                    "version": __version__, "checks": []}
    try:
        import jinja2
        import numpy
        import pandas
        import scipy.interpolate
        import scipy.optimize
        import scipy.sparse
        import yaml
        from openpyxl import Workbook, load_workbook
        from packaging.tags import sys_tags
        from streamlit.testing.v1 import AppTest

        assert next(sys_tags()).interpreter == "cp311", "The packaged scientific runtime targets Python 3.11"
        assert scipy.optimize.minimize_scalar(lambda x: (x - 2) ** 2).success
        assert scipy.sparse.eye(2).shape == (2, 2)
        assert yaml.safe_load("ready: true")["ready"]
        assert jinja2.Template("{{ ready }}").render(ready="ok") == "ok"
        assert pandas.DataFrame({"value": numpy.arange(2)}).shape == (2, 1)
        report["checks"].append("bundled scientific dependencies")

        from aptaswitch_core.thermodynamics.setup import LIBRARY_BASENAME, validate_runtime
        bundle = Path(getattr(sys, "_MEIPASS", Path(__file__).resolve().parents[2]))
        choices = [bundle / "thermodynamics" / LIBRARY_BASENAME,
                   bundle.parent / "Resources/thermodynamics" / LIBRARY_BASENAME]
        if not getattr(sys, "frozen", False):
            choices.append(bundle / ".local_deps/thermodynamics" / LIBRARY_BASENAME)
        native_library = next((path for path in choices if path.is_file()), None)
        assert native_library is not None, "The packaged native hybrid library is missing"
        native = validate_runtime(native_library)
        assert native.ok, native.message
        report["thermodynamics"] = native.provenance
        report["checks"].append("bundled hybrid runtime: RNA/DNA MFE, probabilities and joint samples")
        # The smoke process must never validate an unrelated remembered engine.
        os.environ["APTASWITCH_HYBRID_NATIVE_LIBRARY"] = str(native_library)

        report["rbs_calculator"] = _check_rbs_runtime()
        report["checks"].append("bundled OSTIR/ViennaRNA: finite OFF and constrained ON, exact occupancy mask and saved provenance")

        from aptaswitch_core.oxdna_setup import detect_oxdna_runtime

        ox_runtime = detect_oxdna_runtime()
        assert ox_runtime is not None, "The bundled oxDNA runtime is missing"
        report["oxdna"] = {key: ox_runtime[key] for key in ("source_revision", "oxpy_version", "license")}
        report["checks"].append("bundled oxDNA/oxpy modelNA")

        from aptaswitch_studio.web_runtime import svg_to_png_bytes

        png = svg_to_png_bytes('<svg xmlns="http://www.w3.org/2000/svg" width="24" height="24">'
                               '<circle cx="12" cy="12" r="10" fill="#1258DC"/></svg>')
        assert png.startswith(b"\x89PNG\r\n\x1a\n")
        workbook = Workbook()
        workbook.active.append(["sequence", "ACGU"])
        stream = BytesIO()
        workbook.save(stream)
        stream.seek(0)
        assert load_workbook(stream).active["B1"].value == "ACGU"
        report["checks"].append("PNG and Excel exports")

        from aptaswitch_core.architecture import architecture_layout, default_architecture_for_trigger
        from aptaswitch_core.preview import build_switch_model

        classic = default_architecture_for_trigger(30)
        layout = architecture_layout(classic)
        preview = build_switch_model(classic, "ACGT" * 7 + "AC", reporter_sequence="CGUAAAGGC")
        linker_start = layout["aug_start"] + 3 + 9
        assert preview.switch_sequence[linker_start:layout["reporter_start"]] == "AACCUGGCGGCAGCGCAAAAG"
        assert preview.switch_sequence[layout["reporter_start"]:] == "CGUAAAGGC"
        assert preview.off_structure[linker_start:layout["reporter_start"]] == "." * 21
        assert layout["reporter_start"] - layout["aug_start"] == 33
        assert next(domain for domain in preview.domains if domain.label == "Linker").sequence == "AACCUGGCGGCAGCGCAAAAG"
        report["checks"].append("classic 30-nt activator: exact 21-nt linker after the 9-nt downstream stem, reporter in frame")

        app_path = Path(__file__).with_name("web_app.py")
        assets = app_path.parent / "assets"
        for filename in ("logo.png", "app_icon.png", "igem_SU_logo.png"):
            assert (assets / filename).is_file(), filename
        from aptaswitch_studio.molecular_viewer_assets import bundled_viewer_javascript

        viewer_assets = assets / "3dmol"
        viewer_manifest = json.loads((viewer_assets / "manifest.json").read_text(encoding="utf-8"))
        for filename, entry in viewer_manifest["files"].items():
            assert hashlib.sha256((viewer_assets / filename).read_bytes()).hexdigest() == entry["sha256"], filename
        assert len(bundled_viewer_javascript()) > 500_000
        report["checks"].append("embedded offline oxDNA molecular viewer, vendor checksums and licenses")
        from aptaswitch_studio.design_library import design_output_dir, load_design

        # Exercise onboarding instead of bypassing the named-design gate. All
        # manifests and autosaved inputs stay in the disposable smoke library.
        with tempfile.TemporaryDirectory(prefix="aptaswitch-design-smoke-") as directory:
            library = Path(directory).resolve() / "designs"
            with _isolated_ui_settings(Path(directory)):
                app = AppTest.from_file(str(app_path), default_timeout=45)
                app.session_state["dependency_detection_done"] = True
                app.session_state["ui_language"] = "en"
                app.run()
                assert not app.exception, [str(error.value) for error in app.exception]
                assert app.radio(key="nav").value == "Mes designs"
                assert not app.session_state.filtered_state.get("active_design_path")
                assert not any(area.key == "aptamer_input" for area in app.text_area)
                app.text_input(key="design_new_name").set_value("Smoke design")
                app.button(key="create_design").click().run()
                assert not app.exception, [str(error.value) for error in app.exception]
                assert app.radio(key="nav").value == "Conception"
                project = load_design(Path(app.session_state["active_design_path"]))
                assert project.path.parent == library and project.name == "Smoke design"
                assert Path(app.session_state["output_dir"]) == design_output_dir(project)
                app.text_area(key="aptamer_input").set_value("ACGTACGT")
                app.text_area(key="trigger_input").set_value("ACGT" * 6).run()
                assert not app.exception, [str(error.value) for error in app.exception]
                app.button(key="design_home").click().run()
                saved = load_design(project.path)
                assert saved.state["form"]["aptamer_input"] == "ACGTACGT"
                assert saved.state["original_sequences"]["trigger"] == "ACGT" * 6
                app.selectbox(key="design_to_open").set_value(str(project.path))
                app.button(key="open_design").click().run()
                assert not app.exception, [str(error.value) for error in app.exception]
                assert app.text_area(key="aptamer_input").value == "ACGTACGT"
                assert app.text_area(key="trigger_input").value == "ACGT" * 6
                report["checks"].append("named local designs: empty home, creation, isolated autosave and reopen")

                for language in ("en", "fr"):
                    app.session_state["ui_language"] = language
                    for page in ("Mes designs", "Conception", "Résultats", "Toehold switches enregistrés", "Réglages"):
                        app.session_state["nav"] = page
                        app.run()
                        failures = [str(error.value) for error in app.exception]
                        assert not failures, f"{language}/{page}: {failures}"
                        assert app.radio(key="nav").value == page
                        assert len(app.selectbox) + len(app.button) > 0
                        if page == "Conception":
                            app.radio(key="design_section").set_value("switch").run()
                            assert not app.exception, [str(error.value) for error in app.exception]
                            for key in ("t7_leader", "rbs_loop", "frame_linker"):
                                assert not app.text_input(key=key).disabled, key
                            removed = {"customize_tswitch", "toehold_length", "rbs_aug_distance",
                                       "binding_start", "binding_end"}
                            assert not removed.intersection(
                                widget.key for controls in (app.checkbox, app.slider, app.number_input)
                                for widget in controls
                            )
                        if page == "Réglages":
                            from aptaswitch_core.localization import language_context, tr
                            with language_context(language):
                                assert any(tab.label == tr("Licences et remerciements") for tab in app.tabs)
                report["checks"].append("My designs, Design, Results, Saved toehold switches and Licences in English and French")
                report["checks"].append("standard architecture only; editable leader, RBS and linker in English and French")
        report["ok"] = True
    except Exception:
        report["error"] = traceback.format_exc()
    Path(output_path).write_text(json.dumps(report, indent=2), encoding="utf-8")
    return 0 if report["ok"] else 1
