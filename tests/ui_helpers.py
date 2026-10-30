"""Explicit French sessions for UI assertions written before localization.

New localization tests use Streamlit's AppTest directly to exercise the English
application default. Existing French UI assertions continue to cover that locale.
"""

from contextlib import ExitStack
import importlib
from pathlib import Path
import tempfile
import weakref
from unittest.mock import patch

from streamlit.testing.v1 import AppTest as StreamlitAppTest


def prepare_test_design(app, *, project=None):
    """Give a full-app fixture a real, isolated design without testing onboarding.

    New landing-page tests use raw AppTest instead. All run entry points,
    including widget.run(), receive the temporary library root. The project's
    lifetime follows the AppTest instance and never uses a user's design tree.
    """
    from aptaswitch_studio.design_library import create_design, design_output_dir

    if project is None:
        owner = tempfile.TemporaryDirectory(prefix="aptaswitch-ui-design-")
        app._test_design_temporary = owner
        app._test_design_cleanup = weakref.finalize(app, owner.cleanup)
        project = create_design(Path(owner.name), "Test design")
    app._test_design_project = project
    app._test_design_library_root = Path(project.path).parent
    app.session_state["active_design_path"] = str(project.path)
    app.session_state["active_design_name"] = project.name
    app.session_state["output_dir"] = str(design_output_dir(project))
    app.session_state["nav"] = "Conception"
    if not getattr(app, "_test_design_run_wrapped", False):
        original_run = app._run

        def isolated_run(*args, **kwargs):
            root = app._test_design_library_root
            with ExitStack() as context:
                for name in ("aptaswitch_studio.design_library", "aptaswitch_studio.design_library_ui"):
                    module = importlib.import_module(name)
                    if hasattr(module, "default_design_root"):
                        context.enter_context(patch.object(module, "default_design_root", return_value=root))
                return original_run(*args, **kwargs)

        app._run = isolated_run
        app._test_design_run_wrapped = True
    return design_output_dir(project)


def candidate_tables(app):
    """Ranking tables, independently of optimizer or provenance summaries."""
    return [table for table in app.get("table") if "Score" in table.value.columns]


class FrenchAppTest(StreamlitAppTest):
    @staticmethod
    def _french(app):
        app.session_state["ui_language"] = "fr"
        return app

    @classmethod
    def from_file(cls, *args, **kwargs):
        app = cls._french(super().from_file(*args, **kwargs))
        path = args[0] if args else kwargs.get("script_path", "")
        if Path(path).name == "web_app.py":
            prepare_test_design(app)
        return app

    @classmethod
    def from_string(cls, *args, **kwargs):
        return cls._french(super().from_string(*args, **kwargs))

    @classmethod
    def from_function(cls, *args, **kwargs):
        return cls._french(super().from_function(*args, **kwargs))

    @classmethod
    def _from_string(cls, script, **kwargs):
        # Component-only apps do not call web_app.main(), so activate their
        # locale inside the Streamlit script thread as well as session state.
        prefix = "from aptaswitch_core.localization import set_language\nset_language('fr')\n"
        return super()._from_string(prefix + script, **kwargs)
