"""Explicit French sessions for UI assertions written before localization.

New localization tests use Streamlit's AppTest directly to exercise the English
application default. Existing French UI assertions continue to cover that locale.
"""

from streamlit.testing.v1 import AppTest as StreamlitAppTest


class FrenchAppTest(StreamlitAppTest):
    @staticmethod
    def _french(app):
        app.session_state["ui_language"] = "fr"
        return app

    @classmethod
    def from_file(cls, *args, **kwargs):
        return cls._french(super().from_file(*args, **kwargs))

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
