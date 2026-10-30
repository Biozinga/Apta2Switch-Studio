"""Architecture selector with genuinely disabled, forthcoming choices."""

from __future__ import annotations

from html import escape

import streamlit as st

from aptaswitch_core.localization import tr


def render_architecture_selector() -> None:
    """Show the supported architecture and unavailable alternatives.

    Streamlit's selectbox only supports disabling the entire control. A native
    select preserves keyboard accessibility and disables individual options.
    Since only the activator is selectable, no client-to-server state bridge is
    needed; add one when another architecture is actually implemented.
    """
    st.markdown(tr("Architecture du switch"),
                help=tr("Les architectures répressives ne sont pas encore disponibles."))
    markup = """
    <div class="switch-architecture-field">
      <style>
        .switch-architecture-field { width: 100%; }
        .switch-architecture-field select {
          box-sizing: border-box; width: 100%; min-height: 2.5rem;
          padding: 0.5rem 0.75rem; border: 1px solid transparent;
          border-radius: 0.5rem; font: inherit; font-size: 1rem;
          color: var(--text-color, #15233a);
          background-color: var(--secondary-background-color, #eaf1f8);
          cursor: pointer;
        }
        .switch-architecture-field select:focus-visible {
          outline: 2px solid var(--primary-color, #1258dc);
          outline-offset: 2px;
        }
        .switch-architecture-field option:disabled { color: #8c95a3; }
      </style>
      <select id="switch-architecture-choice" name="switch-architecture"
              aria-label="Architecture du switch">
        <option value="activator" selected>Toehold activateur</option>
        <option value="repressor" disabled>Toehold répresseur</option>
        <option value="repressor_3wj" disabled>Toehold répresseur 3WJ</option>
      </select>
    </div>
    """
    for source in (
        "Architecture du switch", "Toehold activateur", "Toehold répresseur 3WJ",
        "Toehold répresseur",
    ):
        markup = markup.replace(source, escape(tr(source)))
    st.html(markup)
