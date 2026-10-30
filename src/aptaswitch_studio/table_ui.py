"""Static result tables: no cell selection, editing or library toolbar."""

import streamlit as st


def render_readonly_table(data, *, height: int | None = None, key: str | None = None) -> None:
    """Keep table presentation uniform; navigation and exports live outside it.

    Accept the same data and pandas Styler objects as ``st.table``. Callers own
    scientific precision, labels and any explicit preview/subset of their data.
    A bounded height enables scrolling without dropping rows.
    """
    with st.container(key=key):
        st.table(data, hide_index=True, border="horizontal", width="stretch",
                 height=height if height is not None else "content")
