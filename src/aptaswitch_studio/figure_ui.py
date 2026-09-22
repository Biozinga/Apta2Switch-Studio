"""One display and export component for all molecular diagrams.

The established switch view is the reference: fitted 2D figures, natural-size
scrollable linear maps, and the same SVG/PNG download controls throughout.
"""

from __future__ import annotations

import base64
from collections.abc import Callable, Sequence
from dataclasses import dataclass
from html import escape

import streamlit as st

from aptaswitch_core.localization import tr

from .web_runtime import svg_to_png_bytes

MOLECULAR_VIEW_LABELS = ("Structure 2D", "Séquence linéaire")
Diagram = str | tuple[str, int, int]


@dataclass(frozen=True)
class ExportOption:
    label: str
    data: str | bytes
    filename: str
    mime: str | None = None
    key: str | None = None


def render_export_menu(options: Sequence[ExportOption], *, label: str = "Exporter") -> None:
    """Keep every format/file choice behind one consistent export control."""
    if not options:
        return
    with st.popover(tr(label), icon=":material/download:"):
        for option in options:
            st.download_button(
                tr(option.label), option.data, option.filename, option.mime,
                key=option.key, width="stretch", on_click="ignore",
            )


def svg_export_options(
    svg: str, *, filename: str, key: str = "",
    download_labels: tuple[str, str] = ("Télécharger la structure SVG", "Télécharger la structure PNG"),
) -> tuple[ExportOption, ExportOption]:
    return (
        ExportOption(download_labels[0], svg.encode("utf-8"), f"{filename}.svg", "image/svg+xml", f"{key}_svg" if key else None),
        ExportOption(download_labels[1], svg_to_png_bytes(svg), f"{filename}.png", "image/png", f"{key}_png" if key else None),
    )


def render_svg_exports(svg: str, *, filename: str, key: str = "",
                       download_labels: tuple[str, str] = ("Télécharger la structure SVG", "Télécharger la structure PNG")) -> None:
    render_export_menu(svg_export_options(svg, filename=filename, key=key, download_labels=download_labels))


def render_structure_copy(structure: str) -> None:
    """Streamlit's code block provides the native clipboard action."""
    with st.popover("Structure", icon=":material/content_copy:", help=tr("Afficher et copier la notation dot-bracket")):
        st.code(structure, language=None, wrap_lines=True)


def render_svg(svg: str, *, height: int, natural_width: int | None = None, scroll: bool = False) -> None:
    # An image data URI avoids Streamlit's inline-SVG sanitization and preserves
    # the original switch viewer's sizing and horizontal scrolling behavior.
    encoded = base64.b64encode(svg.encode("utf-8")).decode("ascii")
    image_css = (
        f"width:{max(1, int(natural_width))}px;height:auto;max-width:none;max-height:100%;object-fit:contain;"
        if natural_width
        else "width:100%;height:100%;max-width:100%;max-height:100%;object-fit:contain;"
    )
    overflow_x = "auto" if scroll else "hidden"
    st.html(f"""
    <div style="width:100%;height:{max(1, int(height))}px;overflow-x:{overflow_x};overflow-y:hidden;
      border-radius:12px;background:#fff;display:flex;align-items:center;">
      <img src="data:image/svg+xml;base64,{encoded}" alt="{escape(tr('Aperçu moléculaire'), quote=True)}"
        style="{image_css}display:block;margin:0 auto;" />
    </div>
    """)


def render_diagram(
    diagram: Callable[..., Diagram],
    *,
    height: int = 510,
    legend: Sequence[tuple[str, str]] = (),
    caption: str = "",
    filename: str | None = None,
    key: str = "",
    download_labels: tuple[str, str] = ("Télécharger la structure SVG", "Télécharger la structure PNG"),
    structure: str | None = None,
) -> None:
    """Render one molecular figure in live and export modes using the same function.

    A renderer is usually a ``functools.partial`` binding the scientific input.
    No folding calculation is repeated: only the SVG layout is drawn twice.
    """
    live = diagram(show_annotations=False)
    if isinstance(live, tuple):
        svg, width, linear_height = live
        render_svg(svg, height=max(250, linear_height + 24), natural_width=width, scroll=True)
    else:
        svg = live
        render_svg(svg, height=height)
    if legend:
        chips = "".join(
            f'<span class="legend-chip"><span class="legend-dot" style="background:{escape(color, quote=True)}"></span>{escape(tr(label))}</span>'
            for label, color in legend
        )
        st.markdown(chips, unsafe_allow_html=True)
    if caption:
        st.caption(tr(caption))
    if structure or filename:
        with st.container(horizontal=True, vertical_alignment="center"):
            if structure:
                render_structure_copy(structure)
            if filename:
                exported = diagram(show_annotations=True)
                export_svg = exported[0] if isinstance(exported, tuple) else exported
                render_svg_exports(export_svg, filename=filename, key=key, download_labels=download_labels)
