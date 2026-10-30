"""Live 2D secondary-structure preview for toehold switches.

This module builds a switch design straight from the runtime architecture,
without a thermodynamics engine. Bases selected for optimization are shown as
``X`` placeholders. It exposes:

* :func:`build_switch_model` - domains, sequence, colors and OFF/ON structures.
* :func:`radial_layout` - a pure-Python radial layout of a dot-bracket string.
* :func:`render_structure_svg` - a standalone SVG string for one structure.
* :func:`switch_preview_svg` - convenience wrapper used by the web UI.

The layout is a classic loop/helix radial drawing: every loop becomes a regular
polygon and every helix a straight ladder. It handles multi-strand structures
written with a ``+`` strand break (used for the trigger + switch ON complex).
"""

from __future__ import annotations

from .localization import tr

import math
from dataclasses import dataclass
from html import escape
from typing import List, Optional, Sequence, Tuple

from .architecture import architecture_layout, target_switch_structures, trigger_binding_slices
from .models import SwitchArchitecture
from .sequences import OPTIMIZABLE_BASE, SequenceValidationError, reading_frame_report, reverse_complement_rna_pattern

Point = Tuple[float, float]
Pair = Tuple[int, int]

# Domain colors preserve the original switch figures.
DOMAIN_COLORS = {
    "T7": "#cfd8dc",
    "Toehold": "#33b5e5",
    "Stem L": "#1976d2",
    "U2 L": "#7e57c2",
    "Bulge": "#ffb74d",
    "U L": "#5e35b1",
    "Boucle RBS": "#ffb74d",
    "RBS": "#ef5350",
    "U R": "#5e35b1",
    "AUG": "#00a676",
    "U2 R": "#7e57c2",
    "Stem R": "#1565c0",
    "Linker": "#8bc34a",
    "Reporter": "#1db954",
    "Spacer": "#90a4ae",
}
TRIGGER_COLOR = "#ffd54f"
APTAMER_COLOR = DOMAIN_COLORS["Stem L"]
EXTENSION_COLOR = DOMAIN_COLORS["AUG"]
EXTENSION_3PRIME_COLOR = "#8b5cf6"
LIGAND_COLOR = "#e85d9e"
PLACEHOLDER_BASE = OPTIMIZABLE_BASE

_DARK_FILLS = {
    "#1976d2",
    "#7e57c2",
    "#5e35b1",
    "#ef5350",
    "#00a676",
    "#1565c0",
    "#1db954",
    "#33b5e5",
}


@dataclass(frozen=True)
class SwitchDomain:
    """One contiguous stretch of the switch with a semantic role."""

    label: str
    sequence: str
    color: str
    optimizable: bool


@dataclass(frozen=True)
class SwitchModel:
    """Everything needed to draw the OFF and ON views of a design."""

    domains: List[SwitchDomain]
    switch_sequence: str
    switch_colors: List[str]
    off_structure: str
    on_structure: str
    trigger_sequence: str
    reporter_label: str

    @property
    def on_sequence(self) -> str:
        return self.trigger_sequence + self.switch_sequence

    @property
    def on_colors(self) -> List[str]:
        return [TRIGGER_COLOR] * len(self.trigger_sequence) + self.switch_colors


def _frame_linker(architecture: SwitchArchitecture) -> str:
    return architecture.resolved_frame_linker()


def build_switch_model(
    architecture: SwitchArchitecture,
    trigger_dna: Optional[str] = None,
    reporter_sequence: str = "CGUAAAGGCGAGGAGCUGUUC",
    reporter_label: str = "sfGFP",
) -> SwitchModel:
    """Build the switch design (with ``X`` placeholders) from the architecture.

    The recognized domains are derived from the reverse complement of the
    trigger when one is supplied; otherwise they are shown as ``X`` too. Every
    domain selected for optimization is always rendered as ``X``. The switch
    is RNA and the trigger is DNA.
    """

    arch = architecture
    arch.validate_sequence_patterns()
    toe_len = arch.toehold_length
    stem_len = arch.stem_length

    binding: Optional[str] = None
    if trigger_dna:
        binding = reverse_complement_rna_pattern(trigger_dna)

    binding_slices = trigger_binding_slices(arch)
    paired_region = max(section.stop for section in binding_slices.values())
    if (binding is not None and arch.binding_layout == "full_left_arm"
            and len(binding) != paired_region):
        raise ValueError(tr('La région reconnue doit correspondre à la longueur du trigger ({v0} != {v1}).',
                            v0=paired_region, v1=len(binding)))
    if binding is not None and len(binding) < paired_region:
        binding = None

    def binding_domain(name: str, label: str, length: int) -> SwitchDomain:
        section = binding_slices.get(name)
        fixed = binding is not None and section is not None
        sequence = binding[section] if fixed else PLACEHOLDER_BASE * length
        return SwitchDomain(label, sequence, DOMAIN_COLORS[label], not fixed or PLACEHOLDER_BASE in sequence)

    upper = arch.upper_stem_length
    upper2 = arch.upper_stem2_length
    bulge = arch.bulge_length
    rbs_prefix_length = arch.rbs_prefix_length
    aug_spacer = max(0, int(arch.aug_spacer))
    t7_seq = arch.t7_leader.upper().replace("T", "U")
    rbs_seq = arch.rbs_loop.upper().replace("T", "U")
    linker_seq = _frame_linker(arch).replace("T", "U")
    reporter_seq = (reporter_sequence or "").upper().replace("T", "U")
    aug_seq = "AUG"

    n = PLACEHOLDER_BASE
    domains = [
        SwitchDomain("Leader", t7_seq, DOMAIN_COLORS["T7"], False),
        binding_domain("toehold", "Toehold", toe_len),
        binding_domain("stem_l", "Stem L", stem_len),
    ]
    if upper2:
        domains.append(binding_domain("upper2_l", "U2 L", upper2))
    domains += [
        binding_domain("bulge_l", "Bulge", bulge),
        binding_domain("upper_l", "U L", upper),
    ]
    if rbs_prefix_length:
        domains.append(SwitchDomain("Boucle RBS", n * rbs_prefix_length, DOMAIN_COLORS["Boucle RBS"], True))
    domains += [
        SwitchDomain("RBS", rbs_seq, DOMAIN_COLORS["RBS"], False),
        SwitchDomain("U R", n * upper, DOMAIN_COLORS["U R"], True),
    ]
    if aug_spacer:
        domains.append(SwitchDomain("Spacer", n * aug_spacer, DOMAIN_COLORS["Spacer"], True))
    domains.append(SwitchDomain("AUG", aug_seq, DOMAIN_COLORS["AUG"], False))
    if upper2:
        domains.append(SwitchDomain("U2 R", n * upper2, DOMAIN_COLORS["U2 R"], True))
    domains.append(SwitchDomain("Stem R", n * stem_len, DOMAIN_COLORS["Stem R"], True))
    if linker_seq:
        domains.append(SwitchDomain("Linker", linker_seq, DOMAIN_COLORS["Linker"], PLACEHOLDER_BASE in linker_seq))
    domains.append(SwitchDomain(reporter_label, reporter_seq, DOMAIN_COLORS["Reporter"], False))

    switch_sequence = "".join(domain.sequence for domain in domains)
    switch_colors: List[str] = []
    for domain in domains:
        switch_colors.extend([domain.color] * len(domain.sequence))

    reporter_len = len(reporter_seq)
    off_structure, on_structure = target_switch_structures(arch, reporter_len)

    trigger_sequence = (trigger_dna or PLACEHOLDER_BASE * paired_region).upper().replace("U", "T")

    return SwitchModel(
        domains=domains,
        switch_sequence=switch_sequence,
        switch_colors=switch_colors,
        off_structure=off_structure,
        on_structure=on_structure,
        trigger_sequence=trigger_sequence,
        reporter_label=reporter_label,
    )


# ---------------------------------------------------------------------------
# Structure parsing and radial layout
# ---------------------------------------------------------------------------


def parse_structure(structure: str) -> Tuple[List[Pair], int, set]:
    """Return base pairs, base count and strand-break positions.

    ``cuts`` holds 0-based indices ``i`` with no backbone bond between ``i`` and
    ``i + 1`` (a ``+`` strand break in the dot-bracket notation).
    """

    stacks = {"(": [], "[": [], "{": [], "<": []}
    closers = {")": "(", "]": "[", "}": "{", ">": "<"}
    pairs: List[Pair] = []
    cuts: set = set()
    index = 0

    for char in structure:
        if char in " \t\n":
            continue
        if char == "+":
            if index > 0:
                cuts.add(index - 1)
            continue
        if char == ".":
            index += 1
            continue
        if char in stacks:
            stacks[char].append(index)
            index += 1
            continue
        if char in closers:
            opener = closers[char]
            if not stacks[opener]:
                raise ValueError(tr('Structure déséquilibrée : {v0}', v0=structure))
            pairs.append((stacks[opener].pop(), index))
            index += 1
            continue
        raise ValueError(tr('Caractère de structure non pris en charge : {v0!r}', v0=char))

    if any(stacks.values()):
        raise ValueError(tr('Structure déséquilibrée : {v0}', v0=structure))
    return pairs, index, cuts


def _norm(vx: float, vy: float) -> Point:
    length = math.hypot(vx, vy)
    if length < 1e-9:
        return (0.0, 1.0)
    return (vx / length, vy / length)


def _wrap(angle: float) -> float:
    while angle > math.pi:
        angle -= 2 * math.pi
    while angle <= -math.pi:
        angle += 2 * math.pi
    return angle


def radial_layout(structure: str, spacing: float = 1.0) -> List[Point]:
    """Compute 2D coordinates for every base of a dot-bracket structure."""

    pairs, n, _ = parse_structure(structure)
    if n == 0:
        return []

    partner = [-1] * n
    for i, j in pairs:
        partner[i] = j
        partner[j] = i
    coords: List[Optional[Point]] = [None] * n

    def place_helix(a: int, b: int, out_dir: Point) -> None:
        ax, ay = coords[a]  # type: ignore[misc]
        bx, by = coords[b]  # type: ignore[misc]
        i, j = a, b
        while i + 1 < j - 1 and partner[i + 1] == j - 1:
            ax += out_dir[0] * spacing
            ay += out_dir[1] * spacing
            bx += out_dir[0] * spacing
            by += out_dir[1] * spacing
            coords[i + 1] = (ax, ay)
            coords[j - 1] = (bx, by)
            i += 1
            j -= 1
        place_loop(i, j, out_dir)

    def place_loop(i: int, j: int, into_dir: Point) -> None:
        verts = [i]
        kids: List[Pair] = []
        k = i + 1
        while k < j:
            if partner[k] == -1:
                verts.append(k)
                k += 1
            else:
                b = partner[k]
                verts.append(k)
                verts.append(b)
                kids.append((k, b))
                k = b + 1
        verts.append(j)
        m = len(verts)
        if m <= 2:
            return

        pi = coords[i]  # type: ignore[assignment]
        pj = coords[j]  # type: ignore[assignment]
        radius = spacing / (2 * math.sin(math.pi / m))
        apothem = math.sqrt(max(radius * radius - (spacing / 2) ** 2, 0.0))
        mid = ((pi[0] + pj[0]) / 2, (pi[1] + pj[1]) / 2)
        nd = _norm(into_dir[0], into_dir[1])
        center = (mid[0] + apothem * nd[0], mid[1] + apothem * nd[1])

        theta0 = math.atan2(pi[1] - center[1], pi[0] - center[0])
        phi = math.atan2(pj[1] - center[1], pj[0] - center[0])
        sign = -1.0 if _wrap(phi - theta0) > 0 else 1.0
        step = 2 * math.pi / m
        for t, v in enumerate(verts):
            angle = theta0 + sign * t * step
            coords[v] = (
                center[0] + radius * math.cos(angle),
                center[1] + radius * math.sin(angle),
            )

        for a, b in kids:
            pa = coords[a]  # type: ignore[assignment]
            pb = coords[b]  # type: ignore[assignment]
            mid_ab = ((pa[0] + pb[0]) / 2, (pa[1] + pb[1]) / 2)
            out_dir = _norm(mid_ab[0] - center[0], mid_ab[1] - center[1])
            place_helix(a, b, out_dir)

    # Exterior loop: place top-level elements on a circle with one open slot.
    exterior: List[int] = []
    ext_children: List[Pair] = []
    k = 0
    while k < n:
        if partner[k] == -1:
            exterior.append(k)
            k += 1
        else:
            b = partner[k]
            exterior.append(k)
            exterior.append(b)
            ext_children.append((k, b))
            k = b + 1

    m = max(len(exterior) + 1, 3)
    radius = spacing / (2 * math.sin(math.pi / m))
    center = (0.0, 0.0)
    step = 2 * math.pi / m
    start_angle = math.pi / 2 + step / 2
    for t, v in enumerate(exterior):
        angle = start_angle + t * step
        coords[v] = (
            center[0] + radius * math.cos(angle),
            center[1] + radius * math.sin(angle),
        )
    for a, b in ext_children:
        pa = coords[a]  # type: ignore[assignment]
        pb = coords[b]  # type: ignore[assignment]
        mid_ab = ((pa[0] + pb[0]) / 2, (pa[1] + pb[1]) / 2)
        out_dir = _norm(mid_ab[0] - center[0], mid_ab[1] - center[1])
        place_helix(a, b, out_dir)

    return [pt if pt is not None else (0.0, 0.0) for pt in coords]


# ---------------------------------------------------------------------------
# SVG rendering
# ---------------------------------------------------------------------------


def _fit_transform(
    points: Sequence[Point], box: Tuple[float, float, float, float], padding: float
) -> Tuple[float, float, float]:
    """Compute the (scale, dx, dy) that fits ``points`` into ``box``."""
    x0, y0, width, height = box
    xs = [p[0] for p in points]
    ys = [p[1] for p in points]
    min_x, max_x = min(xs), max(xs)
    min_y, max_y = min(ys), max(ys)
    raw_w = max(max_x - min_x, 1e-6)
    raw_h = max(max_y - min_y, 1e-6)
    scale = min((width - 2 * padding) / raw_w, (height - 2 * padding) / raw_h)
    dx = x0 + (width - raw_w * scale) / 2 - min_x * scale
    dy = y0 + (height - raw_h * scale) / 2 - min_y * scale
    return scale, dx, dy


def _align_long_axis_horizontal(points: Sequence[Point]) -> List[Point]:
    """Rotate a layout so its principal axis follows the canvas width.

    Hairpins with long leaders/reporters otherwise inherit an arbitrary
    diagonal angle from the radial layout.  Aligning their principal axis
    before fitting makes better use of a wide preview without changing any
    base order, pairing, or relative distance.
    """

    if len(points) < 2:
        return list(points)
    cx = sum(x for x, _ in points) / len(points)
    cy = sum(y for _, y in points) / len(points)
    xx = sum((x - cx) ** 2 for x, _ in points)
    yy = sum((y - cy) ** 2 for _, y in points)
    xy = sum((x - cx) * (y - cy) for x, y in points)
    if xx + yy <= 1e-12:
        return list(points)

    principal_angle = 0.5 * math.atan2(2.0 * xy, xx - yy)
    cosine = math.cos(principal_angle)
    sine = math.sin(principal_angle)
    rotated = [
        (
            (x - cx) * cosine + (y - cy) * sine,
            -(x - cx) * sine + (y - cy) * cosine,
        )
        for x, y in points
    ]
    # Keep the molecular 5' end on the left for a stable, natural reading
    # direction instead of allowing the PCA axis sign to flip between inputs.
    start_x = (points[0][0] - cx) * cosine + (points[0][1] - cy) * sine
    end_x = (points[-1][0] - cx) * cosine + (points[-1][1] - cy) * sine
    if end_x < start_x:
        rotated = [(-x, -y) for x, y in rotated]
    return rotated


def _apply_transform(points: Sequence[Point], scale: float, dx: float, dy: float) -> List[Point]:
    return [(x * scale + dx, y * scale + dy) for x, y in points]


def _fit(
    points: Sequence[Point], box: Tuple[float, float, float, float], padding: float
) -> Tuple[List[Point], float]:
    scale, dx, dy = _fit_transform(points, box, padding)
    return _apply_transform(points, scale, dx, dy), scale


def _text_color(fill: str) -> str:
    return "#ffffff" if fill in _DARK_FILLS else "#0f172a"


def _strand_ends(n: int, cuts: set) -> Tuple[set, set]:
    starts = {0}
    ends = {n - 1}
    for cut in cuts:
        ends.add(cut)
        starts.add(cut + 1)
    return starts, ends


def _consecutive_runs(indices: Sequence[int]) -> List[List[int]]:
    """Split a sorted list of indices into runs of consecutive integers."""
    runs: List[List[int]] = []
    for idx in indices:
        if runs and idx == runs[-1][-1] + 1:
            runs[-1].append(idx)
        else:
            runs.append([idx])
    return runs


def _default_pair_style(cuts: set):
    def style(i: int, j: int) -> Tuple[str, float]:
        cross = any(i <= cut < j for cut in cuts)
        return ("#94a3b8" if not cross else "#475569", 0.55 if not cross else 0.8)

    return style


def _draw_strand_body(
    svg: List[str],
    *,
    points: Sequence[Point],
    sequence: str,
    colors: Sequence[str],
    pairs: Sequence[Pair],
    cuts: set,
    radius: float,
) -> None:
    """Append base pairs, backbone and labeled bases for molecular 2D views."""
    style = _default_pair_style(cuts)

    # Base pairs (drawn first, below the backbone and bases).
    for i, j in pairs:
        x1, y1 = points[i]
        x2, y2 = points[j]
        stroke, opacity = style(i, j)
        svg.append(
            f'<line x1="{x1:.1f}" y1="{y1:.1f}" x2="{x2:.1f}" y2="{y2:.1f}" '
            f'stroke="{stroke}" stroke-width="1.4" opacity="{opacity:.2f}"/>'
        )

    # Backbone (skip strand breaks).
    n = len(points)
    for idx in range(n - 1):
        if idx in cuts:
            continue
        x1, y1 = points[idx]
        x2, y2 = points[idx + 1]
        svg.append(
            f'<line x1="{x1:.1f}" y1="{y1:.1f}" x2="{x2:.1f}" y2="{y2:.1f}" '
            f'stroke="#64748b" stroke-width="1.5" opacity="0.85"/>'
        )

    # Bases.
    for idx in range(n):
        x, y = points[idx]
        fill = colors[idx] if idx < len(colors) else "#e2e8f0"
        base = sequence[idx] if idx < len(sequence) else PLACEHOLDER_BASE
        dash = ' stroke-dasharray="2 2"' if base == PLACEHOLDER_BASE else ""
        svg.append(
            f'<circle cx="{x:.1f}" cy="{y:.1f}" r="{radius:.1f}" fill="{fill}" '
            f'stroke="#ffffff" stroke-width="1.1"{dash}/>'
        )
        svg.append(
            f'<text x="{x:.1f}" y="{y + radius * 0.35:.1f}" font-size="{radius * 1.15:.1f}" '
            f'font-weight="700" text-anchor="middle" fill="{_text_color(fill)}">{escape(base)}</text>'
        )


def render_structure_svg(
    sequence: str,
    structure: str,
    colors: Sequence[str],
    *,
    width: int = 700,
    height: int = 540,
    title: str = "",
    subtitle: str = "",
    highlight: Optional[set] = None,
    orient_horizontal: bool = False,
    base_outlines: Optional[dict[int, str]] = None,
    legend: Sequence[Tuple[str, str]] = (),
    probability_scale: Sequence[str] = (),
    footer_notes: Sequence[Tuple[str, str]] = (),
    description: str = "",
    stretch_structure_label: bool = True,
    show_annotations: bool = True,
) -> str:
    """Render one secondary structure to a standalone SVG string.

    ``highlight`` is a set of 0-based base indices; each run of consecutive
    indices is framed by a clean rounded red box (used to flag premature stop
    codons).
    Optional legend entries are ``(label, color)`` pairs. A probability scale
    supplies colors sampled uniformly from 0 to 1; footer notes are
    ``(text, color)`` pairs. These annotations share this renderer's canvas and
    typography rather than adding a second SVG around the structure.
    ``stretch_structure_label=False`` keeps short dot-bracket strings at their
    natural glyph width, while still fitting long strings to the canvas.
    Set ``show_annotations=False`` for an app view without export captions or
    dot-bracket notation. Base labels, strand ends, domain legends, probability
    scales and positional markings remain part of the molecular drawing.
    """

    if not show_annotations:
        title, subtitle, description, footer_notes = "", "", "", ()
    highlight = highlight or set()
    pairs, n, cuts = parse_structure(structure)
    if n == 0:
        return _empty_svg(width, height, title, tr('Aucune séquence à afficher'))
    if len(sequence) < n:
        sequence = sequence + PLACEHOLDER_BASE * (n - len(sequence))

    raw = radial_layout(structure)
    if orient_horizontal:
        raw = _align_long_axis_horizontal(raw)
    top = 54 if title else 24
    # Long architecture legends wrap so every domain remains legible in exports.
    legend_positions, legend_x, legend_row = [], 18., 0
    for label, color in legend:
        item_width = 34 + len(tr(label)) * 6.5
        if legend_x > 18 and legend_x + item_width > width - 18:
            legend_x, legend_row = 18., legend_row + 1
        legend_positions.append((legend_x, legend_row, label, color))
        legend_x += item_width
    legend_height = 28 * (legend_row + 1) if legend else 0
    footer_height = legend_height + (40 if probability_scale else 0) + 16 * len(footer_notes)
    drawing_height = height - footer_height
    bottom = 46 if show_annotations else 24
    points, scale = _fit(raw, (0, top, width, drawing_height - top - bottom), padding=34)
    radius = max(4.0, min(9.0, scale * 0.34))

    svg: List[str] = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" '
        f'viewBox="0 0 {width} {height}" font-family="Arial,Helvetica,sans-serif">',
        f'<rect width="{width}" height="{height}" fill="#f8fafc"/>',
    ]
    if description:
        svg.append(f'<title>{escape(tr(title))}</title><desc>{escape(tr(description))}</desc>')
    if title:
        svg.append(
            f'<text x="18" y="30" font-size="17" font-weight="700" fill="#0f172a">{escape(tr(title))}</text>'
        )
    if subtitle:
        svg.append(
            f'<text x="18" y="47" font-size="12" fill="#64748b">{escape(tr(subtitle))}</text>'
        )

    _draw_strand_body(
        svg,
        points=points,
        sequence=sequence,
        colors=colors,
        pairs=pairs,
        cuts=cuts,
        radius=radius,
    )

    for index, color in (base_outlines or {}).items():
        if 0 <= index < n:
            x, y = points[index]
            svg.append(f'<circle class="base-annotation" data-position="{index + 1}" cx="{x:.2f}" cy="{y:.2f}" r="{radius + 3:.2f}" fill="none" stroke="{escape(color)}" stroke-width="2.5"/>')

    # Highlight frame (e.g. premature stop codon): a single clean rounded box
    # around each run of consecutive flagged bases, instead of overlapping rings.
    for run in _consecutive_runs(sorted(i for i in highlight if 0 <= i < n)):
        xs = [points[i][0] for i in run]
        ys = [points[i][1] for i in run]
        x0, y0 = points[run[0]]
        x1, y1 = points[run[-1]]
        cx, cy = (x0 + x1) / 2.0, (y0 + y1) / 2.0
        pad = radius * 0.85
        angle = math.degrees(math.atan2(y1 - y0, x1 - x0)) if len(run) > 1 else 0.0
        length = math.hypot(x1 - x0, y1 - y0)
        # Extra spread of the bases off the run axis (curved regions).
        spread = 0.0
        if len(run) > 2:
            nx, ny = math.sin(math.radians(angle)), -math.cos(math.radians(angle))
            spread = max(abs((px - cx) * nx + (py - cy) * ny) for px, py in zip(xs, ys))
        half_h = radius + pad + spread
        box_w = length + 2.0 * (radius + pad)
        box_h = 2.0 * half_h
        svg.append(
            f'<rect x="{cx - box_w / 2.0:.1f}" y="{cy - box_h / 2.0:.1f}" '
            f'width="{box_w:.1f}" height="{box_h:.1f}" rx="{half_h:.1f}" ry="{half_h:.1f}" '
            f'fill="none" stroke="#dc2626" stroke-width="2.4" '
            f'transform="rotate({angle:.2f} {cx:.1f} {cy:.1f})"/>'
        )

    # 5' / 3' strand-end markers.
    starts, ends = _strand_ends(n, cuts)
    for idx in sorted(starts):
        x, y = points[idx]
        svg.append(
            f'<text x="{x:.1f}" y="{y - radius - 4:.1f}" font-size="11" font-weight="700" '
            f'text-anchor="middle" fill="#0f172a">5\u2032</text>'
        )
    for idx in sorted(ends):
        x, y = points[idx]
        svg.append(
            f'<text x="{x:.1f}" y="{y + radius + 12:.1f}" font-size="11" font-weight="700" '
            f'text-anchor="middle" fill="#0f172a">3\u2032</text>'
        )

    # Dot-bracket notation belongs in exported figures; the app offers copying
    # the original string separately from the drawing.
    if show_annotations:
        dot = structure.replace("+", " ")
        label_fit = (f'textLength="{width - 36}" lengthAdjust="spacingAndGlyphs"'
                     if stretch_structure_label or len(dot) * 6.6 > width - 36 else '')
        svg.append(
            f'<text x="18" y="{drawing_height - 16}" font-size="11" fill="#475569" '
            f'font-family="SFMono-Regular,Consolas,Menlo,monospace" '
            f'{label_fit}>{escape(dot)}</text>'
        )

    footer_y = drawing_height + 12
    if legend:
        for x, row, label, color in legend_positions:
            y = footer_y + 28 * row
            svg.append(f'<circle cx="{x + 6:.1f}" cy="{y - 4}" r="5" fill="{escape(color)}"/>')
            svg.append(f'<text x="{x + 18:.1f}" y="{y}" font-size="11" fill="#334155">{escape(tr(label))}</text>')
        footer_y += legend_height
    if probability_scale:
        bar_width = width - 36
        for index, color in enumerate(probability_scale):
            x = 18 + index * bar_width / len(probability_scale)
            svg.append(f'<rect x="{x:.2f}" y="{footer_y - 12}" width="{bar_width / len(probability_scale) + 0.1:.2f}" height="10" fill="{escape(color)}"/>')
        svg.append(f'<text x="18" y="{footer_y + 14}" font-size="11" fill="#334155">0</text>')
        svg.append(f'<text x="{width - 18}" y="{footer_y + 14}" font-size="11" text-anchor="end" fill="#334155">1</text>')
        footer_y += 40
    for note, color in footer_notes:
        svg.append(f'<text x="18" y="{footer_y}" font-size="11" fill="{escape(color)}">{escape(tr(note))}</text>')
        footer_y += 16

    svg.append("</svg>")
    return "".join(svg)


def _empty_svg(width: int, height: int, title: str, message: str) -> str:
    return (
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" '
        f'viewBox="0 0 {width} {height}" font-family="Arial,Helvetica,sans-serif">'
        f'<rect width="{width}" height="{height}" fill="#f8fafc"/>'
        f'<text x="18" y="30" font-size="17" font-weight="700" fill="#0f172a">{escape(tr(title))}</text>'
        f'<text x="{width / 2:.0f}" y="{height / 2:.0f}" font-size="13" fill="#94a3b8" '
        f'text-anchor="middle">{escape(tr(message))}</text>'
        "</svg>"
    )


def switch_preview_svg(
    architecture: SwitchArchitecture,
    trigger_dna: Optional[str] = None,
    reporter_sequence: str = "CGUAAAGGCGAGGAGCUGUUC",
    reporter_label: str = "sfGFP",
    state: str = "off",
    *,
    width: int = 700,
    height: int = 540,
    show_annotations: bool = True,
) -> str:
    """Build the model and render the OFF or ON 2D structure as SVG."""

    model = build_switch_model(
        architecture, trigger_dna, reporter_sequence, reporter_label
    )
    if state.lower() == "on":
        return render_structure_svg(
            model.on_sequence,
            model.on_structure,
            model.on_colors,
            width=width,
            height=height,
            title="Switch ON",
            subtitle=tr('Trigger + Switch ARN ({v0})', v0=reporter_label),
            orient_horizontal=True,
            show_annotations=show_annotations,
        )
    return render_structure_svg(
        model.switch_sequence,
        model.off_structure,
        model.switch_colors,
        width=width,
        height=height,
        title="Switch OFF",
        subtitle=tr('Switch ARN ({v0}) - X = base à optimiser', v0=reporter_label),
        orient_horizontal=True,
        show_annotations=show_annotations,
    )


def _reading_frame_layout(architecture: SwitchArchitecture) -> Tuple[int, int]:
    """Return exact (aug_start, reporter_start) nucleotide offsets."""
    layout = architecture_layout(architecture)
    return layout["aug_start"], layout["reporter_start"]


def premature_stop_positions(architecture: SwitchArchitecture, switch_seq: str) -> List[int]:
    """Base indices (in the switch) of premature in-frame stop codons."""
    aug_start, reporter_start = _reading_frame_layout(architecture)
    try:
        report = reading_frame_report(switch_seq, aug_start, reporter_start=reporter_start)
    except SequenceValidationError:
        return []
    positions: List[int] = []
    for stop in report["stops_found"]:
        if stop["position"] < reporter_start:
            positions.extend([stop["position"], stop["position"] + 1, stop["position"] + 2])
    return positions


def candidate_display_structure(structure: str) -> str:
    """Validate calculated notation without inventing missing structural data."""
    structure = "".join((structure or "").split())
    if not structure or any(not strand for strand in structure.split("+")):
        raise ValueError(tr('Structure calculée absente ou invalide : aucun schéma ne peut être affiché.'))
    try:
        parse_structure(structure)
    except ValueError as exc:
        raise ValueError(tr('Structure calculée invalide : aucun schéma de remplacement n’est généré.')) from exc
    return structure


def candidate_structure_svg(
    architecture: SwitchArchitecture,
    trigger_seq: str,
    switch_seq: str,
    structure: str,
    *,
    state: str = "off",
    reporter_sequence: str = "CGUAAAGGCGAGGAGCUGUUC",
    reporter_label: str = "sfGFP",
    width: int = 700,
    height: int = 520,
    show_annotations: bool = True,
) -> str:
    """Render the OFF or ON secondary structure of a finished design candidate.

    Uses the candidate's own sequences and recorded computed structure, colors the
    bases by domain, and rings any premature stop codon in red.
    """
    model = build_switch_model(
        architecture, trigger_seq, reporter_sequence, reporter_label
    )
    switch_seq = (switch_seq or "").upper()
    trigger_seq = (trigger_seq or "").upper()
    if not switch_seq or not trigger_seq or set(switch_seq) - set("ACGU") or set(trigger_seq) - set("ACGT"):
        raise ValueError(tr('Les séquences calculées du candidat sont absentes ou invalides.'))
    if len(switch_seq) != len(model.switch_sequence):
        raise ValueError(tr('La séquence calculée ne correspond pas à l’architecture du switch.'))
    stops = premature_stop_positions(architecture, switch_seq)
    structure = candidate_display_structure(structure)
    expected_lengths = [len(trigger_seq), len(switch_seq)] if state.lower() == "on" else [len(switch_seq)]
    if [len(strand) for strand in structure.split("+")] != expected_lengths:
        raise ValueError(tr('La structure calculée ne correspond pas aux longueurs des séquences du candidat.'))

    if state.lower() == "on":
        sequence = trigger_seq + switch_seq
        colors = [TRIGGER_COLOR] * len(trigger_seq) + model.switch_colors
        highlight = {len(trigger_seq) + pos for pos in stops}
        title = "Switch ON"
        subtitle = tr('Trigger + Switch ARN ({v0})', v0=reporter_label)
    else:
        sequence = switch_seq
        colors = model.switch_colors
        highlight = set(stops)
        title = "Switch OFF"
        subtitle = tr('Switch ARN ({v0})', v0=reporter_label)

    if stops:
        subtitle += tr('  -  ⚠ codon STOP prématuré (cercle rouge)')

    return render_structure_svg(
        sequence,
        structure,
        colors,
        width=width,
        height=height,
        title=title,
        subtitle=subtitle,
        highlight=highlight,
        orient_horizontal=True,
        show_annotations=show_annotations,
    )


@dataclass(frozen=True)
class LinearBase:
    """A displayed nucleotide (or an explicitly abbreviated domain token)."""

    base: str
    color: str
    position: int | str | None = None
    domain: str = ""
    tooltip: str = ""


@dataclass(frozen=True)
class LinearDomain:
    """A labeled bar covering display positions within one strand."""

    label: str
    start: int
    length: int
    color: str


@dataclass(frozen=True)
class LinearStrand:
    """A linear row in displayed order, including its 5′/3′ orientation."""

    label: str
    bases: Sequence[LinearBase]
    offset: int = 0
    reverse: bool = False
    domains: Sequence[LinearDomain] = ()
    dash_placeholder: bool = True


def render_linear_svg(
    rows: Sequence[LinearStrand],
    *,
    pairs: Sequence[Tuple[int, int, int, int]] = (),
    pair_color: str = "#b0bec5",
    number_positions: bool = False,
    title: str = "",
    subtitle: str = "",
    show_annotations: bool = True,
) -> Tuple[str, int, int]:
    """Draw nucleotide rows using the original tSwitch linear visual grammar.

    Returns ``(svg, width, height)`` at its natural, horizontally scrollable
    size. Pair entries contain ``(row_a, base_a, row_b, base_b)`` using zero-
    based displayed positions. Only explicitly supplied pairs are drawn.
    Intramolecular pairs use gray arcs above their strand; intermolecular pairs
    connect the two rows. Optional numbering uses each base's original position.
    ``LinearBase`` metadata never changes the base-box geometry.
    ``show_annotations=False`` hides export captions while retaining strand
    names, direction markers and domain labels needed to read the sequence.
    """
    if not show_annotations:
        title, subtitle = "", ""
    step, box_w, box_h = 22.0, 20.0, 26.0
    x0 = 118.0
    header = 50 if title or subtitle else 0
    arc_heights = {}
    arc_clearance = [0.0] * len(rows)
    for ra, ia, rb, ib in pairs:
        if not (0 <= ra < len(rows) and 0 <= rb < len(rows)
                and 0 <= ia < len(rows[ra].bases) and 0 <= ib < len(rows[rb].bases)):
            raise ValueError(tr('Une paire de bases du schéma linéaire fait référence à un nucléotide absent'))
        if ra == rb:
            # Logarithmic height keeps long full-reporter hairpins readable
            # without producing an SVG thousands of pixels tall.
            arc_heights[ra, ia, rb, ib] = 18.0 + 8.0 * math.log2(1 + abs(ib - ia))
            arc_clearance[ra] = max(arc_clearance[ra], arc_heights[ra, ia, rb, ib])
    row_y = []
    for index, clearance in enumerate(arc_clearance):
        y = (row_y[-1] + 90.0 if row_y else 56.0 + header)
        row_y.append(y + max(0.0, clearance - (26.0 if index else 38.0)))
    row_x = [x0 + row.offset * step for row in rows]
    width = int(x0 + max((row.offset + len(row.bases) for row in rows), default=1) * step + 48)
    height = int((row_y[-1] if row_y else 56 + header) + box_h + 78)
    if title or subtitle:
        width = max(700, width, int(36 + len(title) * 9), int(36 + len(subtitle) * 6.5))
    svg: List[str] = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" '
        f'viewBox="0 0 {width} {height}" font-family="Arial,Helvetica,sans-serif">',
        f'<rect width="{width}" height="{height}" fill="#ffffff"/>',
    ]
    if title:
        svg.append(f'<title>{escape(tr(title))}</title>')
        svg.append(f'<text x="18" y="30" font-size="17" font-weight="700" fill="#0f172a">{escape(tr(title))}</text>')
    if subtitle:
        svg.append(f'<text x="18" y="47" font-size="12" fill="#64748b">{escape(tr(subtitle))}</text>')

    for ra, ia, rb, ib in pairs:
        if ra > rb:
            ra, ia, rb, ib = rb, ib, ra, ia
        tx, sx = row_x[ra] + ia * step + box_w / 2, row_x[rb] + ib * step + box_w / 2
        if ra == rb:
            y = row_y[ra]
            top = y - arc_heights[ra, ia, rb, ib]
            svg.append(
                f'<path d="M {tx:.1f},{y:.1f} C {tx:.1f},{top:.1f} {sx:.1f},{top:.1f} {sx:.1f},{y:.1f}" '
                f'fill="none" stroke="#94a3b8" stroke-width="1" opacity="0.85"/>'
            )
            continue
        svg.append(
            f'<line x1="{tx:.1f}" y1="{row_y[ra] + box_h:.1f}" x2="{sx:.1f}" '
            f'y2="{row_y[rb]:.1f}" stroke="{escape(pair_color)}" stroke-width="1" opacity="0.85"/>'
        )

    for row, x_start, y in zip(rows, row_x, row_y):
        for idx, nucleotide in enumerate(row.bases):
            x = x_start + idx * step
            base, fill = nucleotide.base, nucleotide.color
            metadata = bool(nucleotide.domain or nucleotide.position is not None or nucleotide.tooltip)
            if metadata:
                position = "" if nucleotide.position is None else str(nucleotide.position)
                svg.append(f'<g data-domain="{escape(nucleotide.domain)}" data-position="{escape(position)}">')
                if nucleotide.tooltip:
                    svg.append(f'<title>{escape(tr(nucleotide.tooltip))}</title>')
            dash = ' stroke-dasharray="2 2"' if row.dash_placeholder and base == PLACEHOLDER_BASE else ""
            svg.append(
                f'<rect x="{x:.1f}" y="{y:.1f}" width="{box_w:.1f}" height="{box_h:.1f}" '
                f'rx="3" fill="{escape(fill)}" stroke="#ffffff"{dash}/>'
            )
            svg.append(
                f'<text x="{x + box_w / 2:.1f}" y="{y + box_h / 2 + 4:.1f}" font-size="13" '
                f'font-weight="700" text-anchor="middle" fill="{_text_color(fill)}">{escape(base)}</text>'
            )
            if metadata:
                svg.append('</g>')
            if number_positions and isinstance(nucleotide.position, int) and nucleotide.position % 10 == 0:
                svg.append(
                    f'<text x="{x + box_w / 2:.1f}" y="{y - 6:.1f}" font-size="9" '
                    f'text-anchor="middle" fill="#475569">{nucleotide.position}</text>'
                )

    for row, y in zip(rows, row_y):
        svg.append(
            f'<text x="{x0 - 12:.1f}" y="{y + box_h / 2 + 4:.1f}" font-size="12" '
            f'font-weight="700" text-anchor="end" fill="#263238">{escape(tr(row.label))}</text>'
        )
    for row, x, y in zip(rows, row_x, row_y):
        start, end = ("3′", "5′") if row.reverse else ("5′", "3′")
        svg.append(
            f'<text x="{x - 10:.1f}" y="{y - 6:.1f}" font-size="11" '
            f'text-anchor="middle" fill="#37474f">{start}</text>'
        )
        svg.append(
            f'<text x="{x + len(row.bases) * step + 2:.1f}" y="{y - 6:.1f}" '
            f'font-size="11" text-anchor="middle" fill="#37474f">{end}</text>'
        )
    for row, x, y in zip(rows, row_x, row_y):
        bar_y = y + box_h + 12
        for domain in row.domains:
            if domain.length <= 0:
                continue
            sx = x + domain.start * step
            ex = x + (domain.start + domain.length) * step - (step - box_w)
            svg.append(
                f'<line x1="{sx:.1f}" y1="{bar_y:.1f}" x2="{ex:.1f}" y2="{bar_y:.1f}" '
                f'stroke="{escape(domain.color)}" stroke-width="5" stroke-linecap="round"/>'
            )
            svg.append(
                f'<text x="{(sx + ex) / 2:.1f}" y="{bar_y + 18:.1f}" font-size="10" '
                f'text-anchor="middle" fill="#37474f">{escape(tr(domain.label))}</text>'
            )
    svg.append('</svg>')
    return "".join(svg), width, height


def linear_sequence_svg(
    architecture: SwitchArchitecture,
    trigger_dna: Optional[str] = None,
    reporter_sequence: str = "CGUAAAGGCGAGGAGCUGUUC",
    reporter_label: str = "sfGFP",
    *,
    switch_sequence: Optional[str] = None,
    show_annotations: bool = True,
) -> Tuple[str, int, int]:
    """Adapt a switch architecture to the common linear renderer.

    The trigger remains reverse-aligned with the toehold/stem, as in the
    original tSwitch diagram. A finished candidate may supply its measured
    ``switch_sequence`` to replace the live preview's optimizable ``X`` bases.
    """
    model = build_switch_model(
        architecture, trigger_dna, reporter_sequence, reporter_label
    )
    switch_seq = model.switch_sequence
    if switch_sequence is not None:
        switch_seq = "".join(switch_sequence.split()).upper()
        if len(switch_seq) != len(model.switch_sequence):
            raise ValueError(tr('La longueur de la séquence du switch ne correspond pas à l’architecture'))
        if set(switch_seq) - set("ACGU"):
            raise ValueError(tr('Un switch calculé ne doit contenir que des bases nucléotidiques déterminées'))
    aligned_start = len(architecture.t7_leader)
    domains, start = [], 0
    for domain in model.domains:
        domains.append(LinearDomain(domain.label, start, len(domain.sequence), domain.color))
        start += len(domain.sequence)
    rows = (
        LinearStrand("Trigger", [LinearBase(b, TRIGGER_COLOR) for b in model.trigger_sequence[::-1]],
                     offset=aligned_start, reverse=True, dash_placeholder=False),
        LinearStrand("Switch ARN", [LinearBase(b, c) for b, c in zip(switch_seq, model.switch_colors)], domains=domains),
    )
    pairs = [(0, i, 1, aligned_start + i) for i in range(len(model.trigger_sequence)) if aligned_start + i < len(switch_seq)]
    return render_linear_svg(rows, pairs=pairs, show_annotations=show_annotations)
