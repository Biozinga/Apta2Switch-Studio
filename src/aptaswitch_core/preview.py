"""Live 2D secondary-structure preview for toehold switches.

This module builds a switch design straight from the runtime architecture,
without NUPACK or ViennaRNA. Bases that NUPACK would optimize are shown as
``N`` placeholders. It exposes:

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

from .architecture import architecture_layout
from .models import SwitchArchitecture
from .sequences import SequenceValidationError, reading_frame_report, reverse_complement_rna

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
PLACEHOLDER_BASE = "N"

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
    switch_material: str = "rna",
    trigger_material: str = "dna",
) -> SwitchModel:
    """Build the switch design (with ``N`` placeholders) from the architecture.

    The toehold and lower stem are derived from the reverse complement of the
    trigger when one is supplied; otherwise they are shown as ``N`` too. Every
    domain that NUPACK optimizes is always rendered as ``N``. ``switch_material``
    / ``trigger_material`` control whether letters are shown as RNA (U) or DNA (T).
    """

    def _chem(seq: str, material: str) -> str:
        return seq.replace("U", "T") if material == "dna" else seq.replace("T", "U")

    arch = architecture
    toe_len = arch.toehold_length
    stem_len = arch.stem_length

    binding: Optional[str] = None
    if trigger_dna:
        try:
            binding = reverse_complement_rna(trigger_dna)
        except SequenceValidationError:
            binding = None

    if binding is not None and len(binding) >= toe_len + stem_len:
        toehold_seq = _chem(binding[:toe_len], switch_material)
        stem_l_seq = _chem(binding[toe_len : toe_len + stem_len], switch_material)
    else:
        toehold_seq = PLACEHOLDER_BASE * toe_len
        stem_l_seq = PLACEHOLDER_BASE * stem_len

    upper = arch.upper_stem_length
    upper2 = arch.upper_stem2_length
    bulge = arch.bulge_length
    # Older run objects can survive a Streamlit reload with their original
    # architecture class, which predates the optimized RBS-loop prefix.
    rbs_prefix_length = getattr(arch, "rbs_prefix_length", 0)
    aug_spacer = max(0, int(getattr(arch, "aug_spacer", 0)))
    t7_seq = _chem(arch.t7_leader.upper(), switch_material)
    rbs_seq = _chem(arch.rbs_loop.upper(), switch_material)
    linker_seq = _chem(_frame_linker(arch), switch_material)
    reporter_seq = _chem((reporter_sequence or "").upper(), switch_material)
    aug_seq = _chem("AUG", switch_material)

    n = PLACEHOLDER_BASE
    domains = [
        SwitchDomain("Leader", t7_seq, DOMAIN_COLORS["T7"], False),
        SwitchDomain("Toehold", toehold_seq, DOMAIN_COLORS["Toehold"], binding is None),
        SwitchDomain("Stem L", stem_l_seq, DOMAIN_COLORS["Stem L"], binding is None),
    ]
    if upper2:
        domains.append(SwitchDomain("U2 L", n * upper2, DOMAIN_COLORS["U2 L"], True))
    domains += [
        SwitchDomain("Bulge", n * bulge, DOMAIN_COLORS["Bulge"], True),
        SwitchDomain("U L", n * upper, DOMAIN_COLORS["U L"], True),
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

    rbs_len = len(rbs_seq)
    t7_len = len(t7_seq)
    linker_len = len(linker_seq)
    reporter_len = len(reporter_seq)

    off_structure = "".join(
        [
            "." * t7_len,
            "." * toe_len,
            "(" * stem_len,
            "(" * upper2,
            "." * bulge,
            "(" * upper,
            "." * rbs_prefix_length,
            "." * rbs_len,
            ")" * upper,
            "." * aug_spacer,
            "." * 3,
            ")" * upper2,
            ")" * stem_len,
            "." * linker_len,
            "." * reporter_len,
        ]
    )

    paired_region = toe_len + stem_len
    tail = len(switch_sequence) - t7_len - paired_region
    on_switch = "." * t7_len + ")" * paired_region + "." * tail
    on_structure = "(" * paired_region + "+" + on_switch

    trigger_sequence = _chem((trigger_dna or PLACEHOLDER_BASE * paired_region).upper(), trigger_material)

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


def _align_long_axis_horizontal(
    points: Sequence[Point], *, reference: Optional[Sequence[Point]] = None,
) -> List[Point]:
    """Rotate a layout so its principal axis follows the canvas width.

    Hairpins with long leaders/reporters otherwise inherit an arbitrary
    diagonal angle from the radial layout.  Aligning their principal axis
    before fitting makes better use of a wide preview without changing any
    base order, pairing, or relative distance. An optional reference fixes one
    common rotation for an animation instead of reorienting each frame.
    """

    reference = points if reference is None else reference
    if len(reference) < 2:
        return list(points)
    cx = sum(x for x, _ in reference) / len(reference)
    cy = sum(y for _, y in reference) / len(reference)
    xx = sum((x - cx) ** 2 for x, _ in reference)
    yy = sum((y - cy) ** 2 for _, y in reference)
    xy = sum((x - cx) * (y - cy) for x, y in reference)
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
    reference_start_x = (reference[0][0] - cx) * cosine + (reference[0][1] - cy) * sine
    reference_end_x = (reference[-1][0] - cx) * cosine + (reference[-1][1] - cy) * sine
    if reference_end_x < reference_start_x:
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
    pair_style=None,
) -> None:
    """Append base-pair lines, backbone lines and base circles/labels to ``svg``.

    Shared by the static OFF/ON/candidate renderer and the trigger-binding
    animation frames, so both use the exact same visual grammar.
    """
    style = pair_style or _default_pair_style(cuts)

    # Base pairs (drawn first, below the backbone and bases).
    for i, j in pairs:
        x1, y1 = points[i]
        x2, y2 = points[j]
        stroke, opacity = style(i, j)
        if opacity <= 0:
            continue
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
    footer_height = (28 if legend else 0) + (40 if probability_scale else 0) + 16 * len(footer_notes)
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
        for index, (label, color) in enumerate(legend):
            x = 18 + index * (width - 36) / len(legend)
            svg.append(f'<circle cx="{x + 6:.1f}" cy="{footer_y - 4}" r="5" fill="{escape(color)}"/>')
            svg.append(f'<text x="{x + 18:.1f}" y="{footer_y}" font-size="11" fill="#334155">{escape(tr(label))}</text>')
        footer_y += 28
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
    switch_material: str = "rna",
    trigger_material: str = "dna",
    width: int = 700,
    height: int = 540,
    show_annotations: bool = True,
) -> str:
    """Build the model and render the OFF or ON 2D structure as SVG."""

    model = build_switch_model(
        architecture, trigger_dna, reporter_sequence, reporter_label, switch_material, trigger_material
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
        subtitle=tr('Switch ARN ({v0}) - N = base a optimiser (NUPACK)', v0=reporter_label),
        orient_horizontal=True,
        show_annotations=show_annotations,
    )


# ---------------------------------------------------------------------------
# Trigger-binding animation ("Play animation" in the live 2D section)
# ---------------------------------------------------------------------------
#
# Rendered in the exact same radial 2D style as the static OFF/ON view (same
# colors, circles, backbone via :func:`_draw_strand_body`), so the animation
# looks like a moving version of the "Live 2D structure" panel above it.
#
# The one thing the generic radial layout gets wrong for an *animation* is the
# trigger: bases not yet paired to the switch get coiled by the algorithm into
# a loop (it treats any unpaired stretch as part of a polygon). Real toehold
# switches show the trigger as a straight incoming strand, so every frame
# below explicitly straightens the still-unbound part of the trigger into a
# line that continues the already-correct duplex ladder, instead of trusting
# the generic loop placement for it.


def _cap_structure(architecture: SwitchArchitecture) -> Tuple[str, str]:
    """Dot-bracket for the nested upper-stem/bulge/RBS cap, closed and open."""
    arch = architecture
    closed = (
        "(" * arch.upper_stem2_length
        + "." * arch.bulge_length
        + "(" * arch.upper_stem_length
        + "." * getattr(arch, "rbs_prefix_length", 0)
        + "." * len(arch.rbs_loop)
        + ")" * arch.upper_stem_length
        + "." * max(0, int(getattr(arch, "aug_spacer", 0)))
        + "." * 3
        + ")" * arch.upper_stem2_length
    )
    return closed, "." * len(closed)


def _invasion_structures(architecture: SwitchArchitecture) -> List[Tuple[str, str]]:
    """(trigger_structure, switch_structure) pairs for each displacement step.

    Step ``k`` (0 <= k <= stem_length) models the trigger having invaded the
    first ``k`` base pairs of the lower stem after landing on the toehold:
    ``k = 0`` is "just landed" (toehold bound, stem/cap still closed) and
    ``k = stem_length`` is the fully open ON state. This mirrors exactly how
    :func:`build_switch_model` already derives the OFF/ON structures, just at
    intermediate resolution.
    """
    arch = architecture
    toe_len = arch.toehold_length
    stem_len = arch.stem_length
    t7_len = len(arch.t7_leader)
    cap_closed, cap_open = _cap_structure(arch)

    steps: List[Tuple[str, str]] = []
    for k in range(stem_len + 1):
        invaded = toe_len + k
        remaining = stem_len - k
        trigger_struct = "(" * invaded + "." * remaining
        if remaining > 0:
            switch_struct = "." * t7_len + ")" * invaded + "(" * remaining + cap_closed + ")" * remaining
        else:
            switch_struct = "." * t7_len + ")" * invaded + cap_open
        steps.append((trigger_struct, switch_struct))
    return steps


def _lerp_points(a: Sequence[Point], b: Sequence[Point], t: float) -> List[Point]:
    return [(ax + (bx - ax) * t, ay + (by - ay) * t) for (ax, ay), (bx, by) in zip(a, b)]


def _straighten_dangling_trigger(raw: List[Point], bound_count: int, trigger_len: int) -> List[Point]:
    """Replace the still-unbound tail of the trigger with a straight line.

    ``radial_layout`` already places the *bound* trigger bases (indices
    ``0..bound_count-1``) as a perfectly straight duplex ladder (that is how
    helices are drawn). This continues that exact same line, at the exact
    same spacing, for the remaining ``bound_count..trigger_len-1`` bases
    instead of letting the generic algorithm coil them into a loop.
    """
    if bound_count >= trigger_len:
        return raw
    points = list(raw)
    if bound_count >= 2:
        step_x = points[bound_count - 1][0] - points[bound_count - 2][0]
        step_y = points[bound_count - 1][1] - points[bound_count - 2][1]
    else:
        # Too little bound to infer a direction yet: point away from the
        # switch centroid as a reasonable fallback.
        switch_pts = points[trigger_len:]
        cx = sum(x for x, _ in switch_pts) / max(1, len(switch_pts))
        cy = sum(y for _, y in switch_pts) / max(1, len(switch_pts))
        anchor = points[bound_count - 1] if bound_count else points[0]
        dx, dy = _norm(anchor[0] - cx, anchor[1] - cy)
        step_x, step_y = dx, dy
    ax, ay = points[max(0, bound_count - 1)]
    for i in range(bound_count, trigger_len):
        offset = i - (bound_count - 1) if bound_count else i + 1
        points[i] = (ax + step_x * offset, ay + step_y * offset)
    return points


def switch_animation_frames(
    architecture: SwitchArchitecture,
    trigger_dna: Optional[str] = None,
    reporter_sequence: str = "CGUAAAGGCGAGGAGCUGUUC",
    reporter_label: str = "sfGFP",
    *,
    switch_material: str = "rna",
    trigger_material: str = "dna",
    width: int = 760,
    height: int = 560,
    approach_steps: int = 14,
    tween_steps: int = 6,
    hold_start: int = 5,
    hold_end: int = 10,
    show_annotations: bool = True,
) -> List[str]:
    """Build an ordered list of SVG frames: the trigger (drawn as a straight
    strand throughout) slides onto the toehold, then progressively unzips the
    stem (strand displacement) until the switch is fully open (ON).

    Reuses the exact same primitives as the static OFF/ON 2D view
    (:func:`parse_structure`, :func:`radial_layout`, :func:`_draw_strand_body`)
    so the animation shares its visual grammar (colors, circles, backbone)
    with the rest of the app. Returned as a flip-book of standalone SVG
    strings meant to be played back with a timer in the UI.
    """
    model = build_switch_model(
        architecture, trigger_dna, reporter_sequence, reporter_label, switch_material, trigger_material
    )
    switch_seq = model.switch_sequence
    trigger_seq = model.trigger_sequence
    combined_sequence = trigger_seq + switch_seq
    combined_colors = model.on_colors
    trigger_len = len(trigger_seq)

    arch = architecture
    toe_len = arch.toehold_length
    stem_len = max(1, arch.stem_length)
    t7_len = len(arch.t7_leader)
    tail_len = (
        len(switch_seq)
        - t7_len
        - arch.toehold_length
        - arch.stem_length
        - len(_cap_structure(arch)[0])
    )
    tail = "." * max(0, tail_len)

    # One real, valid two-strand structure per displacement step, with the
    # still-unbound trigger tail straightened out (see docstring above).
    keyframes: List[Tuple[List[Pair], List[Point], set]] = []
    for k, (trigger_struct, switch_struct) in enumerate(_invasion_structures(arch)):
        structure = f"{trigger_struct}+{switch_struct}{tail}"
        pairs, _n, cuts = parse_structure(structure)
        raw = radial_layout(structure)
        raw = _straighten_dangling_trigger(raw, bound_count=toe_len + k, trigger_len=trigger_len)
        keyframes.append((pairs, raw, cuts))

    # Orient the final ON state as in the static view and apply that same
    # rotation to the whole clip, avoiding a rotating camera during binding.
    on_reference = keyframes[-1][1]
    keyframes = [
        (pairs, _align_long_axis_horizontal(raw, reference=on_reference), cuts)
        for pairs, raw, cuts in keyframes
    ]

    # One fixed camera for the whole clip (consistent scale for every frame),
    # THEN pin every keyframe so the toehold's outer contact point - switch
    # index ``t7_len``, where the trigger's tip always touches first - sits at
    # the exact same pixel every frame. Without this, the generic radial
    # layout can subtly re-arrange the whole molecule as the stem shortens,
    # which reads as "the switch sliding" instead of "the trigger grabbing a
    # fixed toehold and the stem opening around it".
    top = 54 if show_annotations else 24
    bottom = 46 if show_annotations else 24
    all_raw = [pt for _pairs, raw, _cuts in keyframes for pt in raw]
    scale, dx, dy = _fit_transform(all_raw, (0, top, width, height - top - bottom), padding=34)

    anchor_index = trigger_len + t7_len
    target_anchor = (
        keyframes[0][1][anchor_index][0] * scale + dx,
        keyframes[0][1][anchor_index][1] * scale + dy,
    )
    pinned: List[List[Point]] = []
    for _pairs, raw, _cuts in keyframes:
        scaled = [(x * scale, y * scale) for x, y in raw]
        ax, ay = scaled[anchor_index]
        shift = (target_anchor[0] - ax, target_anchor[1] - ay)
        pinned.append([(x + shift[0], y + shift[1]) for x, y in scaled])

    # Anchoring every keyframe to the same point can push a differently-shaped
    # keyframe (e.g. the fully-open ON state, which reorganizes into a single
    # big loop once nothing is closed anymore) outside the canvas. Re-fit the
    # union of every pinned frame as one uniform correction (same scale/shift
    # applied to all of them alike, so the anchor stays exactly as fixed as
    # before - just moved together onto a spot that fits everything).
    all_pinned = [pt for frame in pinned for pt in frame]
    scale2, dx2, dy2 = _fit_transform(all_pinned, (0, top, width, height - top - bottom), padding=34)
    frames_xy: List[List[Point]] = [_apply_transform(frame, scale2, dx2, dy2) for frame in pinned]
    radius = max(4.0, min(9.0, scale * scale2 * 0.34))

    def render(
        points: Sequence[Point],
        pairs: Sequence[Pair],
        cuts: set,
        subtitle: str,
        done: bool = False,
        pulse: Optional[Tuple[float, float, float]] = None,
    ) -> str:
        svg: List[str] = [
            f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" '
            f'viewBox="0 0 {width} {height}" font-family="Arial,Helvetica,sans-serif">',
            f'<rect width="{width}" height="{height}" fill="#f8fafc"/>',
        ]
        if show_annotations:
            svg.extend((
                f'<text x="18" y="30" font-size="17" font-weight="700" fill="#0f172a">{escape(tr("Fixation du trigger"))}</text>',
                f'<text x="18" y="47" font-size="12" fill="#64748b">{escape(tr(subtitle))}</text>',
            ))
        _draw_strand_body(
            svg,
            points=points,
            sequence=combined_sequence,
            colors=combined_colors,
            pairs=pairs,
            cuts=cuts,
            radius=radius,
        )
        if pulse is not None:
            px, py, pop = pulse
            svg.append(
                f'<circle cx="{px:.1f}" cy="{py:.1f}" r="{radius * 2.4:.1f}" fill="none" '
                f'stroke="#f59e0b" stroke-width="2" opacity="{pop:.2f}"/>'
            )
        if done and show_annotations:
            svg.append(
                f'<text x="{width - 16}" y="30" font-size="14" font-weight="700" '
                f'text-anchor="end" fill="#059669">\u2713 Switch ON</text>'
            )
        svg.append("</svg>")
        return "".join(svg)

    frames: List[str] = []

    # --- Approach: trigger (already a straight line, see keyframe 0) slides
    # in along its own axis and lands on the toehold.
    pairs0, cuts0 = keyframes[0][0], keyframes[0][2]
    landed_trigger = frames_xy[0][:trigger_len]
    switch_fixed = frames_xy[0][trigger_len:]
    if trigger_len >= 2:
        axis = _norm(
            landed_trigger[1][0] - landed_trigger[0][0],
            landed_trigger[1][1] - landed_trigger[0][1],
        )
    else:
        cx = sum(x for x, _ in switch_fixed) / len(switch_fixed)
        cy = sum(y for _, y in switch_fixed) / len(switch_fixed)
        axis = _norm(landed_trigger[0][0] - cx, landed_trigger[0][1] - cy)
    offset_distance = 130.0
    start_trigger = [(x - axis[0] * offset_distance, y - axis[1] * offset_distance) for x, y in landed_trigger]
    unbound_pairs = [(i, j) for i, j in pairs0 if not (i < trigger_len <= j)]

    for step in range(approach_steps):
        t = step / max(1, approach_steps - 1)
        eased = 1.0 - (1.0 - t) ** 3  # ease-out: fast start, gentle landing
        trigger_pts = _lerp_points(start_trigger, landed_trigger, eased)
        frames.append(
            render(
                list(trigger_pts) + list(switch_fixed),
                unbound_pairs,
                cuts0,
                tr("Le trigger s'approche du toehold..."),
            )
        )

    # --- Beat: hold on "just landed, stem still closed", with a brief pulse
    # ring at the toehold to mark the exact moment/spot it grabs on.
    anchor_px = frames_xy[0][anchor_index]
    for i in range(hold_start):
        pulse_t = i / max(1, hold_start - 1)
        pulse = (anchor_px[0], anchor_px[1], max(0.0, 1.0 - pulse_t))
        frames.append(render(frames_xy[0], pairs0, cuts0, tr("Le trigger s'accroche au toehold."), pulse=pulse))

    # --- Invasion: unzip the stem one base pair at a time, with smooth tweening.
    for k in range(stem_len):
        pairs_a, cuts_a = keyframes[k][0], keyframes[k][2]
        pairs_b, cuts_b = keyframes[k + 1][0], keyframes[k + 1][2]
        pts_a, pts_b = frames_xy[k], frames_xy[k + 1]
        for step in range(1, tween_steps + 1):
            t = step / tween_steps
            pts_t = _lerp_points(pts_a, pts_b, t)

            def style(i: int, j: int, _a=set(pairs_a), _b=set(pairs_b), _t=t, _cuts=cuts_a):
                in_a, in_b = (i, j) in _a, (i, j) in _b
                base_stroke, base_op = _default_pair_style(_cuts)(i, j)
                if in_a and in_b:
                    return base_stroke, base_op
                if in_a and not in_b:
                    return base_stroke, base_op * (1.0 - _t)
                return base_stroke, base_op * _t

            fading_pairs = list(set(pairs_a) | set(pairs_b))
            frames.append(
                render(pts_t, fading_pairs, cuts_a if t < 0.5 else cuts_b, tr("Le stem s'ouvre progressivement..."))
            )

    # --- Hold on the final, fully-open ON state.
    on_svg = render(
        frames_xy[-1], keyframes[-1][0], keyframes[-1][2], tr('Switch ON : trigger + switch appariés.'), done=True
    )
    frames.extend([on_svg] * hold_end)

    return frames


def _reading_frame_layout(architecture: SwitchArchitecture) -> Tuple[int, int]:
    """Return (aug_start, reporter_start_codon) for stop-codon scanning."""
    layout = architecture_layout(architecture)
    return layout["aug_start"], layout["reporter_start_codon"]


def premature_stop_positions(architecture: SwitchArchitecture, switch_seq: str) -> List[int]:
    """Base indices (in the switch) of premature in-frame stop codons."""
    aug_start, reporter_start_codon = _reading_frame_layout(architecture)
    try:
        report = reading_frame_report(switch_seq, aug_start, reporter_start_codon)
    except SequenceValidationError:
        return []
    positions: List[int] = []
    for stop in report["stops_found"]:
        if stop["codon_num"] < reporter_start_codon:
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
    switch_material: str = "rna",
    trigger_material: str = "dna",
    width: int = 700,
    height: int = 520,
    show_annotations: bool = True,
) -> str:
    """Render the OFF or ON secondary structure of a finished design candidate.

    Uses the candidate's own sequences and NUPACK-computed structure, colors the
    bases by domain, and rings any premature stop codon in red.
    """
    model = build_switch_model(
        architecture, trigger_seq, reporter_sequence, reporter_label, switch_material, trigger_material
    )
    switch_seq = (switch_seq or "").upper()
    trigger_seq = (trigger_seq or "").upper()
    if not switch_seq or not trigger_seq or set(switch_seq + trigger_seq) - set("ACGTU"):
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
    title: str = "",
    subtitle: str = "",
    show_annotations: bool = True,
) -> Tuple[str, int, int]:
    """Draw nucleotide rows using the original tSwitch linear visual grammar.

    Returns ``(svg, width, height)`` at its natural, horizontally scrollable
    size. Pair entries contain ``(row_a, base_a, row_b, base_b)`` using zero-
    based displayed positions. Only explicitly supplied pairs are drawn.
    ``LinearBase`` metadata never changes the base-box geometry.
    ``show_annotations=False`` hides export captions while retaining strand
    names, direction markers and domain labels needed to read the sequence.
    """
    if not show_annotations:
        title, subtitle = "", ""
    step, box_w, box_h = 22.0, 20.0, 26.0
    x0 = 118.0
    header = 50 if title or subtitle else 0
    row_y = [56.0 + header + index * 90 for index in range(len(rows))]
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
        if not (0 <= ra < len(rows) and 0 <= rb < len(rows)
                and 0 <= ia < len(rows[ra].bases) and 0 <= ib < len(rows[rb].bases)):
            raise ValueError(tr('Une paire de bases du schéma linéaire fait référence à un nucléotide absent'))
        tx, sx = row_x[ra] + ia * step + box_w / 2, row_x[rb] + ib * step + box_w / 2
        svg.append(
            f'<line x1="{tx:.1f}" y1="{row_y[ra] + box_h:.1f}" x2="{sx:.1f}" '
            f'y2="{row_y[rb]:.1f}" stroke="#b0bec5" stroke-width="1" opacity="0.85"/>'
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
    switch_material: str = "rna",
    trigger_material: str = "dna",
    *,
    switch_sequence: Optional[str] = None,
    show_annotations: bool = True,
) -> Tuple[str, int, int]:
    """Adapt a switch architecture to the common linear renderer.

    The trigger remains reverse-aligned with the toehold/stem, as in the
    original tSwitch diagram. A finished candidate may supply its measured
    ``switch_sequence`` to replace the live preview's optimizable ``N`` bases.
    """
    model = build_switch_model(
        architecture, trigger_dna, reporter_sequence, reporter_label, switch_material, trigger_material
    )
    switch_seq = model.switch_sequence
    if switch_sequence is not None:
        switch_seq = "".join(switch_sequence.split()).upper()
        if len(switch_seq) != len(model.switch_sequence):
            raise ValueError(tr('La longueur de la séquence du switch ne correspond pas à l’architecture'))
        if set(switch_seq) - set("ACGTU"):
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
