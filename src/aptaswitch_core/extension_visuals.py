"""Standalone SVG figures from measured trigger-extension NUPACK results.

No folding, extrapolation or placeholder data is performed here. Structures
must contain both strands in aptamer + trigger order, and probability figures
use the reference-pair probabilities returned by the extension engine.
"""

from __future__ import annotations

from .localization import tr

import math
from collections.abc import Mapping, Sequence
from dataclasses import asdict, is_dataclass
from html import escape
from pathlib import Path
from typing import Any

from . import preview
from .preview import (APTAMER_COLOR, TRIGGER_COLOR, EXTENSION_COLOR,
                      EXTENSION_3PRIME_COLOR, LIGAND_COLOR, parse_structure)
from .extension_preview import ligand_positions
_SCORE_COLORS = ("#2166ac", "#35a573", "#fee08b", "#f46d43", "#b2182b", "#762a83")
_PROBABILITY_COLORS = ("#440154", "#3b52b4", "#00b4cf", "#36d46d", "#f2df36", "#e65022", "#a91016")


def _row(value: Any) -> dict[str, Any]:
    return dict(value) if isinstance(value, Mapping) else asdict(value) if is_dataclass(value) else vars(value)


def _number(value: Any, field: str) -> float:
    number = float(value)
    if not math.isfinite(number):
        raise ValueError(tr('Valeur non finie pour {v0}', v0=field))
    return number


def _probability(value: Any) -> float:
    number = _number(value, "probability")
    if not -1e-12 <= number <= 1 + 1e-12:
        raise ValueError(tr('Les probabilités doivent être comprises entre 0 et 1'))
    return min(1.0, max(0.0, number))


def _color(value: float, palette: Sequence[str]) -> str:
    location = max(0.0, min(1.0, value)) * (len(palette) - 1)
    index = min(int(location), len(palette) - 2)
    fraction = location - index
    start, end = palette[index], palette[index + 1]
    channels = [round(int(start[i:i + 2], 16) * (1 - fraction) + int(end[i:i + 2], 16) * fraction) for i in (1, 3, 5)]
    return "#" + "".join(f"{channel:02x}" for channel in channels)


def _text(x: float, y: float, value: Any, *, size: int = 13, anchor: str = "start", fill: str = "#334155", extra: str = "") -> str:
    return f'<text x="{x:.2f}" y="{y:.2f}" font-size="{size}" text-anchor="{anchor}" fill="{fill}" {extra}>{escape(tr(str(value)))}</text>'


def _start(width: int, height: int, title: str, description: str) -> list[str]:
    return [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}" font-family="Arial,Helvetica,sans-serif" role="img">',
        f"<title>{escape(tr(title))}</title><desc>{escape(tr(description))}</desc>",
        f'<rect width="{width}" height="{height}" fill="#ffffff"/>',
        _text(24, 31, title, size=20, fill="#0f172a", extra='font-weight="700"'),
    ]


def _line(x1: float, y1: float, x2: float, y2: float, *, color: str = "#cbd5e1", width: float = 1) -> str:
    return f'<line x1="{x1:.2f}" y1="{y1:.2f}" x2="{x2:.2f}" y2="{y2:.2f}" stroke="{color}" stroke-width="{width}"/>'


def _data_bounds(values: Sequence[float], *, lower: float | None = None,
                 upper: float | None = None) -> tuple[float, float]:
    """Fit measured extrema with a small margin, even for very small spreads.

    A fallback span is used only for identical values. There is no minimum
    span for distinct measurements and no forced zero on positive-only data.
    """
    low, high = min(values), max(values)
    spread = high - low
    padding = spread * 0.08 if spread else (abs(low) * 0.08 or 1e-6)
    # Keep a representable interval when distinct values are a few ULPs apart.
    padding = max(padding, math.ulp(low), math.ulp(high))
    low, high = low - padding, high + padding
    if lower is not None:
        low = max(lower, low)
    if upper is not None:
        high = min(upper, high)
    return low, high


def _tick_labels(values: Sequence[float], minimum_precision: int = 3) -> list[str]:
    """Use enough significant digits to distinguish the actual tick values."""
    distinct = len(set(values))
    span = max(values) - min(values)
    for precision in range(minimum_precision, 18):
        labels = [f"{value:.{precision}g}" if value else "0" for value in values]
        accurate = all(abs(float(label) - value) <= max(
            span * 0.005 if span else abs(value) * 0.00005, math.ulp(value),
        ) for label, value in zip(labels, values))
        if len(set(labels)) == distinct and accurate:
            return labels
    return labels


def _axis_labels(values: Sequence[float], minimum_precision: int = 3) -> tuple[list[str], str]:
    """Compact scientific graduations with an explicit offset when needed."""
    low, high = min(values), max(values)
    span, center = high - low, (low + high) / 2
    offset = float(f"{center:.6g}") if span and abs(center) > span * 10000 else 0.0
    residuals = [value - offset for value in values]
    magnitude = max(map(abs, residuals))
    exponent = math.floor(math.log10(magnitude)) if magnitude and (offset or magnitude < 0.001 or magnitude >= 10000) else 0
    factor = 10.0 ** exponent
    labels = _tick_labels([value / factor for value in residuals], minimum_precision)
    power = str(exponent).translate(str.maketrans("-0123456789", "⁻⁰¹²³⁴⁵⁶⁷⁸⁹"))
    note = f"× 10{power}" if exponent else ""
    if offset:
        note = tr('{v0:.6g} + graduation {v1}', v0=offset, v1=note)
    return labels, note


def render_extension_structure_svg(
    aptamer_sequence: str,
    trigger_sequence: str,
    structure: str,
    *,
    extension_length: int = 0,
    extension_3prime_length: int = 0,
    title: str = "Complexe aptamère + trigger",
    subtitle: str = "",
    nucleotide_probabilities: Sequence[float] | None = None,
    ligand_aptamer_positions: Sequence[int] = (),
    width: int = 700,
    height: int = 540,
    show_annotations: bool = True,
) -> str:
    """Draw the actual two-strand MFE structure, without filling missing bases.

    Optional probabilities give, at every nucleotide, the probability of its
    MFE pair (or its unpaired probability when unpaired in the MFE structure).
    They are ordered exactly like ``aptamer_sequence + trigger_sequence``.
    """
    aptamer = "".join(aptamer_sequence.split()).upper()
    trigger = "".join(trigger_sequence.split()).upper()
    dot = "".join(structure.split())
    if not aptamer or not trigger or set(aptamer + trigger) - set("ACGTU"):
        raise ValueError(tr('Les séquences complètes de l’aptamère et du trigger sont nécessaires'))
    if set(dot) - set(".()+") or len(dot.split("+")) != 2:
        raise ValueError(tr('Une structure dot-bracket aptamère + trigger est attendue'))
    if list(map(len, dot.split("+"))) != [len(aptamer), len(trigger)]:
        raise ValueError(tr('Les longueurs de brins de la structure ne correspondent pas aux séquences'))
    parse_structure(dot)  # Reject unbalanced input before invoking the layout.
    if not isinstance(extension_length, int) or not 0 <= extension_length <= len(trigger):
        raise ValueError(tr('Longueur totale d’extension invalide'))
    if not isinstance(extension_3prime_length, int) or not 0 <= extension_3prime_length <= extension_length:
        raise ValueError(tr('Longueur d’extension 3′ invalide'))
    prefix_length = extension_length - extension_3prime_length
    if width < 320 or height < 320:
        raise ValueError(tr('Le canevas de structure doit mesurer au moins 320 × 320'))

    colors = ([APTAMER_COLOR] * len(aptamer) + [EXTENSION_COLOR] * prefix_length
              + [TRIGGER_COLOR] * (len(trigger) - extension_length)
              + [EXTENSION_3PRIME_COLOR] * extension_3prime_length)
    probabilities = None
    if nucleotide_probabilities is not None:
        probabilities = [_probability(p) for p in nucleotide_probabilities]
        if len(probabilities) != len(aptamer) + len(trigger):
            raise ValueError(tr('Une probabilité par nucléotide est attendue'))
        colors = [_color(p, _PROBABILITY_COLORS) for p in probabilities]

    subtitle = subtitle or tr('Aptamère : {v0} nt · Trigger : {v1} nt · Ajouts : 5′ {v2} nt / 3′ {v3} nt', v0=len(aptamer), v1=len(trigger), v2=prefix_length, v3=extension_3prime_length)
    selected = ligand_positions(ligand_aptamer_positions, len(aptamer))
    if probabilities is None:
        for position in selected:
            colors[position - 1] = LIGAND_COLOR
    legend = []
    probability_scale = []
    notes = []
    if probabilities is None:
        legend = [("Aptamère", APTAMER_COLOR), ("Trigger initial", TRIGGER_COLOR)]
        if selected:
            legend.append((tr('Bases liées au ligand'), LIGAND_COLOR))
        if prefix_length:
            legend.append(("Extension 5′", EXTENSION_COLOR))
        if extension_3prime_length:
            legend.append(("Extension 3′", EXTENSION_3PRIME_COLOR))
    else:
        probability_scale = [_color(index / 99, _PROBABILITY_COLORS) for index in range(100)]
        notes.append((tr('P(paire MFE), ou P(non appariée) pour une base libre.'), "#64748b"))
    notes.append((tr('Ajouts : 5′ {v0} nt / 3′ {v1} nt · Structure NUPACK.', v0=prefix_length, v1=extension_3prime_length), "#64748b"))
    if selected:
        notes.append((tr('Cercles roses : bases de liaison au ligand (annotation utilisateur).'), LIGAND_COLOR))
    return preview.render_structure_svg(
        aptamer + trigger, dot, colors, width=width, height=height,
        title=title, subtitle=subtitle, description=tr('Structure MFE NUPACK : {v0}', v0=dot),
        base_outlines={p - 1: LIGAND_COLOR for p in selected},
        legend=legend, probability_scale=probability_scale, footer_notes=notes,
        stretch_structure_label=False,
        show_annotations=show_annotations,
    )


def extension_selection_landscape_svg(
    candidates: Sequence[Any], *, title: str = "", width: int = 1000, height: int = 640,
    fraction: float = 0.1, limit: int | None = None,
) -> str:
    """Plot the best decile of available candidates, with axes fitted to it."""
    if width < 640 or height < 400 or (limit is not None and limit < 1) or not 0 < fraction <= 1:
        raise ValueError(tr('Dimensions du graphique ou sélection des candidats invalides'))
    ranked = sorted((_row(c) for c in candidates), key=lambda row: _number(row["final_score"], "final_score"))
    count = math.ceil(len(ranked) * fraction)
    rows = ranked[:min(count, limit) if limit is not None else count]
    title = tr('{v0} — {v1} / {v2} candidats', v0=title or tr('Paysage de sélection'), v1=len(rows), v2=len(ranked))
    if fraction < 1:
        title += tr(' (meilleurs {v0:.0%})', v0=fraction)
    svg = _start(width, height, title, tr('Chaque point représente un candidat calculé. La couleur représente son score final ; un score plus bas est meilleur.'))
    if not rows:
        return "".join(svg + [_text(width / 2, height / 2, tr('Aucun candidat évalué'), anchor="middle"), "</svg>"])
    values = [(_number(row["extension_pair_probability_sum"], "extension_pair_probability_sum") + _number(row["extension_to_aptamer_probability_sum"], "extension_to_aptamer_probability_sum"), _number(row["reference_pair_probability_loss"], "reference_pair_probability_loss"), _number(row["final_score"], "final_score")) for row in rows]
    left, top, right, bottom = 115.0, 85.0, width - 150.0, height - 140.0
    xs, ys, scores = zip(*values)
    xmin, xmax = _data_bounds(xs, lower=0.0 if min(xs) >= 0 else None)
    ymin, ymax = _data_bounds(ys, lower=0.0 if min(ys) >= 0 else None)
    transform_x = lambda value: left + (value - xmin) / (xmax - xmin) * (right - left)
    transform_y = lambda value: bottom - (value - ymin) / (ymax - ymin) * (bottom - top)
    score_min, score_max = min(scores), max(scores)
    score_span = score_max - score_min or 1.0
    xticks = [xmin + step / 5 * (xmax - xmin) for step in range(6)]
    yticks = [ymin + step / 5 * (ymax - ymin) for step in range(6)]
    (xlabels, xnote), (ylabels, ynote) = _axis_labels(xticks), _axis_labels(yticks)
    if ynote:
        svg.append(_text(left, top - 16, ynote, size=11, extra='class="y-scale"'))
    if xnote:
        svg.append(_text((left + right) / 2, bottom + 42, xnote, anchor="middle", size=11, extra='class="x-scale"'))
    for step in range(6):
        fraction = step / 5
        x, y = left + (right - left) * fraction, bottom - (bottom - top) * fraction
        svg.extend((_line(x, top, x, bottom, color="#edf2f7"), _line(left, y, right, y, color="#edf2f7")))
        svg.extend((_text(x, bottom + 23, xlabels[step], anchor="middle", size=12,
                          extra=f'class="x-tick" data-value="{xticks[step]:.17g}"'),
                    _text(left - 12, y + 4, ylabels[step], anchor="end", size=12,
                          extra=f'class="y-tick" data-value="{yticks[step]:.17g}"')))
    svg.extend((_line(left, top, left, bottom, color="#475569"), _line(left, bottom, right, bottom, color="#475569")))
    for row, (x, y, score) in zip(rows, values):
        label = tr('Rang {v0} · {v1} nt · 5′ {v2} · 3′ {v3} · Score {v4:.8g} · Isolation {v5:.8g} · Perte {v6:.8g}', v0=row.get('rank', '?'), v1=row.get('trigger_length', '?'), v2=row.get('extension_5prime') or '—', v3=row.get('extension_3prime') or '—', v4=score, v5=x, v6=y)
        svg.append(f'<circle class="candidate-point" cx="{transform_x(x):.2f}" cy="{transform_y(y):.2f}" r="5.5" fill="{_color((score - score_min) / score_span, _SCORE_COLORS)}" stroke="#ffffff" stroke-width="1" data-x="{x:.15g}" data-y="{y:.15g}" data-score="{score:.15g}"><title>{escape(label)}</title></circle>')
    bx, by = transform_x(values[0][0]), transform_y(values[0][1])
    svg.append(f'<path class="best-candidate" d="M {bx:.2f} {by - 11:.2f} L {bx + 11:.2f} {by:.2f} L {bx:.2f} {by + 11:.2f} L {bx - 11:.2f} {by:.2f} Z" fill="none" stroke="#b91c1c" stroke-width="2.5"/>')
    label_x = bx + 18 if bx < (left + right) / 2 else bx - 18
    svg.append(_text(label_x, max(top + 15, by - 16), tr('Meilleur score'), fill="#b91c1c", size=12, anchor="start" if label_x > bx else "end"))
    bar_x, bar_width = right + 28, 18
    for index in range(100):
        fraction = index / 99 if score_max > score_min else 0.0
        svg.append(f'<rect x="{bar_x}" y="{bottom - (index + 1) * (bottom - top) / 100:.2f}" width="{bar_width}" height="{(bottom - top) / 100 + 0.1:.2f}" fill="{_color(fraction, _SCORE_COLORS)}"/>')
    fractions = (0, 0.25, 0.5, 0.75, 1) if score_max > score_min else (0.5,)
    score_ticks = [score_min + fraction * (score_max - score_min) for fraction in fractions]
    score_labels, score_note = _axis_labels(score_ticks, minimum_precision=4)
    for fraction, value, label in zip(fractions, score_ticks, score_labels):
        svg.append(_text(bar_x + bar_width + 8, bottom - fraction * (bottom - top) + 4,
                         label, size=11, extra=f'class="score-tick" data-value="{value:.17g}"'))
    svg.append(_text(bar_x, top - 38 if score_note else top - 16, tr('Score final'), size=11))
    if score_note:
        # A narrow color legend needs two short lines for an offset and scale.
        for index, part in enumerate(score_note.replace("graduation", "grad.").split(" + ")):
            svg.append(_text(bar_x, top - 23 + index * 13, part + (" +" if " + " in score_note and index == 0 else ""), size=10, extra='class="score-scale"'))
    svg.append(_text(bar_x, bottom + 42, tr('Plus bas = meilleur'), size=10))
    svg.extend((_text((left + right) / 2, bottom + 66, tr('Pénalité d’isolation de l’extension'), anchor="middle", size=15), _text((left + right) / 2, bottom + 87, tr('Σ P(extension appariée) + Σ P(extension ↔ aptamère)'), anchor="middle", size=12)))
    svg.append(_text(26, (top + bottom) / 2, tr('Σ perte de probabilité des paires de référence'), anchor="middle", size=14, extra=f'transform="rotate(-90 26 {(top + bottom) / 2:.2f})"'))
    svg.append(_text(width / 2, height - 31, tr('Axes linéaires ajustés aux valeurs des candidats affichés · Losange : meilleur score.'), anchor="middle", size=11, fill="#64748b"))
    zero_count = sum(y == 0 for y in ys)
    note = (tr('{v0} candidat(s) ont une perte exactement nulle ; leurs points sont alignés sur zéro.', v0=zero_count)
            if zero_count > 1 else tr('Valeurs issues de l’analyse d’ensemble NUPACK ; les points ne sont pas décalés artificiellement.'))
    svg.append(_text(width / 2, height - 14, note, anchor="middle", size=11, fill="#64748b"))
    svg.append("</svg>")
    return "".join(svg)


def _pair_label(row: Mapping[str, Any]) -> str:
    if row.get("aptamer_index") is not None and row.get("trigger_index") is not None:
        return f"A{row['aptamer_index']} – T{row['trigger_index']}"
    return f"{row.get('reference_position_a', '?')} – {row.get('reference_position_b', '?')}"


def reference_pair_probability_svg(
    reference_pairs: Sequence[Any], *, candidate_label: str = "Trigger étendu", width: int = 1000, height: int = 620
) -> str:
    """Compare every reference MFE pair probability, with candidate − reference."""
    if width < 640 or height < 360:
        raise ValueError(tr('Dimensions du graphique de probabilités invalides'))
    rows = [_row(pair) for pair in reference_pairs]
    values = [(row, _probability(row["reference_probability"]), _probability(row["candidate_probability"])) for row in rows]
    values.sort(key=lambda value: value[1], reverse=True)
    height = max(height, len(values) * 22 + 210)
    title = tr('Probabilités des {v0} paires de référence', v0=len(values))
    svg = _start(width, height, title, tr('Paires de la structure MFE originale : probabilités de référence, probabilités du candidat, et différences signées. A désigne l’aptamère et T le trigger initial.'))
    if not values:
        return "".join(svg + [_text(width / 2, height / 2, tr('La structure de référence ne contient aucune paire'), anchor="middle"), "</svg>"])
    top, bottom = 100.0, height - 118.0
    left, right = 118.0, width * 0.63
    dleft, dright = width * 0.73, width - 45.0
    all_probabilities = [p for _, a, b in values for p in (a, b)]
    pmin, pmax = _data_bounds(all_probabilities, lower=0.0, upper=1.0)
    delta_max = max(abs(b - a) for _, a, b in values)
    delta_bound = delta_max * 1.15 if delta_max else 1e-6
    px = lambda p: left + (p - pmin) / (pmax - pmin) * (right - left)
    dx = lambda delta: dleft + (delta + delta_bound) / (2 * delta_bound) * (dright - dleft)
    svg.extend((_text(left, 63, tr('Référence · {v0}', v0=tr(candidate_label)), size=12), _text((dleft + dright) / 2, 63, tr('|Δ| maximal = {v0:.4g}', v0=delta_max), anchor="middle", size=12)))
    probability_ticks = [pmin + step / 4 * (pmax - pmin) for step in range(5)]
    probability_labels, probability_note = _axis_labels(probability_ticks)
    if probability_note:
        svg.append(_text(left, top - 17, probability_note, size=10, extra='class="probability-scale"'))
    for value, label in zip(probability_ticks, probability_labels):
        x = px(value)
        svg.extend((_line(x, top - 10, x, bottom + 9, color="#edf2f7"), _text(x, bottom + 30, label, anchor="middle", size=11, extra='class="probability-tick"')))
    delta_ticks = (-delta_bound, 0, delta_bound)
    delta_labels, delta_note = _axis_labels(delta_ticks)
    if delta_note:
        svg.append(_text((dleft + dright) / 2, top - 17, delta_note, anchor="middle", size=10, extra='class="delta-scale"'))
    for delta, label in zip(delta_ticks, delta_labels):
        x = dx(delta)
        svg.extend((_line(x, top - 10, x, bottom + 9, color="#94a3b8" if delta == 0 else "#edf2f7"), _text(x, bottom + 30, label, anchor="middle", size=11, extra='class="delta-tick"')))
    for index, (row, reference, candidate) in enumerate(values):
        y = (top + bottom) / 2 if len(values) == 1 else top + index / (len(values) - 1) * (bottom - top)
        delta, label = candidate - reference, _pair_label(row)
        svg.extend((_text(left - 12, y + 4, label, size=11, anchor="end"), _line(px(reference), y, px(candidate), y, color="#cbd5e1", width=2), _line(dx(0), y, dx(delta), y, color="#f4c8c2", width=2)))
        for kind, probability, fill in (("reference", reference, "#475569"), ("candidate", candidate, "#c9392b")):
            svg.append(f'<circle class="{kind}-probability" cx="{px(probability):.2f}" cy="{y:.2f}" r="4.2" fill="{fill}" data-probability="{probability:.15g}"><title>{escape(label)} · {escape(tr(kind))} : {probability:.8g}</title></circle>')
        svg.append(f'<circle class="probability-delta" cx="{dx(delta):.2f}" cy="{y:.2f}" r="4.2" fill="#c9392b" data-delta="{delta:.15g}"><title>{escape(label)} · Δ = {delta:.8g}</title></circle>')
    svg.extend((_line(left, bottom + 9, right, bottom + 9, color="#475569"), _line(dleft, bottom + 9, dright, bottom + 9, color="#475569")))
    svg.extend((_text((left + right) / 2, bottom + 57, tr('Probabilité d’appariement d’ensemble'), anchor="middle", size=13), _text((dleft + dright) / 2, bottom + 57, tr('Δ (candidat − référence)'), anchor="middle", size=12)))
    legend_y = height - 34
    for x, fill, label in ((left, "#475569", tr('Référence originale')), (left + 200, "#c9392b", candidate_label)):
        svg.extend((f'<circle cx="{x}" cy="{legend_y - 4}" r="4.2" fill="{fill}"/>', _text(x + 12, legend_y, label, size=11)))
    svg.append(_text(width / 2, height - 12, tr('A = aptamère ; T = trigger initial. Positions locales à partir de 1 ; paires intrabrin : positions globales.'), anchor="middle", size=10, fill="#64748b"))
    svg.append("</svg>")
    return "".join(svg)


def extension_structure_comparison_svg(
    result: Mapping[str, Any], *, panel_width: int = 600, panel_height: int = 540
) -> str:
    """Place the reference and best candidate for each target side by side."""
    inputs, reference = result["input"], result["reference"]
    if not reference:
        raise ValueError(tr('Une structure de référence calculée est nécessaire'))
    aptamer, trigger = inputs["aptamer_dna"], inputs["trigger_dna"]
    best_by_length: dict[int, dict] = {}
    for item in sorted((_row(c) for c in result["candidates"]), key=lambda row: row["final_score"]):
        best_by_length.setdefault(int(item["trigger_length"]), item)
    panels = [render_extension_structure_svg(
        aptamer, trigger, reference["mfe_structure"],
        title=tr('Référence originale — {v0} nt', v0=len(trigger)),
        subtitle=f"ΔG MFE = {reference['mfe_energy_kcal_mol']:.3f} kcal/mol",
        nucleotide_probabilities=reference.get("nucleotide_probabilities"),
        ligand_aptamer_positions=inputs.get("ligand_aptamer_positions", ()),
        width=panel_width, height=panel_height,
    )]
    for length, candidate in sorted(best_by_length.items()):
        panels.append(render_extension_structure_svg(
            aptamer, candidate["extended_trigger"], candidate["mfe_structure"],
            extension_length=candidate["extension_length"],
            extension_3prime_length=len(candidate.get("extension_3prime", "")),
            title=tr('Trigger étendu — {v0} nt', v0=length),
            subtitle=f"ΔG MFE = {candidate['mfe_energy_kcal_mol']:.3f} kcal/mol · Score = {candidate['final_score']:.4g}",
            nucleotide_probabilities=candidate.get("nucleotide_probabilities"),
            ligand_aptamer_positions=inputs.get("ligand_aptamer_positions", ()),
            width=panel_width, height=panel_height,
        ))
    width, height = panel_width * len(panels), panel_height + 78
    svg = _start(width, height, tr('Structures 2D : référence et extensions retenues'), tr('Structures MFE NUPACK : complexe original puis meilleur score final pour chaque longueur cible. Chaque panneau est ajusté à sa propre échelle spatiale.'))
    svg.append(_text(24, 53, tr('Complexe original et meilleur score par longueur · Chaque structure est ajustée à son panneau.'), size=12, fill="#64748b"))
    for index, panel in enumerate(panels):
        svg.append(f'<g class="structure-panel" transform="translate({index * panel_width}, 68)">{panel}</g>')
        if index:
            svg.append(_line(index * panel_width, 75, index * panel_width, height - 18, color="#e2e8f0"))
    svg.append("</svg>")
    return "".join(svg)


def export_extension_visuals(result: Mapping[str, Any], output_dir: str | Path) -> dict[str, str]:
    """Export reference, the best structure per length, and scientific figures."""
    destination = Path(output_dir)
    destination.mkdir(parents=True, exist_ok=True)
    inputs, reference = result["input"], result["reference"]
    if not reference:
        return {}
    aptamer, trigger = inputs["aptamer_dna"], inputs["trigger_dna"]
    exports: dict[str, str] = {}

    def save(key: str, svg: str) -> None:
        path = destination / f"{key}.svg"
        path.write_text(svg, encoding="utf-8")
        exports[key] = str(path)

    save("svg_reference", render_extension_structure_svg(aptamer, trigger, reference["mfe_structure"], title=tr('Référence originale — {v0} nt', v0=len(trigger)), subtitle=f"ΔG MFE = {reference['mfe_energy_kcal_mol']:.3f} kcal/mol", nucleotide_probabilities=reference.get("nucleotide_probabilities"), ligand_aptamer_positions=inputs.get("ligand_aptamer_positions", ())))
    candidates = sorted((_row(candidate) for candidate in result["candidates"]), key=lambda row: row["final_score"])
    save("svg_landscape", extension_selection_landscape_svg(candidates))
    save("svg_comparison", extension_structure_comparison_svg(result))
    seen: set[int] = set()
    for candidate in candidates:
        length = int(candidate["trigger_length"])
        if length in seen:
            continue
        seen.add(length)
        group = [item for item in candidates if int(item["trigger_length"]) == length]
        save(f"svg_landscape_{length}", extension_selection_landscape_svg(group, title=tr('Paysage de sélection · {v0} nt', v0=length)))
        save(f"svg_candidate_{length}", render_extension_structure_svg(aptamer, candidate["extended_trigger"], candidate["mfe_structure"], extension_length=candidate["extension_length"], extension_3prime_length=len(candidate.get("extension_3prime", "")), title=tr('Trigger étendu — {v0} nt', v0=length), subtitle=f"ΔG MFE = {candidate['mfe_energy_kcal_mol']:.3f} kcal/mol · Score = {candidate['final_score']:.4g}", nucleotide_probabilities=candidate.get("nucleotide_probabilities"), ligand_aptamer_positions=inputs.get("ligand_aptamer_positions", ())))
        save(f"svg_probabilities_{length}", reference_pair_probability_svg(candidate["reference_pairs"], candidate_label=tr('Trigger étendu · {v0} nt', v0=length)))
    return exports
