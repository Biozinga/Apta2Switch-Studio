"""Draw saved full-reporter MFE results using the shared molecular renderers."""
from __future__ import annotations

from collections import Counter
from dataclasses import dataclass, replace

from .localization import tr
from .models import SwitchArchitecture
from .preview import (
    APTAMER_COLOR, LinearBase, LinearDomain, LinearStrand, TRIGGER_COLOR, build_switch_model,
    candidate_display_structure, parse_structure, premature_stop_positions, render_linear_svg,
    render_structure_svg,
)

# Historical analyses remain viewable at their original size; the limit for
# launching a new calculation is deliberately independent of saved rendering.
_MAX_SAVED_COMPLEX_SIZE = 5


@dataclass(frozen=True)
class ContextDiagram:
    sequence: str
    structure: str
    colors: list[str]
    legend: list[tuple[str, str]]
    rows: list[LinearStrand]
    highlight: set[int]
    title: str

    def structure_svg(self, *, show_annotations=True):
        return render_structure_svg(
            self.sequence, self.structure, self.colors, width=1000, height=650,
            title=self.title, highlight=self.highlight, legend=self.legend,
            orient_horizontal=True, show_annotations=show_annotations,
        )

    def linear_svg(self, *, show_annotations=True):
        # Dot-bracket indices omit '+' separators. Map them to each actual
        # strand copy, never to species names or architectural target pairs.
        positions = [(row, base) for row, strand in enumerate(self.rows) for base in range(len(strand.bases))]
        pairs = [(*positions[left], *positions[right]) for left, right in parse_structure(self.structure)[0]]
        return render_linear_svg(
            self.rows, pairs=pairs, pair_color="#2563eb", number_positions=True,
            title=self.title, show_annotations=show_annotations,
        )


def _switch_row(snapshot: dict, context: dict) -> tuple[LinearStrand, list[tuple[str, str]], set[int]]:
    """Share the saved architecture annotations between every switch copy."""
    architecture = SwitchArchitecture(**snapshot["request"]["architecture"])
    model = build_switch_model(
        architecture, context["trigger_sequence"], context["reporter_sequence_after_start"],
        snapshot["request"]["reporter"]["label"],
    )
    if len(context["switch_sequence"]) != len(model.switch_sequence):
        raise ValueError(tr("La structure enregistrée ne correspond pas aux séquences du contexte complet."))

    domains, base_domains, legend = [], [], []
    start = 0
    for domain in model.domains:
        if not domain.sequence:
            continue
        domains.append(LinearDomain(domain.label, start, len(domain.sequence), domain.color))
        base_domains.extend([domain.label] * len(domain.sequence))
        start += len(domain.sequence)
        entry = (domain.label, domain.color)
        if entry not in legend:
            legend.append(entry)
    row = LinearStrand("Switch", [
        LinearBase(base, color, position=index + 1, domain=domain)
        for index, (base, color, domain) in enumerate(zip(context["switch_sequence"], model.switch_colors, base_domains))
    ], domains=domains)
    return row, legend, set(premature_stop_positions(architecture, context["switch_sequence"]))


def context_complex_diagram(snapshot: dict, context: dict, mfe: dict) -> ContextDiagram:
    """Draw an exact saved complex, including unexpected stoichiometries.

    Every strand copy follows the saved thermodynamic complex order, sequence and chemistry.
    The immutable candidate provides domain annotations only, never target pairs.
    """
    names = mfe.get("strand_order")
    if (not isinstance(names, list) or not 1 <= len(names) <= _MAX_SAVED_COMPLEX_SIZE
            or any(not isinstance(name, str) or name not in {"Switch", "Trigger", "Aptamer"} for name in names)):
        raise ValueError(tr("La structure enregistrée ne correspond pas aux séquences du contexte complet."))
    sequences = []
    for name in names:
        sequence = context.get(f"{name.lower()}_sequence")
        material = context.get(f"{name.lower()}_material")
        expected_material = "rna" if name == "Switch" else "dna"
        alphabet = "ACGU" if name == "Switch" else "ACGT"
        if material != expected_material or not isinstance(sequence, str) or not sequence or set(sequence) - set(alphabet):
            raise ValueError(tr("La structure enregistrée ne correspond pas aux séquences du contexte complet."))
        sequences.append(sequence)
    structure = candidate_display_structure(mfe.get("structure"))
    if (mfe.get("sequences") != sequences
            or [len(part) for part in structure.split("+")] != list(map(len, sequences))):
        raise ValueError(tr("La structure enregistrée ne correspond pas aux séquences du contexte complet."))

    templates, legend, switch_stops = {}, [], set()
    if "Switch" in names:
        templates["Switch"], legend, switch_stops = _switch_row(snapshot, context)
    for name, color in (("Trigger", TRIGGER_COLOR), ("Aptamer", APTAMER_COLOR)):
        if name in names:
            sequence = context[f"{name.lower()}_sequence"]
            templates[name] = LinearStrand(name, [
                LinearBase(base, color, position=index + 1, domain=name)
                for index, base in enumerate(sequence)
            ], domains=[LinearDomain(name, 0, len(sequence), color)])
            legend.append((name, color))

    totals, occurrences = Counter(names), Counter()
    rows, colors, highlight = [], [], set()
    offset = 0
    for name, sequence in zip(names, sequences):
        occurrences[name] += 1
        label = f"{name} {occurrences[name]}" if totals[name] > 1 else name
        row = replace(templates[name], label=label)
        rows.append(row)
        colors.extend(base.color for base in row.bases)
        if name == "Switch":
            highlight.update(offset + index for index in switch_stops)
        offset += len(sequence)
    return ContextDiagram(
        sequence="".join(sequences), structure=structure,
        colors=colors,
        legend=legend, rows=rows,
        highlight=highlight,
        title=tr("MFE {complex} · rapporteur complet", complex=" + ".join(row.label for row in rows)),
    )
