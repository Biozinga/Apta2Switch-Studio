"""Data adapters for the common linear sequence renderer; no predicted folds."""

from collections.abc import Sequence

from .localization import tr

from . import preview
from .preview import (
    APTAMER_COLOR, EXTENSION_COLOR, EXTENSION_3PRIME_COLOR, LIGAND_COLOR,
    TRIGGER_COLOR, LinearBase, LinearDomain, LinearStrand,
)
from .sequences import normalize_dna
from .trigger_extension import extension_split


def ligand_positions(positions: Sequence[int], aptamer_length: int) -> tuple[int, ...]:
    """Validate human-facing, one-based aptamer positions (never trigger indices)."""
    if any(type(p) is not int or not 1 <= p <= aptamer_length for p in positions):
        raise ValueError(tr('Les positions du ligand doivent appartenir à l’aptamère (numérotation à partir de 1).'))
    return tuple(sorted(set(positions)))


def _resolved_sequence(sequence: str) -> str:
    sequence = "".join(sequence.split()).upper()
    if not sequence or set(sequence) - set("ACGTU"):
        raise ValueError(tr('Une séquence nucléotidique complète est nécessaire'))
    return sequence


def _sequence_bases(sequence: str, color: str, domain: str) -> list[LinearBase]:
    return [LinearBase(base, color, i, domain, f"{domain} · {i}") for i, base in enumerate(sequence, 1)]


def _sequence_rows(aptamer: str, core: str, prefix: Sequence[LinearBase], suffix: Sequence[LinearBase],
                   selected_positions: Sequence[int]) -> tuple[LinearStrand, LinearStrand]:
    selected = set(ligand_positions(selected_positions, len(aptamer)))
    aptamer_bases = [LinearBase(base, LIGAND_COLOR if i in selected else APTAMER_COLOR, i,
                                "ligand" if i in selected else "aptamère", tr('Aptamère · {v0}', v0=i) + (tr(' · liaison au ligand') if i in selected else ""))
                     for i, base in enumerate(aptamer, 1)]
    domains = []
    if prefix:
        domains.append(LinearDomain("5′", 0, len(prefix), EXTENSION_COLOR))
    domains.append(LinearDomain("Trigger initial", len(prefix), len(core), TRIGGER_COLOR))
    if suffix:
        domains.append(LinearDomain("3′", len(prefix) + len(core), len(suffix), EXTENSION_3PRIME_COLOR))
    return (
        LinearStrand("Aptamère", aptamer_bases, domains=[LinearDomain("Aptamère", 0, len(aptamer), APTAMER_COLOR)]),
        LinearStrand("Trigger", list(prefix) + _sequence_bases(core, TRIGGER_COLOR, "trigger-initial") + list(suffix), domains=domains),
    )


def aptamer_trigger_linear_svg(
    aptamer_sequence: str,
    trigger_sequence: str,
    *,
    extension_length: int = 0,
    extension_3prime_length: int = 0,
    ligand_aptamer_positions: Sequence[int] = (),
    title: str = "",
    subtitle: str = "",
    show_annotations: bool = True,
) -> tuple[str, int, int]:
    """Draw actual aptamer/trigger sequences with the same boxes as tSwitch.

    Both strands are in their given 5′→3′ order. Extension lengths describe
    the existing prefix/suffix in ``trigger_sequence``; no bases or base pairs
    are inferred. Position annotations on the initial trigger stay local to
    that initial sequence, independent of the chosen extension direction.
    """
    aptamer = _resolved_sequence(aptamer_sequence)
    trigger = _resolved_sequence(trigger_sequence)
    if type(extension_length) is not int or not 0 <= extension_length < len(trigger):
        raise ValueError(tr('Longueur totale d’extension invalide'))
    if type(extension_3prime_length) is not int or not 0 <= extension_3prime_length <= extension_length:
        raise ValueError(tr('Longueur d’extension 3′ invalide'))
    prefix_length = extension_length - extension_3prime_length
    core_end = len(trigger) - extension_3prime_length
    prefix = _sequence_bases(trigger[:prefix_length], EXTENSION_COLOR, "extension-5prime")
    suffix = _sequence_bases(trigger[core_end:], EXTENSION_3PRIME_COLOR, "extension-3prime")
    rows = _sequence_rows(aptamer, trigger[prefix_length:core_end], prefix, suffix, ligand_aptamer_positions)
    return preview.render_linear_svg(rows, title=title, subtitle=subtitle, show_annotations=show_annotations)


def extension_preview_svg(aptamer: str, trigger: str, target_length: int, side: str,
                          selected_positions: Sequence[int] = (), *, rna: bool = False,
                          show_annotations: bool = True) -> tuple[str, int, int]:
    """Adapt fixed strands and unknown extensions to the common linear view.

    Long unknown domains are explicitly abbreviated without allocating their
    whole requested sequence. The requested counts remain in the figure;
    actual input nucleotides and ligand positions are never abbreviated.
    """
    aptamer = normalize_dna(aptamer, field_name="aptamère")
    trigger = normalize_dna(trigger, field_name="trigger")
    prefix_length, suffix_length = extension_split(target_length - len(trigger), side)
    if rna:
        aptamer, trigger = aptamer.replace("T", "U"), trigger.replace("T", "U")

    def unknown(count: int, domain: str, color: str) -> list[LinearBase]:
        if count <= 28:
            return [LinearBase("N", color, i, domain, tr('{v0} · {v1} · à optimiser', v0=domain, v1=i)) for i in range(1, count + 1)]
        return ([LinearBase("N", color, i, domain) for i in range(1, 4)]
                + [LinearBase("…", color, f"{count} nt", domain, tr('Domaine de {v0} bases à optimiser · affichage abrégé', v0=count))]
                + [LinearBase("N", color, i, domain) for i in range(count - 2, count + 1)])

    rows = _sequence_rows(aptamer, trigger,
                          unknown(prefix_length, "extension-5prime", EXTENSION_COLOR),
                          unknown(suffix_length, "extension-3prime", EXTENSION_3PRIME_COLOR), selected_positions)
    return preview.render_linear_svg(
        rows,
        title=tr('Extension · trigger final de {v0} nt', v0=target_length),
        subtitle=tr('5′ : +{v0} nt · 3′ : +{v1} nt · N : base à optimiser avec NUPACK · Rose : liaison au ligand', v0=prefix_length, v1=suffix_length),
        show_annotations=show_annotations,
    )
