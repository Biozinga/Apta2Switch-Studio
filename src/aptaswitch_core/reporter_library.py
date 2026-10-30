"""User-extensible reporter library.

The app ships with five built-in reporters (sfGFP, mNeonGreen, mCherry, mRFP1,
and LacZ, see
:data:`aptaswitch_core.models.REPORTERS`). This module lets the user add as
many custom reporters as they want (name + complete coding sequence), persisted
with the internal design prefix in the local application settings file, so they
show up alongside the built-in ones in the reporter picker on every future run.
"""

from __future__ import annotations

from .localization import tr

import re
from typing import Dict

from .models import REPORTERS, ReporterContext
from .settings import load_settings, save_settings
from .sequences import normalize_rna

__all__ = [
    "all_reporters",
    "load_custom_reporters",
    "add_custom_reporter",
    "remove_custom_reporter",
    "is_builtin_reporter",
    "reporter_from_full_sequence",
]


def is_builtin_reporter(key: str) -> bool:
    return key in REPORTERS


def _slugify(label: str) -> str:
    slug = re.sub(r"[^a-z0-9]+", "-", label.strip().lower()).strip("-")
    return slug or "reporter"


def load_custom_reporters() -> Dict[str, ReporterContext]:
    """Custom reporters saved by the user, keyed by their unique slug."""
    raw = load_settings().get("reporters", [])
    reporters: Dict[str, ReporterContext] = {}
    if not isinstance(raw, list):
        return reporters
    for entry in raw:
        if not isinstance(entry, dict) or not entry.get("full_sequence_after_start_rna"):
            continue
        key = entry.get("key")
        sequence = entry.get("sequence")
        if not isinstance(key, str) or not key.strip() or not isinstance(sequence, str) or not sequence:
            continue
        label = entry.get("label") or key
        reporter = ReporterContext(
            key=key, label=label, sequence_after_start_rna=sequence,
            full_sequence_after_start_rna=entry.get("full_sequence_after_start_rna"),
            source_accession=entry.get("source_accession"), source_url=entry.get("source_url"),
        )
        try:
            reporter.validate()
        except ValueError:
            continue
        reporters[key] = reporter
    return reporters


def _save_custom_reporters(reporters: Dict[str, ReporterContext]) -> None:
    settings = load_settings()
    settings["reporters"] = [
        {"key": reporter.key, "label": reporter.label, "sequence": reporter.sequence_after_start_rna,
         "full_sequence_after_start_rna": reporter.full_sequence_after_start_rna,
         "source_accession": reporter.source_accession, "source_url": reporter.source_url}
        for reporter in reporters.values()
    ]
    save_settings(settings)


def all_reporters() -> Dict[str, ReporterContext]:
    """Built-in reporters plus every custom reporter the user has added."""
    merged: Dict[str, ReporterContext] = dict(REPORTERS)
    merged.update(load_custom_reporters())
    return merged


def reporter_from_full_sequence(
    key: str, label: str, sequence: str, *, includes_start_codon: bool = True,
) -> ReporterContext:
    """Record an explicit full CDS, removing exactly its initial start codon.

    A terminal STOP is accepted and retained. A missing terminal STOP can be
    intentional for a fusion; completeness is the user's explicit declaration,
    not an inference from length or a guessed extension of a historical prefix.
    """
    normalized = normalize_rna(sequence, field_name="reporter sequence")
    if includes_start_codon:
        if not normalized.startswith("AUG"):
            raise ValueError(tr('La séquence complète du rapporteur doit commencer par AUG ou ATG.'))
        normalized = normalized[3:]
    if len(normalized) < 30 or len(normalized) % 3:
        raise ValueError(tr('Le rapporteur complet doit contenir au moins 30 bases après AUG, dans un nombre entier de codons.'))
    if any(normalized[i:i + 3] in {"UAA", "UAG", "UGA"} for i in range(0, len(normalized) - 3, 3)):
        raise ValueError(tr('La séquence du rapporteur contient un codon STOP interne dans le cadre de lecture.'))
    reporter = ReporterContext(
        key=key, label=label, sequence_after_start_rna=normalized[:30],
        full_sequence_after_start_rna=normalized,
    )
    reporter.validate()
    return reporter


def add_custom_reporter(
    label: str, sequence: str, *, includes_start_codon: bool = True,
) -> ReporterContext:
    """Validate and persist a complete reporter and its 30-base design prefix.

    Raises :class:`aptaswitch_core.sequences.SequenceValidationError` if the
    sequence is not valid RNA/DNA, or ``ValueError`` if the name is empty.
    """
    label = label.strip()
    if not label:
        raise ValueError(tr('Le nom du rapporteur ne peut pas être vide.'))
    custom = load_custom_reporters()
    base_slug = _slugify(label)
    key = base_slug
    suffix = 2
    taken = set(REPORTERS) | set(custom)
    while key in taken:
        key = f"{base_slug}-{suffix}"
        suffix += 1

    reporter = reporter_from_full_sequence(key, label, sequence, includes_start_codon=includes_start_codon)
    custom[key] = reporter
    _save_custom_reporters(custom)
    return reporter


def remove_custom_reporter(key: str) -> None:
    if is_builtin_reporter(key):
        raise ValueError(tr('Les rapporteurs intégrés ne peuvent pas être supprimés.'))
    custom = load_custom_reporters()
    if key in custom:
        del custom[key]
        _save_custom_reporters(custom)
