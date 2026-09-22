"""User-extensible reporter library.

The app ships with two built-in reporters (sfGFP, mNeonGreen, see
:data:`aptaswitch_core.models.REPORTERS`). This module lets the user add as
many custom reporters as they want (name + RNA sequence right after the start
codon), persisted in the same local settings file used for NUPACK
paths, so they show up alongside the built-in ones in the reporter picker on
every future run.
"""

from __future__ import annotations

from .localization import tr

import re
from typing import Dict

from .models import REPORTERS, ReporterContext
from .nupack_setup import load_settings, save_settings
from .sequences import normalize_rna

__all__ = [
    "all_reporters",
    "load_custom_reporters",
    "add_custom_reporter",
    "remove_custom_reporter",
    "is_builtin_reporter",
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
        if not isinstance(entry, dict):
            continue
        key = str(entry.get("key") or "").strip()
        sequence = str(entry.get("sequence") or "").strip()
        if not key or not sequence:
            continue
        label = str(entry.get("label") or key).strip()
        reporters[key] = ReporterContext(key=key, label=label, sequence_after_start_rna=sequence)
    return reporters


def _save_custom_reporters(reporters: Dict[str, ReporterContext]) -> None:
    settings = load_settings()
    settings["reporters"] = [
        {"key": reporter.key, "label": reporter.label, "sequence": reporter.sequence_after_start_rna}
        for reporter in reporters.values()
    ]
    save_settings(settings)


def all_reporters() -> Dict[str, ReporterContext]:
    """Built-in reporters plus every custom reporter the user has added."""
    merged: Dict[str, ReporterContext] = dict(REPORTERS)
    merged.update(load_custom_reporters())
    return merged


def add_custom_reporter(label: str, sequence: str) -> ReporterContext:
    """Validate, persist, and return a new custom reporter.

    Raises :class:`aptaswitch_core.sequences.SequenceValidationError` if the
    sequence is not valid RNA/DNA, or ``ValueError`` if the name is empty.
    """
    label = label.strip()
    if not label:
        raise ValueError(tr('Le nom du rapporteur ne peut pas être vide.'))
    normalized = normalize_rna(sequence, field_name="reporter sequence")

    custom = load_custom_reporters()
    base_slug = _slugify(label)
    key = base_slug
    suffix = 2
    taken = set(REPORTERS) | set(custom)
    while key in taken:
        key = f"{base_slug}-{suffix}"
        suffix += 1

    reporter = ReporterContext(key=key, label=label, sequence_after_start_rna=normalized)
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
