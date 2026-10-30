"""Sequence utilities for DNA/RNA inputs."""

from __future__ import annotations

from .localization import tr

DNA_ALPHABET = frozenset("ACGT")
RNA_ALPHABET = frozenset("ACGU")
OPTIMIZABLE_BASE = "X"
RNA_COMPLEMENT = str.maketrans("ACGU", "UGCA")
STOP_CODONS = {"UAA", "UAG", "UGA"}


class SequenceValidationError(ValueError):
    """Raised when a sequence contains unsupported bases."""


def strip_sequence(sequence: str) -> str:
    """Remove whitespace from user sequence input."""
    return "".join(str(sequence or "").split()).upper()


def normalize_dna(sequence: str, *, field_name: str = "sequence") -> str:
    """Normalize DNA/RNA user input to DNA."""
    cleaned = strip_sequence(sequence).replace("U", "T")
    if not cleaned:
        raise SequenceValidationError(tr('{v0} est requis.', v0=tr(field_name)))
    invalid = sorted(set(cleaned) - DNA_ALPHABET)
    if invalid:
        joined = ", ".join(invalid)
        raise SequenceValidationError(tr('{v0} contient des bases non prises en charge : {v1}.', v0=tr(field_name), v1=joined))
    return cleaned


def normalize_rna(sequence: str, *, field_name: str = "sequence") -> str:
    """Normalize DNA/RNA user input to RNA."""
    cleaned = strip_sequence(sequence).replace("T", "U")
    if not cleaned:
        raise SequenceValidationError(tr('{v0} est requis.', v0=tr(field_name)))
    invalid = sorted(set(cleaned) - RNA_ALPHABET)
    if invalid:
        joined = ", ".join(invalid)
        raise SequenceValidationError(tr('{v0} contient des bases non prises en charge : {v1}.', v0=tr(field_name), v1=joined))
    return cleaned


def reverse_complement_rna(sequence: str) -> str:
    return normalize_rna(sequence).translate(RNA_COMPLEMENT)[::-1]


def reverse_complement_rna_pattern(sequence: str) -> str:
    """Complement an internal design pattern; ordinary user inputs stay strict."""
    pattern = strip_sequence(sequence).replace("T", "U")
    if not pattern or set(pattern) - (RNA_ALPHABET | {OPTIMIZABLE_BASE}):
        raise SequenceValidationError(tr("Le motif de conception ne peut contenir que A, C, G, T, U ou X."))
    return pattern.translate(RNA_COMPLEMENT)[::-1]


def reading_frame_report(switch_rna: str, aug_start: int, *, reporter_start: int) -> dict:
    """Scan in the AUG frame using the reporter's exact nucleotide boundary.

    A STOP whose first base precedes the reporter is premature, including a
    codon crossing that boundary in an independently supplied sequence.
    """
    switch = normalize_rna(switch_rna, field_name="switch")
    if (type(aug_start) is not int or type(reporter_start) is not int
            or not 0 <= aug_start < reporter_start <= len(switch)):
        raise ValueError("The AUG and reporter require ordered nucleotide offsets within the switch")
    coding_region = switch[aug_start:]
    stops = []
    for idx in range(0, len(coding_region) - 2, 3):
        codon = coding_region[idx : idx + 3]
        if codon in STOP_CODONS:
            codon_number = (idx // 3) + 1
            stops.append(
                {
                    "codon_num": codon_number,
                    "codon": codon,
                    "position": aug_start + idx,
                }
            )
    first_stop = stops[0]["codon_num"] if stops else None
    return {
        "aug_position": aug_start,
        "stops_found": stops,
        "has_premature_stop": bool(stops) and stops[0]["position"] < reporter_start,
        "first_stop": first_stop,
    }
