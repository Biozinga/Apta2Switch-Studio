"""Sequence utilities for DNA/RNA inputs."""

from __future__ import annotations

from .localization import tr

DNA_ALPHABET = frozenset("ACGT")
RNA_ALPHABET = frozenset("ACGU")
DNA_COMPLEMENT = str.maketrans("ACGT", "TGCA")
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


def to_rna(sequence: str) -> str:
    return normalize_dna(sequence).replace("T", "U")


def to_dna(sequence: str) -> str:
    return normalize_dna(sequence)


def reverse_complement_dna(sequence: str) -> str:
    return normalize_dna(sequence).translate(DNA_COMPLEMENT)[::-1]


def reverse_complement_rna(sequence: str) -> str:
    return normalize_rna(sequence).translate(RNA_COMPLEMENT)[::-1]


def as_nupack_dna(sequence: str) -> str:
    return f"d{normalize_dna(sequence)}"


def as_nupack_rna(sequence: str) -> str:
    return f"r{normalize_rna(sequence)}"


def reading_frame_report(switch_rna: str, aug_start: int, reporter_start_codon: int) -> dict:
    """Scan the coding region for premature stop codons."""
    switch = normalize_rna(switch_rna, field_name="switch")
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
        "has_premature_stop": first_stop is not None and first_stop < reporter_start_codon,
        "first_stop": first_stop,
    }
