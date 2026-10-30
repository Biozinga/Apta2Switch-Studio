"""Geometry and configuration checks for sequence-guided oxNA preparation."""
from __future__ import annotations

import json
import math
from pathlib import Path
import re

import numpy as np


class OxDNAPreparationCancelled(RuntimeError):
    """The user cancelled preparation before a configuration was accepted."""


def _number(value, label: str, *, positive=False) -> float:
    if isinstance(value, bool):
        raise ValueError(f"{label} must be finite")
    try:
        result = float(value)
    except (ValueError, TypeError) as exc:
        raise ValueError(f"{label} must be finite") from exc
    if not math.isfinite(result) or (positive and result <= 0):
        raise ValueError(f"{label} must be finite" + (" and positive" if positive else ""))
    return result


def _write_report(path: Path, report: dict) -> None:
    temporary = path.with_suffix(".json.tmp")
    temporary.write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")
    temporary.replace(path)


def _read_configuration(path: Path, count: int, *, box_size: float) -> tuple[np.ndarray, np.ndarray]:
    """Read one native checkpoint without changing positions or orientations."""
    return _parse_configuration(path.read_text(encoding="utf-8"), count, box_size=box_size)


def _parse_configuration(frame: str, count: int, *, box_size=None):
    lines = [line.strip() for line in frame.splitlines() if line.strip()]
    if len(lines) != count + 3:
        raise ValueError("oxDNA configuration has an unexpected number of nucleotides or frames")
    header = []
    for name, line, size in zip(("t", "b", "E"), lines[:3], (1, 3, 3)):
        match = re.fullmatch(rf"{name}\s*=\s*(.*)", line)
        if not match:
            raise ValueError("Malformed oxDNA configuration header")
        values = np.asarray([float(value) for value in match[1].split()])
        if len(values) != size or not np.isfinite(values).all():
            raise ValueError("Non-finite or malformed oxDNA configuration header")
        header.append(values)
    box = header[1]
    if np.any(box <= 0) or (box_size is not None and not np.allclose(box, box_size, rtol=0, atol=1e-5)):
        raise ValueError("The oxDNA simulation box changed unexpectedly")
    rows = np.asarray([[float(value) for value in line.split()] for line in lines[3:]], dtype=float)
    if rows.shape != (count, 15) or not np.isfinite(rows).all():
        raise ValueError("oxDNA configuration requires 15 finite values per nucleotide")
    a1, a3 = rows[:, 3:6], rows[:, 6:9]
    if (not np.allclose(np.linalg.norm(a1, axis=1), 1, atol=1e-3, rtol=0)
            or not np.allclose(np.linalg.norm(a3, axis=1), 1, atol=1e-3, rtol=0)
            or not np.allclose(np.sum(a1 * a3, axis=1), 0, atol=1e-3, rtol=0)):
        raise ValueError("oxDNA nucleotide orientation vectors are not orthonormal")
    return rows, box


def configuration_geometry(frame: str, strand_lengths, chemistries) -> dict:
    """Screen a prepared oxNA model, not production FFS trajectories.

    Constants are the bundled oxDNA commit 75ae96f's oxDNA2/oxRNA2/DRH
    defaults (model.h, rna_model.h, drh_model.h). R* is the start of the
    strongly repulsive LJ core, not a hard-sphere radius. Ordinary soft
    contacts outside R* are allowed. This is a local geometry check, not
    proof of an equilibrated ensemble or absence of chain entanglement.

    Covalent chains must already be unwrapped, as required by native NA.
    For nonbonded pairs, image the COM displacement before adding sites.
    """
    lengths, acids = tuple(strand_lengths), tuple(chemistries)
    if (not lengths or len(lengths) != len(acids)
            or any(type(n) is not int or n < 1 for n in lengths)
            or any(acid not in ('rna', 'dna') for acid in acids)):
        raise ValueError('Geometry checks require nonempty RNA/DNA strands.')
    rows, box = _parse_configuration(frame, sum(lengths))
    a1, a3 = rows[:, 3:6], rows[:, 6:9]
    offsets = np.empty((len(rows), 3))
    bonded = np.zeros((len(rows), len(rows)), dtype=bool)
    broken, stretches = [], []
    start = 0
    for length, acid in zip(lengths, acids):
        end = start + length
        offsets[start:end] = (-.4*a1[start:end] + .2*a3[start:end] if acid == 'rna'
                             else -.34*a1[start:end] + .3408*np.cross(a3[start:end], a1[start:end]))
        idx = np.arange(start, end-1)
        bonded[idx, idx+1] = bonded[idx+1, idx] = True
        backs = rows[start:end, :3] + offsets[start:end]
        distance = np.linalg.norm(np.diff(backs, axis=0), axis=1)
        rest = .761070781051 if acid == 'rna' else .7564
        stretches.extend(np.abs(distance-rest).tolist())
        broken.extend([int(start+k+1), int(start+k+2)] for k in np.flatnonzero(np.abs(distance-rest) >= .25))
        start = end
    raw_delta = rows[None, :, :3] - rows[:, None, :3]
    delta = raw_delta - box*np.floor(raw_delta/box + .5)
    delta[bonded] = raw_delta[bonded]
    upper = np.triu(np.ones(bonded.shape, dtype=bool), 1)
    bases = .4*a1
    contacts, minima, count = [], {}, 0
    for label, left, right, radius, mask in (
            ('backbone-backbone', offsets, offsets, .675, upper & ~bonded),
            ('base-base', bases, bases, .32, upper),
            ('backbone-base', offsets, bases, .50, upper),
            ('base-backbone', bases, offsets, .50, upper)):
        distance = np.linalg.norm(delta + right[None, :, :] - left[:, None, :], axis=2)
        minima[label] = float(distance[mask].min()) if mask.any() else None
        indices = np.argwhere(mask & (distance < radius))
        count += len(indices)
        # Keep diagnostics bounded even for a badly overlapping seed.
        for i, j in indices[:max(0, 12-len(contacts))]:
            contacts.append({'nucleotides': [int(i+1), int(j+1)], 'sites': label,
                             'distance': float(distance[i, j]), 'core_boundary': radius})
    return {'passed': not broken and count == 0, 'model': 'NA (oxDNA2/oxRNA2)',
            'length_unit': 'oxDNA', 'minimum_site_distances': minima,
            'repulsive_core_contacts': count, 'contact_examples': contacts,
            'invalid_backbone_bonds': broken, 'maximum_backbone_extension': max(stretches, default=0.),
            'equilibrium_validated': False, 'topology_validated': False}


def _place_strands(first: np.ndarray, second: np.ndarray, separation: float, box: float,
                   *, periodic_margin: float = 4., auto_size: bool = False,
                   reorient: bool = False, reset_velocities: bool = True):
    """Place monomers with a periodic gap; automatic resizing is preparation only."""
    first, second = first.copy(), second.copy()
    for rows in (first, second):
        if reorient and len(rows) > 1:
            # Proper rigid rotation: put the long principal axis along Z and
            # the short axis along the separation direction X. Never reflect
            # a chiral nucleotide or deform its intramolecular geometry.
            centered = rows[:, :3] - rows[:, :3].mean(axis=0)
            _, rotation = np.linalg.eigh(centered.T @ centered)
            if np.linalg.det(rotation) < 0:
                rotation[:, 0] *= -1
            for start in range(0, 15, 3):
                rows[:, start:start+3] = rows[:, start:start+3] @ rotation
        rows[:, :3] -= (rows[:, :3].min(axis=0) + rows[:, :3].max(axis=0)) / 2
    second[:, 0] += first[:, 0].max() - second[:, 0].min() + separation
    joined = np.vstack([first, second])
    span = np.ptp(joined[:, :3], axis=0)
    minimum_box = max(span[0] + max(separation, periodic_margin),
                      span[1] + periodic_margin, span[2] + periodic_margin)
    if auto_size:
        box = max(box, float(math.ceil(minimum_box)))
    if box < minimum_box - 1e-8:
        raise ValueError(f"Box size {box:g} is too small for the separated input structures; minimum {minimum_box:g}")
    joined[:, :3] += box / 2 - (joined[:, :3].min(axis=0) + joined[:, :3].max(axis=0)) / 2
    if reset_velocities:
        joined[:, 9:] = 0  # Preparation stages explicitly refresh velocities.
    return joined, box
