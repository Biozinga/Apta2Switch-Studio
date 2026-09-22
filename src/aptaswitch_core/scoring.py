"""Selection scoring utilities."""

from __future__ import annotations

from .models import ScoringWeights


def compute_selection_score(
    *,
    on_yield_pct: float,
    leak_pct: float,
    defect: float,
    ddg_activation: float,
    delta_g_rbs_linker: float,
    weights: ScoringWeights,
) -> float:
    """Return the weighted comparison score without erasing ranking differences.

    This score is not a percentage. Clipping it to 0–100 would hide ON and
    secondary leak differences when the weighted totals exceed that range.
    """
    raw_score = (
        (weights.on_yield * on_yield_pct)
        - (weights.leak * leak_pct)
        - (weights.defect * defect * 100.0)
        - (weights.ddg * ddg_activation)
        + (weights.rbs * delta_g_rbs_linker)
    )
    return round(raw_score, 4)


def sort_candidates(candidates):
    # Exclusion is a separate state, not a score penalty. Even a zero-score
    # evaluated candidate precedes every excluded sequence.
    return sorted(candidates, key=lambda item: (
        item.excluded,
        0 if item.excluded else (-item.score if item.score is not None else float("inf")),
        0 if item.excluded else (item.defect if item.defect is not None else float("inf")),
        item.trial,
    ))
