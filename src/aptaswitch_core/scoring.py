"""Selection scoring utilities."""

from __future__ import annotations

from math import atan, isfinite, pi

from .models import ScoringWeights


def selection_score_100(raw_score: float) -> float:
    """Map a finite comparison score monotonically onto the 0–100 scale.

    The fixed scale preserves comparisons across runs and avoids clipping
    distinct good or poor candidates to the same threshold. This is a score,
    not a percentage yield or an experimentally calibrated success rate.
    """
    if isinstance(raw_score, bool) or not isfinite(raw_score):
        raise ValueError("A selection score must be a finite number.")
    return round(50.0 + (100.0 / pi) * atan(raw_score / 100.0), 6)


def candidate_score_100(candidate, score_scale: str = "bounded_v1") -> float | None:
    """Validate the recorded bounded display score."""
    if score_scale != "bounded_v1":
        raise ValueError("Unknown selection score scale.")
    if candidate.excluded or candidate.score is None:
        return None
    score = candidate.score
    if isinstance(score, bool) or not isfinite(score) or not 0 <= score <= 100:
        raise ValueError("The bounded selection score must lie between zero and 100.")
    return round(float(score), 6)


def compute_selection_score(
    *,
    on_yield_pct: float,
    leak_pct: float,
    ddg_activation: float,
    delta_g_rbs_linker: float,
    weights: ScoringWeights,
    functional_defect: float,
    trigger_linearity_defect: float | None = None,
    rbs_calculator_merit: float | None = None,
) -> float:
    """Calculate the functional comparison sum retained as ``score_raw``.

    The engine applies ``selection_score_100`` for the public bounded score.
    Only inverse-extension callers provide the isolated-trigger defect.
    """
    if (isinstance(functional_defect, bool) or functional_defect is None
            or not isfinite(functional_defect) or not 0 <= functional_defect <= 1):
        raise ValueError("Functional scoring requires a finite OFF/ON defect between 0 and 1.")
    trigger_penalty = 0.0
    rbs_bonus = 0.0
    if rbs_calculator_merit is not None:
        if (isinstance(rbs_calculator_merit, bool) or not isfinite(rbs_calculator_merit)
                or not 0 <= rbs_calculator_merit <= 1):
            raise ValueError("RBS comparison merit must be finite and between zero and one.")
        if isinstance(weights.rbs_calculator, bool) or not isfinite(weights.rbs_calculator) or weights.rbs_calculator < 0:
            raise ValueError("RBS comparison weight must be finite and nonnegative.")
        rbs_bonus = weights.rbs_calculator * 100.0 * rbs_calculator_merit
    if trigger_linearity_defect is not None:
        if (isinstance(trigger_linearity_defect, bool)
                or not isfinite(trigger_linearity_defect)
                or not 0 <= trigger_linearity_defect <= 1):
            raise ValueError("Trigger linearity scoring requires a finite defect between zero and one.")
        weight = weights.trigger_linearity
        if isinstance(weight, bool) or not isfinite(weight) or weight < 0:
            raise ValueError("The trigger linearity weight must be finite and nonnegative.")
        trigger_penalty = weight * 100.0 * trigger_linearity_defect
    raw_score = (
        (weights.on_yield * on_yield_pct)
        - (weights.leak * leak_pct)
        - (weights.defect * functional_defect * 100.0)
        - (weights.ddg * ddg_activation)
        + (weights.rbs * delta_g_rbs_linker)
        - trigger_penalty
        + rbs_bonus
    )
    if not isfinite(raw_score):
        raise ValueError("A selection score must be a finite number.")
    raw_score = round(raw_score, 4)
    return raw_score


def sort_candidates(candidates, *, defect_weight: float | None = None):
    def comparison_score(item):
        # A bounded display can round nearby raw values to the same number.
        # Retain their actual order before applying secondary criteria.
        raw = getattr(item, "score_raw", None)
        score = item.score if raw is None else raw
        return -score if score is not None else float("inf")

    def tie_defect(item):
        if defect_weight == 0:
            return 0
        value = item.regional_defects.get("balanced")
        return value if value is not None else float("inf")

    # Exclusion is a separate state, not a score penalty. Even a zero-score
    # evaluated candidate precedes every excluded sequence.
    return sorted(candidates, key=lambda item: (
        item.excluded,
        0 if item.excluded else comparison_score(item),
        0 if item.excluded else tie_defect(item),
        item.trial,
    ))
