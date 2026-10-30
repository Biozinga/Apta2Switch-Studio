"""Engine-independent validation and reporting of aptamer/trigger populations.

Complex orders are grouped only after equilibrium has been solved. Fractions
of complex particles and fractions of incorporated strands remain distinct.
"""
from __future__ import annotations

from dataclasses import dataclass
import math
from numbers import Integral

from ..localization import tr


AVOGADRO_CONSTANT = 6.02214076e23
MASS_BALANCE_RELATIVE_TOLERANCE = 1e-5
# Only roundoff smaller than one trillionth of the most concentrated input
# species may be clipped to zero. All larger negative values are errors.
NEGATIVE_ROUNDOFF_RELATIVE_TOLERANCE = 1e-12


@dataclass(frozen=True)
class TubeAnalysisConditions:
    aptamer_concentration_m: float = 5e-6
    trigger_concentration_m: float = 5e-6
    volume_l: float = 20e-6
    max_complex_size: int = 3

    def __post_init__(self) -> None:
        for name in ("aptamer_concentration_m", "trigger_concentration_m", "volume_l"):
            value = getattr(self, name)
            try:
                numeric = float(value)
            except (TypeError, ValueError, OverflowError):
                numeric = math.nan
            if isinstance(value, bool) or not math.isfinite(numeric) or numeric <= 0:
                raise ValueError(tr('Les concentrations des brins et le volume du tube doivent être finis et strictement positifs.'))
            object.__setattr__(self, name, numeric)
        if (isinstance(self.max_complex_size, bool)
                or not isinstance(self.max_complex_size, Integral)
                or not 1 <= self.max_complex_size <= 3):
            raise ValueError(tr('La taille maximale des complexes doit être un entier entre 1 et 3.'))
        object.__setattr__(self, "max_complex_size", int(self.max_complex_size))



def _population(complex_concentrations, conditions: TubeAnalysisConditions, *, mixed: bool) -> dict:
    aptamer_total = conditions.aptamer_concentration_m if mixed else 0.0
    trigger_total = conditions.trigger_concentration_m
    negative_tolerance = max(aptamer_total, trigger_total) * NEGATIVE_ROUNDOFF_RELATIVE_TOLERANCE
    grouped: dict[tuple[int, int], list[float]] = {}
    for complex_obj, value in complex_concentrations.items():
        try:
            concentration = float(value)
        except (TypeError, ValueError, OverflowError) as error:
            raise RuntimeError(tr('Le moteur thermodynamique a retourné une concentration de complexe invalide.')) from error
        if not math.isfinite(concentration) or concentration < -negative_tolerance:
            raise RuntimeError(tr('Le moteur thermodynamique a retourné une concentration de complexe invalide.'))
        concentration = max(0.0, concentration)
        names = [strand.name for strand in complex_obj.strands]
        if (not names or len(names) > conditions.max_complex_size
                or any(name not in {"A", "T"} for name in names)
                or (not mixed and "A" in names)):
            raise RuntimeError(tr('Le moteur thermodynamique a retourné un complexe incompatible avec le tube demandé.'))
        key = (names.count("A"), names.count("T"))
        grouped.setdefault(key, []).append(concentration)

    concentrations = {key: math.fsum(values) for key, values in grouped.items()}
    recovered_aptamer = math.fsum(n_a * concentration for (n_a, _), concentration in concentrations.items())
    recovered_trigger = math.fsum(n_t * concentration for (_, n_t), concentration in concentrations.items())
    for recovered, expected in ((recovered_aptamer, aptamer_total), (recovered_trigger, trigger_total)):
        if not math.isclose(recovered, expected, rel_tol=MASS_BALANCE_RELATIVE_TOLERANCE, abs_tol=0.0):
            raise RuntimeError(tr('Les concentrations calculées par le moteur thermodynamique ne conservent pas les quantités de brins du tube.'))

    total_complex_concentration = math.fsum(concentrations.values())
    molecule_factor = conditions.volume_l * AVOGADRO_CONSTANT
    rows = []
    for (n_a, n_t), concentration in concentrations.items():
        strand_count = n_a + n_t
        category = "free" if strand_count == 1 else "target" if (n_a, n_t) == (1, 1) else "parasite"
        label = " + ".join(
            name if count == 1 else f"{count} {name}"
            for name, count in (("A", n_a), ("T", n_t)) if count
        )
        count = concentration * molecule_factor
        rows.append({
            "aptamer_count": n_a,
            "trigger_count": n_t,
            "strand_count": strand_count,
            "label": label,
            "category": category,
            "concentration_m": concentration,
            "complex_molecule_count": count,
            "trigger_molecule_count": n_t * count,
            "aptamer_molecule_count": n_a * count,
            "trigger_fraction_pct": 100.0 * (n_t * concentration / trigger_total),
            "aptamer_fraction_pct": 100.0 * (n_a * concentration / aptamer_total) if aptamer_total else 0.0,
            "particle_fraction_pct": 100.0 * (concentration / total_complex_concentration),
        })
    rows.sort(key=lambda row: (-row["concentration_m"], row["strand_count"], row["label"]))
    summary = {
        f"{species}_{category}_pct": math.fsum(
            row[f"{species}_fraction_pct"] for row in rows if row["category"] == category
        )
        for species in ("trigger", "aptamer")
        for category in ("free", "target", "parasite")
    }
    summary.update({
        "parasite_particle_fraction_pct": math.fsum(row["particle_fraction_pct"] for row in rows if row["category"] == "parasite"),
        "parasite_complex_molecule_count": math.fsum(row["complex_molecule_count"] for row in rows if row["category"] == "parasite"),
        "total_complex_molecule_count": math.fsum(row["complex_molecule_count"] for row in rows),
    })
    if not all(math.isfinite(value) for row in rows for value in row.values() if isinstance(value, (int, float))):
        raise RuntimeError(tr('Les populations calculées dépassent les limites numériques.'))
    return {
        "total_strand_concentrations_m": {"aptamer": aptamer_total, "trigger": trigger_total},
        "complexes": rows,
        "summary": summary,
    }
