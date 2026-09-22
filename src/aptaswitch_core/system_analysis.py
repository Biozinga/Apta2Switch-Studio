"""MFE analysis of the supplied aptamer and trigger before any extension."""

from __future__ import annotations

from .localization import tr

import math

from .models import ThermoConditions
from .nupack_engine import nupack_import_context
from .sequences import normalize_dna, normalize_rna
from .trigger_extension import parse_pairs


def analyze_system(
    import_path: str,
    aptamer_dna: str,
    trigger_dna: str,
    thermo: ThermoConditions,
) -> dict:
    """Compute one real two-strand MFE using the requested model and salts.

    Both strands use the same chemistry (``dna04`` or ``rna06``), as in the
    extension controls. Ligand annotations do not constrain the prediction:
    this function describes only the supplied nucleotide strands.
    """
    if not import_path:
        raise RuntimeError(tr('Configurez NUPACK pour calculer la structure du système initial.'))
    if thermo.material not in {"dna04", "rna06"}:
        raise ValueError(tr('Le système initial nécessite un modèle ADN dna04 ou ARN rna06.'))
    chemistry = "dna" if thermo.material == "dna04" else "rna"
    normalize = normalize_dna if chemistry == "dna" else normalize_rna
    aptamer = normalize(aptamer_dna, field_name="aptamer")
    trigger = normalize(trigger_dna, field_name="trigger")
    model_settings = {
        "material": thermo.material,
        "temperature_c": float(thermo.temperature_c),
        "sodium_m": float(thermo.sodium_m),
        "magnesium_m": float(thermo.magnesium_m),
        "ensemble": "stacking",
        "aptamer_material": chemistry,
        "trigger_material": chemistry,
    }
    if not all(math.isfinite(model_settings[key]) for key in ("temperature_c", "sodium_m", "magnesium_m")):
        raise ValueError(tr('Les conditions thermodynamiques doivent être finies.'))
    if model_settings["sodium_m"] < 0 or model_settings["magnesium_m"] < 0:
        raise ValueError(tr('Les concentrations ioniques doivent être positives ou nulles.'))

    prefix = "d" if chemistry == "dna" else "r"
    with nupack_import_context(import_path) as nupack:
        model = nupack.Model(
            material=model_settings["material"],
            celsius=model_settings["temperature_c"],
            sodium=model_settings["sodium_m"],
            magnesium=model_settings["magnesium_m"],
            ensemble=model_settings["ensemble"],
        )
        records = nupack.mfe(strands=[prefix + aptamer, prefix + trigger], model=model)
        if not records:
            raise RuntimeError(tr('NUPACK n’a retourné aucune structure MFE pour le système initial.'))
        structure = str(records[0].structure)
        energy = float(records[0].energy)

    strands = structure.split("+")
    if len(strands) != 2 or [len(part) for part in strands] != [len(aptamer), len(trigger)]:
        raise RuntimeError(tr('NUPACK a retourné une structure incompatible avec les séquences fournies.'))
    if not math.isfinite(energy):
        raise RuntimeError(tr('NUPACK a retourné une énergie MFE non finie.'))
    pairs = parse_pairs(structure)
    paired_trigger_indices = sorted(right - len(aptamer) for left, right in pairs if left < len(aptamer) <= right)
    return {
        "aptamer_sequence": aptamer,
        "trigger_sequence": trigger,
        "mfe_structure": structure,
        "mfe_energy_kcal_mol": energy,
        "paired_trigger_indices": paired_trigger_indices,
        "model": model_settings,
    }
