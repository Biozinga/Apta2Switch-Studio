"""Shared definitions for recorded association measurements, without rescoring."""
from aptaswitch_core.localization import tr
from .web_runtime import concentration_from_m, preferred_concentration_unit


ASSOCIATION_WITH_APTAMER_LABEL = "Trigger–switch avec aptamère (%)"
ASSOCIATION_WITH_APTAMER_WEIGHT_LABEL = "Trigger–switch avec aptamère · malus"


def association_with_aptamer_help(thermo, *, weight=None):
    def concentration(value):
        unit = preferred_concentration_unit(value)
        return f"{concentration_from_m(value, unit):g} {unit}"

    trigger = thermo.resolved_trigger_concentration_m
    switch = thermo.resolved_switch_concentration_m
    text = tr(
        "**Association trigger–switch en présence d’aptamère.** Pourcentage des switches dans un complexe contenant exactement 1 trigger et 1 switch, à l’équilibre : 100 × [complexe trigger–switch] / [switch] initial. Le dénominateur comprend tous les switches introduits, libres ou liés. Le tube contient aptamère ({trigger_value}), trigger ({trigger_value}) et switch ({switch_value}), soit un ratio aptamère:trigger:switch de {ratio}:{ratio}:1. La concentration d’aptamère est égale à celle du trigger. Les complexes concurrents jusqu’à {max_size} brins participent au calcul d’équilibre, mais seuls les dimères trigger–switch sont comptés au numérateur. Dans votre protocole, les aptamères sont dans un autre tube : c’est un scénario de comparaison, pas une simulation du transfert entre les tubes. Ce n’est pas une fuite de traduction en OFF ni un pourcentage de protéine produite ; le ligand n’est pas simulé.",
        trigger_value=concentration(trigger), switch_value=concentration(switch),
        ratio=f"{trigger / switch:g}" if switch > 0 else tr("indéfini"),
        max_size=thermo.validation_max_complex_size,
    )
    if weight == 0:
        text += " " + tr("Poids 0 : cette mesure reste informative et n’intervient pas dans le score.")
    elif weight is not None:
        text += " " + tr("Poids enregistré : {weight:g}. Le score brut retire poids × pourcentage ; réduire cette association d’un point de pourcentage apporte {weight:g} points bruts. Plus faible est meilleur pour ce seul critère.", weight=weight)
    return text
