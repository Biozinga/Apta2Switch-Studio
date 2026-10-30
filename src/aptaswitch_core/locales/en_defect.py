"""English text for regional OFF/ON and global defect diagnostics."""

TRANSLATIONS = {
    "Les poids des quatre zones doivent être des nombres finis positifs ou nuls.": "The four region weights must be finite nonnegative numbers.",
    "Au moins un poids de zone doit être strictement positif. Pour désactiver ce critère, mettez le poids global du défaut à zéro.": "At least one region weight must be strictly positive. To disable this criterion, set the overall defect weight to zero.",
    "Les poids des quatre zones sont invalides.": "The four region weights are invalid.",
    "Défaut OFF/ON · malus": "OFF/ON defect · penalty",
    "Defect global": "Global defect",
    "Défaut structure (global)": "Structural defect (global)",
    "Défaut concentration (global)": "Concentration defect (global)",
    "Défaut OFF": "OFF defect",
    "Défaut ON": "ON defect",
    "Défaut OFF/ON": "OFF/ON defect",
    "Les défauts sont des fractions de 0 à 1 ; plus faible est meilleur. Le défaut de structure global et le défaut de concentration global s’additionnent pour donner le Defect global. Le compromis régional est une mesure indépendante. « — » signifie non calculé, ou état désactivé lorsque ses deux poids sont nuls.": "Defects are fractions from 0 to 1; lower is better. Global structural and concentration defects add up to the global defect. The regional compromise is an independent measure. “—” means not calculated, or a disabled state when both its weights are zero.",
    "Détail des défauts OFF/ON": "OFF/ON defect details",
    "Ce résultat ne contient pas de défauts régionaux enregistrés ; leurs valeurs restent « — » dans le tableau.": "This result has no saved regional defects; their values remain “—” in the table.",
    "État": "State",
    "Zone": "Region",
    "Défaut d’ensemble (0–1)": "Ensemble defect (0–1)",
    "Nucléotides": "Nucleotides",
    "Nucléotides incorrects attendus": "Expected incorrect nucleotides",
    "Part du défaut du switch (%)": "Share of switch defect (%)",
    "Positions du switch (à partir de 1)": "Switch positions (starting at 1)",
    "La part du défaut rapporte les erreurs de la zone aux erreurs de tout le switch dans le même état, rapporteur inclus. Elle ne représente pas son poids dans le score ; elle vaut « — » lorsque le défaut total du switch est nul.": "The defect share compares errors in the region with errors in the entire switch in the same state, including the reporter. It is not the region’s weight in the score; it is “—” when total switch defect is zero.",
    "Priorités par défaut : accessibilité RBS–linker et compromis structural OFF/ON, puis association ON et ΔΔG. La fuite a un poids très faible pour départager les candidats proches, les aptamères étant dans un autre tube. Les poids restent modifiables ; 0 désactive un critère.": "Default priorities: RBS–linker accessibility and the structural OFF/ON compromise, followed by ON association and ΔΔG. Leak has a very low weight to separate close candidates, with aptamers in another tube. Weights remain adjustable; 0 disables a criterion.",
    "#### Répartition du défaut OFF/ON": "#### OFF/ON defect weights",
    "Poids relatif de cette zone. 0 retire sa contribution au compromis ; son défaut local reste calculé.": "Relative weight of this region. 0 removes its contribution to the compromise; its local defect is still calculated.",
    "Part du compromis : {part}": "Compromise share: {part}",
    "Part du compromis (%)": "Compromise share (%)",
    "Au moins un poids de zone doit être supérieur à 0. Pour désactiver tout le malus, réglez « Défaut OFF/ON · malus » à 0 en conservant des poids de zones valides.": "At least one region weight must be greater than 0. To disable the entire penalty, set “OFF/ON defect · penalty” to 0 while keeping valid region weights.",
    "Poids effectifs : OFF {off:g} % · ON {on:g} %.": "Effective weights: OFF {off:g}% · ON {on:g}%.",
    "Poids effectifs enregistrés : OFF {off:g} % · ON {on:g} %.": "Saved effective weights: OFF {off:g}% · ON {on:g}%.",
    "{state} est désactivé dans le compromis (poids nul) : son défaut agrégé vaut « — ». Les défauts locaux restent calculés.": "{state} is disabled in the compromise (zero weight): its aggregate defect is “—”. Local defects are still calculated.",
}

# Functional OFF/ON metric v2; preserve the translations for saved regional results.
TRANSLATIONS.update({
    'Le critère de défaut doit être fonctionnel, régional ou global.': 'The defect metric must be functional, regional, or global.',
    'Toehold libre': 'Unpaired toehold',
    'Tiges avant/après AUG': 'Stems before/after AUG',
    'Préfixe RBS → AUG libre': 'Unpaired RBS loop prefix → AUG',
    'Préfixe RBS → fin linker libre': 'Unpaired RBS loop prefix → end of linker',
    'Par défaut : 26,7 % pour chacun des trois premiers critères et 20 % pour RBS → fin du linker (poids 1 / 1 / 1 / 0,75). Le coefficient « Défaut OFF/ON · malus » règle séparément l’amplitude totale du malus.': 'Defaults: 26.7% for each of the first three criteria and 20% for RBS → end of linker (weights 1 / 1 / 1 / 0.75). The “OFF/ON defect · penalty” coefficient separately controls the overall penalty strength.',
    'Les deux critères ON se recouvrent : la région préfixe RBS → AUG contribue aux deux. Cela renforce son importance ; leurs parts du défaut ne doivent pas être additionnées.': 'The two ON criteria overlap: the RBS loop prefix → AUG region contributes to both. This increases its importance; their shares of the defect must not be added together.',
    'Le malus utilise quatre défauts locaux : toehold libre et tiges encadrant AUG en OFF ; préfixe RBS → fin AUG et préfixe RBS → fin linker libres en ON. Les poids sont ceux enregistrés avec ce calcul.': 'The penalty uses four local defects: an unpaired toehold and maintained stems flanking AUG in OFF; unpaired regions from the RBS loop prefix to the end of AUG and to the end of the linker in ON. Weights are those saved with this run.',
    'Détail des tiges OFF : 6 paires avant AUG et 9 après AUG dans le classique 30 nt. Les longueurs suivent l’architecture personnalisée.': 'OFF stem details: 6 pairs before AUG and 9 after AUG in the classic 30-nt architecture. Lengths follow the custom architecture.',
    'Le défaut des tiges est la moyenne sur leurs paires cibles : contributions 6/15 et 9/15 dans le classique. AUG et le bulge restent non appariés dans la cible OFF ; le maintien des tiges mesure leur séquestration indirecte.': 'Stem defect is the average over their target pairs: contributions of 6/15 and 9/15 in the classic architecture. AUG and the bulge remain unpaired in the OFF target; maintained stems measure their indirect sequestration.',
    'Tige OFF': 'OFF stem',
    'Paires cibles': 'Target pairs',
    'Avant AUG': 'Before AUG',
    'Après AUG': 'After AUG',
})

TRANSLATIONS.update({
})
