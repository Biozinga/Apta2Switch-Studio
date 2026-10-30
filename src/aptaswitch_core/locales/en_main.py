"""English translations for the main interface."""

TRANSLATIONS = {
    'OSTIR ne trouve aucun site d’initiation calculable autour de l’AUG choisi ({state}). Aucun score nul n’est inventé. Vérifiez le RBS personnalisé ou mettez son poids RBS ON/OFF à 0 pour concevoir sans ce critère.': 'OSTIR found no evaluable initiation site around the selected AUG ({state}). No zero score is fabricated. Check the custom RBS or set its RBS ON/OFF weight to 0 to design without this criterion.',

    'Initiation de traduction · OSTIR': 'Translation initiation · OSTIR',
    'Comparaison exploratoire : ARN seul OFF / ARN contraint ON. Les valeurs ne sont pas des mesures de production de protéine.': 'Exploratory comparison: unconstrained OFF RNA / constrained ON RNA. These values are not measurements of protein production.',
    'OFF · référence ARN': 'OFF · RNA reference',
    'ON · ARN contraint': 'ON · constrained RNA',
    'Rapport ON/OFF exploratoire': 'Exploratory ON/OFF ratio',
    'Indice RBS (/100)': 'RBS index (/100)',
    'Plus élevé est meilleur. L’indice combine l’activité ON et le contraste ON/OFF ; un ratio élevé entre deux états presque inactifs reste peu favorisé. Ce n’est pas un pourcentage de fonctionnement.': 'Higher is better. The index combines ON activity and ON/OFF contrast; a high ratio between two nearly inactive states receives little reward. This is not a percentage of function.',
    'Hypothèses et paramètres OSTIR': 'OSTIR assumptions and parameters',
    'L’ARN du candidat est conservé. En ON, seules les bases censées être occupées par le trigger sont empêchées de former d’autres paires. Le RBS reste disponible pour le ribosome. Le calcul suppose le trigger fixé : il ne calcule ni l’énergie de liaison ADN–ARN, ni la concentration du complexe.': 'The candidate RNA is preserved. In ON, only bases assumed to be occupied by the trigger are prevented from forming other pairs. The RBS remains available to the ribosome. The calculation assumes a bound trigger: it computes neither DNA–RNA binding energy nor complex concentration.',
    'Référence distincte des conditions du run : E. coli · 37 °C · Na⁺ 1,021 M · OSTIR 1.1.3 / ViennaRNA 2.7.2. Taux en unités arbitraires OSTIR.': 'Reference separate from run conditions: E. coli · 37 °C · Na⁺ 1.021 M · OSTIR 1.1.3 / ViennaRNA 2.7.2. Rates in OSTIR arbitrary units.',
    'Fenêtre du modèle : bases {start} à {end} du switch, autour de l’AUG en position {aug}. Les interactions avec les bases hors de cette fenêtre ne sont pas évaluées.': 'Model window: switch bases {start} to {end}, around AUG at position {aug}. Interactions with bases outside this window are not evaluated.',
    'K est un repère de normalisation du modèle, pas un objectif expérimental. Le score brut reçoit +100 × poids RBS ON/OFF × cet indice entre 0 et 1.': 'K is a model normalization pivot, not an experimental target. The raw score receives +100 × RBS ON/OFF weight × this index between 0 and 1.',
    'Indicateur OSTIR : + 100 × wRBS ON/OFF × indice exploratoire.': 'OSTIR indicator: +100 × wRBS ON/OFF × exploratory index.',
    'log10 ON/OFF · OSTIR': 'log10 ON/OFF · OSTIR',
    'L’affinage combine le défaut global, la stabilité RBS–linker, l’indice OSTIR et, en extension inverse, la linéarité du trigger, selon leurs poids. Les limites des défauts locaux restent prioritaires. Ce coût n’est pas une probabilité.': 'Refinement combines global defect, RBS–linker stability, the OSTIR index and, for inverse extension, trigger linearity, according to their weights. Local-defect limits retain priority. This cost is not a probability.',
    'Affinage du score': 'Score refinement',

    'Les bases fixées sont incompatibles avec les contraintes d’appariement.': 'Fixed bases are incompatible with the pairing constraints.',
    'Les bases fixées et les contraintes d’appariement imposent un STOP prématuré au codon {codon}.': 'Fixed bases and pairing constraints force a premature STOP at codon {codon}.',
    'Les bases fixées et les contraintes d’appariement ne permettent aucune séquence sans STOP prématuré avant le rapporteur.': 'Fixed bases and pairing constraints permit no sequence without a premature STOP before the reporter.',
    'Le poids RBS–linker doit être un nombre fini positif ou nul.': 'The RBS–linker weight must be a finite nonnegative number.',
    'Ce checkpoint utilise une autre méthode de conception. Ses candidats restent consultables ; lancez un nouveau design pour utiliser la méthode actuelle.': 'This checkpoint uses a different design method. Its candidates remain available; start a new design to use the current method.',
    ' · coût d’affinage {cost:.4f}': ' · refinement cost {cost:.4f}',
    'Coût d’affinage': 'Refinement cost',
    'L’affinage combine le défaut global et la stabilité RBS–linker selon son poids. Ce coût n’est pas une probabilité ; la cible globale de recherche porte sur ce coût combiné.': 'Refinement combines global defect and RBS–linker stability using its weight. This cost is not a probability; the global search target applies to this combined cost.',
    'La génération et les mutations évitent les STOP avant le rapporteur. Si les bases fixées imposent un STOP, la conception s’arrête avant les analyses thermodynamiques.': 'Generation and mutations avoid STOP codons before the reporter. If fixed bases force a STOP, design stops before thermodynamic analyses.',
    'Le linker personnalisé place le rapporteur hors du cadre de lecture : la longueur après AUG avant le rapporteur doit être un multiple de 3.': 'The custom linker puts the reporter out of frame: the length after AUG and before the reporter must be a multiple of 3.',
    'Défaut global et stabilité RBS–linker': 'Global defect and RBS–linker stability',
    'Ouverture simultanée du RBS → AUG en ON': 'Simultaneous RBS → AUG opening in ON',
    'Région entièrement libre (%)': 'Entire region unpaired (%)',
    'Coût d’ouverture (kcal/mol)': 'Opening cost (kcal/mol)',
    'Bases {start} à {end} du switch, du préfixe RBS à la fin d’AUG.': 'Switch bases {start} to {end}, from the RBS prefix through AUG.',
    'Probabilité que toutes ces bases soient libres simultanément, dans le complexe ADN–ARN déjà formé. Ce n’est ni le produit des probabilités par base, ni le rendement d’association, ni une prédiction d’expression. Mesure informative, sans malus supplémentaire dans le score.': 'Probability that all these bases are simultaneously unpaired in the already formed DNA–RNA complex. This is not a product of per-base probabilities, association yield, or expression prediction. Informative measurement, with no additional score penalty.',
    'Aucune structure du modèle ne satisfait cette ouverture dans le complexe formé ; son coût est infini.': 'No structure in the model satisfies this opening in the formed complex; the opening cost is infinite.',
    'Probabilité trop petite pour être représentée numériquement ; le coût d’ouverture reste calculé à partir de son logarithme.': 'Probability is too small to represent numerically; the opening cost is still calculated from its logarithm.',
    'RBS → AUG entièrement libre ON (%)': 'Entire RBS → AUG unpaired ON (%)',
    'ΔG ouverture RBS → AUG ON': 'RBS → AUG opening ΔG ON',
    'ARN isolé : algorithmes ViennaRNA 2.7.2. Complexes ADN–ARN : moteur hybride expérimental Apta2Switch. Les paramètres ARN et la correction saline sont communs ; ce ne sont pas les réglages par défaut de RNAfold.':
        'Isolated RNA: ViennaRNA 2.7.2 algorithms. DNA–RNA complexes: Apta2Switch experimental hybrid engine. RNA parameters and salt correction are shared; these are not the default RNAfold settings.',
    'Aperçu abrégé du journal. Le fichier complet est conservé dans les exports.':
        'Abbreviated log preview. The complete file is preserved in the exports.',
    'Console du moteur': 'Engine log',
    '{run} · {engine} · {ranked} classés · {excluded} exclus':
        '{run} · {engine} · {ranked} ranked · {excluded} excluded',
    'Calcul arrêté : les candidats terminés restent consultables.':
        'Calculation stopped: completed candidates remain available.',
    'Vue des résultats': 'Results view',
    'Classement': 'Ranking',
    'Comparaison de l’extension': 'Extension comparison',
    'Comprendre le classement': 'Understanding the ranking',
    'Points à vérifier ({count})': 'Points to review ({count})',
    'Tableau complet des candidats': 'Complete candidate table',
    'Explorer ce candidat': 'Explore this candidate',
    'Agrandir le tableau complet': 'Expand the complete table',
    'Rang {rank} · essai {trial} · score {score}': 'Rank {rank} · trial {trial} · score {score}',
    'Optimiser le linker (XXX…)': 'Optimize linker (XXX…)',
    'Remplace toutes les bases du linker par X sans changer sa longueur. Sa séquence sera choisie pendant la conception.':
        'Replaces all linker bases with X without changing its length. Its sequence will be chosen during design.',
    'Rétablir les valeurs par défaut': 'Restore defaults',
    'Rétablit le leader GGG, le RBS, le linker et l’architecture standard adaptée au trigger.':
        'Restores the GGG leader, RBS, linker and standard architecture for the trigger.',
    'Cliquez dans la séquence pour remplacer le linker par défaut. A, C, G, U ou T imposent une base ; X laisse une base à optimiser. Une séquence vide supprime le linker.':
        'Click the sequence to replace the default linker. A, C, G, U or T fixes a base; X leaves a base to optimize. An empty sequence removes the linker.',
    'Linker : {length} nt · après les {stem_length} nt de tige qui suivent AUG, avant le rapporteur.':
        'Linker: {length} nt · after the {stem_length} stem nucleotides following AUG, before the reporter.',
    'indéfini':
        'undefined',
    'Navigation':
        'Navigation',
    'Calcul en cours · {current}/{total}':
        'Calculation running · {current}/{total}',
    'Le toehold activateur standard nécessite un trigger d’au moins {trigger_region_length} nt. Activez l’extension du trigger.':
        'The standard toehold activator requires a trigger of at least {trigger_region_length} nt. Enable trigger extension.',
    'La région du trigger utilisée pour le switch doit contenir au moins 2 bases.':
        'The trigger region used for the switch must contain at least 2 bases.',
    'Renseignez les séquences disponibles, les conditions du système et, si elles sont connues, les bases impliquées dans la liaison au ligand.':
        'Enter the available sequences, system conditions and, if known, the bases involved in ligand binding.',
    '#### Séquences du système':
        '#### System sequences',
    'Molécule cible':
        'Target molecule',
    "Séquence de l'aptamère":
        'Aptamer sequence',
    'Séquence ADN — espaces et retours à la ligne acceptés':
        'DNA sequence — spaces and line breaks accepted',
    'Séquence du trigger':
        'Trigger sequence',
    'Aptamère : {value} nt · Trigger : {value1} nt':
        'Aptamer: {value} nt · Trigger: {value1} nt',
    'Attention : votre trigger mesure {value} nt. Sa longueur peut être trop faible ; vous pouvez activer l’extension du trigger.':
        'Warning: your trigger is {value} nt long. It may be too short; you can enable trigger extension.',
    'Bases de l’aptamère impliquées dans la liaison au ligand (si connues)':
        'Aptamer bases involved in ligand binding (if known)',
    'Leader 5′':
        '5′ leader',
    'RBS':
        'RBS',
    'Chimie des brins':
        'Strand chemistry',
    'Linker automatique':
        'Automatic linker',
    'Linker personnalisé':
        'Custom linker',
    'Le modèle standard nécessite un trigger d’au moins {trigger_region_length} nt. Activez l’extension du trigger.':
        'The standard model requires a trigger of at least {trigger_region_length} nt. Enable trigger extension.',
    'État du switch':
        'Switch state',
    'Conditions':
        'Conditions',
    'Scoring':
        'Scoring',
    'Sortie':
        'Output',
    'Essais':
        'Trials',
    'Température (°C)':
        'Temperature (°C)',
    'Reporter':
        'Reporter',
    '##### Concentrations':
        '##### Concentrations',
    'Trigger':
        'Trigger',
    'Toehold switch':
        'Toehold switch',
    'Sodium':
        'Sodium',
    'Magnésium':
        'Magnesium',
    'Le rendement ON est calculé par rapport à la concentration initiale du toehold switch.':
        'ON yield is calculated relative to the initial toehold switch concentration.',
    'Voir la formule':
        'Show the formula',
    'Exclure les candidats avec un STOP prématuré':
        'Exclude candidates with a premature STOP',
    'Exclure les candidats sans linéarité parfaite du RBS–linker en ON':
        'Exclude candidates without a fully unpaired RBS–linker region in ON',
    'Criblage de la structure ON après le filtre STOP, avant les autres analyses. Les candidats exclus restent en fin de tableau, sans score ni calcul de rendement ou de fuite.':
        'The ON structure is screened after the STOP filter, before the other analyses. Excluded candidates remain at the bottom of the table, without a score, yield calculation or leak calculation.',
    'Dossier local des runs':
        'Local runs folder',
    'Exports : JSON, CSV, XLSX, SVG, journal et projet .aptaswitch.json.':
        'Exports: JSON, CSV, XLSX, SVG, log and .aptaswitch.json project.',
    'Je confirme ce calcul long':
        'I confirm this long calculation',
    'Lancer la conception':
        'Start design',
    'Progression':
        'Progress',
    '{current} / {total}':
        '{current} / {total}',
    'Temps restant estimé':
        'Estimated time remaining',
    'Le temps restant est recalculé après chaque candidat terminé.':
        'The remaining time is recalculated after each completed candidate.',
    'Annuler le calcul':
        'Cancel calculation',
    'Rang':
        'Rank',
    'Score':
        'Score',
    'Defect':
        'Defect',
    'ON %':
        'ON %',
    'Leak %':
        'Leak %',
    'ΔG OFF':
        'ΔG OFF',
    'ΔG ON':
        'ΔG ON',
    'ΔΔG':
        'ΔΔG',
    'ΔG RBS–linker':
        'ΔG RBS–linker',
    'STOP':
        'STOP',
    'Premier STOP':
        'First STOP',
    'Essai':
        'Trial',
    'Type de résultats':
        'Result type',
    'Ouvrez un run enregistré ci-dessus ou lancez un calcul depuis Conception.':
        'Open a saved run above or start a calculation from Design.',
    'Run':
        'Run',
    'Candidats':
        'Candidates',
    'Tous les candidats sont exclus par le criblage. Aucun candidat n’a été retenu pour le classement.':
        'All candidates were excluded by screening. No candidate was retained for ranking.',
    'Copie locale : {output_dir}':
        'Local copy: {output_dir}',
    "Le calcul n'a produit aucun candidat.":
        'The calculation did not produce any candidates.',
    'Candidat à explorer':
        'Candidate to explore',
    'Candidat exclu · {exclusion_reason}':
        'Excluded candidate · {exclusion_reason}',
    'Les calculs de structure, rendement et fuite ont été ignorés pour ce candidat. Sa séquence reste consultable et exportée.':
        'Structure, yield and leak calculations were skipped for this candidate. Its sequence remains available to view and export.',
    '**Trigger 5′ → 3′**':
        '**Trigger 5′ → 3′**',
    '**Switch 5′ → 3′**':
        '**Switch 5′ → 3′**',
    'Structure ON conservée lors du criblage. Les autres analyses ont été ignorées ; la séquence et cette structure restent consultables et exportables.':
        'ON structure retained during screening. The other analyses were skipped; the sequence and this structure remain available to view and export.',
    'Structure':
        'Structure',
    'Structure calculée indisponible : {exc}':
        'Calculated structure unavailable: {exc}',
    'Aperçu indisponible : {exc}':
        'Preview unavailable: {exc}',
    'Séquence linéaire indisponible : {exc}':
        'Linear sequence unavailable: {exc}',
    'Séquences du candidat':
        'Candidate sequences',
    '← Candidat précédent':
        '← Previous candidate',
    'Afficher le candidat précédent dans le classement.':
        'Show the previous candidate in the ranking.',
    'Candidat suivant →':
        'Next candidate →',
    'Afficher le candidat suivant dans le classement.':
        'Show the next candidate in the ranking.',
    'Moteurs externes':
        'External engines',
    "Ces outils ne sont pas distribués avec Apta2Switch-Studio ; leurs licences propres s'appliquent.":
        'These tools are not distributed with Apta2Switch-Studio; their own licenses apply.',
    'Confirmez la responsabilité de licence avant de continuer.':
        'Confirm responsibility for the license before continuing.',
    'Bibliothèque de reporters':
        'Reporter library',
    '← Retour à la conception':
        '← Back to design',
    '**{label}**':
        '**{label}**',
    'Intégré':
        'Built-in',
    'Supprimer':
        'Delete',
    '#### Ajouter un reporter':
        '#### Add a reporter',
    'Nom':
        'Name',
    'ex. mCherry':
        'e.g. mCherry',
    'Ajouter':
        'Add',
    'Reporter {label} ajouté.':
        'Reporter {label} added.',
    'Préférences locales':
        'Local preferences',
    'Enregistrer':
        'Save',
    'Le dossier ne peut pas être vide.':
        'The folder cannot be empty.',
    '#### Réinitialisation':
        '#### Reset',
    'Efface la configuration des moteurs, les reporters personnalisés et le dossier par défaut. Les exports ne sont pas supprimés.':
        'Clears engine configuration, custom reporters and the default folder. Exports are not deleted.',
    'Je confirme la réinitialisation de tous les réglages':
        'I confirm resetting all settings',
    'Réinitialiser les réglages':
        'Reset settings',
    'Conception':
        'Design',
    'Résultats':
        'Results',
    'Réglages':
        'Settings',
    'Langue':
        'Language',
    'Reporters':
        'Reporters',
    'Général':
        'General',
    'Ajouter un reporter…':
        'Add a reporter…',
    'Estimation en cours…':
        'Estimating…',
    "Aucune région libre d'au moins 2 nt n'a été trouvée.":
        'No unbound region of at least 2 nt was found.',
    'Un calcul est déjà en cours. Suivez-le ou annulez-le dans Résultats.':
        'A calculation is already running. Follow or cancel it in Results.',
    "Impossible d'écrire les exports : {exc}":
        'Unable to write exports: {exc}',
    "Le dossier {requested_output} n'était pas inscriptible. Les exports seront enregistrés dans {writable_output}.":
        'The folder {requested_output} was not writable. Exports will be saved in {writable_output}.',
    'ON % · bonus':
        'ON % · bonus',
    'ΔΔG · liaison':
        'ΔΔG · binding',
    'Fuite · départage':
        'Leak · tie-breaker',
    'En attente des premiers messages…':
        'Waiting for the first messages…',
    'Évalué — STOP prématuré':
        'Evaluated — premature STOP',
    'Évalué':
        'Evaluated',
    'critère de criblage':
        'screening criterion',
    'STOP prématuré':
        'premature STOP',
    'Exclu — {reason}':
        'Excluded — {reason}',
    'RBS–linker non linéaire en ON':
        'RBS–linker not fully unpaired in ON',
    'oui':
        'yes',
    'non':
        'no',
    'Annulé':
        'Cancelled',
    'Terminé':
        'Completed',
    'Tous les exports (.zip)':
        'All exports (.zip)',
    'Exporter les résultats':
        'Export results',
    '{status} · essai {trial}':
        '{status} · trial {trial}',
    'Un composant de l’application est manquant. Téléchargez à nouveau l’application complète depuis la page des versions.':
        'An application component is missing. Download the complete application again from the releases page.',
    'Impossible d’ouvrir le navigateur. Ouvrez {url}.':
        'Unable to open the browser. Open {url}.',
    'Le reste de la région reconnue couvre la tige inférieure, le renflement et la tige supérieure.':
        'The rest of the recognized region spans the lower stem, bulge and upper stem.',
    'La disposition de reconnaissance du trigger est inconnue.':
        'The trigger recognition layout is unknown.',
    'La disposition de la région reconnue du trigger est invalide.':
        'The trigger recognition layout is invalid.',
    'La région reconnue doit correspondre à la longueur du trigger ({v0} != {v1}).':
        'The recognized region must match the trigger length ({v0} != {v1}).',
    'Région reconnue du trigger : {trigger_length} nt · toehold : {toehold_length} nt · avant AUG : {upper_length} paires · AUG : 3 nt non appariés · après AUG : {lower_length} paires':
        'Recognized trigger region: {trigger_length} nt · toehold: {toehold_length} nt · before AUG: {upper_length} base pairs · AUG: 3 unpaired nt · after AUG: {lower_length} base pairs',
}

TRANSLATIONS.update({
    'Linéarité du trigger seul · malus': 'Isolated trigger linearity · penalty',
    'L’échelle du score doit être brute ou bornée sur 100.': 'The score scale must be raw or bounded from 0 to 100.',
    'Le poids de linéarité du trigger doit être un nombre fini positif ou nul.': 'The trigger linearity weight must be a finite nonnegative number.',
    'Le score affiché est borné de 0 à 100 par une transformation fixe de la somme pondérée. Elle conserve l’ordre des candidats ; la somme brute départage les notes arrondies identiques. Un score brut nul donne 50/100. Ce n’est pas un pourcentage de fonctionnement. Les poids agissent sur le score brut, avant conversion.': 'The displayed score is bounded from 0 to 100 by a fixed transformation of the weighted sum. It preserves candidate order; the raw sum breaks ties between identical rounded scores. A raw score of zero maps to 50/100. This is not a percentage of function. Weights affect the raw score before conversion.',
    'Score sur 100 : plus élevé est meilleur. La somme pondérée brute reste consultable dans le détail du candidat.': 'Score out of 100: higher is better. The raw weighted sum remains available in the candidate details.',
    'Détail du score /100': 'Score details /100',
    'Score /100': 'Score /100',
    'Score brut': 'Raw score',
    'La note /100 convertit la somme pondérée avec une échelle fixe, indépendante des autres candidats. Elle ne mesure pas une probabilité de fonctionnement.': 'The score out of 100 converts the weighted sum using a fixed scale independent of other candidates. It does not measure the probability of function.',
    'Poids enregistré du malus de linéarité du trigger seul : {weight:g}.': 'Saved weight of the isolated trigger linearity penalty: {weight:g}.',
    'La linéarité du trigger seul n’a pas été enregistrée pour ce candidat.': 'Isolated trigger linearity was not saved for this candidate.',
    'Trigger seul · bases libres': 'Isolated trigger · unpaired bases',
    'Trigger seul · défaut de linéarité': 'Isolated trigger · linearity defect',
    '**Fuite : critère secondaire de départage.** Dans votre protocole, les aptamères sont dans un autre tube. Le calcul ci-dessous réunit les trois espèces : c’est un scénario de comparaison, pas une simulation de la fuite ou du transfert entre vos tubes. Il donne le pourcentage des switches présents dans un complexe 1 trigger + 1 switch à l’équilibre, malgré la présence de l’aptamère censé retenir le trigger. Le tube contient aptamère ({trigger_value}), trigger ({trigger_value}) et switch ({switch_value}), soit un ratio aptamère:trigger:switch de {ratio}:{ratio}:1. La concentration d’aptamère est actuellement égale à celle du trigger. Calcul : 100 × [complexe trigger–switch] / [switch] initial. Plus faible est meilleur dans ce scénario secondaire. Le score brut retire w × fuite ; réduire la fuite d’un point de pourcentage apporte w points au score brut. Par défaut, le poids de 0,05 limite l’écart dû à la fuite à 5 points sur toute la plage 0–100 %, pour départager surtout des candidats très proches. C’est un indicateur d’association, pas une mesure de traduction ; le ligand n’est pas simulé.': '**Leak: a secondary tie-breaker.** In your protocol, the aptamers are in a separate tube. The calculation below combines all three species: it is a comparison scenario, not a simulation of leakage or transfer between your tubes. It gives the percentage of switches in a 1 trigger + 1 switch complex at equilibrium, despite the presence of the aptamer intended to retain the trigger. The tube contains aptamer ({trigger_value}), trigger ({trigger_value}) and switch ({switch_value}), giving an aptamer:trigger:switch ratio of {ratio}:{ratio}:1. The aptamer concentration currently equals the trigger concentration. Calculation: 100 × [trigger–switch complex] / initial [switch]. Lower is better in this secondary scenario. The raw score subtracts w × leak; reducing leak by one percentage point adds w points to the raw score. By default, the weight of 0.05 limits the leak contribution to 5 points across the full 0–100% range, mainly to distinguish closely ranked candidates. This is an association indicator, not a measure of translation; the ligand is not simulated.',
    '**Association ON avec le trigger libre.** Pourcentage des switches présents dans un complexe 1 trigger + 1 switch à l’équilibre. Le tube contient seulement le trigger libre ({trigger_value}) et le switch ({switch_value}), soit un ratio trigger:switch de {ratio}:1, sans aptamère ni ligand. Calcul : 100 × [complexe trigger–switch] / [switch] initial. Plus élevé est meilleur pour favoriser la capture du trigger dans le tube de switch, séparé de celui des aptamères. Ce critère est prioritaire sur la fuite : par défaut, un point de pourcentage de ON pèse 60 fois plus qu’un point de fuite. Le score brut ajoute w × ON ; gagner un point de pourcentage apporte w points au score brut. Ce pourcentage est un indicateur du complexe ON ; il ne mesure pas directement l’activation de la traduction ni la quantité de protéine produite.': '**ON association with free trigger.** Percentage of switches in a 1 trigger + 1 switch complex at equilibrium. The tube contains only free trigger ({trigger_value}) and switch ({switch_value}), giving a trigger:switch ratio of {ratio}:1, without aptamer or ligand. Calculation: 100 × [trigger–switch complex] / initial [switch]. Higher is better to favor trigger capture in the switch tube, separate from the aptamer tube. This criterion takes priority over leak: by default, one percentage point of ON carries 60 times more weight than one point of leak. The raw score adds w × ON; gaining one percentage point adds w points to the raw score. This percentage indicates the ON complex; it does not directly measure translation activation or the amount of protein produced.',
    '**Énergie de formation du complexe trigger–switch.** ΔΔG = ΔG MFE du complexe ON − ΔG MFE du switch seul − ΔG MFE du trigger seul, en kcal/mol, avec les conditions choisies. Une valeur plus négative favorise énergétiquement l’association. Le score brut retire w × ΔΔG : une diminution de 1 kcal/mol apporte w points au score brut. Son poids est faible par défaut car ce critère complète le rendement en tube sans décrire à lui seul le fonctionnement. Ce n’est ni une barrière d’activation ni une vitesse de réaction.': '**Trigger–switch complex formation energy.** ΔΔG = MFE ΔG of the ON complex − MFE ΔG of the switch alone − MFE ΔG of the trigger alone, in kcal/mol, under the chosen conditions. A more negative value energetically favors association. The raw score subtracts w × ΔΔG: a decrease of 1 kcal/mol adds w points to the raw score. Its default weight is low because this criterion complements the tube yield without describing function on its own. It is neither an activation barrier nor a reaction rate.',
})

TRANSLATIONS.update({
    'Système aptamère–trigger': 'Aptamer–trigger system',
    'Architecture du toehold switch': 'Toehold switch architecture',
    'Aperçus du toehold switch': 'Toehold switch previews',
    'Calcul et classement du toehold switch': 'Toehold switch design and ranking',
    '##### Poids du classement': '##### Ranking weights',
    'oxDNA / oxNA · conformations et simulation': 'oxDNA / oxNA · conformations and simulation',
    'Sel monovalent Na⁺': 'Monovalent Na⁺ salt',
    'Page du tableau': 'Table page',
    'Lignes {first} à {last} sur {total}.': 'Rows {first}–{last} of {total}.',
    'Succès aux interfaces': 'Successes at interfaces',
    'Balayer Na⁺': 'Sweep Na⁺',
    'Balayer Mg²⁺': 'Sweep Mg²⁺',
    'Tous les candidats du classement enregistré, y compris les exclus.': 'All candidates in the saved ranking, including excluded candidates.',
    'Candidat à explorer': 'Candidate to explore',
    'Les 5 % meilleurs : {shown} sur {total} candidats classés.': 'Top 5%: {shown} of {total} ranked candidates.',
})

TRANSLATIONS.update({
    'Terminer maintenant': 'Finish now',
    'Le temps restant est recalculé après chaque essai terminé.': 'The remaining time is recalculated after each completed trial.',
    'Reprendre le calcul': 'Resume calculation',
    'Terminer avec les résultats actuels': 'Finish with current results',
    'Calcul en pause': 'Calculation paused',
    'Calcul interrompu': 'Calculation interrupted',
    'en pause': 'paused',
    'terminé par anticipation': 'stopped',
    'Finalisation en cours : comparaison des extensions si prévue, classement et exports. Aucun nouvel essai ne sera lancé.': 'Finalizing: extension comparison if requested, ranking and exports. No new trial will start.',
    'Résultats partiels enregistrés · {completed}/{total} essais terminés. Le calcul se poursuit.': 'Partial results saved · {completed}/{total} trials completed. The calculation continues.',
    '**{state} · {completed}/{total} essais terminés**': '**{state} · {completed}/{total} trials completed**',
    'La reprise conserve les séquences, les conditions et les paramètres de ce calcul. Les changements du formulaire ne s’appliquent pas à cette reprise.': 'Resuming preserves this calculation’s sequences, conditions and settings. Changes to the form do not apply to this continuation.',
    'Terminer avec les résultats actuels finalise la comparaison des extensions prévue et les exports, sans créer de nouvel essai.': 'Finishing with current results completes the requested extension comparison and exports without starting another trial.',
    'Pause demandée : l’essai en cours se termine, puis les résultats et les paramètres sont enregistrés sur disque. Vous pourrez reprendre ce calcul après avoir fermé l’application.': 'Pause requested: the current trial will finish, then results and settings will be saved to disk. You can resume this calculation after closing the application.',
    'Un calcul est déjà en cours. Suivez-le dans Résultats.': 'A calculation is already running. Follow it in Results.',
    'Calcul terminé à votre demande · {completed}/{total} essais terminés.': 'Calculation finished at your request · {completed}/{total} trials completed.',
    'Après une fermeture inattendue, seul un essai non encore enregistré peut être recommencé.': 'After an unexpected shutdown, only a trial that has not yet been saved may need to restart.',
    'Le point de reprise est introuvable. Rouvrez le calcul depuis son dossier d’origine.': 'The checkpoint could not be found. Reopen the calculation from its original folder.',
    'Fin demandée : l’essai en cours se termine, puis les candidats obtenus sont classés, la comparaison des extensions est finalisée si elle est prévue, et les exports sont créés.': 'Finish requested: the current trial will finish, then completed candidates will be ranked, the requested extension comparison will be finalized, and exports will be created.',
    'Ces actions attendent la fin de l’essai en cours. Terminer prépare le classement et les exports ; mettre en pause enregistre le calcul pour le reprendre plus tard.': 'These actions wait for the current trial to finish. Finish prepares ranking and exports; pause saves the calculation to resume later.',
    'Le checkpoint ne se trouve pas dans le dossier de résultats de ce design.': 'The checkpoint is not in this design’s results folder.',
    'Ouvrez le checkpoint de conception à reprendre.': 'Open the design checkpoint to resume.',
    'Ce fichier ne contient pas un checkpoint de conception.': 'This file does not contain a design checkpoint.',
    'Seul un calcul en pause ou interrompu peut être repris.': 'Only a paused or interrupted calculation can be resumed.',
    'Le checkpoint ne correspond pas à son dossier de calcul.': 'The checkpoint does not match its calculation folder.',
    'Ce résultat ne contient pas de checkpoint de conception reprenable.': 'This result does not contain a resumable design checkpoint.',
    'Un calcul de conception est déjà en cours pour ce design.': 'A switch design calculation is already running for this design.',
    'Ce calcul de conception est déjà ouvert dans une autre instance de l’application.': 'This design calculation is already open in another instance of the application.',
    'Un calcul de conception utilise déjà ce checkpoint.': 'A design calculation is already using this checkpoint.',
    'Ce checkpoint est déjà utilisé par une autre instance de l’application.': 'This checkpoint is already in use by another instance of the application.',
    'Pause demandée — l’essai en cours sera terminé et sauvegardé.': 'Pause requested — the current trial will finish and be saved.',
    'Fin demandée — l’essai en cours sera terminé et sauvegardé.': 'Finish requested — the current trial will finish and be saved.',
    'Aucun checkpoint de conception n’est disponible.': 'No design checkpoint is available.',
    'Calcul en pause — les essais terminés sont sauvegardés.': 'Calculation paused — completed trials have been saved.',
    'Calcul terminé avec les essais déjà réalisés.': 'Calculation finished with the completed trials.',
    'Conception en pause · les essais terminés sont sauvegardés.': 'Design paused · completed trials have been saved.',
    'Conception arrêtée à la fin de l’essai · classement des résultats…': 'Design stopped after the current trial · ranking results…',
    'La comparaison des extensions sera disponible après la finalisation de ce calcul.': 'Extension comparison will be available once this calculation is finalized.',
    'Aucun candidat terminé n’est encore disponible dans ce point de reprise.': 'No completed candidate is available in this checkpoint yet.',
    'La comparaison des extensions inverses a échoué. Les essais de conception sont sauvegardés ; reprenez le calcul pour réessayer la comparaison.': 'Inverse extension comparison failed. Design trials have been saved; resume the calculation to retry the comparison.',
    'La reprise nécessite un point de sauvegarde de conception chargé.': 'Resuming design requires a loaded RunResult checkpoint.',
    'Les métadonnées d’optimisation du point de reprise sont invalides.': 'The checkpoint optimization metadata is invalid.',
    'Ce résultat n’est pas un point de reprise de conception pris en charge.': 'This result is not a supported design checkpoint.',
    'Ce calcul ne contient aucun point de reprise inachevé.': 'This design has no unfinished checkpoint to resume.',
    'Les paramètres de conception diffèrent de ceux du point de reprise.': 'The design settings differ from the saved checkpoint.',
    'Le moteur thermodynamique ou ses paramètres diffèrent de ceux du point de reprise.': 'The thermodynamic engine or parameters differ from the saved checkpoint.',
    'Le point de reprise appartient à un autre dossier de résultats.': 'The design checkpoint belongs to another output directory.',
    'Le contenu du point de reprise a échoué au contrôle d’intégrité.': 'The design checkpoint contents failed their integrity check.',
    'Un point de reprise plus récent ou différent est déjà enregistré.': 'A newer or different design checkpoint is already saved.',
    'Les paramètres d’optimisation du point de reprise sont incohérents.': 'The checkpoint optimization settings are inconsistent.',
    'Le nombre d’essais terminés du point de reprise est invalide.': 'The checkpoint completed-trial count is invalid.',
    'L’ordre des essais ou leurs graines aléatoires sont incohérents dans le point de reprise.': 'The checkpoint trial sequence or random seeds are inconsistent.',
    'Les candidats du point de reprise sont incohérents ou dupliqués.': 'The checkpoint candidates are inconsistent or duplicated.',
    'Ce point de reprise omet des essais de conception inachevés.': 'This checkpoint omits unfinished design trials.',
})

TRANSLATIONS.update({
    'Trigger–switch avec aptamère (%)': 'Trigger–switch with aptamer (%)',
    'Trigger–switch avec aptamère · malus': 'Trigger–switch with aptamer · penalty',
    'Poids 0 : cette mesure reste informative et n’intervient pas dans le score.': 'Weight 0: this measurement is informative and does not affect the score.',
    'Poids enregistré : {weight:g}. Le score brut retire poids × pourcentage ; réduire cette association d’un point de pourcentage apporte {weight:g} points bruts. Plus faible est meilleur pour ce seul critère.': 'Recorded weight: {weight:g}. The raw score subtracts weight × percentage; reducing this association by one percentage point adds {weight:g} raw points. Lower is better for this criterion alone.',
    'Somme pondérée : + wON × ON − wA × association en présence d’aptamère − 100 × wdéfaut × défaut OFF/ON − wΔΔG × ΔΔG + wRBS × ΔG RBS–linker.': 'Weighted sum: + wON × ON − wA × association with aptamer − 100 × wdefect × OFF/ON defect − wΔΔG × ΔΔG + wRBS × ΔG RBS–linker.',
    'Les calculs de structure et d’association ont été ignorés pour ce candidat. Sa séquence reste consultable et exportée.': 'Structure and association calculations were skipped for this candidate. Its sequence remains available to view and export.',
    '**Association trigger–switch en présence d’aptamère.** Pourcentage des switches dans un complexe contenant exactement 1 trigger et 1 switch, à l’équilibre : 100 × [complexe trigger–switch] / [switch] initial. Le dénominateur comprend tous les switches introduits, libres ou liés. Le tube contient aptamère ({trigger_value}), trigger ({trigger_value}) et switch ({switch_value}), soit un ratio aptamère:trigger:switch de {ratio}:{ratio}:1. La concentration d’aptamère est égale à celle du trigger. Les complexes concurrents jusqu’à {max_size} brins participent au calcul d’équilibre, mais seuls les dimères trigger–switch sont comptés au numérateur. Dans votre protocole, les aptamères sont dans un autre tube : c’est un scénario de comparaison, pas une simulation du transfert entre les tubes. Ce n’est pas une fuite de traduction en OFF ni un pourcentage de protéine produite ; le ligand n’est pas simulé.': '**Trigger–switch association in the presence of aptamer.** Percentage of switches in a complex containing exactly 1 trigger and 1 switch, at equilibrium: 100 × [trigger–switch complex] / initial [switch]. The denominator includes all switches added, free or bound. The tube contains aptamer ({trigger_value}), trigger ({trigger_value}) and switch ({switch_value}), giving an aptamer:trigger:switch ratio of {ratio}:{ratio}:1. The aptamer concentration equals the trigger concentration. Competing complexes of up to {max_size} strands contribute to the equilibrium calculation, but only trigger–switch dimers count in the numerator. In your protocol, aptamers are in a separate tube: this is a comparison scenario, not a simulation of transfer between tubes. This is neither OFF translation leakage nor a percentage of protein produced; the ligand is not simulated.',
    '**Association ON avec le trigger libre.** Pourcentage des switches présents dans un complexe 1 trigger + 1 switch à l’équilibre. Le tube contient seulement le trigger libre ({trigger_value}) et le switch ({switch_value}), soit un ratio trigger:switch de {ratio}:1, sans aptamère ni ligand. Calcul : 100 × [complexe trigger–switch] / [switch] initial. Plus élevé est meilleur pour favoriser la capture du trigger dans le tube de switch, séparé de celui des aptamères. Par défaut, l’association ON participe au score et le critère d’association en présence d’aptamère reste informatif (poids 0). Le score brut ajoute w × ON ; gagner un point de pourcentage apporte w points au score brut. Ce pourcentage est un indicateur du complexe ON ; il ne mesure pas directement l’activation de la traduction ni la quantité de protéine produite.': '**ON association with free trigger.** Percentage of switches in a 1 trigger + 1 switch complex at equilibrium. The tube contains only free trigger ({trigger_value}) and switch ({switch_value}), giving a trigger:switch ratio of {ratio}:1, without aptamer or ligand. Calculation: 100 × [trigger–switch complex] / initial [switch]. Higher is better to favor trigger capture in the switch tube, separate from the aptamer tube. By default, ON association contributes to the score and association in the presence of aptamer remains informative (weight 0). The raw score adds w × ON; gaining one percentage point adds w raw points. This percentage indicates the ON complex; it does not directly measure translation activation or the amount of protein produced.',
    '**Trigger étendu libre, en mode extension inverse.** Le moteur hybride expérimental calcule les probabilités du trigger seul, avec sa chimie et les conditions du design du switch. Le défaut est la fraction moyenne de bases appariées : 1 − moyenne des probabilités non appariées. Le score brut retire w × 100 × ce défaut ; poids 3 par défaut. Une baisse d’un point de pourcentage de bases appariées apporte donc 3 points bruts avec ce réglage. Ce poids est une priorité de classement réglable, pas un optimum expérimental établi. Un poids nul désactive ce malus. Cette mesure porte sur l’ensemble des conformations, pas seulement sur la MFE ; elle ne garantit pas que toutes les bases soient libres simultanément. Elle ne s’applique pas au design sans extension inverse.': '**Free extended trigger, in inverse-extension mode.** The experimental hybrid engine calculates isolated-trigger probabilities using its chemistry and the switch-design conditions. The defect is the mean paired-base fraction: 1 − mean unpaired probability. The raw score subtracts w × 100 × this defect; the default weight is 3. Reducing the paired-base fraction by one percentage point therefore adds 3 raw points at this setting. This weight is an adjustable ranking priority, not an established experimental optimum. A zero weight disables the penalty. This measurement covers the conformational ensemble, not only the MFE; it does not guarantee that all bases are simultaneously unpaired. It does not apply to design without inverse extension.',
    'Priorités par défaut : accessibilité RBS–linker et compromis structural OFF/ON, puis association ON et ΔΔG. En extension inverse, le malus de repliement du trigger seul a un poids de 3. L’association trigger–switch en présence d’aptamère reste informative avec un poids de 0. Les poids restent modifiables ; 0 désactive un critère.': 'Default priorities: RBS–linker accessibility and the OFF/ON structural compromise, then ON association and ΔΔG. In inverse extension, the isolated-trigger folding penalty has weight 3. Trigger–switch association in the presence of aptamer remains informative with weight 0. All weights are adjustable; 0 disables a criterion.',
})

TRANSLATIONS['A : pourcentage de switches associés au trigger en présence d’aptamère.'] = 'A: percentage of switches bound to the trigger in the presence of aptamer.'

TRANSLATIONS.update({
    'RBS ON/OFF · bonus': 'RBS ON/OFF · bonus',
    'M_RBS : mérite normalisé OSTIR, combinant initiation ON et contraste ON/OFF. Voir l’aide du poids RBS ON/OFF pour le modèle et ses limites.': 'M_RBS: normalized OSTIR merit combining ON initiation and ON/OFF contrast. See the RBS ON/OFF weight help for the model and its limits.',
    '**RBS ON/OFF · prédiction exploratoire.** OSTIR 1.1.3 avec ViennaRNA 2.7.2 estime l’initiation de la traduction pour E. coli à 37 °C, selon sa référence ARN à 1,021 M de sel. Le calcul conserve la séquence du candidat dans une fenêtre de 35 bases de chaque côté de l’AUG voulu. OFF utilise l’ARN sans contrainte. En ON, seules les bases prévues pour le duplex avec le trigger sont rendues indisponibles aux appariements ARN internes ; aucune énergie de liaison ADN–ARN n’est ajoutée au modèle OSTIR. C’est une approximation de l’état occupé, pas une validation du complexe hybride. Le bonus combine un ON élevé et un contraste ON/OFF : qON = TIR_ON / (TIR_ON + K), avec K = exp(7,279194329), puis qcontraste = max(0, (TIR_ON − TIR_OFF) / (TIR_ON + TIR_OFF)) et mérite = qON × qcontraste. TIR_ON et TIR_OFF désignent les initiations relatives prédites par OSTIR, pas les pourcentages d’association en tube. Le score brut ajoute poids × 100 × mérite ; poids 3 par défaut. L’affinage ajoute un coût poids × (1 − mérite) en préservant les priorités locales. Cette normalisation et ce poids sont des choix de classement réglables, pas des paramètres optimaux établis expérimentalement. Les valeurs OSTIR ne sont ni un taux absolu de protéine ni un ratio d’induction validé. Un poids nul désactive ce critère.': '**RBS ON/OFF · exploratory prediction.** OSTIR 1.1.3 with ViennaRNA 2.7.2 estimates translation initiation for E. coli at 37 °C, using its RNA reference at 1.021 M salt. The calculation retains the candidate sequence in a window of 35 bases on either side of the intended AUG. OFF uses unconstrained RNA. In ON, only the bases assigned to the trigger duplex are unavailable for internal RNA pairing; no DNA–RNA binding energy is added to the OSTIR model. This approximates the occupied state; it does not validate the hybrid complex. The bonus combines high ON initiation with ON/OFF contrast: qON = TIR_ON / (TIR_ON + K), with K = exp(7.279194329), then qcontrast = max(0, (TIR_ON − TIR_OFF) / (TIR_ON + TIR_OFF)) and merit = qON × qcontrast. TIR_ON and TIR_OFF denote OSTIR relative initiation predictions, not tube association percentages. The raw score adds weight × 100 × merit; the default weight is 3. Refinement adds weight × (1 − merit) while preserving local priorities. This normalization and weight are adjustable ranking choices, not experimentally established optimal parameters. OSTIR values are neither absolute protein production rates nor a validated induction ratio. A zero weight disables this criterion.',
    'Priorités par défaut : accessibilité RBS–linker et compromis structural OFF/ON, puis association ON et ΔΔG. Le bonus RBS ON/OFF a un poids de 3. En extension inverse, le malus de repliement du trigger seul a un poids de 3. L’association trigger–switch en présence d’aptamère reste informative avec un poids de 0. Les poids restent modifiables ; 0 désactive un critère.': 'Default priorities: RBS–linker accessibility and the OFF/ON structural compromise, then ON association and ΔΔG. The RBS ON/OFF bonus has weight 3. In inverse extension, the isolated-trigger folding penalty has weight 3. Trigger–switch association in the presence of aptamer remains informative with weight 0. All weights are adjustable; 0 disables a criterion.',
})
