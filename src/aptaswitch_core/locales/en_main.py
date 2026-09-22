"""English translations for the main interface."""

TRANSLATIONS = {
    'indéfini':
        'undefined',
    '<div class="apta-subtitle">Concevoir, classer et explorer des toehold switches activés par aptamère.</div>':
        '<div class="apta-subtitle">Design, rank and explore aptamer-activated toehold switches.</div>',
    'Navigation':
        'Navigation',
    'Calcul en cours · {current}/{total}':
        'Calculation running · {current}/{total}',
    'Le toehold activateur standard nécessite un trigger d’au moins {trigger_region_length} nt. Étendez le trigger en section 2 ou activez « Je veux personnaliser mon tSwitch ».':
        'The standard toehold activator requires a trigger of at least {trigger_region_length} nt. Extend the trigger in section 2 or enable “Customize my tSwitch”.',
    'La région du trigger utilisée pour le switch doit contenir au moins 2 bases.':
        'The trigger region used for the switch must contain at least 2 bases.',
    'Configurez NUPACK dans Réglages pour lancer la conception.':
        'Configure NUPACK in Settings to start designing.',
    "Ajoutez un trigger valide pour générer l'animation.":
        'Enter a valid trigger to generate the animation.',
    'Lecture automatique':
        'Autoplay',
    'Étape':
        'Step',
    '1 · Système initial : données et structure 2D':
        '1 · Initial system: data and 2D structure',
    'Renseignez les séquences disponibles, les conditions du système et, si elles sont connues, les bases impliquées dans la liaison au ligand.':
        'Enter the available sequences, system conditions and, if known, the bases involved in ligand binding.',
    '#### Séquences du système':
        '#### System sequences',
    'Molécule cible':
        'Target molecule',
    "Séquence de l'aptamère":
        'Aptamer sequence',
    'Séquence ADN ou ARN — espaces et retours à la ligne acceptés':
        'DNA or RNA sequence — spaces and line breaks accepted',
    'Séquence du trigger':
        'Trigger sequence',
    'Séquence ADN ou ARN':
        'DNA or RNA sequence',
    'Aptamère : {value} nt · Trigger : {value1} nt':
        'Aptamer: {value} nt · Trigger: {value1} nt',
    'Attention : votre trigger mesure {value} nt. Sa longueur peut être trop faible ; vous pouvez l’étendre dans la section 2.':
        'Warning: your trigger is {value} nt long. It may be too short; you can extend it in section 2.',
    'Bases de l’aptamère impliquées dans la liaison au ligand (si connues)':
        'Aptamer bases involved in ligand binding (if known)',
    '#### Région du trigger utilisée pour le switch':
        '#### Trigger region used for the switch',
    'Première base':
        'First base',
    'Dernière base':
        'Last base',
    'Exclure la zone liée':
        'Exclude the bound region',
    'La zone liée est identifiée dans la structure du système initial, calculée en section 1.':
        'The bound region is identified in the initial system structure calculated in section 1.',
    '3 · Architecture du switch':
        '3 · Switch architecture',
    'Leader 5′':
        '5′ leader',
    'RBS':
        'RBS',
    'La boucle cible contient 3 bases libres (NNN), optimisées par NUPACK et non appariées, suivies de votre RBS.':
        'The target loop contains 3 free bases (NNN), optimized by NUPACK to remain unpaired, followed by your RBS.',
    'Chimie / modèle':
        'Thermodynamic model',
    'Je veux personnaliser mon tSwitch':
        'Customize my tSwitch',
    'La personnalisation de l’architecture est déconseillée : elle peut modifier le repliement et le fonctionnement du switch. Conservez les paramètres standards sauf besoin particulier.':
        'Customizing the architecture is not recommended: it may affect switch folding and function. Keep the standard parameters unless you have a specific need.',
    'Longueur du toehold':
        'Toehold length',
    'Le reste de la région du trigger forme le stem inférieur.':
        'The rest of the trigger region forms the lower stem.',
    '\n                <div style="height:42px;display:flex;border-radius:8px;overflow:hidden;color:white;font-weight:700;">\n                  <div style="width:{toe_pct:.2f}%;background:#1258dc;padding:.65rem;white-space:nowrap">Toehold · {toehold_length} nt</div>\n                  <div style="flex:1;background:#0e9f76;padding:.65rem;white-space:nowrap">Stem · {stem_length} nt</div>\n                </div>\n                ':
        '\n                <div style="height:42px;display:flex;border-radius:8px;overflow:hidden;color:white;font-weight:700;">\n                  <div style="width:{toe_pct:.2f}%;background:#1258dc;padding:.65rem;white-space:nowrap">Toehold · {toehold_length} nt</div>\n                  <div style="flex:1;background:#0e9f76;padding:.65rem;white-space:nowrap">Stem · {stem_length} nt</div>\n                </div>\n                ',
    'Distance entre le RBS et AUG (nt)':
        'Distance between RBS and AUG (nt)',
    'Nombre total de bases N entre la fin du RBS et AUG. Ces bases, optimisées par NUPACK, forment le bras droit de la tige supérieure dans la cible OFF.':
        'Total number of N bases between the end of the RBS and AUG. NUPACK optimizes these bases, which form the right arm of the upper stem in the OFF target.',
    'Linker automatique':
        'Automatic linker',
    'La longueur est calculée pour préserver le cadre de lecture (0, 1 ou 2 bases). NUPACK choisit ces bases pour se rapprocher des structures cibles.':
        'The length is calculated to preserve the reading frame (0, 1 or 2 bases). NUPACK chooses these bases to match the target structures.',
    'Linker personnalisé':
        'Custom linker',
    'N laisse une base libre pour NUPACK. A, C, G, U ou T imposent une base si vous souhaitez contraindre le linker.':
        'N leaves a base free for NUPACK. A, C, G, U or T fixes a base if you want to constrain the linker.',
    'Le modèle standard nécessite un trigger d’au moins {trigger_region_length} nt. Étendez-le en section 2 ou personnalisez le tSwitch.':
        'The standard model requires a trigger of at least {trigger_region_length} nt. Extend it in section 2 or customize the tSwitch.',
    '4 · Aperçus du switch en direct':
        '4 · Live switch previews',
    'Structures cibles construites à partir de vos séquences. N indique les bases à optimiser par NUPACK ; les structures prédites seront disponibles après le calcul.':
        'Target structures built from your sequences. N marks bases to be optimized by NUPACK; predicted structures will be available after calculation.',
    'Animation de liaison':
        'Binding animation',
    'État du switch':
        'Switch state',
    "Générer l'animation schématique":
        'Generate the schematic animation',
    'La génération peut prendre quelques secondes pour les séquences longues.':
        'Generation may take a few seconds for long sequences.',
    'Préparation des images…':
        'Preparing frames…',
    "Activez l'animation pour suivre l'accrochage au toehold puis l'ouverture du stem.":
        'Enable the animation to follow toehold binding and stem opening.',
    '**Accessibilité de la région RBS–linker.** Énergie libre de repliement MFE, en kcal/mol, calculée par NUPACK sur cette région isolée : de la première base de la boucle contenant le RBS jusqu’à la base précédant le premier nucléotide du rapporteur. Les trois bases optimisées avant le RBS (NNN dans la cible) sont incluses ; la tige située avant la boucle et le rapporteur sont exclus. La région inclut le RBS, l’espace RBS–AUG, AUG, le bras aval du stem et le linker éventuel. Un ΔG moins négatif, proche de 0, est favorisé : moins de structure secondaire à défaire pour la traduction. Ce poids ajoute w × ΔG au score brut : gagner 1 kcal/mol ajoute w points. [Définition de Green, tableau S3](https://yin.hms.harvard.edu/publications/2014.tswitch1.table.s3.xlsx).':
        '**RBS–linker accessibility.** MFE folding free energy, in kcal/mol, calculated by NUPACK for this isolated region: from the first base of the loop containing the RBS to the base immediately before the first reporter nucleotide. The three optimized bases before the RBS (NNN in the target) are included; the stem before the loop and the reporter are excluded. The region includes the RBS, the RBS–AUG spacer, AUG, the downstream stem arm and any linker. A less negative ΔG, close to 0, is preferred: less secondary structure needs to unfold for translation. This weight adds w × ΔG to the raw score: gaining 1 kcal/mol adds w points. [Green definition, Table S3](https://yin.hms.harvard.edu/publications/2014.tswitch1.table.s3.xlsx).',
    '**Fidélité aux cibles OFF et ON.** Le défaut d’ensemble normalisé retourné par NUPACK mesure l’écart aux structures et aux concentrations cibles des tubes de design : switch seul en OFF ; trigger + switch en ON. Il tient compte des mauvais appariements et du manque de complexes aux concentrations visées. Plus faible est meilleur ; 0 est idéal. Un défaut de 0,05 correspond à 5 %. Le score retire w × (100 × défaut) : réduire le défaut d’un point de pourcentage apporte w points au score brut. Ce n’est pas une mesure d’expression.':
        '**Agreement with the OFF and ON targets.** The normalized ensemble defect returned by NUPACK measures deviation from the target structures and concentrations in the design tubes: switch alone in OFF; trigger + switch in ON. It accounts for incorrect base pairing and insufficient complexes at the intended concentrations. Lower is better; 0 is ideal. A defect of 0.05 corresponds to 5%. The score subtracts w × (100 × defect): reducing the defect by one percentage point adds w points to the raw score. This is not a measure of expression.',
    '**Fuite : critère secondaire de départage.** Dans votre protocole, les aptamères sont dans un autre tube. Le calcul ci-dessous réunit les trois espèces : c’est un scénario de comparaison, pas une simulation de la fuite ou du transfert entre vos tubes. Il donne le pourcentage des switches présents dans un complexe 1 trigger + 1 switch à l’équilibre, malgré la présence de l’aptamère censé retenir le trigger. Le tube contient aptamère ({trigger_value}), trigger ({trigger_value}) et switch ({switch_value}), soit un ratio aptamère:trigger:switch de {ratio}:{ratio}:1. La concentration d’aptamère est actuellement égale à celle du trigger. Calcul : 100 × [complexe trigger–switch] / [switch] initial. Plus faible est meilleur dans ce scénario secondaire. Le score retire w × fuite ; réduire la fuite d’un point de pourcentage apporte w points au score. Par défaut, le poids de 0,05 limite l’écart dû à la fuite à 5 points sur toute la plage 0–100 %, pour départager surtout des candidats très proches. C’est un indicateur d’association, pas une mesure de traduction ; le ligand n’est pas simulé.':
        '**Leak: a secondary tie-breaker.** In your protocol, the aptamers are in a separate tube. The calculation below combines all three species: it is a comparison scenario, not a simulation of leakage or transfer between your tubes. It gives the percentage of switches in a 1 trigger + 1 switch complex at equilibrium, despite the presence of the aptamer intended to retain the trigger. The tube contains aptamer ({trigger_value}), trigger ({trigger_value}) and switch ({switch_value}), giving an aptamer:trigger:switch ratio of {ratio}:{ratio}:1. The aptamer concentration currently equals the trigger concentration. Calculation: 100 × [trigger–switch complex] / initial [switch]. Lower is better in this secondary scenario. The score subtracts w × leak; reducing leak by one percentage point adds w points to the score. By default, the weight of 0.05 limits the leak contribution to 5 points across the full 0–100% range, mainly to distinguish closely ranked candidates. This is an association indicator, not a measure of translation; the ligand is not simulated.',
    '**Association ON avec le trigger libre.** Pourcentage des switches présents dans un complexe 1 trigger + 1 switch à l’équilibre. Le tube contient seulement le trigger libre ({trigger_value}) et le switch ({switch_value}), soit un ratio trigger:switch de {ratio}:1, sans aptamère ni ligand. Calcul : 100 × [complexe trigger–switch] / [switch] initial. Plus élevé est meilleur pour favoriser la capture du trigger dans le tube de switch, séparé de celui des aptamères. Ce critère est prioritaire sur la fuite : par défaut, un point de pourcentage de ON pèse 60 fois plus qu’un point de fuite. Le score ajoute w × ON ; gagner un point de pourcentage apporte w points au score brut. Ce pourcentage est un indicateur du complexe ON ; il ne mesure pas directement l’activation de la traduction ni la quantité de protéine produite.':
        '**ON association with free trigger.** Percentage of switches in a 1 trigger + 1 switch complex at equilibrium. The tube contains only free trigger ({trigger_value}) and switch ({switch_value}), giving a trigger:switch ratio of {ratio}:1, without aptamer or ligand. Calculation: 100 × [trigger–switch complex] / initial [switch]. Higher is better to favor trigger capture in the switch tube, separate from the aptamer tube. This criterion takes priority over leak: by default, one percentage point of ON carries 60 times more weight than one point of leak. The score adds w × ON; gaining one percentage point adds w points to the raw score. This percentage indicates the ON complex; it does not directly measure translation activation or the amount of protein produced.',
    '**Énergie de formation du complexe trigger–switch.** ΔΔG = ΔG MFE du complexe ON − ΔG MFE du switch seul − ΔG MFE du trigger seul, en kcal/mol, avec les conditions choisies. Une valeur plus négative favorise énergétiquement l’association. Le score retire w × ΔΔG : une diminution de 1 kcal/mol apporte w points au score brut. Son poids est faible par défaut car ce critère complète le rendement en tube sans décrire à lui seul le fonctionnement. Ce n’est ni une barrière d’activation ni une vitesse de réaction.':
        '**Trigger–switch complex formation energy.** ΔΔG = MFE ΔG of the ON complex − MFE ΔG of the switch alone − MFE ΔG of the trigger alone, in kcal/mol, under the chosen conditions. A more negative value energetically favors association. The score subtracts w × ΔΔG: a decrease of 1 kcal/mol adds w points to the raw score. Its default weight is low because this criterion complements the tube yield without describing function on its own. It is neither an activation barrier nor a reaction rate.',
    '5 · Calcul et scoring du switch':
        '5 · Switch calculation and scoring',
    'Conditions':
        'Conditions',
    'Scoring':
        'Scoring',
    'Sortie':
        'Output',
    'Moteur de calcul : NUPACK':
        'Calculation engine: NUPACK',
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
    'Priorités par défaut : accessibilité RBS–linker et fidélité aux cibles, puis association ON et ΔΔG. La fuite a un poids très faible pour départager les candidats proches, les aptamères étant dans un autre tube. Les poids restent modifiables ; 0 désactive un critère.':
        'Default priorities: RBS–linker accessibility and target agreement, followed by ON association and ΔΔG. Leak has a very low weight to distinguish closely ranked candidates, as the aptamers are in a separate tube. You can adjust the weights; 0 disables a criterion.',
    'Voir la formule':
        'Show the formula',
    'Le score conserve la somme pondérée, sans plafond ni plancher : plus il est élevé, mieux le candidat est classé. Il peut dépasser 100 ou être négatif ; ce n’est pas un pourcentage. Les poids sont des coefficients de classement, pas des pourcentages d’importance ni des coefficients publiés par Green.':
        'The score retains the weighted sum, without a ceiling or floor: higher scores rank better. It can exceed 100 or be negative; it is not a percentage. The weights are ranking coefficients, not percentages of importance or coefficients published by Green.',
    'Exclure les candidats avec un STOP prématuré':
        'Exclude candidates with a premature STOP',
    'Dès la séquence générée, les candidats avec un STOP avant le reporter sont exclus. Ils restent en fin de tableau, sans score ni analyses thermodynamiques supplémentaires.':
        'As soon as the sequence is generated, candidates with a STOP before the reporter are excluded. They remain at the bottom of the table, without a score or further thermodynamic analyses.',
    'Exclure les candidats sans linéarité parfaite du RBS–linker en ON':
        'Exclude candidates without a fully unpaired RBS–linker region in ON',
    'Toutes les bases, de la première base de la boucle contenant le RBS (y compris les bases optimisées avant le RBS) jusqu’à la dernière base avant le rapporteur, doivent être non appariées dans les structures MFE ON renvoyées par NUPACK pour le complexe trigger–switch. Un seul appariement dans cette région suffit à exclure le candidat. Ce critère porte sur ces structures prédites, pas sur l’ensemble des conformations.':
        'All bases, from the first base of the loop containing the RBS (including the optimized bases before the RBS) to the last base before the reporter, must be unpaired in the ON MFE structures returned by NUPACK for the trigger–switch complex. A single paired base in this region excludes the candidate. This criterion applies to these predicted structures, not to the full conformational ensemble.',
    'Criblage de la structure ON après le filtre STOP, avant les autres analyses. Les candidats exclus restent en fin de tableau, sans score ni calcul de rendement ou de fuite.':
        'The ON structure is screened after the STOP filter, before the other analyses. Excluded candidates remain at the bottom of the table, without a score, yield calculation or leak calculation.',
    'Dossier local des runs':
        'Local runs folder',
    'Les designs et extensions sont rangés dans les sous-dossiers designs/ et extensions/. Vous pouvez les rouvrir dans Résultats.':
        'Designs and extensions are stored in the designs/ and extensions/ subfolders. You can reopen them in Results.',
    'Exports : JSON, CSV, XLSX, SVG, journal et projet .aptaswitch.json.':
        'Exports: JSON, CSV, XLSX, SVG, log and .aptaswitch.json project.',
    'Avec NUPACK, plus de 50 essais peuvent prendre plusieurs heures ou plusieurs jours.':
        'With NUPACK, more than 50 trials may take several hours or several days.',
    'Je confirme ce calcul long':
        'I confirm this long calculation',
    'Lancer la conception':
        'Start design',
    '## Conception':
        '## Design',
    '<div class="apta-subtitle">Définissez les séquences, contrôlez l’architecture puis lancez le classement.</div>':
        '<div class="apta-subtitle">Define the sequences, check the architecture and start ranking.</div>',
    'Moteur':
        'Engine',
    'Progression':
        'Progress',
    '{current} / {total}':
        '{current} / {total}',
    'Temps écoulé':
        'Elapsed time',
    'Temps restant estimé':
        'Estimated time remaining',
    '#### Console du moteur':
        '#### Engine console',
    'Le temps restant est recalculé après chaque candidat terminé.':
        'The remaining time is recalculated after each completed candidate.',
    'Annuler le calcul':
        'Cancel calculation',
    'Rang':
        'Rank',
    'Statut':
        'Status',
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
    '## Résultats':
        '## Results',
    'Type de résultats':
        'Result type',
    '<div class="apta-subtitle">Classement, structures calculées et exports.</div>':
        '<div class="apta-subtitle">Ranking, calculated structures and exports.</div>',
    'Ouvrez un run enregistré ci-dessus ou lancez un calcul depuis Conception.':
        'Open a saved run above or start a calculation from Design.',
    'Run':
        'Run',
    'Candidats':
        'Candidates',
    '{excluded_count} candidat(s) exclu(s) par le criblage, affiché(s) en fin de tableau avec leur motif. — = non calculé ; les exclus ne reçoivent ni rang ni score.':
        '{excluded_count} candidate(s) excluded by screening, shown at the bottom of the table with the reason. — = not calculated; excluded candidates receive no rank or score.',
    'Tous les candidats sont exclus par le criblage. Aucun candidat n’a été retenu pour le classement.':
        'All candidates were excluded by screening. No candidate was retained for ranking.',
    'Tableau de classement en lecture seule.':
        'Read-only ranking table.',
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
    'Codon STOP prématuré — codon {first_stop}':
        'Premature STOP codon — codon {first_stop}',
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
    '#### NUPACK · design et thermodynamique':
        '#### NUPACK · design and thermodynamics',
    'Dossier NUPACK, .local_deps ou wheel':
        'NUPACK folder, .local_deps or wheel',
    '/chemin/vers/nupack':
        '/path/to/nupack',
    "Je confirme être responsable de la licence et de l'autorisation d'utilisation de NUPACK.":
        'I confirm that I am responsible for the NUPACK license and permission to use it.',
    'Valider NUPACK':
        'Validate NUPACK',
    "Indiquez d'abord un chemin NUPACK.":
        'Enter a NUPACK path first.',
    'Confirmez la responsabilité de licence avant de continuer.':
        'Confirm responsibility for the license before continuing.',
    "Validation de l'import NUPACK…":
        'Validating the NUPACK import…',
    'NUPACK activé depuis {nupack_path}':
        'NUPACK enabled from {nupack_path}',
    'Configurez NUPACK pour lancer les calculs de structure, d’extension et de conception.':
        'Configure NUPACK to run structure, extension and design calculations.',
    'Bibliothèque de reporters':
        'Reporter library',
    '← Retour à la conception':
        '← Back to design',
    'La séquence commence juste après le codon de départ AUG du switch.':
        'The sequence starts immediately after the switch AUG start codon.',
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
    'Séquence ARN/ADN après AUG':
        'RNA/DNA sequence after AUG',
    'Ajouter':
        'Add',
    'Reporter {label} ajouté.':
        'Reporter {label} added.',
    'Préférences locales':
        'Local preferences',
    'Par défaut, les runs sont enregistrés dans le dossier de l’application, sous runs/designs ou runs/extensions.':
        'By default, runs are saved in the application folder, under runs/designs or runs/extensions.',
    "Dossier d'export par défaut":
        'Default export folder',
    'Enregistrer':
        'Save',
    'Le dossier ne peut pas être vide.':
        'The folder cannot be empty.',
    'Dossier par défaut enregistré.':
        'Default folder saved.',
    '#### Réinitialisation':
        '#### Reset',
    'Efface la configuration des moteurs, les reporters personnalisés et le dossier par défaut. Les exports ne sont pas supprimés.':
        'Clears engine configuration, custom reporters and the default folder. Exports are not deleted.',
    'Je confirme la réinitialisation de tous les réglages':
        'I confirm resetting all settings',
    'Réinitialiser les réglages':
        'Reset settings',
    '## Réglages':
        '## Settings',
    '<div class="apta-subtitle">Dépendances scientifiques, reporters et préférences enregistrées sur cette machine.</div>':
        '<div class="apta-subtitle">Scientific dependencies, reporters and preferences saved on this machine.</div>',
    'Conception':
        'Design',
    'Résultats':
        'Results',
    'Réglages':
        'Settings',
    'Langue':
        'Language',
    'Dépendances':
        'Dependencies',
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
    'Defect · malus':
        'Defect · penalty',
    'ON % · bonus':
        'ON % · bonus',
    'ΔΔG · liaison':
        'ΔΔG · binding',
    'Fuite · départage':
        'Leak · tie-breaker',
    'En attente des premiers messages…':
        'Waiting for the first messages…',
    'Le calcul a échoué.':
        'The calculation failed.',
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
    'Design de switches':
        'Switch design',
    'Extension de triggers':
        'Trigger extension',
    'Tous les exports (.zip)':
        'All exports (.zip)',
    'Exporter les résultats':
        'Export results',
    '{status} · essai {trial}':
        '{status} · trial {trial}',
    'Rang {rank} · essai {trial} · score {score:.2f}':
        'Rank {rank} · trial {trial} · score {score:.2f}',
    'Apta2Switch-Studio — interface locale de design de toehold switches.':
        'Apta2Switch-Studio — local interface for designing toehold switches.',
    'Switch ARN + trigger ADN (rna-dna06)':
        'RNA switch + DNA trigger (rna-dna06)',
    'Switch ADN + trigger ARN (rna-dna06)':
        'DNA switch + RNA trigger (rna-dna06)',
    'Switch ADN + trigger ADN (dna04)':
        'DNA switch + DNA trigger (dna04)',
    'Switch ARN + trigger ARN (rna06)':
        'RNA switch + RNA trigger (rna06)',
    'Switch ARN + trigger ARN (rna95)':
        'RNA switch + RNA trigger (rna95)',
    'Switch ARN + trigger ARN (rna99)':
        'RNA switch + RNA trigger (rna99)',
    'Impossible d’ouvrir le dossier NUPACK : {error}':
        'Unable to open the NUPACK folder: {error}',
    'Téléchargez NUPACK, décompressez son archive et placez le dossier obtenu dans le dossier ci-dessous. Puis confirmez votre licence et cliquez sur Valider NUPACK. Aucune commande ni installation de Python n’est nécessaire.':
        'Download NUPACK, extract its archive and place the extracted folder in the folder below. Then confirm your license and click Validate NUPACK. No commands or Python installation are needed.',
    'Dossier NUPACK : {path}':
        'NUPACK folder: {path}',
    'Télécharger NUPACK':
        'Download NUPACK',
    'Ouvrir le dossier NUPACK':
        'Open NUPACK folder',
    'Dossier NUPACK ou fichier .whl':
        'NUPACK folder or .whl file',
    'Un composant de l’application est manquant. Téléchargez à nouveau l’application complète depuis la page des versions.':
        'An application component is missing. Download the complete application again from the releases page.',
    'Impossible d’ouvrir le navigateur. Ouvrez {url}.':
        'Unable to open the browser. Open {url}.',
}
