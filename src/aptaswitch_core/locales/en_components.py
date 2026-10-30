"""English translations for reusable interface components."""

TRANSLATIONS = {
    'Modèle de l’extension': 'Extension model',
    'Saisissez l’aptamère pour cocher les bases impliquées dans la liaison au ligand.':
        'Enter the aptamer sequence to select the bases involved in ligand binding.',
    'Cochez les bases connues, de 5′ vers 3′. Les positions commencent à 1 ; plusieurs régions séparées peuvent être sélectionnées.':
        'Select the known bases from 5′ to 3′. Positions start at 1; you can select several separate regions.',
    'Bases sélectionnées : ':
        'Selected bases: ',
    'Aucune base sélectionnée. L’extension reste possible sans annotation du ligand.':
        'No bases selected. You can extend the trigger without ligand annotations.',
    'Effacer la sélection du ligand':
        'Clear ligand selection',
    'Ces annotations repèrent la région de liaison dans les schémas et les résultats ; elles ne modélisent pas le ligand et ne modifient pas le score d’extension.':
        'These annotations highlight the binding region in diagrams and results; they do not model the ligand or change the extension score.',
    '#### Aperçu en direct de l’extension':
        '#### Live extension preview',
    'Aperçu de l’extension SVG':
        'Extension preview SVG',
    'Aperçu de l’extension PNG':
        'Extension preview PNG',
    '#### Conditions du système':
        '#### System conditions',
    'ADN':
        'DNA',
    'ARN':
        'RNA',
    'Température du système (°C)':
        'System temperature (°C)',
    'Choisissez des conditions aussi proches que possible des conditions expérimentales ayant démontré l’efficacité du système aptamère–trigger pour votre cible. Elles seront utilisées pour la structure initiale et l’extension éventuelle.':
        'Choose conditions as close as possible to the experimental conditions that demonstrated the effectiveness of the aptamer–trigger system for your target. They will be used for the initial structure and any extension.',
    '#### Aperçus du système initial':
        '#### Initial system previews',
    'La structure initiale sera actualisée à la fin du calcul en cours.':
        'The initial structure will update when the current calculation finishes.',
    'Saisissez l’aptamère et le trigger pour afficher leur structure 2D.':
        'Enter the aptamer and trigger sequences to display their 2D structure.',
    'La structure initiale n’a pas pu être calculée : ':
        'The initial structure could not be calculated: ',
    'Réessayer l’analyse du système':
        'Retry system analysis',
    'Système initial · aptamère + trigger':
        'Initial system · aptamer + trigger',
    'Structure initiale SVG':
        'Initial structure SVG',
    'Structure initiale PNG':
        'Initial structure PNG',
    'Séquences du système initial':
        'Initial system sequences',
    'Séquence initiale SVG':
        'Initial sequence SVG',
    'Séquence initiale PNG':
        'Initial sequence PNG',
    'Un calcul est déjà en cours. Attendez sa fin ou annulez-le dans Résultats.':
        'A calculation is already running. Wait for it to finish or cancel it in Results.',
    'Choisissez une longueur finale supérieure à celle du trigger original.':
        'Choose a final length greater than the original trigger length.',
    "La longueur cible est limitée à 30 nt. Activez l'option de dépassement pour aller au-delà.":
        'The target length is limited to 30 nt. Enable longer extensions to go beyond this limit.',
    'Extension exhaustive':
        'Exhaustive extension',
    'Testez toutes les combinaisons de bases ajoutées, puis choisissez un trigger étendu pour concevoir le switch.':
        'Test every combination of added bases, then select an extended trigger for switch design.',
    'Extrémité à étendre':
        'End to extend',
    'Répartition égale des bases ajoutées ; si le nombre est impair, une base de plus est ajoutée en 5′.':
        'Added bases are split equally; if the number is odd, one extra base is added at the 5′ end.',
    'Longueurs finales':
        'Final lengths',
    'Deux versions : 24 et 30 nt':
        'Two versions: 24 and 30 nt',
    'Longueur personnalisée':
        'Custom length',
    'Autoriser une longueur supérieure à 30 nt':
        'Allow lengths above 30 nt',
    'Longueur cible (nt)':
        'Target length (nt)',
    'Attention : au-delà de 30 nt, le calcul peut prendre beaucoup de temps. Chaque base ajoutée multiplie la recherche par 4.':
        'Caution: calculations above 30 nt may take a long time. Every added base multiplies the search space by 4.',
    'Cette recherche exhaustive comporte des millions de candidats ou davantage ; elle peut durer plusieurs heures ou jours.':
        'This exhaustive search contains millions of candidates or more and may take several hours or days.',
    'Choisissez une longueur personnalisée supérieure à la longueur actuelle.':
        'Choose a custom length greater than the current length.',
    'Paramètres du calcul d’extension':
        'Extension calculation settings',
    'Candidats conservés par longueur':
        'Candidates retained per length',
    'Maximum de candidats en analyse d’ensemble par longueur':
        'Maximum candidates for ensemble analysis per length',
    'Toutes les extensions passent le criblage MFE. Si cette limite est atteinte, seule une partie des meilleurs ex æquo passe l’analyse d’ensemble ; cette limitation sera indiquée dans les résultats.':
        'All extensions undergo MFE screening. If this limit is reached, only some of the best tied candidates undergo ensemble analysis; this limitation will be recorded in the results.',
    'Étendre le trigger':
        'Extend trigger',
    'Un calcul est en cours. Son suivi et son annulation sont disponibles dans Résultats.':
        'A calculation is running. You can track or cancel it in Results.',
    'Extension annulée · résultats partiels.':
        'Extension cancelled · partial results.',
    'Moteur':
        'Engine',
    'Progression de l’étape':
        'Stage progress',
    'Temps écoulé':
        'Elapsed time',
    'Temps restant de l’étape':
        'Time remaining for this stage',
    '#### Console de l’extension':
        '#### Extension console',
    'Préparation…':
        'Preparing…',
    'Annuler l’extension':
        'Cancel extension',
    'Lancez une extension depuis les séquences de la page Conception.':
        'Start an extension using the sequences on the Design page.',
    'Système personnalisé':
        'Custom system',
    'Bases de l’aptamère liées au ligand (cercles roses) : ':
        'Ligand-binding aptamer bases (pink circles): ',
    'Extensions criblées':
        'Extensions screened',
    'Analyses d’ensemble':
        'Ensemble analyses',
    'Candidats conservés':
        'Candidates retained',
    'Un score plus bas est meilleur. Ces calculs évaluent la conservation de structure du complexe à deux brins ; ils ne prédisent pas le rendement du switch complet.':
        'Lower scores are better. These calculations assess structural preservation of the two-strand complex; they do not predict the performance of the complete switch.',
    '#### Aperçus · originale et versions étendues':
        '#### Previews · original and extended versions',
    'Utiliser pour le design du switch':
        'Use for switch design',
    'Le design de switches utilise actuellement un aptamère ADN. Les résultats ARN restent exportables.':
        'Switch design currently uses a DNA aptamer. RNA results can still be exported.',
    'Le graphique montre les 10 % de candidats disponibles ayant les meilleurs scores pour cette longueur (effectif arrondi au supérieur). Le tableau conserve tous les candidats.':
        'The chart shows the top 10% of available candidates for this length by score (rounded up). The table contains all candidates.',
    'Paysage de sélection SVG':
        'Selection landscape SVG',
    'Paysage de sélection PNG':
        'Selection landscape PNG',
    'Probabilités des paires SVG':
        'Pair probabilities SVG',
    'Probabilités des paires PNG':
        'Pair probabilities PNG',
    'La référence ne contient aucune paire MFE à comparer.':
        'The reference contains no MFE pairs to compare.',
    '#### Tableau des triggers étendus':
        '#### Extended trigger table',
    'Rang / longueur':
        'Rank / length',
    'Longueur (nt)':
        'Length (nt)',
    'Trigger étendu 5′ → 3′':
        'Extended trigger 5′ → 3′',
    'RMSE cœur':
        'Core RMSE',
    'Perte des paires':
        'Pair probability loss',
    'Appariement extension':
        'Extension pairing',
    'Extension ↔ aptamère':
        'Extension ↔ aptamer',
    'Paires MFE perdues':
        'Lost MFE pairs',
    'Nouvelles paires du cœur':
        'New core pairs',
    'Structure 2D (dot-bracket)':
        '2D structure (dot-bracket)',
    'Aucun candidat final disponible pour ce calcul.':
        'No final candidates are available for this calculation.',
    'Tous les exports d’extension (.zip)':
        'All extension exports (.zip)',
    'Exporter les résultats d’extension':
        'Export extension results',
    'Design de switches':
        'Switch design',
    'Extension de triggers':
        'Trigger extension',
    'Attendez la fin du calcul en cours avant d’ouvrir un autre run.':
        'Wait for the current calculation to finish before opening another run.',
    'Choisissez un run ou indiquez son chemin.':
        'Choose a run or enter its path.',
    'Ouvrir un run enregistré':
        'Open a saved run',
    'Les nouveaux runs sont rangés dans designs/ ou extensions/. Les anciens dossiers sont aussi recherchés.':
        'New runs are stored in designs/ or extensions/. Older folders are also searched.',
    'Run à ouvrir':
        'Run to open',
    'Aucun run enregistré trouvé':
        'No saved runs found',
    'Ouvrir ce run':
        'Open this run',
    'Actualiser la liste':
        'Refresh list',
    'Afficher le dossier':
        'Show folder',
    'Ouvrir depuis un chemin':
        'Open from a path',
    'Dossier du run ou fichier de résultats':
        'Run folder or results file',
    '/chemin/vers/le/run':
        '/path/to/run',
    'Accepte le dossier d’un run, son fichier *_results.json ou son projet .aptaswitch.json.':
        'Accepts a run folder, its *_results.json file or its .aptaswitch.json project.',
    'Ouvrir':
        'Open',
    'L’ouverture d’un autre run sera disponible à la fin du calcul en cours.':
        'You can open another run when the current calculation finishes.',
    'Paramètres enregistrés et console du run':
        'Saved parameters and run console',
    'paramètres':
        'parameters',
    'entrée':
        'input',
    'modèle':
        'model',
    'criblage':
        'screening',
    'Aucun journal enregistré à côté de ce fichier de résultats.':
        'No log was saved next to this results file.',
    'Structure 2D':
        '2D structure',
    'Séquence linéaire':
        'Linear sequence',
    'Télécharger la structure SVG':
        'Download structure SVG',
    'Télécharger la structure PNG':
        'Download structure PNG',
    'Séquence linéaire SVG':
        'Linear sequence SVG',
    'Séquence linéaire PNG':
        'Linear sequence PNG',
    'Exporter':
        'Export',
    'Afficher et copier la notation dot-bracket':
        'View and copy the dot-bracket notation',
    'Aperçu moléculaire':
        'Molecular preview',
    'Architecture du switch':
        'Switch architecture',
    'Toehold activateur':
        'Toehold activator',
    'Toehold répresseur':
        'Toehold repressor',
    'Toehold répresseur 3WJ':
        '3WJ toehold repressor',
    'Les architectures répressives ne sont pas encore disponibles.':
        'Repressor architectures are not available yet.',
    'Unité':
        'Unit',
    '5′ uniquement':
        '5′ only',
    '3′ uniquement':
        '3′ only',
    '**Trigger final · {length} nt**':
        '**Final trigger · {length} nt**',
    '{paired}/{total} bases du trigger appariées à l’aptamère · {conditions}':
        '{paired}/{total} trigger bases paired to the aptamer · {conditions}',
    'Les exports seront enregistrés dans {output}.':
        'Exports will be saved in {output}.',
    '{length} nt : +{prefix} bases en 5′ et +{suffix} bases en 3′ · {count} extensions à tester.':
        '{length} nt: +{prefix} bases at 5′ and +{suffix} bases at 3′ · {count} extensions to test.',
    '{name} · extensions':
        '{name} · extensions',
    'Candidat à comparer · {length} nt':
        'Candidate to compare · {length} nt',
    'Rang {rank} · {extension} · score {score:.4f}':
        'Rank {rank} · {extension} · score {score:.4f}',
    'Original · {length} nt':
        'Original · {length} nt',
    'Version {length} nt · rang {rank}':
        'Version {length} nt · rank {rank}',
    '#### Analyse des candidats · {length} nt':
        '#### Candidate analysis · {length} nt',
    'Candidat {length} nt · {extension}':
        'Candidate {length} nt · {extension}',
    'Copie locale : {path}':
        'Local copy: {path}',
    'Impossible d’ouvrir le dossier {path} : {error}':
        'Unable to open folder {path}: {error}',
    '{date} · {kind} · {name} · {count} candidat(s) · {run_id}':
        '{date} · {kind} · {name} · {count} candidate(s) · {run_id}',
    'Dossier des runs : {path}':
        'Runs folder: {path}',
    'Fichier sélectionné : {path}':
        'Selected file: {path}',
    'Run enregistré ouvert : {source} · résultats conservés, sans nouveau calcul.':
        'Opened saved run: {source} · original results, no new calculation.',
    'Extrémité choisie pour ce run : {side}':
        'Extension end selected for this run: {side}',
    'Le journal n’est plus accessible : {error}':
        'The log is no longer accessible: {error}',
    'Unité utilisée pour {label}.':
        'Unit used for {label}.',
    'Région reconnue du trigger : {trigger_length} nt · toehold : {toehold_length} nt · tige : {stem_length} nt':
        '{trigger_length} nt trigger region: {toehold_length} nt toehold + {stem_length} nt stem',
    'Ordre des brins : aptamère puis trigger, chacun de 5′ vers 3′. Comparez les paires et la notation ci-dessous : l’orientation du dessin 2D peut être différente. Les conditions du système sont distinctes des conditions du design de switch.':
        'Strand order: aptamer then trigger, each from 5′ to 3′. Compare the pairs and notation below: the 2D drawing orientation may differ. System conditions are separate from switch design conditions.',
    'Séquences FASTA et structure':
        'FASTA sequences and structure',
}

# Extension phase transitions and resumable MFE selection.
TRANSLATIONS.update({
    'La reprise exige les mêmes séquences, longueurs et conditions que le calcul enregistré.':
        'Resuming requires the same sequences, lengths and conditions as the saved calculation.',
    'La référence thermodynamique diffère du calcul enregistré ; reprise interrompue.': 'The thermodynamic reference differs from the saved calculation; resume stopped.',
    'Sauvegarde MFE incohérente : reprise interrompue.':
        'Inconsistent MFE checkpoint: resume stopped.',
    'Récupération des candidats MFE retenus · {done:,}/{total:,}.':
        'Recovering retained MFE candidates · {done:,}/{total:,}.',
    'Les candidats enregistrés n’ont pas pu être reconstitués avec ce modèle thermodynamique.': 'The saved candidates could not be reconstructed with this thermodynamic model.',
    'Criblage MFE arrêté avant la fin : le classement porte uniquement sur les candidats déjà examinés.':
        'MFE screening stopped early: ranking covers only the candidates already examined.',
    'Passage à l’analyse des candidats retenus après le calcul en cours…':
        'Moving to analysis of retained candidates after the current calculation…',
    'Passer à l’analyse des meilleurs candidats':
        'Analyze the best candidates so far',
    'Passer à la suite arrête le criblage MFE et analyse les candidats déjà retenus, dans la limite choisie. La recherche restera indiquée comme partielle.':
        'Moving on stops MFE screening and analyzes the candidates already retained, up to the selected limit. The search will remain marked as partial.',
    'Reprise avec les séquences et conditions enregistrées ; le criblage complet ne sera pas recommencé.':
        'Resuming with the saved sequences and conditions; full screening will not be repeated.',
    'Reprendre l’analyse des candidats retenus':
        'Resume analysis of retained candidates',
    'Les séquences et conditions d’origine seront réutilisées. L’analyse d’ensemble sera relancée sur le groupe MFE sauvegardé ; les anciens calculs sans sauvegarde nécessitent une récupération des premiers ex æquo.':
        'The original sequences and conditions will be reused. Ensemble analysis will restart on the saved MFE group; older runs without a checkpoint require recovery of the first tied candidates.',
    'Calculs interrompus · sauvegardes MFE':
        'Interrupted calculations · MFE checkpoints',
    'Sauvegarde à reprendre':
        'Checkpoint to resume',
    'Analyser les candidats de cette sauvegarde':
        'Analyze candidates from this checkpoint',
    'Soit {billions:.2f} milliards de possibilités. La limite d’analyse d’ensemble ne limite pas ce criblage.':
        'That is {billions:.2f} billion possibilities. The ensemble analysis limit does not limit this screening.',
})
TRANSLATIONS.update({
    'Espace de recherche : {count} possibilités · {percent:.4f} % examinées.':
        'Search space: {count} possibilities · {percent:.4f}% examined.',
})
TRANSLATIONS.update({
    'Reprise du groupe MFE retenu · vérification de la référence thermodynamique…': 'Resuming the retained MFE group · checking the thermodynamic reference…',
})

TRANSLATIONS.update({
    'La longueur cible de l’extension inverse doit dépasser celle du trigger original.':
        'The inverse extension target length must exceed the original trigger length.',
    'Le nombre de switches à comparer doit être un entier positif.':
        'The number of switches to compare must be a positive integer.',
    'La position de l’extension inverse est invalide.':
        'The inverse extension position is invalid.',
    'L’extension inverse conserve le trigger original entier ; une fenêtre partielle n’est pas compatible.':
        'Inverse extension preserves the entire original trigger; a partial binding window is not compatible.',
    'Le linker ne peut contenir que A, C, G, T, U ou X (base à optimiser).':
        'The linker may only contain A, C, G, T, U or X (a base to optimize).',
    'Le leader et le RBS doivent contenir uniquement des bases déterminées A, C, G, T ou U.':
        'The leader and RBS must contain only specified A, C, G, T or U bases.',
    'Le motif de conception ne peut contenir que A, C, G, T, U ou X.':
        'The design pattern may only contain A, C, G, T, U or X.',
    'Le trigger proposé ne respecte pas les bases fixes du motif d’extension.':
        'The proposed trigger does not preserve the fixed bases of the extension pattern.',
    'Le switch proposé n’est pas complémentaire du trigger étendu optimisé.':
        'The proposed switch is not complementary to the optimized extended trigger.',
    'Ce calcul ne contient pas d’extension conçue avec le switch.':
        'This calculation does not contain an extension designed with the switch.',
    'Comparaison des switches enregistrés : leurs séquences ne seront pas redessinées.':
        'Comparing saved switches: their sequences will not be redesigned.',
    'Comparaison des extensions': 'Extension comparison',
    'Conditions de comparaison': 'Comparison conditions',
    'Référence originale aptamère + trigger pour l’extension inverse…': 'Original aptamer + trigger reference for inverse extension…',
    'Comparaison de l’extension inverse {index}/{total}…': 'Comparing inverse extension {index}/{total}…',
    'Comparaison des extensions inverses terminée.': 'Inverse extension comparison complete.',
})
TRANSLATIONS.update({
})

TRANSLATIONS.update({
    'Le modèle thermodynamique diffère du calcul enregistré ; reprise interrompue.':
        'The thermodynamic model differs from the saved calculation; resume stopped.',
})

TRANSLATIONS.update({
    'Annulation demandée — fin du calcul en cours…':
        'Cancellation requested — waiting for the current calculation to finish…',
})
