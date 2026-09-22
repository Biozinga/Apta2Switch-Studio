"""English translations for reusable interface components."""

TRANSLATIONS = {
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
    'N indique une base laissée libre pour l’optimisation NUPACK.':
        'N marks a base left free for NUPACK to optimize.',
    'Aperçu de l’extension SVG':
        'Extension preview SVG',
    'Aperçu de l’extension PNG':
        'Extension preview PNG',
    '#### Conditions du système':
        '#### System conditions',
    'Modèle thermodynamique':
        'Thermodynamic model',
    'ADN · dna04':
        'DNA · dna04',
    'ARN · rna06':
        'RNA · rna06',
    'ADN : U est converti en T. ARN : T est converti en U, pour les deux brins.':
        'DNA: U is converted to T. RNA: T is converted to U, for both strands.',
    'Température du système (°C)':
        'System temperature (°C)',
    'Choisissez des conditions aussi proches que possible des conditions expérimentales ayant démontré l’efficacité du système aptamère–trigger pour votre cible. Elles seront utilisées pour la structure initiale et l’extension éventuelle.':
        'Choose conditions as close as possible to the experimental conditions that demonstrated the effectiveness of the aptamer–trigger system for your target. They will be used for the initial structure and any extension.',
    '#### Aperçus du système initial':
        '#### Initial system previews',
    'La structure initiale sera actualisée à la fin du calcul en cours.':
        'The initial structure will update when the current calculation finishes.',
    'Calcul de la structure initiale avec NUPACK…':
        'Calculating the initial structure with NUPACK…',
    'Saisissez l’aptamère et le trigger pour afficher leur structure 2D.':
        'Enter the aptamer and trigger sequences to display their 2D structure.',
    'NUPACK est requis pour la structure 2D. Configurez-le dans Réglages → Dépendances.':
        'NUPACK is required for the 2D structure. Configure it in Settings → Dependencies.',
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
    'Structure MFE calculée par NUPACK, ensemble stacking. Rose : bases de l’aptamère annotées pour le ligand. Le ligand lui-même n’est pas simulé.':
        'MFE structure calculated by NUPACK using the stacking ensemble. Pink: aptamer bases annotated for ligand binding. The ligand itself is not simulated.',
    'Séquences du système initial':
        'Initial system sequences',
    'Séquence initiale SVG':
        'Initial sequence SVG',
    'Séquence initiale PNG':
        'Initial sequence PNG',
    'Un calcul est déjà en cours. Attendez sa fin ou annulez-le dans Résultats.':
        'A calculation is already running. Wait for it to finish or cancel it in Results.',
    'aptamère':
        'aptamer',
    'L’extension est désactivée pour les triggers de 24 à 30 nt. Vous pouvez passer à la conception du switch.':
        'Extension is disabled for triggers between 24 and 30 nt. You can proceed to switch design.',
    'Choisissez une longueur finale supérieure à celle du trigger original.':
        'Choose a final length greater than the original trigger length.',
    "La longueur cible est limitée à 30 nt. Activez l'option de dépassement pour aller au-delà.":
        'The target length is limited to 30 nt. Enable longer extensions to go beyond this limit.',
    'Configurez NUPACK dans Réglages → Dépendances pour calculer les extensions.':
        'Configure NUPACK in Settings → Dependencies to calculate extensions.',
    '2 · Extension du trigger (facultatif)':
        '2 · Trigger extension (optional)',
    'Conservez votre trigger tel quel ou choisissez une extension. Le calcul utilisera les séquences, les annotations et les conditions de la section 1.':
        'Keep your original trigger or choose an extension. The calculation uses the sequences, annotations and conditions from section 1.',
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
    'NUPACK est requis : configurez-le dans Réglages → Dépendances.':
        'NUPACK is required: configure it in Settings → Dependencies.',
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
    'Suivi du criblage MFE, puis des probabilités d’ensemble pour chaque longueur. L’estimation porte sur l’étape en cours.':
        'Monitoring MFE screening, then ensemble probabilities for each length. The estimate applies to the current stage.',
    'Le calcul a échoué.':
        'The calculation failed.',
    'Extension du trigger · conservation du complexe aptamère–trigger · MFE et probabilités d’ensemble NUPACK':
        'Trigger extension · preserving the aptamer–trigger complex · NUPACK MFE and ensemble probabilities',
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
    'Moitié 5′ / moitié 3′':
        'Half 5′ / half 3′',
    '**Trigger final · {length} nt**':
        '**Final trigger · {length} nt**',
    '{paired}/{total} bases du trigger appariées à l’aptamère · {conditions}':
        '{paired}/{total} trigger bases paired to the aptamer · {conditions}',
    'Les exports seront enregistrés dans {output}.':
        'Exports will be saved in {output}.',
    'Votre trigger mesure {length} nt. L’extension est désactivée entre 24 et 30 nt inclus ; vous pouvez passer à la conception du switch.':
        'Your trigger is {length} nt long. Extension is disabled from 24 to 30 nt inclusive; you can proceed to switch design.',
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
    'Équivalent : {value:.6g} M':
        'Equivalent: {value:.6g} M',
    'Région reconnue du trigger : {trigger_length} nt · toehold : {toehold_length} nt · tige : {stem_length} nt':
        '{trigger_length} nt trigger region: {toehold_length} nt toehold + {stem_length} nt stem',
}
