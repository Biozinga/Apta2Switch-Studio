"""English translations for computational progress, validation and figures."""

TRANSLATIONS = {'Les positions du ligand doivent appartenir à l’aptamère (numérotation à partir de 1).': 'Ligand positions '
                                                                                          'must lie within '
                                                                                          'the aptamer '
                                                                                          '(numbered from '
                                                                                          '1).',
 'Choisissez une extension en 5′, en 3′ ou moitié de chaque côté.': 'Choose a 5′ extension, a 3′ extension, '
                                                                    'or half on each side.',
 'La longueur d’extension doit être un entier positif ou nul.': 'Extension length must be a non-negative '
                                                                'integer.',
 'La longueur du trigger original doit être un entier positif.': 'Original trigger length must be a positive '
                                                                 'integer.',
 'Choisissez au moins une longueur finale de trigger.': 'Choose at least one final trigger length.',
 'Les longueurs finales doivent être des entiers.': 'Final lengths must be integers.',
 'Chaque longueur finale doit dépasser celle du trigger original.': 'Each final length must exceed the '
                                                                    'original trigger length.',
 'NUPACK est obligatoire pour calculer une extension réelle. Configurez son chemin.': 'NUPACK is required to '
                                                                                      'calculate an '
                                                                                      'extension. Configure '
                                                                                      'its path.',
 "Le nombre exporté et la limite de candidats d'ensemble doivent être des entiers positifs.": 'The export '
                                                                                              'count and '
                                                                                              'ensemble '
                                                                                              'candidate '
                                                                                              'limit must be '
                                                                                              'positive '
                                                                                              'integers.',
 'La chimie du trigger doit être DNA ou RNA.': 'Trigger chemistry must be DNA or RNA.',
 'Les conditions thermodynamiques doivent être finies.': 'Thermodynamic conditions must be finite.',
 'Les concentrations ioniques doivent être positives ou nulles.': 'Ion concentrations must be non-negative.',
 'Au-delà de 30 nt, le calcul peut prendre beaucoup de temps ; chaque base ajoutée multiplie la recherche par quatre.': 'Above '
                                                                                                                        '30 '
                                                                                                                        'nt, '
                                                                                                                        'calculations '
                                                                                                                        'can '
                                                                                                                        'take '
                                                                                                                        'a '
                                                                                                                        'long '
                                                                                                                        'time; '
                                                                                                                        'each '
                                                                                                                        'additional '
                                                                                                                        'base '
                                                                                                                        'multiplies '
                                                                                                                        'the '
                                                                                                                        'search '
                                                                                                                        'space '
                                                                                                                        'by '
                                                                                                                        'four.',
 'NUPACK · référence aptamère + trigger original : MFE, ensemble et probabilités de paires…': 'NUPACK · '
                                                                                              'aptamer + '
                                                                                              'original '
                                                                                              'trigger '
                                                                                              'reference: '
                                                                                              'MFE, ensemble '
                                                                                              'and base-pair '
                                                                                              'probabilities…',
 "La structure MFE de référence ne présente aucune paire entre l'aptamère et le trigger ; la conservation de cette référence ne démontre pas leur liaison.": 'The '
                                                                                                                                                             'reference '
                                                                                                                                                             'MFE '
                                                                                                                                                             'structure '
                                                                                                                                                             'has '
                                                                                                                                                             'no '
                                                                                                                                                             'base '
                                                                                                                                                             'pairs '
                                                                                                                                                             'between '
                                                                                                                                                             'the '
                                                                                                                                                             'aptamer '
                                                                                                                                                             'and '
                                                                                                                                                             'trigger; '
                                                                                                                                                             'preserving '
                                                                                                                                                             'this '
                                                                                                                                                             'reference '
                                                                                                                                                             'does '
                                                                                                                                                             'not '
                                                                                                                                                             'demonstrate '
                                                                                                                                                             'their '
                                                                                                                                                             'binding.',
 'Calcul annulé : les résultats éventuellement présents sont partiels.': 'Calculation cancelled: any '
                                                                         'available results are partial.',
 "NUPACK ou une de ses dépendances est indisponible. Aucune extension n'a été simulée.": 'NUPACK or one of '
                                                                                         'its dependencies '
                                                                                         'is unavailable. No '
                                                                                         'extension was '
                                                                                         'simulated.',
 'annulée': 'cancelled',
 'terminée': 'completed',
 'Criblage RBS–linker : bornes de région invalides.': 'RBS–linker screening: invalid region boundaries.',
 'Criblage RBS–linker : aucune structure MFE ON renvoyée.': 'RBS–linker screening: no ON MFE structure '
                                                            'returned.',
 'Criblage RBS–linker : les longueurs de la structure MFE ON ne correspondent pas aux brins.': 'RBS–linker '
                                                                                               'screening: '
                                                                                               'ON MFE '
                                                                                               'structure '
                                                                                               'lengths do '
                                                                                               'not match '
                                                                                               'the strands.',
 'Criblage RBS–linker : structure MFE ON mal formée.': 'RBS–linker screening: malformed ON MFE structure.',
 'NUPACK · préparation du modèle et des tubes…': 'NUPACK · preparing the model and tubes…',
 'NUPACK n’a renvoyé aucune structure MFE ON.': 'NUPACK returned no ON MFE structure.',
 'criblage RBS–linker dans les structures MFE ON': 'screening RBS–linker in ON MFE structures',
 'analyses thermodynamiques': 'thermodynamic analyses',
 'Configurez NUPACK pour calculer la structure du système initial.': 'Configure NUPACK to calculate the '
                                                                     'initial system structure.',
 'Le système initial nécessite un modèle ADN dna04 ou ARN rna06.': 'The initial system requires a DNA dna04 '
                                                                   'or RNA rna06 model.',
 'NUPACK n’a retourné aucune structure MFE pour le système initial.': 'NUPACK returned no MFE structure for '
                                                                      'the initial system.',
 'NUPACK a retourné une structure incompatible avec les séquences fournies.': 'NUPACK returned a structure '
                                                                              'incompatible with the '
                                                                              'supplied sequences.',
 'NUPACK a retourné une énergie MFE non finie.': 'NUPACK returned a non-finite MFE energy.',
 'Structure calculée absente ou invalide : aucun schéma ne peut être affiché.': 'Calculated structure '
                                                                                'missing or invalid: no '
                                                                                'diagram can be displayed.',
 'Structure calculée invalide : aucun schéma de remplacement n’est généré.': 'Invalid calculated structure: '
                                                                             'no replacement diagram is '
                                                                             'generated.',
 'Les séquences calculées du candidat sont absentes ou invalides.': 'The calculated candidate sequences are '
                                                                    'missing or invalid.',
 'La séquence calculée ne correspond pas à l’architecture du switch.': 'The calculated sequence does not '
                                                                       'match the switch architecture.',
 'La structure calculée ne correspond pas aux longueurs des séquences du candidat.': 'The calculated '
                                                                                     'structure does not '
                                                                                     'match the candidate '
                                                                                     'sequence lengths.',
 'Aucune séquence à afficher': 'No sequence to display',
 "Le trigger s'approche du toehold...": 'The trigger approaches the toehold...',
 "Le trigger s'accroche au toehold.": 'The trigger binds to the toehold.',
 "Le stem s'ouvre progressivement...": 'The stem opens progressively...',
 'Switch ON : trigger + switch appariés.': 'Switch ON: trigger + switch paired.',
 '  -  ⚠ codon STOP prématuré (cercle rouge)': '  -  ⚠ premature STOP codon (red circle)',
 'Seul le moteur NUPACK est pris en charge. Les runs à blanc ont été supprimés.': 'Only the NUPACK engine is '
                                                                                  'supported. Dry runs have '
                                                                                  'been removed.',
 'Écriture des exports…': 'Writing exports…',
 'Exports terminés.': 'Exports complete.',
 'Bases liées au ligand': 'Ligand-binding bases',
 'P(paire MFE), ou P(non appariée) pour une base libre.': 'P(MFE pair), or P(unpaired) for an unpaired base.',
 'Cercles roses : bases de liaison au ligand (annotation utilisateur).': 'Pink circles: ligand-binding bases '
                                                                         '(user annotation).',
 'Chaque point représente un candidat calculé. La couleur représente son score final ; un score plus bas est meilleur.': 'Each '
                                                                                                                         'point '
                                                                                                                         'represents '
                                                                                                                         'a '
                                                                                                                         'calculated '
                                                                                                                         'candidate. '
                                                                                                                         'Color '
                                                                                                                         'indicates '
                                                                                                                         'its '
                                                                                                                         'final '
                                                                                                                         'score; '
                                                                                                                         'lower '
                                                                                                                         'is '
                                                                                                                         'better.',
 'Aucun candidat évalué': 'No evaluated candidates',
 'Meilleur score': 'Best score',
 'Score final': 'Final score',
 'Plus bas = meilleur': 'Lower = better',
 'Pénalité d’isolation de l’extension': 'Extension isolation penalty',
 'Σ P(extension appariée) + Σ P(extension ↔ aptamère)': 'Σ P(paired extension) + Σ P(extension ↔ aptamer)',
 'Σ perte de probabilité des paires de référence': 'Σ reference base-pair probability loss',
 'Axes linéaires ajustés aux valeurs des candidats affichés · Losange : meilleur score.': 'Linear axes '
                                                                                          'fitted to the '
                                                                                          'displayed '
                                                                                          'candidates · '
                                                                                          'Diamond: best '
                                                                                          'score.',
 'Valeurs issues de l’analyse d’ensemble NUPACK ; les points ne sont pas décalés artificiellement.': 'Values '
                                                                                                     'from '
                                                                                                     'NUPACK '
                                                                                                     'ensemble '
                                                                                                     'analysis; '
                                                                                                     'points '
                                                                                                     'are '
                                                                                                     'not '
                                                                                                     'artificially '
                                                                                                     'shifted.',
 'Paires de la structure MFE originale : probabilités de référence, probabilités du candidat, et différences signées. A désigne l’aptamère et T le trigger initial.': 'Base '
                                                                                                                                                                      'pairs '
                                                                                                                                                                      'from '
                                                                                                                                                                      'the '
                                                                                                                                                                      'original '
                                                                                                                                                                      'MFE '
                                                                                                                                                                      'structure: '
                                                                                                                                                                      'reference '
                                                                                                                                                                      'and '
                                                                                                                                                                      'candidate '
                                                                                                                                                                      'probabilities, '
                                                                                                                                                                      'and '
                                                                                                                                                                      'signed '
                                                                                                                                                                      'differences. '
                                                                                                                                                                      'A '
                                                                                                                                                                      'denotes '
                                                                                                                                                                      'the '
                                                                                                                                                                      'aptamer '
                                                                                                                                                                      'and '
                                                                                                                                                                      'T '
                                                                                                                                                                      'the '
                                                                                                                                                                      'original '
                                                                                                                                                                      'trigger.',
 'La structure de référence ne contient aucune paire': 'The reference structure contains no base pairs',
 'Probabilité d’appariement d’ensemble': 'Ensemble base-pair probability',
 'Δ (candidat − référence)': 'Δ (candidate − reference)',
 'Référence originale': 'Original reference',
 'A = aptamère ; T = trigger initial. Positions locales à partir de 1 ; paires intrabrin : positions globales.': 'A '
                                                                                                                 '= '
                                                                                                                 'aptamer; '
                                                                                                                 'T '
                                                                                                                 '= '
                                                                                                                 'original '
                                                                                                                 'trigger. '
                                                                                                                 'Local '
                                                                                                                 'positions '
                                                                                                                 'start '
                                                                                                                 'at '
                                                                                                                 '1; '
                                                                                                                 'intrastrand '
                                                                                                                 'pairs '
                                                                                                                 'use '
                                                                                                                 'global '
                                                                                                                 'positions.',
 'Structures 2D : référence et extensions retenues': '2D structures: reference and selected extensions',
 'Structures MFE NUPACK : complexe original puis meilleur score final pour chaque longueur cible. Chaque panneau est ajusté à sa propre échelle spatiale.': 'NUPACK '
                                                                                                                                                            'MFE '
                                                                                                                                                            'structures: '
                                                                                                                                                            'original '
                                                                                                                                                            'complex '
                                                                                                                                                            'followed '
                                                                                                                                                            'by '
                                                                                                                                                            'the '
                                                                                                                                                            'best '
                                                                                                                                                            'final '
                                                                                                                                                            'score '
                                                                                                                                                            'for '
                                                                                                                                                            'each '
                                                                                                                                                            'target '
                                                                                                                                                            'length. '
                                                                                                                                                            'Each '
                                                                                                                                                            'panel '
                                                                                                                                                            'has '
                                                                                                                                                            'its '
                                                                                                                                                            'own '
                                                                                                                                                            'spatial '
                                                                                                                                                            'scale.',
 'Complexe original et meilleur score par longueur · Chaque structure est ajustée à son panneau.': 'Original '
                                                                                                   'complex '
                                                                                                   'and best '
                                                                                                   'score '
                                                                                                   'per '
                                                                                                   'length · '
                                                                                                   'Each '
                                                                                                   'structure '
                                                                                                   'is '
                                                                                                   'fitted '
                                                                                                   'to its '
                                                                                                   'panel.',
 ' · liaison au ligand': ' · ligand binding',
 'Les longueurs du toehold et du stem doivent être positives.': 'Toehold and stem lengths must be positive.',
 'La longueur du stem supérieur doit être positive.': 'Upper stem length must be positive.',
 'La longueur du stem supérieur facultatif ne peut pas être négative.': 'Optional upper stem length cannot '
                                                                        'be negative.',
 'La longueur du bulge doit être positive.': 'Bulge length must be positive.',
 'Les longueurs du préfixe RBS et de l’espace avant AUG ne peuvent pas être négatives.': 'RBS prefix and AUG '
                                                                                         'spacer lengths '
                                                                                         'cannot be '
                                                                                         'negative.',
 'Le nombre d’essais doit être positif.': 'Trials must be positive.',
 'La concentration du trigger doit être positive.': 'Trigger concentration must be positive.',
 'La concentration du toehold switch doit être positive.': 'Toehold-switch concentration must be positive.',
 'NUPACK n’est pas configuré. Sélectionnez d’abord votre chemin NUPACK.': 'NUPACK is not configured. Select '
                                                                          'a user-provided NUPACK path '
                                                                          'first.',
 'Une structure de référence à deux brins est attendue': 'Expected a two-strand reference structure',
 'NUPACK a renvoyé une structure dont les longueurs de brins sont inattendues.': 'NUPACK returned a '
                                                                                 'structure with unexpected '
                                                                                 'strand lengths.',
 'NUPACK n’a renvoyé aucune structure MFE.': 'NUPACK returned no MFE structure.',
 'Une paire de bases du schéma linéaire fait référence à un nucléotide absent': 'A linear base pair '
                                                                                'references a missing '
                                                                                'nucleotide',
 'La longueur de la séquence du switch ne correspond pas à l’architecture': 'Switch sequence length does not '
                                                                            'match the architecture',
 'Un switch calculé ne doit contenir que des bases nucléotidiques déterminées': 'A completed switch must '
                                                                                'contain only resolved '
                                                                                'nucleotide bases',
 'Une séquence nucléotidique complète est nécessaire': 'A complete nucleotide sequence is required',
 'Longueur totale d’extension invalide': 'Invalid total extension length',
 'Longueur d’extension 3′ invalide': 'Invalid 3-prime extension length',
 'Les probabilités doivent être comprises entre 0 et 1': 'Probabilities must lie between 0 and 1',
 'Les séquences complètes de l’aptamère et du trigger sont nécessaires': 'A complete aptamer and trigger '
                                                                         'sequence is required',
 'Une structure dot-bracket aptamère + trigger est attendue': 'Expected an aptamer + trigger dot-bracket '
                                                              'structure',
 'Les longueurs de brins de la structure ne correspondent pas aux séquences': 'Structure strand lengths do '
                                                                              'not match the sequences',
 'Le canevas de structure doit mesurer au moins 320 × 320': 'Structure canvas must be at least 320 by 320',
 'Une probabilité par nucléotide est attendue': 'Expected one probability per nucleotide',
 'Dimensions du graphique ou sélection des candidats invalides': 'Invalid landscape dimensions or candidate '
                                                                 'selection',
 'Dimensions du graphique de probabilités invalides': 'Invalid probability figure dimensions',
 'Une structure de référence calculée est nécessaire': 'A measured reference structure is required',
 'Impossible d’importer NUPACK.': 'Unable to import nupack.',
 'Ce chemin n’existe pas.': 'Path does not exist.',
 'Sélectionnez un dossier, un dossier .local_deps ou un fichier .whl.': 'Select a folder, .local_deps '
                                                                        'folder, or .whl file.',
 'Aucun paquet NUPACK importable trouvé dans ce dossier.': 'No importable nupack package found in this '
                                                           'folder.',
 'Un fichier .whl est attendu.': 'Expected a .whl file.',
 'Le nom du fichier wheel ne semble pas correspondre à NUPACK.': 'The wheel filename does not look like '
                                                                 'NUPACK.',
 'L’exécutable Python n’est pas disponible.': 'Python executable not available.',
 'L’installation de validation du fichier wheel a échoué.': 'Wheel validation install failed.',
 'Le nom du rapporteur ne peut pas être vide.': 'Reporter name cannot be empty.',
 'Les rapporteurs intégrés ne peuvent pas être supprimés.': 'Built-in reporters cannot be removed.',
 'Candidats AptaSwitch': 'AptaSwitch candidates',
 'Avertissements :': 'Warnings:',
 '{v0} est requis.': '{v0} is required.',
 '{v0} contient des bases non prises en charge : {v1}.': '{v0} contains unsupported bases: {v1}.',
 'La longueur toehold + stem doit correspondre à celle du trigger ({v0} != {v1}).': 'Toehold + stem length '
                                                                                    'must match the trigger '
                                                                                    'length ({v0} != {v1}).',
 'NUPACK · essai {v0}/{v1} : optimisation…': 'NUPACK · trial {v0}/{v1}: optimizing…',
 'NUPACK · essai {v0}/{v1} : {v2}…': 'NUPACK · trial {v0}/{v1}: {v2}…',
 'NUPACK · essai {v0}/{v1} terminé.': 'NUPACK · trial {v0}/{v1} complete.',
 'NUPACK · essai {v0}/{v1} exclu : STOP prématuré au codon {v2}. MFE et analyses en tube ignorées.': 'NUPACK '
                                                                                                     '· '
                                                                                                     'trial '
                                                                                                     '{v0}/{v1} '
                                                                                                     'excluded: '
                                                                                                     'premature '
                                                                                                     'STOP '
                                                                                                     'at '
                                                                                                     'codon '
                                                                                                     '{v2}. '
                                                                                                     'MFE '
                                                                                                     'and '
                                                                                                     'tube '
                                                                                                     'analyses '
                                                                                                     'skipped.',
 'NUPACK · essai {v0}/{v1} : criblage RBS–linker validé, analyses thermodynamiques…': 'NUPACK · trial '
                                                                                      '{v0}/{v1}: RBS–linker '
                                                                                      'screening passed, '
                                                                                      'thermodynamic '
                                                                                      'analyses…',
 'NUPACK · essai {v0}/{v1} exclu : RBS–linker non linéaire en ON. Autres MFE et analyses en tube ignorées.': 'NUPACK '
                                                                                                             '· '
                                                                                                             'trial '
                                                                                                             '{v0}/{v1} '
                                                                                                             'excluded: '
                                                                                                             'RBS–linker '
                                                                                                             'not '
                                                                                                             'linear '
                                                                                                             'in '
                                                                                                             'ON. '
                                                                                                             'Other '
                                                                                                             'MFE '
                                                                                                             'and '
                                                                                                             'tube '
                                                                                                             'analyses '
                                                                                                             'skipped.',
 'NUPACK · extension {v0} : {v1:,} extensions à examiner, recherche exhaustive.': 'NUPACK · {v0} extension: '
                                                                                  '{v1:,} extensions to '
                                                                                  'examine, exhaustive '
                                                                                  'search.',
 'Structure dot-bracket déséquilibrée : {v0}': 'Unbalanced dot-bracket structure: {v0}',
 'Structure déséquilibrée : {v0}': 'Unbalanced structure: {v0}',
 'Caractère de structure non pris en charge : {v0!r}': 'Unsupported structure character: {v0!r}',
 'Référence calculée : MFE {v0:.6f} kcal/mol ; début du criblage MFE.': 'Reference calculated: MFE {v0:.6f} '
                                                                        'kcal/mol; starting MFE screening.',
 'Extension {v0} : {v1:,}/{v2:,} MFE, {v3:,} ensembles analysés, {v4:,} candidats exportés.': 'Extension '
                                                                                              '{v0}: '
                                                                                              '{v1:,}/{v2:,} '
                                                                                              'MFE, {v3:,} '
                                                                                              'ensembles '
                                                                                              'analyzed, '
                                                                                              '{v4:,} '
                                                                                              'candidates '
                                                                                              'exported.',
 'Cible {v0} nt · ajout de {v1} bases en 5′ et {v2} bases en 3′.': 'Target {v0} nt · adding {v1} bases at 5′ '
                                                                   'and {v2} bases at 3′.',
 'Cible {v0} nt : aucune extension examinée ne conserve exactement la structure MFE de référence avec une extension non appariée (meilleure clé {v1}). Les candidats présentés sont les meilleurs compromis calculés.': 'Target '
                                                                                                                                                                                                                        '{v0} '
                                                                                                                                                                                                                        'nt: '
                                                                                                                                                                                                                        'no '
                                                                                                                                                                                                                        'examined '
                                                                                                                                                                                                                        'extension '
                                                                                                                                                                                                                        'exactly '
                                                                                                                                                                                                                        'preserves '
                                                                                                                                                                                                                        'the '
                                                                                                                                                                                                                        'reference '
                                                                                                                                                                                                                        'MFE '
                                                                                                                                                                                                                        'structure '
                                                                                                                                                                                                                        'with '
                                                                                                                                                                                                                        'an '
                                                                                                                                                                                                                        'unpaired '
                                                                                                                                                                                                                        'extension '
                                                                                                                                                                                                                        '(best '
                                                                                                                                                                                                                        'key '
                                                                                                                                                                                                                        '{v1}). '
                                                                                                                                                                                                                        'The '
                                                                                                                                                                                                                        'displayed '
                                                                                                                                                                                                                        'candidates '
                                                                                                                                                                                                                        'are '
                                                                                                                                                                                                                        'the '
                                                                                                                                                                                                                        'best '
                                                                                                                                                                                                                        'calculated '
                                                                                                                                                                                                                        'compromises.',
 "Cible {v0} nt : {v1:,} ex æquo MFE ; seuls les {v2:,} premiers en ordre alphabétique passent à l'ensemble. Le classement d'ensemble est partiel.": 'Target '
                                                                                                                                                     '{v0} '
                                                                                                                                                     'nt: '
                                                                                                                                                     '{v1:,} '
                                                                                                                                                     'MFE '
                                                                                                                                                     'ties; '
                                                                                                                                                     'only '
                                                                                                                                                     'the '
                                                                                                                                                     'first '
                                                                                                                                                     '{v2:,} '
                                                                                                                                                     'in '
                                                                                                                                                     'alphabetical '
                                                                                                                                                     'order '
                                                                                                                                                     'proceed '
                                                                                                                                                     'to '
                                                                                                                                                     'ensemble '
                                                                                                                                                     'analysis. '
                                                                                                                                                     'Ensemble '
                                                                                                                                                     'ranking '
                                                                                                                                                     'is '
                                                                                                                                                     'partial.',
 'Ensemble NUPACK · candidat {v0:,}/{v1:,} ; 5′ {v2} ; 3′ {v3}…': 'NUPACK ensemble · candidate '
                                                                  '{v0:,}/{v1:,}; 5′ {v2}; 3′ {v3}…',
 'MFE · cible {v0} nt : {v1:,}/{v2:,} ; meilleur groupe {v3:,}, clé {v4}.': 'MFE · target {v0} nt: '
                                                                            '{v1:,}/{v2:,}; best group '
                                                                            '{v3:,}, key {v4}.',
 'Trigger + Switch ARN ({v0})': 'Trigger + RNA switch ({v0})',
 'Switch ARN ({v0})': 'RNA switch ({v0})',
 'Switch ARN ({v0}) - N = base a optimiser (NUPACK)': 'RNA switch ({v0}) - N = base to optimize (NUPACK)',
 '{v0} — {v1} / {v2} candidats': '{v0} — {v1} / {v2} candidates',
 'Probabilités des {v0} paires de référence': 'Probabilities of the {v0} reference base pairs',
 '{v0:.6g} + graduation {v1}': '{v0:.6g} + tick {v1}',
 'Aptamère : {v0} nt · Trigger : {v1} nt · Ajouts : 5′ {v2} nt / 3′ {v3} nt': 'Aptamer: {v0} nt · Trigger: '
                                                                              '{v1} nt · Additions: 5′ {v2} '
                                                                              'nt / 3′ {v3} nt',
 ' (meilleurs {v0:.0%})': ' (best {v0:.0%})',
 'Rang {v0} · {v1} nt · 5′ {v2} · 3′ {v3} · Score {v4:.8g} · Isolation {v5:.8g} · Perte {v6:.8g}': 'Rank '
                                                                                                   '{v0} · '
                                                                                                   '{v1} nt '
                                                                                                   '· 5′ '
                                                                                                   '{v2} · '
                                                                                                   '3′ {v3} '
                                                                                                   '· Score '
                                                                                                   '{v4:.8g} '
                                                                                                   '· '
                                                                                                   'Isolation '
                                                                                                   '{v5:.8g} '
                                                                                                   '· Loss '
                                                                                                   '{v6:.8g}',
 '{v0} candidat(s) ont une perte exactement nulle ; leurs points sont alignés sur zéro.': '{v0} candidate(s) '
                                                                                          'have exactly zero '
                                                                                          'loss; their '
                                                                                          'points align at '
                                                                                          'zero.',
 'Valeur non finie pour {v0}': 'Non-finite value for {v0}',
 'Ajouts : 5′ {v0} nt / 3′ {v1} nt · Structure NUPACK.': 'Additions: 5′ {v0} nt / 3′ {v1} nt · NUPACK '
                                                         'structure.',
 'Structure MFE NUPACK : {v0}': 'NUPACK MFE structure: {v0}',
 'Référence · {v0}': 'Reference · {v0}',
 '|Δ| maximal = {v0:.4g}': 'Maximum |Δ| = {v0:.4g}',
 'Référence originale — {v0} nt': 'Original reference — {v0} nt',
 'Trigger étendu — {v0} nt': 'Extended trigger — {v0} nt',
 'Paysage de sélection · {v0} nt': 'Selection landscape · {v0} nt',
 'Trigger étendu · {v0} nt': 'Extended trigger · {v0} nt',
 'Extension · trigger final de {v0} nt': 'Extension · final trigger of {v0} nt',
 '5′ : +{v0} nt · 3′ : +{v1} nt · N : base à optimiser avec NUPACK · Rose : liaison au ligand': '5′: +{v0} '
                                                                                                'nt · 3′: '
                                                                                                '+{v1} nt · '
                                                                                                'N: base to '
                                                                                                'optimize '
                                                                                                'with NUPACK '
                                                                                                '· Pink: '
                                                                                                'ligand '
                                                                                                'binding',
 'Aptamère · {v0}': 'Aptamer · {v0}',
 '{v0} · {v1} · à optimiser': '{v0} · {v1} · to optimize',
 'Domaine de {v0} bases à optimiser · affichage abrégé': 'Domain of {v0} bases to optimize · abbreviated '
                                                         'display',
 'Import NUPACK validé depuis {v0}.': 'NUPACK import validated from {v0}.',
 'Fichier wheel NUPACK validé : {v0}.': 'NUPACK wheel validated: {v0}.',
 'Identifiant du run : {v0}': 'Run ID: {v0}',
 'Statut : {v0}': 'Status: {v0}',
 'Début : {v0}': 'Started: {v0}',
 'Fin : {v0}': 'Completed: {v0}',
 'Moteur : {v0}': 'Engine: {v0}',
 'Essais demandés : {v0}': 'Trials requested: {v0}',
 'Candidats : {v0}': 'Candidates: {v0}',
 'Exclure les candidats avec un STOP prématuré : {v0}': 'Exclude premature STOP candidates: {v0}',
 'Exclure les candidats au RBS–linker non linéaire : {v0}': 'Exclude non-linear RBS-linker candidates: {v0}',
 'Candidats exclus (analyses restantes ignorées) : {v0}': 'Excluded candidates (remaining analyses skipped): '
                                                          '{v0}',
 'Concentration du trigger (M) : {v0}': 'Trigger concentration (M): {v0}',
 'Concentration du toehold switch (M) : {v0}': 'Toehold-switch concentration (M): {v0}',
 'Aptamère': 'Aptamer',
 'aptamère': 'aptamer',
 'séquence': 'sequence',
 'séquence du rapporteur': 'reporter sequence',
 'Complexe aptamère + trigger': 'Aptamer + trigger complex',
 'Trigger initial': 'Original trigger',
 'Trigger étendu': 'Extended trigger',
 'Boucle RBS': 'RBS loop',
 'Switch ARN': 'RNA switch',
 'Paysage de sélection': 'Selection landscape',
 'Moitié 5′ / moitié 3′': 'Half 5′ / half 3′',
 'Fixation du trigger': 'Trigger binding',
 'Codon STOP prématuré — codon {first_stop}': 'Premature STOP codon — codon {first_stop}',
 'RBS–linker non linéaire en ON — au moins une base appariée dans une structure MFE renvoyée': 'RBS–linker '
                                                                                               'not linear '
                                                                                               'in ON — at '
                                                                                               'least one '
                                                                                               'paired base '
                                                                                               'in a '
                                                                                               'returned MFE '
                                                                                               'structure',
 'L’exclusion des STOP prématurés est désactivée : ces candidats seront analysés et pourront être bien classés, même si leur traduction s’arrête avant le reporter.': 'Premature '
                                                                                                                                                                      'STOP '
                                                                                                                                                                      'exclusion '
                                                                                                                                                                      'is '
                                                                                                                                                                      'disabled: '
                                                                                                                                                                      'these '
                                                                                                                                                                      'candidates '
                                                                                                                                                                      'will '
                                                                                                                                                                      'be '
                                                                                                                                                                      'analyzed '
                                                                                                                                                                      'and '
                                                                                                                                                                      'may '
                                                                                                                                                                      'rank '
                                                                                                                                                                      'highly, '
                                                                                                                                                                      'even '
                                                                                                                                                                      'if '
                                                                                                                                                                      'their '
                                                                                                                                                                      'translation '
                                                                                                                                                                      'stops '
                                                                                                                                                                      'before '
                                                                                                                                                                      'the '
                                                                                                                                                                      'reporter.',
 'candidat': 'candidate',
 'référence': 'reference',
 'Conception Apta2Switch-Studio la mieux classée': 'Apta2Switch-Studio top design',
 'Trigger ADN 5′ vers 3′': 'Trigger DNA 5′ to 3′',
 'Toehold switch ARN 5′ vers 3′': 'Toehold switch RNA 5′ to 3′',
 'Ce SVG est un aperçu compact du résultat, pas une simulation de structure.': 'This SVG is a compact result '
                                                                               'preview, not a structural '
                                                                               'simulation.',
 'Score {score:.4f} | ON {on:.2f}% | fuite {leak:.2f}% | ΔG RBS–linker {dg:.2f} kcal/mol': 'Score '
                                                                                           '{score:.4f} | ON '
                                                                                           '{on:.2f}% | leak '
                                                                                           '{leak:.2f}% | ΔG '
                                                                                           'RBS–linker '
                                                                                           '{dg:.2f} '
                                                                                           'kcal/mol',
 'Type de run inconnu : {kind}': 'Unknown run type: {kind}',
 'terminé': 'completed',
 'annulé': 'cancelled',
 'en cours': 'running',
 'La conversion de la structure en PNG a échoué.': 'Converting the structure to PNG failed.',
 "Aucun dossier d'export disponible": 'No export directory available',
 'Préparation du calcul…': 'Preparing calculation…',
 'Annulation demandée…': 'Cancellation requested…',
 'Le calcul a échoué.': 'Calculation failed.',
 'Calcul annulé — résultats partiels exportés.': 'Calculation cancelled — partial results exported.',
 'Calcul terminé.': 'Calculation complete.',
 'Candidats étendus': 'Extended candidates',
 'Paramètres': 'Parameters',
 'Référence': 'Reference',
 'Recherche': 'Search',
 'Champ': 'Field',
 'Valeur': 'Value',
 'Données': 'Data',
 'Progression MFE': 'MFE progress',
 'Progression ensemble': 'Ensemble progress',
 'Classement (courbe)': 'Ranking (curve)',
 'Courbe': 'Curve',
 'Aucune donnée': 'No data',
 'Exécution': 'Execution',
 'Données longues': 'Long data',
 'Feuille': 'Sheet',
 'Cellule': 'Cell',
 'Partie': 'Part',
 'Écriture des séquences, métriques, structures 2D et courbes…': 'Writing sequences, metrics, 2D structures '
                                                                 'and plots…',
 'Extension annulée — résultats partiels exportés.': 'Extension cancelled — partial results exported.',
 'Extension terminée — résultats exportés.': 'Extension complete — results exported.',
 'Préparation de l’extension NUPACK…': 'Preparing NUPACK extension…',
 'Annulation demandée — fin du calcul NUPACK en cours…': 'Cancellation requested — waiting for the current '
                                                         'NUPACK calculation…',
 'L’extension a échoué.': 'Extension failed.',
 'Extension terminée.': 'Extension complete.',
 'Fichier de run': 'Run file',
 'Fichier de résultats du projet': 'Project results file',
 'Ce dossier ne contient pas de run. Sélectionnez le dossier d’un run ou son fichier de résultats JSON.': 'This '
                                                                                                          'directory '
                                                                                                          'contains '
                                                                                                          'no '
                                                                                                          'runs. '
                                                                                                          'Select '
                                                                                                          'a '
                                                                                                          'run '
                                                                                                          'directory '
                                                                                                          'or '
                                                                                                          'its '
                                                                                                          'JSON '
                                                                                                          'results '
                                                                                                          'file.',
 'Ce dossier contient plusieurs runs. Sélectionnez le fichier de résultats JSON du run à ouvrir.': 'This '
                                                                                                   'directory '
                                                                                                   'contains '
                                                                                                   'several '
                                                                                                   'runs. '
                                                                                                   'Select '
                                                                                                   'the JSON '
                                                                                                   'results '
                                                                                                   'file of '
                                                                                                   'the run '
                                                                                                   'to open.',
 'Le manifeste doit désigner un fichier de résultats dans le même dossier.': 'The manifest must point to a '
                                                                             'results file in the same '
                                                                             'directory.',
 'Le manifeste ne désigne pas un fichier de résultats.': 'The manifest does not point to a results file.',
 'Avertissements': 'Warnings',
 'Avertissement': 'Warning',
 'Paramètres du design': 'Design parameters',
 'Ce run n’est pas un calcul NUPACK. Les anciens runs de démonstration ne peuvent pas être ouverts comme résultats scientifiques.': 'This '
                                                                                                                                    'run '
                                                                                                                                    'is '
                                                                                                                                    'not '
                                                                                                                                    'a '
                                                                                                                                    'NUPACK '
                                                                                                                                    'calculation. '
                                                                                                                                    'Old '
                                                                                                                                    'demonstration '
                                                                                                                                    'runs '
                                                                                                                                    'cannot '
                                                                                                                                    'be '
                                                                                                                                    'opened '
                                                                                                                                    'as '
                                                                                                                                    'scientific '
                                                                                                                                    'results.',
 'Séquences': 'Sequences',
 'Nom du système': 'System name',
 'Début de la région reconnue': 'Start of recognized region',
 'Longueur de la région reconnue': 'Length of recognized region',
 'Linker : séquence invalide.': 'Linker: invalid sequence.',
 'Rapporteur': 'Reporter',
 'Séquence du rapporteur': 'Reporter sequence',
 'Pondérations': 'Weights',
 'Nombre d’essais': 'Number of trials',
 'Candidats du design': 'Design candidates',
 'Candidat du design': 'Design candidate',
 'Numéro d’essai': 'Trial number',
 'Plusieurs candidats portent le même numéro d’essai.': 'Several candidates have the same trial number.',
 'Séquence du switch': 'Switch sequence',
 'Séquence du trigger candidat': 'Candidate trigger sequence',
 'Date de début': 'Start date',
 'Date de fin': 'End date',
 'Statut': 'Status',
 'Entrées de l’extension': 'Extension inputs',
 'Trigger original': 'Original trigger',
 'L’extrémité choisie pour l’extension est inconnue.': 'The selected extension end is unknown.',
 'Position liée au ligand': 'Ligand-binding position',
 'Une base liée au ligand dépasse la longueur de l’aptamère.': 'A ligand-binding base exceeds the aptamer '
                                                               'length.',
 'La structure de référence du run d’extension est absente.': 'The extension run reference structure is '
                                                              'missing.',
 'Candidat étendu': 'Extended candidate',
 'Le trigger étendu ne conserve pas le trigger original entre ses extensions.': 'The extended trigger does '
                                                                                'not preserve the original '
                                                                                'trigger between its '
                                                                                'extensions.',
 'Les longueurs enregistrées ne correspondent pas au trigger étendu.': 'The recorded lengths do not match '
                                                                       'the extended trigger.',
 'Paires de référence': 'Reference base pairs',
 'Paire de référence': 'Reference base pair',
 'Identifiant du run': 'Run ID',
 'L’identifiant du run contient un séparateur de dossier.': 'The run ID contains a directory separator.',
 'Ce JSON n’est pas un résultat de design de switches ou d’extension de triggers. Ouvrez le fichier *_results.json du run.': 'This '
                                                                                                                             'JSON '
                                                                                                                             'is '
                                                                                                                             'not '
                                                                                                                             'a '
                                                                                                                             'switch '
                                                                                                                             'design '
                                                                                                                             'or '
                                                                                                                             'trigger '
                                                                                                                             'extension '
                                                                                                                             'result. '
                                                                                                                             'Open '
                                                                                                                             'the '
                                                                                                                             'run’s '
                                                                                                                             '*_results.json '
                                                                                                                             'file.',
 'Unité de concentration inconnue : {v0}': 'Unknown concentration unit: {v0}',
 'ERREUR · {v0}': 'ERROR · {v0}',
 'Voir Données longues · {v0}!{v1} · {v2} parties': 'See Long data · {v0}!{v1} · {v2} parts',
 'Recherche NUPACK · {v0:,} extensions à examiner.': 'NUPACK search · {v0:,} extensions to examine.',
 '{v0} : un objet JSON est attendu.': '{v0}: a JSON object is expected.',
 '{v0} : une liste est attendue.': '{v0}: a list is expected.',
 '{v0} : texte absent ou invalide.': '{v0}: missing or invalid text.',
 '{v0} : valeur numérique absente ou invalide.': '{v0}: missing or invalid numeric value.',
 '{v0} : entier attendu (minimum {v1}).': '{v0}: an integer is expected (minimum {v1}).',
 '{v0} : séquence de bases invalide.': '{v0}: invalid base sequence.',
 'Impossible de lire le run « {v0} » : {v1}': 'Unable to read run “{v0}”: {v1}',
 'Le fichier « {v0} » n’est pas un JSON valide (ligne {v1}).': 'File “{v0}” is not valid JSON (line {v1}).',
 'Le fichier « {v0} » contient une imbrication JSON invalide.': 'File “{v0}” contains invalid JSON nesting.',
 '{v0} : données manquantes ou invalides.': '{v0}: missing or invalid data.',
 'Pondération · {v0}': 'Weight · {v0}',
 'Candidat · {v0} : booléen attendu.': 'Candidate · {v0}: a boolean is expected.',
 'Candidat · {v0}': 'Candidate · {v0}',
 '{v0} : la structure ne correspond pas aux deux séquences.': '{v0}: the structure does not match the two '
                                                              'sequences.',
 '{v0} : structure non équilibrée.': '{v0}: unbalanced structure.',
 '{v0} · énergie MFE': '{v0} · MFE energy',
 '{v0} · probabilités': '{v0} · probabilities',
 '{v0} : il faut une probabilité par base.': '{v0}: one probability per base is required.',
 '{v0} · probabilité': '{v0} · probability',
 '{v0} : probabilité en dehors de [0, 1].': '{v0}: probability outside [0, 1].',
 'Criblage · {v0}': 'Screening · {v0}',
 'Impossible d’ouvrir ce chemin de run : {v0}': 'Unable to open this run path: {v0}',
 'Moteur : NUPACK — recherche exhaustive d’extension ({v0})': 'Engine: NUPACK — exhaustive extension search '
                                                              '({v0})',
 'Candidats exportés : {v0}': 'Candidates exported: {v0}'}

TRANSLATIONS.update({'Aucun fichier NUPACK compatible avec cette application (Python {python}, {system}, {machine}). Téléchargez le paquet NUPACK pour votre système et déposez son dossier décompressé ici.': 'No '
                                                                                                                                                                                           'NUPACK '
                                                                                                                                                                                           'file '
                                                                                                                                                                                           'is '
                                                                                                                                                                                           'compatible '
                                                                                                                                                                                           'with '
                                                                                                                                                                                           'this '
                                                                                                                                                                                           'app '
                                                                                                                                                                                           '(Python '
                                                                                                                                                                                           '{python}, '
                                                                                                                                                                                           '{system}, '
                                                                                                                                                                                           '{machine}). '
                                                                                                                                                                                           'Download '
                                                                                                                                                                                           'NUPACK '
                                                                                                                                                                                           'for '
                                                                                                                                                                                           'your '
                                                                                                                                                                                           'operating '
                                                                                                                                                                                           'system '
                                                                                                                                                                                           'and '
                                                                                                                                                                                           'place '
                                                                                                                                                                                           'its '
                                                                                                                                                                                           'unzipped '
                                                                                                                                                                                           'folder '
                                                                                                                                                                                           'here.',
 'NUPACK ne prend pas en charge Windows en natif. Les calculs de cette bêta sont disponibles sur macOS et Linux.': 'NUPACK '
                                                                                                                   'does '
                                                                                                                   'not '
                                                                                                                   'support '
                                                                                                                   'native '
                                                                                                                   'Windows. '
                                                                                                                   'Calculations '
                                                                                                                   'in '
                                                                                                                   'this '
                                                                                                                   'beta '
                                                                                                                   'are '
                                                                                                                   'available '
                                                                                                                   'on '
                                                                                                                   'macOS '
                                                                                                                   'and '
                                                                                                                   'Linux.',
 'Le nom du fichier wheel NUPACK est invalide.': 'The NUPACK wheel filename is invalid.',
 'Le fichier wheel contient un chemin non valide.': 'The wheel contains an invalid path.',
 'Ce fichier wheel ne contient pas le paquet NUPACK.': 'This wheel does not contain the NUPACK package.'})
