"""Consented BLASTn screening, taxonomy lookup and portable search reports."""

TRANSLATIONS = {
    "Triticum aestivum (blé), Hordeum vulgare (orge), Fusarium…": "Triticum aestivum (wheat), Hordeum vulgare (barley), Fusarium…",
    "Séquence": "Sequence",
    "Génomes de référence": "Reference genomes",
    "Transcrits ARN de référence": "Reference RNA transcripts",
    "Télécharger le contrôle (JSON)": "Download screening report (JSON)",
    "Télécharger les correspondances (CSV)": "Download matches (CSV)",
    "Contrôle du trigger dans l’échantillon · BLASTn": "Trigger screening in the sample · BLASTn",
    "Avant l’extension ou la conception du toehold switch, recherchez les séquences ressemblantes dans les organismes présents dans l’échantillon.":
        "Before extension or toehold switch design, search for similar sequences in organisms present in the sample.",
    "Ce que ce contrôle peut établir": "What this screening can establish",
    "**Génomes** : recherche dans les génomes de référence, donc l’ADN génomique pour les plantes et champignons. Certains virus possèdent un génome ARN.\n\n**Transcrits ARN** : recherche dans les transcrits de référence, parfois prédits.\n\n**En solution** : BLAST ne mesure ni l’expression, ni la quantité, ni la libération de ce brin. Cela dépend notamment du tissu, de la lyse, de l’extraction et de la dégradation. Un alignement seul ne démontre pas un faux positif.":
        "**Genomes**: search reference genomes, which are genomic DNA for plants and fungi. Some viruses have RNA genomes.\n\n**RNA transcripts**: search reference transcripts, some of which are predicted.\n\n**In solution**: BLAST does not measure expression, quantity or release of this strand. These depend on factors including tissue, lysis, extraction and degradation. An alignment alone does not demonstrate a false positive.",
    "Documentation : [BLAST NCBI](https://blast.ncbi.nlm.nih.gov/doc/blast-help/urlapi.html) · [Références génomiques et transcrits RefSeq](https://www.ncbi.nlm.nih.gov/refseq/about/)":
        "Documentation: [NCBI BLAST](https://blast.ncbi.nlm.nih.gov/doc/blast-help/urlapi.html) · [RefSeq genomes and transcripts](https://www.ncbi.nlm.nih.gov/refseq/about/)",
    "Rechercher un organisme (nom ou TaxID)": "Find an organism (name or TaxID)",
    "Rechercher dans NCBI Taxonomy": "Search NCBI Taxonomy",
    "Saisissez un organisme à rechercher.": "Enter an organism to search for.",
    "Aucun organisme trouvé. Essayez son nom scientifique ou son TaxID NCBI.":
        "No organism found. Try its scientific name or NCBI TaxID.",
    "Organisme à ajouter": "Organism to add",
    "Ajouter à l’échantillon": "Add to sample",
    "Retirer": "Remove",
    "Ajoutez au moins un organisme. Choisissez une espèce précise lorsque vous la connaissez.":
        "Add at least one organism. Choose a specific species when you know it.",
    "Sélectionnez au maximum 20 organismes par contrôle.": "Select up to 20 organisms per screening.",
    "Séquences à rechercher": "Sequences to search",
    "Maximum de séquences retournées par organisme et catégorie": "Maximum sequences returned per organism and category",
    "Trigger à contrôler : {length} nt": "Trigger to screen: {length} nt",
    "Un trigger très court peut correspondre à de nombreuses séquences. L’absence de résultat BLAST ne suffit pas à établir sa spécificité.":
        "A very short trigger may match many sequences. No BLAST results is not enough to establish its specificity.",
    "J’autorise l’envoi de ce trigger et de cette sélection d’organismes à NCBI pour ce contrôle.":
        "I authorize sending this trigger and the selected organisms to NCBI for this screening.",
    "Contrôler le trigger avec BLASTn": "Screen trigger with BLASTn",
    "Confirmez l’envoi de ce trigger à NCBI avant de lancer le contrôle.":
        "Confirm sending this trigger to NCBI before starting the screening.",
    "Attendez la fin du calcul en cours avant de lancer BLASTn.": "Wait for the current calculation to finish before starting BLASTn.",
    "Ouvrez un design avant de lancer BLASTn.": "Open a design before starting BLASTn.",
    "NCBI traite les recherches à distance. Cela peut prendre plusieurs minutes ; l’application respecte les délais de consultation du service.":
        "NCBI processes searches remotely. This may take several minutes; the app follows the service’s polling intervals.",
    "Annuler le contrôle BLASTn": "Cancel BLASTn screening",
    "Dernier contrôle enregistré — une nouvelle recherche est en cours.": "Last saved screening — a new search is running.",
    "Ce contrôle concerne un ancien trigger ou une autre sélection d’organismes. Relancez-le pour les paramètres actuels.":
        "This screening used an earlier trigger or a different organism selection. Run it again for the current settings.",
    "Le contrôle est incomplet. Les recherches en erreur ou annulées ne permettent pas de conclure à une absence de correspondance.":
        "Screening is incomplete. Failed or cancelled searches cannot establish an absence of matches.",
    "Trigger contrôlé : {sequence}": "Screened trigger: {sequence}",
    "Contrôle enregistré : {date}": "Screening saved: {date}",
    "Cette catégorie n’a pas été recherchée.": "This category was not searched.",
    "Base : {database} · Requête NCBI : {rid}": "Database: {database} · NCBI request: {rid}",
    "Recherche non terminée : {error}": "Search not completed: {error}",
    "Aucune correspondance rapportée par cette recherche. Cela ne valide pas l’absence de faux positifs.":
        "No matches reported by this search. This does not establish the absence of false positives.",
    "Accession": "Accession",
    "Organisme": "Organism",
    "Non renseigné par NCBI": "Not provided by NCBI",
    "Identité (%)": "Identity (%)",
    "Trigger couvert (%)": "Trigger coverage (%)",
    "Longueur alignée": "Alignment length",
    "Mésappariements": "Mismatches",
    "Gaps": "Gaps",
    "Brin aligné": "Aligned strand",
    "Relation au trigger": "Relationship to trigger",
    "Même sens": "Same orientation",
    "Complément inversé": "Reverse complement",
    "Non précisée": "Not specified",
    "Description": "Description",
    "Identité et couverture décrivent l’alignement, pas une probabilité d’activation du switch. La liste retournée peut être limitée.":
        "Identity and coverage describe the alignment, not the probability of switch activation. The returned list may be limited.",
    "Voir les alignements et les séquences candidates": "View alignments and candidate sequences",
    "Alignement": "Alignment",
    "Ouvrir la séquence dans NCBI": "Open sequence in NCBI",
    "Le brin complémentaire peut présenter un risque d’hybridation ; il n’est pas automatiquement équivalent à un trigger activateur dans le sens direct.":
        "The complementary strand may pose a hybridization risk; it is not automatically equivalent to an activating trigger in the forward orientation.",
    "Transcrit prédit par annotation : sa présence dans votre échantillon n’est pas démontrée.":
        "Transcript predicted by annotation: its presence in your sample has not been demonstrated.",
    "Référence de transcrit : l’expression dans votre tissu et sa libération en solution restent à vérifier.":
        "Reference transcript: expression in your tissue and release into solution still need to be verified.",
    "Le dossier des recherches BLAST ne peut pas être un lien symbolique.": "The BLAST searches folder cannot be a symbolic link.",
    "Préparation du contrôle BLAST…": "Preparing BLAST screening…",
    "Le service BLAST a renvoyé un état de recherche invalide.": "The BLAST service returned an invalid search status.",
    "Contrôle BLAST terminé.": "BLAST screening completed.",
    "Contrôle BLAST incomplet : certaines recherches ont échoué.": "BLAST screening incomplete: some searches failed.",
    "Le contrôle BLAST a échoué.": "BLAST screening failed.",
    "Contrôle BLAST annulé.": "BLAST screening cancelled.",
    "Le rapport BLAST enregistré est absent ou invalide.": "The saved BLAST report is missing or invalid.",
    "Annulation BLAST demandée ; aucune nouvelle requête ne sera envoyée.": "BLAST cancellation requested; no further requests will be sent.",
    "L’identifiant taxonomique NCBI doit être un entier positif valide.": "The NCBI taxonomy ID must be a valid positive integer.",
    "Sélectionnez un organisme identifié dans NCBI Taxonomy.": "Select an organism identified in NCBI Taxonomy.",
    "La séquence du trigger doit contenir des bases A, C, G, T ou U.": "The trigger sequence must contain A, C, G, T or U bases.",
    "BLASTn nécessite au moins sept bases non ambiguës pour le trigger.": "BLASTn requires at least seven unambiguous trigger bases.",
    "Sélectionnez entre un et vingt organismes pour BLASTn.": "Select between one and twenty organisms for BLASTn.",
    "Sélectionnez les génomes de référence, les transcrits ARN, ou les deux.": "Select reference genomes, RNA transcripts, or both.",
    "Le nombre maximal de séquences BLAST doit être compris entre 1 et 500.": "The maximum number of BLAST sequences must be between 1 and 500.",
    "Indiquez un nom d’organisme ou un identifiant NCBI Taxonomy.": "Enter an organism name or NCBI Taxonomy ID.",
    "Le nombre de résultats taxonomiques doit être compris entre 1 et 50.": "The number of taxonomy results must be between 1 and 50.",
    "NCBI Taxonomy a renvoyé une réponse non valide.": "NCBI Taxonomy returned an invalid response.",
    "NCBI a renvoyé une erreur HTTP {code}.": "NCBI returned HTTP error {code}.",
    "La connexion à NCBI a échoué : {error}": "The connection to NCBI failed: {error}",
    "La réponse NCBI dépasse la taille maximale autorisée.": "The NCBI response exceeds the maximum allowed size.",
    "La recherche BLAST ne correspond pas à la sélection autorisée.": "The BLAST search does not match the authorized selection.",
    "Les délais de connexion et d’attente BLAST doivent être strictement positifs.":
        "The BLAST connection and waiting timeouts must be strictly positive.",
    "Le délai d’attente BLAST est dépassé. La recherche distante peut encore être en cours.":
        "The BLAST waiting period has expired. The remote search may still be running.",
    "Le nombre maximal de séquences retournées est atteint ; d’autres correspondances peuvent ne pas être affichées.":
        "The maximum number of returned sequences has been reached; additional matches may not be shown.",
    "Le rapport BLAST reçu est invalide ou ne correspond pas à la recherche : {error}":
        "The received BLAST report is invalid or does not match the search: {error}",
    "NCBI indique que la recherche BLAST a échoué ou a expiré : {status}.": "NCBI reports that the BLAST search failed or expired: {status}.",
    "NCBI n’a pas confirmé le filtre taxonomique dans le rapport ; la portée des résultats reste à vérifier.":
        "NCBI did not confirm the taxonomy filter in the report; the scope of the results still needs to be checked.",
    "NCBI n’a pas fourni d’identifiant pour la recherche BLAST.": "NCBI did not provide an ID for the BLAST search.",
    "BLASTn · envoi pour {organism} · {scope}": "BLASTn · submitting for {organism} · {scope}",
    "BLASTn · attente NCBI · {organism} · RID {rid}": "BLASTn · waiting for NCBI · {organism} · RID {rid}",
    "BLASTn · recherche terminée pour {organism} · {scope}": "BLASTn · search completed for {organism} · {scope}",
}
