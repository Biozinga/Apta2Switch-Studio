"""Named design library, portable project folders and saved design inputs."""

TRANSLATIONS = {
    "Design actif": "Active design",
    "Mes designs": "My designs",
    "Nouveau design": "New design",
    "Importer un design": "Import a design",
    "Calculs de ce design": "Calculations in this design",
    "Retrouvez un design, créez-en un nouveau ou importez un dossier complet. Séquences, extensions, calculs et structures 3D restent réunis.":
        "Open a design, create a new one or import a complete folder. Sequences, extensions, calculations and 3D structures stay together.",
    "Un calcul est en cours. Vous pourrez changer de design lorsqu’il sera terminé.":
        "A calculation is running. You can switch designs when it finishes.",
    "Revenir au design en cours": "Return to the current design",
    "Designs enregistrés": "Saved designs",
    "Choisir un design": "Choose a design",
    "Modifié le {date} · Aptamère : {aptamer} nt · Trigger : {trigger} nt":
        "Updated {date} · Aptamer: {aptamer} nt · Trigger: {trigger} nt",
    "Ouvrir ce design": "Open this design",
    "Aucun design enregistré pour le moment. Créez votre premier design ci-dessous.":
        "No saved designs yet. Create your first design below.",
    "Nom du design": "Design name",
    "Ex. : Adénosine — essai 1": "E.g.: Adenosine — trial 1",
    "Créer le design": "Create design",
    "Nom après import (facultatif)": "Name after import (optional)",
    "Archive du design (.zip)": "Design archive (.zip)",
    "Importer l’archive": "Import archive",
    "Importer depuis un dossier local": "Import from a local folder",
    "Dossier du design, design.json ou archive ZIP": "Design folder, design.json or ZIP archive",
    "Importer ce design": "Import this design",
    "Attention, ce design existe déjà. Voulez-vous l’écraser ? Une copie de l’ancien dossier sera conservée.":
        "This design already exists. Do you want to overwrite it? A copy of the previous folder will be kept.",
    "Oui, écraser ce design": "Yes, overwrite this design",
    "Annuler": "Cancel",
    "Dossier des designs : {path}": "Designs folder: {path}",
    "Sauvegarde automatique dans le dossier de ce design.": "Automatically saved in this design’s folder.",
    "Dossier et export": "Folder and export",
    "Préparer l’export du design": "Prepare design export",
    "Télécharger le design (.zip)": "Download design (.zip)",
    "Le design a changé. Préparez un nouvel export.": "The design has changed. Prepare a new export.",
    "Attendez la fin du calcul avant de changer de design.": "Wait for the calculation to finish before switching designs.",
    "Attendez la fin du calcul en cours avant de changer de design.": "Wait for the current calculation to finish before switching designs.",
    "Indiquez le dossier ou l’archive du design à importer.": "Enter the folder or archive of the design to import.",
    "Aptamère original": "Original aptamer",
    "Séquences d’origine": "Original sequences",
    "Revenir aux séquences d’origine": "Restore original sequences",
    "Vous utilisez un trigger étendu. Les séquences d’origine restent conservées dans ce design.":
        "You are using an extended trigger. The original sequences remain saved in this design.",
    "Chaque design possède son dossier avec ses séquences, ses calculs et ses structures 3D.":
        "Each design has its own folder containing its sequences, calculations and 3D structures.",
    "Dossier des designs": "Designs folder",
    "Dossier des designs enregistré. Les designs existants restent à leur emplacement actuel.":
        "Designs folder saved. Existing designs remain in their current location.",
    "Le design n’a pas pu être sauvegardé : {error}": "The design could not be saved: {error}",
    "Tous les calculs sont enregistrés dans le dossier du design actif.":
        "All calculations are saved in the active design’s folder.",
    "Choisissez un calcul de ce design. Les anciens calculs externes s’importent depuis Mes designs.":
        "Choose a calculation from this design. Earlier calculations from other folders can be imported from My designs.",
    "Importez ce calcul depuis Mes designs pour le rattacher à un design.":
        "Import this calculation from My designs to add it to a design.",
    "Retrouvez les calculs et extensions de ce design.": "Find this design’s calculations and extensions.",
    "Les paramètres enregistrés du design sont invalides : {field}.":
        "The saved design settings are invalid: {field}.",
    "Le dossier de séquences ne peut pas être un lien symbolique.": "The sequences folder cannot be a symbolic link.",
    "Le design « {name} » existe déjà.": "The design “{name}” already exists.",
    "Le nom du design doit être un texte.": "The design name must be text.",
    "Choisissez un nom de design sans séparateur de dossier ni caractère spécial.":
        "Choose a design name without folder separators or special characters.",
    "Ce nom de design est réservé au système.": "This design name is reserved by the system.",
    "Les données du design doivent être un objet JSON.": "The design data must be a JSON object.",
    "Les données du design sont trop volumineuses.": "The design data is too large.",
    "Les données du design ne peuvent pas être enregistrées en JSON.": "The design data cannot be saved as JSON.",
    "Le dossier ne contient pas de manifeste design.json valide.": "The folder does not contain a valid design.json manifest.",
    "Le manifeste du design est trop volumineux.": "The design manifest is too large.",
    "Le format ou la version de ce design est inconnu.": "This design’s format or version is unknown.",
    "Les dates du design sont invalides.": "The design dates are invalid.",
    "Impossible de lire le manifeste du design : {error}": "Could not read the design manifest: {error}",
    "Le manifeste du design ne peut pas être un lien symbolique.": "The design manifest cannot be a symbolic link.",
    "Sélectionnez un dossier de design valide.": "Select a valid design folder.",
    "Plusieurs dossiers portent ce nom de design. Renommez-les avant de continuer.":
        "Several folders have this design name. Rename them before continuing.",
    "Un fichier ou un lien porte déjà ce nom. Choisissez un autre nom de design.":
        "A file or link already has this name. Choose another design name.",
    "Le dossier de sauvegarde ne peut pas être un lien symbolique.": "The backup folder cannot be a symbolic link.",
    "Le design dépasse la taille maximale autorisée pour un import ou un export.":
        "The design exceeds the maximum size allowed for import or export.",
    "L’archive contient trop de fichiers.": "The archive contains too many files.",
    "L’archive contient un chemin de fichier invalide.": "The archive contains an invalid file path.",
    "Les liens symboliques et fichiers spéciaux ne sont pas acceptés dans un design.":
        "Symbolic links and special files are not allowed in a design.",
    "Les archives chiffrées ne sont pas acceptées.": "Encrypted archives are not supported.",
    "L’archive contient plusieurs fichiers portant le même nom.": "The archive contains several files with the same name.",
    "Un fichier de l’archive est trop volumineux.": "A file in the archive is too large.",
    "L’archive décompressée est trop volumineuse.": "The uncompressed archive is too large.",
    "L’archive contient un conflit entre un fichier et un dossier.": "The archive contains a conflict between a file and a folder.",
    "L’archive ne contient pas un seul dossier de design valide.": "The archive must contain exactly one valid design folder.",
    "L’archive doit contenir uniquement le dossier du design.": "The archive must contain only the design folder.",
    "L’archive du design est trop volumineuse.": "The design archive is too large.",
    "La taille d’un fichier de l’archive est invalide.": "A file in the archive has an invalid size.",
    "Impossible d’ouvrir cette archive de design : {error}": "Could not open this design archive: {error}",
    "Importez ce design dans une bibliothèque située hors de son dossier.":
        "Import this design into a library outside its own folder.",
}
