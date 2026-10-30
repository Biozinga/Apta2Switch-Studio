"""Credits and primary references for the bundled scientific tools."""
from __future__ import annotations

import streamlit as st

from aptaswitch_core.localization import tr


# Licences apply to their individual components. A scientific citation is not
# a software licence, and the app's own licence does not replace either one.
CREDITS = (
    ("OSTIR · laboratoire Barrick et contributeurs",
     "Estimation locale de l’initiation de traduction ; comparaison exploratoire entre l’ARN OFF et une projection ARN de l’état ON.",
     "OSTIR 1.1.3 : GPLv3 ou ultérieure. Les sources correspondantes et la licence accompagnent le paquet. Les conditions GPL s’appliquent à la distribution de l’intégration.",
     (("Projet et auteurs", "https://github.com/barricklab/ostir"),
      ("Publication", "https://doi.org/10.21105/joss.03362"),
      ("Licence", "https://github.com/barricklab/ostir/blob/master/LICENSE"))),
    ("ViennaRNA · Université de Vienne",
     "Repliement des brins ARN et paramètres thermodynamiques du moteur expérimental.",
     "Licence ViennaRNA : attribution obligatoire ; redistribution payante soumise à l’accord des auteurs.",
     (("Projet et auteurs", "https://www.tbi.univie.ac.at/RNA/"),
      ("Publication", "https://doi.org/10.1186/1748-7188-6-26"),
      ("Licence", "https://github.com/ViennaRNA/ViennaRNA/blob/v2.7.2/COPYING"))),
    ("oxDNA / oxNA · équipe oxDNA",
     "Simulation moléculaire en gros grains de l’hybridation ADN–ARN et de l’ouverture du switch.",
     "Licence GPLv3. Les notices et les sources correspondantes accompagnent la distribution du moteur.",
     (("Projet et auteurs", "https://github.com/lorenzo-rovigatti/oxDNA"),
      ("Publication", "https://doi.org/10.21105/joss.04693"),
      ("Licence", "https://github.com/lorenzo-rovigatti/oxDNA/blob/master/LICENSE"))),
    ("Sugimoto et collaborateurs · Biopython",
     "Paramètres thermodynamiques d’empilement des duplexes ADN–ARN, publiés en 1995 ; table issue de Biopython.",
     "Table Biopython sous licence BSD à trois clauses. La publication est citée comme source scientifique.",
     (("Publication", "https://doi.org/10.1021/bi00035a029"),
      ("Biopython et licence", "https://github.com/biopython/biopython/blob/master/LICENSE.rst"))),
    ("3Dmol.js · Université de Pittsburgh et contributeurs",
     "Affichage interactif des conformations et des trajectoires moléculaires.",
     "Licence BSD à trois clauses ; notices des composants incluses avec le visualiseur.",
     (("Projet et auteurs", "https://3dmol.org/"),
      ("Licence", "https://github.com/3dmol/3Dmol.js/blob/master/LICENSE"))),
)


def render_licenses() -> None:
    st.subheader(tr("Licences et remerciements"), anchor=False)
    st.markdown(tr("Merci aux équipes qui développent et partagent ces outils et leurs travaux scientifiques."))
    st.caption(tr("Les moteurs et le modèle sont intégrés à l’application. Aucune installation ni validation manuelle n’est nécessaire."))
    for title, purpose, license_text, links in CREDITS:
        with st.container(border=True):
            st.markdown(f"**{tr(title)}**")
            st.write(tr(purpose))
            st.caption(tr(license_text))
            st.markdown(" · ".join(f"[{tr(label)}]({url})" for label, url in links))
    st.caption(tr("Apta2Switch-Studio · équipe iGEM Sorbonne Université. Le code propre de cette livraison utilise la licence Apta2Switch-Studio : redistribution gratuite et attribution, y compris pour les dérivés. Les droits MIT des versions antérieures restent acquis. OSTIR conserve la GPLv3 ou ultérieure, oxDNA la GPLv3 et ViennaRNA sa licence propre."))
