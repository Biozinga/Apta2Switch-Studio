# Third-party notices

Apta2Switch-Studio's original code in this source release is licensed under
**Apta2Switch-Studio Source License** (`LicenseRef-Apta2Switch-Source`); see
[LICENSE](LICENSE). This custom license permits research, educational and
commercial use and modification; redistribution of the covered material and
derived versions must be free, apart from actual physical-media costs. It
requires retained attribution and notices. It is not an OSI-approved license.
The MIT notice from earlier releases remains in [LICENSES/MIT-legacy.txt](LICENSES/MIT-legacy.txt). Existing
third-party copyright notices and licenses are retained and are not replaced
by the application's custom license.

## Material included in the application source

| Material | License / notice | Location |
| --- | --- | --- |
| 3Dmol.js and incorporated visualization components | BSD-3-Clause and the component notices retained by upstream | `src/aptaswitch_studio/assets/3dmol/LICENSE.txt`, `3Dmol-min.js.LICENSE.txt` and `manifest.json` |
| Small Sugimoto RNA/DNA thermodynamic table reproduced from Biopython | BSD-3-Clause; Sebastian Bassi, Markus Piotrowski and the Biopython Contributors | `src/aptaswitch_core/thermodynamics/hybrid_motifs.py` and `LICENSES/Biopython-table-BSD-3-Clause.txt` |
| Apta2Switch-Studio code received under the earlier MIT release | Original MIT copyright and permission notice retained | `LICENSES/MIT-legacy.txt` |

## Scientific sources used by local builds and source archives

The GitHub repository records the following upstream sources in
[UPSTREAM-SOURCES.json](UPSTREAM-SOURCES.json), including exact versions,
download URLs, SHA-256 values and licence paths. Local builders download and
verify these sources. The optional generated source aggregate keeps their
original archives separate and under their own licences, with a manifest for
that generated output. These works are not relicensed under the application's
licence.

| Component | Version / revision | License and upstream source |
| --- | --- | --- |
| ViennaRNA / RNAlib | 2.7.2 | ViennaRNA's own license, **not GPL**: [`COPYING`](https://github.com/ViennaRNA/ViennaRNA/blob/v2.7.2/COPYING), also reproduced verbatim in `LICENSES/ViennaRNA-2.7.2.txt`. Its `src/ViennaRNA/plotting/naview.c` has an additional copyright notice retained inside the upstream archive. |
| OSTIR | 1.1.3 | GPL version 3 or later, as identified in the exact release's Python package metadata; full `LICENSE` inside the upstream archive. [Upstream project](https://github.com/barricklab/ostir). |
| oxDNA / oxpy | `75ae96f063954b2f29524bb897c274aee7ba4466` | GNU GPL version 3; [`LICENSE` at the pinned revision](https://github.com/lorenzo-rovigatti/oxDNA/blob/75ae96f063954b2f29524bb897c274aee7ba4466/LICENSE). Other copyright and component notices remain in that archive. |

Credit for ViennaRNA belongs to its authors and the Institute for Theoretical
Chemistry of the University of Vienna. The ViennaRNA license permits research,
education and commercial use and modification, subject to its own conditions;
it restricts redistribution for a fee beyond media costs and asks that its
authors be contacted about inclusion in a commercial product. The exact text
is authoritative. Apta2Switch-Studio uses its own custom license with a
similar fee restriction; it does not apply the ViennaRNA license to our own
code or grant additional rights to ViennaRNA.

## Source release and locally built applications

The GitHub source delivery contains application sources and upstream-source
references. The optional aggregate additionally carries separate upstream
archives; neither contains a frozen or prelinked combined application.
Build commands compile an application locally on the recipient's own machine.
Python, the compiler, operating-system libraries and the remaining Python
dependencies are obtained separately by the build process and retain their
respective licenses.

OSTIR uses ViennaRNA. ViennaRNA's redistribution restriction and GPL
requirements are not made compatible merely by adding source files or by
selecting a custom ViennaRNA-style license for Apta2Switch-Studio. Our own
restriction on paid redistribution also cannot override OSTIR or oxDNA GPL
requirements. This release does not assert permission to redistribute the resulting combined binary. Such distribution requires
resolving the applicable licenses with the relevant rightsholders. The same
review applies to any other combination of our code with third-party
components. No exception for upstream code is invented here.

The local freezing scripts retain scientific sources, manifests and component
licenses beside the generated runtimes. Those collected files must accompany
any subsequently authorized binary release. A list of names alone is not a
replacement for complete notices or corresponding source.

## Other dependencies

Direct dependencies and build dependencies are declared in `pyproject.toml`.
They include Streamlit, openpyxl, resvg_py, NumPy, SciPy, pandas, PyYAML, Jinja2,
packaging, PyInstaller and Pillow, plus their transitive dependencies. They
are not relicensed by this repository. Their exact installed versions,
copyright notices and license files must be preserved when distributing a
runtime; PyInstaller's own license includes a special bundling exception and
does not replace the licenses of bundled programs.


## Runtime source retention

Local builders retain the exact upstream archives and complete notices alongside generated scientific runtimes. Thermodynamics outputs include the library manifest, ViennaRNA copying terms and the incorporated Biopython attribution. The oxDNA runtime retains its pinned source, worker/helper sources, build instructions and dependency notices; its runtime manifest records payload hashes and parameters. OSTIR source delivery includes the exact OSTIR/ViennaRNA archives, full licences, application-source snapshot and hash manifest. These corresponding sources and complete notices must accompany any subsequently authorised binary redistribution. Python packages and platform libraries retain their own licences.

## Reporter sequences

Built-in coding sequences in [reporter_sequences.py](src/aptaswitch_core/reporter_sequences.py) retain versioned NCBI accessions, CDS coordinates and SHA-256 values. sfGFP and mNeonGreen were retrieved on 27 September 2026; the other entries on 8 October 2026. No reverse translation or codon optimisation was applied.

| Reporter | NCBI record | CDS coordinates, inclusive | Complete CDS | After initial ATG |
| --- | --- | --- | ---: | ---: |
| sfGFP | [HQ873313.1](https://www.ncbi.nlm.nih.gov/nuccore/HQ873313.1) | 1–717, forward | 717 nt | 714 nt |
| mNeonGreen | [KC295282.1](https://www.ncbi.nlm.nih.gov/nuccore/KC295282.1) | 1–711, forward | 711 nt | 708 nt |
| mCherry | [AY678264.1](https://www.ncbi.nlm.nih.gov/nuccore/AY678264.1) | 1–711, forward | 711 nt | 708 nt |
| mRFP1 | [AF506027.1](https://www.ncbi.nlm.nih.gov/nuccore/AF506027.1) | 1–678, forward | 678 nt | 675 nt |
| E. coli K-12 MG1655 LacZ | [NC_000913.3](https://www.ncbi.nlm.nih.gov/nuccore/NC_000913.3?from=363231&to=366305&strand=2) | complement of 363231–366305 | 3075 nt | 3072 nt |

Every retained CDS includes its terminal stop codon. The flanking TGA after sfGFP's annotated CDS is excluded. LacZ is the complete enzyme CDS, not an alpha-complementation fragment. Assembly removes the initial ATG because the switch supplies its own AUG, then converts the reporter transcript to RNA. The usual design fragment is the first 30 bases after the start; complete sequence and provenance remain separate. Custom codon variants must be recorded as distinct custom reporters.

[Reporter integrity tests](tests/test_reporter_full_sequences.py) check hashes, start/stop codons, internal stops and translations against stored records. Sequence storage and these checks do not validate full-transcript analysis or biological switch performance; see [model limits](README.md#model-and-limits). The [NCBI Molecular Data Usage policy](https://www.ncbi.nlm.nih.gov/home/about/policies/) does not transfer submitters' intellectual-property rights. Credit the accessions and associated authors/publications; the application's code licence does not relicense biological sequences or reagents.
