# Apta2Switch-Studio — Beta 2

Local application for designing and analysing **DNA aptamer–trigger / RNA toehold-switch systems**, developed by the iGEM Sorbonne University team. The interface is available in English and French. Version: **0.2.1b2**.

[Download the sources](https://github.com/Biozinga/Apta2Switch-Studio/archive/refs/heads/main.zip) · [Scientific methods](https://2026.igem.wiki/sorbonne-university/model/software)

## Run from source

Use **CPython 3.11** on macOS or Linux, with a C/C++ compiler, Make, CMake and pkg-config. On macOS, install Xcode Command Line Tools first. Build scripts download pinned scientific sources and verify their SHA-256 values; internet access is required for missing dependencies.

```sh
git clone https://github.com/Biozinga/Apta2Switch-Studio.git
cd Apta2Switch-Studio
python3.11 -m venv .venv
.venv/bin/python -m pip install -e ".[build]"
.venv/bin/python scripts/prepare_thermodynamics_runtime.py
.venv/bin/aptaswitch-studio
```

For coarse-grained simulations, prepare the additional worker with `.venv/bin/python scripts/prepare_bundled_science.py`. Scientific runtimes are discovered under the ignored `.local_deps/` directory. Designs, calculations, settings and credentials remain local.

Desktop builds use the same sources:

```sh
# macOS: local, unsigned application
APTASWITCH_PYTHON=python3.11 APTASWITCH_MACOS_SIGNING=local bash scripts/build_macos.sh

# Linux x86_64
APTASWITCH_PYTHON=python3.11 bash scripts/build_linux.sh
```

Add `--check-env` to check prerequisites without building. Outputs are written to `dist/`; keep the complete extracted application folder together. A local macOS build is not notarised. Windows launcher/WSL2 infrastructure is supplied, but this delivery has no newly validated Windows release. [Dependency constraints](requirements-constraints.txt) record one development environment, rather than a universal platform lockfile.

## Use the application

Create or reopen a named design, supply the aptamer, trigger and switch, then generate candidates or use **Reading** mode to analyse supplied sequences. The aptamer and trigger remain DNA; the switch and reporter transcript remain RNA. In templates, **X** marks a nucleotide to optimise.

Compare recorded metrics and structures, save selected candidates and export their sequences or results as JSON, CSV, XLSX or portable design archives. Saved results retain their recorded sequences, conditions and model identity. A previously stored calculation is not a new measurement.

## Model and limits

The **experimental native hybrid engine** combines ViennaRNA homogeneous motifs, Sugimoto RNA/DNA stacks and an explicit hybrid-loop approximation. Each strand retains its actual chemistry. Supported new calculations contain **up to three strands**, at **25–55 °C**, **Na⁺ 0.12–1 M** and **Mg²⁺ 0**. The hybrid-loop approximation is an Apta2Switch implementation, not an official ViennaRNA mixed-chemistry model.

- **oxDNA FFS remains unfinished:** no complete end-to-end run was completed for this delivery.
- **Full-reporter context analysis is non-functional and untested.** Storing complete reporter sequences does not establish full-transcript analysis.
- Implementation tests, docking, short simulations and model scores do not establish experimental performance.

## Licences and sources

Original code uses the [Apta2Switch-Studio Source License](LICENSE), `LicenseRef-Apta2Switch-Source`. It permits research, educational and commercial use and modification; redistribution must be free apart from actual physical-media costs. Preserve attribution, notices and identification of modifications. This custom licence is not OSI-approved; earlier [MIT notices](LICENSES/MIT-legacy.txt) remain valid.

Dependencies retain their own terms. [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md) records attribution, reporter provenance and distribution conditions. [UPSTREAM-SOURCES.json](UPSTREAM-SOURCES.json) records exact scientific-source versions, URLs and hashes. No combined binary release is announced; compatibility of the ViennaRNA/OSTIR combination remains to be clarified before binary redistribution.
