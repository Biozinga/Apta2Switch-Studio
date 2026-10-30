"""Build an isolated experimental hybrid library against audited Vienna sources.

Nothing is downloaded or installed. The supplied source tree must contain the
already built static library. Build provenance records all actual source and
library hashes; it does not establish scientific validation by itself.
"""
import argparse
import ast
import hashlib
import json
from pathlib import Path
import platform
import re
import subprocess
import tempfile


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def parameter_header(path):
    tree = ast.parse(path.read_text())
    table = next(ast.literal_eval(node.value) for node in tree.body
                 if isinstance(node, ast.Assign)
                 and any(isinstance(target, ast.Name) and target.id == "_SUGIMOTO_HS" for target in node.targets))
    # Keep the source attribution/license in the generated file and avoid a
    # second independently maintained scientific table.
    comment = ast.get_docstring(tree).replace("*/", "* /")
    complement = str.maketrans("ACGT", "TGCA")
    values = [table[a+b+"/"+(a+b).translate(complement)] for a in "ACGT" for b in "ACGT"]
    rows = ",\n".join("{%r, %r}" % tuple(pair) for pair in values)
    return ("/*\n" + comment + "\n*/\n"
            + "static constexpr double hybrid_hs_init[2]={%r,%r};\n" % tuple(table["init"])
            + "static constexpr double hybrid_hs_stacks[16][2]={\n" + rows + "\n};\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--vienna-source", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--compiler", default="clang++")
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    native = root / "src/aptaswitch_core/thermodynamics/native"
    oracle = native.parent / "hybrid_motifs.py"
    vienna = args.vienna_source.resolve()
    static = vienna / "src/ViennaRNA/libRNA.a"
    params = vienna / "src/ViennaRNA/params/params.c"
    config = vienna / "src/ViennaRNA/vrna_config.h"
    versions = re.findall(r'^\s*#\s*define\s+VRNA_VERSION\s+"([^"]+)"',
                          config.read_text(), flags=re.MULTILINE)
    if versions != ["2.7.2"]:
        raise ValueError("Vienna headers must declare the audited RNAlib version 2.7.2")
    known = {"a6da5b5038e07c6b4b0b843429115db84b6bcb8b48e51074f3180d2842358e7a",
             "0163c59586c5389acc0b1284a97fda7fd0652716b1d26548466e5df735956cf9"}
    params_hash, config_hash = digest(params), digest(config)
    if params_hash not in known:
        raise ValueError("Vienna parameter source is not an audited 2.7.2 version")
    sources = [native / name for name in ("connected_dp.cpp", "hybrid_provider.cpp", "hybrid_api.cpp")]
    headers = [native / name for name in ("connected_dp.hpp", "hybrid_provider.hpp", "hybrid_api.hpp")]
    source_hashes = {str(p.relative_to(root)):digest(p) for p in sources+headers+[oracle]}
    static_hash = digest(static)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    generated = parameter_header(oracle)
    with tempfile.TemporaryDirectory(prefix="aptaswitch-hybrid-build-") as temp:
        generated_path = Path(temp) / "hybrid_parameters.hpp"
        generated_path.write_text(generated)
        binary = Path(temp) / args.output.name
        exports = Path(temp) / "exports.txt"
        if platform.system() == "Darwin":
            exports.write_text("\n".join("_"+name for name in (
                "aptadp_abi_version", "aptadp_create_v1", "aptadp_create_hybrid_v1", "aptadp_create_rna_v1", "aptadp_destroy",
                "aptadp_log_partition", "aptadp_pair_probabilities", "aptadp_minimum", "aptadp_minimum_pair_energies",
                "aptadp_minimum_structures", "aptadp_sample", "aptadp_memory_estimate"))+"\n")
            link_options = ["-Wl,-exported_symbols_list,"+str(exports)]
        elif platform.system() == "Linux":
            exports.write_text("{ global: aptadp_*; local: *; };\n")
            link_options = ["-Wl,--version-script="+str(exports), "-Wl,--exclude-libs,ALL"]
        else:
            raise ValueError("This isolated build script currently supports macOS and Linux toolchains only")
        command = [args.compiler,"-std=c++17","-O2","-Wall","-Wextra","-Werror","-fPIC",
                   "-dynamiclib" if platform.system() == "Darwin" else "-shared",
                   "-I"+str(vienna/"src"),"-I"+str(vienna),"-I"+temp,
                   *(str(source) for source in sources), str(static), "-pthread", "-lm",
                   *link_options, "-o", str(binary)]
        subprocess.run(command,check=True)
        if (source_hashes != {str(p.relative_to(root)):digest(p) for p in sources+headers+[oracle]}
                or static_hash != digest(static) or params_hash != digest(params)
                or config_hash != digest(config)):
            raise RuntimeError("Build input changed during compilation; rerun after sources are stable")
        # Replace a complete file, preserving any old inode still used by a
        # validation process. Never overwrite a loaded dynamic library in place.
        import shutil
        with tempfile.NamedTemporaryFile(dir=args.output.parent, delete=False) as out:
            staged = Path(out.name)
            with binary.open("rb") as source:
                shutil.copyfileobj(source, out)
        staged.replace(args.output)
    manifest = dict(
        abi=1, status="built", model="experimental-Sugimoto95-Vienna04-loopmean",
        platform=platform.platform(), machine=platform.machine(),
        compiler=subprocess.check_output([args.compiler,"--version"],text=True).splitlines()[0],
        sources=source_hashes,
        vienna_version=versions[0],
        single_rna_solver="rnalib2-pf-bpp-sample-mfe-v1",
        single_rna_minimum_convention="turner2004-centikcal-salt-rounded-v1",
        connected_dp_storage="lazy-inside-v1",
        generated_parameters_sha256=hashlib.sha256(generated.encode()).hexdigest(),
        vienna_params_c_sha256=params_hash,vienna_config_h_sha256=config_hash,
        vienna_static_library_sha256=static_hash,
        native_library_sha256=digest(args.output),
    )
    args.output.with_suffix(args.output.suffix+".json").write_text(json.dumps(manifest,indent=2)+"\n")
    print(args.output.resolve())


if __name__ == "__main__":
    main()
