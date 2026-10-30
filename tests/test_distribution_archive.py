from __future__ import annotations

import importlib.util
import ast
import io
import json
import hashlib
from pathlib import Path
import shutil
import stat
import subprocess
import sys
import tarfile
import tempfile
import unittest
from zipfile import ZipFile, ZipInfo


ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("bundle_checker", ROOT / "scripts/check_bundled_app.py")
checker = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(checker)
WORKER_FILES = (
    "oxdna_ffs_worker.py", "oxdna_kinetics_worker.py", "oxdna_native_ffs.py",
    "oxdna_kinetics.py", "oxdna_initialization.py", "oxdna_guided_worker.py",
    "oxdna_geometry.py", "oxdna_loop_geometry.py", "oxdna_box.py",
    "oxdna_free_md_worker.py", "oxdna_free_md_metrics.py",
)
SCIENCE_FILES = {
    "oxdna_runtime": (
        "worker/aptaswitch-oxdna-worker", "worker/_internal/bin/oxDNA", "LICENSE.txt",
        "parameters/oxDNA2_sequence_dependent_parameters.txt",
        "parameters/rna_sequence_dependent_parameters.txt",
        "parameters/NA_sequence_dependent_parameters.txt", "sources/oxDNA-source.tar.gz",
        "sources/BUILDING.txt", "sources/prepare_bundled_science.py",
        "sources/prepare_oxdna_ffs_runtime.py", "sources/run_oxdna_worker.py",
        "sources/rebuild_oxdna_worker.py",
        *("sources/" + name for name in WORKER_FILES),
    ),
}


def science_manifest(directory):
    manifest = {"schema_version": 1, "kind": "oxdna",
                "files": {name: hashlib.sha256(b"fixture").hexdigest() for name in SCIENCE_FILES[directory]},
                "source_revision": "fixture"}
    if directory == "oxdna_runtime":
        manifest.update(worker="worker/aptaswitch-oxdna-worker", binary="worker/_internal/bin/oxDNA",
                        parameters=dict(zip(("DNA", "RNA", "NA"), SCIENCE_FILES[directory][3:6])),
                        license="GPL-3.0")
    return manifest


class DistributionArchiveTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.work = Path(temporary.name)

    def entries(self, kind):
        root = "Apta2Switch-Studio"
        if kind == "zip":
            application = root + "/Apta2Switch-Studio.app/Contents"
            executable = application + "/MacOS/Apta2Switch-Studio"
            base = application + "/Resources/"
        else:
            executable = root + "/Apta2Switch-Studio"
            base = root + "/_internal/"
        result = {root + "/" + name: True for name in ("runs/designs", "runs/extensions")}
        result.update({name: False for name in (executable, root + "/START-HERE.txt", root + "/LICENSE")})
        viewer = base + "aptaswitch_studio/assets/3dmol/"
        result.update({viewer + name: False for name in (
            "3Dmol-min.js", "LICENSE.txt", "3Dmol-min.js.LICENSE.txt", "manifest.json",
        )})
        library = "libaptaswitch_hybrid.dylib" if kind == "zip" else "libaptaswitch_hybrid.so"
        result.update({base + "thermodynamics/" + name: False for name in (
            library, library + ".json", "ViennaRNA-COPYING.txt", "THIRD-PARTY-NOTICES.txt")})
        result.update({base + "aptaswitch_core/" + name: False for name in WORKER_FILES})
        result.update({base + name: False for name in (
            "ostir/rna_turner1999.par", "ostir/rna_turner2004.par",
            "OSTIR-1.1.3.dist-info/METADATA", "OSTIR-1.1.3.dist-info/LICENSE",
            "viennarna-2.7.2.dist-info/METADATA", "viennarna-2.7.2.dist-info/licenses/COPYING",
            "RNA/_RNA.cpython-311.so", "ostir_sources/manifest.json",
            *("ostir_sources/" + name for name in checker.OSTIR_SOURCE_FILES),
        )})
        for directory, files in SCIENCE_FILES.items():
            result.update({base + directory + "/" + name: False for name in ("runtime.json", *files)})
        return result, base

    def write(self, kind, entries, contents=None):
        file = self.work / ("app.zip" if kind == "zip" else "app.tar.gz")
        def data(name):
            if contents and name in contents:
                return contents[name]
            if "/thermodynamics/libaptaswitch_hybrid." in name and name.endswith(".json"):
                return json.dumps({"status": "built", "abi": 1,
                    "model": "experimental-Sugimoto95-Vienna04-loopmean",
                    "native_library_sha256": hashlib.sha256(b"fixture").hexdigest()}).encode()
            if name.endswith("/ostir_sources/manifest.json"):
                return json.dumps({"ostir_version": "1.1.3", "viennarna_version": "2.7.2",
                    "license": "GPL-3.0-or-later", "files": {
                        filename: hashlib.sha256(b"fixture").hexdigest()
                        for filename in checker.OSTIR_SOURCE_FILES}}).encode()
            for directory in SCIENCE_FILES:
                if name.endswith("/" + directory + "/runtime.json"):
                    return json.dumps(science_manifest(directory)).encode()
            return b"fixture"
        if kind == "zip":
            with ZipFile(file, "w") as archive:
                for name, directory in entries.items():
                    if isinstance(directory, str):
                        member = ZipInfo(name)
                        member.create_system = 3
                        member.external_attr = (stat.S_IFLNK | 0o777) << 16
                        archive.writestr(member, directory)
                    else:
                        archive.writestr(name + "/" if directory else name, b"" if directory else data(name))
        else:
            with tarfile.open(file, "w:gz") as archive:
                for name, directory in entries.items():
                    member = tarfile.TarInfo(name)
                    if isinstance(directory, str):
                        member.type = tarfile.SYMTYPE
                        member.linkname = directory
                        archive.addfile(member)
                    elif directory:
                        member.type = tarfile.DIRTYPE
                        archive.addfile(member)
                    else:
                        member.size = len(data(name))
                        archive.addfile(member, io.BytesIO(data(name)))
        return file

    def test_mac_zip_and_linux_tar_include_runtime_and_empty_user_folders(self):
        for kind in ("zip", "tar"):
            entries, _ = self.entries(kind)
            with self.subTest(kind=kind):
                checker.check_distribution_archive(self.write(kind, entries))

    def test_bundled_science_payload_files_and_notices_cannot_be_dropped(self):
        for kind in ("zip", "tar"):
            for directory, files in SCIENCE_FILES.items():
                for missing in ("runtime.json", *files):
                    entries, base = self.entries(kind)
                    entries.pop(base + directory + "/" + missing)
                    with self.subTest(kind=kind, missing=missing), self.assertRaises(RuntimeError):
                        checker.check_distribution_archive(self.write(kind, entries))

    def test_ostir_tables_native_module_and_corresponding_source_are_required(self):
        for kind in ("zip", "tar"):
            for relative in ("ostir/rna_turner2004.par", "RNA/_RNA.cpython-311.so",
                             "OSTIR-1.1.3.dist-info/LICENSE", "ostir_sources/Apta2Switch-source.tar.gz",
                             "ostir_sources/ostir-1.1.3.tar.gz", "ostir_sources/manifest.json"):
                entries, base = self.entries(kind)
                entries.pop(base + relative)
                with self.subTest(kind=kind, missing=relative), self.assertRaises(RuntimeError):
                    checker.check_distribution_archive(self.write(kind, entries))

    def test_ostir_source_hashes_are_checked(self):
        for kind in ("zip", "tar"):
            entries, base = self.entries(kind)
            with self.subTest(kind=kind), self.assertRaisesRegex(RuntimeError, "OSTIR corresponding source"):
                checker.check_distribution_archive(self.write(kind, entries, {
                    base + "ostir_sources/Apta2Switch-source.tar.gz": b"changed source"}))

    def test_bundled_science_payload_hashes_are_checked(self):
        for kind in ("zip", "tar"):
            entries, base = self.entries(kind)
            for directory in SCIENCE_FILES:
                changed = base + directory + "/" + SCIENCE_FILES[directory][0]
                with self.subTest(kind=kind, directory=directory), self.assertRaises(RuntimeError):
                    checker.check_distribution_archive(self.write(kind, entries, {changed: b"changed"}))

    def test_framework_file_and_directory_links_resolve_before_hashing(self):
        for kind in ("zip", "tar"):
            entries, base = self.entries(kind)
            runtime = base + "oxdna_runtime/"
            framework = "worker/_internal/Python.framework/"
            concrete = framework + "Versions/3.11/Python"
            entries[runtime + concrete] = False
            entries[runtime + framework + "Versions/Current"] = "3.11"
            entries[runtime + framework + "Python"] = "Versions/Current/Python"
            manifest = science_manifest("oxdna_runtime")
            for relative in (concrete, framework + "Python", framework + "Versions/Current/Python"):
                manifest["files"][relative] = hashlib.sha256(b"fixture").hexdigest()
            contents = {runtime + "runtime.json": json.dumps(manifest).encode()}
            with self.subTest(kind=kind):
                checker.check_distribution_archive(self.write(kind, entries, contents))
                with self.assertRaisesRegex(RuntimeError, "differs from its manifest"):
                    checker.check_distribution_archive(self.write(kind, entries, {**contents, runtime + concrete: b"tampered"}))

    def test_framework_directory_links_cannot_cycle_or_escape(self):
        for kind in ("zip", "tar"):
            for target, error in (("Current", "Cyclic"), ("../../../../../../../../../../outside", "escapes"),
                                  ("/Apta2Switch-Studio/elsewhere", "escapes")):
                entries, base = self.entries(kind)
                runtime = base + "oxdna_runtime/"
                framework = "worker/_internal/Python.framework/"
                entries[runtime + framework + "Versions/Current"] = target
                entries[runtime + framework + "Python"] = "Versions/Current/Python"
                manifest = science_manifest("oxdna_runtime")
                manifest["files"][framework + "Python"] = hashlib.sha256(b"fixture").hexdigest()
                contents = {runtime + "runtime.json": json.dumps(manifest).encode()}
                with self.subTest(kind=kind, target=target), self.assertRaisesRegex(RuntimeError, error):
                    checker.check_distribution_archive(self.write(kind, entries, contents))

    def test_only_the_declared_model_can_be_bundled(self):
        for kind in ("zip", "tar"):
            entries, base = self.entries(kind)
            entries[base + "rna_score_runtime/another-model.h5"] = False
            with self.subTest(kind=kind), self.assertRaises(RuntimeError):
                checker.check_distribution_archive(self.write(kind, entries))

    def test_science_manifests_cannot_authorize_paths_outside_the_runtime(self):
        for kind in ("zip", "tar"):
            entries, base = self.entries(kind)
            for bad in ("../private.h5", "/private.h5", "dir\\private.h5"):
                manifest = science_manifest("oxdna_runtime")
                manifest["files"][bad] = hashlib.sha256(b"fixture").hexdigest()
                contents = {base + "oxdna_runtime/runtime.json": json.dumps(manifest).encode()}
                with self.subTest(kind=kind, path=bad), self.assertRaises(RuntimeError):
                    checker.check_distribution_archive(self.write(kind, entries, contents))

    def test_native_hybrid_binary_manifest_and_notices_are_required(self):
        for kind in ("zip", "tar"):
            entries, runtime = self.entries(kind)
            prefix = runtime + "thermodynamics/"
            for missing in [name for name in entries if name.startswith(prefix)]:
                with self.subTest(kind=kind, missing=missing), self.assertRaisesRegex(RuntimeError, "thermodynamics"):
                    checker.check_distribution_archive(self.write(kind, {name: value for name, value in entries.items() if name != missing}))

    def test_native_hybrid_hash_is_verified_after_distribution_signing(self):
        for kind in ("zip", "tar"):
            entries, runtime = self.entries(kind)
            prefix = runtime + "thermodynamics/"
            library = prefix + ("libaptaswitch_hybrid.dylib" if kind == "zip" else "libaptaswitch_hybrid.so")
            with self.subTest(kind=kind), self.assertRaisesRegex(RuntimeError, "differs from its manifest"):
                checker.check_distribution_archive(self.write(kind, entries, {library: b"modified after signing"}))

    def test_native_hybrid_model_and_nonempty_notices_are_checked(self):
        for kind in ("zip", "tar"):
            entries, runtime = self.entries(kind)
            prefix = runtime + "thermodynamics/"
            manifest = next(name for name in entries if name.startswith(prefix) and name.endswith(".json"))
            for contents in ({manifest: b"{}"}, {prefix + "ViennaRNA-COPYING.txt": b""}):
                with self.subTest(kind=kind, contents=contents), self.assertRaises(RuntimeError):
                    checker.check_distribution_archive(self.write(kind, entries, contents))

    def test_archive_without_runtime_is_rejected(self):
        for kind in ("zip", "tar"):
            entries, runtime = self.entries(kind)
            entries = {name: directory for name, directory in entries.items() if not name.startswith(runtime)}
            with self.subTest(kind=kind), self.assertRaises(RuntimeError):
                checker.check_distribution_archive(self.write(kind, entries))

    def test_molecular_viewer_and_its_license_notices_must_be_embedded(self):
        for kind in ("zip", "tar"):
            for missing in ("3Dmol-min.js", "LICENSE.txt", "3Dmol-min.js.LICENSE.txt", "manifest.json"):
                entries, runtime = self.entries(kind)
                viewer = runtime + "aptaswitch_studio/assets/3dmol/"
                entries.pop(viewer + missing)
                with self.subTest(kind=kind, missing=missing), self.assertRaisesRegex(RuntimeError, "molecular viewer"):
                    checker.check_distribution_archive(self.write(kind, entries))

    def test_unmanifested_external_engines_and_models_must_not_ship(self):
        for kind in ("zip", "tar"):
            for unwanted in ("rna_score/pretrained_MLP.h5",
                             "oxdna/build/bin/oxDNA", "_internal/oxpy/__init__.py"):
                entries, _ = self.entries(kind)
                entries["Apta2Switch-Studio/" + unwanted] = False
                with self.subTest(kind=kind, unwanted=unwanted), self.assertRaisesRegex(RuntimeError, "must not be bundled"):
                    checker.check_distribution_archive(self.write(kind, entries))


    def test_model_weights_and_local_configuration_cannot_hide_inside_the_application(self):
        for kind in ("zip", "tar"):
            for unwanted in (".env", ".env.local", ".aws/credentials", ".git/config", ".local_deps/secret",
                             "models/pretrained.h5", "models/network.pt", "models/network.onnx"):
                entries, runtime = self.entries(kind)
                entries[runtime + unwanted] = False
                with self.subTest(kind=kind, unwanted=unwanted), self.assertRaisesRegex(RuntimeError, "must not be bundled"):
                    checker.check_distribution_archive(self.write(kind, entries))

    def test_every_external_worker_dependency_must_be_shipped_as_source(self):
        for kind in ("zip", "tar"):
            for missing in WORKER_FILES:
                entries, runtime = self.entries(kind)
                entries.pop(runtime + "aptaswitch_core/" + missing)
                with self.subTest(kind=kind, missing=missing), self.assertRaisesRegex(RuntimeError, "oxDNA workers"):
                    checker.check_distribution_archive(self.write(kind, entries))

    def test_spec_worker_sources_import_without_the_application_or_checkout(self):
        # Discover actual worker source files from the spec and verify that
        # they remain independent from the application's frozen modules.
        tree = ast.parse((ROOT / "packaging/aptaswitch-studio.spec").read_text())
        sources = []
        for node in ast.walk(tree):
            if isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute) and node.func.attr == "join":
                values = [arg.value for arg in node.args[1:] if isinstance(arg, ast.Constant)]
                if len(values) == 3 and values[:2] == ["src", "aptaswitch_core"] and values[2].endswith(".py"):
                    sources.append(values[2])
        self.assertTrue(sources)
        for name in sources:
            shutil.copyfile(ROOT / "src/aptaswitch_core" / name, self.work / name)
        code = (
            "import sys; sys.path.insert(0, sys.argv[1]); "
            "import oxdna_ffs_worker, oxdna_kinetics_worker, oxdna_guided_worker, oxdna_native_ffs, oxdna_free_md_worker, oxdna_free_md_metrics; "
            "assert 'aptaswitch_core' not in sys.modules"
        )
        result = subprocess.run([sys.executable, "-I", "-B", "-c", code, str(self.work)],
                                cwd=self.work, capture_output=True, text=True, timeout=30)
        self.assertEqual(result.returncode, 0, result.stderr)

    def test_archive_links_must_remain_inside_the_application(self):
        for kind in ("zip", "tar"):
            for target in ("/private/native", "../../native", "..\\external"):
                entries, _ = self.entries(kind)
                entries["Apta2Switch-Studio/library-link"] = target
                with self.subTest(kind=kind, target=target), self.assertRaisesRegex(RuntimeError, "link escapes"):
                    checker.check_distribution_archive(self.write(kind, entries))
            entries, runtime = self.entries(kind)
            entries[runtime + "library-alias"] = "thermodynamics/" + ("libaptaswitch_hybrid.dylib" if kind == "zip" else "libaptaswitch_hybrid.so")
            checker.check_distribution_archive(self.write(kind, entries))

    def test_user_results_must_not_ship(self):
        for kind in ("zip", "tar"):
            for unwanted in ("runs/extensions/private-result.json", "designs/My design/design.json"):
                entries, _ = self.entries(kind)
                entries["Apta2Switch-Studio/" + unwanted] = False
                with self.subTest(kind=kind, unwanted=unwanted), self.assertRaisesRegex(RuntimeError, "User runs"):
                    checker.check_distribution_archive(self.write(kind, entries))

    def test_empty_user_folders_must_remain_directories(self):
        for kind in ("zip", "tar"):
            entries, _ = self.entries(kind)
            entries["Apta2Switch-Studio/runs/designs"] = False
            with self.subTest(kind=kind), self.assertRaisesRegex(RuntimeError, "supporting directory"):
                checker.check_distribution_archive(self.write(kind, entries))

    def test_foreign_top_level_or_traversing_path_is_rejected(self):
        for kind in ("zip", "tar"):
            for bad_path in ("foreign/private.txt", "Apta2Switch-Studio/../../private.txt", "/absolute/file",
                             "Apta2Switch-Studio/dir\\private.txt"):
                entries, _ = self.entries(kind)
                entries[bad_path] = False
                with self.subTest(kind=kind, bad_path=bad_path), self.assertRaisesRegex(RuntimeError, "one Apta2Switch"):
                    checker.check_distribution_archive(self.write(kind, entries))


if __name__ == "__main__":
    unittest.main()
