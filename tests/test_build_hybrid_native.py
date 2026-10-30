"""Build metadata and publication checks without compiling or loading a library."""
import hashlib
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch


ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "build_hybrid_native", ROOT / "scripts/build_hybrid_native.py")
builder = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(builder)
AUDITED_PARAMS = "0163c59586c5389acc0b1284a97fda7fd0652716b1d26548466e5df735956cf9"


class HybridBuildManifestTests(unittest.TestCase):
    def setUp(self):
        self.root = Path(self.enterContext(tempfile.TemporaryDirectory())).resolve()
        self.vienna = self.root / "ViennaRNA-2.7.2"
        self.params = self.vienna / "src/ViennaRNA/params/params.c"
        self.params.parent.mkdir(parents=True)
        self.params.write_bytes(b"audited parameter fixture")
        self.config = self.params.parent.parent / "vrna_config.h"
        self.config.write_text('#define VRNA_VERSION  "2.7.2"\n')
        self.static = self.params.parent.parent / "libRNA.a"
        self.static.write_bytes(b"static library fixture")
        self.output = self.root / "runtime/libaptaswitch_hybrid.dylib"
        self.commands, self.exports = [], []
        self.mutate = lambda: None
        original_digest = builder.digest

        def fixture_digest(path):
            if path == self.params and path.read_bytes() == b"audited parameter fixture":
                return AUDITED_PARAMS
            return original_digest(path)

        def compile_fixture(command, **kwargs):
            self.assertTrue(kwargs["check"])
            self.commands.append(command)
            target = Path(command[command.index("-o") + 1])
            self.exports.append((target.parent / "exports.txt").read_text())
            target.write_bytes(b"complete compiled library fixture")
            self.mutate()

        self.enterContext(patch.object(builder, "digest", side_effect=fixture_digest))
        self.enterContext(patch.object(builder.subprocess, "run", side_effect=compile_fixture))
        self.enterContext(patch.object(builder.subprocess, "check_output", return_value="fixture compiler 1\n"))
        self.enterContext(patch.object(builder.platform, "system", return_value="Darwin"))
        self.enterContext(patch.object(builder.platform, "platform", return_value="fixture-Darwin"))
        self.enterContext(patch("sys.argv", ["build_hybrid_native.py", "--vienna-source", str(self.vienna),
                                             "--output", str(self.output)]))
        self.enterContext(patch("builtins.print"))

    def test_build_records_off_solver_header_identity_and_macos_export(self):
        builder.main()
        manifest = json.loads(self.output.with_suffix(self.output.suffix + ".json").read_text())
        self.assertEqual(manifest["vienna_version"], "2.7.2")
        self.assertEqual(manifest["single_rna_solver"], "rnalib2-pf-bpp-sample-mfe-v1")
        self.assertEqual(manifest["single_rna_minimum_convention"], "turner2004-centikcal-salt-rounded-v1")
        self.assertEqual(manifest["connected_dp_storage"], "lazy-inside-v1")
        self.assertEqual(manifest["vienna_config_h_sha256"], hashlib.sha256(self.config.read_bytes()).hexdigest())
        self.assertEqual(manifest["vienna_params_c_sha256"], AUDITED_PARAMS)
        self.assertEqual(manifest["native_library_sha256"], hashlib.sha256(self.output.read_bytes()).hexdigest())
        self.assertIn("_aptadp_create_rna_v1\n", self.exports[0])
        self.assertIn("_aptadp_create_hybrid_v1\n", self.exports[0])

    def test_wrong_missing_or_ambiguous_header_version_fails_before_compile(self):
        for header in ('#define VRNA_VERSION "2.7.1"\n', "/* no version */\n",
                       '#define VRNA_VERSION "2.7.2"\n#define VRNA_VERSION "2.7.3"\n'):
            with self.subTest(header=header):
                self.config.write_text(header)
                with self.assertRaisesRegex(ValueError, "RNAlib version 2.7.2"):
                    builder.main()
                self.assertFalse(self.output.exists())
        self.assertEqual(self.commands, [])

    def test_header_or_parameter_change_during_compile_preserves_existing_payload(self):
        self.output.parent.mkdir()
        self.output.write_bytes(b"previous library")
        manifest_path = self.output.with_suffix(self.output.suffix + ".json")
        manifest_path.write_bytes(b"previous manifest")
        for changed in (self.config, self.params):
            with self.subTest(changed=changed.name):
                original = changed.read_bytes()
                self.mutate = lambda: changed.write_bytes(original + b"changed during compilation")
                with self.assertRaisesRegex(RuntimeError, "Build input changed"):
                    builder.main()
                self.assertEqual(self.output.read_bytes(), b"previous library")
                self.assertEqual(manifest_path.read_bytes(), b"previous manifest")
                changed.write_bytes(original)


if __name__ == "__main__":
    unittest.main()
