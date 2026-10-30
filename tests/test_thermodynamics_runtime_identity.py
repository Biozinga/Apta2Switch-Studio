"""Binary provenance is independent of the unchanged upstream version string."""
import hashlib
import json
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest

from aptaswitch_core.thermodynamics.runtime_identity import identify_vienna_runtime


class RuntimeIdentityTests(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.binary = Path(self.directory.name, "_RNA.so")
        self.binary.write_bytes(b"independent synthetic native build A")
        self.native = SimpleNamespace(__version__="2.7.2",
                                      _RNA=SimpleNamespace(__file__=str(self.binary)))

    def manifest(self, **updates):
        values = {
            "upstream_version": "2.7.2", "patch_id": "aptaswitch-salt-loop-31-v1",
            "upstream_archive_sha256": "1ab5f4a4f76fc85a2243546088e45f5d85f2d7a56cc656e969b005cce9bfab5f",
            "original_params_c_sha256": "a6da5b5038e07c6b4b0b843429115db84b6bcb8b48e51074f3180d2842358e7a",
            "patch_sha256": "700a7e70eb0aa2bc1b1e7778f827bd46ecaabba9ce4ab299b1919485ea424a7e",
            "expected_patched_params_c_sha256": "0163c59586c5389acc0b1284a97fda7fd0652716b1d26548466e5df735956cf9",
            "status": "built", "native_module_sha256": hashlib.sha256(self.binary.read_bytes()).hexdigest(),
        }
        values.update(updates)
        path = Path(self.directory.name, "build.json")
        path.write_text(json.dumps(values))
        return path

    def test_absent_manifest_does_not_claim_unmodified_release(self):
        identity = identify_vienna_runtime(self.native)
        self.assertIsNone(identity.patch_id)
        self.assertIn("patch=unrecorded", identity.version_label)
        self.assertEqual(identity.module_sha256, hashlib.sha256(self.binary.read_bytes()).hexdigest())

    def test_two_binaries_with_same_version_have_distinct_identity(self):
        before = identify_vienna_runtime(self.native)
        self.binary.write_bytes(b"independent synthetic native build B")
        after = identify_vienna_runtime(self.native)
        self.assertEqual(before.upstream_version, after.upstream_version)
        self.assertNotEqual(before.version_label, after.version_label)

    def test_audited_manifest_binds_recipe_to_loaded_binary(self):
        path = self.manifest()
        identity = identify_vienna_runtime(self.native, path)
        self.assertEqual(identity.patch_id, "aptaswitch-salt-loop-31-v1")
        self.assertEqual(identity.manifest_sha256, hashlib.sha256(path.read_bytes()).hexdigest())
        self.assertIn(identity.module_sha256, identity.version_label)

    def test_manifest_for_another_binary_is_rejected(self):
        path = self.manifest(native_module_sha256="0" * 64)
        with self.assertRaisesRegex(ValueError, "loaded ViennaRNA binary"):
            identify_vienna_runtime(self.native, path)

    def test_changed_recipe_or_unfinished_build_is_rejected(self):
        for updates in ({"patch_id": "some-other-fix"}, {"status": "prepared"},
                        {"patch_sha256": "0" * 64}, {"upstream_version": "2.7.3"}):
            with self.subTest(updates=updates), self.assertRaises(ValueError):
                identify_vienna_runtime(self.native, self.manifest(**updates))

    def test_invalid_manifests_and_missing_binaries_fail_explicitly(self):
        path = Path(self.directory.name, "invalid.json")
        for data in (b"[1,2,3]", b"not json", b"null"):
            path.write_bytes(data)
            with self.assertRaises(ValueError):
                identify_vienna_runtime(self.native, path)
        self.binary.unlink()
        with self.assertRaises(RuntimeError):
            identify_vienna_runtime(self.native)
        with self.assertRaises(RuntimeError):
            identify_vienna_runtime(SimpleNamespace(__version__="2.7.2"))
        with self.assertRaises(RuntimeError):
            identify_vienna_runtime(SimpleNamespace(__version__="2.7.3"))


if __name__ == "__main__":
    unittest.main()
