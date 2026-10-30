"""Automatic bundled oxDNA discovery retains integrity and physics checks."""
from __future__ import annotations

import hashlib
import os
from pathlib import Path
import subprocess
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

from aptaswitch_core import oxdna_setup as setup


class OxDNASetupTests(unittest.TestCase):
    def setUp(self):
        self.root = Path(self.enterContext(tempfile.TemporaryDirectory())).resolve()
        self.enterContext(patch.object(setup, "runtime_directory", return_value=self.root))
        self.probe = self.enterContext(patch.object(setup, "_probe_oxpy", return_value="3.7"))
        self.native_probe = self.enterContext(patch.object(setup, "_probe_native_ffs", return_value={
            "force_equivalence": True, "native_first_crossing": True}))
        setup._checked_runtime.cache_clear()

    def installation(self):
        import json
        files = {"worker/engine": b"executable", "worker/bin/oxDNA": b"native", "LICENSE.txt": b"GPL3"}
        files.update({"parameters/" + filename: kind.encode() for kind, filename in setup.PARAMETERS.items()})
        for relative, content in files.items():
            path = self.root / relative
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(content)
            path.chmod(0o755)
        metadata = {"schema_version": 1, "kind": "oxdna", "license": "GPL-3.0", "source_revision": "pinned",
                    "worker": "worker/engine", "binary": "worker/bin/oxDNA",
                    "parameters": {kind: "parameters/" + filename for kind, filename in setup.PARAMETERS.items()},
                    "files": {name: hashlib.sha256(content).hexdigest() for name, content in files.items()}}
        (self.root / "runtime.json").write_text(json.dumps(metadata))
        return metadata

    def test_bundled_manifest_preserves_executable_and_chemistry_provenance(self):
        metadata = self.installation()
        runtime = setup.detect_oxdna_runtime()
        self.assertTrue(runtime["bundled"])
        self.assertEqual(runtime["python"], str(self.root / metadata["worker"]))
        self.assertEqual(runtime["binary"], str(self.root / metadata["binary"]))
        self.assertEqual(runtime["oxpy_version"], "3.7")
        for kind, relative in metadata["parameters"].items():
            self.assertEqual(runtime["seq_dep_file_" + kind], str(self.root / relative))
            self.assertEqual(runtime["parameter_sha256"][kind], metadata["files"][relative])
        setup.detect_oxdna_runtime()
        self.probe.assert_called_once()

    def test_user_environment_does_not_redirect_the_bundled_engine(self):
        self.installation()
        with patch.dict(os.environ, {"APTASWITCH_OXDNA_HOME": "/native/user/engine"}):
            self.assertEqual(setup.detect_oxdna_runtime()["home"], str(self.root))

    def test_missing_payload_is_unavailable_but_corrupt_payload_raises(self):
        self.assertIsNone(setup.detect_oxdna_runtime())
        self.installation()
        (self.root / "worker/engine").write_text("corrupt")
        with self.assertRaisesRegex(ValueError, "checksum mismatch"):
            setup.detect_oxdna_runtime()
        self.probe.assert_not_called()

    def test_missing_and_nonexecutable_files_are_rejected_before_import(self):
        metadata = self.installation()
        binary = self.root / metadata["binary"]
        binary.chmod(0o644)
        with self.assertRaises(ValueError):
            setup.validate_oxdna_runtime(str(self.root))
        binary.chmod(0o755)
        (self.root / metadata["parameters"]["NA"]).unlink()
        with self.assertRaises(ValueError):
            setup.validate_oxdna_runtime(str(self.root))
        self.probe.assert_not_called()

    def test_escaping_manifest_path_is_rejected(self):
        import json
        metadata = self.installation()
        metadata["files"]["../private"] = "a" * 64
        (self.root / "runtime.json").write_text(json.dumps(metadata))
        with self.assertRaisesRegex(ValueError, "Unsafe"):
            setup.validate_oxdna_runtime(str(self.root))



class OxpyProbeTests(unittest.TestCase):
    def test_native_ffs_probe_rejects_force_mismatch_and_propagates_verified_capability(self):
        payload = '{"force_equivalence": true, "native_first_crossing": true}'
        with patch.object(setup.subprocess, "run", return_value=SimpleNamespace(returncode=0, stdout=payload)) as run:
            self.assertTrue(setup._probe_native_ffs({"python": "/native/python"})["native_first_crossing"])
        self.assertEqual(run.call_args.args[0][0], "/native/python")
        self.assertIn("--probe", run.call_args.args[0])
        for result in (SimpleNamespace(returncode=1, stdout="", stderr="forces differ"),
                       SimpleNamespace(returncode=0, stdout='{"force_equivalence":false}')):
            with patch.object(setup.subprocess, "run", return_value=result), self.assertRaises(ValueError):
                setup._probe_native_ffs({"python": "/native/python"})

    def test_probe_uses_external_python_environment_and_required_oxpy_api(self):
        with patch.object(setup, "external_process_environment", return_value={"LANG": "C"}), \
             patch.object(setup.subprocess, "run", return_value=SimpleNamespace(returncode=0, stdout="3.7\n")) as run:
            self.assertEqual(setup._probe_oxpy("/local/.venv/bin/python"), "3.7")
        args, options = run.call_args
        self.assertEqual(args[0][:2], ["/local/.venv/bin/python", "-c"])
        self.assertIn("OxpyManager", args[0][2])
        self.assertEqual(options["env"], {"LANG": "C"})
        self.assertEqual(options["timeout"], 15)

    def test_import_failure_and_timeout_never_validate(self):
        with patch.object(setup.subprocess, "run", return_value=SimpleNamespace(returncode=1, stdout="")):
            with self.assertRaises(ValueError):
                setup._probe_oxpy("/missing/python")
        with patch.object(setup.subprocess, "run", side_effect=subprocess.TimeoutExpired("python", 15)):
            with self.assertRaises(ValueError):
                setup._probe_oxpy("/timeout/python")


if __name__ == "__main__":
    unittest.main()
