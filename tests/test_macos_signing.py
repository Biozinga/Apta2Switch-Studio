from __future__ import annotations

import importlib.util
import json
from pathlib import Path
import plistlib
import subprocess
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch


ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("sign_macos_app", ROOT / "scripts/sign_macos_app.py")
signer = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = signer
SPEC.loader.exec_module(signer)
IDENTITY = "Developer ID Application: Example Maintainer (EXAMPLE123)"
RELEASE = signer.SigningConfig("release", IDENTITY, "aptaswitch-notary")
LOCAL = signer.SigningConfig("local", "-", "")
MACHO = b"\xcf\xfa\xed\xfe\x0c\x00\x00\x01" + b"\x00" * 56


class MacSigningTests(unittest.TestCase):
    def setUp(self):
        self.root = Path(self.enterContext(tempfile.TemporaryDirectory())).resolve()
        self.app = self.root / "Apta2Switch-Studio.app"
        self.main = self.write("Contents/MacOS/Apta2Switch-Studio", MACHO)
        self.framework = self.app / "Contents/Frameworks/Python.framework"
        self.framework_binary = self.write("Contents/Frameworks/Python.framework/Versions/3.11/Python", MACHO)
        self.extension = self.write("Contents/Frameworks/numpy/core.so", MACHO)
        self.commands = []
        self.sealed_manifests = None
        self.entitlements = {}
        self.status = "Accepted"
        self.failure = None
        self.enterContext(patch("platform.system", return_value="Darwin"))
        self.enterContext(patch("platform.machine", return_value="arm64"))

        runtime = self.app / "Contents/Resources/oxdna_runtime"
        self.worker = self.write(str((runtime / "worker/aptaswitch-oxdna-worker").relative_to(self.app)), MACHO)
        self.engine = self.write(str((runtime / "worker/_internal/bin/oxDNA").relative_to(self.app)), MACHO)
        self.non_native = self.write(str((runtime / "sources/oxDNA-source.tar.gz").relative_to(self.app)), b"fixture source")
        for name in ("LICENSE.txt", "parameters/dna.txt", "parameters/rna.txt", "parameters/na.txt",
                     "sources/BUILDING.txt", "sources/prepare_bundled_science.py",
                     "sources/prepare_oxdna_ffs_runtime.py", "sources/run_oxdna_worker.py"):
            self.write(str((runtime / name).relative_to(self.app)), b"fixture")
        self.manifest(runtime, {"schema_version": 1, "kind": "oxdna", "files": self.hashes(runtime),
                                "worker": "worker/aptaswitch-oxdna-worker", "binary": "worker/_internal/bin/oxDNA",
                                "parameters": {kind: "parameters/" + kind.lower() + ".txt" for kind in ("DNA", "RNA", "NA")},
                                "source_revision": "fixture", "license": "GPL-3.0"})
        hybrid = self.app / "Contents/Resources/thermodynamics"
        self.hybrid = self.write(str((hybrid / "libaptaswitch_hybrid.dylib").relative_to(self.app)), MACHO)
        self.hybrid_manifest = hybrid / "libaptaswitch_hybrid.dylib.json"
        self.hybrid_manifest.write_text(json.dumps({"status": "built", "abi": 1,
            "model": "experimental-Sugimoto95-Vienna04-loopmean",
            "native_library_sha256": signer.sha256(self.hybrid)}))
        for name in ("ViennaRNA-COPYING.txt", "THIRD-PARTY-NOTICES.txt"):
            self.write(str((hybrid / name).relative_to(self.app)), b"fixture license notice")

    def write(self, name, data):
        path = self.app / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(data)
        path.chmod(0o755)
        return path

    def manifest(self, root, data):
        (root / "runtime.json").write_text(json.dumps(data))

    def hashes(self, root):
        return {path.relative_to(root).as_posix(): signer.sha256(path) for path in root.rglob("*") if path.is_file()}

    def run_command(self, command, **kwargs):
        self.assertTrue(kwargs["check"])
        self.commands.append(command)
        if self.failure and self.failure(command):
            raise subprocess.CalledProcessError(1, command)
        if command[:2] == ["codesign", "--force"]:
            target = Path(command[-1])
            if "--entitlements" in command:
                self.entitlements[target] = plistlib.loads(Path(command[command.index("--entitlements") + 1]).read_bytes())
            if target.is_file():
                target.write_bytes(target.read_bytes() + b"fake signature")
            elif target == self.app:
                self.sealed_manifests = {path: path.read_bytes() for path in self.app.rglob("runtime.json")}
                self.sealed_manifests[self.hybrid_manifest] = self.hybrid_manifest.read_bytes()
        if command[:3] == ["xcrun", "notarytool", "submit"]:
            self.assertTrue(Path(command[3]).parent.is_dir())
            return SimpleNamespace(stdout=json.dumps({"status": self.status, "id": "submission-fixture"}))
        return SimpleNamespace(stdout=f'1) ABCDEF "{IDENTITY}"')

    def test_release_requires_explicit_developer_identity_and_profile(self):
        for environment in ({}, {"APTASWITCH_MACOS_SIGN_IDENTITY": "-"},
                            {"APTASWITCH_MACOS_SIGN_IDENTITY": IDENTITY}):
            with self.subTest(environment=environment), self.assertRaisesRegex(RuntimeError, "release requires"):
                signer.configuration(environment)
        self.assertEqual(signer.configuration({"APTASWITCH_MACOS_SIGNING": "local"}), LOCAL)
        self.assertEqual(signer.configuration({"APTASWITCH_MACOS_SIGN_IDENTITY": IDENTITY,
                                              "APTASWITCH_MACOS_NOTARY_PROFILE": "aptaswitch-notary"}), RELEASE)

    def test_preflight_checks_identity_without_uploading(self):
        with patch.object(signer.subprocess, "run", side_effect=self.run_command):
            signer.preflight(RELEASE)
        self.assertEqual([command[0] for command in self.commands], ["security", "xcrun", "xcrun"])
        self.assertTrue(all("submit" not in command for command in self.commands))
        with patch.object(signer.subprocess, "run", return_value=SimpleNamespace(stdout="0 valid identities")):
            with self.assertRaisesRegex(RuntimeError, "not usable"):
                signer.preflight(RELEASE)

    def test_all_native_files_signed_before_containers_and_manifests_sealed_last(self):
        original = {root: data for root, data, _ in signer.runtime_manifests(self.app)}
        binaries = signer.native_files(self.app)
        with patch.object(signer.subprocess, "run", side_effect=self.run_command):
            signer.sign_app(self.app, RELEASE)
        signed = [Path(command[-1]) for command in self.commands if command[:2] == ["codesign", "--force"]]
        self.assertEqual(signed[:len(binaries)], binaries)
        self.assertLess(signed.index(self.framework_binary), signed.index(self.framework))
        self.assertEqual(signed[-1], self.app)
        self.assertNotIn(self.non_native, signed)
        self.assertIn(self.worker, signed)
        self.assertIn(self.engine, signed)
        for command in self.commands[:-1]:
            self.assertNotIn("--deep", command)
            self.assertIn("runtime", command)
            self.assertIn("--timestamp", command)
        self.assertIn(signer.DEVELOPER_ID_REQUIREMENT, self.commands[-1])
        self.assertEqual(self.entitlements, {})
        for path, content in self.sealed_manifests.items():
            self.assertEqual(path.read_bytes(), content)
        for root, data, field in signer.runtime_manifests(self.app):
            self.assertEqual(data["unsigned_" + field], original[root][field])
        science = json.loads((self.worker.parent.parent / "runtime.json").read_text())
        self.assertEqual(science["files"]["worker/aptaswitch-oxdna-worker"], signer.sha256(self.worker))
        self.assertNotEqual(science["files"]["worker/aptaswitch-oxdna-worker"],
                            science["unsigned_files"]["worker/aptaswitch-oxdna-worker"])
        self.assertEqual(science["files"]["sources/oxDNA-source.tar.gz"],
                         science["unsigned_files"]["sources/oxDNA-source.tar.gz"])
        hybrid = json.loads(self.hybrid_manifest.read_text())
        self.assertEqual(hybrid["native_library_sha256"], signer.sha256(self.hybrid))
        self.assertNotEqual(hybrid["native_library_sha256"], hybrid["unsigned_native_library_sha256"])
        self.assertEqual(hybrid["code_signing"]["identity"], IDENTITY)
        self.assertIn(self.hybrid, signed)
        self.assertEqual(self.sealed_manifests[self.hybrid_manifest], self.hybrid_manifest.read_bytes())

    def test_local_mode_signs_without_claiming_notarization(self):
        with patch.object(signer.subprocess, "run", side_effect=self.run_command):
            signer.sign_app(self.app, LOCAL)
            with self.assertRaisesRegex(RuntimeError, "cannot be notarized"):
                signer.notarize_app(self.app, LOCAL)
        self.assertFalse(any(command[0] == "xcrun" for command in self.commands))
        self.assertFalse(any("--timestamp" in command for command in self.commands))
        self.assertEqual(self.entitlements, {})

    def test_changed_payload_is_rejected_before_any_signing(self):
        for path in (self.worker, self.engine):
            original = path.read_bytes()
            path.write_bytes(original + b"tampered")
            try:
                with self.subTest(payload=path.name), patch.object(signer.subprocess, "run") as command:
                    with self.assertRaisesRegex(ValueError, "checksum mismatch"):
                        signer.sign_app(self.app, RELEASE)
                command.assert_not_called()
            finally:
                path.write_bytes(original)

    def test_changed_native_hybrid_payload_is_rejected_before_any_signing(self):
        self.hybrid.write_bytes(MACHO + b"tampered")
        with patch.object(signer.subprocess, "run") as command:
            with self.assertRaisesRegex(RuntimeError, "does not match"):
                signer.sign_app(self.app, RELEASE)
        command.assert_not_called()

    def test_symlinks_never_sign_external_files_or_duplicate_native_files(self):
        alias = self.app / "Contents/Resources/native-alias"
        try:
            alias.symlink_to(self.extension)
        except (OSError, NotImplementedError):
            self.skipTest("The test account cannot create symlinks")
        self.assertNotIn(alias, signer.native_files(self.app))
        outside = self.root / "outside"
        outside.write_bytes(MACHO)
        alias.unlink()
        alias.symlink_to(outside)
        with self.assertRaisesRegex(RuntimeError, "escapes"):
            signer.native_files(self.app)

    def test_signature_failure_stops_before_outer_seal(self):
        self.failure = lambda command: command[-1] == str(self.extension)
        with patch.object(signer.subprocess, "run", side_effect=self.run_command):
            with self.assertRaises(subprocess.CalledProcessError):
                signer.sign_app(self.app, RELEASE)
        self.assertIsNone(self.sealed_manifests)

    def test_notarization_acceptance_precedes_ticket_validation_and_gatekeeper(self):
        with patch.object(signer.subprocess, "run", side_effect=self.run_command):
            signer.notarize_app(self.app, RELEASE)
        operations = [command[:3] for command in self.commands]
        self.assertEqual(operations, [
            ["codesign", "--verify", "--deep"], ["ditto", "-c", "-k"],
            ["xcrun", "notarytool", "submit"], ["xcrun", "stapler", "staple"],
            ["xcrun", "stapler", "validate"], ["codesign", "--verify", "--deep"],
            ["spctl", "--assess", "--type"],
        ])

    def test_ci_notarization_uses_its_dedicated_keychain(self):
        config = signer.SigningConfig("release", IDENTITY, "aptaswitch-notary", "/tmp/fixture.keychain-db")
        with patch.object(signer.subprocess, "run", side_effect=self.run_command):
            signer.notarize_app(self.app, config)
        command = next(command for command in self.commands if command[:3] == ["xcrun", "notarytool", "submit"])
        self.assertEqual(command[command.index("--keychain") + 1], config.keychain)

    def test_rejected_notarization_never_staples_or_assesses(self):
        for status in ("Invalid", "In Progress", None):
            with self.subTest(status=status):
                self.commands.clear()
                self.status = status
                with patch.object(signer.subprocess, "run", side_effect=self.run_command):
                    with self.assertRaisesRegex(RuntimeError, "not accepted"):
                        signer.notarize_app(self.app, RELEASE)
                self.assertFalse(any("stapler" in command or "spctl" in command for command in self.commands))

    def test_malformed_notarization_reply_never_staples(self):
        for output in ("Not JSON", "[]", "null"):
            def command(args, **kwargs):
                result = self.run_command(args, **kwargs)
                if args[:3] == ["xcrun", "notarytool", "submit"]:
                    return SimpleNamespace(stdout=output)
                return result
            with self.subTest(output=output):
                self.commands.clear()
                with patch.object(signer.subprocess, "run", side_effect=command):
                    with self.assertRaisesRegex(RuntimeError, "valid notarization result"):
                        signer.notarize_app(self.app, RELEASE)
                self.assertFalse(any("stapler" in args for args in self.commands))

    def test_ticket_or_gatekeeper_failure_prevents_success(self):
        for failed in (["xcrun", "stapler", "validate"], ["spctl", "--assess", "--type"]):
            with self.subTest(failed=failed):
                self.failure = lambda command: command[:3] == failed
                with patch.object(signer.subprocess, "run", side_effect=self.run_command):
                    with self.assertRaises(subprocess.CalledProcessError):
                        signer.notarize_app(self.app, RELEASE)


if __name__ == "__main__":
    unittest.main()
