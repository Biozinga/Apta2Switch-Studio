from __future__ import annotations

import base64
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import tempfile
from types import SimpleNamespace
import traceback
import unittest
from unittest.mock import patch


ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("macos_ci_keychain", ROOT / "scripts/macos_ci_keychain.py")
keychains = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(keychains)
IDENTITY = "Developer ID Application: Fixture Owner (ABC123)"


@unittest.skipIf(os.name == "nt", "macOS keychain permissions require a POSIX filesystem")
class MacCIKeychainTests(unittest.TestCase):
    def setUp(self):
        self.root = Path(self.enterContext(tempfile.TemporaryDirectory())).resolve()
        self.environment = {
            "RUNNER_TEMP": str(self.root), "GITHUB_ENV": str(self.root / "github-env"),
            "APPLE_CERTIFICATE_P12_BASE64": base64.b64encode(b"fixture-private-certificate").decode(),
            "APPLE_CERTIFICATE_PASSWORD": "PRIVATE-P12-PASSWORD", "APPLE_ID": "fixture@example.invalid",
            "APPLE_TEAM_ID": "TEAMFIXTURE", "APPLE_APP_PASSWORD": "PRIVATE-APPLE-PASSWORD",
            "APTASWITCH_MACOS_SIGN_IDENTITY": IDENTITY,
            "APTASWITCH_MACOS_NOTARY_PROFILE": "aptaswitch-notary",
        }
        self.previous = ["/Users/runner/Library/Keychains/login.keychain-db", "/Library/Keychains/System With Spaces.keychain"]
        self.commands = []
        self.fail = None
        self.identities = f'1) ABCDEF "{IDENTITY}"'
        self.state_path = self.root / keychains.STATE_NAME
        self.enterContext(patch.object(keychains.subprocess, "run", side_effect=self.run_command))

    def run_command(self, command, **kwargs):
        self.commands.append(command)
        self.assertEqual(kwargs, {"check": True, "capture_output": True, "text": True})
        if self.fail and self.fail(command):
            raise subprocess.CalledProcessError(1, command, output="PRIVATE-APPLE-PASSWORD", stderr="PRIVATE-P12-PASSWORD")
        if command == ["security", "list-keychains", "-d", "user"]:
            return SimpleNamespace(stdout="\n".join(f'    "{path}"' for path in self.previous))
        if command[:2] == ["security", "create-keychain"]:
            self.assertTrue(self.state_path.is_file(), "Recovery state must precede any mutation")
            self.assertEqual(self.state_path.stat().st_mode & 0o777, 0o600)
            Path(command[-1]).write_bytes(b"fixture keychain")
        if command[:2] == ["security", "import"]:
            source = Path(command[2])
            self.assertEqual(source.stat().st_mode & 0o777, 0o600)
            self.assertEqual(source.parent.stat().st_mode & 0o777, 0o700)
            self.assertEqual(source.read_bytes(), b"fixture-private-certificate")
        if command[:2] == ["security", "delete-keychain"]:
            Path(command[-1]).unlink()
        if command[:2] == ["security", "find-identity"]:
            return SimpleNamespace(stdout=self.identities)
        return SimpleNamespace(stdout="")

    def test_prepare_orders_calls_and_preserves_all_existing_keychains(self):
        path = keychains.prepare(self.environment)
        self.assertEqual([command[:2] for command in self.commands], [
            ["security", "list-keychains"], ["security", "create-keychain"],
            ["security", "set-keychain-settings"], ["security", "unlock-keychain"],
            ["security", "import"], ["security", "set-key-partition-list"],
            ["security", "list-keychains"], ["security", "find-identity"], ["xcrun", "notarytool"],
        ])
        self.assertEqual(self.commands[6], ["security", "list-keychains", "-d", "user", "-s", *self.previous, str(path)])
        self.assertEqual(self.commands[-1][self.commands[-1].index("--keychain") + 1], str(path))
        self.assertEqual(Path(self.environment["GITHUB_ENV"]).read_text(), f"APTASWITCH_MACOS_NOTARY_KEYCHAIN={path}\n")
        self.assertEqual(json.loads(self.state_path.read_text())["previous_search_list"], self.previous)
        for value in ("PRIVATE-P12-PASSWORD", "PRIVATE-APPLE-PASSWORD", "fixture-private-certificate"):
            self.assertNotIn(value, self.state_path.read_text())
            self.assertNotIn(value, Path(self.environment["GITHUB_ENV"]).read_text())
        self.assertFalse(list(self.root.glob("aptaswitch-certificate-*")))

    def test_cleanup_restores_exact_search_list_and_deletes_only_owned_keychain(self):
        path = keychains.prepare(self.environment)
        self.commands.clear()
        keychains.cleanup({"RUNNER_TEMP": str(self.root)})
        self.assertEqual(self.commands, [
            ["security", "list-keychains", "-d", "user", "-s", *self.previous],
            ["security", "delete-keychain", str(path)],
        ])
        self.assertFalse(path.exists())
        self.assertFalse(self.state_path.exists())
        self.commands.clear()
        keychains.cleanup({"RUNNER_TEMP": str(self.root)})
        self.assertEqual(self.commands, [])

    def test_import_failure_retains_recovery_state_and_never_exposes_credentials(self):
        self.fail = lambda command: command[:2] == ["security", "import"]
        with self.assertRaisesRegex(RuntimeError, "importing the Developer ID"):
            try:
                keychains.prepare(self.environment)
            except RuntimeError as exc:
                self.assertNotIn("PRIVATE-P12-PASSWORD", str(exc))
                self.assertNotIn("PRIVATE-APPLE-PASSWORD", traceback.format_exc())
                raise
        self.assertTrue(self.state_path.is_file())
        self.assertFalse(list(self.root.glob("aptaswitch-certificate-*")))
        self.fail = None
        keychains.cleanup(self.environment)
        self.assertFalse(self.state_path.exists())
        self.assertFalse(list(self.root.glob("*.keychain-db")))

    def test_missing_or_invalid_identity_never_records_notary_credentials(self):
        del self.environment["APPLE_ID"]
        with self.assertRaisesRegex(RuntimeError, "Missing signing configuration: APPLE_ID"):
            keychains.prepare(self.environment)
        self.assertEqual(self.commands, [])
        self.environment["APPLE_ID"] = "fixture@example.invalid"
        self.identities = "0 valid identities found"
        with self.assertRaisesRegex(RuntimeError, "does not provide"):
            keychains.prepare(self.environment)
        self.assertFalse(any(command[0] == "xcrun" for command in self.commands))
        keychains.cleanup(self.environment)

    def test_invalid_base64_creates_no_keychain_or_state(self):
        self.environment["APPLE_CERTIFICATE_P12_BASE64"] = "PRIVATE invalid %%%"
        with self.assertRaisesRegex(RuntimeError, "not a valid"):
            keychains.prepare(self.environment)
        self.assertFalse(self.state_path.exists())
        self.assertEqual(self.commands, [])

    def test_existing_state_cannot_be_overwritten(self):
        path = keychains.prepare(self.environment)
        before = self.state_path.read_bytes()
        self.commands.clear()
        with self.assertRaisesRegex(RuntimeError, "already exists"):
            keychains.prepare(self.environment)
        self.assertEqual(self.state_path.read_bytes(), before)
        self.assertTrue(path.is_file())
        self.assertFalse(any(command[1] == "create-keychain" for command in self.commands))

    def test_cleanup_refuses_foreign_keychain_path_or_symlink(self):
        path = keychains.prepare(self.environment)
        original = json.loads(self.state_path.read_text())
        invalid = dict(original, keychain="/Users/runner/Library/Keychains/login.keychain-db")
        self.state_path.write_text(json.dumps(invalid))
        self.commands.clear()
        with self.assertRaisesRegex(RuntimeError, "invalid"):
            keychains.cleanup(self.environment)
        self.assertEqual(self.commands, [])
        self.state_path.write_text(json.dumps(original))
        path.unlink()
        path.symlink_to(self.root / "foreign")
        with self.assertRaisesRegex(RuntimeError, "invalid"):
            keychains.cleanup(self.environment)
        self.assertEqual(self.commands, [])

    def test_cleanup_failure_keeps_record_for_retry(self):
        path = keychains.prepare(self.environment)
        self.fail = lambda command: command[:2] == ["security", "delete-keychain"]
        with self.assertRaisesRegex(RuntimeError, "removing"):
            keychains.cleanup(self.environment)
        self.assertTrue(path.exists())
        self.assertTrue(self.state_path.exists())
        self.fail = None
        keychains.cleanup(self.environment)
        self.assertFalse(path.exists())
        self.assertFalse(self.state_path.exists())


if __name__ == "__main__":
    unittest.main()
