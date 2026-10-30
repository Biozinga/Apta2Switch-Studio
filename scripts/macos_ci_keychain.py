"""Prepare and remove a dedicated CI keychain without changing existing identities."""
from __future__ import annotations

import argparse
import base64
import binascii
import json
import os
from pathlib import Path
import re
import secrets
import shlex
import subprocess
import sys
import tempfile
import uuid


STATE_NAME = "aptaswitch-keychain-state.json"
REQUIRED = ("APPLE_CERTIFICATE_P12_BASE64", "APPLE_CERTIFICATE_PASSWORD", "APPLE_ID",
            "APPLE_TEAM_ID", "APPLE_APP_PASSWORD", "APTASWITCH_MACOS_SIGN_IDENTITY")


def run(command: list[str], phase: str) -> str:
    """Security/notarytool errors can repeat secret argv; never relay their output."""
    try:
        result = subprocess.run(command, check=True, capture_output=True, text=True)
    except (OSError, subprocess.CalledProcessError):
        raise RuntimeError(f"macOS signing setup failed during {phase}; sensitive command output was suppressed") from None
    return result.stdout


def runner_directory(environment: dict[str, str]) -> Path:
    raw = environment.get("RUNNER_TEMP", "")
    if not raw or "\n" in raw or "\r" in raw:
        raise RuntimeError("RUNNER_TEMP must identify the CI runner's temporary directory")
    root = Path(raw).resolve()
    if not root.is_dir():
        raise RuntimeError("RUNNER_TEMP does not exist")
    return root


def prepare(environment: dict[str, str]) -> Path:
    root = runner_directory(environment)
    missing = [name for name in REQUIRED if not environment.get(name)]
    if missing:
        raise RuntimeError("Missing signing configuration: " + ", ".join(missing))
    if not environment["APTASWITCH_MACOS_SIGN_IDENTITY"].startswith("Developer ID Application: "):
        raise RuntimeError("The signing identity must be a Developer ID Application certificate")
    try:
        certificate = base64.b64decode(environment["APPLE_CERTIFICATE_P12_BASE64"], validate=True)
        if not certificate:
            raise ValueError("empty")
    except (ValueError, binascii.Error):
        raise RuntimeError("APPLE_CERTIFICATE_P12_BASE64 is not a valid nonempty base64 certificate") from None
    profile = environment.get("APTASWITCH_MACOS_NOTARY_PROFILE", "aptaswitch-notary")
    if not profile:
        raise RuntimeError("APTASWITCH_MACOS_NOTARY_PROFILE must not be empty")
    previous = shlex.split(run(["security", "list-keychains", "-d", "user"], "reading the keychain search list"))
    token = uuid.uuid4().hex
    keychain = root / f"aptaswitch-signing-{token}.keychain-db"
    state_path = root / STATE_NAME
    state = {"schema_version": 1, "token": token, "keychain": str(keychain), "previous_search_list": previous}
    try:
        descriptor = os.open(state_path, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o600)
    except FileExistsError:
        raise RuntimeError("A signing keychain state already exists; run cleanup before preparing another") from None
    # Save recovery information before create-keychain or any search-list change.
    with os.fdopen(descriptor, "w") as stream:
        json.dump(state, stream)
        stream.flush()
        os.fsync(stream.fileno())
    password = secrets.token_urlsafe(32)
    run(["security", "create-keychain", "-p", password, str(keychain)], "creating the temporary keychain")
    run(["security", "set-keychain-settings", "-lut", "21600", str(keychain)], "configuring the temporary keychain")
    run(["security", "unlock-keychain", "-p", password, str(keychain)], "unlocking the temporary keychain")
    with tempfile.TemporaryDirectory(prefix="aptaswitch-certificate-", dir=root) as folder:
        source = Path(folder) / "identity.p12"
        descriptor = os.open(source, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o600)
        with os.fdopen(descriptor, "wb") as stream:
            stream.write(certificate)
        run(["security", "import", str(source), "-k", str(keychain),
             "-P", environment["APPLE_CERTIFICATE_PASSWORD"], "-T", "/usr/bin/codesign", "-T", "/usr/bin/security"],
            "importing the Developer ID certificate")
    run(["security", "set-key-partition-list", "-S", "apple-tool:,apple:,codesign:", "-s", "-k", password, str(keychain)],
        "authorizing the signing tools")
    run(["security", "list-keychains", "-d", "user", "-s", *previous, str(keychain)], "adding the temporary keychain")
    identities = run(["security", "find-identity", "-v", "-p", "codesigning", str(keychain)], "checking the imported identity")
    if f'"{environment["APTASWITCH_MACOS_SIGN_IDENTITY"]}"' not in identities:
        raise RuntimeError("The imported certificate does not provide the requested Developer ID signing identity")
    run(["xcrun", "notarytool", "store-credentials", profile, "--keychain", str(keychain),
         "--apple-id", environment["APPLE_ID"], "--team-id", environment["APPLE_TEAM_ID"],
         "--password", environment["APPLE_APP_PASSWORD"]], "storing the notarization profile")
    if environment.get("GITHUB_ENV"):
        with Path(environment["GITHUB_ENV"]).open("a") as stream:
            stream.write(f"APTASWITCH_MACOS_NOTARY_KEYCHAIN={keychain}\n")
    return keychain


def cleanup(environment: dict[str, str]) -> None:
    root = runner_directory(environment)
    state_path = root / STATE_NAME
    if not state_path.exists() and not state_path.is_symlink():
        return
    if state_path.is_symlink() or not state_path.is_file():
        raise RuntimeError("Refusing an unsafe signing keychain state")
    metadata = state_path.stat()
    if metadata.st_mode & 0o077 or (hasattr(os, "getuid") and metadata.st_uid != os.getuid()):
        raise RuntimeError("The signing keychain state must be private and owned by the current user")
    try:
        state = json.loads(state_path.read_text())
        token = state["token"]
        if state["schema_version"] != 1 or not isinstance(token, str) or not re.fullmatch(r"[0-9a-f]{32}", token):
            raise ValueError("invalid token")
        keychain = root / f"aptaswitch-signing-{token}.keychain-db"
        previous = state["previous_search_list"]
        if (state["keychain"] != str(keychain) or keychain.is_symlink()
                or not isinstance(previous, list)
                or any(not isinstance(item, str) or not item or "\x00" in item for item in previous)):
            raise ValueError("unsafe path")
    except (KeyError, TypeError, ValueError):
        raise RuntimeError("Refusing an invalid signing keychain state") from None
    errors = []
    try:
        run(["security", "list-keychains", "-d", "user", "-s", *previous], "restoring the keychain search list")
    except RuntimeError as exc:
        errors.append(str(exc))
    if keychain.exists():
        try:
            run(["security", "delete-keychain", str(keychain)], "removing the temporary keychain")
        except RuntimeError as exc:
            errors.append(str(exc))
    if errors:
        # Keep the recovery record so the always() cleanup can be retried.
        raise RuntimeError("; ".join(errors))
    state_path.unlink()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("prepare", "cleanup"))
    args = parser.parse_args()
    try:
        (prepare if args.action == "prepare" else cleanup)(dict(os.environ))
    except (RuntimeError, OSError) as exc:
        # RuntimeError text is constructed above and contains no credential or
        # captured tool output. Filesystem errors are deliberately generic.
        message = str(exc) if isinstance(exc, RuntimeError) else "Filesystem access failed while managing the signing keychain"
        print(message, file=sys.stderr)
        raise SystemExit(1) from None
    print(f"macOS CI signing keychain: {args.action} completed")


if __name__ == "__main__":
    main()
