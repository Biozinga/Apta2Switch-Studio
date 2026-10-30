"""Sign the assembled app, then notarize and verify it before distribution.

Developer credentials stay in the maintainer's keychain. The explicit local
mode only provides the ad-hoc signatures required by Apple Silicon; it does
not make a downloaded app trusted by Gatekeeper.
"""
from __future__ import annotations

import argparse
from dataclasses import dataclass
import hashlib
import json
import os
from pathlib import Path
import plistlib
import runpy
import subprocess
import tempfile


ROOT = Path(__file__).resolve().parents[1]
MACHO_MAGIC = {b"\xfe\xed\xfa\xce", b"\xce\xfa\xed\xfe", b"\xfe\xed\xfa\xcf",
               b"\xcf\xfa\xed\xfe", b"\xca\xfe\xba\xbe", b"\xbe\xba\xfe\xca",
               b"\xca\xfe\xba\xbf", b"\xbf\xba\xfe\xca"}
DEVELOPER_ID_REQUIREMENT = (
    "anchor apple generic and certificate leaf[field.1.2.840.113635.100.6.1.13] exists"
)


@dataclass(frozen=True)
class SigningConfig:
    mode: str
    identity: str
    profile: str
    keychain: str = ""


def configuration(environment: dict[str, str] | None = None) -> SigningConfig:
    environment = os.environ if environment is None else environment
    mode = environment.get("APTASWITCH_MACOS_SIGNING", "release")
    if mode == "local":
        return SigningConfig(mode, "-", "")
    if mode != "release":
        raise RuntimeError("APTASWITCH_MACOS_SIGNING must be release or local")
    identity = environment.get("APTASWITCH_MACOS_SIGN_IDENTITY", "").strip()
    profile = environment.get("APTASWITCH_MACOS_NOTARY_PROFILE", "").strip()
    if not identity.startswith("Developer ID Application: ") or not profile:
        raise RuntimeError(
            "A release requires APTASWITCH_MACOS_SIGN_IDENTITY (Developer ID Application) "
            "and APTASWITCH_MACOS_NOTARY_PROFILE (notarytool keychain profile). "
            "For a private local build only, set APTASWITCH_MACOS_SIGNING=local."
        )
    return SigningConfig(mode, identity, profile, environment.get("APTASWITCH_MACOS_NOTARY_KEYCHAIN", "").strip())


def preflight(config: SigningConfig) -> None:
    """Fail before downloading/building anything when release credentials are absent."""
    if config.mode == "release":
        result = subprocess.run(["security", "find-identity", "-v", "-p", "codesigning"],
                                check=True, capture_output=True, text=True)
        if f'"{config.identity}"' not in result.stdout:
            raise RuntimeError("The requested Developer ID Application identity is not usable in the keychain")
        for tool in ("notarytool", "stapler"):
            subprocess.run(["xcrun", "--find", tool], check=True, capture_output=True, text=True)


def sha256(path: Path) -> str:
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def native_files(app: Path) -> list[Path]:
    """Inspect file magic, including extensions and binaries collected as data.

    PyInstaller's Resources/Frameworks links refer to the same files. Sign only
    real files, never twice through aliases or outside the app boundary.
    """
    files = []
    for path in app.rglob("*"):
        if path.is_symlink():
            if not path.resolve().is_relative_to(app.resolve()):
                raise RuntimeError(f"App symlink escapes the bundle: {path}")
            continue
        if path.is_file():
            with path.open("rb") as stream:
                if stream.read(4) in MACHO_MAGIC:
                    files.append(path)
    return sorted(files, key=lambda path: (-len(path.parts), str(path)))


def runtime_manifests(app: Path) -> list[tuple[Path, dict, str]]:
    """Validate immutable payloads before any signature changes their hashes."""
    resources = app / "Contents/Resources"
    manifests = []
    for directory in ("oxdna_runtime",):
        root = resources / directory
        validator = runpy.run_path(str(ROOT / "scripts/prepare_bundled_science.py"))["validate_payload"]
        manifests.append((root, validator(root), "files"))
    root = resources / "thermodynamics"
    validator = runpy.run_path(str(ROOT / "scripts/prepare_thermodynamics_runtime.py"))["validate_payload"]
    manifests.append((root, validator(root), "native_library_sha256"))
    return manifests


def refresh_manifests(manifests: list[tuple[Path, dict, str]], config: SigningConfig) -> None:
    for root, manifest, field in manifests:
        manifest.setdefault("unsigned_" + field, manifest[field])
        if field == "native_library_sha256":
            manifest[field] = sha256(root / "libaptaswitch_hybrid.dylib")
        else:
            manifest[field] = {name: sha256(root / name) for name in manifest[field]}
        manifest["code_signing"] = {
            "mode": config.mode, "identity": config.identity,
            "hardened_runtime": config.mode == "release",
        }
        name = "libaptaswitch_hybrid.dylib.json" if field == "native_library_sha256" else "runtime.json"
        (root / name).write_text(json.dumps(manifest, indent=2) + "\n")


def verify_signature(app: Path, config: SigningConfig) -> None:
    command = ["codesign", "--verify", "--deep", "--strict", "--verbose=2"]
    if config.mode == "release":
        command += ["-R", DEVELOPER_ID_REQUIREMENT]
    subprocess.run([*command, str(app)], check=True)


def sign_app(app: Path, config: SigningConfig) -> None:
    app = app.resolve()
    main = app / "Contents/MacOS/Apta2Switch-Studio"
    if app.suffix != ".app" or not main.is_file():
        raise RuntimeError("Expected the assembled Apta2Switch-Studio.app bundle")
    manifests = runtime_manifests(app)
    binaries = native_files(app)
    if main not in binaries:
        raise RuntimeError("The application launcher is not a Mach-O executable")
    base = ["codesign", "--force", "--sign", config.identity]
    if config.mode == "release":
        base += ["--options", "runtime", "--timestamp"]
    # Every native dependency is bundled and signed with this identity.
    # Keep hardened-runtime library validation enabled for release builds.
    def sign(path: Path) -> None:
        subprocess.run([*base, str(path)], check=True)

    for binary in binaries:
        sign(binary)
    # Manifest bytes are resources sealed by the following container
    # signatures. Nothing in the app is rewritten after its outer seal.
    refresh_manifests(manifests, config)
    containers = [path for path in app.rglob("*") if path.is_dir() and not path.is_symlink()
                  and path.suffix in {".framework", ".app", ".xpc", ".bundle"}
                  and any(binary.is_relative_to(path) for binary in binaries)]
    for container in sorted(containers, key=lambda path: (-len(path.parts), str(path))):
        sign(container)
    sign(app)
    verify_signature(app, config)


def notarize_app(app: Path, config: SigningConfig) -> None:
    """Accepted submission, stapled ticket and Gatekeeper acceptance are mandatory."""
    if config.mode != "release":
        raise RuntimeError("A local ad-hoc build cannot be notarized as a release")
    app = app.resolve()
    verify_signature(app, config)
    with tempfile.TemporaryDirectory(prefix="aptaswitch-notarize-") as folder:
        archive = Path(folder) / "Apta2Switch-Studio.zip"
        subprocess.run(["ditto", "-c", "-k", "--sequesterRsrc", "--keepParent", str(app), str(archive)], check=True)
        result = subprocess.run([
            "xcrun", "notarytool", "submit", str(archive), "--keychain-profile", config.profile,
            *(["--keychain", config.keychain] if config.keychain else []),
            "--wait", "--output-format", "json",
        ], check=True, capture_output=True, text=True)
        try:
            report = json.loads(result.stdout)
            if not isinstance(report, dict):
                raise ValueError("Expected a notarization result object")
        except (ValueError, TypeError) as exc:
            raise RuntimeError("Apple did not return a valid notarization result; no release archive is produced") from exc
        if report.get("status") != "Accepted":
            raise RuntimeError(
                f"Apple notarization was not accepted (status={report.get('status')}, id={report.get('id')}). "
                "Inspect the submission with xcrun notarytool log before retrying."
            )
    subprocess.run(["xcrun", "stapler", "staple", str(app)], check=True)
    subprocess.run(["xcrun", "stapler", "validate", str(app)], check=True)
    verify_signature(app, config)
    subprocess.run(["spctl", "--assess", "--type", "execute", "--verbose=2", str(app)], check=True)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("app", nargs="?", type=Path)
    action = parser.add_mutually_exclusive_group()
    action.add_argument("--check-config", action="store_true")
    action.add_argument("--notarize", action="store_true")
    args = parser.parse_args()
    config = configuration()
    if args.check_config:
        preflight(config)
        print(f"macOS signing mode: {config.mode}")
    elif args.app is None:
        parser.error("app is required for signing or notarization")
    elif args.notarize:
        notarize_app(args.app, config)
    else:
        sign_app(args.app, config)


if __name__ == "__main__":
    main()
