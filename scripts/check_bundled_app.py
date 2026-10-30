"""Check release archives or launch a binary outside the checkout."""

from __future__ import annotations

import json
import hashlib
import os
from contextlib import contextmanager
from html.parser import HTMLParser
from pathlib import Path, PurePosixPath
import posixpath
import shutil
import stat
import subprocess
import sys
import tarfile
import tempfile
import time
import urllib.error
import urllib.parse
import urllib.request
from zipfile import ZipFile


EXTERNAL_WORKER_FILES = (
    "oxdna_ffs_worker.py", "oxdna_kinetics_worker.py", "oxdna_native_ffs.py",
    "oxdna_kinetics.py", "oxdna_initialization.py", "oxdna_guided_worker.py",
    "oxdna_geometry.py", "oxdna_loop_geometry.py", "oxdna_box.py",
    "oxdna_free_md_worker.py", "oxdna_free_md_metrics.py",
)

OSTIR_SOURCE_FILES = (
    "OSTIR-LICENSE.txt", "NOTICE.txt", "ostir-1.1.3.tar.gz",
    "ViennaRNA-2.7.2.tar.gz", "Apta2Switch-source.tar.gz",
)


def _check_ostir_archive(archive_path: Path, base: str, all_files: set[str]) -> None:
    """Do not ship frozen OSTIR without tables, notices and exact source files."""
    required = {
        base + "ostir/rna_turner1999.par", base + "ostir/rna_turner2004.par",
        base + "OSTIR-1.1.3.dist-info/METADATA", base + "OSTIR-1.1.3.dist-info/LICENSE",
        base + "viennarna-2.7.2.dist-info/METADATA",
        base + "viennarna-2.7.2.dist-info/licenses/COPYING",
        base + "ostir_sources/manifest.json",
        *(base + "ostir_sources/" + name for name in OSTIR_SOURCE_FILES),
    }
    if required - all_files:
        raise RuntimeError(f"Incomplete OSTIR package or corresponding source: {sorted(required - all_files)}")
    if not any(name.startswith(base + "RNA/_RNA") and name.endswith((".so", ".pyd")) for name in all_files):
        raise RuntimeError("The native RNA module required by OSTIR is missing")
    archive_type = ZipFile if archive_path.suffix == ".zip" else tarfile.open
    with archive_type(archive_path) as archive:
        def open_member(name):
            return archive.open(name) if isinstance(archive, ZipFile) else archive.extractfile(name)
        prefix = base + "ostir_sources/"
        with open_member(prefix + "manifest.json") as stream:
            manifest = json.load(stream)
        if (manifest.get("ostir_version") != "1.1.3" or manifest.get("viennarna_version") != "2.7.2"
                or manifest.get("license") != "GPL-3.0-or-later"
                or not set(OSTIR_SOURCE_FILES).issubset(manifest.get("files", {}))):
            raise RuntimeError("OSTIR source/version/licence manifest is incomplete")
        for filename in OSTIR_SOURCE_FILES:
            with open_member(prefix + filename) as stream:
                actual = hashlib.file_digest(stream, "sha256").hexdigest()
            if actual != manifest["files"][filename]:
                raise RuntimeError(f"OSTIR corresponding source differs from manifest: {filename}")


def _check_thermodynamics_archive(archive_path: Path, base: str, all_files: set[str]) -> None:
    runtime = base + "thermodynamics/"
    library = "libaptaswitch_hybrid.dylib" if archive_path.suffix == ".zip" else "libaptaswitch_hybrid.so"
    expected = {runtime + name for name in (library, library + ".json", "ViennaRNA-COPYING.txt", "THIRD-PARTY-NOTICES.txt")}
    if expected - all_files:
        raise RuntimeError(f"Incomplete native thermodynamics runtime or notices: {sorted(expected - all_files)}")
    archive_type = ZipFile if archive_path.suffix == ".zip" else tarfile.open
    with archive_type(archive_path) as archive:
        def open_member(name):
            return archive.open(name) if isinstance(archive, ZipFile) else archive.extractfile(name)
        with open_member(runtime + library + ".json") as stream:
            raw = stream.read(1024 * 1024 + 1)
        if len(raw) > 1024 * 1024:
            raise RuntimeError("The thermodynamics manifest is too large")
        try:
            manifest = json.loads(raw)
        except (UnicodeError, ValueError) as exc:
            raise RuntimeError("Invalid thermodynamics build manifest") from exc
        if (not isinstance(manifest, dict) or manifest.get("status") != "built"
                or manifest.get("model") != "experimental-Sugimoto95-Vienna04-loopmean"
                or type(manifest.get("abi")) is not int or manifest["abi"] != 1):
            raise RuntimeError("The archive does not identify the supported native thermodynamic model")
        digest = hashlib.sha256()
        with open_member(runtime + library) as stream:
            for chunk in iter(lambda: stream.read(1024 * 1024), b""):
                digest.update(chunk)
        if digest.hexdigest() != manifest.get("native_library_sha256"):
            raise RuntimeError("The bundled thermodynamics binary differs from its manifest")
        for notice in ("ViennaRNA-COPYING.txt", "THIRD-PARTY-NOTICES.txt"):
            with open_member(runtime + notice) as stream:
                if not stream.read(1):
                    raise RuntimeError("The thermodynamics runtime has an empty license notice")


@contextmanager
def _open_random_access_archive(archive_path: Path):
    """Decompress tar once before verifying thousands of members out of order."""
    if archive_path.suffix == ".zip":
        with ZipFile(archive_path) as archive:
            yield archive
        return
    # Seeking backwards through gzip otherwise decompresses the same payload
    # for every checksum. Use disk rather than holding the whole app in memory.
    with tarfile.open(archive_path) as source, tempfile.TemporaryFile() as unpacked:
        source.fileobj.seek(0)
        shutil.copyfileobj(source.fileobj, unpacked, length=1024 * 1024)
        unpacked.seek(0)
        with tarfile.open(fileobj=unpacked, mode="r:") as archive:
            yield archive


def _check_science_archive(archive_path: Path, base: str, all_files: set[str]) -> set[str]:
    """Verify each bundled engine, legal material and exact payload hashes."""
    admitted = set()
    with _open_random_access_archive(archive_path) as archive:
        if isinstance(archive, ZipFile):
            members = {item.filename.rstrip("/"): item for item in archive.infolist()}
            links = {name: (archive.read(item).decode("utf-8"), False)
                     for name, item in members.items() if stat.S_ISLNK(item.external_attr >> 16)}
        else:
            members = {item.name.rstrip("/"): item for item in archive.getmembers()}
            links = {name: (item.linkname, item.islnk()) for name, item in members.items()
                     if item.issym() or item.islnk()}

        def read(name):
            seen = set()
            while True:
                if name in seen or len(seen) >= 64:
                    raise RuntimeError("Cyclic scientific runtime archive link")
                seen.add(name)
                parts = PurePosixPath(name).parts
                # Resolve directory components before looking up the final file:
                # Python.framework/Python -> Versions/Current/Python, with
                # Current itself a directory symlink to the concrete version.
                for index in range(1, len(parts) + 1):
                    prefix = "/".join(parts[:index])
                    if prefix not in links:
                        continue
                    target, hardlink = links[prefix]
                    destination = posixpath.normpath(posixpath.join(
                        "" if hardlink else posixpath.dirname(prefix), target, *parts[index:]))
                    if ("\\" in target or PurePosixPath(target).is_absolute()
                            or not destination.startswith("Apta2Switch-Studio/")):
                        raise RuntimeError("Archive link escapes the application folder")
                    name = destination
                    break
                else:
                    member = members[name]
                    if (member.is_dir() if isinstance(archive, ZipFile) else not member.isfile()):
                        raise RuntimeError("Scientific runtime manifest refers to a non-file")
                    stream = archive.open(member) if isinstance(archive, ZipFile) else archive.extractfile(member)
                    with stream:
                        return stream.read()
        for directory, kind in (("oxdna_runtime", "oxdna"),):
            prefix = base + directory + "/"
            try:
                metadata = json.loads(read(prefix + "runtime.json"))
            except (KeyError, ValueError, UnicodeError) as exc:
                raise RuntimeError(f"Missing or invalid bundled {kind} manifest") from exc
            if (not isinstance(metadata, dict) or metadata.get("schema_version") != 1
                    or metadata.get("kind") != kind or not isinstance(metadata.get("files"), dict)
                    or not metadata["files"]):
                raise RuntimeError(f"Invalid bundled {kind} manifest")
            required = {"LICENSE.txt"}
            if kind == "oxdna":
                required.update(("sources/oxDNA-source.tar.gz", "sources/BUILDING.txt",
                                 "sources/prepare_bundled_science.py", "sources/prepare_oxdna_ffs_runtime.py",
                                 "sources/run_oxdna_worker.py", "sources/rebuild_oxdna_worker.py"))
                required.update("sources/" + name for name in EXTERNAL_WORKER_FILES)
                if (metadata.get("license") != "GPL-3.0" or not metadata.get("source_revision")
                        or set(metadata.get("parameters", {})) != {"DNA", "RNA", "NA"}):
                    raise RuntimeError("Missing oxDNA chemistry or GPL source provenance")
                required.update(metadata["parameters"].values())
                required.update((metadata.get("worker"), metadata.get("binary")))
            if None in required or not required.issubset(metadata["files"]):
                raise RuntimeError(f"Incomplete bundled {kind} runtime, sources or notices")
            for relative, digest in metadata["files"].items():
                parts = PurePosixPath(relative).parts
                if not parts or PurePosixPath(relative).is_absolute() or ".." in parts or "\\" in relative:
                    raise RuntimeError("Unsafe scientific runtime manifest path")
                name = prefix + relative
                try:
                    actual = hashlib.sha256(read(name)).hexdigest()
                except KeyError as exc:
                    raise RuntimeError(f"Missing bundled scientific file: {name}") from exc
                if actual != digest:
                    raise RuntimeError(f"Bundled scientific file differs from its manifest: {name}")
                admitted.add(name)
    return admitted


def check_distribution_archive(archive_path: Path) -> None:
    """Require complete offline engines while excluding user data and obsolete engines."""
    root = "Apta2Switch-Studio"
    if archive_path.suffix == ".zip":
        executable = f"{root}/Apta2Switch-Studio.app/Contents/MacOS/Apta2Switch-Studio"
        with ZipFile(archive_path) as archive:
            entries = [(item.filename.rstrip("/"), item.is_dir()) for item in archive.infolist()]
            links = [(item.filename.rstrip("/"), archive.read(item).decode("utf-8"), False)
                     for item in archive.infolist() if stat.S_ISLNK(item.external_attr >> 16)]
    else:
        executable = f"{root}/Apta2Switch-Studio"
        with tarfile.open(archive_path, "r:gz") as archive:
            entries = [(item.name.rstrip("/"), item.isdir()) for item in archive]
            links = [(item.name.rstrip("/"), item.linkname, item.islnk())
                     for item in archive.getmembers() if item.issym() or item.islnk()]
    directories = {f"{root}/{name}" for name in ("runs/designs", "runs/extensions")}
    required = directories | {executable, f"{root}/START-HERE.txt", f"{root}/LICENSE"}
    missing = required - {name for name, _ in entries}
    if missing:
        raise RuntimeError(f"Incomplete application archive: {sorted(missing)}")
    all_files = {name for name, is_dir in entries if not is_dir}
    if required - directories - all_files:
        raise RuntimeError("The application, start instructions and license must be files")
    manifests = [name for name in all_files if name.endswith("/oxdna_runtime/runtime.json")]
    if not manifests:
        raise RuntimeError("The bundled oxDNA runtime is missing from the application archive")
    admitted = set()
    for manifest in manifests:
        base = manifest.removesuffix("oxdna_runtime/runtime.json")
        _check_thermodynamics_archive(archive_path, base, all_files)
        _check_ostir_archive(archive_path, base, all_files)
        admitted.update(_check_science_archive(archive_path, base, all_files))
        required_workers = {base + "aptaswitch_core/" + name for name in EXTERNAL_WORKER_FILES}
        viewer = base + "aptaswitch_studio/assets/3dmol/"
        required_viewer = {viewer + name for name in ("3Dmol-min.js", "LICENSE.txt", "3Dmol-min.js.LICENSE.txt", "manifest.json")}
        if (required_workers | required_viewer) - all_files:
            raise RuntimeError("Incomplete bundled oxDNA workers or offline molecular viewer")
    for name, is_directory in entries:
        parts = PurePosixPath(name).parts
        if not parts or parts[0] not in {root, "__MACOSX"} or ".." in parts or "\\" in name:
            raise RuntimeError("The archive must contain one Apta2Switch-Studio folder")
        if name in directories and not is_directory:
            raise RuntimeError(f"Expected an empty supporting directory: {name}")
        if {part.lower() for part in parts} & {"3drna", "3ddna", "nsp", "structure3d_runtime", "amber_runtime", "rnaalign_runtime", "rna_score_runtime", "rna_score"}:
            raise RuntimeError(f"Obsolete external engines must not be bundled: {name}")
        if ({part.lower() for part in parts} & {"oxdna", "oxpy"}
                and name not in admitted
                and not (is_directory and any(path.startswith(name + "/") for path in admitted))):
            raise RuntimeError(f"External engines outside authenticated payloads must not be bundled: {name}")
        if (any(part in {".git", ".aws", ".ssh", ".codex", ".venv", ".local_deps"}
                or part == ".env" or part.startswith(".env.") for part in parts)
                or (PurePosixPath(name).suffix.lower() in {".h5", ".hdf5", ".keras", ".pt", ".pth", ".onnx"}
                    and name not in admitted)):
            raise RuntimeError(f"Local configuration or unverified model weights must not be bundled: {name}")
        if name.startswith((f"{root}/runs/", f"{root}/designs/")) and not is_directory:
            raise RuntimeError(f"User runs must not be bundled: {name}")
    for name, target, hardlink in links:
        destination = posixpath.normpath(target if hardlink else posixpath.join(posixpath.dirname(name), target))
        if "\\" in target or not destination.startswith(root + "/"):
            raise RuntimeError(f"Archive link escapes the application folder: {name}")


class _AppDocument(HTMLParser):
    """Read the app mount point and browser entry points from shipped HTML."""

    def __init__(self) -> None:
        super().__init__(convert_charrefs=True)
        self.has_root = False
        self.assets: list[tuple[str, str]] = []

    def handle_starttag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
        attributes = dict(attrs)
        if attributes.get("id") == "root":
            self.has_root = True
        if tag == "script" and attributes.get("src"):
            self.assets.append(("javascript", attributes["src"]))
        if tag == "link" and "stylesheet" in (attributes.get("rel") or "").split() and attributes.get("href"):
            self.assets.append(("css", attributes["href"]))


def _get(opener, url: str, limit: int = 16 * 1024 * 1024) -> tuple[bytes, str]:
    try:
        response = opener.open(url, timeout=3)
    except urllib.error.HTTPError as exc:
        exc.close()
        raise
    with response:
        if response.status != 200:
            raise RuntimeError(f"HTTP {response.status}: {url}")
        body = response.read(limit + 1)
        if not body or len(body) > limit:
            raise RuntimeError(f"Empty or unexpectedly large HTTP response: {url}")
        return body, response.headers.get_content_type()


def _wait_for_server(process: subprocess.Popen, state_path: Path, opener, timeout: float = 60) -> str:
    deadline = time.monotonic() + timeout
    last_error = "The application did not write its server state."
    while time.monotonic() < deadline:
        if process.poll() is not None:
            raise RuntimeError(f"Application exited before becoming ready (exit {process.returncode}).")
        try:
            state = json.loads(state_path.read_text(encoding="utf-8"))
            port = int(state["port"])
            if not 1 <= port <= 65535:
                raise ValueError("Invalid server port")
            url = f"http://127.0.0.1:{port}/"
            body, _ = _get(opener, url + "_stcore/health", limit=1024)
            if body.strip() == b"ok":
                return url
            last_error = "The health endpoint did not return 'ok'."
        except (OSError, ValueError, KeyError, TypeError, RuntimeError) as exc:
            last_error = str(exc)
        time.sleep(0.1)
    raise RuntimeError(f"Application startup timed out: {last_error}")


def _check_http_page(opener, url: str) -> list[str]:
    html, content_type = _get(opener, url, limit=2 * 1024 * 1024)
    if content_type != "text/html":
        raise RuntimeError(f"The application returned {content_type}, not HTML: {url}")
    document = _AppDocument()
    document.feed(html.decode("utf-8"))
    if not document.has_root:
        raise RuntimeError("The application HTML is missing the Streamlit root element.")
    kinds = {kind for kind, _ in document.assets}
    if not {"javascript", "css"}.issubset(kinds):
        raise RuntimeError("The application HTML is missing its JavaScript or stylesheet entry point.")
    checked = []
    for kind, reference in dict.fromkeys(document.assets):
        asset_url = urllib.parse.urljoin(url, reference)
        parsed = urllib.parse.urlsplit(asset_url)
        origin = urllib.parse.urlsplit(url)
        if (parsed.scheme, parsed.netloc) != (origin.scheme, origin.netloc):
            raise RuntimeError(f"The standalone app depends on a remote entry point: {reference}")
        _, asset_type = _get(opener, asset_url)
        # A missing asset can fall back to index.html with HTTP 200. Browsers
        # cannot run that response as JavaScript or apply it as a stylesheet.
        expected = {"text/css"} if kind == "css" else {
            "application/javascript", "text/javascript", "application/x-javascript",
            "application/ecmascript", "text/ecmascript",
        }
        if asset_type not in expected:
            raise RuntimeError(f"Invalid {kind} response ({asset_type}): {reference}")
        checked.append(parsed.path)
    return checked


def _stop_child(process: subprocess.Popen) -> None:
    """Stop only the process created by this check, including on HTTP failures."""
    if process.poll() is None:
        try:
            process.terminate()
        except ProcessLookupError:
            pass
    try:
        process.wait(timeout=10)
    except subprocess.TimeoutExpired:
        process.kill()
        process.wait(timeout=5)


def check_http_launch(executable: Path, folder: Path, environment: dict[str, str]) -> dict:
    state_dir = folder / "server-state"
    child_environment = {
        **environment,
        "APTASWITCH_STATE_DIR": str(state_dir),
        "APTASWITCH_NO_BROWSER": "1",
    }
    log_path = folder / "server.log"
    # Ignore machine-level proxy configuration for this loopback-only check.
    opener = urllib.request.build_opener(urllib.request.ProxyHandler({}))
    try:
        with log_path.open("wb") as log:
            process = subprocess.Popen(
                [str(executable)], cwd=folder, env=child_environment,
                stdout=log, stderr=subprocess.STDOUT,
            )
            try:
                url = _wait_for_server(process, state_dir / "web_server.json", opener)
                assets = _check_http_page(opener, url)
                if process.poll() is not None:
                    raise RuntimeError(f"Application stopped while serving the page (exit {process.returncode}).")
            finally:
                _stop_child(process)
    except Exception as exc:
        server_log = log_path.read_text(encoding="utf-8", errors="replace")[-20000:] if log_path.is_file() else ""
        raise RuntimeError(f"Packaged HTTP launch failed: {exc}\nApplication log:\n{server_log}") from exc
    return {"ok": True, "checks": ["HTTP health", "HTML root", "JavaScript and CSS"], "assets": assets}


def main() -> int:
    if len(sys.argv) == 3 and sys.argv[1] == "--archive":
        check_distribution_archive(Path(sys.argv[2]))
        print("Complete offline scientific payloads and archive structure verified")
        return 0
    executable = Path(sys.argv[1]).resolve(strict=True)
    environment = os.environ.copy()
    for key in ("PYTHONHOME", "PYTHONPATH", "VIRTUAL_ENV", "APTASWITCH_HYBRID_NATIVE_LIBRARY"):
        environment.pop(key, None)
    with tempfile.TemporaryDirectory(prefix="aptaswitch-binary-check-") as folder:
        result_path = Path(folder) / "result.json"
        completed = subprocess.run(
            [str(executable), "--smoke-test", str(result_path)],
            cwd=folder, env=environment, capture_output=True, text=True, timeout=180,
        )
        if not result_path.is_file():
            raise SystemExit(f"Packaged app failed to produce a report ({completed.returncode}):\n"
                             f"{completed.stdout}\n{completed.stderr}")
        report = json.loads(result_path.read_text(encoding="utf-8"))
        print(json.dumps(report, indent=2))
        if completed.returncode or not report.get("ok") or not report.get("frozen"):
            raise SystemExit("Packaged application check failed")
        http_report = check_http_launch(executable, Path(folder), environment)
        print(json.dumps({"http_launch": http_report}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
