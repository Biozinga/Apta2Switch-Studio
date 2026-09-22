"""Launch a release binary outside the checkout, without a Python environment."""

from __future__ import annotations

import json
import os
from html.parser import HTMLParser
from pathlib import Path
import subprocess
import sys
import tempfile
import time
import urllib.error
import urllib.parse
import urllib.request


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
    executable = Path(sys.argv[1]).resolve(strict=True)
    environment = os.environ.copy()
    for key in ("PYTHONHOME", "PYTHONPATH", "VIRTUAL_ENV"):
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
