"""Launch the local Streamlit application in the user's browser."""

from __future__ import annotations

from contextvars import copy_context
import json
import os
import socket
import sys
import threading
import time
import urllib.error
import urllib.request
from pathlib import Path

from aptaswitch_core.localization import set_language, tr
from aptaswitch_core.nupack_setup import load_settings
from aptaswitch_studio import __version__
from aptaswitch_studio.desktop import open_browser


STATE_DIR = Path(os.environ.get("APTASWITCH_STATE_DIR", Path.home() / ".aptaswitch-studio")).expanduser()
SERVER_STATE_PATH = STATE_DIR / "web_server.json"


def _server_identity() -> dict[str, str]:
    installation = Path(sys.executable if getattr(sys, "frozen", False) else __file__).resolve()
    return {"app_id": "aptaswitch-studio", "version": __version__, "installation": str(installation)}


def _server_is_healthy(port: int) -> bool:
    try:
        with urllib.request.urlopen(
            f"http://127.0.0.1:{port}/_stcore/health",
            timeout=0.4,
        ) as response:
            if response.status != 200 or response.read(16).strip() != b"ok":
                return False
        # A surviving development server can still answer health checks after
        # its installation has been moved/deleted, while the browser gets 500.
        with urllib.request.urlopen(f"http://127.0.0.1:{port}/", timeout=0.8) as response:
            return (
                response.status == 200
                and response.headers.get("Content-Type", "").lower().startswith("text/html")
                and b'id="root"' in response.read(65536)
            )
    except (OSError, urllib.error.URLError):
        return False


def _running_port() -> int | None:
    try:
        data = json.loads(SERVER_STATE_PATH.read_text(encoding="utf-8"))
        port = int(data["port"])
        pid = int(data["pid"])
        if not 1 <= port <= 65535 or pid <= 0:
            return None
        if any(data.get(key) != value for key, value in _server_identity().items()):
            return None
        os.kill(pid, 0)
    except (OSError, KeyError, TypeError, ValueError, json.JSONDecodeError):
        return None
    return port if _server_is_healthy(port) else None


def _available_port() -> int:
    for port in range(8765, 8800):
        with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as probe:
            try:
                probe.bind(("127.0.0.1", port))
            except OSError:
                continue
            return port
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as probe:
        probe.bind(("127.0.0.1", 0))
        return int(probe.getsockname()[1])


def _save_server_state(port: int) -> None:
    STATE_DIR.mkdir(parents=True, exist_ok=True)
    SERVER_STATE_PATH.write_text(
        json.dumps({**_server_identity(), "pid": os.getpid(), "port": port}),
        encoding="utf-8",
    )


def _clear_own_server_state() -> None:
    try:
        data = json.loads(SERVER_STATE_PATH.read_text(encoding="utf-8"))
        if int(data.get("pid", -1)) == os.getpid():
            SERVER_STATE_PATH.unlink(missing_ok=True)
    except (OSError, TypeError, ValueError, json.JSONDecodeError):
        pass


def _open_app_browser(port: int) -> None:
    if os.environ.get("APTASWITCH_NO_BROWSER") == "1":
        return
    url = f"http://127.0.0.1:{port}"
    if not open_browser(url) and sys.stderr is not None:
        print(tr("Impossible d’ouvrir le navigateur. Ouvrez {url}.", url=url), file=sys.stderr)


def _open_browser_when_ready(port: int, stopped: threading.Event) -> None:
    deadline = time.monotonic() + 30.0
    while not stopped.is_set() and time.monotonic() < deadline:
        if _server_is_healthy(port):
            _open_app_browser(port)
            return
        stopped.wait(0.1)


def main() -> int:
    set_language(load_settings().get("ui_language", "en"))
    try:
        from streamlit.web import cli as streamlit_cli
    except ImportError as exc:
        print(
            tr("Un composant de l’application est manquant. Téléchargez à nouveau l’application complète depuis la page des versions."),
            file=sys.stderr,
        )
        raise SystemExit(1) from exc

    running_port = _running_port()
    if running_port is not None:
        _open_app_browser(running_port)
        return 0

    port = _available_port()
    _save_server_state(port)
    app_path = Path(__file__).resolve().with_name("web_app.py")
    sys.argv = [
        "streamlit",
        "run",
        str(app_path),
        "--global.developmentMode=false",
        "--server.address=127.0.0.1",
        f"--server.port={port}",
        "--server.headless=true",
        "--browser.serverAddress=127.0.0.1",
        f"--browser.serverPort={port}",
        "--browser.gatherUsageStats=false",
        "--theme.primaryColor=#1258DC",
        "--theme.backgroundColor=#F6F9FD",
        "--theme.secondaryBackgroundColor=#EAF1F8",
        "--theme.textColor=#15233A",
    ]
    stopped = threading.Event()
    browser_context = copy_context()
    browser_thread = threading.Thread(
        target=browser_context.run, args=(_open_browser_when_ready, port, stopped),
        daemon=True, name="aptaswitch-open-browser",
    )
    browser_thread.start()
    try:
        return int(streamlit_cli.main() or 0)
    finally:
        stopped.set()
        _clear_own_server_state()


if __name__ == "__main__":
    raise SystemExit(main())
