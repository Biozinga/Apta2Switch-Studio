from __future__ import annotations

import json
import os
import tempfile
import threading
import unittest
import urllib.error
from pathlib import Path
from unittest.mock import MagicMock, patch
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import aptaswitch_studio.__main__ as launcher
from aptaswitch_studio import desktop


class LauncherTests(unittest.TestCase):
    @staticmethod
    def response(body: bytes, content_type: str = "text/plain"):
        response = MagicMock()
        response.__enter__.return_value = response
        response.status = 200
        response.headers = {"Content-Type": content_type}
        response.read.return_value = body
        return response

    def test_health_endpoint_ok_does_not_accept_broken_browser_page(self):
        # Regression: the discarded source installation's process was alive
        # and returned health=ok, but opening the actual page returned 500.
        with patch.object(launcher.urllib.request, "urlopen", side_effect=[
            self.response(b"ok"),
            urllib.error.HTTPError("http://127.0.0.1:8765/", 500, "Internal Server Error", {}, None),
        ]) as request:
            self.assertFalse(launcher._server_is_healthy(8765))
        self.assertEqual(request.call_args_list[-1].args[0], "http://127.0.0.1:8765/")

    def test_health_requires_actual_streamlit_html(self):
        for body, content_type, expected in [
            (b'<html><div id="root"></div></html>', "text/html; charset=utf-8", True),
            (b"ok", "text/plain", False),
            (b"<html>Installation missing</html>", "text/html", False),
        ]:
            with self.subTest(body=body), patch.object(launcher.urllib.request, "urlopen", side_effect=[
                self.response(b"ok"), self.response(body, content_type),
            ]):
                self.assertEqual(launcher._server_is_healthy(8765), expected)

    def test_reopening_uses_safe_browser_helper_without_starting_a_server(self):
        with (
            patch.object(launcher, "load_settings", return_value={}),
            patch.object(launcher, "_running_port", return_value=8765),
            patch.object(launcher, "open_browser", return_value=True) as browser,
            patch("streamlit.web.cli.main") as cli,
        ):
            self.assertEqual(launcher.main(), 0)
        browser.assert_called_once_with("http://127.0.0.1:8765")
        cli.assert_not_called()

    def test_first_launch_disables_streamlit_browser_and_starts_health_waiter(self):
        stopped = threading.Event()
        with (
            patch.object(launcher, "load_settings", return_value={}),
            patch.object(launcher, "_running_port", return_value=None),
            patch.object(launcher, "_available_port", return_value=8791),
            patch.object(launcher, "_save_server_state") as save_state,
            patch.object(launcher, "_clear_own_server_state") as clear_state,
            patch.object(launcher.threading, "Event", return_value=stopped),
            patch.object(launcher.threading, "Thread") as thread,
            patch.object(launcher, "open_browser") as browser,
            patch.object(sys, "argv", []),
            patch("streamlit.web.cli.main", return_value=0) as cli,
        ):
            self.assertEqual(launcher.main(), 0)
            self.assertIn("--server.headless=true", sys.argv)
        save_state.assert_called_once_with(8791)
        thread.return_value.start.assert_called_once()
        self.assertEqual(thread.call_args.kwargs["args"][1:], (8791, stopped))
        self.assertTrue(stopped.is_set())
        browser.assert_not_called()
        cli.assert_called_once()
        clear_state.assert_called_once()

    def test_first_launch_opens_only_after_server_is_healthy(self):
        stopped = MagicMock(spec=threading.Event)
        stopped.is_set.return_value = False
        with (
            patch.object(launcher, "_server_is_healthy", side_effect=[False, False, True]) as health,
            patch.object(launcher, "open_browser", return_value=True) as browser,
        ):
            launcher._open_browser_when_ready(8791, stopped)
        self.assertEqual(health.call_count, 3)
        self.assertEqual(stopped.wait.call_count, 2)
        browser.assert_called_once_with("http://127.0.0.1:8791")

    def test_waiter_stops_without_opening_browser_if_server_never_starts(self):
        with (
            patch.object(launcher.time, "monotonic", side_effect=[0.0, 31.0]),
            patch.object(launcher, "_server_is_healthy") as health,
            patch.object(launcher, "open_browser") as browser,
        ):
            launcher._open_browser_when_ready(8791, threading.Event())
        health.assert_not_called()
        browser.assert_not_called()

    def test_running_port_reuses_healthy_saved_server(self):
        with tempfile.TemporaryDirectory() as tmp:
            state_path = Path(tmp) / "web_server.json"
            state_path.write_text(json.dumps({**launcher._server_identity(), "pid": os.getpid(), "port": 8765}), encoding="utf-8")
            with (
                patch.object(launcher, "SERVER_STATE_PATH", state_path),
                patch.object(launcher, "_server_is_healthy", return_value=True),
            ):
                self.assertEqual(launcher._running_port(), 8765)

    def test_old_versions_other_installations_and_legacy_servers_are_not_reused(self):
        identity = launcher._server_identity()
        for change in [
            {key: None for key in identity},
            {"version": "0.2.0b1"},
            {"installation": "/discarded/installation/app.py"},
            {"app_id": "different-app"},
        ]:
            with self.subTest(change=change), tempfile.TemporaryDirectory() as tmp:
                state_path = Path(tmp) / "web_server.json"
                state_path.write_text(json.dumps({**identity, "pid": os.getpid(), "port": 8765, **change}))
                with patch.object(launcher, "SERVER_STATE_PATH", state_path), patch.object(launcher, "_server_is_healthy") as health:
                    self.assertIsNone(launcher._running_port())
                    health.assert_not_called()

    def test_dead_process_cannot_be_reused_even_if_another_server_took_its_port(self):
        with tempfile.TemporaryDirectory() as tmp:
            state_path = Path(tmp) / "web_server.json"
            state_path.write_text(json.dumps({**launcher._server_identity(), "pid": 42, "port": 8765}))
            with (
                patch.object(launcher, "SERVER_STATE_PATH", state_path),
                patch.object(launcher.os, "kill", side_effect=ProcessLookupError),
                patch.object(launcher, "_server_is_healthy") as health,
            ):
                self.assertIsNone(launcher._running_port())
                health.assert_not_called()

    def test_new_state_records_identity_and_cleanup_preserves_another_process(self):
        with (
            tempfile.TemporaryDirectory() as tmp,
            patch.object(launcher, "STATE_DIR", Path(tmp)),
            patch.object(launcher, "SERVER_STATE_PATH", Path(tmp) / "web_server.json"),
        ):
            launcher._save_server_state(8767)
            state = json.loads(launcher.SERVER_STATE_PATH.read_text())
            self.assertEqual(state, {**launcher._server_identity(), "pid": os.getpid(), "port": 8767})
            state["pid"] = os.getpid() + 1
            launcher.SERVER_STATE_PATH.write_text(json.dumps(state))
            launcher._clear_own_server_state()
            self.assertTrue(launcher.SERVER_STATE_PATH.exists())

    def test_automated_server_check_does_not_open_a_browser(self):
        with patch.dict(os.environ, {"APTASWITCH_NO_BROWSER": "1"}), patch.object(launcher, "open_browser") as browser:
            launcher._open_app_browser(8767)
            browser.assert_not_called()

    def test_stale_or_invalid_server_state_is_ignored(self):
        with tempfile.TemporaryDirectory() as tmp:
            state_path = Path(tmp) / "web_server.json"
            state_path.write_text("not-json", encoding="utf-8")
            with patch.object(launcher, "SERVER_STATE_PATH", state_path):
                self.assertIsNone(launcher._running_port())


class DesktopTests(unittest.TestCase):
    def test_linux_frozen_restores_original_search_path_only_in_child_environment(self):
        inherited = {"LD_LIBRARY_PATH": "/bundle/_internal:/custom/lib",
                     "LD_LIBRARY_PATH_ORIG": "/custom/lib", "PATH": "/usr/bin"}
        with (
            patch.dict(os.environ, inherited, clear=True),
            patch.object(sys, "platform", "linux"),
            patch.object(sys, "frozen", True, create=True),
        ):
            child = desktop.external_process_environment()
            self.assertEqual(child["LD_LIBRARY_PATH"], "/custom/lib")
            self.assertEqual(dict(os.environ), inherited)

    def test_linux_frozen_removes_injected_path_without_original(self):
        with (
            patch.dict(os.environ, {"LD_LIBRARY_PATH": "/bundle/_internal"}, clear=True),
            patch.object(sys, "platform", "linux"),
            patch.object(sys, "frozen", True, create=True),
        ):
            self.assertNotIn("LD_LIBRARY_PATH", desktop.external_process_environment())
            self.assertEqual(os.environ["LD_LIBRARY_PATH"], "/bundle/_internal")

    def test_source_launch_preserves_ordinary_environment(self):
        with (
            patch.dict(os.environ, {"LD_LIBRARY_PATH": "/custom/lib"}, clear=True),
            patch.object(sys, "platform", "linux"),
            patch.object(sys, "frozen", False, create=True),
        ):
            self.assertEqual(desktop.external_process_environment(), dict(os.environ))

    def test_linux_browser_inherits_clean_copy_and_gio_fallback_does_too(self):
        with (
            patch.dict(os.environ, {"LD_LIBRARY_PATH": "/bundle/_internal"}, clear=True),
            patch.object(sys, "platform", "linux"),
            patch.object(sys, "frozen", True, create=True),
            patch.object(desktop.subprocess, "Popen", side_effect=[FileNotFoundError(), MagicMock()]) as process,
        ):
            self.assertTrue(desktop.open_browser("http://127.0.0.1:8791"))
            self.assertEqual(os.environ["LD_LIBRARY_PATH"], "/bundle/_internal")
        self.assertEqual(process.call_args_list[0].args[0], ["xdg-open", "http://127.0.0.1:8791"])
        self.assertEqual(process.call_args_list[1].args[0], ["gio", "open", "http://127.0.0.1:8791"])
        for call in process.call_args_list:
            self.assertNotIn("LD_LIBRARY_PATH", call.kwargs["env"])

    def test_linux_without_desktop_opener_reports_failure(self):
        with (
            patch.object(sys, "platform", "linux"),
            patch.object(desktop.subprocess, "Popen", side_effect=FileNotFoundError()),
        ):
            self.assertFalse(desktop.open_browser("http://127.0.0.1:8791"))

    def test_macos_keeps_native_browser_opening(self):
        with (
            patch.object(sys, "platform", "darwin"),
            patch.object(desktop.webbrowser, "open", return_value=True) as browser,
        ):
            self.assertTrue(desktop.open_browser("http://127.0.0.1:8791"))
        browser.assert_called_once_with("http://127.0.0.1:8791")


if __name__ == "__main__":
    unittest.main()
