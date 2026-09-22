from __future__ import annotations

from contextlib import contextmanager
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
import threading
import unittest
import urllib.error
import urllib.request

from scripts.check_bundled_app import _check_http_page, _get


@contextmanager
def local_app(*, broken_html: bool = False, missing_script: bool = False, empty_css: bool = False):
    """A real HTTP regression fixture: health may be green while the app fails."""
    document = (
        b'<!DOCTYPE html><html><head><script src="./static/app.js"></script>'
        b'<link rel="stylesheet" href="./static/app.css"></head>'
        b'<body><div id="root"></div></body></html>'
    )

    class Handler(BaseHTTPRequestHandler):
        def do_GET(self):
            status, mime, body = 200, "text/html", document
            if self.path == "/_stcore/health":
                mime, body = "text/plain", b"ok"
            elif self.path == "/" and broken_html:
                status, body = 500, b"Internal Server Error"
            elif self.path == "/static/app.js" and not missing_script:
                mime, body = "application/javascript", b"window.appReady = true;"
            elif self.path == "/static/app.css":
                mime, body = "text/css", b"" if empty_css else b"body { margin: 0 }"
            self.send_response(status)
            self.send_header("Content-Type", mime)
            self.send_header("Content-Length", str(len(body)))
            self.end_headers()
            self.wfile.write(body)

        def log_message(self, *args):
            pass

    server = ThreadingHTTPServer(("127.0.0.1", 0), Handler)
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    try:
        yield f"http://127.0.0.1:{server.server_port}/"
    finally:
        server.shutdown()
        server.server_close()
        thread.join(timeout=3)


class PackagedHTTPTests(unittest.TestCase):
    def setUp(self):
        self.opener = urllib.request.build_opener(urllib.request.ProxyHandler({}))

    def test_reads_real_html_javascript_and_css(self):
        with local_app() as url:
            self.assertEqual(_check_http_page(self.opener, url), ["/static/app.js", "/static/app.css"])

    def test_healthy_server_with_broken_html_is_not_a_success(self):
        with local_app(broken_html=True) as url:
            self.assertEqual(_get(self.opener, url + "_stcore/health")[0], b"ok")
            with self.assertRaises(urllib.error.HTTPError) as raised:
                _check_http_page(self.opener, url)
            self.assertEqual(raised.exception.code, 500)

    def test_missing_javascript_disguised_by_html_fallback_is_rejected(self):
        with local_app(missing_script=True) as url:
            with self.assertRaisesRegex(RuntimeError, "Invalid javascript response"):
                _check_http_page(self.opener, url)

    def test_empty_stylesheet_is_rejected(self):
        with local_app(empty_css=True) as url:
            with self.assertRaisesRegex(RuntimeError, "Empty"):
                _check_http_page(self.opener, url)


if __name__ == "__main__":
    unittest.main()
