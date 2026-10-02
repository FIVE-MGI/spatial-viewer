#!/usr/bin/env python
r"""Serve an export folder over HTTP with CORS, for testing before upload.

    python static_serve.py <export root> 8765

then open the page with ?static=http://localhost:8765 - that path reads the files
directly and skips the password and the blob token entirely.

The README used to point at a copy of this that only existed in a scratch folder;
it lives here now so the local workflow keeps working.
"""
import sys
from functools import partial
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer


class CORSHandler(SimpleHTTPRequestHandler):
    def end_headers(self):
        self.send_header("Access-Control-Allow-Origin", "*")
        self.send_header("Cache-Control", "no-store")   # re-exports show up without a hard reload
        super().end_headers()

    def log_message(self, fmt, *args):
        if "404" in (fmt % args):
            super().log_message(fmt, *args)


if __name__ == "__main__":
    root = sys.argv[1] if len(sys.argv) > 1 else "."
    port = int(sys.argv[2]) if len(sys.argv) > 2 else 8765
    handler = partial(CORSHandler, directory=root)
    print(f"serving {root} on http://localhost:{port}  (Ctrl+C to stop)")
    ThreadingHTTPServer(("127.0.0.1", port), handler).serve_forever()
