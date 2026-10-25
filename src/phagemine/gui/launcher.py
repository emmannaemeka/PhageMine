"""Launch the bundled local GUI and open a browser after startup."""
from __future__ import annotations

import argparse
import importlib.util
import socket
import sys
import threading
import time
import urllib.error
import urllib.request
import webbrowser
from pathlib import Path


def _port(value: str) -> int:
    number = int(value)
    if not 1 <= number <= 65535:
        raise argparse.ArgumentTypeError("port must be between 1 and 65535")
    return number


def available_port(requested: int | None) -> int:
    """Keep explicit ports fixed; find a free default when another app is open."""
    candidates = [requested] if requested is not None else range(8501, 8521)
    for port in candidates:
        try:
            with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as server:
                server.bind(("127.0.0.1", port))
            return port
        except OSError:
            continue
    raise OSError(f"Port {requested} is unavailable. Choose another with --port." if requested else
                  "No free GUI port between 8501 and 8520. Choose one with --port.")


def open_when_ready(url: str, stopped: threading.Event, timeout: float = 30) -> None:
    """Do not open a blank tab while the server is still starting."""
    opener = urllib.request.build_opener(urllib.request.ProxyHandler({}))
    deadline = time.monotonic() + timeout
    while not stopped.is_set() and time.monotonic() < deadline:
        try:
            with opener.open(url + "_stcore/health", timeout=1) as response:
                ready = response.status == 200
            if ready and not stopped.is_set():
                try:
                    if not webbrowser.open(url, new=2):
                        print(f"Open PhageMine in your browser: {url}", file=sys.stderr)
                except (OSError, webbrowser.Error):
                    print(f"Open PhageMine in your browser: {url}", file=sys.stderr)
                return
        except (OSError, urllib.error.URLError):
            pass
        stopped.wait(0.2)
    if not stopped.is_set():
        print(f"Browser startup timed out. Once the server is ready, open {url}", file=sys.stderr)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Open the PhageMine desktop browser interface.")
    parser.add_argument("--port", "--server.port", type=_port, default=None, help="Local port (default: first free port from 8501)")
    parser.add_argument("--no-browser", action="store_true", help="Start the server without opening a browser")
    parser.add_argument("--server.headless", choices=("true", "false"), dest="server_headless", default=None, help=argparse.SUPPRESS)
    options, passthrough = parser.parse_known_args(sys.argv[1:] if argv is None else argv)
    # Address is deliberately local. This application can access local files and
    # launch analysis commands; deploying it publicly requires an access layer.
    if any(arg.startswith(("--server.address", "--browser.serverAddress", "--server.baseUrlPath")) for arg in passthrough):
        parser.error("The GUI binds to 127.0.0.1. Use a local browser or an SSH tunnel for remote access.")
    if importlib.util.find_spec("streamlit") is None:
        print('PhageMine GUI dependencies are missing. Repair the installation with: python -m pip install --upgrade "phagemine[gui]"', file=sys.stderr)
        return 2
    try:
        port = available_port(options.port)
    except OSError as exc:
        print(f"Could not start PhageMine GUI: {exc}", file=sys.stderr)
        return 2
    from streamlit.web import cli as stcli

    url = f"http://127.0.0.1:{port}/"
    print(f"PhageMine GUI: {url}\nKeep this terminal open. Press Ctrl+C to stop.", flush=True)
    stopped = threading.Event()
    worker = None
    if not options.no_browser and options.server_headless != "true":
        worker = threading.Thread(target=open_when_ready, args=(url, stopped), daemon=True)
        worker.start()
    original_argv = sys.argv
    app = Path(__file__).with_name("app.py")
    sys.argv = ["streamlit", "run", str(app), "--browser.gatherUsageStats=false", "--server.headless=true",
                "--server.address=127.0.0.1", f"--server.port={port}", "--client.toolbarMode=minimal", *passthrough]
    try:
        return int(stcli.main() or 0)
    finally:
        stopped.set()
        sys.argv = original_argv
