"""Console launcher for the optional Streamlit interface."""
from __future__ import annotations

import importlib.util
import socket
import sys
import threading
import time
import urllib.request
import webbrowser
from pathlib import Path


LOOPBACK = "127.0.0.1"


def available_port() -> int:
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as listener:
        listener.bind((LOOPBACK, 0))
        return int(listener.getsockname()[1])


def streamlit_arguments(app: Path, passthrough: list[str], port: int) -> list[str]:
    return ["streamlit", "run", "--global.developmentMode=false",
            "--browser.gatherUsageStats=false", "--server.headless=true",
            f"--server.address={LOOPBACK}",
            f"--browser.serverAddress={LOOPBACK}", f"--server.port={port}", *passthrough, str(app)]


def _open_when_ready(url: str, timeout: float = 30.0) -> None:
    deadline = time.monotonic() + timeout
    health = url + "/_stcore/health"
    while time.monotonic() < deadline:
        try:
            with urllib.request.urlopen(health, timeout=1) as response:
                if response.status == 200:
                    webbrowser.open(url, new=1)
                    return
        except OSError:
            time.sleep(0.2)


def main(argv: list[str] | None = None) -> int:
    if importlib.util.find_spec("streamlit") is None:
        print('PhageMine GUI requires optional dependencies. Install them with: pip install "phagemine[gui]"', file=sys.stderr)
        return 2
    from streamlit.web import cli as stcli

    app = Path(__file__).with_name("app.py")
    passthrough = list(sys.argv[1:] if argv is None else argv)
    no_browser = "--no-browser" in passthrough or "--help" in passthrough
    passthrough = [arg for arg in passthrough if arg != "--no-browser"]
    port = available_port()
    url = f"http://{LOOPBACK}:{port}"
    if not no_browser:
        threading.Thread(target=_open_when_ready, args=(url,), daemon=True).start()
    sys.argv = streamlit_arguments(app, passthrough, port)
    return int(stcli.main() or 0)
