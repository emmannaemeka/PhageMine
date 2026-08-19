"""Console launcher for the optional Streamlit interface."""
from __future__ import annotations

import importlib.util
import sys
from pathlib import Path


def main(argv: list[str] | None = None) -> int:
    if importlib.util.find_spec("streamlit") is None:
        print('PhageMine GUI requires optional dependencies. Install them with: pip install "phagemine[gui]"', file=sys.stderr)
        return 2
    from streamlit.web import cli as stcli

    app = Path(__file__).with_name("app.py")
    passthrough = sys.argv[1:] if argv is None else argv
    sys.argv = ["streamlit", "run", "--browser.gatherUsageStats=false", "--server.headless=true", *passthrough, str(app)]
    return int(stcli.main() or 0)
