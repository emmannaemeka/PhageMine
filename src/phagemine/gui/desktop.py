"""Windowless desktop entry point for frozen PhageMine applications."""
from __future__ import annotations

import os
import sys


def main() -> int:
    if "--phagemine-cli" in sys.argv:
        # PyInstaller's Windows windowed mode supplies no standard streams.
        # The CLI must still be usable internally by the GUI and smoke tests.
        if sys.stdout is None:
            sys.stdout = open(os.devnull, "w")
        if sys.stderr is None:
            sys.stderr = open(os.devnull, "w")
        from phagemine.cli import main as cli_main
        index = sys.argv.index("--phagemine-cli")
        return cli_main(sys.argv[index + 1:])
    os.environ["PHAGEMINE_DESKTOP"] = "1"
    from phagemine.gui.launcher import main as gui_main
    return gui_main(sys.argv[1:])


if __name__ == "__main__":
    raise SystemExit(main())
