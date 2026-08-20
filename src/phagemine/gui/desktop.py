"""Windowless desktop entry point for frozen PhageMine applications."""
from __future__ import annotations

import os
import sys
from pathlib import Path


def _configure_cli_streams() -> None:
    """Give a Windows windowed executable an explicit diagnostic stream."""
    option = "--phagemine-cli-output"
    if option in sys.argv:
        index = sys.argv.index(option)
        try:
            destination = Path(sys.argv[index + 1])
        except IndexError as exc:
            raise SystemExit(f"{option} requires a path") from exc
        del sys.argv[index:index + 2]
        destination.parent.mkdir(parents=True, exist_ok=True)
        stream = destination.open("w", encoding="utf-8")
        sys.stdout = stream
        sys.stderr = stream
    elif getattr(sys, "frozen", False) and sys.platform == "win32":
        # Windowed PyInstaller processes have no reliable CRT console stream.
        sys.stdout = open(os.devnull, "w")
        sys.stderr = open(os.devnull, "w")


def main() -> int:
    if "--phagemine-cli" in sys.argv:
        _configure_cli_streams()
        from phagemine.cli import main as cli_main
        index = sys.argv.index("--phagemine-cli")
        return cli_main(sys.argv[index + 1:])
    os.environ["PHAGEMINE_DESKTOP"] = "1"
    from phagemine.gui.launcher import main as gui_main
    return gui_main(sys.argv[1:])


if __name__ == "__main__":
    raise SystemExit(main())
