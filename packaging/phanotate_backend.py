# SPDX-License-Identifier: GPL-3.0-or-later
"""Frozen command wrapper for the unmodified upstream PHANOTATE program.

This wrapper is distributed as a separate process under GPL-3.0-or-later. It
loads the PHANOTATE 1.6.7 script bundled by the desktop build and does not
change its gene-calling algorithm or defaults.
"""
from __future__ import annotations

import runpy
import sys
from pathlib import Path


def main() -> int:
    root = Path(getattr(sys, "_MEIPASS", Path(__file__).resolve().parent))
    program = root / "phanotate.py"
    if not program.is_file():
        print("Bundled PHANOTATE program is missing.", file=sys.stderr)
        return 2
    sys.argv = [str(program), *sys.argv[1:]]
    runpy.run_path(str(program), run_name="__main__")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
