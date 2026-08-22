"""Safe subprocess boundary for invoking the unchanged PhageMine CLI."""
from __future__ import annotations

import json
import shlex
import subprocess
import sys
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Sequence

from phagemine.gui import GUI_VERSION

MODES = ("annotate", "discover", "both")
EVIDENCE_LEVELS = ("core", "standard", "full")


def build_cli_args(input_path: str | Path, output: str | Path, *, mode: str = "annotate",
                   evidence: str = "core", threads: int = 1, cohort: bool | None = None,
                   resume: bool = False) -> list[str]:
    source, destination = Path(input_path), Path(output)
    if mode not in MODES:
        raise ValueError(f"unsupported mode: {mode}")
    if evidence not in EVIDENCE_LEVELS:
        raise ValueError(f"unsupported evidence level: {evidence}")
    if int(threads) < 1:
        raise ValueError("threads must be at least 1")
    is_cohort = source.is_dir() if cohort is None else cohort
    if not is_cohort:
        # Single-genome CLI has no evidence-profile switch. Core maps exactly to
        # its default adapters; richer profiles are available through batch.
        if evidence != "core":
            raise ValueError("standard/full evidence profiles require cohort (batch) input")
        command = {"annotate": "annotate", "discover": "mine", "both": "run"}[mode]
        return ["phagemine", command, str(source), "--output", str(destination), "--threads", str(int(threads))]
    args = ["phagemine", "batch", str(source), "--output", str(destination),
            "--mode", mode, "--evidence", evidence, "--threads", str(int(threads))]
    if resume:
        args.append("--resume-existing")
    return args


def display_command(args: Sequence[str]) -> str:
    return shlex.join(list(args))


@dataclass
class RunResult:
    command: list[str]
    stdout: str
    stderr: str
    exit_code: int
    output_directory: Path


def execute(args: Sequence[str]) -> RunResult:
    """Run via the installed CLI semantics, without a shell."""
    command = list(args)
    executable = [sys.executable, "-m", "phagemine", *command[1:]] if command and command[0] == "phagemine" else command
    completed = subprocess.run(executable, capture_output=True, text=True, check=False, shell=False)
    try:
        output = Path(command[command.index("--output") + 1])
    except (ValueError, IndexError):
        output = Path(".")
    result = RunResult(command, completed.stdout, completed.stderr, completed.returncode, output)
    if completed.returncode == 0:
        write_gui_metadata(output, result)
    return result


def write_gui_metadata(output: Path, result: RunResult) -> None:
    output.mkdir(parents=True, exist_ok=True)
    payload = {"interface": "PhageMine GUI", "gui_version": GUI_VERSION,
               "executed_command": result.command, "exit_code": result.exit_code,
               "recorded_at": datetime.now(timezone.utc).isoformat()}
    (output / "gui_metadata.json").write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
