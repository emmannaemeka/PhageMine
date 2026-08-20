"""Background process management without duplicating PhageMine checkpoints."""
from __future__ import annotations

import os
import shutil
import subprocess
import tempfile
import uuid
from dataclasses import dataclass
from pathlib import Path

from .execution import RunResult, execution_argv, output_from_command, write_gui_metadata


@dataclass
class RunHandle:
    run_id: str
    command: list[str]
    process: subprocess.Popen
    stdout_path: Path
    stderr_path: Path
    stdout_handle: object
    stderr_handle: object
    output_directory: Path
    input_cleanup: Path | None = None
    finalized: bool = False


RUNS: dict[str, RunHandle] = {}


def start_run(args: list[str], *, input_cleanup: Path | None = None) -> RunHandle:
    command = list(args)
    log_root = Path(tempfile.mkdtemp(prefix="phagemine-run-"))
    stdout_path, stderr_path = log_root / "stdout.log", log_root / "stderr.log"
    stdout_handle, stderr_handle = stdout_path.open("w"), stderr_path.open("w")
    kwargs = {}
    if os.name == "nt":
        kwargs["creationflags"] = getattr(subprocess, "CREATE_NO_WINDOW", 0)
    process = subprocess.Popen(execution_argv(command), stdout=stdout_handle, stderr=stderr_handle,
                               text=True, shell=False, **kwargs)
    handle = RunHandle(uuid.uuid4().hex, command, process, stdout_path, stderr_path,
                       stdout_handle, stderr_handle, output_from_command(command), input_cleanup)
    RUNS[handle.run_id] = handle
    return handle


def poll_run(run_id: str) -> dict:
    handle = RUNS.get(run_id)
    if handle is None:
        return {"state": "UNKNOWN", "exit_code": None, "stdout": "", "stderr": ""}
    exit_code = handle.process.poll()
    if exit_code is not None and not handle.finalized:
        handle.stdout_handle.close(); handle.stderr_handle.close()
        if exit_code == 0:
            write_gui_metadata(handle.output_directory, RunResult(handle.command, "", "", exit_code, handle.output_directory))
        if handle.input_cleanup:
            shutil.rmtree(handle.input_cleanup, ignore_errors=True)
        handle.finalized = True
    stdout = handle.stdout_path.read_text(errors="replace") if handle.stdout_path.is_file() else ""
    stderr = handle.stderr_path.read_text(errors="replace") if handle.stderr_path.is_file() else ""
    return {"state": "RUNNING" if exit_code is None else ("COMPLETE" if exit_code == 0 else "FAILED"),
            "exit_code": exit_code, "stdout": stdout, "stderr": stderr,
            "command": handle.command, "output_directory": str(handle.output_directory)}
