"""Job records and a detached worker process.

The chat page never runs analyses itself. A job is a JSON record plus a
separate ``python -m phagemine.assistant.runner`` process, so analyses survive
browser reloads. Only the runner sets COMPLETE, and only after every command
exited 0 *and* the expected artifacts validated.
"""
from __future__ import annotations

import os
import re
import signal
import subprocess
import sys
import uuid
from contextlib import contextmanager
from datetime import datetime, timezone
from pathlib import Path

try:
    import fcntl
except ImportError:  # pragma: no cover - Windows
    fcntl = None

from .store import Project, now

FINAL = {"COMPLETE", "FAILED", "CANCELLED"}
_PROGRESS = re.compile(r"\[PhageMine (\d+)/(\d+)\] (RUNNING|DONE|SKIPPED): ([^(;\n]+)")


def new_job_id(kind: str) -> str:
    prefix = re.sub(r"[^a-z]", "", kind.lower())[:12] or "job"
    return f"{prefix}-{uuid.uuid4().hex[:8]}"


def create_job(project: Project, kind: str, steps: list[list[str]], *, output: str,
               depends_on: list[str] | None = None, meta: dict | None = None,
               launch: bool = True) -> dict:
    """Record a job and start its worker. ``steps`` are argv lists run without a shell."""
    if not steps or any(not isinstance(step, list) or not step for step in steps):
        raise ValueError("A job needs at least one command")
    job_id = new_job_id(kind)
    # "{job_id}" in an output path or argument becomes this job's ID.
    output = output.replace("{job_id}", job_id)
    steps = [[part.replace("{job_id}", job_id) for part in step] for step in steps]
    for key, value in list((meta or {}).items()):
        if key == "directions":
            meta[key] = [{**d, "output": d["output"].replace("{job_id}", job_id)} for d in value]
    record = {"id": job_id, "kind": kind, "status": "QUEUED", "created": now(), "steps": steps,
              "output": output, "depends_on": list(depends_on or []), "log": f"jobs/{job_id}.log",
              "exit_codes": [], "error": None, **(meta or {})}
    project.save_job(record)
    if launch:
        start_worker(project, job_id)
    return project.job(job_id)


@contextmanager
def _locked(project: Project, job_id: str):
    """Serialise read-modify-write of one job record across processes."""
    path = project.job_path(job_id).with_suffix(".lock")
    with path.open("a") as handle:
        if fcntl:
            fcntl.flock(handle, fcntl.LOCK_EX)
        try:
            yield
        finally:
            if fcntl:
                fcntl.flock(handle, fcntl.LOCK_UN)


def transition(project: Project, job_id: str, **fields) -> dict:
    """Update a job atomically. A final status (COMPLETE, FAILED, CANCELLED) is never
    replaced; other fields (for example ``reported``) may still be recorded."""
    with _locked(project, job_id):
        record = project.job(job_id)
        if record.get("status") in FINAL and "status" in fields:
            fields = {k: v for k, v in fields.items() if k not in {"status", "error", "finished"}}
        record.update(fields)
        project.save_job(record)
        return record


def start_worker(project: Project, job_id: str) -> None:
    log = project.root / "jobs" / f"{job_id}.worker.log"
    with log.open("a") as handle:
        process = subprocess.Popen([sys.executable, "-m", "phagemine.assistant.runner", str(project.root), job_id],
                                   stdout=handle, stderr=subprocess.STDOUT, stdin=subprocess.DEVNULL,
                                   start_new_session=True, close_fds=True)
    transition(project, job_id, worker_pid=process.pid)


def _alive(pid: int | None) -> bool:
    if not pid:
        return False
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        return True
    try:  # a finished child that was never reaped is a zombie, not alive
        with open(f"/proc/{pid}/stat") as handle:
            return handle.read().split(")")[-1].split()[0] != "Z"
    except OSError:
        return True


def _age_seconds(stamp: str | None) -> float:
    try:
        return (datetime.now(timezone.utc) - datetime.fromisoformat(stamp)).total_seconds()
    except (TypeError, ValueError):
        return 0.0


def status(project: Project, job_id: str, log_lines: int = 12) -> dict:
    """Job state as recorded by the worker, plus progress parsed from its log."""
    record = project.job(job_id)
    state = record.get("status")
    if state not in FINAL and record.get("worker_pid") and not _alive(record["worker_pid"]):
        if state != "QUEUED" or _age_seconds(record.get("created")) > 30:
            record = transition(project, job_id, status="FAILED", finished=now(),
                                error="The worker process stopped unexpectedly")
    log_path = project.root / record.get("log", f"jobs/{job_id}.log")
    text = log_path.read_text(errors="replace") if log_path.is_file() else ""
    stages = _PROGRESS.findall(text)
    progress = None
    if stages:
        done, total, state, stage = stages[-1]
        progress = {"completed": int(done), "total": int(total), "state": state, "stage": stage.strip()}
    record["progress"] = progress
    record["log_tail"] = [line for line in text.splitlines() if line.strip()][-log_lines:]
    return record


def cancel(project: Project, job_id: str) -> dict:
    record = project.job(job_id)
    if record.get("status") in FINAL:
        return record
    for pid in (record.get("child_pid"), record.get("worker_pid")):
        if pid and _alive(pid):
            try:
                os.killpg(pid, signal.SIGTERM)
            except (ProcessLookupError, PermissionError):
                pass
    return transition(project, job_id, status="CANCELLED", finished=now(), error="Cancelled by the researcher")


def cli_argv(step: list[str]) -> list[str]:
    """Run ``phagemine ...`` through the current interpreter, never through a shell."""
    return [sys.executable, "-m", "phagemine", *step[1:]] if step[0] == "phagemine" else list(step)
