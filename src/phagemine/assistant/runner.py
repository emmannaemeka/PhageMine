"""Detached worker: ``python -m phagemine.assistant.runner <project_root> <job_id>``.

Waits for dependencies, runs each command, then validates the expected
artifacts. COMPLETE is written only when both checks pass.
"""
from __future__ import annotations

import json
import os
import subprocess
import sys
import time
from pathlib import Path

from .jobs import FINAL, cli_argv
from .store import Project, now


def validate_outputs(project: Project, record: dict) -> None:
    """Raise ValueError unless the job's expected artifacts exist and are well formed."""
    kind, output = record["kind"], project.contained(record["output"])
    if kind == "annotation":
        from phagemine.artifact_ops import completed_run
        completed_run(output)
        manifest = json.loads((output / "run_manifest.json").read_text())
        for name in ("functional_classification.json", "genomic_context.json"):
            if not (output / name).is_file():
                raise ValueError(f"missing {name}")
        record["pipeline_version"] = manifest.get("pipeline_version")
        record["input_sha256"] = manifest.get("input_sha256")
    elif kind == "comparison":
        for direction in record.get("directions", []):
            folder = project.contained(direction["output"])
            for name in ("comparative_context.json", "ortholog_links.tsv"):
                if not (folder / name).is_file():
                    raise ValueError(f"missing {direction['query']} comparison output {name}")
    elif kind == "ncbi_revision":
        if not (output / "submission_v2" / "genbank_submission" / "validation.json").is_file():
            raise ValueError("revised submission package was not written")
        if not (output / "submission_correction_audit.json").is_file():
            raise ValueError("correction audit was not written")
    else:
        raise ValueError(f"no validator for job kind {kind}")


def _wait_for_dependencies(project: Project, record: dict) -> str | None:
    from .jobs import status
    pending = list(record.get("depends_on", []))
    while pending:
        for job_id in list(pending):
            state = status(project, job_id)["status"]
            if state == "COMPLETE":
                pending.remove(job_id)
            elif state in FINAL:
                return f"prerequisite job {job_id} ended as {state}"
        if pending:
            if project.job(record["id"]).get("status") == "CANCELLED":
                return "cancelled"
            time.sleep(2)
    return None


def run(project_root: str, job_id: str) -> int:
    project = Project(Path(project_root))
    record = project.job(job_id)
    if record.get("status") in FINAL:
        return 0
    record.update(status="WAITING" if record.get("depends_on") else "RUNNING", started=now())
    project.save_job(record)
    problem = _wait_for_dependencies(project, record)
    record = project.job(job_id)
    if record.get("status") == "CANCELLED":
        return 1
    if problem:
        record.update(status="FAILED", error=problem, finished=now())
        project.save_job(record)
        return 1
    record["status"] = "RUNNING"
    project.save_job(record)
    log_path = project.root / record["log"]
    project.contained(record["output"]).parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("a") as log:
        for step in record["steps"]:
            log.write(f"$ {' '.join(step)}\n")
            log.flush()
            process = subprocess.Popen(cli_argv(step), stdout=log, stderr=subprocess.STDOUT,
                                       stdin=subprocess.DEVNULL, cwd=project.root, shell=False,
                                       start_new_session=True, env={**os.environ, "PYTHONUNBUFFERED": "1"})
            record["child_pid"] = process.pid
            project.save_job(record)
            code = process.wait()
            record = project.job(job_id)
            record.setdefault("exit_codes", []).append(code)
            if record.get("status") == "CANCELLED":
                return 1
            if code != 0:
                record.update(status="FAILED", error=f"command exited with code {code}", finished=now())
                project.save_job(record)
                return 1
    try:
        validate_outputs(project, record)
    except (ValueError, OSError) as exc:
        record.update(status="FAILED", error=f"commands finished but outputs did not validate: {exc}", finished=now())
        project.save_job(record)
        return 1
    record.update(status="COMPLETE", error=None, finished=now())
    project.save_job(record)
    return 0


if __name__ == "__main__":
    sys.exit(run(sys.argv[1], sys.argv[2]))
