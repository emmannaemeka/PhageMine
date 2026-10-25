"""Detached worker: ``python -m phagemine.assistant.runner <project_root> <job_id>``.

Waits for dependencies, runs each command, then validates the outputs.
COMPLETE is written only when both pass. All state changes go through
``jobs.transition``, which never replaces a final status.
"""
from __future__ import annotations

import csv
import json
import os
import subprocess
import sys
import time
from pathlib import Path

from .jobs import FINAL, cli_argv, transition
from .store import Project, now


def _json(path: Path):
    try:
        return json.loads(path.read_text())
    except (OSError, ValueError) as exc:
        raise ValueError(f"{path.name} is missing or not valid JSON: {exc}") from exc


def _run_proteins(run: Path) -> set[str]:
    return {row["protein_id"] for row in _json(run / "evidence.json")}


def validate_comparison_direction(project: Project, record: dict, direction: dict) -> dict:
    """Parse one direction's outputs and check every reference against its source runs."""
    folder = project.contained(direction["output"])
    rows = _json(folder / "comparative_context.json")
    manifest = _json(folder / "comparison_manifest.json")
    if not isinstance(rows, list) or any(not isinstance(r, dict) or "query_protein_id" not in r for r in rows):
        raise ValueError(f"{direction['query']}: comparative_context.json has an unexpected structure")
    runs = {name: project.contained(path) for name, path in record["run_dirs"].items()}
    query_run = runs[direction["query"]]
    query_ids = _run_proteins(query_run)
    if {r["query_protein_id"] for r in rows} != query_ids:
        raise ValueError(f"{direction['query']}: comparison rows do not match the annotated proteins")
    by_path = {str(path): name for name, path in runs.items()}
    references = {name: _run_proteins(path) for name, path in runs.items() if name != direction["query"]}
    with (folder / "ortholog_links.tsv").open(newline="") as handle:
        links = list(csv.DictReader(handle, delimiter="\t"))
    for link in links:
        target = by_path.get(str(Path(link["reference_genome_id"]).resolve()))
        if target is None or target == direction["query"]:
            raise ValueError(f"{direction['query']}: link points to an unknown run {link['reference_genome_id']}")
        if link["query_protein_id"] not in query_ids or link["reference_protein_id"] not in references[target]:
            raise ValueError(f"{direction['query']}: link references a protein that is not in its run")
    hashes = {name: _json(path / "run_manifest.json").get("input_sha256") for name, path in runs.items()}
    recorded = [manifest["query_run"], *manifest["reference_runs"]]
    for entry in recorded:
        name = by_path.get(str(Path(entry["path"]).resolve()))
        if name is None or entry.get("input_sha256") != hashes[name]:
            raise ValueError(f"{direction['query']}: comparison manifest does not match the source runs")
    return manifest


def validate_outputs(project: Project, record: dict) -> dict:
    """Raise ValueError unless the job's outputs are complete and consistent.
    Returns facts to store on the job record."""
    kind, output = record["kind"], project.contained(record["output"])
    if kind == "annotation":
        from phagemine.operational_validation import validate_run
        validate_run(output)  # protein IDs, sequences and coordinates agree across all outputs
        for name in ("functional_classification.json", "genomic_context.json"):
            _json(output / name)
        manifest = _json(output / "run_manifest.json")
        return {"pipeline_version": manifest.get("pipeline_version"), "input_sha256": manifest.get("input_sha256")}
    if kind == "comparison":
        mmseqs = {"searches_attempted": 0, "searches_succeeded": 0, "searches_failed": 0}
        for direction in record.get("directions", []):
            manifest = validate_comparison_direction(project, record, direction)
            for key in mmseqs:
                mmseqs[key] += int(manifest["mmseqs"].get(key) or 0)
        executed = mmseqs["searches_succeeded"] > 0 and mmseqs["searches_failed"] == 0
        return {"mmseqs2_executed": executed, "mmseqs_searches": mmseqs}
    if kind == "ncbi_revision":
        validation = _json(output / "submission_v2" / "genbank_submission" / "validation.json")
        _json(output / "submission_correction_audit.json")
        # Job completion is not submission readiness: record both separately.
        return {"submission_readiness": validation.get("submission_readiness")}
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
    record = transition(project, job_id, status="WAITING" if record.get("depends_on") else "RUNNING", started=now())
    problem = _wait_for_dependencies(project, record)
    if problem:
        transition(project, job_id, status="FAILED", error=problem, finished=now())
        return 1
    record = transition(project, job_id, status="RUNNING")
    if record["status"] in FINAL:
        return 1
    log_path = project.root / record["log"]
    project.contained(record["output"]).parent.mkdir(parents=True, exist_ok=True)
    codes = []
    with log_path.open("a") as log:
        for step in record["steps"]:
            log.write(f"$ {' '.join(step)}\n")
            log.flush()
            process = subprocess.Popen(cli_argv(step), stdout=log, stderr=subprocess.STDOUT,
                                       stdin=subprocess.DEVNULL, cwd=project.root, shell=False,
                                       start_new_session=True, env={**os.environ, "PYTHONUNBUFFERED": "1"})
            transition(project, job_id, child_pid=process.pid)
            code = process.wait()
            codes.append(code)
            record = transition(project, job_id, exit_codes=list(codes))
            if record["status"] in FINAL:  # cancelled while running
                return 1
            if code != 0:
                transition(project, job_id, status="FAILED", error=f"command exited with code {code}", finished=now())
                return 1
    try:
        facts = validate_outputs(project, project.job(job_id))
    except (ValueError, OSError, KeyError) as exc:
        transition(project, job_id, status="FAILED", error=f"commands finished but outputs did not validate: {exc}", finished=now())
        return 1
    final = transition(project, job_id, status="COMPLETE", error=None, finished=now(), **facts)
    return 0 if final["status"] == "COMPLETE" else 1


if __name__ == "__main__":
    sys.exit(run(sys.argv[1], sys.argv[2]))
