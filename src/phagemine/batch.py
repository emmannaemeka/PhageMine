"""Sequential, checkpointed orchestration of the single-genome workflow."""
from __future__ import annotations

import csv
import hashlib
import json
import os
import re
import traceback
from datetime import datetime, timezone
from pathlib import Path

from .gene_prediction import create_predictor
from .pipeline import run
from .progress import ProgressReporter
from .resume import recover_evidence_complete
from .reporting import prefix_checkpoint_artifacts

FASTA_EXTENSIONS = {".fasta", ".fa", ".fna"}
REQUIRED_RUN = ("run_manifest.json", "evidence.json", "functional_classification.json",
                "genomic_context.json", "modules.json")
SUMMARY_FIELDS = ("sample_id", "input_file", "status", "error_message", "genome_length",
                  "predicted_proteins", "KNOWN_FUNCTION", "PROBABLE_FUNCTION",
                  "FUNCTIONAL_CLASS_ONLY", "CONSERVED_UNKNOWN", "CONFLICTING_EVIDENCE",
                  "UNRESOLVED", "module_count", "output_directory")


def _sample_id(path: Path) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", path.stem).strip("._-") or "sample"


def discover_inputs(input_dir: str | Path, recursive: bool = False) -> list[Path]:
    root = Path(input_dir)
    if not root.is_dir():
        raise ValueError(f"Batch input directory does not exist: {root}")
    paths = root.rglob("*") if recursive else root.iterdir()
    return sorted((p for p in paths if p.is_file() and p.suffix.lower() in FASTA_EXTENSIONS),
                  key=lambda p: str(p.relative_to(root)).lower())


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _valid_output(path: Path, input_path: Path | None = None) -> bool:
    if not all((path / name).is_file() for name in REQUIRED_RUN):
        return False
    if input_path is None:
        return True
    status_path = path / "sample_status.json"
    if not status_path.is_file():
        return False
    try:
        status = json.loads(status_path.read_text())
        return status.get("input_sha256") == _sha256(input_path) and status.get("status") in {"SUCCESS", "REUSED"}
    except (OSError, ValueError, TypeError):
        return False


def _atomic_json(path: Path, value) -> None:
    temp = path.with_name(path.name + ".tmp")
    temp.write_text(json.dumps(value, indent=2, sort_keys=True))
    os.replace(temp, path)


def batch(input_dir: str | Path, output: str | Path, recursive=False, resume_existing=False,
          fail_fast=False, gene_predictor="phanotate", phanotate=None, progress=None,
          reconcile_orfs=False, prodigal=None) -> list[dict]:
    root = Path(input_dir).resolve()
    project = Path(output).resolve()
    project.mkdir(parents=True, exist_ok=True)
    inputs = discover_inputs(root, recursive)
    progress = progress or ProgressReporter(quiet=True)
    started = datetime.now(timezone.utc).isoformat()
    rows: list[dict] = []
    manifest = {"batch_schema_version": "1.2", "pipeline": "PhageMine", "started_at": started,
                "ended_at": None, "input_directory": str(root), "recursive": recursive,
                "inputs": [str(p) for p in inputs], "samples": [], "failures": [],
                "configuration": {"gene_predictor": gene_predictor, "resume_existing": resume_existing,
                                   "fail_fast": fail_fast}}

    def write_summary() -> None:
        _atomic_json(project / "batch_summary.json", rows)
        temp = project / "batch_summary.tsv.tmp"
        with temp.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=SUMMARY_FIELDS, delimiter="\t")
            writer.writeheader()
            writer.writerows({key: row.get(key) for key in SUMMARY_FIELDS} for row in rows)
        os.replace(temp, project / "batch_summary.tsv")

    def persist() -> None:
        manifest["samples"] = rows
        manifest["failures"] = [row for row in rows if row["status"] == "FAILED"]
        _atomic_json(project / "batch_manifest.json", manifest)
        write_summary()

    _atomic_json(project / "batch_manifest.json", manifest)
    write_summary()
    used_ids: set[str] = set()

    for index, path in enumerate(inputs, 1):
        sid = _sample_id(path)
        base_sid = sid
        suffix = 2
        while sid in used_ids:
            sid = f"{base_sid}_{suffix}"
            suffix += 1
        used_ids.add(sid)
        destination = project / sid
        preexisting = destination.exists()
        can_reuse = resume_existing and _valid_output(destination, path)
        previous_status = None
        if preexisting and (destination / "sample_status.json").is_file():
            try:
                previous_status = json.loads((destination / "sample_status.json").read_text())
            except (OSError, ValueError):
                previous_status = None
        destination.mkdir(parents=True, exist_ok=True)
        for folder in ("checkpoints", "logs", ".working"):
            (destination / folder).mkdir(exist_ok=True)
        status_path = destination / "sample_status.json"
        input_hash = _sha256(path)
        status = {"sample_id": sid, "input_file": str(path), "input_sha256": input_hash,
                  "status": "RUNNING", "failed_stage": None, "current_stage": None,
                  "checkpoints": {}, "started_at": datetime.now(timezone.utc).isoformat()}
        if previous_status and resume_existing and previous_status.get("input_sha256") == input_hash:
            status["checkpoints"] = previous_status.get("checkpoints", {})
            if previous_status.get("status") == "FAILED" and status["checkpoints"].get("evidence_integration", {}).get("status") == "COMPLETE":
                status["audit_history"] = previous_status.get("audit_history", []) + [{"status": previous_status.get("status"), "failed_stage": previous_status.get("failed_stage"), "error_message": previous_status.get("error_message")}]
        _atomic_json(status_path, status)
        row = {"sample_id": sid, "input_file": str(path), "status": "FAILED", "error_message": None,
               "genome_length": None, "predicted_proteins": None, "module_count": None,
               "output_directory": str(destination)}
        for state in ("KNOWN_FUNCTION", "PROBABLE_FUNCTION", "FUNCTIONAL_CLASS_ONLY",
                      "CONSERVED_UNKNOWN", "CONFLICTING_EVIDENCE", "UNRESOLVED"):
            row[state] = None
        rows.append(row)
        log = destination / "logs" / "run.log"
        log.write_text(f"START {sid} {path}\n")

        def checkpoint(state_name, stage, detail):
            key = stage.replace("/", "_").replace(" ", "_")
            now = datetime.now(timezone.utc).isoformat()
            status["checkpoints"][key] = {"stage": stage, "status": state_name, "time": now, "detail": detail}
            status["current_stage"] = stage if state_name == "RUNNING" else None
            if state_name == "COMPLETE":
                prefix_checkpoint_artifacts(destination, sid)
            _atomic_json(status_path, status)
            with log.open("a") as handle:
                handle.write(f"{state_name} {stage}\n")

        sample_progress = ProgressReporter(stream=progress.stream, quiet=progress.quiet,
                                            no_progress=progress.no_progress, callback=checkpoint)
        progress._write(f"[Batch {index}/{len(inputs)}] RUNNING: {sid}")
        recovering = False
        try:
            if can_reuse:
                row["status"] = "REUSED"
                status["status"] = "REUSED"
                status["reused"] = True
                _atomic_json(status_path, status)
            elif resume_existing and previous_status and previous_status.get("input_sha256") == input_hash and status.get("audit_history") and status["checkpoints"].get("evidence_integration", {}).get("status") == "COMPLETE":
                recovering = True
                recover_evidence_complete(destination, sample_progress)
                row["status"] = "SUCCESS"
                status = json.loads(status_path.read_text())
            elif resume_existing and preexisting and status.get("checkpoints"):
                raise ValueError("checkpointed partial run cannot be resumed by the monolithic single-genome workflow")
            elif preexisting and any(destination.iterdir()):
                raise ValueError(f"output directory exists but is not a valid completed run: {destination}")
            else:
                run(path, destination, command="run", predictor=create_predictor(gene_predictor, phanotate),
                    reconcile_orfs=reconcile_orfs, prodigal=prodigal,
                    progress=sample_progress)
                row["status"] = "SUCCESS"
                status["status"] = "SUCCESS"
                _atomic_json(status_path, status)
            classification = json.loads((destination / "functional_classification.json").read_text())
            modules = json.loads((destination / "modules.json").read_text())
            fasta_lines = [line for line in (destination / "analysis_genome.fasta").read_text().splitlines() if not line.startswith(">")]
            row["genome_length"] = len("".join(fasta_lines))
            row["predicted_proteins"] = len(classification)
            row["module_count"] = len(modules)
            for state_name in ("KNOWN_FUNCTION", "PROBABLE_FUNCTION", "FUNCTIONAL_CLASS_ONLY",
                               "CONSERVED_UNKNOWN", "CONFLICTING_EVIDENCE", "UNRESOLVED"):
                row[state_name] = sum(item.get("functional_state") == state_name for item in classification)
            with log.open("a") as handle:
                handle.write("Completed single-genome workflow\n")
        except Exception as exc:
            row["error_message"] = str(exc)
            traceback_text = traceback.format_exc()
            failed_stage = status.get("current_stage")
            if recovering:
                failed_stage = "resume validation"
            if "fusion.py" in traceback_text or "classify_protein" in traceback_text:
                failed_stage = "functional classification / evidence fusion"
            status.update({"status": "FAILED", "failed_stage": failed_stage,
                           "error_type": type(exc).__name__, "error_message": str(exc),
                           "finished_at": datetime.now(timezone.utc).isoformat()})
            traceback_path = destination / "logs" / "traceback.log"
            traceback_path.write_text(traceback_text)
            status["traceback_path"] = str(traceback_path)
            _atomic_json(status_path, status)
            with log.open("a") as handle:
                handle.write(f"FAILED {status.get('failed_stage')}: {exc}\n")
            progress._write(f"[Batch {index}/{len(inputs)}] FAILED: {sid}")
            progress._write(f"Stage: {status.get('failed_stage') or 'unknown'}")
            progress._write(f"Reason: {exc}")
            progress._write(f"Resume available from: {destination}")
            persist()
            if fail_fast:
                break
            continue
        progress._write(f"[Batch {index}/{len(inputs)}] DONE: {sid}")
        persist()

    manifest["ended_at"] = datetime.now(timezone.utc).isoformat()
    persist()
    return rows
