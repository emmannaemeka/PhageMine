"""Per-project storage: read-only inputs, job records, conversation, outputs.

Layout::

    <projects root>/<project_id>/
        project.json        name, created time
        inputs/             copies of uploaded files (never modified)
        jobs/<job_id>.json  job records written by the worker
        runs/<job_id>/      PhageMine outputs, one new folder per job
        figures/  summaries/  ncbi/
        conversation.json   model-format messages
        transcript.json     what the chat page shows
        outbound.jsonl      audit log of payloads sent to the model
"""
from __future__ import annotations

import hashlib
import json
import os
import re
import stat
import uuid
from datetime import datetime, timezone
from pathlib import Path

FASTA_SUFFIXES = {".fa", ".fasta", ".fna", ".fas", ".fsa"}
_SAFE_NAME = re.compile(r"[^A-Za-z0-9._-]+")


def now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def safe_name(name: str) -> str:
    cleaned = _SAFE_NAME.sub("_", Path(name).name).strip("._")
    if not cleaned:
        raise ValueError("File name is empty after removing unsafe characters")
    return cleaned


def projects_root() -> Path:
    root = os.environ.get("PHAGEMINE_PROJECTS") or str(Path.home() / "PhageMine" / "projects")
    path = Path(root).expanduser().resolve()
    path.mkdir(parents=True, exist_ok=True)
    return path


def _read_json(path: Path, default):
    try:
        return json.loads(path.read_text())
    except (OSError, ValueError):
        return default


def _write_json(path: Path, payload) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{uuid.uuid4().hex}.tmp")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True, default=str) + "\n")
    os.replace(temporary, path)


class Project:
    """One isolated analysis project. Every path a tool touches is checked against it."""

    def __init__(self, root: Path):
        self.root = Path(root).resolve()
        if not (self.root / "project.json").is_file():
            raise ValueError(f"Not a PhageMine project: {self.root}")
        self.meta = _read_json(self.root / "project.json", {})

    # ---- creation and listing -------------------------------------------------
    @classmethod
    def create(cls, name: str, base: Path | None = None) -> "Project":
        base = Path(base) if base else projects_root()
        slug = safe_name(name.lower().replace(" ", "-"))[:40] or "project"
        project_id = f"{slug}-{uuid.uuid4().hex[:6]}"
        root = base / project_id
        for folder in ("inputs", "jobs", "runs", "figures", "summaries", "ncbi"):
            (root / folder).mkdir(parents=True, exist_ok=True)
        _write_json(root / "project.json", {"id": project_id, "name": name, "created": now()})
        return cls(root)

    @staticmethod
    def list(base: Path | None = None) -> list["Project"]:
        base = Path(base) if base else projects_root()
        found = [Project(path) for path in sorted(base.iterdir()) if (path / "project.json").is_file()]
        return sorted(found, key=lambda p: p.meta.get("created", ""), reverse=True)

    @property
    def id(self) -> str:
        return self.meta.get("id", self.root.name)

    @property
    def name(self) -> str:
        return self.meta.get("name", self.id)

    # ---- path safety ----------------------------------------------------------
    def contained(self, path: str | Path) -> Path:
        """Resolve a path and refuse anything outside this project."""
        resolved = (self.root / path).resolve() if not Path(path).is_absolute() else Path(path).resolve()
        if resolved != self.root and self.root not in resolved.parents:
            raise ValueError("Path is outside this project")
        return resolved

    # ---- inputs ---------------------------------------------------------------
    def add_input(self, name: str, data: bytes) -> dict:
        """Copy an uploaded file into inputs/ and make it read-only."""
        target_name = safe_name(name)
        target = self.root / "inputs" / target_name
        if target.exists():
            if hashlib.sha256(data).hexdigest() == sha256(target):
                return self.input_record(target)
            raise ValueError(f"An input called {target_name} already exists with different content; rename the file")
        target.write_bytes(data)
        target.chmod(stat.S_IRUSR | stat.S_IRGRP | stat.S_IROTH)
        return self.input_record(target)

    def input_record(self, path: Path) -> dict:
        return {"file": path.name, "bytes": path.stat().st_size, "sha256": sha256(path),
                "kind": "fasta" if path.suffix.lower() in FASTA_SUFFIXES else "other"}

    def inputs(self) -> list[dict]:
        return [self.input_record(path) for path in sorted((self.root / "inputs").iterdir()) if path.is_file()]

    def input_path(self, name: str) -> Path:
        path = self.contained(Path("inputs") / safe_name(name))
        if not path.is_file():
            raise ValueError(f"No uploaded file called {name}")
        return path

    # ---- jobs -----------------------------------------------------------------
    def job_path(self, job_id: str) -> Path:
        if not re.fullmatch(r"[a-z]+-[0-9a-f]{8}", job_id or ""):
            raise ValueError(f"Unknown job ID: {job_id}")
        return self.root / "jobs" / f"{job_id}.json"

    def job(self, job_id: str) -> dict:
        path = self.job_path(job_id)
        if not path.is_file():
            raise ValueError(f"Unknown job ID: {job_id}")
        return _read_json(path, {})

    def jobs(self) -> list[dict]:
        records = [_read_json(path, {}) for path in (self.root / "jobs").glob("*.json")]
        return sorted((r for r in records if r), key=lambda r: r.get("created", ""))

    def save_job(self, record: dict) -> None:
        _write_json(self.job_path(record["id"]), record)

    def annotation_runs(self) -> dict[str, dict]:
        """Completed annotation jobs keyed by sample name (latest wins)."""
        runs = {}
        for job in self.jobs():
            if job.get("kind") == "annotation" and job.get("status") == "COMPLETE":
                runs[job["sample"]] = job
        return runs

    def run_dir(self, sample: str) -> Path:
        runs = self.annotation_runs()
        if sample not in runs:
            known = ", ".join(sorted(runs)) or "none yet"
            raise ValueError(f"No completed annotation for '{sample}'. Completed samples: {known}")
        return self.contained(runs[sample]["output"])

    # ---- conversation ---------------------------------------------------------
    def messages(self) -> list[dict]:
        return _read_json(self.root / "conversation.json", [])

    def save_messages(self, messages: list[dict]) -> None:
        _write_json(self.root / "conversation.json", messages)

    def transcript(self) -> list[dict]:
        return _read_json(self.root / "transcript.json", [])

    def append_transcript(self, entry: dict) -> None:
        items = self.transcript()
        items.append({**entry, "time": now()})
        _write_json(self.root / "transcript.json", items)

    def user_turns(self) -> int:
        return sum(1 for item in self.transcript() if item.get("role") == "user")

    def log_outbound(self, payload: dict) -> None:
        with (self.root / "outbound.jsonl").open("a") as handle:
            handle.write(json.dumps({"time": now(), **payload}, default=str) + "\n")

    # ---- generic small JSON state --------------------------------------------
    def read_state(self, relative: str, default):
        return _read_json(self.contained(relative), default)

    def write_state(self, relative: str, payload) -> None:
        _write_json(self.contained(relative), payload)

    def outputs(self) -> list[Path]:
        """Files the researcher can download: figures, summaries, NCBI previews."""
        found = []
        for folder in ("figures", "summaries"):
            found += sorted(p for p in (self.root / folder).glob("*") if p.is_file())
        return found
