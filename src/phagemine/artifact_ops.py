"""Integrity checks and atomic derived results for offline annotation workflows."""
from __future__ import annotations

from contextlib import contextmanager
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import shutil
import tempfile


def checksum(path: str | Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_json(path: str | Path, expected_type=dict):
    path = Path(path)
    try:
        value = json.loads(path.read_text())
    except (OSError, ValueError) as exc:
        raise ValueError(f"Cannot read {path.name}: {exc}") from exc
    if not isinstance(value, expected_type):
        raise ValueError(f"Invalid schema in {path.name}: expected {expected_type.__name__}")
    return value


def completed_run(root: str | Path) -> Path:
    root = Path(root).expanduser().resolve()
    for name in ("run_manifest.json", "analysis_genome.fasta", "annotation.tsv", "evidence.json", "genes.gff3"):
        if not (root / name).is_file():
            raise ValueError(f"Completed single-genome run required; missing {name}")
    manifest = read_json(root / "run_manifest.json")
    if "stage_timings_seconds" not in manifest:
        raise ValueError("Run has not completed; the final manifest is missing")
    status = root / "sample_status.json"
    if status.is_file() and read_json(status).get("status") not in {"SUCCESS", "REUSED"}:
        raise ValueError("Sample did not complete successfully")
    return root


@contextmanager
def derived_run(source: str | Path, output: str | Path, operation: str):
    """Write a new run atomically, retaining its parent and raw evidence."""
    source = completed_run(source)
    output = Path(output).expanduser().resolve()
    if output == source or source in output.parents or output in source.parents:
        raise ValueError("Output must be separate from the source result directory")
    if output.exists():
        raise ValueError("Output already exists; choose a new result directory")
    if any(path.is_symlink() for path in source.rglob("*")):
        raise ValueError("Result contains symbolic links; use a self-contained result directory")
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = Path(tempfile.mkdtemp(prefix=".phagemine-derived-", dir=output.parent))
    stage = temporary / "result"
    try:
        shutil.copytree(source, stage)
        yield stage
        manifest = read_json(stage / "run_manifest.json")
        manifest.setdefault("derived_results", []).append({
            "operation": operation, "source_manifest_sha256": checksum(source / "run_manifest.json"),
            "source_annotation_sha256": checksum(source / "annotation.tsv"),
            "recorded_at": datetime.now(timezone.utc).isoformat(),
            "gene_coordinates_unchanged": True, "raw_evidence_unchanged": True,
        })
        if checksum(stage / "evidence.json") != checksum(source / "evidence.json"):
            raise RuntimeError("Offline processing changed the original evidence")
        (stage / "run_manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
        # Cached comparisons, figures and original reports retain their original
        # meaning; invalidate archive checksums after changing any artifact.
        for name in ("reproducibility_manifest.json", "SHA256SUMS"):
            (stage / name).unlink(missing_ok=True)
        if output.exists():
            raise ValueError("Output appeared during processing; refusing to overwrite it")
        stage.rename(output)
    finally:
        shutil.rmtree(temporary, ignore_errors=True)


def record_sidecar(root: Path, name: str, payload: dict) -> None:
    manifest = read_json(root / "run_manifest.json")
    manifest[name] = payload
    (root / "run_manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
