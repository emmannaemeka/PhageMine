"""Secure staging and validation for researcher-supplied FASTA uploads."""
from __future__ import annotations

import re
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

from phagemine.io import read_fasta

FASTA_EXTENSIONS = {".fa", ".fasta", ".fna", ".fas"}


@dataclass(frozen=True)
class UploadedGenome:
    name: str
    data: bytes


@dataclass(frozen=True)
class StagedInput:
    path: Path
    root: Path
    files: tuple[Path, ...]
    cohort: bool


def sanitize_filename(name: str) -> str:
    """Return a portable basename while preserving a supported FASTA suffix."""
    basename = Path(str(name).replace("\\", "/")).name
    stem, suffix = Path(basename).stem, Path(basename).suffix.lower()
    if suffix not in FASTA_EXTENSIONS:
        raise ValueError(f"{basename or 'File'} is not a supported FASTA file")
    safe_stem = re.sub(r"[^A-Za-z0-9._-]+", "_", stem).strip("._-") or "genome"
    return safe_stem[:120] + suffix


def validate_fasta_bytes(name: str, data: bytes) -> tuple[str, int]:
    safe_name = sanitize_filename(name)
    if not data:
        raise ValueError(f"{safe_name} is empty")
    try:
        data.decode("utf-8")
    except UnicodeDecodeError as exc:
        raise ValueError(f"{safe_name} is not UTF-8 text FASTA") from exc
    with tempfile.TemporaryDirectory(prefix="phagemine-validate-") as temporary:
        path = Path(temporary) / safe_name
        path.write_bytes(data)
        try:
            sequence_id, sequence = read_fasta(path)
        except (OSError, ValueError) as exc:
            raise ValueError(f"{safe_name}: {exc}") from exc
    return sequence_id, len(sequence)


def validate_uploads(uploads: Iterable[UploadedGenome]) -> list[dict]:
    rows = []
    for upload in uploads:
        sequence_id, length = validate_fasta_bytes(upload.name, upload.data)
        rows.append({"filename": sanitize_filename(upload.name), "sequence_id": sequence_id, "bases": length})
    return rows


def stage_uploads(uploads: Iterable[UploadedGenome]) -> StagedInput:
    items = list(uploads)
    if not items:
        raise ValueError("Add at least one FASTA file before starting analysis")
    root = Path(tempfile.mkdtemp(prefix="phagemine-upload-"))
    paths: list[Path] = []
    used: set[str] = set()
    try:
        for upload in items:
            validate_fasta_bytes(upload.name, upload.data)
            name = sanitize_filename(upload.name)
            candidate, counter = name, 2
            while candidate.casefold() in used:
                path = Path(name)
                candidate = f"{path.stem}_{counter}{path.suffix}"
                counter += 1
            used.add(candidate.casefold())
            destination = root / candidate
            destination.write_bytes(upload.data)
            paths.append(destination)
    except Exception:
        import shutil
        shutil.rmtree(root, ignore_errors=True)
        raise
    cohort = len(paths) > 1
    return StagedInput(paths[0] if not cohort else root, root, tuple(paths), cohort)


def project_name(value: str) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9._-]+", "_", value.strip()).strip("._-")
    if not cleaned:
        raise ValueError("Enter a project name")
    return cleaned[:100]


def default_results_root() -> Path:
    documents = Path.home() / "Documents"
    return (documents if documents.is_dir() else Path.home()) / "PhageMine Results"
