"""Readers for native PhageMine annotation/discovery artifacts."""
from __future__ import annotations

import csv
import json
from pathlib import Path


def _json(path: Path, default):
    if not path.is_file():
        return default
    try:
        return json.loads(path.read_text())
    except (OSError, ValueError, TypeError) as exc:
        raise ValueError(f"Invalid PhageMine output {path}: {exc}") from exc


def _tsv(path: Path) -> list[dict]:
    if not path.is_file():
        return []
    try:
        with path.open(newline="") as handle:
            return list(csv.DictReader(handle, delimiter="\t"))
    except (OSError, csv.Error) as exc:
        raise ValueError(f"Invalid PhageMine table {path}: {exc}") from exc


def annotation_records(root: str | Path) -> list[dict]:
    root = Path(root)
    sample_roots = [root] if (root / "annotation.tsv").is_file() else sorted(p.parent for p in root.rglob("annotation.tsv"))
    records = []
    for sample in sample_roots:
        base = {r.get("protein_id"): r for r in _tsv(sample / "annotation.tsv")}
        classes = {r.get("protein_id"): r for r in _json(sample / "functional_classification.json", [])}
        evidence = {r.get("protein_id"): r for r in _json(sample / "evidence.json", [])}
        contexts = {r.get("protein_id"): r for r in _json(sample / "genomic_context.json", [])}
        for protein_id in sorted(set(base) | set(classes)):
            row = {**base.get(protein_id, {}), **classes.get(protein_id, {})}
            row["protein_id"] = protein_id; row["sample"] = sample.name
            row["evidence"] = evidence.get(protein_id, {}).get("evidence", [])
            row["sequence"] = evidence.get(protein_id, {}).get("sequence")
            row["cds"] = evidence.get(protein_id, {}).get("cds")
            row["genomic_context"] = contexts.get(protein_id)
            records.append(row)
    return records


def discovery_tables(root: str | Path) -> dict[str, list[dict]]:
    root = Path(root)
    if not root.is_dir():
        raise ValueError(f"Result directory does not exist: {root}")
    return {name: _tsv(root / filename) for name, filename in {
        "families": "pmf_families.tsv", "members": "pmf_members.tsv",
        "recurrence": "family_recurrence.tsv", "unknown_proteome": "cohort_unknown_proteome.tsv",
        "neighbourhoods": "conserved_neighbourhoods.tsv", "validation": "pmfdb_validation.tsv",
        "ranking": "discovery_ranking.tsv", "nearest_phages": "inphared_nearest_phages.tsv"}.items()}


def figures(root: str | Path) -> list[Path]:
    root = Path(root)
    return sorted(p for folder in root.rglob("figures") for p in folder.iterdir() if p.is_file() and p.suffix.lower() in {".png", ".jpg", ".jpeg", ".svg"})


def downloadable_files(root: str | Path) -> list[Path]:
    root = Path(root)
    allowed = {".tsv", ".json", ".faa", ".fna", ".gff3", ".png", ".jpg", ".jpeg", ".svg", ".pdf", ".md", ".html", ".gb", ".gbk", ".sqn", ".tbl"}
    return sorted(p for p in root.rglob("*") if p.is_file() and p.suffix.lower() in allowed)
