"""Readers for native PhageMine annotation/discovery artifacts."""
from __future__ import annotations

import csv
import json
from pathlib import Path


def _json(path: Path, default):
    if not path.is_file():
        return default
    try:
        data = json.loads(path.read_text())
        if isinstance(default, list) and (not isinstance(data, list) or any(not isinstance(row, dict) for row in data)):
            raise ValueError("Expected an array of records")
        if isinstance(default, dict) and not isinstance(data, dict):
            raise ValueError("Expected an object")
        return data
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
    root = Path(root).expanduser()
    # BOTH copies annotation artifacts into discovery; show each sample once.
    if (root / "annotation" / "batch_manifest.json").is_file():
        root = root / "annotation"
    sample_roots = [root] if (root / "annotation.tsv").is_file() else sorted(p.parent for p in root.rglob("annotation.tsv"))
    records = []
    for sample in sample_roots:
        base = {r.get("protein_id"): r for r in _tsv(sample / "annotation.tsv")}
        classes = {r.get("protein_id"): r for r in _json(sample / "functional_classification.json", [])}
        evidence = {r.get("protein_id"): r for r in _json(sample / "evidence.json", [])}
        contexts = {r.get("protein_id"): r for r in _json(sample / "genomic_context.json", [])}
        curated = {r.get("protein_id"): r for r in _json(sample / "curated_annotations.json", [])}
        external = _json(sample / "external_evidence.json", {})
        conflicts = _json(sample / "annotation_conflicts.json", {})
        for protein_id in sorted(set(base) | set(classes)):
            row = {**base.get(protein_id, {}), **classes.get(protein_id, {})}
            row.update(curated.get(protein_id, {}))
            row["protein_id"] = protein_id
            row["sample"] = str(sample.relative_to(root)) if sample != root else sample.name
            row["record_id"] = f"{row['sample']} / {protein_id}"
            row["run_directory"] = str(sample.resolve())
            row["external_evidence"] = [record for record in external.get("records", []) if record.get("protein_id") == protein_id]
            row["conflict_review"] = [record for record in conflicts.get("records", []) if record.get("protein_id") == protein_id]
            row["evidence"] = evidence.get(protein_id, {}).get("evidence", [])
            row["sequence"] = evidence.get(protein_id, {}).get("sequence")
            row["cds"] = evidence.get(protein_id, {}).get("cds")
            row["genomic_context"] = contexts.get(protein_id)
            records.append(row)
    return records


def discovery_tables(root: str | Path) -> dict[str, list[dict]]:
    root = Path(root).expanduser()
    if not root.is_dir():
        raise ValueError(f"Result directory does not exist: {root}")
    root = discovery_root(root)
    return {name: _tsv(root / filename) for name, filename in {
        "families": "pmf_families.tsv", "members": "pmf_members.tsv",
        "recurrence": "family_recurrence.tsv", "unknown_proteome": "cohort_unknown_proteome.tsv",
        "neighbourhoods": "conserved_neighbourhoods.tsv", "validation": "pmfdb_validation.tsv",
        "ranking": "discovery_ranking.tsv", "nearest_phages": "inphared_nearest_phages.tsv"}.items()}


def discovery_root(root: str | Path) -> Path:
    """Resolve native single, discovery, and BOTH output layouts."""
    root = Path(root).expanduser()
    for candidate in (root, root / "discovery", root / "comparative"):
        if (candidate / "pmf_families.tsv").is_file():
            return candidate
    return root


def figures(root: str | Path) -> list[Path]:
    root = Path(root).expanduser()
    return sorted(p for folder in root.rglob("figures") for p in folder.iterdir() if p.is_file() and p.suffix.lower() in {".png", ".jpg", ".jpeg", ".svg"})


def downloadable_files(root: str | Path) -> list[Path]:
    root = Path(root).expanduser()
    allowed = {".tsv", ".json", ".fa", ".fasta", ".fsa", ".faa", ".fna", ".gff", ".gff3", ".png", ".jpg", ".jpeg", ".svg", ".pdf", ".md", ".html", ".gb", ".gbk", ".sqn", ".tbl", ".txt"}
    return [p for p in export_files(root) if p.suffix.lower() in allowed]


def export_files(root: str | Path) -> list[Path]:
    """Export regular result files without following links outside the run."""
    root = Path(root).expanduser().resolve()
    return sorted(p for p in root.rglob("*") if p.is_file() and not p.is_symlink()
                  and root in p.resolve().parents)
