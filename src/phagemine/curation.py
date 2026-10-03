"""Explicit, audited annotation edits without rerunning prediction or searches."""
from __future__ import annotations

import csv
from datetime import datetime, timezone
import hashlib
import html
import json
from pathlib import Path
import shutil
from urllib.parse import quote, unquote

from .artifact_ops import checksum, completed_run, derived_run, read_json, record_sidecar

FIELDS = {"product", "gene", "note"}


def read_annotations(root: Path) -> tuple[list[str], list[dict]]:
    with (root / "annotation.tsv").open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        columns, rows = list(reader.fieldnames or []), list(reader)
    ids = [row.get("protein_id") for row in rows]
    if not ids or any(not value for value in ids) or len(set(ids)) != len(ids):
        raise ValueError("Annotation table must contain unique protein IDs")
    return columns, rows


def template(run: str | Path) -> dict:
    root = completed_run(run)
    _, rows = read_annotations(root)
    return {"schema_version": 1, "source_annotation_sha256": checksum(root / "annotation.tsv"),
            "reviewer": "", "annotations": {},
            "available_proteins": [{"protein_id": row["protein_id"], "product": row.get("product", "")}
                                   for row in rows]}


def _safe_text(value, name: str, required: bool = True) -> str:
    if not isinstance(value, str) or any(ord(char) < 32 for char in value) or (required and not value.strip()):
        raise ValueError(f"{name} must be a nonempty single-line string")
    return value.strip()


def _write_gff(root: Path, annotations: dict) -> None:
    lines = []
    for line in (root / "genes.gff3").read_text().splitlines():
        if not line or line.startswith("#"):
            lines.append(line); continue
        fields = line.split("\t")
        if len(fields) != 9:
            raise ValueError("Invalid GFF3 in completed run")
        attrs = {unquote(k): unquote(v) for k, v in (item.split("=", 1) for item in fields[8].split(";") if "=" in item)}
        row = annotations.get(attrs.get("ID"))
        if fields[2] == "CDS" and row:
            attrs.update({"product": row["product"], "gene": row.get("gene", ""),
                          "curation_status": row.get("curation_state", "AUTOMATED"),
                          "curation_note": row.get("curation_note", "")})
            fields[8] = ";".join(f"{quote(str(k), safe='')}={quote(str(v), safe='')}" for k, v in attrs.items())
        lines.append("\t".join(fields))
    (root / "genes.gff3").write_text("\n".join(lines) + "\n")


def _submission(root: Path, annotations: dict) -> None:
    """Regenerate submission outputs; never retain stale official validation."""
    from .genbank import write_package
    from .io import read_fasta
    from .models import Evidence, EvidenceLevel, Protein, SubmissionMetadata
    from .sequencing_provenance import SequencingProvenance
    old = root / "genbank_submission"
    metadata = read_json(old / "submission_metadata.json") if (old / "submission_metadata.json").is_file() else {}
    proteins = []
    for row in read_json(root / "evidence.json", list):
        fields = {key: value for key, value in row.items() if key in Protein.__dataclass_fields__}
        fields["annotation_level"] = EvidenceLevel(fields["annotation_level"])
        fields["evidence"] = [Evidence(**{**e, "level": EvidenceLevel(e["level"])}) for e in fields.get("evidence", [])]
        proteins.append(Protein(**fields))
    genome_id, sequence = read_fasta(root / "analysis_genome.fasta")
    provenance = read_json(root / "run_manifest.json")
    sequencing = read_json(root / "sequencing_provenance.json") if (root / "sequencing_provenance.json").is_file() else {}
    if old.exists():
        shutil.rmtree(old)
    overrides = {pid: {"product": row["product"], "note": row.get("curation_note") or "User-reviewed computational annotation"}
                 for pid, row in annotations.items() if row.get("curation_state") == "MANUALLY_REVIEWED"}
    write_package(root, genome_id, sequence, proteins, provenance,
                  metadata=SubmissionMetadata.from_dict(metadata),
                  sequencing_provenance=SequencingProvenance.from_dict(sequencing), feature_overrides=overrides)
    if (root / "rna_features.json").is_file():
        from .rna_features import append_submission_features
        append_submission_features(root, read_json(root / "rna_features.json")["features"])


def apply(run: str | Path, changes: dict | str | Path, output: str | Path) -> dict:
    source = completed_run(run)
    changes = read_json(changes) if not isinstance(changes, dict) else changes
    if changes.get("schema_version") != 1:
        raise ValueError("Unsupported curation schema version")
    if changes.get("source_annotation_sha256") != checksum(source / "annotation.tsv"):
        raise ValueError("Annotation table changed since review; create a fresh curation template")
    reviewer = _safe_text(changes.get("reviewer"), "reviewer")
    edits = changes.get("annotations")
    if not isinstance(edits, dict) or not edits:
        raise ValueError("Supply at least one explicit annotation edit")
    columns, rows = read_annotations(source)
    by_id = {row["protein_id"]: row for row in rows}
    if set(edits) - set(by_id):
        raise ValueError("Unknown protein ID in annotation edits")
    for pid, edit in edits.items():
        if not isinstance(edit, dict) or set(edit) - FIELDS - {"reason", "evidence_reference"}:
            raise ValueError(f"Unsupported edits for {pid}; coordinates and sequences cannot be edited")
        _safe_text(edit.get("reason"), "reason")
        if not FIELDS.intersection(edit):
            raise ValueError(f"No annotation fields supplied for {pid}")
        for field in FIELDS.intersection(edit):
            _safe_text(edit[field], field, required=field == "product")
        if "evidence_reference" in edit:
            _safe_text(edit["evidence_reference"], "evidence_reference")
        if edit.get("product", "hypothetical protein").lower() not in {"hypothetical protein", "unknown function", "uncharacterized protein"} and not edit.get("evidence_reference"):
            raise ValueError("Named products require an evidence reference in the review record")
    audit = []
    now = datetime.now(timezone.utc).isoformat()
    audit_path = source / "curation_audit.jsonl"
    previous = "0" * 64
    if audit_path.is_file():
        previous = verify_audit(audit_path)
    for pid, edit in edits.items():
        row = by_id[pid]
        for field in sorted(FIELDS.intersection(edit)):
            table_field = "curation_note" if field == "note" else field
            event = {"protein_id": pid, "field": field, "before": row.get(table_field, ""),
                     "after": edit[field], "reviewer": reviewer, "reason": edit["reason"],
                     "evidence_reference": edit.get("evidence_reference"), "timestamp": now,
                     "source_annotation_sha256": changes["source_annotation_sha256"], "previous_sha256": previous}
            previous = hashlib.sha256(json.dumps(event, sort_keys=True, separators=(",", ":")).encode()).hexdigest()
            event["sha256"] = previous; audit.append(event)
            row[table_field] = edit[field]
            if field == "product": row["proposed_function"] = edit[field]
        row.update({"curation_state": "MANUALLY_REVIEWED", "curation_reviewer": reviewer,
                    "curation_timestamp": now})
    extra = ["curation_state", "curation_note", "curation_reviewer", "curation_timestamp"]
    columns += [column for column in extra if column not in columns]
    with derived_run(source, output, "manual_annotation_curation") as root:
        if not (root / "automated_annotation.tsv").exists():
            shutil.copyfile(root / "annotation.tsv", root / "automated_annotation.tsv")
        with (root / "annotation.tsv").open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t")
            writer.writeheader(); writer.writerows(rows)
        with (root / "curation_audit.jsonl").open("a") as handle:
            for event in audit: handle.write(json.dumps(event, sort_keys=True) + "\n")
        (root / "curated_annotations.json").write_text(json.dumps(rows, indent=2, sort_keys=True) + "\n")
        _write_gff(root, by_id)
        sequences = {row["protein_id"]: row["sequence"] for row in read_json(root / "evidence.json", list)}
        (root / "annotated_proteins.faa").write_text("".join(
            f'>{row["protein_id"]} product="{row["product"].replace(chr(34), chr(39))}" curation={row.get("curation_state", "AUTOMATED")}\n{sequences[row["protein_id"]]}\n'
            for row in rows))
        _submission(root, by_id)
        record_sidecar(root, "curation", {"status": "MANUALLY_REVIEWED", "events": len(audit),
                        "audit_sha256": previous, "reviewer": reviewer, "confidence_calibrated": False})
        # Preserve original figures/reports, and make the current entry point
        # display the edited products rather than stale automated assertions.
        if not (root / "automated_report.html").exists():
            shutil.copyfile(root / "report.html", root / "automated_report.html")
        if (root / "report.md").is_file() and not (root / "automated_report.md").exists():
            shutil.copyfile(root / "report.md", root / "automated_report.md")
        def markdown_cell(value):
            return str(value).replace("|", "\\|").replace("\n", " ")
        (root / "report.md").write_text(
            "# Reviewed annotation\n\nHuman-reviewed computational annotations. Coordinates and raw evidence are unchanged.\n\n"
            "[Original analysis](automated_report.html) · [Annotation table](annotation.tsv) · [Audit trail](curation_audit.jsonl)\n\n"
            "| Protein | Product | Review |\n| --- | --- | --- |\n" + "".join(
                f"| {markdown_cell(row['protein_id'])} | {markdown_cell(row['product'])} | {markdown_cell(row.get('curation_state', 'AUTOMATED'))} |\n"
                for row in rows))
        table = "".join(f"<tr><td>{html.escape(row['protein_id'])}</td><td>{html.escape(row['product'])}</td><td>{html.escape(row.get('curation_state', 'AUTOMATED'))}</td></tr>" for row in rows)
        (root / "report.html").write_text("<!doctype html><html><meta charset='utf-8'><title>Curated PhageMine annotation</title><body><h1>Reviewed annotation</h1><p>Human-reviewed computational annotations; confidence is not empirically calibrated. Gene coordinates and raw evidence are unchanged.</p><p><a href='automated_report.html'>Original analysis</a> · <a href='annotation.tsv'>Current annotation table</a> · <a href='curation_audit.jsonl'>Audit trail</a></p><table><tr><th>Protein</th><th>Product</th><th>Review</th></tr>" + table + "</table></body></html>")
    return {"status": "CURATED", "edited_proteins": len(edits), "audit_events": len(audit), "output": str(output)}


def verify_audit(path: str | Path) -> str:
    previous = "0" * 64
    for line in Path(path).read_text().splitlines():
        event = json.loads(line); expected = event.pop("sha256")
        if event.get("previous_sha256") != previous:
            raise ValueError("Curation audit chain is broken")
        actual = hashlib.sha256(json.dumps(event, sort_keys=True, separators=(",", ":")).encode()).hexdigest()
        if actual != expected:
            raise ValueError("Curation audit checksum mismatch")
        previous = actual
    return previous
