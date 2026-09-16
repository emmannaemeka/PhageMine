"""Assertion-level records for the prospective validation study.

The exporter is additive and never changes gene calls or functional decisions.
It provides a strict, provenance-complete record, a recursively blinded view,
review/adjudication storage, and calibration summaries that exclude unreviewed
cases.
"""
from __future__ import annotations

import csv
import hashlib
import json
import platform
import re
from dataclasses import asdict
from pathlib import Path
from typing import Any, Iterable

from .models import Protein

SCHEMA_VERSION = "validation-record-1.1"
REVIEW_LABELS = frozenset({"correct", "partially_correct", "incorrect", "unsupported", "abstain", "unresolved", "appropriate_abstention", "inappropriate_abstention"})
ASSERTION_STATES = frozenset({"supported_specific", "weak_specific", "unsupported_specific", "hypothetical_unknown", "abstained", "conflict_review"})
_IDENTITY_TERMS = re.compile(r"(?i)(phagemine|pharokka|prokka|phanotate|prodigal|hmmscan|hmmsearch|mmseqs|diamond|blast|pyhmmer|pfam|vogdb|phrogs|swiss[- ]?prot|inphared|/opt/|/usr/|/Users/|\\)")


def _json(value: Any) -> Any:
    if hasattr(value, "value"):
        return value.value
    if isinstance(value, dict):
        return {str(k): _json(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json(v) for v in value]
    return value


def _evidence_record(evidence: Any) -> dict[str, Any]:
    row = _json(asdict(evidence) if hasattr(evidence, "__dataclass_fields__") else dict(evidence))
    return {"source": row.get("source"), "evidence_type": row.get("modality"), "statement": row.get("statement"), "database_tool": row.get("source"), "database_version": row.get("source_version"), "identifier": row.get("identifier"), "metrics": row.get("metrics") or {}, "strength": row.get("evidence_strength") or row.get("level"), "supports": bool(row.get("supports", True)), "status": row.get("status"), "threshold": row.get("threshold") or {}, "coordinates": row.get("coordinates"), "provenance": row.get("provenance") or {}}


def _support_state(product: Any, evidence: list[dict[str, Any]], classification: dict[str, Any], threshold: int = 1) -> str:
    explicit = classification.get("validation_assertion_state")
    if explicit in ASSERTION_STATES:
        return explicit
    conflict = classification.get("review_flag") == "REVIEW_REQUIRED" or classification.get("functional_state") in {"CONFLICTING_EVIDENCE", "UNRESOLVED"}
    if conflict and product:
        return "conflict_review"
    if conflict:
        return "conflict_review"
    if not product:
        return "abstained" if classification.get("review_flag") == "REVIEW_REQUIRED" or classification.get("functional_state") == "UNRESOLVED" else "hypothetical_unknown"
    positive = [e for e in evidence if e.get("supports") and str(e.get("strength", "")).upper() not in {"WEAK", "WEAK INFERENCE", "NONE"}]
    weak = [e for e in evidence if e.get("supports")]
    if len(positive) >= threshold:
        return "supported_specific"
    if weak:
        return "weak_specific"
    return "unsupported_specific"


def build_validation_record(protein: Protein, classification: dict[str, Any] | None = None, manifest: dict[str, Any] | None = None, *, support_threshold: int = 1) -> dict[str, Any]:
    classification, manifest = classification or {}, manifest or {}
    genome = str(manifest.get("genome_id") or protein.genome_id)
    evidence = [_evidence_record(item) for item in protein.evidence]
    product = classification.get("display_product") or classification.get("proposed_function")
    state = _support_state(product, evidence, classification, support_threshold)
    confidence = "WEAK" if state == "weak_specific" else (classification.get("confidence") or protein.functional_confidence or "NONE")
    case_payload = f"{genome}|{protein.protein_id}|{protein.start}|{protein.end}|{protein.strand}"
    case_id = "CASE_" + hashlib.sha256(case_payload.encode()).hexdigest()[:20]
    caller = manifest.get("gene_caller") or {}
    databases = manifest.get("databases") or []
    if not databases:
        for adapter in manifest.get("evidence_adapters") or []:
            prov = adapter.get("provenance") or {}
            databases.append({"name": adapter.get("adapter"), "release": prov.get("database_version") or prov.get("version"), "checksum": prov.get("database_checksum"), "path": prov.get("path")})
    structural = manifest.get("structural_evidence") or {"exact_match": None, "tolerant_3nt": None, "tolerant_9nt": None, "tolerant_30nt": None, "independent_caller_support": [], "truth_source": None, "models": [], "matching_relationship": "predicted"}
    location_type = (protein.gene_call_parameters or {}).get("location_type") or "simple"
    return {
        "schema_version": SCHEMA_VERSION, "study_id": manifest.get("study_id"), "case_id": case_id,
        "genome": {"accession_version": genome, "length_nt": manifest.get("length_nt"), "fasta_sha256": manifest.get("input_sha256"), "input_sha256": manifest.get("input_sha256"), "panel_split": manifest.get("panel_split")},
        "tool": {"name": manifest.get("tool") or "PhageMine", "version": manifest.get("pipeline_version") or manifest.get("version"), "source_commit": manifest.get("source_commit"), "command": manifest.get("command"), "os": manifest.get("os") or platform.system(), "architecture": manifest.get("architecture") or platform.machine(), "executable": manifest.get("executable")},
        "databases": _json(databases),
        "cds": {"cds_id": protein.protein_id, "start": protein.start, "end": protein.end, "strand": protein.strand, "location_type": location_type, "protein_sha256": hashlib.sha256(protein.sequence.encode()).hexdigest(), "structural_status": "predicted", "gene_call_source": protein.gene_call_source, "gene_call_parameters": _json(protein.gene_call_parameters), "competing_models": manifest.get("competing_models") or []},
        "structural_evidence": _json(structural), "functional_evidence": evidence,
        "assertion": {"assigned_function": classification.get("normalized_function") or classification.get("proposed_function"), "assigned_product": product, "specific_product": bool(product), "hypothetical": state == "hypothetical_unknown", "evidence_sources": sorted({e.get("source") for e in evidence if e.get("source")}), "evidence_strength": classification.get("evidence_tier_label"), "confidence": confidence, "uncertainty": "REVIEW_REQUIRED" if state == "conflict_review" else ("ABSTAINED" if state == "abstained" else ("INSUFFICIENT_EVIDENCE" if state in {"hypothetical_unknown", "unsupported_specific"} else ("WEAK" if state == "weak_specific" else str(confidence).upper()))), "abstention": state in {"abstained", "hypothetical_unknown"} or (state == "conflict_review" and not product), "review_required": state == "conflict_review" or classification.get("review_flag") == "REVIEW_REQUIRED", "unsupported_risk": "unsafe" if state == "unsupported_specific" else ("possible" if state in {"weak_specific", "conflict_review"} else "none"), "system_evidence_state": state},
        "conflict": {"present": bool(classification.get("conflicting_evidence_ids") or classification.get("conflict_descriptions")), "sources": classification.get("conflicting_sources") or [], "evidence_ids": classification.get("conflicting_evidence_ids") or [], "descriptions": classification.get("conflict_descriptions") or []},
        "provenance": {"input_checksum": manifest.get("input_sha256"), "raw_output_checksum": manifest.get("raw_output_sha256"), "command": manifest.get("command"), "executable": manifest.get("executable"), "tool_version": manifest.get("pipeline_version") or manifest.get("version"), "structural_caller": {"name": caller.get("name"), "version": caller.get("version"), "parameters": caller.get("parameters") or {}}, "database_manifests": _json(databases), "evidence_adapters": manifest.get("evidence_adapters") or [], "parameters": manifest.get("parameters") or {}, "execution": manifest.get("execution") or {}, "resource_metrics": manifest.get("resource_metrics") or {"wall_seconds": None, "cpu_seconds": None, "peak_rss_bytes": None}, "platform": {"os": manifest.get("os") or platform.system(), "architecture": manifest.get("architecture") or platform.machine()}},
        "review": None, "reviews": [], "adjudication": None,
    }


def build_validation_records(proteins: Iterable[Protein], classifications: Iterable[dict[str, Any]] | None = None, manifest: dict[str, Any] | None = None, *, support_threshold: int = 1) -> list[dict[str, Any]]:
    by_id = {str(row.get("protein_id")): row for row in (classifications or [])}
    rows = [build_validation_record(p, by_id.get(p.protein_id), manifest, support_threshold=support_threshold) for p in proteins]
    for row in rows:
        validate_validation_record(row)
    return sorted(rows, key=lambda r: (r["genome"]["accession_version"] or "", r["cds"]["start"], r["cds"]["end"], r["cds"]["cds_id"]))


def validate_validation_record(record: dict[str, Any]) -> None:
    """Validate the required schema envelope; unavailable values may be null."""
    required = {"schema_version", "study_id", "case_id", "genome", "tool", "databases", "cds", "structural_evidence", "functional_evidence", "assertion", "review", "provenance"}
    missing = sorted(required - set(record))
    if missing:
        raise ValueError("validation record missing fields: " + ", ".join(missing))
    for group, fields in {"genome": ("accession_version", "length_nt", "fasta_sha256", "panel_split"), "tool": ("name", "version", "source_commit", "command", "os", "architecture"), "cds": ("cds_id", "start", "end", "strand", "location_type", "protein_sha256", "structural_status"), "assertion": ("assigned_function", "specific_product", "hypothetical", "uncertainty", "abstention", "unsupported_risk"), "provenance": ("input_checksum", "raw_output_checksum", "resource_metrics")}.items():
        absent = [field for field in fields if field not in record[group]]
        if absent:
            raise ValueError(f"{group} missing fields: {', '.join(absent)}")


def _reviewer_token(review: dict[str, Any]) -> str | None:
    return review.get("reviewer_id") or review.get("reviewer_token")


def validate_review_record(review: dict[str, Any]) -> None:
    required = {"case_id", "initial_judgement", "correctness", "evidence_adequacy", "unsupported_specificity", "abstention_appropriateness", "reviewer_confidence", "disagreement", "adjudicated_final_label", "review_time_seconds"}
    missing = sorted(required - set(review))
    if missing:
        raise ValueError("review record missing fields: " + ", ".join(missing))
    if not _reviewer_token(review):
        raise ValueError("reviewer_id must be a non-empty blinded token")
    if review["initial_judgement"] not in REVIEW_LABELS or review["adjudicated_final_label"] not in REVIEW_LABELS:
        raise ValueError("invalid review judgement label")
    if not isinstance(review["disagreement"], bool):
        raise ValueError("disagreement must be boolean")
    if review["review_time_seconds"] is not None and float(review["review_time_seconds"]) < 0:
        raise ValueError("review_time_seconds cannot be negative")


def attach_review(records: list[dict[str, Any]], review: dict[str, Any]) -> list[dict[str, Any]]:
    validate_review_record(review)
    for row in records:
        if row["case_id"] == review["case_id"]:
            updated = dict(row); reviews = list(updated.get("reviews") or []); normalized = dict(review); normalized["reviewer_id"] = _reviewer_token(review); normalized.pop("reviewer_token", None); reviews.append(normalized); updated["reviews"] = reviews; updated["review"] = reviews[0]
            return [updated if x["case_id"] == review["case_id"] else x for x in records]
    raise KeyError(f"unknown case_id: {review['case_id']}")


def adjudicate_reviews(records: list[dict[str, Any]], case_id: str, final_label: str, adjudicator_id: str, rationale: str = "") -> list[dict[str, Any]]:
    if final_label not in REVIEW_LABELS or not adjudicator_id:
        raise ValueError("invalid adjudication label or adjudicator token")
    for row in records:
        if row["case_id"] == case_id:
            updated = dict(row); updated["adjudication"] = {"adjudicated_final_label": final_label, "adjudicator_id": adjudicator_id, "rationale": rationale, "review_count": len(row.get("reviews") or []), "provenance": {"method": "independent_review_adjudication"}}
            return [updated if x["case_id"] == case_id else x for x in records]
    raise KeyError(f"unknown case_id: {case_id}")


def _blind_value(value: Any) -> Any:
    if isinstance(value, dict):
        out = {}
        for key, val in value.items():
            low = key.lower()
            # Nested evidence/provenance dictionaries may expose workflow
            # identity through arbitrary key names (for example ``hmmscan``
            # or a database acronym), not only through the known schema keys.
            if _IDENTITY_TERMS.search(str(key)):
                continue
            if low in {"tool", "tool_name", "executable", "command", "path", "source_commit", "gene_call_source", "gene_call_parameters", "execution", "platform", "database_tool", "database_version", "source_version", "source", "identifier", "evidence_adapters", "database_manifests", "databases", "structural_caller", "competing_models", "models", "caller", "provider"}:
                continue
            out[key] = _blind_value(val)
        return out
    if isinstance(value, list):
        return [_blind_value(v) for v in value]
    if isinstance(value, str):
        return _IDENTITY_TERMS.sub("[REDACTED]", value)
    return value


def blind_records(records: Iterable[dict[str, Any]], salt: str) -> list[dict[str, Any]]:
    if not salt:
        raise ValueError("salt is required for blinded export")
    output = []
    for record in records:
        row = json.loads(json.dumps(record, sort_keys=True, default=str)); case_id = row["case_id"]
        row = _blind_value(row); row["blinded_case_id"] = "BLIND_" + hashlib.sha256((salt + case_id).encode()).hexdigest()[:20]; row.pop("case_id", None)
        if "cds" in row: row["cds"]["cds_id"] = row["blinded_case_id"]
        if "review" in row: row.pop("review", None)
        row.pop("reviews", None); row.pop("adjudication", None)
        output.append(row)
    return sorted(output, key=lambda r: r["blinded_case_id"])


def write_validation_export(root: str | Path, proteins: Iterable[Protein], classifications: Iterable[dict[str, Any]] | None = None, manifest: dict[str, Any] | None = None) -> list[dict[str, Any]]:
    root = Path(root); root.mkdir(parents=True, exist_ok=True); records = build_validation_records(proteins, classifications, manifest)
    (root / "validation_records.json").write_text(json.dumps(records, indent=2, sort_keys=True, default=str) + "\n")
    columns = ["schema_version", "study_id", "case_id", "genome", "cds_id", "start", "end", "strand", "location_type", "assigned_function", "assigned_product", "system_evidence_state", "confidence", "uncertainty", "abstention", "review_required", "conflict_present", "input_checksum", "raw_output_checksum"]
    with (root / "validation_records.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t"); writer.writeheader()
        for r in records:
            writer.writerow({"schema_version": r["schema_version"], "study_id": r["study_id"] or "", "case_id": r["case_id"], "genome": r["genome"]["accession_version"], "cds_id": r["cds"]["cds_id"], "start": r["cds"]["start"], "end": r["cds"]["end"], "strand": r["cds"]["strand"], "location_type": r["cds"]["location_type"], "assigned_function": r["assertion"]["assigned_function"] or "", "assigned_product": r["assertion"]["assigned_product"] or "", "system_evidence_state": r["assertion"]["system_evidence_state"], "confidence": r["assertion"]["confidence"], "uncertainty": r["assertion"]["uncertainty"], "abstention": str(r["assertion"]["abstention"]).lower(), "review_required": str(r["assertion"]["review_required"]).lower(), "conflict_present": str(r["conflict"]["present"]).lower(), "input_checksum": r["provenance"]["input_checksum"] or "", "raw_output_checksum": r["provenance"]["raw_output_checksum"] or ""})
    return records


def calibration_bins(records: Iterable[dict[str, Any]], reviews: Iterable[dict[str, Any]], bins: tuple[str, ...] = ("HIGH", "MEDIUM", "LOW", "WEAK", "NONE")) -> list[dict[str, Any]]:
    labels = {r["case_id"]: r for r in records}; reviewed = {r["case_id"]: r for r in reviews}; out = []
    for confidence in bins:
        subset = [x for cid, x in labels.items() if str(x["assertion"].get("confidence", "NONE")).upper() == confidence and cid in reviewed]
        correct = sum(reviewed[x["case_id"]].get("correctness") in {True, "correct", "CORRECT"} or reviewed[x["case_id"]].get("adjudicated_final_label") == "correct" for x in subset)
        out.append({"confidence": confidence, "reviewed_n": len(subset), "correct_n": correct, "observed_accuracy": correct / len(subset) if subset else None, "calibration_status": "MEASURED" if subset else "NO_REVIEW_DATA"})
    return out
