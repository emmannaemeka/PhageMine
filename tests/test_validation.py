import json

import pytest

from phagemine.models import Evidence, EvidenceLevel, Protein
from phagemine.validation import (
    adjudicate_reviews,
    build_validation_records,
    blind_records,
    calibration_bins,
    validate_review_record,
    validate_validation_record,
    write_validation_export,
)


def protein():
    p = Protein("NC_TEST.1", "PM_000001", 10, 99, "+", "ATG" * 30, "M" * 30, "PHANOTATE")
    p.evidence = [Evidence("domain", "terminase domain", EvidenceLevel.COMPUTATIONAL, "Pfam", "38.2", status="REAL", identifier="PF123", evidence_strength="STRONG", metrics={"evalue": 1e-20}, provenance={"query_protein_id": "PM_000001"})]
    return p


def test_validation_record_serializes_assertion_evidence_and_uncertainty():
    rows = build_validation_records([protein()], [{"protein_id": "PM_000001", "functional_state": "PROBABLE_FUNCTION", "display_product": "terminase large subunit", "confidence": "HIGH", "review_flag": "NONE"}], {"genome_id": "NC_TEST.1", "input_sha256": "abc", "pipeline_version": "1.1.0"})
    row = rows[0]
    assert row["genome"]["input_sha256"] == "abc"
    assert row["assertion"]["assigned_product"] == "terminase large subunit"
    assert row["functional_evidence"][0]["identifier"] == "PF123"
    assert row["assertion"]["uncertainty"] == "HIGH"


def test_conflict_and_abstention_are_explicit():
    p = protein()
    rows = build_validation_records([p], [{"protein_id": p.protein_id, "functional_state": "CONFLICTING_EVIDENCE", "confidence": "LOW", "review_flag": "REVIEW_REQUIRED", "conflicting_evidence_ids": ["E2"], "conflicting_sources": ["VOGDB"]}], {})
    assert rows[0]["assertion"]["uncertainty"] == "REVIEW_REQUIRED"
    assert rows[0]["assertion"]["abstention"] is True
    assert rows[0]["conflict"]["evidence_ids"] == ["E2"]


def test_export_is_deterministic_and_review_validation_blinds_identity(tmp_path):
    p = protein(); cls = [{"protein_id": p.protein_id, "functional_state": "UNRESOLVED", "confidence": "NONE"}]
    first = write_validation_export(tmp_path / "a", [p], cls, {"genome_id": p.genome_id, "input_sha256": "abc"})
    second = write_validation_export(tmp_path / "b", [p], cls, {"genome_id": p.genome_id, "input_sha256": "abc"})
    assert json.loads((tmp_path / "a" / "validation_records.json").read_text()) == json.loads((tmp_path / "b" / "validation_records.json").read_text())
    blinded = blind_records(first, "study-salt")
    assert "tool" not in blinded[0] and blinded[0]["blinded_case_id"].startswith("BLIND_")


def test_review_schema_and_calibration_require_actual_labels():
    with pytest.raises(ValueError):
        validate_review_record({"case_id": "x"})
    review = {"case_id": "x", "reviewer_token": "R1", "initial_judgement": "abstain", "correctness": "correct", "evidence_adequacy": "adequate", "unsupported_specificity": "none", "abstention_appropriateness": "appropriate", "reviewer_confidence": "HIGH", "disagreement": False, "adjudicated_final_label": "appropriate_abstention", "review_time_seconds": 2}
    validate_review_record(review)
    assert calibration_bins([], [])[0]["calibration_status"] == "NO_REVIEW_DATA"


def test_schema_required_fields_are_present_or_explicit_null():
    row = build_validation_records([protein()], manifest={"study_id": "S1", "genome_id": "NC_TEST.1", "input_sha256": "sha", "command": None})[0]
    assert {"study_id", "genome", "tool", "databases", "cds", "structural_evidence", "functional_evidence", "assertion", "review", "provenance"} <= set(row)
    assert row["study_id"] == "S1" and row["tool"]["command"] is None and row["provenance"]["raw_output_checksum"] is None
    assert {"length_nt", "fasta_sha256", "panel_split"} <= set(row["genome"])
    validate_validation_record(row)


def test_blinding_removes_nested_identity_and_preserves_biology():
    p = protein(); p.gene_call_parameters = {"tool": "PHANOTATE", "path": "/Users/private"}
    rows = build_validation_records([p], manifest={"tool": "PhageMine", "genome_id": p.genome_id, "evidence_adapters": [{"adapter": "PhageMine", "provenance": {"path": "/opt/db", "database_version": "1"}}]})
    text = json.dumps(blind_records(rows, "salt"), sort_keys=True)
    for term in ("PhageMine", "PHANOTATE", "Pharokka", "Prokka", "hmmscan", "/opt/", "/Users/"):
        assert term.lower() not in text.lower()
    blinded = blind_records(rows, "salt")[0]
    assert blinded["genome"]["accession_version"] == p.genome_id
    assert blinded["cds"]["start"] == p.start


def test_support_states_and_multi_reviewer_adjudication():
    p = protein(); base = {"genome_id": p.genome_id}
    weak = build_validation_records([p], [{"protein_id": p.protein_id, "display_product": "weak product", "validation_assertion_state": "weak_specific"}], base)[0]
    unsupported = build_validation_records([p], [{"protein_id": p.protein_id, "display_product": "unsupported product", "validation_assertion_state": "unsupported_specific"}], base)[0]
    assert weak["assertion"]["system_evidence_state"] == "weak_specific"
    assert unsupported["assertion"]["system_evidence_state"] == "unsupported_specific"
    reviews = []
    for token in ("R1", "R2"):
        reviews.append({"case_id": weak["case_id"], "reviewer_id": token, "initial_judgement": "correct", "correctness": "correct", "evidence_adequacy": "adequate", "unsupported_specificity": "none", "abstention_appropriateness": "not_applicable", "reviewer_confidence": "HIGH", "disagreement": token == "R2", "adjudicated_final_label": "correct", "review_time_seconds": 4})
    reviewed = [weak]
    for review in reviews: reviewed = __import__("phagemine.validation", fromlist=["attach_review"]).attach_review(reviewed, review)
    reviewed = adjudicate_reviews(reviewed, weak["case_id"], "correct", "ADJ", "consensus")
    assert len(reviewed[0]["reviews"]) == 2 and reviewed[0]["adjudication"]["review_count"] == 2
    assert next(x for x in calibration_bins(reviewed, reviews) if x["confidence"] == "WEAK")["reviewed_n"] == 1
