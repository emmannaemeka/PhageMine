"""Deterministic, non-calibrated review-priority summaries."""
from __future__ import annotations

NO_REVIEW = "NO_REVIEW"
LOW_PRIORITY_REVIEW = "LOW_PRIORITY_REVIEW"
MODERATE_REVIEW = "MODERATE_REVIEW"
HIGH_PRIORITY_REVIEW = "HIGH_PRIORITY_REVIEW"


def review_priority(locus, decision=None) -> str:
    """Map existing reconciliation/adjudication states to review tiers.

    This is a presentation tier only. It introduces no score or biological
    threshold; the original conflict and ``review_required`` fields remain
    authoritative audit data.
    """
    cls = locus.reconciliation_class
    decision_class = getattr(decision, "decision_class", None)
    if locus.strand_status == "DISCORDANT" or cls == "STRAND_DISCORDANCE":
        return HIGH_PRIORITY_REVIEW
    if cls in {"SPLIT_MODEL", "MERGED_MODEL", "COMPLEX_CONFLICT"} or decision_class in {
        "STRAND_CONFLICT_UNRESOLVED", "SPLIT_MERGE_UNRESOLVED", "COMPLEX_CONFLICT",
        "CONFLICTING_EVIDENCE",
    }:
        return HIGH_PRIORITY_REVIEW
    if cls == "COMMON_STOP_ALTERNATE_START":
        return LOW_PRIORITY_REVIEW
    if cls in {"COMMON_START_ALTERNATE_STOP", "BOUNDARY_DISCORDANCE", "CALLER_SPECIFIC"}:
        return MODERATE_REVIEW
    if getattr(locus, "review_required", False):
        return MODERATE_REVIEW
    return NO_REVIEW
