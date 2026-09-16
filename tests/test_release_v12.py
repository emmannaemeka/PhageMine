from types import SimpleNamespace

from phagemine.io import normalize_sequence
from phagemine.review_priority import (
    HIGH_PRIORITY_REVIEW, LOW_PRIORITY_REVIEW, MODERATE_REVIEW, NO_REVIEW,
    review_priority,
)


def locus(cls, strand_status="CONCORDANT", review_required=True):
    return SimpleNamespace(reconciliation_class=cls, strand_status=strand_status,
                           review_required=review_required)


def test_sequence_normalization_is_auditable():
    normalized, audit = normalize_sequence(" acgRyu? ", record_id="g", provider_id="pyrodigal")
    assert normalized == "ACGNNNN"
    assert audit["replacement_count"] == 4
    assert audit["replaced_symbol_counts"] == {"?": 1, "R": 1, "U": 1, "Y": 1}


def test_sequence_normalization_rejects_empty_or_non_nucleotide():
    import pytest
    with pytest.raises(ValueError, match="empty"):
        normalize_sequence(" \n\t", record_id="empty")
    with pytest.raises(ValueError, match="no usable"):
        normalize_sequence("!?", record_id="bad")


def test_review_priority_mapping_is_deterministic():
    assert review_priority(locus("EXACT_CONCORDANCE", review_required=False)) == NO_REVIEW
    assert review_priority(locus("COMMON_STOP_ALTERNATE_START")) == LOW_PRIORITY_REVIEW
    assert review_priority(locus("BOUNDARY_DISCORDANCE")) == MODERATE_REVIEW
    assert review_priority(locus("CALLER_SPECIFIC")) == MODERATE_REVIEW
    assert review_priority(locus("STRAND_DISCORDANCE", "DISCORDANT")) == HIGH_PRIORITY_REVIEW
    assert review_priority(locus("MERGED_MODEL")) == HIGH_PRIORITY_REVIEW
