"""Shared, serializable gene-model representation.

Step 1 deliberately keeps model selection out of this module.  The model
records preserve both normalized coordinates and the caller's raw values so a
final CDS can be traced back to the exact prediction that produced it.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass, field
from typing import Any


@dataclass(frozen=True)
class GeneModel:
    """A caller-produced CDS model using 1-based inclusive coordinates."""

    # The first five fields retain the v1.1 positional constructor API.
    caller: str
    identifier: str
    start: int
    end: int
    strand: str
    sequence: str = ""
    frame: int | None = None
    caller_version: str | None = None
    command: list[str] | None = None
    options: dict[str, Any] | None = None
    genome_id: str | None = None
    segment_id: str | None = None
    start_codon: str | None = None
    stop_codon: str | None = None
    cds_sequence: str = ""
    protein_sequence: str = ""
    input_sequence_sha256: str | None = None
    raw_start: int | None = None
    raw_end: int | None = None
    raw_strand: str | None = None
    coordinate_system: str = "1-based-inclusive"
    source_record: str | None = None
    source_file: str | None = None
    method_family: str | None = None
    method_lineage: str | None = None

    @property
    def raw_identifier(self) -> str:
        """Canonical v1.2 name for the caller's identifier."""
        return self.identifier

    @property
    def length_nt(self) -> int:
        return self.end - self.start + 1

    @property
    def length_aa(self) -> int:
        return len(self.protein_sequence or self.sequence)

    def to_dict(self) -> dict[str, Any]:
        record = asdict(self)
        record["raw_identifier"] = self.raw_identifier
        record["length_nt"] = self.length_nt
        record["length_aa"] = self.length_aa
        return record


@dataclass
class ReconciledLocus:
    """A locus containing all overlapping candidate models, without selection."""

    locus_id: str
    genome_id: str | None
    segment_id: str | None
    start: int
    end: int
    strand_status: str
    candidate_models: list[GeneModel] = field(default_factory=list)
    supporting_providers: list[str] = field(default_factory=list)
    provider_count: int = 0
    method_family_count: int = 0
    method_lineage_count: int = 0
    supporting_method_families: list[str] = field(default_factory=list)
    supporting_method_lineages: list[str] = field(default_factory=list)
    candidate_count: int = 0
    reconciliation_class: str = "CALLER_SPECIFIC"
    exact_coordinate_groups: list[dict[str, Any]] = field(default_factory=list)
    start_groups: dict[str, list[str]] = field(default_factory=dict)
    stop_groups: dict[str, list[str]] = field(default_factory=dict)
    strand_groups: dict[str, list[str]] = field(default_factory=dict)
    review_required: bool = False
    notes: list[str] = field(default_factory=list)

    def to_dict(self) -> dict[str, Any]:
        return {
            "locus_id": self.locus_id, "genome_id": self.genome_id,
            "segment_id": self.segment_id, "start": self.start, "end": self.end,
            "strand_status": self.strand_status,
            "candidate_models": [model.to_dict() for model in self.candidate_models],
            "supporting_providers": list(self.supporting_providers),
            "provider_count": self.provider_count, "candidate_count": self.candidate_count,
            "method_family_count": self.method_family_count, "method_lineage_count": self.method_lineage_count,
            "supporting_method_families": list(self.supporting_method_families),
            "supporting_method_lineages": list(self.supporting_method_lineages),
            "reconciliation_class": self.reconciliation_class,
            "exact_coordinate_groups": self.exact_coordinate_groups,
            "start_groups": self.start_groups, "stop_groups": self.stop_groups,
            "strand_groups": self.strand_groups, "review_required": self.review_required,
            "notes": list(self.notes),
        }
