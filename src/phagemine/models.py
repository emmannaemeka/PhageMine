from __future__ import annotations

from dataclasses import asdict, dataclass, field
from enum import Enum
from typing import Any


class EvidenceLevel(str, Enum):
    EXPERIMENTAL = "experimentally established"
    CURATED = "curated database-supported"
    COMPUTATIONAL = "computational prediction"
    WEAK = "weak inference"
    HYPOTHESIS = "hypothesis requiring experimental validation"


@dataclass
class SubmissionMetadata:
    """User-supplied source and submitter information; absent fields are never inferred."""
    sequence_id: str | None = None
    organism: str | None = None
    phage_name: str | None = None
    host: str | None = None
    isolation_source: str | None = None
    collection_date: str | None = None
    geographic_location: str | None = None
    isolate: str | None = None
    strain: str | None = None
    sequencing_technology: str | None = None
    submitter_name: str | None = None
    submitter_email: str | None = None
    authors: list[str] = field(default_factory=list)
    bioproject_accession: str | None = None
    biosample_accession: str | None = None

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> "SubmissionMetadata":
        allowed = {item.name for item in __import__("dataclasses").fields(cls)}
        return cls(**{key: value for key, value in data.items() if key in allowed})

    def missing_required(self) -> list[str]:
        required = ("organism", "phage_name", "host", "isolation_source", "collection_date", "geographic_location", "sequencing_technology", "submitter_name", "submitter_email")
        missing = [name for name in required if not getattr(self, name)]
        if not self.authors:
            missing.append("authors")
        return missing


@dataclass
class Evidence:
    modality: str
    statement: str
    level: EvidenceLevel
    source: str
    source_version: str
    status: str = "mock"
    supports: bool = True
    metrics: dict[str, Any] = field(default_factory=dict)
    identifier: str | None = None
    threshold: dict[str, Any] = field(default_factory=dict)
    coordinates: dict[str, int] | None = None
    provenance: dict[str, Any] = field(default_factory=dict)
    evidence_strength: str | None = None
    family_name: str | None = None
    description: str | None = None


@dataclass
class Protein:
    genome_id: str
    protein_id: str
    start: int
    end: int
    strand: str
    cds: str
    sequence: str
    gene_call_source: str
    locus_tag: str | None = None
    gene_call_parameters: dict[str, Any] = field(default_factory=dict)
    start_codon: str | None = None
    stop_codon: str | None = None
    annotation: str = "Hypothetical protein"
    annotation_level: EvidenceLevel = EvidenceLevel.HYPOTHESIS
    evidence: list[Evidence] = field(default_factory=list)
    alternatives: list[str] = field(default_factory=list)
    missing_evidence: list[str] = field(default_factory=list)
    functional_confidence: str = "None"
    biological_interest: int = 0
    evidence_diversity: str = "None"
    score_components: dict[str, int] = field(default_factory=dict)

    @property
    def length(self) -> int:
        return len(self.sequence)

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)
