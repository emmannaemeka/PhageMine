"""Explicit, provenance-preserving views of a genome assembly.

The original assembly remains authoritative. Transformations are only made through
explicit method calls and every analysis coordinate refers to ``analysis_sequence``.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass, field
from enum import Enum
from typing import Any
import hashlib


class Topology(str, Enum):
    LINEAR = "LINEAR"
    CIRCULAR = "CIRCULAR"
    UNKNOWN = "UNKNOWN"


class Orientation(str, Enum):
    ORIGINAL = "ORIGINAL"
    REVERSE_COMPLEMENT = "REVERSE_COMPLEMENT"
    UNRESOLVED = "UNRESOLVED"


class Rotation(str, Enum):
    NONE = "NONE"
    ROTATED = "ROTATED"
    UNRESOLVED = "UNRESOLVED"


@dataclass(frozen=True)
class SegmentRecord:
    """A segment-local sequence record; segments are never concatenated."""
    segment_id: str
    sequence: str
    source_record: str | None = None

    @property
    def length(self) -> int:
        return len(self.sequence)

    @property
    def sha256(self) -> str:
        return hashlib.sha256(self.sequence.encode()).hexdigest()


@dataclass(frozen=True)
class GenomeRecord:
    """Biological genome container preserving independent segment identity."""
    genome_id: str
    molecule_type: str
    segments: tuple[SegmentRecord, ...]

    def __post_init__(self):
        if self.molecule_type.lower() not in {"dna", "rna"}:
            raise ValueError("molecule_type must be dna or rna")
        if not self.segments:
            raise ValueError("GenomeRecord requires at least one segment")


def reverse_complement(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


@dataclass(frozen=True)
class TransformEvent:
    operation: str
    parameters: dict[str, Any]
    rationale: str | None = None
    evidence: list[dict[str, Any]] = field(default_factory=list)


@dataclass(frozen=True)
class GenomeRepresentation:
    original_sequence_id: str
    analysis_sequence_id: str
    original_sequence: str
    analysis_sequence: str
    topology: Topology = Topology.UNKNOWN
    orientation: Orientation = Orientation.ORIGINAL
    rotation: Rotation = Rotation.NONE
    transform_history: tuple[TransformEvent, ...] = ()
    evidence: tuple[dict[str, Any], ...] = ()
    reference: dict[str, Any] | None = None

    @classmethod
    def original(cls, sequence_id: str, sequence: str, topology: Topology = Topology.UNKNOWN, evidence: list[dict[str, Any]] | None = None, reference: dict[str, Any] | None = None) -> "GenomeRepresentation":
        return cls(sequence_id, sequence_id, sequence, sequence, topology=topology, evidence=tuple(evidence or ()), reference=reference)

    def with_reverse_complement(self, rationale: str, evidence: list[dict[str, Any]] | None = None, reference: dict[str, Any] | None = None) -> "GenomeRepresentation":
        if not rationale:
            raise ValueError("An explicit rationale is required before reverse-complementing an analysis representation.")
        return GenomeRepresentation(self.original_sequence_id, f"{self.original_sequence_id}__rc", self.original_sequence, reverse_complement(self.analysis_sequence), self.topology, Orientation.REVERSE_COMPLEMENT if self.orientation == Orientation.ORIGINAL else Orientation.ORIGINAL, self.rotation, self.transform_history + (TransformEvent("reverse_complement", {}, rationale, evidence or []),), self.evidence + tuple(evidence or ()), reference if reference is not None else self.reference)

    def with_rotation(self, analysis_origin: int, rationale: str, evidence: list[dict[str, Any]] | None = None, reference: dict[str, Any] | None = None) -> "GenomeRepresentation":
        if self.topology != Topology.CIRCULAR:
            raise ValueError("Circular rotation requires topology=CIRCULAR; do not rotate an UNKNOWN or LINEAR representation.")
        if not 1 <= analysis_origin <= len(self.analysis_sequence):
            raise ValueError("Rotation origin must be a one-based coordinate within the current analysis representation.")
        if not rationale:
            raise ValueError("An explicit rationale is required before rotating an analysis representation.")
        offset = analysis_origin - 1
        sequence = self.analysis_sequence[offset:] + self.analysis_sequence[:offset]
        return GenomeRepresentation(self.original_sequence_id, f"{self.original_sequence_id}__rot{analysis_origin}", self.original_sequence, sequence, self.topology, self.orientation, Rotation.ROTATED, self.transform_history + (TransformEvent("rotate", {"analysis_origin": analysis_origin, "coordinate_scope": "pre-rotation analysis representation"}, rationale, evidence or []),), self.evidence + tuple(evidence or ()), reference if reference is not None else self.reference)

    def manifest(self) -> dict[str, Any]:
        return {"original_sequence_id": self.original_sequence_id, "analysis_sequence_id": self.analysis_sequence_id, "topology": self.topology.value, "orientation": self.orientation.value, "rotation": self.rotation.value, "transform_history": [asdict(event) for event in self.transform_history], "evidence": list(self.evidence), "reference": self.reference, "coordinate_scope": "All PhageMine gene, protein, neighborhood, and annotation coordinates refer to analysis_sequence_id."}
