"""Sequencing/assembly provenance, intentionally separate from genome representation."""
from __future__ import annotations

from dataclasses import asdict, dataclass, field
from enum import Enum
from typing import Any


class SequencingPlatform(str, Enum):
    ILLUMINA = "ILLUMINA"
    OXFORD_NANOPORE = "OXFORD_NANOPORE"
    PACBIO = "PACBIO"
    HYBRID = "HYBRID"
    SANGER = "SANGER"
    OTHER = "OTHER"
    UNKNOWN = "UNKNOWN"


@dataclass(frozen=True)
class SequencingProvenance:
    """User-provided generation metadata; no field is inferred from an assembled FASTA."""
    sequencing_platform: SequencingPlatform = SequencingPlatform.UNKNOWN
    library_type: str | None = None
    read_type: str | None = None
    assembler: str | None = None
    assembler_version: str | None = None
    polishing_method: str | None = None
    polishing_version: str | None = None
    assembly_method: str | None = None
    assembly_version: str | None = None
    source_information: dict[str, Any] = field(default_factory=dict)
    raw_reads_available: bool | None = None
    metadata_source: str | None = None

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> "SequencingProvenance":
        values = dict(data)
        platform = values.get("sequencing_platform", SequencingPlatform.UNKNOWN)
        values["sequencing_platform"] = platform if isinstance(platform, SequencingPlatform) else SequencingPlatform(str(platform).upper())
        allowed = {item.name for item in __import__("dataclasses").fields(cls)}
        return cls(**{key: value for key, value in values.items() if key in allowed})

    def manifest(self) -> dict[str, Any]:
        data = asdict(self)
        data["sequencing_platform"] = self.sequencing_platform.value
        data["interpretation"] = "Sequencing provenance describes how an assembly was generated. It does not determine topology, orientation, rotation, or biological interpretation."
        return data
