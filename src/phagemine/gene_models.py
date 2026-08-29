"""Shared, serializable gene-model representation.

Step 1 deliberately keeps model selection out of this module.  The model
records preserve both normalized coordinates and the caller's raw values so a
final CDS can be traced back to the exact prediction that produced it.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
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
