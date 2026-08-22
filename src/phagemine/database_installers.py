"""Manifest-driven database installers.

Downloads are intentionally opt-in; ``--dry-run`` never touches the network.
The release specifications document the transactional preparation contract used
by production installers and provide a stable plan for automation.
"""
from __future__ import annotations
from dataclasses import dataclass

@dataclass(frozen=True)
class DatabaseSpec:
    name: str
    compressed_gib: float
    preparation: str

SPECS = (
    DatabaseSpec("pfam", 0.25, "decompress; hmmpress; checksum; transactional install; register; manifest; readiness validation"),
    DatabaseSpec("vogdb", 0.35, "concatenate HMMs; hmmpress; retain annotations; checksum; transactional install; register; manifest; readiness validation"),
    DatabaseSpec("swissprot", 0.90, "build DIAMOND database; retain DAT metadata; checksum; transactional install; register; manifest; readiness validation"),
    DatabaseSpec("phrogs", 0.40, "install prepared MMseqs2 database and annotations; checksum; transactional install; register; manifest; readiness validation"),
    DatabaseSpec("pmfdb", 0.60, "build INPHARED-derived protein reference and MMseqs2 index; retain predicted annotations as computational evidence; checksum; transactional install; register; manifest; readiness validation"),
    DatabaseSpec("inphared", 0.90, "install whole-genome reference; mash sketch -i per genome; checksum; transactional install; register; manifest; readiness validation"),
)

def plan(selected: str | None = None) -> list[dict[str, object]]:
    specs = SPECS if selected in (None, "all") else tuple(s for s in SPECS if s.name == selected)
    if selected not in (None, "all") and not specs:
        raise ValueError(f"unknown database: {selected}")
    return [{"database": s.name, "compressed_gib": s.compressed_gib, "preparation": s.preparation} for s in specs]

def dry_run(selected: str | None = "all") -> dict[str, object]:
    items = plan(selected)
    return {"databases": items, "count": len(items), "compressed_download_gib": round(sum(float(i["compressed_gib"]) for i in items), 1), "dry_run": True}
