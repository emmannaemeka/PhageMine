from __future__ import annotations

import hashlib
from pathlib import Path


def read_fasta(path: str | Path) -> tuple[str, str]:
    records: list[tuple[str, str]] = []
    header: str | None = None
    seq: list[str] = []
    for raw in Path(path).read_text().splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                records.append((header, "".join(seq)))
            header, seq = line[1:].split()[0], []
        else:
            seq.append(line.upper())
    if header is not None:
        records.append((header, "".join(seq)))
    if len(records) != 1:
        raise ValueError("MVP accepts exactly one FASTA record per run")
    genome_id, sequence = records[0]
    # Accept standard IUPAC DNA ambiguity codes; preserve the sequence exactly
    # for provenance while allowing annotated public cohorts with ambiguous
    # bases (e.g. Y) to proceed through validation.
    if not sequence or set(sequence) - set("ACGTRYSWKMBDHVN"):
        raise ValueError("FASTA sequence must contain standard IUPAC DNA bases")
    return genome_id, sequence


def checksum(path: str | Path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()
