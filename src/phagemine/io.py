from __future__ import annotations

import hashlib
from pathlib import Path

IUPAC_NUCLEOTIDES = frozenset("ACGTRYSWKMBDHVN")
USABLE_NUCLEOTIDES = frozenset("ACGTN")


def normalize_sequence(sequence: str, *, record_id: str = "", provider_id: str = "") -> tuple[str, dict]:
    """Return a provider-safe sequence and an auditable normalization record.

    IUPAC ambiguity symbols and isolated unsupported characters are replaced by
    ``N``.  Empty/non-nucleotide records are rejected; otherwise the original
    sequence is never silently discarded and the replacement details are
    returned for provenance.
    """
    compact = "".join(str(sequence).split()).upper()
    if not compact:
        raise ValueError(f"FASTA record {record_id or '<unknown>'} is empty")
    usable = [base for base in compact if base in IUPAC_NUCLEOTIDES]
    if not usable:
        raise ValueError(f"FASTA record {record_id or '<unknown>'} contains no usable nucleotide symbols")
    replacements = {}
    normalized_chars = []
    for base in compact:
        if base in USABLE_NUCLEOTIDES:
            normalized_chars.append(base)
        else:
            normalized_chars.append("N")
            replacements[base] = replacements.get(base, 0) + 1
    normalized = "".join(normalized_chars)
    audit = {
        "record_id": record_id,
        "provider_id": provider_id,
        "original_sequence": compact,
        "original_length": len(compact),
        "normalized_length": len(normalized),
        "replaced_symbol_counts": dict(sorted(replacements.items())),
        "replacement_count": sum(replacements.values()),
        "changed": normalized != compact,
        "original_sequence_sha256": hashlib.sha256(compact.encode()).hexdigest(),
        "normalized_sequence_sha256": hashlib.sha256(normalized.encode()).hexdigest(),
    }
    return normalized, audit


def read_fasta_records(path: str | Path) -> list[tuple[str, str]]:
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
            fields = line[1:].split()
            if not fields:
                raise ValueError("FASTA record identifier is missing")
            header, seq = fields[0], []
        else:
            if header is None:
                raise ValueError("FASTA sequence data appeared before the first record identifier")
            seq.append("".join(line.split()).upper())
    if header is not None:
        records.append((header, "".join(seq)))
    for identifier, sequence in records:
        # Validation at the input boundary distinguishes empty/non-nucleotide
        # records from provider-incompatible but otherwise usable symbols.
        normalize_sequence(sequence, record_id=identifier)
    return records


def read_fasta(path: str | Path) -> tuple[str, str]:
    records = read_fasta_records(path)
    if len(records) != 1:
        raise ValueError("MVP accepts exactly one FASTA record per run")
    genome_id, sequence = records[0]
    # Accept standard IUPAC DNA ambiguity codes; preserve the sequence exactly
    # for provenance while allowing annotated public cohorts with ambiguous
    # bases (e.g. Y) to proceed through validation.
    return genome_id, sequence


def checksum(path: str | Path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()
