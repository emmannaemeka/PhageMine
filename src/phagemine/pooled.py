"""Exact-sequence pan-proteome utilities for pooled evidence execution."""
from __future__ import annotations

import hashlib
from collections import defaultdict
from dataclasses import asdict
from typing import Iterable

from .models import Protein


def sequence_sha256(sequence: str) -> str:
    return hashlib.sha256(sequence.encode("ascii")).hexdigest()


def deduplicate_proteins(proteins: Iterable[Protein]) -> tuple[list[Protein], dict[str, list[dict]]]:
    """Return first-occurrence representatives and a complete occurrence map."""
    representatives: list[Protein] = []
    occurrences: dict[str, list[dict]] = defaultdict(list)
    by_sequence: dict[str, Protein] = {}
    for protein in proteins:
        digest = sequence_sha256(protein.sequence)
        occurrences[digest].append({"genome_id": protein.genome_id, "protein_id": protein.protein_id,
                                     "start": protein.start, "end": protein.end, "strand": protein.strand,
                                     "sequence_sha256": digest})
        if digest not in by_sequence:
            by_sequence[digest] = protein
            representatives.append(protein)
    return representatives, dict(occurrences)


def remap_evidence(evidence_by_representative: dict[str, list], representatives: list[Protein], occurrences: dict[str, list[dict]]) -> dict[str, list]:
    """Copy pooled evidence losslessly to every exact-sequence occurrence."""
    digest_by_id = {p.protein_id: sequence_sha256(p.sequence) for p in representatives}
    result: dict[str, list] = {}
    for representative_id, evidence in evidence_by_representative.items():
        digest = digest_by_id[representative_id]
        for occurrence in occurrences.get(digest, []):
            result[occurrence["protein_id"]] = [dict(item, provenance={**item.get("provenance", {}), "pooled": True,
                    "pooled_representative": representative_id, "sequence_sha256": digest}) if isinstance(item, dict) else item for item in evidence]
    return result
