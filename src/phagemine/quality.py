from __future__ import annotations

from .models import Protein


def assess(genome: str, proteins: list[Protein]) -> dict:
    """Non-destructive pre-submission QC; warnings are retained for downstream review."""
    warnings: list[dict] = []
    if "N" in genome:
        warnings.append({"code": "ambiguous_bases", "message": f"Genome contains {genome.count('N')} ambiguous bases."})
    ordered = sorted(proteins, key=lambda protein: protein.start)
    for left, right in zip(ordered, ordered[1:]):
        if right.start <= left.end:
            warnings.append({"code": "overlapping_cds", "message": "Predicted CDS coordinates overlap; review gene calls before submission.", "proteins": [left.protein_id, right.protein_id]})
    return {"sequence_length": len(genome), "predicted_cds": len(proteins), "status": "review_required" if warnings else "passed", "warnings": warnings}
