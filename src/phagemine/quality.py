from __future__ import annotations

from .genome import translate
from .models import Protein


def _reverse_complement(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


def _overlap_type(left: Protein, right: Protein) -> tuple[str, int]:
    overlap = left.end - right.start + 1
    if overlap <= 0:
        return "none", 0
    if right.start == left.end:
        return "shared_boundary", 1
    if right.end <= left.end or left.start >= right.start and right.end <= left.end:
        return "nested_overlap", overlap
    return "partial_overlap", overlap


def assess(genome: str, proteins: list[Protein]) -> dict:
    """Non-destructive QC with explicit coordinate/translation auditing.

    CDS coordinates are interpreted as 1-based, inclusive coordinates, matching
    GFF3 and the PHANOTATE textual output. Overlaps are retained because phage
    genomes can legitimately encode overlapping genes; the report distinguishes
    shared-boundary, partial, and nested overlaps instead of treating every
    contact as a coordinate error.
    """
    errors: list[dict] = []
    warnings: list[dict] = []
    infos: list[dict] = []
    if "N" in genome:
        warnings.append({"code": "ambiguous_bases", "message": f"Genome contains {genome.count('N')} ambiguous bases."})

    ordered = sorted(proteins, key=lambda protein: (protein.start, protein.end, protein.strand))
    for protein in proteins:
        prefix = {"protein_id": protein.protein_id, "start": protein.start, "end": protein.end, "strand": protein.strand}
        if not (1 <= protein.start <= protein.end <= len(genome)):
            errors.append({**prefix, "code": "coordinate_out_of_range", "message": "CDS coordinates are outside the analysis sequence."})
            continue
        if protein.strand not in {"+", "-"}:
            errors.append({**prefix, "code": "invalid_strand", "message": "Strand must be + or -."})
            continue
        observed = genome[protein.start - 1:protein.end]
        if protein.strand == "-":
            observed = _reverse_complement(observed)
        if observed != protein.cds:
            errors.append({**prefix, "code": "cds_sequence_mismatch", "message": "Stored CDS sequence does not match the analysis FASTA at the stated coordinates."})
        if len(protein.cds) % 3 != 0:
            errors.append({**prefix, "code": "cds_length_not_multiple_of_three", "message": "CDS length is not divisible by three."})
        expected_translation = translate(protein.cds)
        if protein.sequence != expected_translation:
            errors.append({**prefix, "code": "translation_mismatch", "message": "Stored protein sequence does not match translation of the stored CDS."})
        if "*" in expected_translation:
            errors.append({**prefix, "code": "internal_stop_codon", "message": "Translation contains an internal stop codon."})
        touches_left_boundary = protein.start == 1
        touches_right_boundary = protein.end == len(genome)
        if touches_left_boundary or touches_right_boundary:
            infos.append({**prefix, "code": "boundary_adjacent_cds", "message": "CDS touches a sequence boundary; absence of a canonical start/stop at the truncated end is not treated as a gene-model error."})
        if len(protein.cds) >= 6:
            start_codon = protein.cds[:3]
            stop_codon = protein.cds[-3:]
            if start_codon not in {"ATG", "GTG", "TTG"} and not touches_right_boundary:
                warnings.append({**prefix, "code": "noncanonical_start_codon", "message": f"Predicted CDS begins with {start_codon}; review whether this is supported by the gene caller."})
            if stop_codon not in {"TAA", "TAG", "TGA"} and not touches_left_boundary:
                warnings.append({**prefix, "code": "missing_terminal_stop_codon", "message": f"Predicted CDS ends with {stop_codon}, not a standard stop codon."})
        else:
            warnings.append({**prefix, "code": "very_short_cds", "message": "CDS is shorter than six nucleotides and cannot support start/stop auditing."})

    for index, left in enumerate(ordered):
        for right in ordered[index + 1:]:
            if right.start > left.end:
                break
            overlap_type, overlap = _overlap_type(left, right)
            if overlap_type == "shared_boundary":
                infos.append({"code": "shared_cds_boundary", "message": "Adjacent CDS calls share one endpoint base in 1-based inclusive coordinates; retained as a gene-caller observation, not treated as a coordinate error.", "proteins": [left.protein_id, right.protein_id], "overlap_bases": overlap, "strand_pair": [left.strand, right.strand]})
            elif overlap_type in {"partial_overlap", "nested_overlap"}:
                warnings.append({"code": "overlapping_cds", "message": f"{overlap_type.replace('_', ' ')} retained from the gene caller; phage genomes can contain legitimate overlapping CDSs, but review is recommended before submission.", "proteins": [left.protein_id, right.protein_id], "overlap_bases": overlap, "strand_pair": [left.strand, right.strand], "same_strand": left.strand == right.strand})

    status = "failed" if errors else ("review_required" if warnings else "passed")
    return {"sequence_length": len(genome), "predicted_cds": len(proteins), "status": status, "errors": errors, "warnings": warnings, "info": infos, "coordinate_system": "1-based-inclusive", "coordinate_scope": "analysis_sequence"}
