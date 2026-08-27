from __future__ import annotations

from .models import Evidence, EvidenceLevel, Protein


def _category(annotation: str) -> str:
    lowered = annotation.lower()
    if "structural" in lowered:
        return "structural"
    if "dna" in lowered:
        return "DNA metabolism"
    return "uncharacterized"


def mine(proteins: list[Protein], mock: bool = False) -> None:
    """Rank candidates without inventing unavailable evidence; mock signals are fixture-only."""
    for i, protein in enumerate(proteins):
        previous = proteins[i - 1] if i else None
        following = proteins[i + 1] if i + 1 < len(proteins) else None
        unknown = protein.annotation_level in {EvidenceLevel.HYPOTHESIS, EvidenceLevel.WEAK}
        components = {"novelty": 0, "conservation": 0, "context": 0, "functional_signal": 0, "convergence": 0, "quality_penalty": 0}
        if unknown:
            if mock:
                components["novelty"] = 16
                protein.evidence.append(Evidence("novelty", "Mock search found no close characterized homolog; this is unusualness, not functional evidence.", EvidenceLevel.HYPOTHESIS, "mock-phage-evidence", "demo-1", metrics={"characterized_homologs": 0, "cluster_status": "mock lineage-restricted"}))
                components["conservation"] = 13
                protein.evidence.append(Evidence("conservation", "Mock protein-family signal indicates conservation in related phages only.", EvidenceLevel.WEAK, "mock-phage-evidence", "demo-1", metrics={"homologous_phages": 7, "distribution": "lineage-restricted"}))
            neighbors = [p for p in (previous, following) if p]
            categories = [_category(p.annotation) for p in neighbors]
            if mock and any(c != "uncharacterized" for c in categories):
                components["context"] = 14
                protein.evidence.append(Evidence("genomic_context", f"Mock neighborhood contains {', '.join(sorted(set(categories)))} annotation(s); context suggests, but does not prove, association.", EvidenceLevel.WEAK, "mock-context", "demo-1", metrics={"upstream": previous.protein_id if previous else None, "downstream": following.protein_id if following else None, "orientation": protein.strand}))
            if any(e.modality == "domain" and e.supports and e.evidence_strength == "STRONG" for e in protein.evidence):
                components["functional_signal"] = 12
            elif mock:
                # weak sequence signal is only generated for unknown candidates
                components["functional_signal"] = 6
                protein.evidence.append(Evidence("functional_signature", "Mock residue-pattern signal is weak and non-specific.", EvidenceLevel.WEAK, "mock-signatures", "demo-1"))
            protein.alternatives = ["Accessory protein with no currently recognizable family.", "Spurious or miscalled ORF if gene boundaries are inaccurate."]
            protein.missing_evidence = ["External similarity and profile searches against versioned databases.", "Comparative genomic context across related phages.", "Experimental phenotype or biochemical assay."]
        if protein.length < 40:
            components["quality_penalty"] = 15
            has_short_orf_qc = any(
                evidence.modality == "quality"
                and evidence.source == "phagemine-qc"
                and evidence.metrics.get("protein_length") == protein.length
                for evidence in protein.evidence
            )
            if not has_short_orf_qc:
                protein.evidence.append(Evidence(
                    "quality",
                    "Short predicted ORF increases risk of a spurious gene call.",
                    EvidenceLevel.WEAK,
                    "phagemine-qc",
                    "0.1.0",
                    status="real",
                    supports=False,
                    metrics={"protein_length": protein.length},
                ))
        support_modalities = {e.modality for e in protein.evidence if e.supports}
        if unknown:
            components["convergence"] = min(16, max(0, len(support_modalities) - 1) * 4)
        protein.score_components = components
        protein.biological_interest = max(0, min(100, sum(v for k, v in components.items() if k != "quality_penalty") - components["quality_penalty"]))
        n_modalities = len(support_modalities)
        protein.evidence_diversity = "High" if n_modalities >= 4 else "Moderate" if n_modalities >= 2 else "Low" if n_modalities else "None"
        if protein.annotation_level == EvidenceLevel.CURATED:
            protein.functional_confidence = "High"
        elif protein.annotation_level == EvidenceLevel.COMPUTATIONAL:
            protein.functional_confidence = "Moderate"
        elif protein.annotation_level == EvidenceLevel.WEAK:
            protein.functional_confidence = "Low"
        else:
            protein.functional_confidence = "None"


def ranked_candidates(proteins: list[Protein]) -> list[Protein]:
    return sorted((p for p in proteins if p.annotation_level in {EvidenceLevel.HYPOTHESIS, EvidenceLevel.WEAK}), key=lambda p: (-p.biological_interest, p.protein_id))
