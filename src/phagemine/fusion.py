"""Deterministic, conservative fusion of acquired evidence."""
from __future__ import annotations

import csv
import json
import re
from pathlib import Path
from typing import Any

from .models import Evidence, Protein

FUSION_RULES_VERSION = "1.0"
MODALITY_FAMILIES = {
    "Pfam": "DOMAIN", "domain": "DOMAIN",
    "Swiss-Prot": "CURATED_SEQUENCE_HOMOLOGY", "sequence_homology": "CURATED_SEQUENCE_HOMOLOGY",
    "VOGDB": "VIRAL_ORTHOLOGY", "viral_orthology": "VIRAL_ORTHOLOGY",
    "PHROGs": "PHAGE_ORTHOLOGY", "phage_orthology": "PHAGE_ORTHOLOGY",
}
INCOMPATIBLE_FUNCTION_GROUPS = (
    {"integrase", "recombinase", "excisionase"},
    {"capsid", "portal protein", "terminase", "tail protein"},
    {"holin", "endolysin", "lysis"},
)
UNKNOWN_LABELS = {
    "", "hypothetical", "hypothetical protein", "unknown protein", "unknown function",
    "protein of unknown function", "uncharacterized protein", "uncharacterised protein",
    "conserved hypothetical protein", "conserved protein of unknown function",
    "putative uncharacterized protein", "no annotation", "none", "null",
}


def normalize_function(description: str | None) -> str | None:
    if description is None:
        return None
    value = re.sub(r"\s+", " ", description.strip().lower())
    value = re.sub(r"^refseq\s+", "", value)
    value = re.sub(r"^sp\|[^|]+\|", "", value)
    value = re.sub(r"\s*\{eco:[^}]+\}\s*$", "", value, flags=re.IGNORECASE).strip()
    return None if value in UNKNOWN_LABELS else value


def _category(evidence: Evidence) -> str | None:
    value = evidence.metrics.get("functional_category") or evidence.metrics.get("category")
    if value is None or str(value).strip().lower() in UNKNOWN_LABELS | {"xu", "unknown"}:
        return None
    return str(value).strip()


def _conservation(evidence: list[Evidence]) -> str:
    accepted = [e for e in evidence if e.supports and e.source in {"PHROGs", "VOGDB"}]
    strong = [e for e in accepted if e.evidence_strength == "STRONG"]
    if strong:
        return "STRONGLY_CONSERVED"
    if accepted:
        return "CONSERVED"
    return "NOT_ESTABLISHED"


def classify_protein(protein: Protein) -> dict[str, Any]:
    accepted = [e for e in protein.evidence if e.supports]
    informative = [(i, e, normalize_function(e.description)) for i, e in enumerate(accepted) if normalize_function(e.description)]
    strong_info = [(i, e, label) for i, e, label in informative if e.evidence_strength == "STRONG"]
    categories = sorted({_category(e) for e in accepted if _category(e)})
    support_ids = [f"{e.source}:{e.identifier or e.family_name or i}" for i, e in enumerate(accepted)]
    supporting_sources = sorted({e.source for _, e, _ in informative})
    supporting_modalities = sorted({MODALITY_FAMILIES.get(e.source, MODALITY_FAMILIES.get(e.modality, e.modality)) for _, e, _ in informative})
    conflicting = []
    labels_by_source: dict[str, set[str]] = {}
    for _, e, label in strong_info:
        labels_by_source.setdefault(e.source, set()).add(label)
    strong_labels = {label for _, _, label in strong_info}
    compatible = all(a == b or a in b or b in a for a in strong_labels for b in strong_labels)
    explicit_incompatibility = any(any(any(term in label for term in group) for label in strong_labels) and
                                   any(any(term in label for term in other) for label in strong_labels)
                                   for index, group in enumerate(INCOMPATIBLE_FUNCTION_GROUPS)
                                   for other in INCOMPATIBLE_FUNCTION_GROUPS[index + 1:])
    ambiguity = []
    if len(strong_labels) > 1 and not compatible and not explicit_incompatibility:
        ambiguity.append("distinct strong labels could not be confidently established as biologically incompatible")
    if len(strong_labels) > 1 and len(labels_by_source) > 1 and explicit_incompatibility:
        conflicting = [(i, e, label) for i, e, label in strong_info]
    conflict_ids = [f"{e.source}:{e.identifier or e.family_name or i}" for i, e, _ in conflicting]
    conflict_sources = sorted({e.source for _, e, _ in conflicting})
    conflict_descriptions = [e.description for _, e, _ in conflicting if e.description]
    labels = sorted({label for _, _, label in informative})
    proposed = labels[0] if len(labels) == 1 else None
    independent_strong = len({e.source for _, e, _ in strong_info}) >= 2
    swiss_strong = any(e.source == "Swiss-Prot" for _, e, _ in strong_info)
    strong_families = {MODALITY_FAMILIES.get(e.source, MODALITY_FAMILIES.get(e.modality, e.modality)) for _, e, _ in strong_info}
    orthology_only = strong_families and strong_families <= {"PHAGE_ORTHOLOGY", "VIRAL_ORTHOLOGY"}
    agreeing_sources = len({e.source for _, e, label in informative if proposed and label == proposed})
    if conflicting:
        state, confidence, reasons = "CONFLICTING_EVIDENCE", "LOW", ["accepted STRONG informative evidence supports incompatible normalized functions"]
    elif strong_info and (independent_strong or (swiss_strong and agreeing_sources >= 2)) and len(strong_families) >= 2 and not orthology_only and not ambiguity:
        state, confidence = "KNOWN_FUNCTION", "HIGH"
        reasons = [f"STRONG {e.source} informative evidence" for _, e, _ in strong_info]
        if independent_strong:
            reasons.append("independent agreeing source corroboration")
        reasons.append("no accepted strong conflict")
    elif strong_info or informative:
        state, confidence = "PROBABLE_FUNCTION", "MODERATE" if strong_info else "LOW"
        reasons = [f"accepted {e.evidence_strength or 'informative'} {e.source} evidence supports a functional interpretation"]
    elif categories:
        state, confidence = "FUNCTIONAL_CLASS_ONLY", "LOW"
        reasons = ["accepted evidence supports a broad functional category without an informative specific function"]
    elif _conservation(protein.evidence) in {"STRONGLY_CONSERVED", "CONSERVED"}:
        state, confidence = "CONSERVED_UNKNOWN", "LOW"
        reasons = ["accepted PHROGs/VOGDB orthology supports conservation but no accepted informative functional evidence is present"]
    else:
        state, confidence, reasons = "UNRESOLVED", "NONE", ["no accepted evidence establishes a function, category, or meaningful conservation"]
    return {
        "protein_id": protein.protein_id, "functional_state": state,
        "proposed_function": proposed, "normalized_function": proposed,
        "functional_category": categories[0] if len(categories) == 1 else ("; ".join(categories) if categories else None),
        "confidence": confidence, "confidence_reasons": reasons,
        "conservation_status": _conservation(protein.evidence),
        "supporting_sources": supporting_sources, "supporting_modalities": supporting_modalities,
        "supporting_source_count": len(supporting_sources), "supporting_modality_count": len(supporting_modalities), "supporting_record_count": len(informative),
        "supporting_evidence_ids": support_ids,
        "conflicting_sources": conflict_sources, "conflicting_evidence_ids": conflict_ids,
        "conflict_descriptions": conflict_descriptions, "ambiguity_flags": ambiguity,
        "reasoning_summary": "; ".join(reasons) + ".", "fusion_rules_version": FUSION_RULES_VERSION,
    }


def classify_proteins(proteins: list[Protein]) -> list[dict[str, Any]]:
    return [classify_protein(protein) for protein in proteins]


def write_classification(root: str | Path, proteins: list[Protein], results: list[dict[str, Any]] | None = None) -> None:
    results = results if results is not None else classify_proteins(proteins)
    by_id = {p.protein_id: p for p in proteins}
    path = Path(root)
    (path / "functional_classification.json").write_text(json.dumps(results, indent=2, sort_keys=True))
    columns = ["protein_id", "start", "end", "strand", "length_aa", "functional_state", "proposed_function", "functional_category", "confidence", "conservation_status", "supporting_sources", "supporting_source_count", "supporting_record_count", "conflict", "reasoning_summary"]
    with (path / "functional_classification.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t")
        writer.writeheader()
        for result in results:
            protein = by_id[result["protein_id"]]
            writer.writerow({**{key: result.get(key) for key in columns}, "start": protein.start, "end": protein.end, "strand": protein.strand, "length_aa": protein.length, "supporting_sources": ";".join(result["supporting_sources"]), "conflict": bool(result["conflicting_evidence_ids"])})
