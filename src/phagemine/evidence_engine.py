"""Evidence-aware adjudication for the next-generation annotation engine.

This module is intentionally policy-first.  It consumes already-produced
evidence records and caller models, preserves every raw record, and emits
auditable interpretations.  PHANOTATE remains the final CDS model; secondary
callers are observational only.  No score is treated as a probability and no
missing evidence is converted into biological absence.
"""
from __future__ import annotations

import csv
import hashlib
import json
import re
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

STRUCTURAL_CONFIDENCE = ("HIGH", "MODERATE", "LOW", "UNRESOLVED")
STRUCTURAL_RELATIONSHIPS = (
    "EXACT_MATCH", "SAME_STOP_DIFFERENT_START", "NEAR_BOUNDARY_MATCH",
    "PHANOTATE_ONLY", "PYRODIGAL_GV_ONLY", "CONFLICTING_ORF",
    "OVERLAPPING_ALTERNATIVE",
)
FUNCTIONAL_CLASSES = (
    "HIGH_CONFIDENCE_SPECIFIC", "SUPPORTED_SPECIFIC", "SUPPORTED_GENERAL",
    "DOMAIN_ONLY", "CONSERVED_UNKNOWN", "HYPOTHETICAL", "CONFLICTING",
)
MODULE_STATUSES = (
    "ESTABLISHED", "STRONGLY_SUPPORTED", "SUPPORTED", "CANDIDATE",
    "EVIDENCE_PRESENT_REVIEW_REQUIRED", "NOT_ESTABLISHED", "NOT_ASSESSABLE",
)
ARCHITECTURES = (
    "TAILED_DSDNA_LIKE", "FILAMENTOUS_PHAGE_LIKE", "SMALL_SSDNA_PHAGE_LIKE",
    "OTHER_PHAGE_ARCHITECTURE", "MIXED_OR_CONFLICTING",
    "ARCHITECTURE_UNRESOLVED", "NOT_ASSESSABLE",
)
COMPARATIVE_STATES = (
    "SIGNIFICANT_CLOSE_REFERENCE", "SIGNIFICANT_DISTANT_REFERENCE",
    "WEAK_REFERENCE", "NO_SIGNIFICANT_REFERENCE", "INSUFFICIENT_DATA",
    "DATABASE_UNAVAILABLE", "TOOL_UNAVAILABLE", "FAILED",
)
_UNKNOWN = re.compile(r"\b(hypothetical|unknown|uncharacteri[sz]ed|unannotated|predicted protein)\b", re.I)


@dataclass(frozen=True)
class StructuralAdjudication:
    locus_id: str
    relationship: str
    structural_cds_confidence: str
    phanotate_coordinates: tuple[int, int, str] | None
    pyrodigal_gv_coordinates: tuple[int, int, str] | None
    reasoning: str
    provenance: dict[str, Any] = field(default_factory=dict)


@dataclass(frozen=True)
class FunctionalAdjudication:
    locus_id: str
    functional_class: str
    functional_assignment_confidence: str
    annotation_specificity: str
    accepted_evidence: tuple[dict[str, Any], ...]
    weak_evidence: tuple[dict[str, Any], ...]
    rejected_evidence: tuple[dict[str, Any], ...]
    conflicting_evidence: tuple[dict[str, Any], ...]
    reasoning: str


def _coord(model: Any) -> tuple[int, int, str]:
    if isinstance(model, Mapping):
        return (int(model["start"]), int(model["end"]), str(model.get("strand", "+")))
    return (int(getattr(model, "start")), int(getattr(model, "end")), str(getattr(model, "strand", "+")))


def _length(coord: tuple[int, int, str]) -> int:
    return abs(coord[1] - coord[0]) + 1


def _overlap(a: tuple[int, int, str], b: tuple[int, int, str]) -> int:
    if a[2] != b[2]:
        return 0
    return max(0, min(a[1], b[1]) - max(a[0], b[0]) + 1)


def _relationship(primary: tuple[int, int, str], secondary: tuple[int, int, str] | None) -> str:
    if secondary is None:
        return "PHANOTATE_ONLY"
    if primary == secondary:
        return "EXACT_MATCH"
    if primary[1] == secondary[1] and primary[0] != secondary[0]:
        return "SAME_STOP_DIFFERENT_START"
    overlap = _overlap(primary, secondary)
    if overlap and overlap / min(_length(primary), _length(secondary)) >= 0.8:
        return "NEAR_BOUNDARY_MATCH"
    if overlap:
        return "OVERLAPPING_ALTERNATIVE"
    return "CONFLICTING_ORF"


def _support_metrics(evidence: Iterable[Mapping[str, Any]]) -> tuple[int, set[str], bool, bool]:
    records = list(evidence)
    sources: set[str] = set()
    strong = False
    independent = False
    for item in records:
        source = str(item.get("source_database") or item.get("source") or item.get("database") or "").strip()
        if source:
            sources.add(source)
        status = str(item.get("status", "")).upper()
        role = str(item.get("role", item.get("supports", ""))).upper()
        if status in {"REJECTED", "CONFLICTING"} or role in {"CONFLICT", "CONFLICTING", "REJECTED"}:
            continue
        strength = str(item.get("evidence_strength", item.get("strength", ""))).upper()
        if strength in {"STRONG", "HIGH", "EXPERIMENTAL", "CURATED"} or item.get("qualified") is True:
            strong = True
        if source and source.lower() not in {s.lower() for s in sources if s != source}:
            independent = independent or source.lower() in {"swiss-prot", "swissprot", "phrogs", "vogdb", "comparative", "synteny", "context"}
    return len(records), sources, strong, independent


def adjudicate_structure(locus_id: str, phanotate: Any, pyrodigal_gv: Any | None,
                         independent_support: Sequence[Mapping[str, Any]] = (), *,
                         provenance: Mapping[str, Any] | None = None) -> StructuralAdjudication:
    """Adjudicate a locus without replacing the PHANOTATE coordinates."""
    p = _coord(phanotate); g = _coord(pyrodigal_gv) if pyrodigal_gv is not None else None
    relationship = _relationship(p, g)
    n, sources, strong, independent = _support_metrics(independent_support)
    if relationship == "EXACT_MATCH":
        confidence, reason = "HIGH", "PHANOTATE and Pyrodigal-gv agree exactly."
    elif strong and (independent or n >= 2):
        confidence, reason = "MODERATE", "Independent biological support corroborates the PHANOTATE model despite caller disagreement."
    elif relationship == "PHANOTATE_ONLY" and not independent_support:
        confidence, reason = "UNRESOLVED", "PHANOTATE is retained as the legacy model but lacks independent structural/biological support."
    elif relationship in {"CONFLICTING_ORF", "OVERLAPPING_ALTERNATIVE"}:
        confidence, reason = "LOW", "Alternative ORF is retained observationally; coordinate conflict requires review."
    else:
        confidence, reason = "LOW", "Boundary disagreement is retained and length alone was not used to reject the PHANOTATE model."
    return StructuralAdjudication(locus_id, relationship, confidence, p, g, reason, dict(provenance or {}, phanotate_final=True, secondary_caller="prodigal_gv", evidence_sources=sorted(sources)))


def adjudicate_function(locus_id: str, product: str | None, evidence: Sequence[Mapping[str, Any]], *, structural_confidence: str = "UNRESOLVED") -> FunctionalAdjudication:
    accepted, weak, rejected, conflicting = [], [], [], []
    for item in evidence:
        status = str(item.get("status", "ACCEPTED")).upper()
        strength = str(item.get("evidence_strength", item.get("strength", ""))).upper()
        if status in {"REJECTED", "UNACCEPTED"}:
            rejected.append(dict(item)); continue
        if status in {"CONFLICTING", "CONFLICT"} or str(item.get("role", "")).upper() in {"CONFLICT", "CONFLICTING"}:
            conflicting.append(dict(item)); continue
        (accepted if strength in {"STRONG", "HIGH", "EXPERIMENTAL", "CURATED"} else weak).append(dict(item))
    text = product or ""
    if conflicting and not accepted:
        cls, confidence, specificity, reason = "CONFLICTING", "LOW", "CONFLICTING", "Conflicting evidence was preserved and no concordant qualified support resolved it."
    elif not text or _UNKNOWN.search(text):
        cls, confidence, specificity, reason = "HYPOTHETICAL", "LOW", "UNKNOWN", "No specific product was asserted."
    elif accepted and any(str(x.get("evidence_type", x.get("type", ""))).lower() == "domain" for x in accepted) and not any(str(x.get("evidence_type", x.get("type", ""))).lower() in {"sequence_similarity", "curated_homology", "experimental"} for x in accepted):
        cls, confidence, specificity, reason = "DOMAIN_ONLY", "LOW", "DOMAIN", "Domain evidence is retained but cannot establish a whole-protein specific function."
    elif accepted and all(str(x.get("evidence_type", x.get("type", ""))).lower() in {"family", "orthology", "profile"} for x in accepted):
        cls, confidence, specificity, reason = "SUPPORTED_SPECIFIC", "MODERATE", "SPECIFIC", "Family/orthology evidence supports a candidate function but does not receive high specific-function confidence by itself."
    elif accepted and structural_confidence == "HIGH":
        cls, confidence, specificity, reason = "HIGH_CONFIDENCE_SPECIFIC", "HIGH", "SPECIFIC", "Qualified concordant evidence supports the asserted specificity."
    elif accepted:
        cls, confidence, specificity, reason = "SUPPORTED_SPECIFIC", "MODERATE", "SPECIFIC", "Qualified evidence supports a functional assertion; family evidence was not treated as sufficient alone."
    elif weak:
        cls, confidence, specificity, reason = "SUPPORTED_GENERAL", "LOW", "GENERAL", "Only weak or partial evidence is available; a broad interpretation is safer."
    else:
        cls, confidence, specificity, reason = "CONSERVED_UNKNOWN", "LOW", "UNKNOWN", "No qualified functional evidence was available."
    return FunctionalAdjudication(locus_id, cls, confidence, specificity, tuple(accepted), tuple(weak), tuple(rejected), tuple(conflicting), reason)


def aggregate_module_evidence(loci: Sequence[Mapping[str, Any]], *, modules: Sequence[str] | None = None) -> list[dict[str, Any]]:
    """Aggregate qualified evidence without counting repeated database hits as votes."""
    names = tuple(modules or ("DNA_PACKAGING", "HEAD_CAPSID", "PORTAL", "TAIL", "TAPE_MEASURE", "BASEPLATE", "HOST_RECOGNITION", "LYSIS", "DNA_REPLICATION"))
    rows=[]
    for module in names:
        candidates=[]; conflicts=[]
        for locus in loci:
            evidence = locus.get("evidence", [])
            for item in evidence:
                label = str(item.get("module", item.get("functional_category", ""))).upper()
                if module.replace("_", " ") in label or module in label:
                    role = str(item.get("role", "SUPPORT")).upper()
                    (conflicts if role in {"CONFLICT", "CONFLICTING"} else candidates).append(item)
        qualified = [x for x in candidates if str(x.get("evidence_strength", x.get("strength", ""))).upper() in {"STRONG", "HIGH", "CURATED", "EXPERIMENTAL"}]
        if conflicts and not qualified: status = "EVIDENCE_PRESENT_REVIEW_REQUIRED"
        elif len({str(x.get("source_database", x.get("source", ""))) for x in qualified}) >= 2: status = "STRONGLY_SUPPORTED"
        elif qualified: status = "SUPPORTED"
        elif candidates: status = "CANDIDATE"
        else: status = "NOT_ESTABLISHED"
        rows.append({"module": module, "status": status, "locus_ids": ";".join(sorted({str(x.get("locus_id", "")) for x in candidates if x.get("locus_id")})), "qualified_evidence_count": len(qualified), "evidence_count": len(candidates), "conflicting_evidence_count": len(conflicts), "interpretation": "NOT_ESTABLISHED is not biological absence." if status == "NOT_ESTABLISHED" else "Qualified module-level evidence; inspect locus provenance.", "provenance": candidates + conflicts})
    return rows


def assess_architecture(module_rows: Sequence[Mapping[str, Any]], *, context: Sequence[Mapping[str, Any]] = (), declared_type: str | None = None) -> dict[str, Any]:
    """Infer an architecture hypothesis from coherent module/context evidence."""
    scores = {x: 0 for x in ARCHITECTURES}
    established = {str(x.get("module")): str(x.get("status")) in {"ESTABLISHED", "STRONGLY_SUPPORTED", "SUPPORTED"} for x in module_rows}
    if established.get("TAIL") or established.get("BASEPLATE") or established.get("PORTAL"):
        scores["TAILED_DSDNA_LIKE"] += 2
    filamentous_keys = {"REP_REPLICATION", "ZOT_EXTRUSION", "COAT_VIRION", "MEMBRANE_STRUCTURAL", "INTEGRATION", "REGULATION"}
    filamentous = sum(1 for x in module_rows if str(x.get("module")) in filamentous_keys and str(x.get("status")) in {"ESTABLISHED", "STRONGLY_SUPPORTED", "SUPPORTED"})
    if filamentous >= 2:
        scores["FILAMENTOUS_PHAGE_LIKE"] += filamentous
    if established.get("DNA_REPLICATION") and not established.get("PORTAL") and not established.get("TAIL"):
        scores["SMALL_SSDNA_PHAGE_LIKE"] += 1
    best = max(scores, key=scores.get); top = scores[best]
    if top == 0: best = "ARCHITECTURE_UNRESOLVED"
    elif list(scores.values()).count(top) > 1: best = "MIXED_OR_CONFLICTING"
    return {"architecture_hypothesis": best, "scores": scores, "declared_genome_type": declared_type or "UNKNOWN", "genome_strategy": "UNKNOWN", "taxonomy": "NOT_INFERRED", "reasoning": "Architecture is a hypothesis from coherent module/context evidence; taxonomy and genome strategy remain separate.", "context_evidence": [dict(x) for x in context]}


def architecture_hallmarks(architecture: str, module_rows: Sequence[Mapping[str, Any]]) -> list[dict[str, Any]]:
    """Return architecture-specific hallmark states without forcing absence."""
    groups = {
        "FILAMENTOUS_PHAGE_LIKE": ("REP_REPLICATION", "ZOT_EXTRUSION", "COAT_VIRION", "MEMBRANE_STRUCTURAL", "INTEGRATION", "REGULATION", "GENOME_ORGANIZATION"),
        "TAILED_DSDNA_LIKE": ("DNA_PACKAGING", "PORTAL", "HEAD_CAPSID", "TAIL", "TAPE_MEASURE", "BASEPLATE", "LYSIS"),
    }
    names = groups.get(architecture, ())
    if architecture == "FILAMENTOUS_PHAGE_LIKE":
        names = names + ("DNA_PACKAGING", "PORTAL", "TAPE_MEASURE", "TAIL", "BASEPLATE")
    by={str(row.get("module")): row for row in module_rows}; rows=[]
    for name in names:
        row=by.get(name); status=str(row.get("status")) if row else "NOT_ESTABLISHED"
        if architecture == "FILAMENTOUS_PHAGE_LIKE" and name in {"DNA_PACKAGING", "PORTAL", "TAPE_MEASURE", "TAIL", "BASEPLATE"}:
            status="NOT_APPLICABLE_TO_ARCHITECTURE"
        rows.append({"architecture":architecture,"hallmark":name,"status":status,"interpretation":"NOT_ESTABLISHED is not biological absence; architecture-specific non-applicable components are not penalized.","provenance":row.get("provenance",[]) if row else []})
    return rows


def assess_comparative(*, mash_distance: float | None = None, shared_hashes: int | None = None, protein_conservation: float | None = None, gene_sharing: float | None = None, database_status: str = "READY", tool_status: str = "READY", alignment_significant: bool | None = None) -> dict[str, Any]:
    """Classify comparative significance; nearest-by-distance is never enough."""
    if database_status not in {"READY", "AVAILABLE"}: return {"state": "DATABASE_UNAVAILABLE", "significant": False, "nearest_candidate": True if mash_distance is not None else False, "interpretation": "Comparative database unavailable."}
    if tool_status not in {"READY", "AVAILABLE"}: return {"state": "TOOL_UNAVAILABLE", "significant": False, "nearest_candidate": mash_distance is not None, "interpretation": "Comparative tool unavailable."}
    if alignment_significant is False or ((shared_hashes == 0 or mash_distance == 1.0) and not (protein_conservation and protein_conservation >= 0.25)):
        return {"state": "NO_SIGNIFICANT_REFERENCE", "significant": False, "nearest_candidate": mash_distance is not None, "interpretation": "Mathematical nearest candidate lacks biologically significant support; no genus/species claim is made."}
    if alignment_significant is None and protein_conservation is None and gene_sharing is None: return {"state": "INSUFFICIENT_DATA", "significant": False, "nearest_candidate": mash_distance is not None, "interpretation": "No sufficient protein/content/synteny evidence."}
    score = max(float(protein_conservation or 0), float(gene_sharing or 0))
    state = "SIGNIFICANT_CLOSE_REFERENCE" if score >= 0.7 else ("SIGNIFICANT_DISTANT_REFERENCE" if score >= 0.25 else "WEAK_REFERENCE")
    return {"state": state, "significant": state.startswith("SIGNIFICANT"), "nearest_candidate": mash_distance is not None, "protein_conservation": protein_conservation, "gene_sharing": gene_sharing, "interpretation": "Comparative significance uses protein/content evidence; nucleotide distance alone is not a biological relative."}


def write_engine_outputs(output_dir: str | Path, *, structural: Sequence[StructuralAdjudication], functional: Sequence[FunctionalAdjudication], modules: Sequence[Mapping[str, Any]], architecture: Mapping[str, Any], comparative: Mapping[str, Any], hallmarks: Sequence[Mapping[str, Any]] | None = None) -> dict[str, str]:
    root=Path(output_dir); root.mkdir(parents=True, exist_ok=True)
    struct_rows=[asdict(x) for x in structural]; func_rows=[asdict(x) for x in functional]
    (root/"architecture_assessment.json").write_text(json.dumps(dict(architecture), indent=2, sort_keys=True)+"\n")
    (root/"comparative_assessment.json").write_text(json.dumps(dict(comparative), indent=2, sort_keys=True)+"\n")
    (root/"module_evidence.json").write_text(json.dumps(list(modules), indent=2, sort_keys=True)+"\n")
    (root/"architecture_hallmarks.json").write_text(json.dumps(list(hallmarks or []), indent=2, sort_keys=True)+"\n")
    (root/"structural_adjudication.json").write_text(json.dumps(struct_rows, indent=2, sort_keys=True)+"\n")
    (root/"functional_adjudication.json").write_text(json.dumps(func_rows, indent=2, sort_keys=True)+"\n")
    def tsv(name, rows):
        p=root/name
        flat=[]
        for row in rows:
            flat.append({k:(json.dumps(v,sort_keys=True) if isinstance(v,(dict,list,tuple)) else v) for k,v in row.items()})
        cols=list(flat[0]) if flat else ["status"]
        with p.open("w",newline="") as h:
            w=csv.DictWriter(h,fieldnames=cols,delimiter="\t"); w.writeheader(); w.writerows(flat)
    tsv("structural_adjudication.tsv",struct_rows); tsv("functional_adjudication.tsv",func_rows); tsv("module_evidence.tsv",modules); tsv("architecture_assessment.tsv",[architecture]); tsv("architecture_hallmarks.tsv",hallmarks or []); tsv("comparative_assessment.tsv",[comparative])
    manifest={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in root.iterdir() if p.is_file() and p.name not in {"engine_manifest.json"}}
    (root/"engine_manifest.json").write_text(json.dumps({"schema_version":"1.0","final_cds_policy":"PHANOTATE","secondary_callers_observational":True,"outputs":manifest},indent=2,sort_keys=True)+"\n")
    return manifest


__all__=["STRUCTURAL_CONFIDENCE","STRUCTURAL_RELATIONSHIPS","FUNCTIONAL_CLASSES","MODULE_STATUSES","ARCHITECTURES","COMPARATIVE_STATES","StructuralAdjudication","FunctionalAdjudication","adjudicate_structure","adjudicate_function","aggregate_module_evidence","assess_architecture","architecture_hallmarks","assess_comparative","write_engine_outputs"]
