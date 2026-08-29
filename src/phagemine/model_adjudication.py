"""Model-specific evidence and observational gene-model recommendations.

This layer never mutates the legacy PHANOTATE final CDS set.  It evaluates
candidate translations attached to a reconciled locus and records a
recommendation separately for later review.
"""
from __future__ import annotations

import csv
import hashlib
import json
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any

from .gene_models import GeneModel, ReconciledLocus
from .models import Protein


@dataclass
class CandidateModel:
    candidate_id: str
    locus_id: str
    provider_id: str
    raw_identifier: str
    start: int
    end: int
    strand: str
    cds_sequence: str
    protein_sequence: str
    protein_sha256: str
    caller_version: str | None
    input_sequence_sha256: str | None
    evidence: list[dict[str, Any]] = field(default_factory=list)
    evidence_status: str = "NOT_SEARCHED"

    @classmethod
    def from_model(cls, locus_id: str, model: GeneModel, candidate_id: str | None = None) -> "CandidateModel":
        protein = model.protein_sequence or model.sequence or ""
        cid = candidate_id or "GM_" + hashlib.sha256(f"{locus_id}|{model.caller}|{model.raw_identifier}|{model.start}|{model.end}|{model.strand}|{protein}".encode()).hexdigest()[:16]
        return cls(cid, locus_id, model.caller, model.raw_identifier, model.start, model.end,
                   model.strand, model.cds_sequence or "", protein,
                   hashlib.sha256(protein.encode()).hexdigest(), model.caller_version,
                   model.input_sequence_sha256)

    def to_dict(self):
        return asdict(self)


@dataclass
class GeneModelDecision:
    locus_id: str
    candidate_ids: list[str]
    recommended_candidate_id: str | None
    recommended_provider_id: str | None
    decision_class: str
    decision_status: str
    decision_reason: str
    caller_agreement_summary: dict[str, Any]
    model_specific_evidence_summary: dict[str, Any]
    evidence_comparability: str
    gene_model_confidence: str
    confidence_calibrated: bool = False
    review_required: bool = True
    alternative_candidates: list[str] = field(default_factory=list)

    def to_dict(self):
        return asdict(self)


def candidates_from_locus(locus: ReconciledLocus) -> list[CandidateModel]:
    return [CandidateModel.from_model(locus.locus_id, model) for model in locus.candidate_models]


def _adapter_identity(adapter) -> dict[str, Any]:
    provenance = adapter.provenance() if callable(getattr(adapter, "provenance", None)) else getattr(adapter, "provenance", {})
    return {"class": adapter.__class__.__name__, "name": getattr(adapter, "name", None), "provenance": provenance or {}, "parameters": getattr(adapter, "parameters", None)}


def evidence_cache_key(candidate: CandidateModel, adapter) -> str:
    payload = {"protein_sha256": candidate.protein_sha256, "adapter": _adapter_identity(adapter)}
    return hashlib.sha256(json.dumps(payload, sort_keys=True, default=str, separators=(",", ":")).encode()).hexdigest()


def acquire_model_specific_evidence(candidates: list[CandidateModel], adapters=(), cache_dir: str | Path | None = None) -> list[CandidateModel]:
    """Search each distinct translation once per adapter and attach evidence.

    Identical proteins share a cache/search result, while evidence is copied to
    every candidate using that translation. Cache identity includes adapter
    implementation/provenance and therefore database/version/parameter state.
    """
    cache = Path(cache_dir) if cache_dir else None
    if cache: cache.mkdir(parents=True, exist_ok=True)
    for adapter in adapters:
        groups: dict[str, list[CandidateModel]] = {}
        for candidate in candidates:
            groups.setdefault(evidence_cache_key(candidate, adapter), []).append(candidate)
        for key, group in groups.items():
            path = cache / f"{key}.json" if cache else None
            records = None
            status = "SEARCH_EXECUTED"
            if path and path.is_file():
                try:
                    payload = json.loads(path.read_text())
                    if isinstance(payload, dict) and payload.get("cache_key") == key:
                        records = payload.get("evidence", [])
                        status = "CACHE_REUSED"
                except (OSError, ValueError):
                    records = None
            if records is None:
                candidate = group[0]
                protein = Protein("candidate", candidate.candidate_id, candidate.start, candidate.end,
                                  candidate.strand, candidate.cds_sequence, candidate.protein_sequence,
                                  candidate.provider_id)
                result = adapter.analyze([protein])
                records = [asdict(item) for item in getattr(result, "evidence", [])]
                if getattr(result, "status", None) == "UNAVAILABLE":
                    status = "RESOURCE_UNAVAILABLE"
                elif not records:
                    status = "SEARCH_EXECUTED_ZERO_HITS"
                if path:
                    path.write_text(json.dumps({"cache_key": key, "evidence": records}, indent=2, sort_keys=True, default=str))
            for candidate in group:
                candidate.evidence.extend(json.loads(json.dumps(records, default=str)))
                candidate.evidence_status = status
    return candidates


def _strong(records):
    return [item for item in records if item.get("supports") and item.get("evidence_strength") in {"STRONG", "EXPERIMENTAL"}]


def _evidence_summary(candidates):
    return {candidate.candidate_id: {"records": len(candidate.evidence), "strong_records": len(_strong(candidate.evidence)), "status": candidate.evidence_status} for candidate in candidates}


def decide_locus(locus: ReconciledLocus, candidates: list[CandidateModel] | None = None) -> GeneModelDecision:
    candidates = candidates or candidates_from_locus(locus)
    ids = [candidate.candidate_id for candidate in candidates]
    providers = sorted({candidate.provider_id for candidate in candidates})
    exact = len({(candidate.start, candidate.end, candidate.strand) for candidate in candidates}) == 1
    caller_summary = {"providers": providers, "provider_count": len(providers), "candidate_count": len(candidates), "exact_coordinate_agreement": exact, "reconciliation_class": locus.reconciliation_class}
    evidence_summary = _evidence_summary(candidates)
    comparable = "COMPARABLE_EVIDENCE_AVAILABLE" if all(candidate.evidence_status in {"SEARCH_EXECUTED", "CACHE_REUSED", "SEARCH_EXECUTED_ZERO_HITS"} for candidate in candidates) and candidates else ("RESOURCE_UNAVAILABLE" if any(candidate.evidence_status == "RESOURCE_UNAVAILABLE" for candidate in candidates) else "PARTIAL_EVIDENCE_AVAILABILITY")
    recommended = None; provider = None; cls = "UNRESOLVED_GENE_MODEL"; status = "OBSERVATIONAL_RECOMMENDATION"; reason = "Evidence did not resolve the competing gene models."; confidence = "LOW"; review = True
    if exact and len(providers) > 1:
        recommended = min(candidates, key=lambda candidate: (candidate.provider_id, candidate.raw_identifier)).candidate_id
        provider = next(candidate.provider_id for candidate in candidates if candidate.candidate_id == recommended)
        cls, reason, confidence, review = "EXACT_CONSENSUS", "All active providers agree on start, end, and strand.", "HIGH", False
    elif exact:
        candidate = candidates[0]
        recommended = candidate.candidate_id if _strong(candidate.evidence) else None
        provider = candidate.provider_id if recommended else None
        cls = "CALLER_SPECIFIC_WITH_STRONG_EVIDENCE" if recommended else "CALLER_SPECIFIC_WITH_WEAK_EVIDENCE"
        reason = "A single-provider model is retained as an observational recommendation; no final CDS mutation occurs."
        confidence = "MODERATE" if recommended else "LOW"
    elif locus.strand_status == "DISCORDANT":
        cls, reason = "STRAND_CONFLICT_UNRESOLVED", "Overlapping opposite-strand models require independent biological review."
    elif locus.reconciliation_class in {"SPLIT_MODEL", "MERGED_MODEL", "COMPLEX_CONFLICT"}:
        cls, reason = "SPLIT_MERGE_UNRESOLVED", "The locus contains a split/merge or complex architecture; no majority model is selected."
    else:
        strong = {candidate.candidate_id: _strong(candidate.evidence) for candidate in candidates}
        winners = [cid for cid, records in strong.items() if records]
        signatures = {cid: {json.dumps({"source": rec.get("source"), "identifier": rec.get("identifier"), "description": rec.get("description")}, sort_keys=True, default=str) for rec in records} for cid, records in strong.items()}
        unique = [cid for cid in winners if signatures[cid] and not any(signatures[cid] & signatures[other] for other in winners if other != cid)]
        if len(unique) == 1:
            recommended = unique[0]; provider = next(candidate.provider_id for candidate in candidates if candidate.candidate_id == recommended)
            cls = "CALLER_SPECIFIC_WITH_STRONG_EVIDENCE" if len(candidates) == 1 else "MODEL_SPECIFIC_EVIDENCE_FAVORS_CANDIDATE"
            reason = "Distinct strong model-specific evidence supports this candidate; recommendation remains computational."
            confidence = "MODERATE"
        elif locus.reconciliation_class == "COMMON_STOP_ALTERNATE_START":
            cls = "CONSENSUS_ALTERNATE_START" if len(winners) > 1 else "BOUNDARY_CONFLICT_UNRESOLVED"
            reason = "Candidates share a stop but differ at the start; evidence is not uniquely discriminating." if cls.endswith("UNRESOLVED") else "Evidence supports the shared-stop alternate-start architecture without changing final CDS selection."
        elif locus.reconciliation_class == "COMMON_START_ALTERNATE_STOP":
            cls = "CONSENSUS_ALTERNATE_STOP" if len(winners) > 1 else "BOUNDARY_CONFLICT_UNRESOLVED"
            reason = "Candidates share a start but differ at the stop; evidence is not uniquely discriminating." if cls.endswith("UNRESOLVED") else "Evidence supports the shared-start alternate-stop architecture without changing final CDS selection."
        elif len(candidates) == 1:
            cls = "CALLER_SPECIFIC_WITH_WEAK_EVIDENCE" if not winners else "CALLER_SPECIFIC_WITH_STRONG_EVIDENCE"
            reason = "A caller-specific candidate is retained as an observational recommendation; no final CDS mutation occurs."
        else:
            cls = "INSUFFICIENT_COMPARATIVE_EVIDENCE"
    alternatives = [cid for cid in ids if cid != recommended]
    return GeneModelDecision(locus.locus_id, ids, recommended, provider, cls, status, reason,
                             caller_summary, evidence_summary, comparable, confidence,
                             False, review, alternatives)


def write_model_adjudication(root: str | Path, candidates: list[CandidateModel], decisions: list[GeneModelDecision]) -> None:
    root = Path(root); root.mkdir(parents=True, exist_ok=True)
    (root / "model_specific_evidence.json").write_text(json.dumps({"candidates": [candidate.to_dict() for candidate in candidates]}, indent=2, sort_keys=True, default=str))
    with (root / "model_specific_evidence.tsv").open("w", newline="") as handle:
        columns = ["candidate_id", "locus_id", "provider_id", "raw_identifier", "start", "end", "strand", "protein_sha256", "evidence_status", "evidence_count"]
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t"); writer.writeheader()
        for candidate in candidates:
            writer.writerow({"candidate_id": candidate.candidate_id, "locus_id": candidate.locus_id, "provider_id": candidate.provider_id, "raw_identifier": candidate.raw_identifier, "start": candidate.start, "end": candidate.end, "strand": candidate.strand, "protein_sha256": candidate.protein_sha256, "evidence_status": candidate.evidence_status, "evidence_count": len(candidate.evidence)})
    (root / "model_decisions.json").write_text(json.dumps({"decisions": [decision.to_dict() for decision in decisions]}, indent=2, sort_keys=True, default=str))
    with (root / "model_decisions.tsv").open("w", newline="") as handle:
        columns = ["locus_id", "candidate_ids", "reconciliation_class", "candidate_count", "provider_count", "recommended_candidate_id", "recommended_provider_id", "decision_class", "decision_status", "evidence_comparability", "gene_model_confidence", "confidence_calibrated", "review_required", "decision_reason"]
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t"); writer.writeheader()
        for decision in decisions:
            writer.writerow({"locus_id": decision.locus_id, "candidate_ids": ",".join(decision.candidate_ids), "reconciliation_class": decision.caller_agreement_summary.get("reconciliation_class", ""), "candidate_count": len(decision.candidate_ids), "provider_count": decision.caller_agreement_summary.get("provider_count", ""), "recommended_candidate_id": decision.recommended_candidate_id or "", "recommended_provider_id": decision.recommended_provider_id or "", "decision_class": decision.decision_class, "decision_status": decision.decision_status, "evidence_comparability": decision.evidence_comparability, "gene_model_confidence": decision.gene_model_confidence, "confidence_calibrated": str(decision.confidence_calibrated).lower(), "review_required": str(decision.review_required).lower(), "decision_reason": decision.decision_reason})
