"""Opt-in, conservative final gene-model selection.

Reconciliation and adjudication describe competing models and evidence.  This
module is the deliberately separate policy layer that decides whether a model
may replace the legacy PHANOTATE model.  The default policy is strictly
backward compatible and never changes the final CDS set.
"""
from __future__ import annotations

import csv
import json
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Iterable

from .gene_models import GeneModel, ReconciledLocus
from .model_adjudication import CandidateModel, GeneModelDecision, candidates_from_locus
from .review_priority import review_priority


PHANOTATE_ONLY = "phanotate-only"
CONSENSUS = "consensus"
CONSENSUS_POLICY_VERSION = "1.0"


@dataclass(frozen=True)
class FinalGeneModelSelection:
    locus_id: str
    selected_candidate_id: str | None
    selected_provider_id: str | None
    legacy_phanotate_candidate_id: str | None
    selection_policy: str
    selection_rule: str
    selection_reason: str
    adjudication_decision_class: str | None
    adjudication_recommended_candidate: str | None
    evidence_comparability: str
    provider_support: int
    method_family_support: int
    method_lineage_support: int
    review_required: bool
    fallback_used: bool
    changed_from_legacy: bool
    confidence_class: str
    confidence_calibrated: bool = False
    review_priority: str = "NO_REVIEW"

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


def _candidate_map(locus: ReconciledLocus) -> dict[str, CandidateModel]:
    return {candidate.candidate_id: candidate for candidate in candidates_from_locus(locus)}


def _model_map(locus: ReconciledLocus) -> dict[str, GeneModel]:
    candidates = candidates_from_locus(locus)
    return {candidate.candidate_id: model for candidate, model in zip(candidates, locus.candidate_models)}


def _is_valid_candidate(candidate: CandidateModel | None) -> bool:
    return bool(candidate and candidate.start > 0 and candidate.end >= candidate.start and candidate.strand in {"+", "-"})


def _selection_base(locus: ReconciledLocus, policy: str, decision: GeneModelDecision | None,
                    *, selected: CandidateModel | None, legacy: CandidateModel | None,
                    rule: str, reason: str, review: bool, fallback: bool,
                    confidence: str = "LOW") -> FinalGeneModelSelection:
    selected_id = selected.candidate_id if selected else None
    legacy_id = legacy.candidate_id if legacy else None
    changed = bool(selected and legacy and (selected.start, selected.end, selected.strand) != (legacy.start, legacy.end, legacy.strand))
    return FinalGeneModelSelection(
        locus_id=locus.locus_id,
        selected_candidate_id=selected_id,
        selected_provider_id=selected.provider_id if selected else None,
        legacy_phanotate_candidate_id=legacy_id,
        selection_policy=policy,
        selection_rule=rule,
        selection_reason=reason,
        adjudication_decision_class=decision.decision_class if decision else None,
        adjudication_recommended_candidate=decision.recommended_candidate_id if decision else None,
        evidence_comparability=decision.evidence_comparability if decision else "NOT_AVAILABLE",
        provider_support=locus.provider_count,
        method_family_support=locus.method_family_count,
        method_lineage_support=locus.method_lineage_count,
        review_required=review,
        fallback_used=fallback,
        changed_from_legacy=changed,
        confidence_class=confidence,
        review_priority=review_priority(locus, decision),
    )


def select_final_gene_models(
    reconciled_loci: Iterable[ReconciledLocus],
    decisions: Iterable[GeneModelDecision] | None = None,
    policy: str = PHANOTATE_ONLY,
    *,
    profile: str = "standard",
) -> list[FinalGeneModelSelection]:
    """Apply a selection policy without running reconciliation or evidence.

    ``phanotate-only`` is authoritative and exact legacy behavior.  Consensus
    is intentionally narrow: it can alter only an alternate boundary when a
    prior adjudication recommendation is explicitly discriminating and
    comparable.  Caller-specific rescue and structural conflicts always remain
    review candidates rather than silently changing annotation.
    """
    if policy not in {PHANOTATE_ONLY, CONSENSUS}:
        raise ValueError(f"unknown gene-model policy: {policy}")
    if profile not in {"standard", "extended"}:
        raise ValueError(f"unknown gene-model profile: {profile}")
    decision_map = {decision.locus_id: decision for decision in (decisions or [])}
    selections: list[FinalGeneModelSelection] = []
    for locus in sorted(reconciled_loci, key=lambda item: (item.segment_id or "", item.start, item.end, item.locus_id)):
        cmap = _candidate_map(locus)
        legacy = next((candidate for candidate in cmap.values() if candidate.provider_id.lower() == "phanotate"), None)
        decision = decision_map.get(locus.locus_id)
        if policy == PHANOTATE_ONLY:
            selections.append(_selection_base(
                locus, policy, decision, selected=legacy, legacy=legacy,
                rule="PHANOTATE_ONLY_AUTHORITATIVE",
                reason="PHANOTATE remains the authoritative legacy final CDS source.",
                review=False, fallback=False, confidence="HIGH" if legacy else "LOW"))
            continue

        # Exact agreement is safe only when all active models share one exact
        # coordinate group.  This is not a provider vote: the topology itself
        # is unambiguous, and PHANOTATE is part of the shared model when present.
        exact = len({(candidate.start, candidate.end, candidate.strand) for candidate in cmap.values()}) == 1
        if exact and len(cmap) > 1:
            selected = legacy or next(iter(cmap.values()))
            selections.append(_selection_base(
                locus, policy, decision, selected=selected, legacy=legacy,
                rule="EXACT_CONSENSUS_SELECTED",
                reason="All active candidate models have identical coordinates and strand.",
                review=False, fallback=False, confidence="HIGH"))
            continue

        cls = locus.reconciliation_class
        recommended = cmap.get(decision.recommended_candidate_id) if decision else None
        discriminating = bool(
            decision and recommended and recommended.provider_id.lower() != "phanotate"
            and decision.decision_class == "MODEL_SPECIFIC_EVIDENCE_FAVORS_CANDIDATE"
            and decision.evidence_comparability == "COMPARABLE_EVIDENCE_AVAILABLE"
            and _is_valid_candidate(recommended)
            and not (locus.strand_status == "DISCORDANT" or cls in {"SPLIT_MODEL", "MERGED_MODEL", "COMPLEX_CONFLICT"})
        )
        if cls in {"COMMON_STOP_ALTERNATE_START", "COMMON_START_ALTERNATE_STOP"} and discriminating:
            rule = "EVIDENCE_RESOLVED_ALTERNATE_START" if cls == "COMMON_STOP_ALTERNATE_START" else "EVIDENCE_RESOLVED_ALTERNATE_STOP"
            selections.append(_selection_base(
                locus, policy, decision, selected=recommended, legacy=legacy,
                rule=rule,
                reason="Comparable, discriminating model-specific evidence supports the alternative boundary.",
                review=False, fallback=False, confidence="MODERATE"))
            continue
        if cls == "CALLER_SPECIFIC" and legacy is None:
            selections.append(_selection_base(
                locus, policy, decision, selected=None, legacy=None,
                rule="RESCUE_CANDIDATE_NOT_AUTOMATICALLY_SELECTED",
                reason="A non-PHANOTATE caller-specific model is retained as a rescue candidate only.",
                review=True, fallback=False))
            continue
        if cls == "STRAND_DISCORDANCE":
            rule, reason = "PHANOTATE_FALLBACK_STRAND_CONFLICT", "Opposite-strand models are unresolved; PHANOTATE is retained pending review."
        elif cls in {"SPLIT_MODEL", "MERGED_MODEL"}:
            rule, reason = "PHANOTATE_FALLBACK_SPLIT_MERGE", "Split/merge architecture is not automatically restructured."
        elif cls == "COMPLEX_CONFLICT":
            rule, reason = "PHANOTATE_FALLBACK_COMPLEX_CONFLICT", "Complex model conflict is unresolved; PHANOTATE is retained."
        elif cls in {"BOUNDARY_DISCORDANCE", "COMMON_STOP_ALTERNATE_START", "COMMON_START_ALTERNATE_STOP"}:
            rule, reason = "PHANOTATE_FALLBACK_BOUNDARY_UNRESOLVED", "Boundary evidence is insufficiently discriminating; PHANOTATE is retained."
        elif legacy is not None:
            rule, reason = "PHANOTATE_RETAINED_CALLER_SPECIFIC", "Absence of another caller does not establish that a PHANOTATE gene is false."
        else:
            rule, reason = "RESCUE_CANDIDATE_NOT_AUTOMATICALLY_SELECTED", "No PHANOTATE model exists and non-PHANOTATE rescue is not automatic."
        selections.append(_selection_base(
            locus, policy, decision, selected=legacy, legacy=legacy, rule=rule,
            reason=reason, review=True, fallback=legacy is not None))
    return selections


def selected_models_for_loci(reconciled_loci: Iterable[ReconciledLocus], selections: Iterable[FinalGeneModelSelection]) -> list[GeneModel]:
    """Return selected models in deterministic locus order for downstream use."""
    reconciled_loci = list(reconciled_loci)
    by_locus = {locus.locus_id: _candidate_map(locus) for locus in reconciled_loci}
    models: list[GeneModel] = []
    for selection in selections:
        candidate = by_locus.get(selection.locus_id, {}).get(selection.selected_candidate_id or "")
        if candidate:
            locus = next(locus for locus in reconciled_loci if locus.locus_id == selection.locus_id)
            models.extend(model for model in locus.candidate_models if CandidateModel.from_model(locus.locus_id, model).candidate_id == candidate.candidate_id)
    return models


def write_selection_outputs(root: str | Path, loci: Iterable[ReconciledLocus], selections: Iterable[FinalGeneModelSelection],
                            decisions: Iterable[GeneModelDecision] | None = None) -> None:
    """Write additive, machine-readable selection records and a summary."""
    root = Path(root); root.mkdir(parents=True, exist_ok=True)
    loci = list(loci); selections = list(selections); decisions = list(decisions or [])
    candidate_by_locus = {locus.locus_id: _candidate_map(locus) for locus in loci}
    columns = ["final_protein_id", "locus_id", "candidate_id", "genome_id", "segment_id", "start", "end", "strand", "provider_source", "supporting_providers", "supporting_method_families", "supporting_method_lineages", "reconciliation_class", "adjudication_class", "selection_policy", "selection_rule", "selection_reason", "changed_from_legacy", "review_required", "review_priority", "gene_model_confidence", "confidence_calibrated"]
    with (root / "final_gene_models.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t"); writer.writeheader()
        for index, selection in enumerate(selections, 1):
            locus = next((item for item in loci if item.locus_id == selection.locus_id), None)
            candidate = candidate_by_locus.get(selection.locus_id, {}).get(selection.selected_candidate_id or "")
            if not locus or not candidate: continue
            writer.writerow({"final_protein_id": f"protein_{index:06d}", "locus_id": selection.locus_id, "candidate_id": candidate.candidate_id, "genome_id": locus.genome_id or "", "segment_id": locus.segment_id or "", "start": candidate.start, "end": candidate.end, "strand": candidate.strand, "provider_source": candidate.provider_id, "supporting_providers": ",".join(locus.supporting_providers), "supporting_method_families": ",".join(locus.supporting_method_families), "supporting_method_lineages": ",".join(locus.supporting_method_lineages), "reconciliation_class": locus.reconciliation_class, "adjudication_class": selection.adjudication_decision_class or "", "selection_policy": selection.selection_policy, "selection_rule": selection.selection_rule, "selection_reason": selection.selection_reason, "changed_from_legacy": str(selection.changed_from_legacy).lower(), "review_required": str(selection.review_required).lower(), "review_priority": selection.review_priority, "gene_model_confidence": selection.confidence_class, "confidence_calibrated": str(selection.confidence_calibrated).lower()})
    counts = {"total_loci": len(selections), "exact_consensus": sum(s.selection_rule == "EXACT_CONSENSUS_SELECTED" for s in selections), "unchanged_phanotate": sum(s.selected_provider_id == "phanotate" and not s.changed_from_legacy for s in selections), "alternate_start_selected": sum(s.selection_rule == "EVIDENCE_RESOLVED_ALTERNATE_START" for s in selections), "alternate_stop_selected": sum(s.selection_rule == "EVIDENCE_RESOLVED_ALTERNATE_STOP" for s in selections), "phanotate_fallback_boundary": sum(s.selection_rule == "PHANOTATE_FALLBACK_BOUNDARY_UNRESOLVED" for s in selections), "phanotate_fallback_strand": sum(s.selection_rule == "PHANOTATE_FALLBACK_STRAND_CONFLICT" for s in selections), "phanotate_fallback_split_merge": sum(s.selection_rule == "PHANOTATE_FALLBACK_SPLIT_MERGE" for s in selections), "rescue_candidates_not_selected": sum(s.selection_rule == "RESCUE_CANDIDATE_NOT_AUTOMATICALLY_SELECTED" for s in selections), "review_required": sum(s.review_required for s in selections), "changed_from_legacy": sum(s.changed_from_legacy for s in selections), "no_review": sum(s.review_priority == "NO_REVIEW" for s in selections), "low_priority_review": sum(s.review_priority == "LOW_PRIORITY_REVIEW" for s in selections), "moderate_review": sum(s.review_priority == "MODERATE_REVIEW" for s in selections), "high_priority_review": sum(s.review_priority == "HIGH_PRIORITY_REVIEW" for s in selections)}
    with (root / "gene_model_selection_summary.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t"); writer.writerow(["metric", "count"]); writer.writerows(counts.items())
    (root / "gene_model_selection.json").write_text(json.dumps({"policy_version": CONSENSUS_POLICY_VERSION, "selections": [s.to_dict() for s in selections]}, indent=2, sort_keys=True))

    # Extend the Step-1 trace without removing its legacy columns.  This is
    # intentionally separate from final CDS generation: the recommendation
    # and policy decision are provenance, not an instruction to mutate legacy
    # output files.
    trace = root / "final_gene_model_trace.tsv"
    if trace.is_file():
        with trace.open(newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
            old_fields = list(rows[0]) if rows else []
        extra = ["observational_recommended_candidate", "observational_decision_class", "selection_policy", "selection_rule", "changed_from_legacy"]
        fields = old_fields + [field for field in extra if field not in old_fields]
        by_locus = {selection.locus_id: selection for selection in selections}
        with trace.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t"); writer.writeheader()
            for row in rows:
                selection = by_locus.get(row.get("locus_id"))
                if selection:
                    row.update({"observational_recommended_candidate": selection.adjudication_recommended_candidate or "", "observational_decision_class": selection.adjudication_decision_class or "", "selection_policy": selection.selection_policy, "selection_rule": selection.selection_rule, "changed_from_legacy": str(selection.changed_from_legacy).lower()})
                writer.writerow(row)
