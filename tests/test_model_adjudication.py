from dataclasses import dataclass

from phagemine.gene_models import GeneModel
from phagemine.model_adjudication import (
    CandidateModel,
    GeneModelDecision,
    acquire_model_specific_evidence,
    candidates_from_locus,
    decide_locus,
)
from phagemine.reconciliation_engine import reconcile_gene_models


def gm(provider, ident, start, end, strand="+", protein="M" * 20):
    return GeneModel(provider, ident, start, end, strand, sequence=protein, protein_sequence=protein,
                     cds_sequence="ATG" + "AAA" * (len(protein) - 1))


def test_exact_three_caller_is_observational_consensus():
    locus = reconcile_gene_models({"a": [gm("a", "a1", 100, 400)], "b": [gm("b", "b1", 100, 400)], "c": [gm("c", "c1", 100, 400)]})[0]
    decision = decide_locus(locus)
    assert isinstance(decision, GeneModelDecision)
    assert decision.decision_class == "EXACT_CONSENSUS"
    assert decision.confidence_calibrated is False
    assert decision.review_required is False


def test_alternate_start_without_discriminating_evidence_remains_unresolved():
    locus = reconcile_gene_models({"a": [gm("a", "a1", 100, 400)], "b": [gm("b", "b1", 120, 400)]})[0]
    decision = decide_locus(locus)
    assert decision.decision_class == "BOUNDARY_CONFLICT_UNRESOLVED"
    assert decision.review_required


def test_caller_specific_strong_evidence_is_only_a_recommendation():
    locus = reconcile_gene_models({"a": [gm("a", "a1", 100, 400)]})[0]
    candidates = candidates_from_locus(locus)
    candidates[0].evidence = [{"supports": True, "evidence_strength": "STRONG", "source": "Pfam", "identifier": "PF1", "description": "specific"}]
    candidates[0].evidence_status = "SEARCH_EXECUTED"
    decision = decide_locus(locus, candidates)
    assert decision.decision_class == "CALLER_SPECIFIC_WITH_STRONG_EVIDENCE"
    assert decision.recommended_candidate_id == candidates[0].candidate_id
    assert decision.review_required


def test_split_merge_is_not_majority_resolved():
    locus = reconcile_gene_models({"a": [gm("a", "a1", 100, 500)], "b": [gm("b", "b1", 100, 320), gm("b", "b2", 280, 500)]})[0]
    assert decide_locus(locus).decision_class == "SPLIT_MERGE_UNRESOLVED"


class FakeResult:
    status = "SUCCESS"
    evidence = []


class CountingAdapter:
    name = "fake"
    calls = 0
    def provenance(self):
        return {"database_version": "1"}
    def analyze(self, proteins):
        self.calls += len(proteins)
        return FakeResult()


def test_identical_translations_are_searched_once_and_cache_is_safe(tmp_path):
    locus = reconcile_gene_models({"a": [gm("a", "a1", 100, 400)], "b": [gm("b", "b1", 120, 420)]})[0]
    candidates = candidates_from_locus(locus)
    adapter = CountingAdapter()
    acquire_model_specific_evidence(candidates, [adapter], tmp_path)
    assert adapter.calls == 1
    assert all(c.evidence_status == "SEARCH_EXECUTED_ZERO_HITS" for c in candidates)
    adapter2 = CountingAdapter()
    acquire_model_specific_evidence(candidates, [adapter2], tmp_path)
    assert adapter2.calls == 0


def test_different_database_identity_invalidates_cache(tmp_path):
    locus = reconcile_gene_models({"a": [gm("a", "a1", 100, 400)]})[0]
    candidates = candidates_from_locus(locus)
    first = CountingAdapter(); acquire_model_specific_evidence(candidates, [first], tmp_path)
    second = CountingAdapter(); second.provenance = lambda: {"database_version": "2"}
    acquire_model_specific_evidence(candidates, [second], tmp_path)
    assert second.calls == 1
