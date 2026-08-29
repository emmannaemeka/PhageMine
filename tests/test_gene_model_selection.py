from phagemine.gene_models import GeneModel
from phagemine.reconciliation_engine import reconcile_gene_models
from phagemine.model_adjudication import CandidateModel, GeneModelDecision
from phagemine.gene_model_selection import select_final_gene_models, PHANOTATE_ONLY, CONSENSUS


def gm(caller, start, end, strand="+", protein="M" * 20, family=None, lineage=None):
    return GeneModel(caller, f"{caller}_{start}_{end}", start, end, strand,
                     sequence=protein, protein_sequence=protein, cds_sequence="ATG" + "A" * (len(protein) * 3 - 3),
                     genome_id="synthetic", method_family=family or caller, method_lineage=lineage or family or caller)


def locus(models):
    return reconcile_gene_models({m.caller: [m] for m in models})[0]


def test_phanotate_only_is_authoritative():
    l = locus([gm("phanotate", 100, 300), gm("pyrodigal", 120, 300, family="prodigal", lineage="prodigal")])
    s = select_final_gene_models([l], policy=PHANOTATE_ONLY)[0]
    assert s.selected_provider_id == "phanotate"
    assert not s.changed_from_legacy
    assert s.selection_rule == "PHANOTATE_ONLY_AUTHORITATIVE"


def test_exact_consensus_does_not_change_legacy():
    l = locus([gm("phanotate", 100, 300), gm("pyrodigal", 100, 300, family="prodigal", lineage="prodigal")])
    s = select_final_gene_models([l], policy=CONSENSUS)[0]
    assert s.selection_rule == "EXACT_CONSENSUS_SELECTED"
    assert s.selected_provider_id == "phanotate"
    assert not s.changed_from_legacy


def test_alternate_start_without_discriminating_evidence_falls_back():
    p = gm("phanotate", 100, 300)
    q = gm("pyrodigal", 130, 300, family="prodigal", lineage="prodigal")
    l = locus([p, q]); c = CandidateModel.from_model(l.locus_id, q)
    d = GeneModelDecision(l.locus_id, [CandidateModel.from_model(l.locus_id, p).candidate_id, c.candidate_id], None, None,
                          "BOUNDARY_CONFLICT_UNRESOLVED", "OBSERVATIONAL_RECOMMENDATION", "insufficient", {}, {},
                          "COMPARABLE_EVIDENCE_AVAILABLE", "LOW", review_required=True)
    s = select_final_gene_models([l], [d], policy=CONSENSUS)[0]
    assert s.selected_provider_id == "phanotate"
    assert not s.changed_from_legacy
    assert s.review_required


def test_discriminating_alternate_start_may_change_boundary():
    p = gm("phanotate", 100, 300)
    q = gm("pyrodigal", 130, 300, family="prodigal", lineage="prodigal")
    l = locus([p, q]); cp = CandidateModel.from_model(l.locus_id, p); cq = CandidateModel.from_model(l.locus_id, q)
    d = GeneModelDecision(l.locus_id, [cp.candidate_id, cq.candidate_id], cq.candidate_id, "pyrodigal",
                          "MODEL_SPECIFIC_EVIDENCE_FAVORS_CANDIDATE", "OBSERVATIONAL_RECOMMENDATION", "discriminating", {}, {},
                          "COMPARABLE_EVIDENCE_AVAILABLE", "MODERATE", review_required=False)
    s = select_final_gene_models([l], [d], policy=CONSENSUS)[0]
    assert s.selected_provider_id == "pyrodigal"
    assert s.changed_from_legacy
    assert s.selection_rule == "EVIDENCE_RESOLVED_ALTERNATE_START"


def test_nonphanotate_caller_specific_is_rescue_only():
    l = locus([gm("pyrodigal", 100, 300, family="prodigal", lineage="prodigal")])
    s = select_final_gene_models([l], policy=CONSENSUS)[0]
    assert s.selected_candidate_id is None
    assert s.selection_rule == "RESCUE_CANDIDATE_NOT_AUTOMATICALLY_SELECTED"


def test_strand_conflict_falls_back_to_phanotate():
    l = locus([gm("phanotate", 100, 300, "+"), gm("pyrodigal", 100, 300, "-", family="prodigal", lineage="prodigal")])
    s = select_final_gene_models([l], policy=CONSENSUS)[0]
    assert s.selected_provider_id == "phanotate"
    assert s.selection_rule == "PHANOTATE_FALLBACK_STRAND_CONFLICT"
    assert s.review_required


def test_lineage_support_is_metadata_not_voting():
    p = gm("phanotate", 100, 300)
    q = gm("pyrodigal", 130, 300, family="prodigal", lineage="prodigal")
    r = gm("prodigal", 130, 300, family="prodigal", lineage="prodigal")
    l = locus([p, q, r]); s = select_final_gene_models([l], policy=CONSENSUS)[0]
    assert s.provider_support == 3
    assert s.method_lineage_support == 2
    assert s.selected_provider_id == "phanotate"


def test_selection_is_order_invariant():
    models = [gm("phanotate", 100, 300), gm("pyrodigal", 100, 300, family="prodigal", lineage="prodigal")]
    a = select_final_gene_models([locus(models)], policy=CONSENSUS)[0]
    b = select_final_gene_models([locus(list(reversed(models)))], policy=CONSENSUS)[0]
    assert a.to_dict() == b.to_dict()
