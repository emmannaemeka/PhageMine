import json
from pathlib import Path

from phagemine.evidence_engine import (
    ARCHITECTURES, COMPARATIVE_STATES, MODULE_STATUSES,
    adjudicate_function, adjudicate_structure, aggregate_module_evidence,
    architecture_hallmarks, assess_architecture, assess_comparative,
)


def test_phanotate_is_final_and_exact_relationship_is_high():
    result=adjudicate_structure("L1", {"start":1,"end":300,"strand":"+"}, {"start":1,"end":300,"strand":"+"})
    assert result.relationship=="EXACT_MATCH"
    assert result.structural_cds_confidence=="HIGH"
    assert result.provenance["phanotate_final"] is True


def test_short_phanotate_only_can_remain_supported_with_biological_evidence():
    result=adjudicate_structure("L1", {"start":1,"end":60,"strand":"+"}, None, [{"source_database":"PHROGs","evidence_strength":"STRONG"},{"source_database":"comparative","evidence_strength":"STRONG"}])
    assert result.relationship=="PHANOTATE_ONLY"
    assert result.structural_cds_confidence=="MODERATE"


def test_family_membership_does_not_promote_specific_function():
    result=adjudicate_function("L1", "specific enzyme", [{"source_database":"VOGDB","evidence_strength":"STRONG","evidence_type":"family"}], structural_confidence="HIGH")
    assert result.functional_class=="SUPPORTED_SPECIFIC"  # qualified family is retained, but not called high confidence
    assert result.functional_assignment_confidence=="MODERATE"


def test_domain_only_and_conflict_are_preserved():
    domain=adjudicate_function("L1", "specific enzyme", [{"source_database":"Pfam","evidence_strength":"STRONG","evidence_type":"domain"}])
    assert domain.functional_class=="DOMAIN_ONLY"
    conflict=adjudicate_function("L2", "enzyme", [{"source_database":"A","status":"CONFLICTING"}])
    assert conflict.functional_class=="CONFLICTING"


def test_module_and_architecture_rules():
    loci=[{"locus_id":"r","evidence":[{"module":"REP_REPLICATION","evidence_strength":"STRONG","source_database":"Swiss-Prot"},{"module":"COAT_VIRION","evidence_strength":"STRONG","source_database":"PHROGs"}]}]
    modules=aggregate_module_evidence(loci, modules=("REP_REPLICATION","COAT_VIRION","LYSIS","PORTAL"))
    arch=assess_architecture(modules)
    assert arch["architecture_hypothesis"]=="FILAMENTOUS_PHAGE_LIKE"
    hall=architecture_hallmarks("FILAMENTOUS_PHAGE_LIKE", modules)
    assert any(x["hallmark"]=="PORTAL" and x["status"]=="NOT_APPLICABLE_TO_ARCHITECTURE" for x in hall)


def test_comparative_nearest_is_not_significant():
    result=assess_comparative(mash_distance=1.0, shared_hashes=0)
    assert result["state"]=="NO_SIGNIFICANT_REFERENCE"
    assert not result["significant"]
    assert set(COMPARATIVE_STATES)>= {result["state"]}


def test_output_schema(tmp_path):
    from phagemine.evidence_engine import write_engine_outputs
    s=adjudicate_structure("L1", {"start":1,"end":3,"strand":"+"}, None)
    f=adjudicate_function("L1", None, [], structural_confidence=s.structural_cds_confidence)
    write_engine_outputs(tmp_path, structural=[s], functional=[f], modules=[], architecture={"architecture_hypothesis":"ARCHITECTURE_UNRESOLVED"}, comparative={"state":"INSUFFICIENT_DATA"})
    assert (tmp_path/"module_evidence.tsv").exists()
    assert (tmp_path/"architecture_assessment.json").exists()
    assert (tmp_path/"engine_manifest.json").exists()
