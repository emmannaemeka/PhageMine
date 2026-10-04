"""Prevent misleading names and qualifiers without using benchmark answers."""
import pytest
from phagemine.fusion import classify_protein, normalize_function
from phagemine.models import Evidence, EvidenceLevel, Protein


def hit(description, source='PHROGs', **metrics):
    return Evidence('domain' if source == 'Pfam' else 'phage_orthology', 'fixture',
                    EvidenceLevel.CURATED, source, 'fixture', status='real',
                    description=description, evidence_strength='STRONG', metrics=metrics)


def classify(*evidence):
    return classify_protein(Protein('fixture', 'P1', 1, 90, '+', 'ATG'*30,
                                   'M'*30, 'fixture', evidence=list(evidence)))


@pytest.mark.parametrize('description', ['NA', 'N/A', 'nan', '-', 'DUF1234 protein',
                                       'UPF0123 family protein', 'domain of unknown function'])
def test_unknown_labels_cannot_become_functions(description):
    result = classify(hit(description, gene='abcA', ec='1.2.3.4'))
    assert result['proposed_function'] is None
    assert result['gene'] is None and result['ec_number'] is None


def test_conflict_clears_gene_and_ec():
    result = classify(hit('integrase', gene='int', ec='1.2.3.4'),
                      hit('major capsid protein', source='Swiss-Prot'))
    assert result['functional_state'] == 'CONFLICTING_EVIDENCE'
    assert result['gene'] is None and result['ec_number'] is None


def test_domain_identifiers_do_not_name_a_whole_protein():
    result = classify(hit('major capsid protein'),
                      hit('ATPase domain', source='Pfam', gene='atpA', ec='1.2.3.4'))
    assert result['proposed_function'] == 'major capsid protein'
    assert result['gene'] is None and result['ec_number'] is None


def test_selected_product_retains_explicit_identifiers_and_raw_evidence():
    evidence = hit('DNA polymerase', gene='polA', ec='2.7.7.7')
    before = evidence.__dict__.copy()
    result = classify(evidence)
    assert result['gene'] == 'polA' and result['ec_number'] == '2.7.7.7'
    assert evidence.__dict__ == before


def test_rejected_taxon_specific_hit_does_not_lend_identifiers():
    result = classify(hit('kinetochore protein', source='Swiss-Prot', gene='abcA', ec='1.2.3.4'),
                      hit('major capsid protein'))
    assert result['proposed_function'] == 'major capsid protein'
    assert result['gene'] is None and result['ec_number'] is None


def test_family_rewrite_does_not_inherit_member_identifiers():
    result = classify(hit('Mitochondrial chaperone BCS1', source='VOGDB', gene='BCS1', ec='1.2.3.4'),
                      hit('ATPase family associated with various cellular activities (AAA)', source='Pfam'))
    assert result['proposed_function'] == 'bcs1-like aaa-family atpase'
    assert result['gene'] is None and result['ec_number'] is None


def test_informative_conserved_names_and_uncertainty_are_preserved():
    assert normalize_function('Conserved DNA helicase') == 'conserved dna helicase'
    assert normalize_function('Putative DNA helicase') == 'putative dna helicase'


def test_partial_profile_cannot_name_a_whole_protein():
    evidence = hit('DNA polymerase', query_coverage=0.95, profile_coverage=0.08,
                   bit_score=120, gene='polA', ec='2.7.7.7')
    result = classify(evidence)
    assert result['proposed_function'] is None
    assert result['gene'] is None and result['ec_number'] is None
    assert evidence.supports  # The similarity remains evidence of conservation.
    assert any('partial reference-profile' in flag for flag in result['ambiguity_flags'])


def test_partial_profile_does_not_lend_qualifiers_to_full_match():
    result = classify(hit('DNA polymerase', query_coverage=0.95, profile_coverage=0.08, gene='otherA'),
                      hit('DNA polymerase', query_coverage=0.95, profile_coverage=0.95, gene='polA'))
    assert result['gene'] == 'polA'


def test_repeated_profiles_do_not_outvote_a_stronger_distinct_function():
    broad = [hit('structural protein', query_coverage=0.9, profile_coverage=0.9,
                 percent_identity=0.25, bit_score=55, evalue=1e-9) for _ in range(8)]
    specific = hit('head decoration protein', query_coverage=1.0, profile_coverage=1.0,
                   percent_identity=0.85, bit_score=350, evalue=1e-70)
    assert classify(specific, *broad)['proposed_function'] == 'head decoration protein'


def test_duplicate_supporting_records_do_not_shrink_the_hypothesis_margin():
    best = [hit('head decoration protein', query_coverage=1.0, profile_coverage=1.0,
                percent_identity=0.85, bit_score=350, evalue=1e-70) for _ in range(2)]
    other = hit('structural protein', query_coverage=0.9, profile_coverage=0.9,
                percent_identity=0.25, bit_score=55, evalue=1e-9)
    assert classify(*best, other)['proposed_function'] == 'head decoration protein'
