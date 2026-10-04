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


@pytest.mark.parametrize('description', ['Protein C', 'Protein GP45.2', 'Orf80',
                                       'gene product 31', 'putative protein ea47', 'gp12 protein'])
def test_identifiers_alone_cannot_name_protein_functions(description):
    result = classify(hit(description, gene='fixtureA'))
    assert result['proposed_function'] is None
    assert result['gene'] is None


def test_functional_roles_with_identifiers_are_preserved():
    assert normalize_function('capsid protein gp7') == 'capsid protein gp7'
    assert normalize_function('protein kinase') == 'protein kinase'
    assert normalize_function('putative DNA polymerase') == 'putative dna polymerase'


@pytest.mark.parametrize('description', [
    'Uncharacterized 7.3 kDa protein in mobB-gp55 intergenic region',
    'Uncharacterised protein dexA.1', 'Conserved hypothetical protein gp7',
])
def test_extended_unknown_descriptions_cannot_override_supported_function(description):
    unknown = hit(description, source='Swiss-Prot', reviewed=True,
                  organism='Bacteriophage T4', percent_identity=1.0,
                  query_coverage=1.0, subject_coverage=1.0, gene='y03F')
    assert classify(unknown)['proposed_function'] is None
    assert classify(unknown)['gene'] is None
    result = classify(unknown, hit('head decoration protein'))
    assert result['proposed_function'] == 'head decoration protein'
    assert result['gene'] is None
    assert unknown.description == description


def test_record_gene_symbol_alone_cannot_override_functional_description():
    identifier = hit('Protein rIIA', source='Swiss-Prot', gene_name='rIIA',
                     reviewed=True, percent_identity=100, query_coverage=1,
                     subject_coverage=1, organism='Bacteriophage T4')
    result = classify(identifier, hit('rIIA lysis inhibitor'))
    assert result['proposed_function'] == 'riia lysis inhibitor'
    assert result['gene'] is None
    assert classify(identifier)['proposed_function'] is None


def test_record_gene_symbol_does_not_remove_informative_role():
    result = classify(hit('Tail tip assembly protein I', source='Swiss-Prot',
                          gene_name='I'))
    assert result['proposed_function'] == 'tail tip assembly protein i'


def test_gene_number_wrapper_does_not_establish_function():
    assert classify(hit('gene 4.3 protein'))['proposed_function'] is None
    assert classify(hit('gene 5.7 protein'))['proposed_function'] is None
    assert normalize_function('DNA polymerase gene 5 protein') == 'dna polymerase gene 5 protein'
    assert classify(hit('lambda prophage-derived protein ea10'))['proposed_function'] is None
    assert normalize_function('prophage-derived DNA polymerase') == 'prophage-derived dna polymerase'


def test_prophage_origin_does_not_make_a_gene_symbol_informative():
    result = classify(hit('P21 prophage-derived protein fixtureA',
                          source='Swiss-Prot', gene_name='fixtureA'),
                      hit('recombination mediator protein'))
    assert result['proposed_function'] == 'recombination mediator protein'


@pytest.mark.parametrize('source,description,metrics', [
    ('Swiss-Prot', 'Protein AbcX', {'entry_name': 'ABCX_BPFIX'}),
    ('VOGDB', 'sp|FIXTURE|ABCX_BPFIX Protein AbcX', {}),
])
def test_entry_symbol_only_description_cannot_override_supported_role(source, description, metrics):
    result = classify(hit(description, source=source, **metrics),
                      hit('DNA-binding protein'))
    assert result['proposed_function'] == 'dna-binding protein'


def test_entry_symbol_matching_an_enzyme_role_remains_informative():
    result = classify(reviewed('Protein kinase', entry_name='KINASE_BPFIX'))
    assert result['proposed_function'] == 'protein kinase'


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


def reviewed(description, **metrics):
    return hit(description, source='Swiss-Prot', reviewed=True,
               organism='Synthetic bacteriophage fixture', percent_identity=1.0,
               query_coverage=1.0, subject_coverage=1.0, evalue=1e-100, **metrics)


def test_domain_rule_cannot_replace_reviewed_whole_protein_name():
    result = classify(reviewed('putative DNA polymerase', gene='fixtureA'),
                      hit('DNA polymerase family A', source='Pfam'))
    assert result['proposed_function'] == 'putative dna polymerase'
    assert result['gene'] == 'fixtureA'
    assert result['selected_by_curated_phage_anchor']


def test_conflicting_reviewed_records_do_not_gain_a_domain_rule_winner():
    result = classify(reviewed('DNA polymerase', gene='fixtureA'),
                      reviewed('DNA primase', gene='fixtureB'),
                      hit('DNA polymerase family A', source='Pfam'))
    assert result['proposed_function'] is None
    assert result['functional_state'] == 'CONFLICTING_EVIDENCE'
    assert result['gene'] is None
    assert result['review_flag'] == 'REVIEW_REQUIRED'


def full_identity(description, **metrics):
    return reviewed(description, query_length=30, subject_length=30,
                    alignment_length=30, **metrics)


def test_full_identity_record_precedes_differently_named_distant_homologue():
    exact = full_identity('probable head decoration protein', gene='fixtureA')
    distant = hit('capsid-associated protein', source='Swiss-Prot', reviewed=True,
                  organism='Synthetic bacteriophage fixture', percent_identity=65,
                  query_coverage=1, subject_coverage=1, evalue=1e-50,
                  query_length=30, subject_length=30, alignment_length=30)
    result = classify(exact, distant, hit('major capsid protein'))
    assert result['proposed_function'] == 'probable head decoration protein'
    assert result['functional_state'] == 'PROBABLE_FUNCTION'
    assert result['gene'] == 'fixtureA'
    assert len(result['product_alternatives']) == 3


def test_distinct_full_identity_records_still_require_review():
    result = classify(full_identity('DNA polymerase', gene='fixtureA'),
                      full_identity('DNA primase', gene='fixtureB'))
    assert result['proposed_function'] is None
    assert result['gene'] is None
    assert result['functional_state'] == 'CONFLICTING_EVIDENCE'


def test_perfect_partial_alignment_cannot_overrule_complete_record():
    partial = reviewed('DNA polymerase', query_length=30, subject_length=60,
                       alignment_length=30)
    result = classify(partial, full_identity('DNA primase'))
    assert result['proposed_function'] == 'dna primase'


def test_missing_lengths_cannot_grant_full_identity_precedence():
    result = classify(reviewed('DNA polymerase'), reviewed('DNA primase'))
    assert result['proposed_function'] is None


def test_unknown_full_identity_match_cannot_veto_informative_homologue():
    result = classify(full_identity('uncharacterized 7.3 kDa protein'),
                      reviewed('head decoration protein'))
    assert result['proposed_function'] == 'head decoration protein'


def test_complete_full_identity_match_disambiguates_related_capsid_members():
    relative = hit('minor capsid protein', source='Swiss-Prot', reviewed=True,
                   organism='Synthetic bacteriophage fixture', percent_identity=90,
                   query_coverage=1, subject_coverage=0.85, evalue=1e-70,
                   query_length=30, subject_length=35, alignment_length=30)
    result = classify(full_identity('major capsid protein'), relative)
    assert result['proposed_function'] == 'major capsid protein'
    assert len(result['product_alternatives']) == 2
