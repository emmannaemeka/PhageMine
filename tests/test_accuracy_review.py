import hashlib
import importlib.util
import json
from pathlib import Path

import pytest


def reviewer():
    path = Path(__file__).resolve().parents[1] / 'evaluation/v1.3_validation/scripts/complete_accuracy_review.py'
    spec = importlib.util.spec_from_file_location('accuracy_review', path)
    module = importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    return module


def test_review_rejects_changed_sequence_record_and_prediction():
    module = reviewer()
    record = {'sequence': {'value': 'MABC'}, 'primaryAccession': 'fixture'}
    case = {'product_a': 'synthetic enzyme', 'product_b': 'hypothetical protein'}
    decision = {**case, 'reviewed_accession': 'fixture', 'protein_sequence_sha256': module.digest_sequence('MABC'),
                'reviewed_record_sha256': hashlib.sha256(json.dumps(record, sort_keys=True, separators=(',', ':')).encode()).hexdigest(),
                'verdict_a': 'CONSISTENT_WITH_RECORD', 'verdict_b': 'ABSTENTION',
                'reviewer': 'Synthetic fixture', 'rationale': 'Fixture only'}
    assert module.validate_decision(case, decision, 'MABC', {'fixture': record}) == record
    with pytest.raises(ValueError, match='sequence changed'):
        module.validate_decision(case, decision, 'MABD', {'fixture': record})
    with pytest.raises(ValueError, match='prediction changed'):
        module.validate_decision({**case, 'product_a': 'other'}, decision, 'MABC', {'fixture': record})
    with pytest.raises(ValueError, match='record changed'):
        module.validate_decision(case, decision, 'MABC', {'fixture': {**record, 'version': 2}})


def test_no_exact_match_cannot_be_scored_as_supported():
    module = reviewer()
    case = {'product_a': 'synthetic enzyme', 'product_b': 'hypothetical protein'}
    decision = {**case, 'reviewed_accession': '', 'protein_sequence_sha256': module.digest_sequence('MABC'),
                'verdict_a': 'CONSISTENT_WITH_RECORD', 'verdict_b': 'ABSTENTION',
                'reviewer': 'Synthetic fixture', 'rationale': 'Fixture only'}
    with pytest.raises(ValueError, match='exact reviewed sequence'):
        module.validate_decision(case, decision, 'MABC', {})


def test_review_reads_full_system_sequences_and_rejects_duplicate_loci(tmp_path):
    module = reviewer()
    directory = tmp_path / 'NC_fixture'; directory.mkdir()
    (directory / 'proteins.faa').write_text('>PM_fixture\nMABC\n')
    (directory / 'functional_classification.tsv').write_text(
        'protein_id\tstart\tend\tstrand\nPM_fixture\t1\t12\t+\n')
    assert module.load_sequences(tmp_path) == {('NC_fixture', 1, 12, '+'): 'MABC'}
    (directory / 'proteins.faa').write_text('>PM_fixture\nMABC\n>PM_fixture\nMABC\n')
    with pytest.raises(ValueError, match='Duplicate predicted locus'):
        module.load_sequences(tmp_path)


def test_naming_equivalences_preserve_family_members_and_uncertainty():
    from phagemine.annotation_labels import classify
    from phagemine.functional_benchmark import load_synonyms, PENDING
    path = Path(__file__).resolve().parents[1] / 'evaluation/v1.3_validation/naming_equivalences.tsv'
    rules = load_synonyms(path)
    assert classify('major head protein', 'major capsid protein', rules) == 'EQUIVALENT_FUNCTION'
    for a, b in [('capsid vertex protein', 'major capsid protein'),
                 ('putative major head protein', 'major capsid protein'),
                 ('terminase large subunit', 'terminase protein'),
                 ('ATP-dependent DNA ligase', 'DNA ligase')]:
        assert classify(a, b, rules) == PENDING
