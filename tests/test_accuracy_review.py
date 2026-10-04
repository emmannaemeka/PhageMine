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
