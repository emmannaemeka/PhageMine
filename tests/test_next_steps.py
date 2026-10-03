from __future__ import annotations
import csv
from dataclasses import replace
import json
from pathlib import Path

import pytest


def test_database_health_versions_age_missing_and_ambiguity(tmp_path, monkeypatch):
    from phagemine.database_health import health
    from phagemine.resources import EvidenceResourceManager
    monkeypatch.setattr('phagemine.resources._tool_available', lambda tool: True)
    database = tmp_path / 'profiles.hmm'; database.write_text('synthetic')
    for suffix in ('.h3f', '.h3i', '.h3m', '.h3p'): Path(str(database)+suffix).touch()
    manager = EvidenceResourceManager(tmp_path / 'registry.json')
    manager.register('pfam', 'PFAM', database, version='test-1', provenance={'release_date': '2000-01-01'})
    result = health(manager=manager, expected_versions={'PFAM': 'test-2'})
    assert result['status'] == 'NEEDS_ATTENTION'
    assert result['resources'][0]['status'] == 'INVALID'
    assert any('AGE_WARNING' in warning for warning in result['resources'][0]['warnings'])
    assert 'PHROGS' in result['missing_annotation_resources']
    manager.register('pfam-second', 'PFAM', database, version='test-1')
    assert health(manager=manager)['ambiguous_resource_types'] == ['PFAM']
    database.unlink()
    assert all(row['status'] == 'UNAVAILABLE' for row in health(manager=manager)['resources'])


def test_persistent_cache_rebinds_ids_and_invalidates_changed_inputs(tmp_path, monkeypatch):
    from phagemine.evidence import EvidenceAdapterResult
    from phagemine.models import Evidence, EvidenceLevel, Protein
    from phagemine.pfam import PfamHMMAdapter
    from phagemine.search_cache import PersistentEvidenceCache
    database = tmp_path / 'db.hmm'; database.write_text('db-one')
    executable = tmp_path / 'hmmscan'; executable.write_text('#!/bin/sh\necho "HMMER 3.4"\n'); executable.chmod(0o755)
    calls = []
    def analyze(adapter, proteins):
        calls.append([p.protein_id for p in proteins])
        hits = [Evidence('domain', 'synthetic', EvidenceLevel.COMPUTATIONAL, 'Pfam', 'test', status='real',
                         provenance={'protein_id': p.protein_id}, metrics={'query_protein_id': p.protein_id}) for p in proteins]
        return EvidenceAdapterResult(adapter.name, 'REAL', hits, adapter.provenance())
    monkeypatch.setattr(PfamHMMAdapter, 'analyze', analyze)
    adapter = PfamHMMAdapter(database, str(executable))
    protein = Protein('g', 'first', 1, 9, '+', 'ATGAAATAG', 'MK', 'test')
    cache_dir = tmp_path / 'cache'
    PersistentEvidenceCache(cache_dir).analyze(adapter, [protein])
    second = PersistentEvidenceCache(cache_dir).analyze(adapter, [replace(protein, protein_id='second')])
    assert calls == [['first']]
    assert second.evidence[0].provenance['protein_id'] == 'second'
    assert second.evidence[0].metrics['query_protein_id'] == 'second'
    database.write_text('db-two')
    PersistentEvidenceCache(cache_dir).analyze(adapter, [protein])
    assert len(calls) == 2
    adapter.coverage_threshold = .9
    PersistentEvidenceCache(cache_dir).analyze(adapter, [protein])
    assert len(calls) == 3
    # A malformed entry is recomputed rather than trusted.
    for path in cache_dir.glob('*.json'): path.write_text('{}')
    result = PersistentEvidenceCache(cache_dir).analyze(adapter, [protein])
    assert len(calls) == 4 and result.provenance['cache']['corrupt_entries_recomputed'] == 1


def test_failed_search_is_not_cached_and_zero_hits_are_cached(tmp_path, monkeypatch):
    from phagemine.evidence import EvidenceAdapterResult
    from phagemine.models import Protein
    from phagemine.pfam import PfamHMMAdapter
    from phagemine.search_cache import PersistentEvidenceCache
    database = tmp_path / 'db.hmm'; database.write_text('synthetic')
    executable = tmp_path / 'hmmscan'; executable.write_text('#!/bin/sh\necho "HMMER 3.4"\n'); executable.chmod(0o755)
    adapter = PfamHMMAdapter(database, str(executable))
    protein = Protein('g', 'p', 1, 9, '+', 'ATGAAATAG', 'MK', 'test')
    directory = tmp_path / 'cache'
    monkeypatch.setattr(PfamHMMAdapter, 'analyze', lambda adapter, proteins: EvidenceAdapterResult(adapter.name, 'UNAVAILABLE'))
    assert PersistentEvidenceCache(directory).analyze(adapter, [protein]).status == 'UNAVAILABLE'
    assert not list(directory.glob('*.json'))
    calls = []
    def empty(adapter, proteins):
        calls.append(len(proteins)); return EvidenceAdapterResult(adapter.name, 'REAL', provenance=adapter.provenance())
    monkeypatch.setattr(PfamHMMAdapter, 'analyze', empty)
    PersistentEvidenceCache(directory).analyze(adapter, [protein])
    result = PersistentEvidenceCache(directory).analyze(adapter, [protein])
    assert calls == [1] and not result.evidence


def _table(path, rows, curated=False):
    columns = ['accession', 'start', 'end', 'strand', 'product']
    if curated: columns += ['curated_by', 'reference_evidence']
    with path.open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter='\t'); writer.writeheader(); writer.writerows(rows)


def test_curated_benchmark_preserves_unresolved_and_accounts_for_missing_models(tmp_path):
    from phagemine.curated_benchmark import evaluate
    predictions, references = tmp_path / 'predictions.tsv', tmp_path / 'reference.tsv'
    base = {'accession': 'synthetic', 'start': 1, 'end': 9, 'strand': '+'}
    _table(predictions, [{**base, 'product': 'synthetic label A'}])
    _table(references, [{**base, 'product': 'synthetic label B', 'curated_by': 'R', 'reference_evidence': 'fixture'},
                       {**base, 'start': 20, 'end': 30, 'product': 'synthetic label C', 'curated_by': 'R', 'reference_evidence': 'fixture'}], True)
    result = evaluate(predictions, references, tmp_path / 'unresolved')
    assert result['status'] == 'PENDING_ADJUDICATION'
    assert result['summary']['strict_functional_precision'] is None
    assert result['summary']['ABSTENTION'] == 1
    assert result['summary']['exact_model_recall'] == .5
    review = tmp_path / 'review.tsv'
    review.write_text('accession\tstart\tend\tstrand\tcategory\treviewer\trationale\nsynthetic\t1\t9\t+\tEQUIVALENT_FUNCTION\tR\tSynthetic synonyms\n')
    final = evaluate(predictions, references, tmp_path / 'scored', adjudications=review)
    assert final['summary']['strict_functional_precision'] == 1
    assert final['summary']['functional_recall'] == .5
    assert final['accuracy_superiority_established'] is False
    with pytest.raises(ValueError, match='new directory'): evaluate(predictions, references, tmp_path / 'scored')


def test_curated_benchmark_requires_provenance_and_exact_accessions(tmp_path):
    from phagemine.curated_benchmark import evaluate
    prediction, reference = tmp_path / 'p.tsv', tmp_path / 'r.tsv'
    base = {'start': 1, 'end': 9, 'strand': '+', 'product': 'synthetic product'}
    _table(prediction, [{**base, 'accession': 'genome_a'}])
    _table(reference, [{**base, 'accession': 'genome_b', 'curated_by': '', 'reference_evidence': ''}], True)
    with pytest.raises(ValueError, match='curator'): evaluate(prediction, reference, tmp_path / 'bad')
    _table(reference, [{**base, 'accession': 'genome_b', 'curated_by': 'R', 'reference_evidence': 'fixture'}], True)
    result = evaluate(prediction, reference, tmp_path / 'valid')
    assert result['summary']['exact_model_matches'] == 0
    assert result['summary']['unmatched_prediction_models'] == 1


def test_conflict_reports_distinguish_engine_conflicts_from_label_differences(tmp_path):
    from phagemine.conflict_review import write_conflicts
    (tmp_path / 'run_manifest.json').write_text('{}')
    (tmp_path / 'annotation.tsv').write_text('protein_id\tproduct\np1\tsynthetic product A\np2\thypothetical protein\n')
    (tmp_path / 'functional_classification.json').write_text(json.dumps([{'protein_id': 'p1', 'functional_state': 'CONFLICTING_EVIDENCE', 'conflicting_evidence_ids': ['id1'], 'conflicting_sources': ['S']}]))
    (tmp_path / 'evidence.json').write_text(json.dumps([{'protein_id': 'p1', 'evidence': []}, {'protein_id': 'p2', 'evidence': []}]))
    (tmp_path / 'external_evidence.json').write_text(json.dumps({'records': [{'protein_id': pid, 'proposed_product': 'synthetic product B', 'source': 'S', 'external_id': 'id2'} for pid in ('p1', 'p2')]}))
    report = write_conflicts(tmp_path)
    assert [r['kind'] for r in report['records']] == ['ENGINE_EVIDENCE_CONFLICT', 'EXTERNAL_LABEL_DIFFERENCE']
    assert report['automatic_product_changes'] is False
    assert 'does not establish biological incompatibility' in report['records'][1]['reason']
