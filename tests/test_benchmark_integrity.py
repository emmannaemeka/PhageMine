from __future__ import annotations
import csv
import json
from pathlib import Path

import pytest


def table(path, rows, reference=False):
    columns = ['accession', 'start', 'end', 'strand', 'product']
    if reference: columns += ['curated_by', 'reference_evidence']
    with path.open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter='\t'); writer.writeheader(); writer.writerows(rows)


def test_genbank_version_and_wrapped_products_on_separate_records(tmp_path):
    from phagemine.benchmark import import_genbank
    path = tmp_path / 'multi.gbk'
    path.write_text('LOCUS       TEST_A\nACCESSION   TEST_A\nVERSION     TEST_A.2\n     CDS             10..21\n                     /product="synthetic long\n                     product"\n                     /translation="MPEP\n                     TIDE"\n//\nLOCUS       TEST_B\nACCESSION   TEST_B\nVERSION     TEST_B.4\n     CDS             complement(10..21)\n                     /product="another synthetic product"\n//\n')
    rows = import_genbank(path, 'fixture')
    assert [(r['genome_id'], r['product']) for r in rows] == [('TEST_A.2', 'synthetic long product'), ('TEST_B.4', 'another synthetic product')]
    assert rows[0]['protein_length'] == 8
    assert [r['genome_id'] for r in import_genbank(path, 'fixture', 'explicit')] == ['explicit', 'explicit']


def test_stale_review_decisions_and_synonym_changes_are_rejected(tmp_path):
    from phagemine.curated_benchmark import evaluate
    base = {'accession': 'held_fixture', 'start': 10, 'end': 21, 'strand': '+'}
    prediction, reference = tmp_path / 'prediction.tsv', tmp_path / 'reference.tsv'
    table(prediction, [{**base, 'product': 'synthetic A'}])
    table(reference, [{**base, 'product': 'synthetic B', 'curated_by': 'Fixture reviewer', 'reference_evidence': 'fixture only'}], True)
    out = tmp_path / 'first'; evaluate(prediction, reference, out)
    with (out / 'review_queue.tsv').open() as handle:
        reader = csv.DictReader(handle, delimiter='\t'); fields = reader.fieldnames; decisions = list(reader)
    decisions[0].update(category='EQUIVALENT_FUNCTION', reviewer='Fixture reviewer', rationale='Explicit synthetic decision')
    review = tmp_path / 'review.tsv'
    with review.open('w', newline='') as handle:
        writer=csv.DictWriter(handle, fieldnames=fields, delimiter='\t'); writer.writeheader(); writer.writerows(decisions)
    assert evaluate(prediction, reference, tmp_path / 'valid', adjudications=review)['status'] == 'SCORED'
    table(prediction, [{**base, 'product': 'changed synthetic assertion'}])
    with pytest.raises(ValueError, match='source changed'): evaluate(prediction, reference, tmp_path / 'stale', adjudications=review)
    assert not (tmp_path / 'stale').exists()
    table(prediction, [{**base, 'product': 'synthetic A'}])
    synonyms = tmp_path / 'synonyms.tsv'; synonyms.write_text('term_a\tterm_b\trule_id\trationale\nx\ty\tfixture\tSynthetic rule\n')
    with pytest.raises(ValueError, match='source changed'): evaluate(prediction, reference, tmp_path / 'changed-rules', adjudications=review, synonyms=synonyms)


def test_compound_reference_locations_are_not_flattened(tmp_path):
    from phagemine.curated_benchmark import read_records
    path = tmp_path / 'reference.tsv'
    path.write_text('accession\tstart\tend\tstrand\tproduct\tcurated_by\treference_evidence\tcompound_location\ng\t1\t40\t+\tsynthetic\tR\tfixture\tTrue\n')
    with pytest.raises(ValueError, match='Compound CDS'): read_records(path, curated=True)


def test_error_analysis_respects_negative_strand_and_accession(tmp_path):
    import importlib.util
    path = Path(__file__).resolve().parents[1] / 'evaluation/v1.3_validation/scripts/run_archived.py'
    spec = importlib.util.spec_from_file_location('archived_validation', path)
    module = importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    ref = {'accession': 'new_fixture', 'start': 10, 'end': 40, 'strand': '-', 'product': 'synthetic'}
    pred = {**ref, 'end': 43}
    assert module.error_cases('fixture', [pred], [ref])[0]['kind'] == 'START_BOUNDARY_DIFFERENCE'
    assert module.error_cases('fixture', [{**pred, 'accession': 'other_fixture'}], [ref])[0]['kind'] == 'NO_OVERLAPPING_MODEL'


def test_archived_reproduction_preserves_historical_status(tmp_path):
    import subprocess, sys
    root=Path(__file__).resolve().parents[1]
    path=root/'evaluation/v1.3_validation/scripts/run_archived.py'
    result=subprocess.run([sys.executable, str(path), '--output', str(tmp_path/'output')], capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    report=json.loads((tmp_path/'output/manifest.json').read_text())
    assert report['fresh_annotation_run'] is False
    assert report['held_out'] is False
    assert report['tool_status']['Phold'] == 'NOT_RUN_NO_OUTPUTS_AVAILABLE'
    assert report['excluded_compound_reference_cds'] == 11
    for tool in ('PhageMine', 'Pharokka', 'Prokka'):
        report=json.loads((tmp_path/'output'/tool/'curated_benchmark.json').read_text())
        assert report['summary']['strict_functional_precision'] is None
        assert report['accuracy_superiority_established'] is False


def _publisher():
    import importlib.util
    path=Path(__file__).resolve().parents[1] / 'tools/publish_verified_release.py'
    spec=importlib.util.spec_from_file_location('verified_release', path)
    module=importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    return module


def test_release_requires_all_latest_main_checks_on_exact_commit():
    module=_publisher(); sha='a'*40
    runs=[{'id': i, 'name': name, 'head_branch': 'main', 'head_sha': sha, 'event': 'push', 'status': 'completed', 'conclusion': 'success'} for i,name in enumerate(module.REQUIRED,1)]
    def get(path):
        return {'object': {'sha': sha}} if '/git/ref/' in path else {'workflow_runs': runs}
    assert module.ready('owner/repo', sha, get)
    runs.append({**runs[0], 'id': 10, 'conclusion': 'failure'})
    assert not module.ready('owner/repo', sha, get)
    runs.pop(); runs[0]['event']='pull_request'
    assert not module.ready('owner/repo', sha, get)
    runs[0]['event']='push'; runs[0]['head_branch']='fork-branch'
    assert not module.ready('owner/repo', sha, get)
    assert not module.ready('owner/repo', 'b'*40, get)


def test_release_rejects_existing_version_tag_on_another_commit():
    module=_publisher(); sha='a'*40
    def existing(path, optional=False): return {'object': {'type': 'commit', 'sha': 'b'*40}}
    with pytest.raises(ValueError, match='different commit'): module.tag_matches('owner/repo', 'v1.3.0', sha, existing)
    module.tag_matches('owner/repo', 'v1.3.0', sha, lambda path, optional=False: None)
