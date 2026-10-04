import importlib.util
from pathlib import Path
import sys
import csv
import json
import hashlib
import pytest

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation
from Bio.SeqRecord import SeqRecord


def runner():
    directory = Path(__file__).resolve().parents[1] / 'evaluation/v1.3_validation/scripts'
    spec = importlib.util.spec_from_file_location('fresh_runner', directory / 'run_fresh.py')
    module = importlib.util.module_from_spec(spec)
    sys.path.insert(0, str(directory))
    try: spec.loader.exec_module(module)
    finally: sys.path.remove(str(directory))
    return module


def test_fresh_inputs_do_not_leak_reference_annotations(tmp_path):
    source = tmp_path / 'reference'; source.mkdir()
    output = tmp_path / 'output'; output.mkdir()
    record = SeqRecord(Seq('ATG' * 30), id='SYNTHETIC.1', name='SYNTHETIC',
                       annotations={'molecule_type': 'DNA'})
    record.features = [SeqFeature(FeatureLocation(0, 30, strand=1), type='CDS',
                                  qualifiers={'product': ['reference answer marker'], 'gene': ['answerA']}),
                       SeqFeature(CompoundLocation([FeatureLocation(35, 40, strand=1),
                                                     FeatureLocation(45, 55, strand=1)]), type='CDS')]
    SeqIO.write(record, source / 'fixture.gb', 'genbank')
    genomes, references, excluded = runner().prepare(source, output)
    fasta = Path(genomes[0]['fasta']).read_text()
    assert 'answer' not in fasta
    assert fasta == '>SYNTHETIC.1\n' + str(record.seq) + '\n'
    assert references[0]['gene'] == 'answerA'
    assert references[0]['start'] == 1 and references[0]['end'] == 30
    assert len(references) == 1 and excluded == 1
    other = tmp_path / 'filtered'; other.mkdir()
    assert runner().prepare(source, other, 'SYNTHETIC.1')[0][0]['accession'] == 'SYNTHETIC.1'


def audit_runner():
    directory = Path(__file__).resolve().parents[1] / 'evaluation/v1.3_validation/scripts'
    spec = importlib.util.spec_from_file_location('fresh_audit', directory / 'evaluate_fresh_results.py')
    module = importlib.util.module_from_spec(spec)
    sys.path.insert(0, str(directory))
    try: spec.loader.exec_module(module)
    finally: sys.path.remove(str(directory))
    return module


def fresh_artifact(root, accession):
    root.mkdir(); (root / 'inputs').mkdir()
    fasta = root / 'inputs' / (accession + '.fasta'); fasta.write_text('>fixture\nATG\n')
    digest = hashlib.sha256(fasta.read_bytes()).hexdigest()
    manifest = {'fresh_annotation_run': True, 'analysis_type': 'FRESH_REFERENCE_CONCORDANCE',
                'source_commit': 'a' * 40,
                'tool_status': {tool: 'FRESH_RUN_COMPLETED' for tool in ('PhageMine', 'Pharokka', 'Prokka')},
                'genomes': [{'accession': accession, 'fasta_sha256': digest}]}
    (root / 'manifest.json').write_text(json.dumps(manifest))
    (root / 'database_checksums.tsv').write_text('fixture database fingerprint\n')
    (root / 'SHA256SUMS.tsv').write_text(f'path\tsha256\ninputs/{accession}.fasta\t{digest}\n')
    return manifest


def test_fresh_audit_rejects_archived_data_changed_inputs_and_mixed_databases(tmp_path):
    module = audit_runner()
    first, second = tmp_path / 'first', tmp_path / 'second'
    manifest = fresh_artifact(first, 'synthetic1')
    fresh_artifact(second, 'synthetic2')
    assert len(module.validate_sources([first, second])) == 2
    (second / 'database_checksums.tsv').write_text('different database\n')
    with pytest.raises(ValueError, match='different databases'): module.validate_sources([first, second])
    (first / 'inputs/synthetic1.fasta').write_text('>fixture\nCCC\n')
    with pytest.raises(ValueError, match='checksum mismatch'): module.validate_sources([first])
    manifest['fresh_annotation_run'] = False
    (first / 'manifest.json').write_text(json.dumps(manifest))
    with pytest.raises(ValueError, match='completed fresh'): module.validate_sources([first])
