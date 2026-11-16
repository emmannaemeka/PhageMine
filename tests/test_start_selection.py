import csv
from dataclasses import replace
from unittest.mock import patch

import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from phagemine.models import Protein
from phagemine.gene_models import GeneModel
from phagemine.genome import translate
from phagemine.start_selection import select_supported_starts, STARTS_11
from phagemine.genbank import write_package


def fixtures(sequence, old, new, strand):
    genome = sequence if strand == '+' else str(Seq(sequence).reverse_complement())
    def coords(start):
        return (start + 1, len(sequence)) if strand == '+' else (1, len(sequence) - start)
    a, b = coords(old)
    p = Protein('synthetic', 'PM_1', a, b, strand, sequence[old:], translate(sequence[old:]), 'fixture')
    a, b = coords(new)
    return genome, p, GeneModel('prodigal_gv', 'GV1', a, b, strand)


@pytest.mark.parametrize('strand', ['+', '-'])
@pytest.mark.parametrize('start', sorted(STARTS_11))
def test_invalid_start_replaced_by_same_stop_candidate(strand, start):
    genome, p, c = fixtures('CCC' + start + 'AAA' * 10 + 'TAA', 0, 3, strand)
    models, changes = select_supported_starts(genome, [p], [c])
    assert models[0].cds[:3] == start and models[0].sequence.startswith('M')
    assert changes[0]['distance_nt'] == 3
    assert (models[0].end if strand == '+' else models[0].start) == (p.end if strand == '+' else p.start)
    assert p.cds.startswith('CCC')


@pytest.mark.parametrize('strand', ['+', '-'])
def test_upstream_sd_candidate_can_replace_unsupported_start(strand):
    seq = 'AGGAGG' + 'CCC' * 2 + 'GTG' + 'CCC' * 10 + 'ATG' + 'CCC' * 10 + 'TAA'
    genome, p, c = fixtures(seq, 45, 12, strand)
    models, changes = select_supported_starts(genome, [p], [c])
    assert models[0].cds.startswith('GTG')
    assert changes[0]['rule'] == 'SAME_STRAND_SAME_STOP_PRODIGAL_GV_START'


def test_prodigal_preferred_even_when_both_have_sd_motifs():
    seq = 'AGGAGGCCCCCCGTG' + 'CCC' * 4 + 'AGGAGGCCCCCCATG' + 'CCC' * 10 + 'TAA'
    genome, p, c = fixtures(seq, 39, 12, '+')
    models, changes = select_supported_starts(genome, [p], [c])
    assert changes[0]['distance_nt'] == -27
    assert models[0].cds.startswith('GTG')


def test_no_stop_crossing_or_frame_changes_or_new_stop():
    genome, p, c = fixtures('ATGTAACCC' + 'AAA' * 10 + 'TAA', 6, 0, '+')
    for candidate in (c, replace(c, end=c.end-3), replace(c, start=2)):
        models, changes = select_supported_starts(genome, [p], [candidate])
        assert not changes
        assert (models[0].start, models[0].end, models[0].cds) == (p.start, p.end, p.cds)


def test_explicit_partial_preserved():
    genome, p, c = fixtures('CCCATGAAATAA', 0, 3, '+')
    p.partial_5prime = True
    assert select_supported_starts(genome, [p], [c]) == ([p], [])


@pytest.mark.parametrize('strand', ['+', '-'])
def test_orf_table_roundtrip_uses_export_models(tmp_path, strand):
    genome, p, c = fixtures('CCCGTGAAAAAATAA', 0, 3, strand)
    models, _ = select_supported_starts(genome, [p], [c])
    with patch('phagemine.genbank.table2asn_status', return_value={'state':'unavailable'}):
        write_package(tmp_path, 'synthetic', genome, models, {'input_sha256':'fixture'})
    root = tmp_path/'genbank_submission'
    rows = list(csv.DictReader((root/'orfs.tsv').open(), delimiter='\t'))
    record = SeqIO.read(root/'genome.gbk', 'genbank')
    features = [f for f in record.features if f.type == 'CDS']
    assert len(rows) == len(features) == 1
    for row, feature in zip(rows, features):
        assert int(row['start']) == int(feature.location.start)+1
        assert int(row['end']) == int(feature.location.end)
        assert row['strand'] == ('+' if feature.location.strand == 1 else '-')
        assert row['nucleotide_sequence'] == str(feature.extract(record.seq))
        assert row['protein_sequence'] == feature.qualifiers['translation'][0]
        assert row['start_codon'] == row['nucleotide_sequence'][:3]
        assert row['stop_codon'] == row['nucleotide_sequence'][-3:]
    assert '\t\t\ttransl_table\t11\n' in (root/'features.tbl').read_text()


@pytest.mark.parametrize('strand', ['+', '-'])
def test_agreed_start_unchanged_and_unmatched_retained(strand):
    genome, p, c = fixtures('ATGAAAAAATAA', 0, 0, strand)
    assert select_supported_starts(genome, [p], [c]) == ([p], [])
    models, changes = select_supported_starts(genome, [p], [])
    assert not changes and models[0].cds == p.cds
    assert 'review' in models[0].gene_call_parameters['start_selection_review']


@pytest.mark.parametrize('start', ['ATT','ATC','ATA','CTG'])
def test_rare_starts_pass_export_and_are_flagged(tmp_path, start):
    from phagemine.genbank import validate
    from importlib.util import spec_from_file_location, module_from_spec
    from pathlib import Path
    genome, p, _ = fixtures(start+'AAAAAATAA', 0, 0, '+')
    with patch('phagemine.genbank.table2asn_status', return_value={'state':'unavailable'}):
        result = write_package(tmp_path, 'synthetic', genome, [p], {'input_sha256':'fixture'})
    root = tmp_path/'genbank_submission'
    row = next(csv.DictReader((root/'orfs.tsv').open(), delimiter='\t'))
    record = SeqIO.read(root/'genome.gbk', 'genbank')
    feature = next(f for f in record.features if f.type=='CDS')
    assert row['start_codon'] == start and row['partial_5prime']=='False'
    assert 'rare start codon - review' in row['review_flag']
    assert row['protein_sequence'].startswith('M')
    assert not any(x['code']=='noncanonical_start_codon' for x in result['phagemine_pre_submission_validation']['warnings'])
    spec=spec_from_file_location('verify_codons_test', Path(__file__).parents[1]/'scripts/verify_codons.py')
    verifier=module_from_spec(spec);spec.loader.exec_module(verifier)
    assert not verifier.check_cds(record, feature, 11, extended_starts=True)['problems']
    assert verifier.check_cds(record, feature, 11, extended_starts=False)['problems']
