"""Boundary repairs exercise both strands and the serialized submission package."""
import json
import tempfile
import unittest
from dataclasses import asdict
from pathlib import Path
from unittest.mock import patch

from Bio import SeqIO
from phagemine.genbank import repair_boundaries, write_package, _reverse_complement
from phagemine.genome import translate
from phagemine.models import Protein


def fixture(sequence, left, right, strand):
    genome = sequence if strand == '+' else _reverse_complement(sequence)
    start, end = (left + 1, right) if strand == '+' else (len(sequence) - right + 1, len(sequence) - left)
    cds = sequence[left:right]
    return genome, Protein('test', 'PM_1', start, end, strand, cds, translate(cds), 'fixture')


class BoundaryRepairTests(unittest.TestCase):
    def test_stop_distance_and_codon_matrix(self):
        for strand in ('+', '-'):
            for distance in (3, 6, 9):
                for stop in ('TAA', 'TAG', 'TGA'):
                    with self.subTest(strand=strand, distance=distance, stop=stop):
                        seq = 'ATGAAA' + 'AAA' * (distance // 3 - 1) + stop + 'TAG'
                        genome, p = fixture(seq, 0, 6, strand)
                        result = repair_boundaries(genome, [p])[0]
                        self.assertEqual(result.cds, seq[:6 + distance])
                        self.assertFalse(result.partial_3prime)
                        self.assertEqual(result.sequence, 'M' + 'K' * (distance // 3))

    def test_start_distance_and_codon_matrix(self):
        for strand in ('+', '-'):
            for delta in (-3, -6, 3, 6):
                for start in ('ATG', 'GTG', 'TTG'):
                    with self.subTest(strand=strand, delta=delta, start=start):
                        seq = list('CCC' * 7 + 'TAA')
                        seq[6 + delta:9 + delta] = start
                        seq = ''.join(seq)
                        genome, p = fixture(seq, 6, len(seq), strand)
                        result = repair_boundaries(genome, [p])[0]
                        self.assertEqual(result.cds, seq[6 + delta:])
                        self.assertFalse(result.partial_5prime)
                        self.assertTrue(result.sequence.startswith('M'))

    def test_partial_roundtrip_and_provenance(self):
        for strand in ('+', '-'):
            for cds, five, three in [('ATGAAA', False, True), ('CCCTAA', True, False), ('CCCCCC', True, True)]:
                with self.subTest(strand=strand, cds=cds), tempfile.TemporaryDirectory() as tmp:
                    genome, p = fixture(cds, 0, len(cds), strand)
                    before = asdict(p)
                    with patch('phagemine.genbank.table2asn_status', return_value={'state': 'unavailable'}):
                        result = write_package(tmp, 'test', genome, [p], {'input_sha256': 'fixture'})
                    root = Path(tmp) / 'genbank_submission'
                    record = SeqIO.read(root / 'genome.gbk', 'genbank')
                    gene, feature = record.features[1:]
                    self.assertEqual(gene.location, feature.location)
                    self.assertEqual(str(feature.extract(record.seq)), cds)
                    location = ('<' if (five if strand == '+' else three) else '') + '1..' + ('>' if (three if strand == '+' else five) else '') + str(len(cds))
                    if strand == '-':
                        location = 'complement(' + location + ')'
                    self.assertIn(location, (root / 'genome.gbk').read_text())
                    model = json.loads((root / 'cds_provenance.json').read_text())['PM_1']
                    self.assertEqual(model['partial_5prime'], five)
                    self.assertEqual(model['partial_3prime'], three)
                    self.assertIn('boundary_repair', model['gene_call_parameters'])
                    self.assertEqual(asdict(p), before)
                    self.assertTrue(result['phagemine_pre_submission_validation']['valid'])

    def test_repair_is_idempotent(self):
        genome, p = fixture('ATGAAAAAATAA', 3, 9, '+')
        once = repair_boundaries(genome, [p])
        self.assertEqual(once, repair_boundaries(genome, once))

    def test_contig_edge_never_wraps_or_uses_partial_codons(self):
        for strand in ('+', '-'):
            for suffix in ('', 'T', 'TA', 'AAATA', 'AAAAAAAATAA'):
                genome, p = fixture('ATGAAA' + suffix, 0, 6, strand)
                result = repair_boundaries(genome, [p])[0]
                self.assertTrue(result.partial_3prime)
                self.assertEqual(result.cds, ('ATGAAA' + suffix) if len(suffix) < 9 else 'ATGAAA')

    def test_upstream_start_cannot_cross_stop(self):
        for strand in ('+', '-'):
            genome, p = fixture('ATGTAACCCAAATAA', 6, 15, strand)
            self.assertTrue(repair_boundaries(genome, [p])[0].partial_5prime)

    def test_invalid_models_remain_validation_errors(self):
        genome, p = fixture('ATGAAATAA', 0, 9, '+')
        for field, value in [('start', 0), ('strand', '?'), ('cds', 'ATGCCCTAA'), ('sequence', 'BAD')]:
            with self.subTest(field=field):
                from dataclasses import replace
                from phagemine.genbank import validate
                result = repair_boundaries(genome, [replace(p, **{field: value})])
                self.assertFalse(validate('test', genome, result, {'input_sha256': 'fixture'})['valid'])

    def test_repaired_translation_roundtrip(self):
        for strand in ('+', '-'):
            genome, p = fixture('GTGAAAAAATAA', 3, 9, strand)
            with tempfile.TemporaryDirectory() as tmp, patch('phagemine.genbank.table2asn_status', return_value={'state': 'unavailable'}):
                write_package(tmp, 'test', genome, [p], {'input_sha256': 'fixture'})
                root = Path(tmp) / 'genbank_submission'
                record = SeqIO.read(root / 'genome.gbk', 'genbank')
                cds = record.features[2]
                self.assertEqual(str(cds.extract(record.seq).translate(table=11, cds=True)), cds.qualifiers['translation'][0])
                self.assertEqual(str(SeqIO.read(root / 'proteins.faa', 'fasta').seq), cds.qualifiers['translation'][0])
