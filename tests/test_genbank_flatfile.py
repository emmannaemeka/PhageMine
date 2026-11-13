import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from Bio import SeqIO
from Bio.SeqFeature import BeforePosition, AfterPosition
from phagemine.genbank import write_package
from phagemine.genbank_flatfile import write_flatfile
from phagemine.models import Protein, EvidenceLevel


class GenBankFlatfileTests(unittest.TestCase):
    def test_package_roundtrip_products_strands_overrides_and_rna(self):
        genome = 'ATGAAATAATTATTTCAT'
        proteins = [Protein('test', 'PM_000001', 1, 9, '+', genome[:9], 'MK', 'test',
                            annotation='major capsid protein', annotation_level=EvidenceLevel.COMPUTATIONAL),
                    Protein('test', 'PM_000002', 10, 18, '-', 'ATGAAATAA', 'MK', 'test',
                            partial_5prime=True, partial_3prime=True)]
        with tempfile.TemporaryDirectory() as tmp, patch('phagemine.genbank.table2asn_status', return_value={'state': 'unavailable'}):
            write_package(tmp, 'test', genome, proteins, {'input_sha256': 'fixture'})
            root = Path(tmp) / 'genbank_submission'
            record = SeqIO.read(root / 'genome.gbk', 'genbank')
            cds = [f for f in record.features if f.type == 'CDS']
            self.assertEqual(str(record.seq), genome)
            self.assertEqual(cds[0].qualifiers['product'], ['major capsid protein'])
            self.assertEqual(cds[1].qualifiers['product'], ['hypothetical protein'])
            self.assertEqual(cds[1].location.strand, -1)
            self.assertIsInstance(cds[1].location.start, BeforePosition)
            self.assertIsInstance(cds[1].location.end, AfterPosition)
            self.assertEqual(str(cds[1].extract(record.seq)), 'ATGAAATAA')
            self.assertEqual(cds[0].qualifiers['translation'], ['MK'])
            write_package(tmp, 'test', genome, proteins, {'input_sha256': 'fixture'},
                          feature_overrides={'PM_000001': {'product': 'portal protein'}})
            with (root / 'features.tbl').open('a') as handle:
                handle.write('2\t8\ttRNA\n\t\t\tproduct\ttRNA-Lys\n')
            write_flatfile(root)
            record = SeqIO.read(root / 'genome.gbk', 'genbank')
            self.assertEqual(record.features[2].qualifiers['product'], ['portal protein'])
            self.assertEqual(record.features[-1].type, 'tRNA')

    def test_mock_annotation_stays_hypothetical(self):
        from phagemine.genbank import product_name
        protein = Protein('test', 'PM_1', 1, 9, '+', 'ATGAAATAA', 'MK', 'test',
                          annotation='mock capsid protein', annotation_level=EvidenceLevel.COMPUTATIONAL)
        self.assertEqual(product_name(protein), 'hypothetical protein')
