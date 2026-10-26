import copy
import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation
from Bio import SeqIO

from phagemine.annotation_reliability import (coding_capacity, caller_review, parse_alignments,
    reference_gene_review, reference_screen, write_reliability_review)
from phagemine.models import Protein
from phagemine.quality import assess


def protein(start=1, end=90, strand='+'):
    return Protein('g', 'p', start, end, strand, 'ATG' + 'AAA' * ((end-start+1)//3-2) + 'TAA', 'M', 'fixture')


class ReliabilityTests(unittest.TestCase):
    def test_union_density_and_empty(self):
        self.assertEqual(coding_capacity(100, [(1, 60), (40, 90), (1, 60)])['coding_density_percent'], 90)
        self.assertIsNone(coding_capacity(0, [])['coding_density_percent'])
        self.assertEqual(coding_capacity(100, [(0, 10), (80, 110)])['coding_bases'], 0)

    def test_low_density_warning_is_not_code_diagnosis(self):
        result = assess('A' * 100, [])
        warning = next(w for w in result['warnings'] if w['code'] == 'low_coding_density')
        self.assertIn('not proof', warning['message'])
        self.assertEqual(result['coding_density_percent'], 0)

    def test_native_code_and_short_candidates_preserve_models(self):
        p = protein(); original = copy.deepcopy(p)
        records = [{'begin': 1, 'end': 90, 'strand': 1, 'translation_table': 15, 'raw_identifier': 'gv1'},
                   {'begin': 100, 'end': 189, 'strand': -1, 'translation_table': 15, 'raw_identifier': 'gv2'}]
        result = caller_review('A'*200, [p], {'records': records})
        self.assertEqual(result['alternative_translation_tables'], [15])
        self.assertEqual(result['short_cds'][0]['assessment'], 'EXACT_SECONDARY_SUPPORT')
        self.assertEqual(result['short_cds'][1]['assessment'], 'POSSIBLE_MISSED_SHORT_CDS')
        self.assertEqual(p, original)
        self.assertEqual(caller_review('A'*200, [p], {})['alternative_code_status'], 'NOT_ASSESSED')

    def test_hsp_overlap_reverse_orientation_and_projection(self):
        text = '1\t9\t9\t1\t50\tATGAAATAA\tATGAAATAA\n1\t9\t9\t1\t20\tATGAAATAA\tAAAAAAAAA\n'
        coverage, identity, hsps = parse_alignments(text, 10, 9)
        self.assertEqual((coverage, identity), (90, 100))
        record = SeqRecord(Seq('TTATTTCAT'), id='reference')
        record.features = [SeqFeature(FeatureLocation(0, 9, strand=-1), type='CDS')]
        rows = reference_gene_review(record, hsps, [protein(1, 9)])
        self.assertEqual(rows[0]['status'], 'EXACT_BOUNDARY_SUPPORT')
        self.assertEqual(reference_gene_review(record, hsps, [])[0]['status'], 'POSSIBLE_MISSED_CDS')
        self.assertEqual(reference_gene_review(record, hsps, [protein(1, 6)])[0]['status'], 'BOUNDARY_DISAGREEMENT')

    def test_gaps_compound_and_malformed_alignments(self):
        _, _, hsps = parse_alignments('1\t10\t1\t9\t50\tATGAAAATAA\tATG-AAATAA', 10, 9)
        record = SeqRecord(Seq('ATGAAATAA'), id='r')
        record.features = [SeqFeature(FeatureLocation(0, 9, strand=1), type='CDS')]
        self.assertEqual(reference_gene_review(record, hsps, [])[0]['status'], 'INDEL_REQUIRES_REVIEW')
        record.features[0].location = CompoundLocation([FeatureLocation(0,3,strand=1), FeatureLocation(6,9,strand=1)])
        self.assertEqual(reference_gene_review(record, hsps, [])[0]['status'], 'UNMAPPABLE_REFERENCE_CDS')
        with self.assertRaises(ValueError):
            parse_alignments('1\t9\t1\t9\t50\tAAA\tAAA', 9, 9)

    def test_reference_screen_and_no_tool(self):
        from types import SimpleNamespace
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp)/'refs.gb'
            record = SeqRecord(Seq('ATGAAATAA'), id='NC_000001.1', annotations={'molecule_type':'DNA'})
            record.features = [SeqFeature(FeatureLocation(0, 9, strand=1), type='CDS')]
            SeqIO.write([record], path, 'genbank')
            with patch('phagemine.annotation_reliability.shutil.which', return_value='/fake/blastn'), patch('phagemine.annotation_reliability.subprocess.run', return_value=SimpleNamespace(returncode=0, stdout='1\t9\t1\t9\t50\tATGAAATAA\tATGAAATAA', stderr='')):
                result = reference_screen('ATGAAATAA', [protein(1,9)], [path], 'blastn')
                self.assertTrue(result['candidates'][0]['selected'])
                self.assertEqual(result['candidates'][0]['ictv_reference_status'], 'NOT_VERIFIED')
                self.assertEqual(result['gene_reviews'][0]['status'], 'EXACT_BOUNDARY_SUPPORT')
            with patch('phagemine.annotation_reliability.shutil.which', return_value=None):
                self.assertEqual(reference_screen('AAA', [], [path], 'missing')['status'], 'UNAVAILABLE')

    def test_report_does_not_consume_stale_secondary_after_failure(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp); raw=root/'gene_calls'/'raw'; raw.mkdir(parents=True)
            (raw/'prodigal_gv.json').write_text(json.dumps({'records':[{'translation_table':15}]}))
            (root/'gene_calls'/'gene_call_manifest.json').write_text(json.dumps({'providers':[{'provider_id':'prodigal_gv','status':'FAILED_PROVIDER'}]}))
            result=write_reliability_review(root, 'A'*100, [])
            self.assertEqual(result['alternative_code_status'], 'NOT_ASSESSED')
            self.assertTrue((root/result['report']).exists())
            self.assertEqual(result['reference_status'], 'NOT_REQUESTED')

    def test_pipeline_exposes_review_without_changing_gene_predictions(self):
        from phagemine.pipeline import run
        from phagemine.gene_prediction import DemoORFPredictor
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp); fasta=root/'g.fa'; out=root/'result'
            fasta.write_text('>g\nATG' + 'AAA'*40 + 'TAA' + 'C'*80 + '\n')
            run(fasta, out, command='annotate', predictor=DemoORFPredictor(), use_mock_evidence=True, rna_annotation='none')
            manifest=json.loads((out/'run_manifest.json').read_text())
            self.assertIn('annotation_reliability', manifest)
            self.assertIn('annotation_reliability/index.html', (out/'report.html').read_text())
            self.assertEqual(manifest['annotation_reliability']['reference_status'], 'NOT_REQUESTED')
            self.assertEqual(manifest['quality_control']['coding_bases'], 126)


if __name__ == '__main__':
    unittest.main()
