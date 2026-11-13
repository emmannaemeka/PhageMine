import importlib.util
import json
import os
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
from types import SimpleNamespace

from phagemine.models import EvidenceLevel, Protein, SubmissionMetadata

import phagemine.genbank as export
MODULE = Path(export.__file__)

def protein(**kw):
    args = dict(genome_id='fixture', protein_id='PM_000001', start=1, end=9, strand='+', cds='ATGAAATAA', sequence='MK', gene_call_source='fixture', annotation='example protein', annotation_level=EvidenceLevel.COMPUTATIONAL)
    args.update(kw)
    return Protein(**args)

def metadata():
    return SubmissionMetadata(organism='Escherichia phage Test', phage_name='Test', host='Escherichia coli', isolation_source='water', collection_date='2025', geographic_location='Nigeria', sequencing_technology='Illumina', submitter_name='Test Author', submitter_email='test@example.org', authors=['Test Author'], locus_tag_prefix='TEST')

class ExportTests(unittest.TestCase):
    def test_final_computational_product_is_preserved(self):
        self.assertEqual(export.product_name(protein()), 'example protein')
        for level in (EvidenceLevel.WEAK, EvidenceLevel.HYPOTHESIS):
            self.assertEqual(export.product_name(protein(annotation_level=level)), 'hypothetical protein')
        self.assertEqual(export.product_name(protein(annotation='example MOCK')), 'hypothetical protein')

    def test_gene_and_cds_share_coordinates_and_tag(self):
        table = export.feature_table('fixture', [protein()], locus_tag_prefix='TEST_')
        self.assertIn('1\t9\tgene\n', table)
        self.assertIn('1\t9\tCDS\n', table)
        self.assertEqual(table.count('locus_tag\tTEST_000001'), 2)
        self.assertIn('product\texample protein', table)
        self.assertNotIn('\tgene\texample protein', table)

    def test_partial_table_columns_on_both_strands(self):
        for strand, coords in [('+', '<1\t>9'), ('-', '<9\t>1')]:
            p = protein(strand=strand, partial_5prime=True, partial_3prime=True, codon_start=2)
            table = export.feature_table('fixture', [p])
            self.assertIn(coords+'\tgene', table)
            self.assertIn(coords+'\tCDS', table)
            self.assertIn('codon_start\t2', table)

    def test_warning_state_is_reported_without_false_failure(self):
        with tempfile.TemporaryDirectory() as tmp:
            with patch.object(export, 'table2asn_status', return_value={'state':'passed', 'official_ncbi_validation':True}):
                result=export.write_package(tmp, 'fixture', 'ATGAAATAA'+'ACG'*100, [protein(annotation_level=EvidenceLevel.HYPOTHESIS)], {'input_sha256':'fixture'}, metadata())
            report=(Path(tmp)/'genbank_submission/NCBI_VALIDATION_REPORT.txt').read_text()
            self.assertEqual(result['submission_readiness'], 'READY_WITH_WARNINGS')
            self.assertIn('OVERALL: READY_WITH_WARNINGS', report)
            self.assertNotIn('NOT READY FOR', report)

    def test_missing_metadata_is_visible_and_blocks(self):
        with tempfile.TemporaryDirectory() as tmp:
            with patch.object(export, 'table2asn_status', return_value={'state':'unavailable'}):
                result=export.write_package(tmp, 'fixture', 'ATGAAATAA'+'ACG'*100, [protein()], {'input_sha256':'fixture'})
            report=(Path(tmp)/'genbank_submission/NCBI_VALIDATION_REPORT.txt').read_text()
            self.assertEqual(result['submission_readiness'], 'INCOMPLETE_METADATA')
            self.assertFalse(result['phagemine_pre_submission_validation']['ncbi_submission_ready'])
            self.assertIn('Missing required metadata: organism', report)

    def test_real_cds_error_still_blocks(self):
        pre=export.validate('fixture', 'ATGAAAAAA', [protein(cds='ATGAAAAAA', sequence='MKK')], {'input_sha256':'fixture'})
        self.assertFalse(pre['valid'])
        self.assertEqual(export.readiness(pre, metadata(), {'state':'passed'})[0], 'VALIDATION_FAILED')

    def test_partial_feature_away_from_boundary_is_blocked(self):
        p=protein(start=2,end=10,partial_5prime=True)
        result=export.validate('fixture', 'AATGAAATAA', [p], {'input_sha256':'fixture'})
        self.assertIn('NCBI_PARTIAL_LOCATION', {item['code'] for item in result['errors']})

    def test_unresolved_complete_start_requires_review(self):
        p=protein(cds='AAAAAATAA',sequence='KK')
        result=export.validate('fixture', p.cds, [p], {'input_sha256':'fixture'})
        self.assertEqual(export.readiness(result,metadata(),{'state':'passed'})[0], 'REQUIRES_REVIEW')

    def test_missing_official_validation_cannot_be_ready(self):
        state,_=export.readiness({'valid':True,'warnings':[]},metadata(),{'state':'unavailable'})
        self.assertEqual(state,'PENDING_NCBI_VALIDATION')

    def test_source_fields_are_exported(self):
        with tempfile.TemporaryDirectory() as tmp:
            with patch.object(export,'table2asn_status',return_value={'state':'unavailable'}):
                export.write_package(tmp,'fixture','ATGAAATAA',[protein()],{'input_sha256':'fixture'},metadata())
            header=(Path(tmp)/'genbank_submission/genome.fsa').read_text().splitlines()[0]
            self.assertIn('[host=Escherichia coli]',header)
            self.assertIn('[geo_loc_name=Nigeria]',header)
            self.assertIn('[collection_date=2025]',header)

    def test_review_item_still_blocks(self):
        pre={'valid':True, 'warnings':[{'severity':'REQUIRES_REVIEW'}]}
        self.assertEqual(export.readiness(pre, metadata(), {'state':'passed'})[0], 'REQUIRES_REVIEW')

    def test_protein_ids_and_products_match_table(self):
        with tempfile.TemporaryDirectory() as tmp:
            with patch.object(export, 'table2asn_status', return_value={'state':'unavailable'}):
                export.write_package(tmp, 'fixture', 'ATGAAATAA'+'ACG'*100, [protein()], {'input_sha256':'fixture'}, metadata())
            root=Path(tmp)/'genbank_submission'
            self.assertIn('>gnl|PhageMine|TEST_000001 example protein', (root/'proteins.faa').read_text())
            self.assertIn('protein_id\tgnl|PhageMine|TEST_000001', (root/'features.tbl').read_text())
            self.assertEqual((root/'features.tbl').read_bytes(), (root/'genome.tbl').read_bytes())

    def test_table2asn_zero_exit_with_errors_is_failure(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp); output=root/'table2asn_output'; output.mkdir()
            (output/'fixture.sqn').write_text('fixture')
            (output/'fixture.val').write_text('ERROR: fixture validation failure\n')
            with patch.object(export.tempfile, 'mkdtemp', return_value=str(output)), patch.object(export.subprocess,'run',return_value=SimpleNamespace(returncode=0,stdout='',stderr='')):
                result=export.table2asn_status(root, __file__)
            self.assertEqual(result['state'],'failed')
            self.assertFalse(result['official_ncbi_validation'])

    def test_rna_attachment_keeps_converter_table_synchronized(self):
        spec=importlib.util.spec_from_file_location('phagemine.rna_export_under_test', MODULE.parent/'rna_features.py')
        rna=importlib.util.module_from_spec(spec); spec.loader.exec_module(rna)
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp); package=root/'genbank_submission'; package.mkdir()
            (package/'features.tbl').write_text('>Feature fixture\n')
            rna.append_submission_features(root, [{'start':1,'end':9,'strand':'+','type':'tRNA','product':'tRNA-X','source':'fixture','source_version':'1','feature_id':'RNA1'}])
            self.assertEqual((package/'features.tbl').read_bytes(), (package/'genome.tbl').read_bytes())
            self.assertIn('1\t9\ttRNA', (package/'genome.tbl').read_text())

    def test_official_template_can_be_supplied(self):
        with tempfile.TemporaryDirectory() as tmp:
            template=Path(tmp)/'official.sbt'; template.write_text('fixture template')
            with patch.dict(os.environ, {'PHAGEMINE_SUBMISSION_TEMPLATE':str(template)}):
                self.assertEqual(export._submission_template(metadata()), 'fixture template')

    def test_fatal_discrepancy_blocks_zero_exit(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp); output=root/'run'; output.mkdir()
            (output/'genome.sqn').write_text('fixture')
            (output/'genome.val').write_text('')
            (output/'genome.dr').write_text('FATAL: fixture discrepancy')
            with patch.object(export.tempfile,'mkdtemp',return_value=str(output)), patch.object(export.subprocess,'run',return_value=SimpleNamespace(returncode=0,stdout='',stderr='')):
                self.assertEqual(export.table2asn_status(root,__file__)['state'],'failed')

    def test_rna_changes_are_revalidated(self):
        spec=importlib.util.spec_from_file_location('phagemine.rna_export_under_test',MODULE.parent/'rna_features.py')
        rna=importlib.util.module_from_spec(spec); spec.loader.exec_module(rna)
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp); package=root/'genbank_submission'; package.mkdir()
            (package/'features.tbl').write_text('>Feature fixture\n')
            (package/'submission_metadata.json').write_text(json.dumps(metadata().__dict__))
            (package/'validation.json').write_text(json.dumps({'submission_readiness':'READY', 'phagemine_pre_submission_validation':{'valid':True,'warnings':[],'errors':[]}}))
            with patch('phagemine.genbank.table2asn_status',return_value={'state':'failed','validation_findings':['Error: fixture RNA failure']}), patch('phagemine.genbank.readiness',side_effect=export.readiness):
                rna.append_submission_features(root,[{'start':1,'end':9,'strand':'+','type':'tRNA','product':'tRNA-X','source':'fixture','source_version':'1','feature_id':'RNA1'}])
            result=json.loads((package/'validation.json').read_text())
            self.assertEqual(result['submission_readiness'],'VALIDATION_FAILED')
            self.assertFalse(result['phagemine_pre_submission_validation']['ncbi_submission_ready'])
            self.assertIn('Error: fixture RNA failure',(package/'NCBI_VALIDATION_REPORT.txt').read_text())

    def test_table2asn_requires_output_artifact(self):
        with tempfile.TemporaryDirectory() as tmp:
            with patch.object(export.subprocess,'run',return_value=SimpleNamespace(returncode=0,stdout='',stderr='')):
                self.assertEqual(export.table2asn_status(Path(tmp), __file__)['state'],'failed')

if __name__=='__main__': unittest.main()
