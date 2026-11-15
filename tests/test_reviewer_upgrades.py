import json
import shutil
from pathlib import Path
from phagemine.start_refinement import attenuate_rbs,coupling_score
from phagemine.submission_metadata import load_metadata
from phagemine.models import Protein
from phagemine.genbank import write_package
from phagemine.reviewer_package import build_reviewer_package,COLUMNS


def test_gc_and_coupling():
    weak={'score':8.,'motif':'GGAG'};strong={'score':12.,'motif':'AGGAGG'}
    assert attenuate_rbs(weak,{'gc_fraction':.64})['score']==2
    assert attenuate_rbs(weak,{'gc_fraction':.58})['score']==8
    assert attenuate_rbs(strong,{'gc_fraction':.64})['score']==12
    assert coupling_score('TAATG',2,-1)==4
    assert coupling_score('CCATG',2,-1)==0
    assert coupling_score('ATGA',0,-4)==4


def test_metadata_and_workbook(tmp_path):
    from openpyxl import load_workbook
    config=tmp_path/'metadata.json';config.write_text(json.dumps({'topology':None,'organism':'test phage','host':'test host','country':'Nigeria','isolation_source':'water'}))
    m=load_metadata(config);assert m.topology=='linear' and m.geographic_location=='Nigeria'
    seq='ATG'+'AAA'*28+'TAA';p=Protein('test','PM_000001',1,len(seq),'+',seq,'M'+'K'*28,'PHANOTATE')
    write_package(tmp_path,'test',seq,[p],{'input_sha256':'test'},m)
    package=tmp_path/'genbank_submission'
    assert '[organism=test phage]' in (package/'genome.fsa').read_text()
    assert '\t\t\tgeo_loc_name\tNigeria' in (package/'features.tbl').read_text()
    (tmp_path/'proteins.faa').write_text('>PM_000001\n'+p.sequence+'\n')
    (tmp_path/'cds.fna').write_text('>PM_000001\n'+seq+'\n')
    (tmp_path/'genes.gff3').write_text('##gff-version 3\n')
    (tmp_path/'phagemine_qc_report.json').write_text('{}');(tmp_path/'phagemine_qc_report.md').write_text('QC fixture')
    out=build_reviewer_package(tmp_path)
    wb=load_workbook(out/'test_manual_verification.xlsx')
    assert wb.sheetnames==['CDS_Verification','tRNA_tmRNA','Review_Queue_Candidates','Genome_Summary_and_QC']
    ws=wb['CDS_Verification'];assert [c.value for c in ws[1]]==COLUMNS
    assert ws.freeze_panes=='B2' and ws.auto_filter.ref
    assert ws.cell(2,COLUMNS.index('Codon_Verifier_Status')+1).value=='PASS'
    assert len(list(out.iterdir()))==9


def test_overlap_family_exception_without_rbs(tmp_path,monkeypatch):
    from types import SimpleNamespace
    from phagemine.models import Evidence,EvidenceLevel
    from phagemine.overlap_rescue import promote_overlaps
    import phagemine.fusion
    genome='C'*31+'ATG'+'AAA'*24+'TAA'+'C'*200
    outer=Protein('g','PM_000001',1,300,'+',genome[:300],'X'*100,'PHANOTATE')
    monkeypatch.setattr(phagemine.fusion,'classify_proteins',lambda ps:[{'protein_id':p.protein_id,'display_product':'hypothetical protein'} for p in ps])
    row={'start':32,'end':109,'strand':'+','score':0,'rbs_score':0,'gap_eligible':False,'accepted':False}
    (tmp_path/'short_orf_candidates.json').write_text(json.dumps([row]))
    def analyze(ps):
        return [SimpleNamespace(adapter='VOGDB',status='SUCCESS',provenance={},evidence=[Evidence('homology','family',EvidenceLevel.COMPUTATIONAL,'VOGDB','test','REAL',True,{'evalue':1e-10,'query_coverage':1},'VOG_TEST',provenance={'protein_id':p.protein_id}) for p in ps])]
    result=promote_overlaps(genome,[outer],tmp_path,analyze)
    assert len(result)==2
    assert result[-1].gene_call_parameters['overlap_rescue']['code']=='SHORT_ORF_DISTINCT_FAMILY_OVERLAP_PROMOTED'
    audit=json.loads((tmp_path/'overlap_rescue.json').read_text())
    assert audit['prefilter']['family_seed_queries']==1
    assert audit['prefilter']['confirmation_queries']==0


def test_yaml_metadata(tmp_path):
    p=tmp_path/'metadata.yaml';p.write_text('organism: test phage\nhost: test host\ncountry: Nigeria\n')
    assert load_metadata(p).geographic_location=='Nigeria'
