import json
from pathlib import Path
from dataclasses import replace
import pytest
from Bio.Seq import Seq
from Bio import SeqIO
from phagemine.models import Protein, Evidence, EvidenceLevel
from phagemine.genbank import repair_boundaries, write_package
from phagemine.genome import translate
from phagemine.start_refinement import cds_problems, rescue_short_orfs, genome_start_context, candidate_scores
from phagemine.overlap_rescue import overlap_decision


def model(cds, start=1, strand='+', name='PM_000001'):
    return Protein('test', name, start, start+len(cds)-1, strand, cds, translate(cds), 'test')


@pytest.mark.parametrize('strand', ['+', '-'])
@pytest.mark.parametrize('padding', [1, 2])
@pytest.mark.parametrize('end', ['five', 'three'])
def test_partial_snapping_preserves_translation_and_frame(strand, padding, end, tmp_path):
    coding = 'CCCAAATAA' if end == 'five' else 'ATGAAACCC'
    oriented = 'G'*padding+coding if end == 'five' else coding+'G'*padding
    genome = oriented if strand == '+' else str(Seq(oriented).reverse_complement())
    left, right = (padding, len(oriented)) if end == 'five' else (0, len(coding))
    a,b = (left+1,right) if strand=='+' else (len(genome)-right+1,len(genome)-left)
    p = model(coding, a, strand); p.end=b
    p.partial_5prime=end=='five'; p.partial_3prime=end=='three'
    q=repair_boundaries(genome,[p])[0]
    assert q.start==1 and q.end==len(genome)
    assert q.codon_start==(padding+1 if end=='five' else 1)
    assert q.sequence==p.sequence and cds_problems(q,genome)==[]
    assert repair_boundaries(genome,[q])[0]==q
    write_package(tmp_path,'test',genome,[q],{'input_sha256':'test'})
    feature=next(f for f in SeqIO.read(tmp_path/'genbank_submission/genome.gbk','genbank').features if f.type=='CDS')
    nt=feature.extract(Seq(genome))[q.codon_start-1:]
    assert str(nt[:len(nt)//3*3].translate(table=11)).rstrip('*')==q.sequence


def test_ambiguity_resolution_and_internal_stops():
    assert translate('ATGGTKTAA')=='MV'
    assert translate('ATGNNNTAA')=='MX'
    assert translate('ATGTAATAA')=='M*'  # do not hide an internal stop
    from phagemine.gene_prediction import reverse_complement
    assert reverse_complement('TTAMACCAT') == 'ATGGTK TAA'.replace(' ', '')
    genome=str(Seq('ATGGTKTAA').reverse_complement())
    p=model('ATGGTKTAA', strand='-')
    assert cds_problems(repair_boundaries(genome,[p])[0],genome)==[]


def test_rna_mask_excludes_strong_rbs_rescue():
    seq='CCC'*10+'AGGAGG'+'A'*7+'ATG'+'AAA'*24+'TAA'
    out,_,rows=rescue_short_orfs(seq,[])
    assert len(out)==1
    mask=[{'start':out[0].start+10,'end':out[0].end-10,'type':'tmRNA'}]
    out,_,rows=rescue_short_orfs(seq,[],rna_features=mask)
    assert not out
    assert any(r['review_reason']=='SHORT_ORF_OVERLAPS_RNA_REJECTED' for r in rows)


def evidence(identifier='phrog_5672', coverage=.9, status='REAL'):
    return Evidence('homology','Rz1 outer membrane spanin',EvidenceLevel.COMPUTATIONAL,
                    'PHROGs','test',status,True,{'evalue':1e-10,'query_coverage':coverage},identifier)


def test_overlap_requires_support_frame_and_rna_clearance():
    outer=model('ATG'+'AAA'*98+'TAA'); outer.annotation='Rz-like spanin'
    inner=model('ATG'+'AAA'*24+'TAA',start=32,name='candidate')
    inner.sequence='MKKLAGC'+'K'*19
    row={'rbs_score':10}
    assert overlap_decision(inner,row,[outer])['product']=='Rz1-like spanin'
    assert overlap_decision(inner,{'rbs_score':0},[outer]) is None
    assert overlap_decision(replace(inner,start=33),row,[outer]) is None
    assert overlap_decision(inner,row,[outer],[{'start':40,'end':50}]) is None
    inner.sequence='M'+'K'*25; inner.evidence=[evidence(identifier='phrog_other')]
    assert overlap_decision(inner,{'rbs_score':0},[outer])['code']=='SHORT_ORF_DISTINCT_FAMILY_OVERLAP_PROMOTED'
    outer.evidence=[evidence(identifier='phrog_other')]
    assert overlap_decision(inner,{'rbs_score':0},[outer]) is None
    outer.evidence=[]; inner.evidence=[evidence(coverage=.69)]
    assert overlap_decision(inner,{'rbs_score':0},[outer]) is None
    inner.evidence=[evidence(status='mock')]
    assert overlap_decision(inner,{'rbs_score':0},[outer]) is None


def test_gc_sd_calibration_needs_anchors_and_preserves_moderate_gc():
    ps=[model('ATG'+'GCC'*25+'TAA',start=50+i*90,name=f'PM_{i:06d}') for i in range(12)]
    high='GCC'*400
    c=genome_start_context(high,ps,ps)
    assert c['sd_sparse'] and c['sd_weight']==0 and c['adaptation']==1
    assert genome_start_context(high,ps,[])['sd_weight']==1
    moderate='ATGC'*300
    assert genome_start_context(moderate,ps,ps)['sd_weight']==1
    assert genome_start_context(moderate,ps,ps)['adaptation']==0
