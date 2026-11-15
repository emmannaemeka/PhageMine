from dataclasses import replace
import json
from pathlib import Path
import pytest
from phagemine.models import Protein, Evidence, EvidenceLevel
from phagemine.spanins import canonical_spanin, normalize_spanins


def pair(strand='+', offset=100, sequence='MKK'+'L'*10+'LAGC'+'A'*30):
    up=Protein('genome','up',1,300,'+','', 'M'+'A'*99,'test',annotation='rz-like spanin')
    down=Protein('genome','down',1+offset,250,'+','',sequence,'test',annotation='spanin')
    if strand=='-':
        up=replace(up,start=1001-up.end,end=1001-up.start,strand='-')
        down=replace(down,start=1001-down.end,end=1001-down.start,strand='-')
    return up,down


def test_capitalization():
    assert canonical_spanin('rz-like spanin')=='Rz-like spanin'
    assert canonical_spanin('rz spanin')=='Rz-like spanin'
    assert canonical_spanin('rz1-like spanin')=='Rz1-like spanin'
    assert canonical_spanin('unimolecular spanin')=='unimolecular spanin'


@pytest.mark.parametrize('strand',['+','-'])
def test_embedded_plus_one_lipobox(strand):
    up,down=pair(strand)
    before=[(p.start,p.end,p.sequence) for p in (up,down)]
    normalize_spanins([down,up])
    assert (up.annotation,down.annotation)==('Rz-like spanin','Rz1-like spanin')
    assert [(p.start,p.end,p.sequence) for p in (up,down)]==before
    assert down.gene_call_parameters['spanin_product_normalization']['relative_frame']==1
    assert normalize_spanins([up,down])==[]


@pytest.mark.parametrize('offset',[99,101])
def test_wrong_frame_not_upgraded(offset):
    up,down=pair(offset=offset);normalize_spanins([up,down]);assert down.annotation=='spanin'


def test_adjacent_supported_spanin_and_distance_limit():
    up,down=pair();down.start=310;down.end=450
    normalize_spanins([up,down]);assert down.annotation=='Rz1-like spanin'
    up,down=pair();down.start=317;down.end=450
    normalize_spanins([up,down]);assert down.annotation=='spanin'


def test_no_support_and_cross_contig_or_strand():
    up,down=pair(sequence='M'+'A'*40);normalize_spanins([up,down]);assert down.annotation=='spanin'
    for field,value in [('genome_id','another'),('strand','-')]:
        up,down=pair();setattr(down,field,value);normalize_spanins([up,down]);assert down.annotation=='spanin'


def test_weak_or_mock_family_not_used():
    for state,ev in [('mock',1e-30),('REAL',.1)]:
        up,down=pair(sequence='M'+'A'*40)
        down.evidence=[Evidence('phage_orthology','spanin',EvidenceLevel.COMPUTATIONAL,'PHROGs','test',status=state,supports=True,identifier='phrog_5672',metrics={'evalue':ev,'query_coverage':.9})]
        normalize_spanins([up,down]);assert down.annotation=='spanin'


def test_pz_real_evidence_and_export(tmp_path):
    from phagemine.fusion import classify_proteins
    from phagemine.genbank import write_package
    from unittest.mock import patch
    from Bio import SeqIO
    root=Path(__file__).parent/'fixtures/genbank_real'
    proteins=[]
    for row in json.loads((root/'PZ583544.1.spanins.json').read_text()):
        row['annotation_level']=EvidenceLevel(row['annotation_level'])
        row['evidence']=[Evidence(**{**e,'level':EvidenceLevel(e['level'])}) for e in row['evidence']]
        proteins.append(Protein(**row))
    results=classify_proteins(proteins)
    assert {r['protein_id']:r['display_product'] for r in results}=={'PM_000009':'Rz-like spanin','PM_000010':'Rz1-like spanin'}
    genome=str(SeqIO.read(root/'PZ583544.1.fasta','fasta').seq)
    with patch('phagemine.genbank.table2asn_status',return_value={'state':'unavailable'}):
        write_package(tmp_path,'PZ583544.1',genome,proteins,{'input_sha256':'test'})
    rec=SeqIO.read(tmp_path/'genbank_submission/genome.gbk','genbank')
    names={f.qualifiers['locus_tag'][0]:f.qualifiers['product'][0] for f in rec.features if f.type=='CDS'}
    assert names=={'PM_000009':'Rz-like spanin','PM_000010':'Rz1-like spanin'}
