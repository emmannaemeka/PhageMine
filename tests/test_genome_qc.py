from phagemine.genome_qc import (sequence_findings, terminal_repeats, frameshift_candidates,
                                sanitize_product, safety_markers, curate_regulatory)
from phagemine.start_refinement import enumerate_starts, intergenic_candidate, RefinementConfig
from types import SimpleNamespace
from Bio.Seq import Seq


def test_ambiguous_intersection_and_homopolymer():
    rows=sequence_findings('ACGNNRY'+'A'*12,[dict(id='p1',start=4,end=7)])
    assert rows[0]['cds_ids']==['p1']
    assert rows[1]['code']=='HOMOPOLYMER'


def test_direct_and_inverted_repeats():
    motif='ACGTAGCTAGCCTAGGATCG'
    assert terminal_repeats(motif+'T'*40+motif)[0]['code']=='DTR_CANDIDATE'
    assert any(r['code']=='ITR_CANDIDATE' for r in terminal_repeats(motif+'T'*40+str(Seq(motif).reverse_complement())))
    assert terminal_repeats('N'*100)==[]
    assert terminal_repeats(motif+'T'*20+motif+'C'*20+motif)[0]['internal_copies']==[41]


def test_frameshift_requires_ordered_same_subject_different_frame():
    genes=[dict(id='a',start=1,end=90,strand='+'),dict(id='b',start=95,end=184,strand='+')]
    hits={'a':[dict(target='X',source='db',subject_start=1,subject_end=30,evalue=1e-10)],'b':[dict(target='X',source='db',subject_start=31,subject_end=60,evalue=1e-10)]}
    assert frameshift_candidates(genes,hits)[0]['code']=='PUTATIVE_FRAMESHIFT'
    hits['b'][0]['target']='Y'
    assert frameshift_candidates(genes,hits)==[]
    hits['b'][0]['target']='X'; genes[1]['start']=94
    assert frameshift_candidates(genes,hits)==[]


def test_gap_floor_and_no_nested_rescue():
    seq='ATG'+'AAA'*24+'TAA' # 78 nt, above the guide's 75-nt floor
    assert enumerate_starts(seq,len(seq),RefinementConfig().rescue_min_nt)==[0]
    assert enumerate_starts(seq,len(seq),90)==[]
    occupied=[SimpleNamespace(start=1,end=100),SimpleNamespace(start=200,end=300)]
    assert intergenic_candidate(101,199,occupied)
    assert intergenic_candidate(97,203,occupied)
    assert not intergenic_candidate(50,140,occupied)
    assert not intergenic_candidate(20,90,occupied)
    assert not intergenic_candidate(1,300,occupied)


def test_products():
    for value in ('putative phage protein','YP_123456.1','gp200','ORF184','Mcp','hypothetical protein SA5_0153'):
        assert sanitize_product(value)=='hypothetical protein'
    assert sanitize_product('DNA polymerase')=='DNA polymerase'
    assert sanitize_product('integrase')=='integrase'
    assert sanitize_product('DNA polymerase',False)=='hypothetical protein'


def test_safety_does_not_claim_absence_or_trust_mock():
    evidence=[dict(status='REAL',supports=True,description='site-specific recombinase',source='db'),dict(status='REAL',supports=True,description='Shiga toxin',source='db'),dict(status='REAL',supports=True,description='beta-lactamase',source='db')]
    report=safety_markers([dict(protein_id='p',evidence=evidence)])
    assert {r['category'] for r in report['findings']}=={'lysogeny','toxin_virulence','amr'}
    assert report['therapeutic_safety']=='NOT_ESTABLISHED'
    assert safety_markers([dict(protein_id='p',evidence=[dict(status='mock',supports=True,description='toxin')])])['findings']==[]


def test_regulatory_overlap_rejected():
    accepted,rejected=curate_regulatory([dict(start=1,end=20),dict(start=25,end=45)],[dict(start=30,end=100)])
    assert len(accepted)==len(rejected)==1
    assert 'computational' in accepted[0]['note']


def test_aragorn_tmrna_parser():
    import pytest
    from phagemine.rna_features import parse_aragorn_tmrna
    rows=parse_aragorn_tmrna('>seq\n1 genes found\n1 tmRNA c[12,310] 90,120\nTAGPEPTIDE\n','test')
    assert rows[0]['strand']=='-' and rows[0]['type']=='tmRNA'
    assert parse_aragorn_tmrna('>seq\n0 genes found\n','test')==[]
    with pytest.raises(ValueError): parse_aragorn_tmrna('>seq\n1 genes found\n','test')


def test_dedicated_safety_requires_identity_and_both_coverages():
    from phagemine.safety_screen import parse_hits
    hit='gnl|PhageMine|PM_000001\ttoxin1\t80\t1\t90\t100\t1\t90\t100\t1e-30'
    assert parse_hits(hit,'toxin_virulence')[0]['id']=='PM_000001'
    assert parse_hits(hit.replace('\t1\t90\t100\t1e-30','\t1\t10\t100\t1e-30'),'toxin_virulence')==[]


def test_assembly_command_and_no_overwrite(tmp_path, monkeypatch):
    from phagemine.assembly import assemble
    import subprocess
    import pytest
    for name in ('r1','r2'): (tmp_path/name).write_text('@read\nATG\n+\nIII\n')
    monkeypatch.setattr('phagemine.assembly.shutil.which',lambda exe:'/bin/shovill')
    def run(command,**kwargs):
        out=tmp_path/'assembly';out.mkdir();(out/'contigs.fasta').write_text('>c\nATG\n')
        assert command[command.index('--depth')+1]=='100'
        return subprocess.CompletedProcess(command,0,'done','')
    monkeypatch.setattr('phagemine.assembly.subprocess.run',run)
    out=assemble(tmp_path/'r1',tmp_path/'r2',tmp_path/'assembly')
    assert out.is_file()
    with pytest.raises(ValueError):assemble(tmp_path/'r1',tmp_path/'r2',tmp_path/'assembly')


def test_benchmark_same_stop_and_strand():
    from phagemine.cds_benchmark import compare
    ref=[dict(id='ref',contig='c',start=1,end=99,strand='+',start_codon='ATT')]
    pred=[dict(id='p',contig='c',start=4,end=99,strand='+',start_codon='ATG')]
    report=compare(pred,ref)
    assert report['stop_only_count']==1 and report['exact_count']==0
    pred[0]['strand']='-'
    report=compare(pred,ref)
    assert len(report['missed'])==len(report['unmatched_predictions'])==1


def test_joined_reference_not_counted_as_exact_spanning_cds():
    from phagemine.cds_benchmark import compare
    ref=[dict(id='r',contig='c',start=1,end=99,strand='+',parts=[[1,30,1],[40,99,1]])]
    pred=[dict(id='p',contig='c',start=1,end=99,strand='+')]
    report=compare(pred,ref)
    assert report['exact_count']==0 and len(report['structure_mismatch'])==1


def test_products_keep_function_and_drop_locus_designations():
    assert sanitize_product('baseplate wedge protein gp25') == 'baseplate wedge protein'
    assert sanitize_product("5'-3' exonuclease gp43") == "5'-3' exonuclease"
    assert sanitize_product('nucleotide kinase gp1.7') == 'nucleotide kinase'
    assert sanitize_product('terminase, small subunit gp18') == 'terminase, small subunit'
    assert sanitize_product('probable major capsid protein gp17') == 'putative major capsid protein'
    assert sanitize_product('tail fiber protein {ECO:0000303|PubMed:1}') == 'tail fiber protein'
    assert sanitize_product('Nrd.1 hypothetical protein') == 'hypothetical protein'
    assert sanitize_product('conserved phage protein of unknown function') == 'hypothetical protein'
    for value in ('gp5.4', 'ORF 12', 'protein gp7', 'phrog_123'):
        assert sanitize_product(value) == 'hypothetical protein'
    assert sanitize_product(sanitize_product('baseplate wedge protein gp25')) == 'baseplate wedge protein'
