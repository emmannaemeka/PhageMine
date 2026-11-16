"""PZ benchmark expectations are confined to tests, never runtime rules."""
import csv
import json
from pathlib import Path
from types import SimpleNamespace
from Bio import SeqIO
from Bio.Seq import Seq
from phagemine.models import Protein
from phagemine.start_selection import select_supported_starts
from phagemine.start_refinement import refine_starts, rescue_short_orfs, cds_problems
from phagemine.genbank import repair_boundaries


def fixture():
    root=Path(__file__).parent/'fixtures/genbank_real'
    genome=str(SeqIO.read(root/'PZ583544.1.fasta','fasta').seq)
    proteins=[]
    for row in csv.DictReader((root/'PZ583544.1.models.tsv').open(),delimiter='\t'):
        a,b=int(row['start']),int(row['end']);nt=genome[a-1:b]
        if row['strand']=='-':nt=str(Seq(nt).reverse_complement())
        proteins.append(Protein('PZ583544.1',row['protein_id'],a,b,row['strand'],nt,row['sequence'],'PHANOTATE'))
    secondary=[SimpleNamespace(**r, caller="prodigal_gv", identifier=f"GV{i}") for i,r in enumerate(json.loads((root/'PZ583544.1.gv_coordinates.json').read_text()))]
    return genome,proteins,secondary


def test_pz_starts_and_preserved_gene_ids():
    genome,raw,secondary=fixture()
    selected,_=select_supported_starts(genome,raw,secondary)
    refined,events,_=refine_starts(genome,selected,secondary)
    final=repair_boundaries(genome,refined)
    targets={'PM_000001':37,'PM_000008':6642,'PM_000009':6860,'PM_000010':7110,'PM_000069':42415}
    for p in final:
        if p.protein_id in targets:
            assert p.start==targets[p.protein_id]
            assert cds_problems(p,genome)==[]
    assert not final[0].partial_5prime
    assert {e['protein_id'] for e in events}=={'PM_000008','PM_000009','PM_000010','PM_000069'}
    assert [p.protein_id for p in final]==[p.protein_id for p in raw]


def test_pz_overlapping_short_gene_is_review_only():
    genome,raw,secondary=fixture()
    selected,_=select_supported_starts(genome,raw,secondary)
    refined,_,_=refine_starts(genome,selected,secondary)
    rescued,events,candidates=rescue_short_orfs(genome,refined)
    target=next(c for c in candidates if (c['start'],c['end'])==(43091,43192))
    assert target['length_nt']==102 and target['start_codon']=='ATG'
    assert not target['gap_eligible'] and not target['accepted']
    assert target['review_reason']=='SHORT_ORF_OVERLAPS_EXISTING_CDS_REVIEW'
    assert not events and len(rescued)==len(raw)
