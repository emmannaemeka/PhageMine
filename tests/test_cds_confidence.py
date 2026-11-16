from types import SimpleNamespace
from phagemine.models import Protein
from phagemine.cds_confidence import assign_confidence


def test_default_retains_low_confidence_and_strict_excludes_only_micro_orf():
    genome='ATG'+'AAA'*28+'TAA'+'CCC'*10+'ATG'+'AAA'*60+'TAA'
    ps=[Protein('g','PM_000001',1,90,'+',genome[:90],'M'+'K'*28,'PHANOTATE'),
        Protein('g','PM_000002',121,len(genome),'+',genome[120:],'M'+'K'*60,'PHANOTATE')]
    kept,audit=assign_confidence(genome,ps)
    assert len(kept)==2 and audit[0]['confidence_tier']=='low'
    kept,audit=assign_confidence(genome,ps,strict_precision=True)
    assert [p.protein_id for p in kept]==['PM_000002'] and audit[0]['excluded']
    # Same-stop support counts even when callers disagree on the start.
    kept,audit=assign_confidence(genome,ps,[SimpleNamespace(start=4,end=90,strand='+')],True)
    assert len(kept)==2 and audit[0]['caller_support']=='Both'
