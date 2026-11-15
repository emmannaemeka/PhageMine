"""Transparent structural confidence tiers, separate from product confidence."""
from dataclasses import replace
from Bio.Seq import Seq
from .start_refinement import oriented_bounds, rbs_score


def assign_confidence(genome, proteins, secondary=(), strict_precision=False):
    """Return retained models and an audit; never use reference annotations.

    Caller support means the same biological stop on the same strand. A start
    disagreement does not turn a jointly detected gene into a single-caller ORF.
    Rescue models are never described as ab initio caller consensus.
    """
    stops={(p.strand, p.end if p.strand=='+' else p.start) for p in secondary}
    retained, audit=[],[]
    oriented={'+':genome.upper(), '-':str(Seq(genome.upper()).reverse_complement())}
    for p in proteins:
        left,_=oriented_bounds(p,len(genome))
        rbs=rbs_score(oriented[p.strand],left)
        supported=any(e.supports and e.status.upper()=='REAL' and
                      e.source.lower() in {'phrogs','vogdb','pfam','swiss-prot','swissprot','cdd','interpro'}
                      for e in p.evidence)
        rescue=p.gene_call_source=='short_orf_rescue'
        both=not rescue and (p.strand,p.end if p.strand=='+' else p.start) in stops
        single=not rescue and not both
        low=single and len(p.cds)<150 and rbs['score']<8 and not supported
        tier='low' if low else ('high' if both or supported else 'medium')
        row={'protein_id':p.protein_id,'confidence_tier':tier,
             'caller_support':'Overlap-Rescue' if rescue and p.gene_call_parameters.get('overlap_rescue') else
                              ('Short-ORF-Rescue' if rescue else ('Both' if both else p.gene_call_source)),
             'rbs_score':rbs['score'],'protein_evidence_support':supported,
             'single_caller_micro_orf':single and len(p.cds)<150,
             'excluded':bool(strict_precision and low),
             'reason':'UNSUPPORTED_SINGLE_CALLER_MICRO_ORF' if low else 'SUPPORTED_OR_NON_MICRO_ORF',
             'start':p.start,'end':p.end,'strand':p.strand}
        audit.append(row)
        if not row['excluded']:
            retained.append(replace(p,gene_call_parameters={**p.gene_call_parameters,'structural_confidence':row}))
    return retained,audit
