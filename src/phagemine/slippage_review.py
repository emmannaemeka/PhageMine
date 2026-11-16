"""Split homology identifies review candidates, not proven recoding events."""
from .genome_qc import frameshift_candidates


def annotate_slippage_review(proteins):
    cds=[{'id':p.protein_id,'start':p.start,'end':p.end,'strand':p.strand} for p in proteins]
    hits={}
    for p in proteins:
        for e in p.evidence:
            m=dict(e.metrics)
            if 'subject_coordinates' in m:
                m.update(subject_start=m['subject_coordinates']['start'],subject_end=m['subject_coordinates']['end'])
            if e.status.upper()=='REAL' and e.supports:
                hits.setdefault(p.protein_id,[]).append({**m,'target':e.identifier,'source':e.source})
    alerts=frameshift_candidates(cds,hits)
    for alert in alerts:
        alert['proposed_qualifier']='ribosomal_slippage'
        alert['status']='REVIEW_REQUIRED_NOT_ASSERTED'
        for p in proteins:
            if p.protein_id in alert['cds_ids']:
                p.gene_call_parameters.setdefault('ribosomal_slippage_review',[]).append(alert)
    return alerts
