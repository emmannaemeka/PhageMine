from __future__ import annotations
import csv,json
from pathlib import Path

def prioritize(enrichment_db, output):
    p=Path(enrichment_db)/'pmf_functional_enrichment.tsv'; rows=list(csv.DictReader(p.open(),delimiter='\t'))
    candidates=[]
    for r in rows:
        if r['family_functional_status']!='UNKNOWN_MULTI_GENOME': continue
        sources=[x for x in r.get('evidence_sources_present','').split(',') if x]
        completeness=1.0 if r.get('member_annotations') else 0.0
        components={'member_count':int(r['member_count']),'genome_count':int(r['genome_count']),'evidence_source_count':len(sources),'metadata_completeness':completeness}
        rationale=f"recurrent across {r['genome_count']} genomes with {len(sources)} accepted evidence sources but no consistent defensible function."
        candidates.append({**r,'accepted_evidence_sources':','.join(sources),'evidence_source_count':len(sources),'swissprot_status':'PRESENT' if 'Swiss-Prot' in sources else 'ABSENT_OR_NO_ACCEPTED_HIT','pfam_status':'PRESENT' if 'Pfam' in sources else 'ABSENT_OR_NO_ACCEPTED_HIT','vogdb_status':'PRESENT' if 'VOGDB' in sources else 'ABSENT_OR_NO_ACCEPTED_HIT','phrogs_status':'PRESENT' if 'PHROGs' in sources else 'ABSENT_OR_NO_ACCEPTED_HIT','metadata_completeness':completeness,'priority_components':json.dumps(components,sort_keys=True),'priority_rationale':rationale})
    candidates.sort(key=lambda r:(-int(r['genome_count']),-int(r['member_count']),-int(r['evidence_source_count']),r['family_id']))
    for i,r in enumerate(candidates,1): r['rank']=i
    out=Path(output);out.mkdir(parents=True,exist_ok=True); cols=['rank','family_id','member_count','genome_count','source_genomes','host_genera','family_functional_status','accepted_evidence_sources','evidence_source_count','swissprot_status','pfam_status','vogdb_status','phrogs_status','metadata_completeness','priority_components','priority_rationale']
    with (out/'pmf_unresolved_priority.tsv').open('w',newline='') as h:w=csv.DictWriter(h,fieldnames=cols,delimiter='\t');w.writeheader();w.writerows({k:r.get(k,'') for k in cols} for r in candidates)
    (out/'pmf_unresolved_priority.json').write_text(json.dumps([{k:r.get(k,'') for k in cols} for r in candidates],indent=2,sort_keys=True))
    return candidates
