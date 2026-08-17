"""Conservative aggregation of existing member annotations/evidence."""
from __future__ import annotations
import csv,json
from collections import defaultdict
from pathlib import Path

def _load_result(root):
    root=Path(root); ev={}; cls={}
    if (root/'evidence.json').exists():
        for p in json.loads((root/'evidence.json').read_text()): ev[p.get('protein_id')]=p
    if (root/'functional_classification.json').exists():
        for p in json.loads((root/'functional_classification.json').read_text()): cls[p.get('protein_id')]=p
    return ev,cls

def enrich_database(database, result_dirs, output):
    root=Path(database); families=json.loads((root/'families.json').read_text())
    members=[]
    with (root/'family_members.tsv').open() as h: members=list(csv.DictReader(h,delimiter='\t'))
    results={}
    for d in result_dirs:
        ev,cl=_load_result(d); name=Path(d).name.upper()
        for token in ('AJP','C3','C6','IJEOMA','ANYANBIMPE'):
            if token in name: name=token; break
        results.update({f'{name}|{k}':(v,cl.get(k,{})) for k,v in ev.items()})
        results.update({k:(ev.get(k,{}),v) for k,v in cl.items()})
    byfam=defaultdict(list)
    for m in members:
        key=m.get('source_protein_id') or m.get('member_id'); p=results.get(key) or results.get(m.get('member_id'))
        byfam[m['family_id']].append((m,p))
    rows=[]; member_rows=[]
    for f in families:
        ms=byfam[f['family_id']]; anns=[]; sources=set(); states=[]; div=[]
        for m,p in ms:
            ev,cl=p or ({},{}); ann=cl.get('proposed_function') or ev.get('annotation')
            state=cl.get('functional_state');
            if ann: anns.append(ann)
            if state: states.append(state)
            accepted=[e for e in ev.get('evidence',[]) if e.get('supports')]
            sources.update(e.get('source') for e in accepted if e.get('source')); div.append(len({e.get('source') for e in accepted if e.get('source')}))
            member_rows.append({'family_id':f['family_id'],'member_id':m.get('member_id'),'source_genome_id':m.get('source_genome_id'),'source_protein_id':m.get('source_protein_id'),'annotation':ann or '','functional_state':state or '','evidence_sources':','.join(sorted({e.get('source') for e in accepted if e.get('source')})),'evidence_diversity':max(div[-1],0)})
        if not p: status='MISSING_ANNOTATION_METADATA'; consensus=''
        elif not anns and not states: status='MISSING_ANNOTATION_METADATA'; consensus=''
        elif anns and len({a.casefold() for a in anns})==1 and not all(('hypothetical' in a.casefold() or 'unknown' in a.casefold() or 'uncharacterized' in a.casefold()) for a in anns): status='KNOWN_CONSISTENT'; consensus=anns[0]
        elif anns and all(('hypothetical' in a.casefold() or 'unknown' in a.casefold() or 'uncharacterized' in a.casefold()) for a in anns):
            distinct={m.get('source_genome_id') for m,_ in ms}
            status='UNKNOWN_SINGLETON' if len(ms)==1 else ('UNKNOWN_MULTI_GENOME' if len(distinct)>1 else 'UNKNOWN_MULTI_MEMBER')
            consensus=''
        elif len({a.casefold() for a in anns})>1: status='CONFLICTING_MEMBER_ANNOTATIONS'; consensus=''
        else: status='MISSING_ANNOTATION_METADATA'; consensus=''
        rows.append({'family_id':f['family_id'],'member_count':len(ms),'genome_count':f.get('genome_count'),'host_genera':','.join(f.get('host_genera',[])),'family_functional_status':status,'consensus_annotation':consensus,'member_annotations':'; '.join(sorted(set(anns))),'evidence_sources_present':','.join(sorted(sources)),'evidence_diversity':max(div or [0]),'source_genomes':','.join(f.get('source_genome_ids',[]))})
    out=Path(output); out.mkdir(parents=True,exist_ok=True)
    (out/'pmf_functional_enrichment.json').write_text(json.dumps(rows,indent=2,sort_keys=True))
    cols=list(rows[0]) if rows else ['family_id'];
    with (out/'pmf_functional_enrichment.tsv').open('w',newline='') as h: w=csv.DictWriter(h,fieldnames=cols,delimiter='\t');w.writeheader();w.writerows(rows)
    with (out/'pmf_member_evidence.tsv').open('w',newline='') as h: w=csv.DictWriter(h,fieldnames=list(member_rows[0]) if member_rows else ['family_id'],delimiter='\t');w.writeheader();w.writerows(member_rows)
    for name,state in [('pmf_unresolved_families.tsv',None),('pmf_conflicting_families.tsv','CONFLICTING_MEMBER_ANNOTATIONS')]:
        with (out/name).open('w') as h: h.write('family_id\n'+'\n'.join(r['family_id'] for r in rows if (r['family_functional_status'].startswith('UNKNOWN_') if state is None else r['family_functional_status']==state))+('\n' if rows else ''))
    return rows
