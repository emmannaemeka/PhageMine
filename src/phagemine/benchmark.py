"""Import-only benchmarking for PhageMine and external annotations."""
from __future__ import annotations
import csv,json,re
from pathlib import Path
from hashlib import sha256

def import_gff(path, tool, genome_id=None):
    records=[]; path=Path(path)
    for line in path.read_text().splitlines():
        if not line or line.startswith('#'): continue
        f=line.split('\t')
        if len(f)<9 or f[2].lower() not in {'cds','gene'}: continue
        attrs=dict((x.split('=',1) for x in f[8].split(';') if '=' in x)); start,end=int(f[3]),int(f[4]); seq=attrs.get('translation','')
        records.append({'genome_id':genome_id or f[0],'tool':tool,'tool_version':None,'locus_id':attrs.get('ID') or attrs.get('locus_tag') or f"{tool}_{len(records)+1}",'start':start,'end':end,'strand':f[6],'protein_length':len(seq) if seq else None,'protein_sequence_checksum':sha256(seq.encode()).hexdigest() if seq else None,'product':attrs.get('product'),'gene':attrs.get('gene'),'source_annotation':attrs,'original_identifier':attrs.get('ID')})
    return records

def import_prokka(path, genome_id=None): return import_gff(path, 'Prokka', genome_id)
def import_pharokka(path, genome_id=None): return import_gff(path, 'Pharokka', genome_id)

def import_phagemine(results):
    root=Path(results); required=[root/'genes.gff3']
    missing=[str(p) for p in required if not p.is_file()]
    if missing: raise ValueError('PhageMine benchmark import missing required files: '+', '.join(missing))
    annotations={}
    if (root/'annotation.tsv').is_file():
        with (root/'annotation.tsv').open() as h:
            for row in csv.DictReader(h,delimiter='\t'): annotations[row.get('protein_id')]=row
    classification={}
    classification_available=(root/'functional_classification.tsv').is_file()
    if classification_available:
        with (root/'functional_classification.tsv').open() as h:
            for row in csv.DictReader(h,delimiter='\t'): classification[row.get('protein_id')]=row
    records=import_gff(root/'genes.gff3','PHAGEMINE')
    for record in records:
        ident=record['locus_id']; a=annotations.get(ident,{}) ; c=classification.get(ident,{})
        record.update({'product':a.get('annotation'),'functional_state':c.get('functional_state') if classification_available else None,'proposed_function':c.get('proposed_function') if classification_available else None,'confidence':c.get('confidence') if classification_available else None,'functional_classification_status':'AVAILABLE' if classification_available else 'LEGACY_OUTPUT_NOT_AVAILABLE','source_annotation':{'annotation':a,'classification':c if classification_available else None}})
    return records

def _relation(a,b):
    ov=max(0,min(a['end'],b['end'])-max(a['start'],b['start'])+1)
    if a['start']==b['start'] and a['end']==b['end'] and a['strand']==b['strand']: kind='EXACT_MATCH'
    elif a['strand']==b['strand'] and a['start']==b['start'] and a['end']!=b['end']: kind='STOP_DIFFERENCE'
    elif a['strand']==b['strand'] and a['end']==b['end'] and a['start']!=b['start']: kind='START_DIFFERENCE'
    # A one-base contact at a shared endpoint is not a meaningful overlap.
    elif ov > 1 and a['strand']!=b['strand']: kind='STRAND_CONFLICT'
    elif ov > 1 and a['strand']==b['strand']: kind='START_AND_STOP_DIFFERENCE'
    else: kind=None
    return kind,ov

def compare_models(a,b):
    rows=[]; used=set(); pending=set(range(len(a)))
    # Resolve globally in deterministic priority order so a weak overlap
    # cannot consume a target before an exact coordinate match is considered.
    for priority in ('EXACT_MATCH','STOP_DIFFERENCE','START_DIFFERENCE','START_AND_STOP_DIFFERENCE','STRAND_CONFLICT'):
        for i in sorted(pending):
            candidates=[(j,y) for j,y in enumerate(b) if j not in used and _relation(a[i],y)[0]==priority]
            if len(candidates)==1:
                j,y=candidates[0]; used.add(j); pending.remove(i); rows.append({'method_a':a[i],'method_b':y,'relationship':priority})
            elif len(candidates)>1 and priority in ('START_AND_STOP_DIFFERENCE','STRAND_CONFLICT'):
                pending.remove(i); rows.append({'method_a':a[i],'method_b':None,'relationship':'COMPLEX'})
    for i in sorted(pending):
        candidates=[j for j,y in enumerate(b) if j not in used and _relation(a[i],y)[1]>1]
        if not candidates: rows.append({'method_a':a[i],'method_b':None,'relationship':'METHOD_A_ONLY'})
        else: rows.append({'method_a':a[i],'method_b':None,'relationship':'COMPLEX'})
    for i,y in enumerate(b):
        if i not in used: rows.append({'method_a':None,'method_b':y,'relationship':'METHOD_B_ONLY'})
    return rows

def metrics(predicted, reference):
    exact=sum(r['relationship']=='EXACT_MATCH' for r in compare_models(predicted,reference)); p=len(predicted); r=len(reference); precision=exact/p if p else 0; recall=exact/r if r else 0
    return {'predicted':p,'reference':r,'exact_matches':exact,'exact_precision':precision,'exact_recall':recall,'exact_f1':2*precision*recall/(precision+recall) if precision+recall else 0}

def benchmark(methods, output, references=None):
    output=Path(output); output.mkdir(parents=True,exist_ok=True); references=references or {}
    comparisons=[]; summary=[]; metric_rows=[]
    names=list(methods)
    for i,a in enumerate(names):
        summary.append({'method':a,'cds_count':len(methods[a]),'specific_annotations':sum(bool(x.get('product')) for x in methods[a]),'hypothetical':sum(not x.get('product') or 'hypothetical' in x.get('product','').lower() for x in methods[a])})
        if a in references: metric_rows.append({'method':a,**metrics(methods[a],references[a])})
        for b in names[i+1:]: comparisons.extend([{**r,'method_a_name':a,'method_b_name':b} for r in compare_models(methods[a],methods[b])])
    for name,data in [('benchmark_summary',summary),('gene_model_comparison',comparisons),('benchmark_metrics',metric_rows)]:
        (output/(name+'.json')).write_text(json.dumps(data,indent=2,sort_keys=True,default=str))
        cols=sorted({k for row in data for k in row}) if data else ['method'];
        with (output/(name+'.tsv')).open('w',newline='') as h:
            w=csv.DictWriter(h,fieldnames=cols,delimiter='\t'); w.writeheader(); w.writerows(row for row in data)
    (output/'hypothetical_resolution.tsv').write_text('method\tcds_count\thypothetical\n'+'\n'.join(f"{x['method']}\t{x['cds_count']}\t{x['hypothetical']}" for x in summary)+'\n')
    (output/'functional_comparison.json').write_text(json.dumps([],indent=2)); (output/'functional_comparison.tsv').write_text('method\tstatus\n')
    legacy=[name for name,data in methods.items() if name == 'PHAGEMINE' and any(x.get('functional_classification_status') == 'LEGACY_OUTPUT_NOT_AVAILABLE' for x in data)]
    note='\n\nLegacy PhageMine output detected: gene-model comparison is supported, but modern functional-classification metrics are NOT_EVALUABLE.\n' if legacy else ''
    (output/'benchmark_report.md').write_text('# PhageMine benchmark\n\nThis report separates gene-model observations from reference-supported metrics.'+note+'\n')
    return summary
