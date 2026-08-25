"""Import-only benchmarking for PhageMine and external annotations."""
from __future__ import annotations
import csv,json,re
from pathlib import Path
from hashlib import sha256
from .hallmarks import HALLMARKS

def import_gff(path, tool, genome_id=None):
    records=[]; path=Path(path); source_checksum=sha256(path.read_bytes()).hexdigest()
    for line in path.read_text().splitlines():
        if not line or line.startswith('#'): continue
        f=line.split('\t')
        if len(f)<9 or f[2].lower() not in {'cds','gene'}: continue
        attrs=dict((x.split('=',1) for x in f[8].split(';') if '=' in x)); start,end=int(f[3]),int(f[4]); seq=attrs.get('translation','')
        records.append({'genome_id':genome_id or f[0],'tool':tool,'tool_version':None,'locus_id':attrs.get('ID') or attrs.get('locus_tag') or f"{tool}_{len(records)+1}",'start':start,'end':end,'strand':f[6],'protein_length':len(seq) if seq else None,'protein_sequence_checksum':sha256(seq.encode()).hexdigest() if seq else None,'product':attrs.get('product'),'gene':attrs.get('gene'),'source_annotation':attrs,'original_identifier':attrs.get('ID'),'source_path':str(path),'source_sha256':source_checksum})
    return records

def import_prokka(path, genome_id=None): return import_gff(path, 'Prokka', genome_id)
def import_pharokka(path, genome_id=None): return import_gff(path, 'Pharokka', genome_id)
def import_multiphate(path, genome_id=None): return import_gff(path, 'multiPhATE2', genome_id)

def import_genbank(path, tool, genome_id=None):
    """Import CDS coordinates and products from a GenBank flat file."""
    path=Path(path); source_checksum=sha256(path.read_bytes()).hexdigest(); records=[]; current=None; qualifier=None
    for line in path.read_text().splitlines():
        match=re.match(r"^\s{5}CDS\s+(complement\()?<?(\d+)\.\.>?(\d+)\)?", line)
        if match:
            if current: records.append(current)
            current={'genome_id':genome_id,'tool':tool,'tool_version':None,'locus_id':None,
                     'start':int(match.group(2)),'end':int(match.group(3)),'strand':'-' if match.group(1) else '+',
                     'protein_length':None,'protein_sequence_checksum':None,'product':None,'gene':None,
                     'source_annotation':{},'original_identifier':None}
            qualifier=None; continue
        if current is None: continue
        q=re.match(r'^\s+/([A-Za-z_]+)="?(.*)$', line)
        if q:
            qualifier=q.group(1); current['source_annotation'][qualifier]=q.group(2).rstrip('"')
        elif qualifier and re.match(r'^\s{21}\S', line):
            current['source_annotation'][qualifier] += line.strip().rstrip('"')
        else: qualifier=None
    if current: records.append(current)
    for index, record in enumerate(records,1):
        attrs=record['source_annotation']; ident=attrs.get('locus_tag') or attrs.get('protein_id') or f"{tool}_{index}"
        sequence=re.sub(r'\s+','',attrs.get('translation',''))
        record.update({'genome_id':record['genome_id'] or 'unknown','locus_id':ident,'original_identifier':ident,
                       'product':attrs.get('product'),'gene':attrs.get('gene'),
                       'protein_length':len(sequence) if sequence else None,
                       'protein_sequence_checksum':sha256(sequence.encode()).hexdigest() if sequence else None,
                       'source_path':str(path),'source_sha256':source_checksum})
    return records

def import_phold(path, genome_id=None): return import_genbank(path, 'Phold', genome_id)

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
        product=a.get('proposed_function') or a.get('annotation') or c.get('display_product') or c.get('proposed_function')
        record.update({'product':product,'functional_state':c.get('functional_state') if classification_available else None,'proposed_function':c.get('proposed_function') if classification_available else None,'confidence':c.get('confidence') if classification_available else None,'functional_classification_status':'AVAILABLE' if classification_available else 'LEGACY_OUTPUT_NOT_AVAILABLE','source_annotation':{'annotation':a,'classification':c if classification_available else None}})
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
        hypothetical=sum(not x.get('product') or any(term in x.get('product','').lower() for term in ('hypothetical','unknown function','uncharacterized','uncharacterised')) for x in methods[a])
        summary.append({'method':a,'cds_count':len(methods[a]),'specific_annotations':len(methods[a])-hypothetical,'hypothetical':hypothetical})
        if a in references: metric_rows.append({'method':a,**metrics(methods[a],references[a])})
        for b in names[i+1:]: comparisons.extend([{**r,'method_a_name':a,'method_b_name':b} for r in compare_models(methods[a],methods[b])])
    for name,data in [('benchmark_summary',summary),('gene_model_comparison',comparisons),('benchmark_metrics',metric_rows)]:
        (output/(name+'.json')).write_text(json.dumps(data,indent=2,sort_keys=True,default=str))
        cols=sorted({k for row in data for k in row}) if data else ['method'];
        with (output/(name+'.tsv')).open('w',newline='') as h:
            w=csv.DictWriter(h,fieldnames=cols,delimiter='\t'); w.writeheader(); w.writerows(row for row in data)
    (output/'hypothetical_resolution.tsv').write_text('method\tcds_count\thypothetical\n'+'\n'.join(f"{x['method']}\t{x['cds_count']}\t{x['hypothetical']}" for x in summary)+'\n')
    hallmark_rows=[]
    for method,records in methods.items():
        for hallmark,keywords in HALLMARKS.items():
            matches=[record for record in records if any(keyword in str(record.get('product') or '').lower() for keyword in keywords)]
            hallmark_rows.append({'method':method,'hallmark':hallmark,'status':'DETECTED' if matches else 'NOT_ESTABLISHED','locus_ids':';'.join(record['locus_id'] for record in matches),'products':'; '.join(str(record.get('product') or '') for record in matches),'interpretation':'Name-based comparison only; NOT_ESTABLISHED is not biological absence.'})
    hallmark_columns=['method','hallmark','status','locus_ids','products','interpretation']
    with (output/'hallmark_comparison.tsv').open('w',newline='') as handle:
        writer=csv.DictWriter(handle,fieldnames=hallmark_columns,delimiter='\t'); writer.writeheader(); writer.writerows(hallmark_rows)
    (output/'hallmark_comparison.json').write_text(json.dumps(hallmark_rows,indent=2,sort_keys=True))
    functional=[]
    def norm(value):
        value=re.sub(r'\b(putative|probable|predicted)\b','',str(value or '').lower())
        return re.sub(r'[^a-z0-9]+',' ',value).strip()
    for i,a in enumerate(names):
        for b in names[i+1:]:
            for relation in compare_models(methods[a],methods[b]):
                if relation['relationship']!='EXACT_MATCH': continue
                left=relation['method_a'].get('product'); right=relation['method_b'].get('product')
                status='ONE_OR_BOTH_UNANNOTATED' if not left or not right else ('PRODUCT_AGREEMENT' if norm(left)==norm(right) else 'PRODUCT_DIFFERENCE_REQUIRES_REVIEW')
                functional.append({'method_a':a,'method_b':b,'method_a_locus':relation['method_a']['locus_id'],'method_b_locus':relation['method_b']['locus_id'],'coordinates':f"{relation['method_a']['start']}..{relation['method_a']['end']}",'method_a_product':left,'method_b_product':right,'status':status})
    (output/'functional_comparison.json').write_text(json.dumps(functional,indent=2,sort_keys=True))
    functional_columns=['method_a','method_b','method_a_locus','method_b_locus','coordinates','method_a_product','method_b_product','status']
    with (output/'functional_comparison.tsv').open('w',newline='') as handle:
        writer=csv.DictWriter(handle,fieldnames=functional_columns,delimiter='\t'); writer.writeheader(); writer.writerows(functional)
    legacy=[name for name,data in methods.items() if name == 'PHAGEMINE' and any(x.get('functional_classification_status') == 'LEGACY_OUTPUT_NOT_AVAILABLE' for x in data)]
    manifest={'methods':{name:{'record_count':len(records),'source_files':sorted({record.get('source_path') for record in records if record.get('source_path')}),'source_sha256':sorted({record.get('source_sha256') for record in records if record.get('source_sha256')})} for name,records in methods.items()},'interpretation':'Independent software agreement is supporting computational evidence, not experimental truth.'}
    (output/'benchmark_manifest.json').write_text(json.dumps(manifest,indent=2,sort_keys=True)+'\n')
    note='\n\nLegacy PhageMine output detected: gene-model comparison is supported, but modern functional-classification metrics are NOT_EVALUABLE.\n' if legacy else ''
    (output/'benchmark_report.md').write_text('# PhageMine annotation comparison\n\nThis report compares independent annotation outputs. Agreement is supporting computational evidence, not experimental validation. Gene-model differences and product-name differences are reported separately; no consensus call is silently substituted.'+note+'\n')
    return summary
