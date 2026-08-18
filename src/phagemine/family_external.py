from __future__ import annotations
import csv,json,shutil,subprocess,tempfile
import hashlib
from datetime import datetime, timezone
from pathlib import Path

def _sha(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def load_pmfdb(path):
    root=Path(path); required=['reference_phage_proteins.faa','reference_metadata.tsv','reference_manifest.json','reference_qc.tsv']
    if not root.is_dir(): raise ValueError('PMFDB path is not a directory')
    missing=[x for x in required if not (root/x).is_file()]
    if missing: raise ValueError('PMFDB missing required files: '+', '.join(missing))
    manifest=json.loads((root/'reference_manifest.json').read_text()); fields=['pmfdb_version','creation_date','source_database','source_release','retrieval_date','filters','genome_count','protein_count','checksums','schema_version']
    absent=[x for x in fields if x not in manifest]
    if absent: raise ValueError('PMFDB manifest missing required fields: '+', '.join(absent))
    if str(manifest['schema_version']) != '1.0': raise ValueError('Unsupported PMFDB schema_version: '+str(manifest['schema_version']))
    for name in required[:1]+required[1:2]+required[3:4]:
        expected=manifest['checksums'].get(name) or manifest['checksums'].get(str(root/name))
        if expected and expected != _sha(root/name): raise ValueError(f'PMFDB checksum mismatch: {name}')
    return root,manifest

def validate_external(database, reference=None, output=None, metadata=None, mmseqs='mmseqs', top_n=None, reference_db=None):
    pmfdb_manifest=None
    if reference_db:
        root,pmfdb_manifest=load_pmfdb(reference_db); reference=root/'reference_phage_proteins.faa'; metadata=root/'reference_metadata.tsv'
        if Path(output).resolve() == root.resolve() or root.resolve() in Path(output).resolve().parents: raise ValueError('Validation output must not be inside PMFDB release')
    elif not reference:
        raise ValueError('Provide --reference-db or --reference-proteins')
    db=Path(database); families=json.loads((db/'families.json').read_text());
    if top_n: families=families[:int(top_n)]
    meta={}
    if metadata:
        mp=Path(metadata)
        if mp.suffix.lower()=='.json': meta=json.loads(mp.read_text())
        else:
            with mp.open() as h:
                for r in csv.DictReader(h,delimiter='\t'):
                    identifier = r.get('external_protein_id') or r.get('protein_id') or r.get('raw_fasta_id') or r.get('id')
                    if identifier: meta[identifier]=r
    exe=mmseqs if Path(mmseqs).exists() else shutil.which(mmseqs)
    if not exe: raise RuntimeError('MMseqs2 unavailable; provide --mmseqs')
    with tempfile.TemporaryDirectory(prefix='phagemine-external-') as t:
        root=Path(t); q=root/'q.faa'; ref=root/'ref.faa'; out=root/'hits'; tmp=root/'tmp'
        q.write_text(''.join(f">{f['family_id']}\n{f['representative_sequence']}\n" for f in families))
        cached_target = None
        cache_manifest = None
        if pmfdb_manifest:
            candidate = Path(reference).parent.parent.parent / 'cache' / 'mmseqs' / Path(reference).parent.name
            manifest_path = candidate / 'mmseqs_index_manifest.json'
            if (candidate/'target_db.dbtype').is_file() and manifest_path.is_file():
                cache_manifest=json.loads(manifest_path.read_text())
                if (cache_manifest.get('pmfdb_version') == pmfdb_manifest.get('pmfdb_version') and
                    cache_manifest.get('pmfdb_reference_fasta_sha256') == _sha(reference) and
                    cache_manifest.get('index_status') == 'SUCCESS'):
                    cached_target=candidate/'target_db'
        if cached_target:
            query_db=root/'query_db'; result_db=root/'result_db'
            commands=[
                [exe,'createdb',str(q),str(query_db)],
                [exe,'search',str(query_db),str(cached_target),str(result_db),str(tmp),'--max-seqs','100'],
                [exe,'convertalis',str(query_db),str(cached_target),str(result_db),str(out),'--format-output','query,target,pident,alnlen,qlen,tlen,evalue,bits,qcov,tcov'],
            ]
            for cmd in commands:
                run=subprocess.run(cmd,capture_output=True,text=True,check=False)
                if run.returncode: raise RuntimeError(run.stderr.strip() or 'MMseqs2 search failed')
            cmd=commands
        else:
            shutil.copyfile(reference,ref)
            cmd=[exe,'easy-search',str(q),str(ref),str(out),str(tmp),'--format-output','query,target,pident,alnlen,qlen,tlen,evalue,bits,qcov,tcov','--max-seqs','100']
            run=subprocess.run(cmd,capture_output=True,text=True,check=False)
            if run.returncode: raise RuntimeError(run.stderr.strip() or 'MMseqs2 search failed')
        hits=[]
        for line in out.read_text().splitlines():
            f=line.split('\t')
            if len(f)!=10: continue
            hits.append({'family_id':f[0],'external_protein_id':f[1],'identity':float(f[2]),'alignment_length':int(f[3]),'query_length':int(f[4]),'target_length':int(f[5]),'evalue':float(f[6]),'bitscore':float(f[7]),'query_coverage':float(f[8]),'target_coverage':float(f[9]),'raw_mmseqs_row':line})
    grouped={f['family_id']:[] for f in families}
    for h in hits: grouped.setdefault(h['family_id'],[]).append(h)
    summaries=[]
    for f in families:
        hs=grouped[f['family_id']]; enriched=[]; genomes=set(); hosts=set(); species=set(); taxa=set(); chars=0; incomplete=0
        for h in hs:
            m=meta.get(h['external_protein_id'],{}) or {}; ann=m.get('annotation') or m.get('product') or m.get('existing_annotation'); st=str(m.get('annotation_status','')).lower(); characterized=str(m.get('characterized','')).lower() in {'1','true','yes'} or st in {'characterized','known','annotated'}; incomplete += not bool(m) or (not ann and not characterized and st not in {'unknown','uncharacterized','hypothetical'})
            h.update({'external_genome_id':m.get('source_genome_id') or m.get('genome_id') or m.get('genome'),'host_genus':m.get('host_genus'),'host_species':m.get('host_species'),'taxon':m.get('phage_taxonomy') or m.get('taxon'),'existing_annotation':ann,'characterized':characterized,'metadata_incomplete':bool(not m),'provenance':{'search_command':cmd,'mmseqs_executed':True}}); enriched.append(h); genomes.add(h.get('external_genome_id')); hosts.add(h.get('host_genus')); species.add(h.get('host_species')); taxa.add(h.get('taxon')); chars+=bool(h['characterized'])
        grouped[f['family_id']]=enriched
        status='NO_EXTERNAL_MATCH' if not hs else ('REFERENCE_METADATA_INCOMPLETE' if incomplete else ('EXTERNAL_CHARACTERIZED_HOMOLOG_FOUND' if chars==len(hs) else ('EXTERNAL_MATCHES_MIXED' if chars else 'EXTERNAL_MATCHES_UNCHARACTERIZED')))
        best=max(hs,key=lambda x:(x['bitscore'],-x['evalue'])) if hs else {}
        summaries.append({'family_id':f['family_id'],'query_representative':f['representative_protein_id'],'external_match_count':len(hs),'distinct_external_genomes':len({x for x in genomes if x}),'distinct_host_genera':len({x for x in hosts if x}),'distinct_host_species':len({x for x in species if x}),'distinct_taxonomic_groups':len({x for x in taxa if x}),'host_genera':','.join(sorted(x for x in hosts if x)),'best_identity':best.get('identity'),'best_query_coverage':best.get('query_coverage'),'best_target_coverage':best.get('target_coverage'),'best_evalue':best.get('evalue'),'best_bitscore':best.get('bitscore'),'characterized_homolog_count':chars,'uncharacterized_homolog_count':len(hs)-chars,'metadata_incomplete_count':incomplete,'external_validation_status':status})
    outdir=Path(output);outdir.mkdir(parents=True,exist_ok=True); cols=list(summaries[0]) if summaries else ['family_id']
    with (outdir/'pmf_external_validation.tsv').open('w',newline='') as h:w=csv.DictWriter(h,fieldnames=cols,delimiter='\t');w.writeheader();w.writerows(summaries)
    prov={'command':cmd,'reference':str(reference),'reference_type':'PMFDB' if pmfdb_manifest else 'CUSTOM_REFERENCE','validation_timestamp':datetime.now(timezone.utc).isoformat(),'prebuilt_target_database':str(cached_target) if cached_target else None,'prebuilt_target_manifest':cache_manifest}
    if pmfdb_manifest:
        prov.update({'phagemine_version':'0.1.0','pmfdb_version':pmfdb_manifest['pmfdb_version'],'pmfdb_schema_version':pmfdb_manifest['schema_version'],'pmfdb_manifest_sha256':_sha(Path(reference).parent/'reference_manifest.json'),'reference_protein_sha256':_sha(reference),'reference_metadata_sha256':_sha(metadata),'reference_qc_sha256':_sha(Path(reference).parent/'reference_qc.tsv')})
    (outdir/'pmf_external_validation.json').write_text(json.dumps({'provenance':prov,'families':summaries},indent=2,sort_keys=True)); (outdir/'pmf_external_validation_manifest.json').write_text(json.dumps(prov,indent=2,sort_keys=True))
    with (outdir/'pmf_external_hits.tsv').open('w',newline='') as h:
        hc=['family_id','external_protein_id','identity','query_coverage','target_coverage','evalue','bitscore','external_genome_id','host_genus','taxon','characterized'];w=csv.DictWriter(h,fieldnames=hc,delimiter='\t');w.writeheader();w.writerows({k:x.get(k) for k in hc} for hs in grouped.values() for x in hs)
    with (outdir/'pmf_external_host_distribution.tsv').open('w',newline='') as h:
        w=csv.writer(h,delimiter='\t'); w.writerow(['family_id','host_genus','external_protein_count','external_genome_count'])
        for fid,hs in grouped.items():
            for host in sorted({x.get('host_genus') for x in hs if x.get('host_genus')}): w.writerow([fid,host,sum(x.get('host_genus')==host for x in hs),len({x.get('external_genome_id') for x in hs if x.get('host_genus')==host and x.get('external_genome_id')})])
    with (outdir/'pmf_external_characterized_hits.tsv').open('w',newline='') as h:
        cols=['family_id','external_protein_id','external_genome_id','existing_annotation','identity','query_coverage','target_coverage','evalue','bitscore']; w=csv.DictWriter(h,fieldnames=cols,delimiter='\t'); w.writeheader(); w.writerows({k:x.get(k) for k in cols} for fid,hs in grouped.items() for x in hs if x.get('characterized'))
    return summaries
