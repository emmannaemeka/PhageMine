"""Offline, opt-in Protein Family (PMF) foundation."""
from __future__ import annotations
import csv, hashlib, json, shutil, subprocess, tempfile
from datetime import datetime, timezone
from pathlib import Path

UNRESOLVED_LABELS = ('hypothetical protein','unknown protein','uncharacterized protein',
                     'uncharacterised protein','protein of unknown function')

def is_unresolved_label(label):
    if not label: return True
    text=str(label).strip().casefold()
    return any(text == x or text.startswith(x+' ') for x in UNRESOLVED_LABELS)

def _fasta(path):
    out={}; ident=None; seq=[]
    for line in Path(path).read_text().splitlines():
        if line.startswith('>'):
            if ident is not None: out[ident]=''.join(seq)
            ident=line[1:].split()[0]; seq=[]
        else: seq.append(line.strip())
    if ident is not None: out[ident]=''.join(seq)
    return out

def _next_id(registry):
    nums=[int(x['family_id'].split('-')[1]) for x in registry.get('families',[]) if str(x.get('family_id','')).startswith('PMF-') and x['family_id'].split('-')[1].isdigit()]
    return f'PMF-{(max(nums or [0])+1):06d}'

def _members(protein_fasta, metadata=None):
    metadata=metadata or {}; out=[]
    for pid,seq in sorted(_fasta(protein_fasta).items()):
        m=metadata.get(pid,{})
        out.append({'member_id':pid,'sequence_sha256':hashlib.sha256(seq.encode()).hexdigest(),
                    'source_genome_id':m.get('genome_id') or m.get('source_genome_id') or 'UNKNOWN',
                    'source_protein_id':m.get('protein_id',pid),'host_genus':m.get('host_genus'),
                    'original_annotation':m.get('annotation'),'sequence':seq})
    return out

class MMseqsFamilyBackend:
    """Small adapter around MMseqs2; higher layers never construct commands."""
    def __init__(self, executable='mmseqs'): self.executable=executable
    def cluster(self, protein_fasta, minimum_identity=0.3, minimum_coverage=0.5, coverage_mode=0, clustering_mode=0):
        exe=self.executable if Path(self.executable).exists() else shutil.which(self.executable)
        if not exe: raise RuntimeError('MMseqs2 unavailable; provide --mmseqs or install mmseqs')
        with tempfile.TemporaryDirectory(prefix='phagemine-family-') as t:
            root=Path(t); out=root/'cluster'; tmp=root/'tmp'
            cmd=[exe,'easy-cluster',str(Path(protein_fasta).resolve()),str(out),str(tmp),'--min-seq-id',str(minimum_identity),'-c',str(minimum_coverage),'--cov-mode',str(coverage_mode),'--cluster-mode',str(clustering_mode)]
            r=subprocess.run(cmd,capture_output=True,text=True,check=False)
            if r.returncode: raise RuntimeError('MMseqs2 family clustering failed: '+(r.stderr.strip() or r.stdout.strip()))
            pairs=Path(str(out)+'_cluster.tsv')
            groups={}
            for line in pairs.read_text().splitlines():
                rep,member=line.split('\t')[:2]; groups.setdefault(rep,[]).append(member)
            return groups, {'backend':'MMSEQS2','executable':exe,'command':cmd}

def build_database(protein_fasta, output, metadata=None, registry=None, database_version='1.0', backend='EXACT_SEQUENCE', config=None, mmseqs='mmseqs'):
    """Build a deterministic local family DB. Exact grouping is the offline foundation backend."""
    members=_members(protein_fasta,metadata); reg=json.loads(Path(registry).read_text()) if registry and Path(registry).exists() else {'families':[]}
    if isinstance(reg,list): reg={'families':reg}
    old_by_sha={f.get('representative_sequence_sha256'):f for f in reg.get('families',[])}; groups={}
    if backend.upper() == 'MMSEQS2':
        cfg=(config or {}).get('family_build', config or {}); cluster_map, backend_provenance=MMseqsFamilyBackend(mmseqs).cluster(protein_fasta,cfg.get('minimum_identity',0.3),cfg.get('minimum_coverage',0.5),cfg.get('coverage_mode',0),cfg.get('clustering_mode',0))
        byid={m['member_id']:m for m in members}
        for rep, ids in cluster_map.items(): groups[byid[rep]['sequence_sha256']]=[byid[i] for i in ids]
    else:
        backend_provenance={'backend':'EXACT_SEQUENCE'}
        for m in members: groups.setdefault(m['sequence_sha256'],[]).append(m)
    families=[]; used=set()
    for sha,ms in sorted(groups.items()):
        old=old_by_sha.get(sha)
        if not old:
            overlaps=[f for f in reg.get('families',[]) if set(f.get('member_protein_ids',[])) & {m['member_id'] for m in ms}]
            old=overlaps[0] if len(overlaps)==1 else None
        fid=old['family_id'] if old else _next_id({'families':families+reg.get('families',[])})
        used.add(fid); ann=[m['original_annotation'] for m in ms if m.get('original_annotation')]
        families.append({'family_id':fid,'representative_protein_id':ms[0]['member_id'],'representative_sequence':ms[0]['sequence'],
          'representative_sequence_sha256':sha,'member_count':len(ms),'genome_count':len({m['source_genome_id'] for m in ms}),
          'member_protein_ids':[m['member_id'] for m in ms],'source_genome_ids':sorted({m['source_genome_id'] for m in ms}),
          'host_genera':sorted({m['host_genus'] for m in ms if m.get('host_genus')}), 'member_annotations':ann,
          'consensus_functional_status':'UNKNOWN' if not ann or all(is_unresolved_label(a) for a in ann) else 'MIXED',
          'created_at':old.get('created_at') if old else datetime.now(timezone.utc).isoformat(), 'updated_at':datetime.now(timezone.utc).isoformat(),
          'family_database_version':database_version,'members':[{k:v for k,v in m.items() if k!='sequence'} for m in ms]})
    out=Path(output); out.mkdir(parents=True,exist_ok=True)
    (out/'families.json').write_text(json.dumps(families,indent=2,sort_keys=True));
    (out/'representatives.faa').write_text(''.join(f">{f['family_id']}\n{f['representative_sequence']}\n" for f in families))
    with (out/'families.tsv').open('w',newline='') as h:
        cols=['family_id','representative_protein_id','member_count','genome_count','host_genera','consensus_functional_status','family_database_version']; w=csv.DictWriter(h,fieldnames=cols,delimiter='\t');w.writeheader();w.writerows({k:(','.join(f['host_genera']) if k=='host_genera' else f.get(k)) for k in cols} for f in families)
    with (out/'family_members.tsv').open('w',newline='') as h:
        cols=['family_id','member_id','sequence_sha256','source_genome_id','source_protein_id','host_genus','original_annotation'];w=csv.DictWriter(h,fieldnames=cols,delimiter='\t');w.writeheader()
        for f in families:
            for m in f['members']:w.writerow({'family_id':f['family_id'],**m})
    (out/'family_database_manifest.json').write_text(json.dumps({'schema_version':'1.0','database_version':database_version,'backend':backend,'backend_provenance':backend_provenance,'config':config or {'family_build':{}},'families':len(families)},indent=2,sort_keys=True))
    return families

def _family_hits(queries, families, executable, config):
    exe=executable if Path(executable).exists() else shutil.which(executable)
    if not exe: raise RuntimeError('MMseqs2 unavailable; provide --mmseqs or install mmseqs')
    with tempfile.TemporaryDirectory(prefix='phagemine-family-assign-') as t:
        root=Path(t); qf=root/'queries.faa'; rf=root/'representatives.faa'; out=root/'hits'; tmp=root/'tmp'
        qf.write_text(''.join(f">{q['protein_id']}\n{q['sequence']}\n" for q in queries)); rf.write_text(''.join(f">{f['family_id']}\n{f['representative_sequence']}\n" for f in families))
        # easy-search internally creates the temporary databases and is kept
        # entirely inside this backend.
        cmd=[exe,'easy-search',str(qf),str(rf),str(out),str(tmp),'--format-output','query,target,pident,alnlen,qlen,tlen,evalue,bits,qcov,tcov','--max-seqs','50']
        r=subprocess.run(cmd,capture_output=True,text=True,check=False)
        if r.returncode: raise RuntimeError('MMseqs2 family assignment failed: '+(r.stderr.strip() or r.stdout.strip()))
        rows=[]
        for line in out.read_text().splitlines():
            f=line.split('\t')
            if len(f)!=10: continue
            try: rows.append({'query_protein_id':f[0],'family_id':f[1],'identity':float(f[2]),'alignment_length':int(f[3]),'query_length':int(f[4]),'target_length':int(f[5]),'evalue':float(f[6]),'bitscore':float(f[7]),'query_coverage':float(f[8]),'target_coverage':float(f[9]),'raw_mmseqs_row':line})
            except ValueError: continue
        return rows, {'search_backend':'MMSEQS2','search_backend_version':subprocess.run([exe,'version'],capture_output=True,text=True,check=False).stdout.strip().splitlines()[0] if subprocess.run([exe,'version'],capture_output=True,text=True,check=False).stdout else 'unknown','command':cmd,'configuration':config}

def assign_protein_families(proteins, family_database, config=None, mmseqs='mmseqs'):
    config=config or {'family_assignment':{'minimum_identity':0.3,'minimum_query_coverage':0.5,'minimum_target_coverage':0.5,'maximum_evalue':1e-3,'ambiguity_margin':0.02}}
    cfg=config.get('family_assignment',config); db=Path(family_database); files=db/'families.json' if db.is_dir() else db
    if not files.exists(): return [{'protein_id':q['protein_id'],'assignment_status':'FAMILY_DATABASE_UNAVAILABLE','candidate_families':[]} for q in proteins]
    families=json.loads(files.read_text())
    try: hits, provenance=_family_hits(proteins,families,mmseqs,cfg)
    except RuntimeError:
        return [{'protein_id':q['protein_id'],'assignment_status':'FAMILY_DATABASE_UNAVAILABLE','candidate_families':[]} for q in proteins]
    byq={q['protein_id']:[] for q in proteins}
    for h in hits:
        if h['identity']>=cfg.get('minimum_identity',0.3) and h['query_coverage']>=cfg.get('minimum_query_coverage',0.5) and h['target_coverage']>=cfg.get('minimum_target_coverage',0.5) and h['evalue']<=cfg.get('maximum_evalue',1e-3): byq.setdefault(h['query_protein_id'],[]).append(h)
    result=[]
    for q in sorted(proteins,key=lambda x:x['protein_id']):
        hs=sorted(byq.get(q['protein_id'],[]),key=lambda h:(-h['bitscore'],h['family_id']))
        base={'protein_id':q['protein_id'],'query_sequence_sha256':hashlib.sha256(q['sequence'].encode()).hexdigest(),'family_database_version':families[0].get('family_database_version') if families else None,'provenance':provenance}
        if not hs: result.append({**base,'assignment_status':'NO_FAMILY_MATCH','candidate_families':[]}); continue
        best=hs[0]; competing=[h for h in hs if best['bitscore']-h['bitscore'] <= cfg.get('ambiguity_margin',0.02)*max(abs(best['bitscore']),1)]
        result.append({**base,**best,'candidate_families':[h['family_id'] for h in competing],'family_id':best['family_id'] if len(competing)==1 else None,'assignment_status':'CONFIDENT_FAMILY_MATCH' if len(competing)==1 else 'AMBIGUOUS_FAMILY_MATCH'})
    return result

def assign_protein_family(sequence, family_database, config=None, mmseqs='mmseqs'):
    return assign_protein_families([{'protein_id':'QUERY','sequence':sequence}],family_database,config,mmseqs)[0]

def write_assignments(assignments, output):
    out=Path(output);out.mkdir(parents=True,exist_ok=True);(out/'protein_family_assignments.json').write_text(json.dumps(assignments,indent=2,sort_keys=True))
    cols=['protein_id','family_id','assignment_status','candidate_families','identity','query_coverage','target_coverage','family_database_version']
    with (out/'protein_family_assignments.tsv').open('w',newline='') as h:
        w=csv.DictWriter(h,fieldnames=cols,delimiter='\t'); w.writeheader()
        for a in assignments:
            row={k:a.get(k) for k in cols}; row['protein_id']=a.get('protein_id') or a.get('query_protein_id'); row['candidate_families']=','.join(a.get('candidate_families',[])); w.writerow(row)
