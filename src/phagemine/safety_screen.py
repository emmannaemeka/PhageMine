"""Optional similarity screens against user-configured, category-specific databases.

These are review candidates, not clinical safety determinations or replacements
for database-specific curated models such as AMRFinderPlus/CARD.
"""
from pathlib import Path
import hashlib
import json
import shutil
import subprocess

CATEGORIES={'lysogeny','toxin_virulence','amr'}


def parse_hits(text, category, identity=40.0, coverage=0.8):
    findings=[]
    for line in text.splitlines():
        fields=line.split('\t')
        if len(fields)!=10: raise ValueError('Malformed safety DIAMOND output')
        query,target=fields[:2]
        pid,qs,qe,ql,ss,se,sl,ev=map(float,fields[2:])
        if not ql>0 or not sl>0: raise ValueError('Invalid safety alignment lengths')
        qcov=(abs(qe-qs)+1)/ql;scov=(abs(se-ss)+1)/sl
        if pid>=identity and min(qcov,scov)>=coverage and ev<=1e-5:
            findings.append({'category':category,'id':query.split('|')[-1],'target':target,'percent_identity':pid,'query_coverage':qcov,'subject_coverage':scov,'evalue':ev,'interpretation':'homology candidate - review'})
    return findings


def screen(output, databases, executable='diamond', threads=1):
    root=Path(output);package=root/'genbank_submission'
    if not package.exists():package=root
    if threads<1 or not databases or set(databases)-CATEGORIES: raise ValueError('Supply known safety categories and positive threads')
    exe=shutil.which(executable)
    if not exe:raise ValueError('DIAMOND is unavailable')
    query=package/'proteins.faa'
    if not query.is_file():raise ValueError('Final proteins.faa is required')
    directory=root/'safety_screen';directory.mkdir(exist_ok=True)
    result={'status':'COMPLETE','therapeutic_safety':'NOT_ESTABLISHED','screens':{},'findings':[], 'query_sha256':hashlib.sha256(query.read_bytes()).hexdigest()}
    for category,database in databases.items():
        db=Path(database).resolve()
        if not db.is_file():raise ValueError(f'Safety database not found: {db}')
        digest=hashlib.sha256()
        with db.open('rb') as handle:
            for chunk in iter(lambda:handle.read(1024*1024),b''):digest.update(chunk)
        out=directory/f'{category}.tsv'
        command=[exe,'blastp','--db',str(db),'--query',str(query),'--out',str(out),'--threads',str(threads),'--max-target-seqs','10','--outfmt','6','qseqid','sseqid','pident','qstart','qend','qlen','sstart','send','slen','evalue']
        done=subprocess.run(command,capture_output=True,text=True,timeout=3600)
        state='COMPLETE' if done.returncode==0 and out.exists() else 'FAILED'
        result['screens'][category]={'status':state,'command':command,'database_sha256':digest.hexdigest(),'stderr':done.stderr,'thresholds':{'identity_percent':40,'both_coverage':0.8,'evalue':1e-5}}
        if state=='FAILED':result['status']='FAILED'
        else:result['findings'].extend(parse_hits(out.read_text(),category))
    result['unconfigured_categories']=sorted(CATEGORIES-set(databases))
    (root/'safety_screen.json').write_text(json.dumps(result,indent=2)+'\n')
    return result
