"""Offline comparative genomic context for completed PhageMine results."""
from __future__ import annotations

import csv, json, shutil, subprocess, tempfile
from pathlib import Path
from typing import Any

REQUIRED = ("proteins.faa", "genes.gff3", "functional_classification.json", "genomic_context.json", "modules.json", "evidence.json", "run_manifest.json")
MMSEQS_FIELDS = "query,target,fident,alnlen,qstart,qend,qlen,tstart,tend,tlen,evalue,bits,qcov,tcov"


def _mmseqs_version(executable: str) -> str:
    result = subprocess.run([executable, "version"], capture_output=True, text=True, check=False)
    return (result.stdout or result.stderr).strip().splitlines()[0] or "unknown"


def _mmseqs_matches(query: dict[str, Any], reference: dict[str, Any], executable: str | None):
    if not executable or not (Path(executable).exists() or shutil.which(executable)):
        return {}, {"mmseqs_executed": False, "mmseqs": executable}
    with tempfile.TemporaryDirectory(prefix="phagemine-compare-mmseqs-") as temp:
        root=Path(temp); qf=root/'query.faa'; rf=root/'reference.faa'; qdb=root/'query'; rdb=root/'reference'; resultdb=root/'result'; out=root/'matches.tsv'; tmp=root/'tmp'
        qf.write_text(''.join(f">{p['protein_id']}\n{p.get('sequence','')}\n" for p in query.values())); rf.write_text(''.join(f">{p['protein_id']}\n{p.get('sequence','')}\n" for p in reference.values()))
        commands=[[executable,'createdb',str(qf),str(qdb)],[executable,'createdb',str(rf),str(rdb)],[executable,'search',str(qdb),str(rdb),str(resultdb),str(tmp),'--max-seqs','20'],[executable,'convertalis',str(qdb),str(rdb),str(resultdb),str(out),'--format-output',MMSEQS_FIELDS]]
        for command in commands:
            result=subprocess.run(command,capture_output=True,text=True,check=False)
            if result.returncode:
                return {}, {"mmseqs_executed": False,"mmseqs":executable,"error":result.stderr.strip(),"commands":commands}
        matches={}
        for line in out.read_text().splitlines():
            fields=line.split('\t')
            if len(fields)!=14: continue
            try:
                q,t=fields[0],fields[1]; row={"query_protein_id":q,"target_protein_id":t,"percent_identity":float(fields[2]),"alignment_length":int(fields[3]),"query_start":int(fields[4]),"query_end":int(fields[5]),"target_start":int(fields[7]),"target_end":int(fields[8]),"evalue":float(fields[10]),"bit_score":float(fields[11]),"query_coverage":float(fields[12]),"target_coverage":float(fields[13]),"raw_mmseqs_row":line}
            except ValueError: continue
            if row['evalue'] <= 1e-3 and row['query_coverage'] >= 0.5 and row['target_coverage'] >= 0.5: matches.setdefault(q,[]).append(row)
        return matches,{"mmseqs_executed":True,"mmseqs":executable,"mmseqs_version":_mmseqs_version(executable),"commands":commands,"output_format":MMSEQS_FIELDS,"thresholds":{"evalue":1e-3,"query_coverage":0.5,"target_coverage":0.5}}


def _fasta(path):
    records=[]; current=None; seq=[]
    for line in path.read_text().splitlines():
        if line.startswith('>'):
            if current is not None: records.append((current,''.join(seq)))
            current=line[1:].split()[0]; seq=[]
        else: seq.append(line.strip())
    if current is not None: records.append((current,''.join(seq)))
    return records


def _load(root: Path):
    missing=[str(root/n) for n in REQUIRED if not (root/n).is_file()]
    if missing: raise ValueError("Comparison validation failed; missing artifacts: " + ", ".join(missing))
    proteins={pid:{"protein_id":pid,"sequence":seq} for pid,seq in _fasta(root/'proteins.faa')}
    context={x['protein_id']:x for x in json.loads((root/'genomic_context.json').read_text())}
    classification={x['protein_id']:x for x in json.loads((root/'functional_classification.json').read_text())}
    evidence=json.loads((root/'evidence.json').read_text())
    for item in evidence:
        pid=item['protein_id']; proteins.setdefault(pid,{}) .update({"start":item.get('start'),"end":item.get('end'),"strand":item.get('strand'),"genome_id":item.get('genome_id'),"sequence":item.get('sequence')})
    if set(proteins)!=set(context) or set(proteins)!=set(classification): raise ValueError(f"Comparison validation failed; protein IDs disagree in {root}")
    for pid in proteins:
        proteins[pid].update(context[pid]); proteins[pid].update({"functional_state":classification[pid].get('functional_state'),"proposed_function":classification[pid].get('proposed_function'),"functional_category":classification[pid].get('functional_category'),"conservation_status":classification[pid].get('conservation_status')})
    phrog={}
    for item in evidence:
        for e in item.get('evidence',[]):
            if e.get('supports') and e.get('source')=='PHROGs' and e.get('identifier'): phrog.setdefault(e['identifier'],[]).append(item['protein_id'])
    return root, proteins, phrog


def compare(roots: list[str | Path], output: str | Path, mmseqs: str | None = None):
    if not 2 <= len(roots) <= 10: raise ValueError("Comparison requires 2-10 result directories")
    loaded=[_load(Path(r).resolve()) for r in roots]; output=Path(output).resolve(); output.mkdir(parents=True, exist_ok=True)
    query_root, query, query_phrog=loaded[0]; rows=[]; links=[]
    for qid, qp in query.items():
        per=[]
        for ref_root, ref, ref_phrog in loaded[1:]:
            candidates=[]; method='EXACT_SEQUENCE_MATCH'
            for phrog, qids in query_phrog.items():
                if qid in qids:
                    for rid in ref_phrog.get(phrog,[]): candidates.append(rid)
                    if candidates: method='PHROG_ANCHOR'
            mm_matches, mm_provenance = ({}, {"mmseqs_executed":False})
            if candidates:
                mm_matches, mm_provenance = _mmseqs_matches({qid:qp}, ref, mmseqs or shutil.which('mmseqs'))
                if mm_provenance.get('mmseqs_executed') and any(row['target_protein_id'] in candidates for row in mm_matches.get(qid, [])):
                    method='HYBRID'
            if not candidates:
                executable=mmseqs or shutil.which('mmseqs')
                mm_matches, mm_provenance = _mmseqs_matches({qid:qp}, ref, executable)
                candidates=[row['target_protein_id'] for row in mm_matches.get(qid, [])]
                method='MMSEQS2' if mm_provenance.get('mmseqs_executed') and candidates else 'EXACT_SEQUENCE_MATCH'
                if not candidates: candidates=[rid for rid,rp in ref.items() if rp.get('sequence')==qp.get('sequence')]
            if len(candidates)>1: method='HYBRID' if method=='PHROG_ANCHOR' else 'EXACT_SEQUENCE_MATCH'
            for rid in candidates:
                rp=ref[rid]; mmrow=next((x for x in mm_matches.get(qid,[]) if x['target_protein_id']==rid),None); links.append({"query_protein_id":qid,"reference_genome_id":str(ref_root),"reference_protein_id":rid,"orthology_method":method,"percent_identity":mmrow['percent_identity'] if mmrow else (1.0 if qp.get('sequence')==rp.get('sequence') else None),"sequence_identity":mmrow['percent_identity'] if mmrow else (1.0 if qp.get('sequence')==rp.get('sequence') else None),"alignment_length":mmrow['alignment_length'] if mmrow else (len(qp.get('sequence','')) if qp.get('sequence')==rp.get('sequence') else None),"evalue":mmrow['evalue'] if mmrow else None,"bit_score":mmrow['bit_score'] if mmrow else None,"query_coverage":mmrow['query_coverage'] if mmrow else (1.0 if qp.get('sequence')==rp.get('sequence') else None),"target_coverage":mmrow['target_coverage'] if mmrow else (1.0 if qp.get('sequence')==rp.get('sequence') else None),"query_start":mmrow['query_start'] if mmrow else None,"query_end":mmrow['query_end'] if mmrow else None,"target_start":mmrow['target_start'] if mmrow else None,"target_end":mmrow['target_end'] if mmrow else None,"raw_mmseqs_row":mmrow['raw_mmseqs_row'] if mmrow else None,"provenance":{**mm_provenance,"query_result":str(query_root),"reference_result":str(ref_root)}})
            per.append((ref_root,ref,candidates))
        present=sum(bool(c) for _,_,c in per); order=sum(1 for ref_root,ref,c in per if c and all(ref[r].get('gene_order_index') is not None for r in c))
        orientation=sum(1 for ref_root,ref,c in per if c and any(ref[r].get('strand')==qp.get('strand') for r in c))
        modules=sum(1 for ref_root,ref,c in per if c and qp.get('functional_module') and any(ref[r].get('functional_module')==qp.get('functional_module') for r in c))
        fraction=present/len(per) if per else 0
        status='STRONGLY_CONSERVED_CONTEXT' if fraction==1 and order==len(per) and orientation==len(per) and modules==len(per) else 'CONSERVED_CONTEXT' if fraction>=0.5 and order>=1 else 'PARTIALLY_CONSERVED_CONTEXT' if present else 'NO_CLEAR_CONTEXT'
        rows.append({"query_protein_id":qid,"genomes_compared":len(loaded),"homologs_found":present,"fraction_present":fraction,"conserved_order_count":order,"conserved_orientation_count":orientation,"conserved_module_count":modules,"synteny_status":status,"functional_state":qp.get('functional_state'),"proposed_function":qp.get('proposed_function'),"functional_module":qp.get('functional_module'),"is_cross_genome_conserved_unknown":qp.get('functional_state') in {'CONSERVED_UNKNOWN','UNRESOLVED'} and status in {'STRONGLY_CONSERVED_CONTEXT','CONSERVED_CONTEXT'},"ambiguous_mapping_count":sum(1 for _,_,c in per if len(c)>1),"provenance":{"query_result":str(query_root),"mmseqs":mmseqs or shutil.which('mmseqs')}})
    (output/'comparative_context.json').write_text(json.dumps(rows,indent=2,sort_keys=True)); (output/'ortholog_links.tsv').write_text('query_protein_id\treference_genome_id\treference_protein_id\torthology_method\tsequence_identity\talignment_length\tquery_coverage\ttarget_coverage\n'+''.join(f"{x['query_protein_id']}\t{x['reference_genome_id']}\t{x['reference_protein_id']}\t{x['orthology_method']}\t{x['sequence_identity']}\t{x['alignment_length']}\t{x['query_coverage']}\t{x['target_coverage']}\n" for x in links))
    cols=['query_protein_id','genomes_compared','homologs_found','fraction_present','conserved_order_count','conserved_orientation_count','conserved_module_count','synteny_status','functional_state','proposed_function','functional_module','is_cross_genome_conserved_unknown','ambiguous_mapping_count']
    with (output/'comparative_context.tsv').open('w',newline='') as h: w=csv.DictWriter(h,fieldnames=cols,delimiter='\t'); w.writeheader(); w.writerows({k:x.get(k) for k in cols} for x in rows)
    return rows,links
