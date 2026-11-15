"""Evaluate frozen runs against external fixtures, retaining compound locations."""
from pathlib import Path
import argparse
import collections
import json
import re
import subprocess
import sys
from Bio import SeqIO
from phagemine.cds_benchmark import read_cds, compare

ACCESSIONS=['PZ583544.1','NC_005282.1','NC_000866.4','NC_001416.1','NC_001604.1','NC_002656.1','OP373729.1','PX787804.1']


def read_json(path, default):
    return json.loads(path.read_text()) if path.is_file() else default


def describe_ref(row, refs):
    note=' '.join(' '.join(row['qualifiers'].get(k,[])) for k in ('note','exception','product'))
    overlaps=[]
    for other in refs:
        if other is row:continue
        size=sum(max(0,min(a[1],b[1])-max(a[0],b[0])+1) for a in row['parts'] for b in other['parts'])
        if size:overlaps.append({'id':other['id'],'overlap_nt':size,'same_strand':other['strand']==row['strand']})
    return {**row,'reference_overlaps':overlaps,'compound':len(row['parts'])>1,
            'frameshift_or_slippage':bool('ribosomal_slippage' in row['qualifiers'] or re.search('frameshift|slippage|bypass',note,re.I)),
            'intron_note':bool(re.search('intron|splic',note,re.I))}


def readable_location(row):
    parts=row['parts']
    if row['strand']=='-': parts=list(reversed(parts))
    locations=[f"{p[0]}..{p[1]}" for p in parts]
    location=locations[0] if len(locations)==1 else 'join('+','.join(locations)+')'
    return 'complement('+location+')' if row['strand']=='-' else location


def evaluate(accession,fixtures,runs,output,repo):
    root=runs/accession;gb=root/'genbank_submission/genome.gbk'
    if not gb.exists():return {'accession':accession,'state':'NO_EXPORT','log_status':read_json(runs.parent/'logs'/f'{accession}.status.json',{})}
    refpath=fixtures/(accession+'.gbk');refs=read_cds(refpath);pred=read_cds(gb)
    for row in refs+pred:
        row['biopython_location']=row['location']
        row['location']=readable_location(row)
    # Reserve exact reference structures before assigning shared-stop alternatives.
    structures={tuple(tuple(x) for x in p['parts']) for p in pred}
    ordered=sorted(refs,key=lambda r:tuple(tuple(x) for x in r['parts']) not in structures)
    report=compare(pred,ordered)
    queue=read_json(root/'gene_calls/start_refinement/short_orf_candidates.json',[])+read_json(root/'gene_calls/start_refinement/start_candidates.json',[])
    for row in report['missed']:
        row.update(describe_ref(row,refs))
        row['manual_review_candidates']=[c for c in queue if len(row['parts'])==1 and c['start']==row['start'] and c['end']==row['end'] and c['strand']==row['strand']]
    reference=SeqIO.read(refpath,'genbank');prediction=SeqIO.read(gb,'genbank')
    rna_ref=[f for f in reference.features if f.type=='tRNA'];rna_pred=[f for f in prediction.features if f.type=='tRNA']
    rna_exact=[f for f in rna_ref if any(f.location==p.location for p in rna_pred)]
    rna_overlap=[f for f in rna_ref if any(f.location.strand==p.location.strand and min(int(f.location.end),int(p.location.end))>max(int(f.location.start),int(p.location.start)) for p in rna_pred)]
    spanins=[]
    diagnostics=[]
    matches=report['exact']+report['stop_only']+report['structure_mismatch']
    for row in refs:
        q=row['qualifiers'];text=str(q)
        if re.search(r'\bspanin\b|\bRz1?\b',text,re.I):
            match=next((m for m in matches if m['reference']['id']==row['id']),None)
            spanins.append({'reference':row,'prediction':match['prediction'] if match else None})
        if len(row['parts'])>1 or re.search(r'frameshift|ribosomal_slippage|translationally bypass',text,re.I) or (accession=='NC_001416.1' and q.get('gene') in [['G'],['T']]):
            match=next((m for m in matches if m['reference']['id']==row['id']),None)
            overlaps=[p for p in pred if p['strand']==row['strand'] and min(p['end'],row['end'])>=max(p['start'],row['start'])]
            diagnostics.append({'reference':describe_ref(row,refs),'matched_prediction':match['prediction'] if match else None,'overlapping_predictions':overlaps})
    verifier=output/'verification';verifier.mkdir(exist_ok=True)
    jsonpath=verifier/(accession+'.codons.json')
    cmd=[sys.executable,str(repo/'scripts/verify_codons.py'),str(gb),'--tbl',str(gb.parent/'features.tbl'),'--extended-starts','--require-table','11','--json',str(jsonpath)]
    result=subprocess.run(cmd,capture_output=True,text=True)
    (verifier/(accession+'.codons.txt')).write_text(result.stdout+result.stderr)
    verification=read_json(jsonpath,{})
    checks=verification.get('cds',[])
    qc=read_json(root/'phagemine_qc_report.json',{})
    validation=read_json(gb.parent/'validation.json',{})
    changes=read_json(root/'gene_calls/start_refinement/changes.json',[])
    extras={p['id']:p for p in report['unmatched_predictions']}
    rescued_extras=[]
    for change in changes:
        if 'RESCUED' not in change.get('code','') or change['protein_id'] not in extras: continue
        prediction_row=extras[change['protein_id']]
        rna_hits=[{'location':str(f.location),'product':f.qualifiers.get('product'),'overlap_nt':min(prediction_row['end'],int(f.location.end))-max(prediction_row['start'],int(f.location.start)+1)+1} for f in rna_ref if min(prediction_row['end'],int(f.location.end))>=max(prediction_row['start'],int(f.location.start)+1)]
        rescued_extras.append({'prediction':prediction_row,'rescue_decision':change,'reference_trna_overlap':rna_hits})
    codons={name:dict(collections.Counter(p['start_codon'] for p in rows)) for name,rows in [('reference',refs),('predicted',pred),('exact_reference',[m['reference'] for m in report['exact']]),('stop_only_reference',[m['reference'] for m in report['stop_only']]),('stop_only_predicted',[m['prediction'] for m in report['stop_only']])]}
    return {'accession':accession,'state':'COMPLETE','description':reference.description,'length_nt':len(reference),'gc_percent':round(100*(str(reference.seq).count('G')+str(reference.seq).count('C'))/len(reference),2),
            **report,'stop_match_count':len(matches),'spanins':spanins,'trna':{'reference':len(rna_ref),'predicted':len(rna_pred),'exact':len(rna_exact),'overlap_recovered':len(rna_overlap),'missed':[str(f.location) for f in rna_ref if f not in rna_overlap]},
            'codon_verifier':{'returncode':result.returncode,'passed':sum(not c['problems'] for c in checks),'total':len(checks),'failed':[c for c in checks if c['problems']],'package_problems':verification.get('package_problems',[])},
            'qc_status':qc.get('status','MISSING'),'qc_criteria':qc.get('criteria',[]),'pipeline_frameshift_findings':[f for c in qc.get('criteria',[]) if c['question']==3 for f in c['details'] if f.get('code')=='PUTATIVE_FRAMESHIFT'],
            'complex_feature_diagnostics':diagnostics,'start_codon_breakdown':codons,'start_changes':changes,'rescued_extras':rescued_extras,
            'submission_readiness':validation.get('submission_readiness'),
            'evidence_adapters':read_json(root/'run_manifest.json',{}).get('evidence_adapters',[])}


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--fixtures',type=Path,required=True);ap.add_argument('--runs',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);a=ap.parse_args()
    repo=Path(__file__).resolve().parents[2];a.output.mkdir(parents=True,exist_ok=True)
    results=[evaluate(x,a.fixtures,a.runs,a.output,repo) for x in ACCESSIONS]
    report={'method':'Frozen structural rules; reference features used only for evaluation in these runs. The original six genomes informed this revision; OP373729.1 and PX787804.1 provide additional regression coverage. Exact recovery requires identical ordered CDS parts, not just outer bounds. Exact matches reserved before shared-stop assignments. Rates use all reference CDSs as denominator. Reported locations are 1-based inclusive; original Biopython locations are retained separately. Stop-only matches can include terminal fragments of joined reference genes and do not establish full complex-gene recovery. Unmatched predictions are benchmark extras, not experimentally proven false genes.','fixtures':read_json(a.fixtures/'manifest.json',[]),'runtime_fingerprint':read_json(a.runs.parent/'source_fingerprint.json',{}),'genomes':results}
    (a.output/'multi_phage_benchmark.json').write_text(json.dumps(report,indent=2)+'\n')
    lines=['# Eight-genome PhageMine regression benchmark','',report['method'],'','| Genome | Ref / predicted CDS | Exact | Stop-only | Missed / extra | tRNA exact/ref (pred) | Codons pass | QC |','|---|---|---|---|---|---|---|---|']
    for r in results:
        if r['state']!='COMPLETE':lines.append(f"| {r['accession']} | {r['state']} | — | — | — | — | — | — |");continue
        t=r['trna'];v=r['codon_verifier'];lines.append(f"| {r['accession']} | {r['reference_count']} / {r['prediction_count']} | {r['exact_count']} ({r['exact_rate']:.1%}) | {r['stop_only_count']} ({r['stop_only_rate']:.1%}) | {len(r['missed'])} / {len(r['unmatched_predictions'])} | {t['exact']}/{t['reference']} ({t['predicted']}) | {v['passed']}/{v['total']} | {r['qc_status']} |")
    for r in results:
        lines+=['',f"## {r['accession']}",'']
        if r['state']!='COMPLETE':lines+=[json.dumps(r)];continue
        lines += [f"GC: {r['gc_percent']}%. INSDC readiness: {r['submission_readiness']}.",f"Start codons: `{json.dumps(r['start_codon_breakdown'])}`",f"14-question statuses: {', '.join(str(c['question'])+'='+c['status'] for c in r['qc_criteria'])}.",'','### Stop-only examples','']
        for m in r['stop_only'][:10]:
            p,q=m['prediction'],m['reference'];lines.append(f"- {q['id']}: reference {q['start']}..{q['end']} ({q['strand']}, {q['start_codon']}); predicted {p['id']} {p['start']}..{p['end']} ({p['start_codon']}).")
        lines+=['','### Missed reference CDSs','']
        for q in r['missed']:lines.append(f"- {q['id']} {q['location']}, {q['length_nt']} nt; product={q['qualifiers'].get('product')}; compound={q['compound']}; frameshift/slippage={q['frameshift_or_slippage']}; overlaps={len(q['reference_overlaps'])}; candidate-review entries={len(q['manual_review_candidates'])}.")
        lines+=['','### Extra predictions','']
        for p in r['unmatched_predictions']:lines.append(f"- {p['id']} {p['location']}, {p['length_nt']} nt, {p['start_codon']}/{p['stop_codon']}; {p['qualifiers'].get('product')}.")
        lines+=['','### New rescue calls absent from reference','']
        for rescue in r['rescued_extras']:
            p=rescue['prediction'];scores=rescue['rescue_decision']['score_components']
            lines.append(f"- {p['id']} {p['location']}: RBS={scores['rbs_score']}, homology bonus={scores['homology_score']}, upstream gap={scores['upstream_gap_nt']} nt; reference tRNA overlaps={rescue['reference_trna_overlap']}.")
        lines+=['','### Spanins','']
        for s in r['spanins']:
            q,p=s['reference'],s['prediction'];lines.append(f"- Reference {q['id']} {q['location']} gene={q['qualifiers'].get('gene')} {q['qualifiers'].get('product')}; prediction: {p['id']+' '+p['location']+' '+str(p['qualifiers'].get('product')) if p else 'MISSED'}.")
        lines+=['','### Complex CDS and frameshift diagnostics','',f"Pipeline split-homology frameshift warnings: {len(r['pipeline_frameshift_findings'])}."]
        for d in r['complex_feature_diagnostics']:
            q=d['reference'];lines.append(f"- {q['id']} {q['location']}; gene={q['qualifiers'].get('gene')}; note={q['qualifiers'].get('note')}; matched={d['matched_prediction']['id'] if d['matched_prediction'] else 'none'}; overlapping calls={[p['id'] for p in d['overlapping_predictions']]}. Joined structures are not inferred from an ordinary overlapping ORF.")
        if r['codon_verifier']['failed']:
            lines+=['','### Strict codon failures (verbatim)','```json',json.dumps(r['codon_verifier']['failed'],indent=2),'```']
    lines+=['','## Interpretation limits','','The supplied Lambda RefSeq record annotates G and T separately and has no compound CDS; it does not supply a joined gpGT feature to score directly. T7, T4 and Bxb1 do contain compound CDSs. Reference notes and genes are retained, because some RefSeq spanin products say Rz-like even when the gene/note identifies Rz1. A zero frameshift-warning count is not proof that programmed recoding is absent. No caller weights, micro-ORF deletion thresholds or accession-specific runtime rules were tuned during this benchmark.']
    (a.output/'multi_phage_benchmark.md').write_text('\n'.join(lines)+'\n')
    print('\n'.join(lines[:12]))

if __name__=='__main__':main()
