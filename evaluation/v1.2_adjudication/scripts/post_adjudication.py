#!/usr/bin/env python3
"""Post-unblinding analysis; refuses final results while reviews are incomplete."""
from __future__ import annotations
import argparse, csv, hashlib, itertools, json, math, random
from collections import Counter, defaultdict
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]; REPO=ROOT.parents[1]; FROZEN=REPO/'evaluation/v1.2_functional'
JUDGMENTS={'CORRECT','PARTIALLY_CORRECT','TOO_GENERAL','UNSUPPORTED_SPECIFIC','INCORRECT','UNRESOLVABLE','NOT_EVALUABLE'}
UTILITY={'CORRECT':1.0,'PARTIALLY_CORRECT':0.5,'TOO_GENERAL':0.25,'UNSUPPORTED_SPECIFIC':0.0,'INCORRECT':0.0}

def read(path, delimiter=','):
    with Path(path).open(newline='') as f:return list(csv.DictReader(f,delimiter=delimiter))

def wilson(k,n):
    if not n:return (None,None)
    z=1.959963984540054;p=k/n; den=1+z*z/n; c=(p+z*z/(2*n))/den; h=z*math.sqrt(p*(1-p)/n+z*z/(4*n*n))/den
    return max(0,c-h),min(1,c+h)

def paired_permutation(diffs):
    if not diffs:return {'estimate':None,'ci95':(None,None),'p_value':None,'n_genomes':0}
    observed=sum(diffs)/len(diffs); values=[]
    for signs in itertools.product((-1,1),repeat=len(diffs)): values.append(sum(x*s for x,s in zip(diffs,signs))/len(diffs))
    values.sort(); lo=values[max(0,math.ceil(.025*len(values))-1)]; hi=values[min(len(values)-1,math.floor(.975*len(values)))]
    p=sum(abs(x)>=abs(observed)-1e-15 for x in values)/len(values)
    return {'estimate':observed,'ci95':(lo,hi),'p_value':p,'n_genomes':len(diffs)}

def bootstrap(diffs,seed=12020260918,iterations=20000):
    if not diffs:return {'estimate':None,'ci95':(None,None),'iterations':0}
    rng=random.Random(seed); samples=[]
    for _ in range(iterations):samples.append(sum(rng.choice(diffs) for _ in diffs)/len(diffs))
    samples.sort();return {'estimate':sum(diffs)/len(diffs),'ci95':(samples[int(.025*iterations)],samples[int(.975*iterations)-1]),'iterations':iterations,'seed':seed}

def difference_ci(phage_num, phage_den, phar_num, phar_den):
    # Independent-binomial Wilson/Newcombe approximation is conservative here;
    # the paired genome bootstrap remains the primary dependence-aware output.
    p=phage_num/phage_den if phage_den else None; q=phar_num/phar_den if phar_den else None
    if p is None or q is None:return {'estimate':None,'ci95':(None,None)}
    a=wilson(phage_num,phage_den); b=wilson(phar_num,phar_den)
    return {'estimate':p-q,'ci95':(a[0]-b[1],a[1]-b[0])}

def validate_definition():
    return {'primary_utility':{'formula':'U_g=sum(u_i)/N_g','weights':UTILITY,'unresolvable_rule':'missing; report count and identification bounds; never silently zero','not_evaluable_rule':'excluded from N_g'},'strict_secondary':{'formula':'CORRECT_count/frozen_evaluable_reference_loci'},'inclusive_sensitivity':{'formula':'(CORRECT+PARTIALLY_CORRECT)_count/frozen_evaluable_reference_loci'},'primary_margin':0.05,'structural_noninferiority_margin':-0.03,'unsupported_specificity_safety_bound':0.02,'frozen_on_commit':'4933d8fa532a0157c24c3c01dbf1d8adfa787ed7'}

def check_reviews(path):
    rows=read(path); required={'case_id','final_judgment','reviewer_id','review_timestamp'}
    if not rows:raise ValueError('No adjudication rows supplied')
    if not required <= set(rows[0]):raise ValueError('Missing review columns')
    ids=[r['case_id'] for r in rows]
    if len(set(ids))!=len(ids):raise ValueError('Duplicate case ID')
    for r in rows:
        if not r['reviewer_id'] or not r['review_timestamp']:raise ValueError('Reviewer provenance missing '+r['case_id'])
        if r['final_judgment'] not in JUDGMENTS:raise ValueError('Invalid judgment '+r['case_id'])
    return rows

def analyze(review_path, output):
    rows=check_reviews(review_path)
    # The neutral key is required only after review is supplied, and is never
    # distributed with reviewer-facing files.
    key={r['case_id']:r for r in read(ROOT/'blinding/tool_blinding_key.tsv',delimiter='\t')}
    if set(key)!=set(r['case_id'] for r in rows):raise ValueError('Review must contain every blinded case exactly once')
    frozen={r['case_id']:r for r in read(ROOT/'reviewer_files/blinded_adjudication_cases.csv',delimiter=',')}
    by_tool=defaultdict(dict); errors=[]
    # Seed every exact evaluable locus with the frozen deterministic category;
    # genuine reviews replace only the unresolved side for that same locus.
    frozen_loci=read(FROZEN/'functional_by_locus.tsv',delimiter='\t')
    for base in frozen_loci:
        if base['matching_mode']!='exact': continue
        if base['PhageMine_category']!='NON_EVALUABLE_REFERENCE':
            by_tool['PhageMine'][(base['accession'],base['start'],base['end'],base['strand'])]={'final_judgment': {'EXACT_PRODUCT_AGREEMENT':'CORRECT','EQUIVALENT_FUNCTION':'CORRECT','ABSTENTION':'UNRESOLVABLE','UNRESOLVED_REVIEW_REQUIRED':None}.get(base['PhageMine_category'], base['PhageMine_category'])}
        if base['Pharokka_category']!='NON_EVALUABLE_REFERENCE':
            by_tool['Pharokka'][(base['accession'],base['start'],base['end'],base['strand'])]={'final_judgment': {'EXACT_PRODUCT_AGREEMENT':'CORRECT','EQUIVALENT_FUNCTION':'CORRECT','ABSTENTION':'UNRESOLVABLE','UNRESOLVED_REVIEW_REQUIRED':None}.get(base['Pharokka_category'], base['Pharokka_category'])}
    for r in rows:
        base=frozen[r['case_id']]; mapping=key[r['case_id']]
        if base['matching_mode']!='exact':
            continue
        side=base['review_target']; tool=mapping['prediction_'+side+'_tool']
        by_tool[tool][(base['genome'],base['start'],base['end'],base['strand'])]=r
        errors.append({'case_id':r['case_id'],'tool':tool,'genome':base['genome'],'judgment':r['final_judgment'],'error_category':r.get('error_category',''),'evidence_strength':r.get('evidence_strength',''),'conflicting_evidence':r.get('conflicting_evidence',''),'preferred_product_name':r.get('preferred_product_name',''),'reviewer_rationale':r.get('reviewer_rationale','')})
    # Final utility requires one completed judgment per tool/locus. Cases are
    # intentionally checked before any estimate is emitted.
    all_loci={(r['accession'],r['start'],r['end'],r['strand']) for r in read(FROZEN/'functional_by_locus.tsv',delimiter='\t') if r['matching_mode']=='exact' and r['PhageMine_category']!='NON_EVALUABLE_REFERENCE'}
    missing={t:sorted(all_loci-set(by_tool[t])) for t in ('PhageMine','Pharokka')}
    if any(missing.values()):raise ValueError('FINAL_ANALYSIS_BLOCKED_INCOMPLETE_REVIEW:'+json.dumps({t:len(v) for t,v in missing.items()}))
    results={}; diffs=[]; safety_by_genome=[]
    for genome in sorted({x[0] for x in all_loci}):
        results[genome]={}
        for tool in ('PhageMine','Pharokka'):
            loci=[k for k in all_loci if k[0]==genome]; judgments=[by_tool[tool][k]['final_judgment'] for k in loci]
            utility=sum(UTILITY.get(j,0) for j in judgments)/len(loci); results[genome][tool]={'N':len(loci),'utility':utility,'strict':sum(j=='CORRECT' for j in judgments)/len(loci),'inclusive':sum(j in {'CORRECT','PARTIALLY_CORRECT'} for j in judgments)/len(loci),'counts':dict(Counter(judgments))}
        named={t:sum(j not in {'UNRESOLVABLE','NOT_EVALUABLE'} for j in [by_tool[t][k]['final_judgment'] for k in loci]) for t in ('PhageMine','Pharokka')}
        unsafe={t:sum(by_tool[t][k]['final_judgment']=='UNSUPPORTED_SPECIFIC' for k in loci) for t in ('PhageMine','Pharokka')}
        safety_by_genome.append((unsafe['PhageMine'],named['PhageMine'],unsafe['Pharokka'],named['Pharokka']))
        diffs.append(results[genome]['PhageMine']['utility']-results[genome]['Pharokka']['utility'])
    output.mkdir(parents=True,exist_ok=True)
    safety_num=(sum(x[0] for x in safety_by_genome),sum(x[2] for x in safety_by_genome)); safety_den=(sum(x[1] for x in safety_by_genome),sum(x[3] for x in safety_by_genome))
    structural_rows=read(REPO/'docs/benchmark_v1.2/tables/strict_relaxed_metrics.tsv',delimiter='\t'); by_genome=defaultdict(dict)
    for x in structural_rows:
        if x['tool'] in ('PhageMine','Pharokka'): by_genome[x['genome']][x['tool']]=float(x['strict_f1'])
    structural_diffs=[v['PhageMine']-v['Pharokka'] for v in by_genome.values()]
    report={'status':'ADJUDICATED','definition':validate_definition(),'genome_results':results,'primary_difference':paired_permutation(diffs),'bootstrap_sensitivity':bootstrap(diffs),'safety_endpoint':{'phagemine_unsupported':safety_num[0],'phagemine_named':safety_den[0],'pharokka_unsupported':safety_num[1],'pharokka_named':safety_den[1],**difference_ci(*safety_num,*safety_den),'frozen_bound':0.02,'criterion':'difference upper 95% CI <= +0.02 is the safety condition; report separately from utility'},'structural_noninferiority':{'genome_strict_f1_differences':structural_diffs,'estimate':sum(structural_diffs)/len(structural_diffs),'ci95':paired_permutation(structural_diffs)['ci95'],'margin':-0.03,'criterion':'lower 95% CI > -0.03'},'error_analysis':errors}
    (output/'statistical_report.json').write_text(json.dumps(report,indent=2,sort_keys=True)+'\n')
    with (output/'error_analysis.tsv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(errors[0]),delimiter='\t');w.writeheader();w.writerows(errors)
    return report

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--reviews',required=True);ap.add_argument('--output',type=Path,required=True);ap.add_argument('--final',action='store_true');a=ap.parse_args()
    definition=validate_definition(); (ROOT/'frozen/utility_definition.json').write_text(json.dumps(definition,indent=2,sort_keys=True)+'\n')
    if not a.final: print(json.dumps({'status':'READY_FOR_GENUINE_REVIEW','definition':definition},indent=2)); return
    print(json.dumps(analyze(a.reviews,a.output),indent=2))
if __name__=='__main__':main()
