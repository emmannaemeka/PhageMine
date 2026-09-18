#!/usr/bin/env python3
"""Create a separate post-unblinding PhageMine error-analysis dataset.

The input must be a genuine completed review and the output is explicitly
development data. It never writes to evaluation/v1.2_functional/.
"""
from __future__ import annotations
import argparse, csv, json
from pathlib import Path
from post_adjudication import check_reviews, ROOT, FROZEN

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--reviews',required=True);ap.add_argument('--output',type=Path,required=True);ap.add_argument('--unblind',action='store_true');a=ap.parse_args()
    if not a.unblind: raise SystemExit('Refusing: pass --unblind only after genuine review import and documented reviewer authorization')
    reviews=check_reviews(a.reviews)
    key={r['case_id']:r for r in csv.DictReader((ROOT/'blinding/tool_blinding_key.tsv').open(),delimiter='\t')}
    cases={r['case_id']:r for r in csv.DictReader((ROOT/'reviewer_files/blinded_adjudication_cases.csv').open())}
    out=[]
    for r in reviews:
        c=cases[r['case_id']]; mapping=key[r['case_id']]
        for side in ('A','B'):
            if mapping['prediction_'+side+'_tool']!='PhageMine': continue
            out.append({'dataset_status':'DEVELOPMENT_ERROR_ANALYSIS_NOT_EVALUATION','case_id':r['case_id'],'genome':c['genome'],'locus':c['reference_locus'],'coordinates':f"{c['start']}..{c['end']}",'strand':c['strand'],'prediction':c['prediction_'+side+'_product'],'reference_product':c['reference_product'],'final_judgment':r['final_judgment'],'error_category':r.get('error_category',''),'evidence_strength':r.get('evidence_strength',''),'conflicting_evidence':r.get('conflicting_evidence',''),'specificity_appropriate':r.get('specificity_appropriate',''),'preferred_product_name':r.get('preferred_product_name',''),'reviewer_rationale':r.get('reviewer_rationale','')})
    a.output.mkdir(parents=True,exist_ok=True)
    with (a.output/'phagemine_error_analysis.tsv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(out[0]) if out else ['dataset_status']);w.writeheader();w.writerows(out)
    (a.output/'README.md').write_text('# Development error-analysis dataset\n\nThis dataset is derived after genuine blinded review and unblinding. It is development data, not the frozen v1.2 evaluation, and must not be used as its own test set. Separate development, validation and untouched external test genomes are required.\n')
if __name__=='__main__':main()
