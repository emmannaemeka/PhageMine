#!/usr/bin/env python3
"""Non-scoring review QC; deliberately does not read the confidential key."""
from __future__ import annotations
import argparse,csv,json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
JUDGMENTS={'CORRECT','PARTIALLY_CORRECT','TOO_GENERAL','UNSUPPORTED_SPECIFIC','INCORRECT','UNRESOLVABLE','NOT_EVALUABLE'}
FIELDS=('final_judgment','evidence_supports_prediction','evidence_strength','specificity_appropriate','conflicting_evidence','error_category','reviewer_id','review_timestamp')

def report(path):
    with Path(path).open(newline='') as f: rows=list(csv.DictReader(f))
    expected={f'FB{i:04d}' for i in range(1,690)}; ids=[r.get('case_id','') for r in rows]
    duplicate=sorted({x for x in ids if ids.count(x)>1}); missing=sorted(expected-set(ids)); invalid=[]; missing_fields=0; completed=0; counts={}
    for r in rows:
        if r.get('final_judgment'):
            if r['final_judgment'] not in JUDGMENTS: invalid.append({'case_id':r.get('case_id'),'field':'final_judgment','value':r['final_judgment']})
            else: counts[r['final_judgment']]=counts.get(r['final_judgment'],0)+1
        if all(r.get(f,'').strip() for f in FIELDS): completed+=1
        missing_fields += sum(1 for f in FIELDS if not r.get(f,'').strip())
    return {'total_blinded_units':len(rows),'completed_rows':completed,'incomplete_rows':len(rows)-completed,'judgment_counts':counts,'missing_required_field_cells':missing_fields,'invalid_categories':invalid,'duplicate_ids':duplicate,'missing_ids':missing,'key_accessed':False,'scoring_performed':False,'status':'ADJUDICATION_READY_FOR_FREEZE' if len(rows)==689 and not duplicate and not missing and not invalid and completed==689 else 'INCOMPLETE_REVIEW'}

def main():
    ap=argparse.ArgumentParser();ap.add_argument('csv');a=ap.parse_args();print(json.dumps(report(a.csv),indent=2,sort_keys=True))
if __name__=='__main__':main()
