#!/usr/bin/env python3
"""Score archived v1.2 outputs without rerunning annotation tools."""
from pathlib import Path
import argparse
import csv
import hashlib
import json
import random
import sys
from urllib.parse import unquote

REPO = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO/'src'))
from phagemine.functional_benchmark import (CATEGORIES, COMPATIBLE, CORRECT, PENDING,
    classify, informative, load_synonyms, locus, match_loci, paired_test, summarize, wilson)


def read(path):
    with path.open() as f:
        return list(csv.DictReader(f, delimiter='\t'))


def write(path, rows, fields=None):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields or list(rows[0]), delimiter='\t', lineterminator='\n')
        w.writeheader()
        w.writerows({k: 'NA' if v is None else v for k, v in r.items()} for r in rows)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def gff(path, accession):
    records = []
    for line in path.read_text().splitlines():
        f = line.split('\t')
        if len(f) < 9 or f[2] != 'CDS':
            continue
        a = {k: unquote(v) for k,v in (x.split('=',1) for x in f[8].split(';') if '=' in x)}
        records.append(dict(accession=accession, start=int(f[3]), end=int(f[4]), strand=f[6],
                            locus_id=a['ID'], product=a.get('product',''), evidence=json.dumps(a,sort_keys=True),
                            confidence='', functional_state=''))
    return records


def load_inputs(root):
    from Bio import SeqIO
    references, predictions, provenance = [], {'PhageMine': [], 'Pharokka': []}, []
    published = read(REPO/'docs/benchmark_v1.2/tables/reference_provenance.tsv')
    structural = read(REPO/'docs/benchmark_v1.2/tables/strict_relaxed_metrics.tsv')
    for item in published:
        a = item['accession']
        path = root/'reference'/f'{a}.gb'
        assert digest(path) == item['sha256'], f'Reference checksum mismatch: {a}'
        record = SeqIO.read(path, 'genbank')
        assert record.id == a, (record.id,a)
        cds = [f for f in record.features if f.type == 'CDS']
        expected = next(r for r in structural if r['accession']==a and r['tool']=='PhageMine')
        assert len(cds) == int(expected['reference_cds'])
        provenance.append(dict(accession=a, retrieval_source='Frozen original benchmark inputs/'+a+'.gb',
            original_retrieval_date='NOT_RECORDED', archive_date='2026-09-18', sha256=digest(path),
            reference_cds_count=len(cds), published_checksum_match=True))
        for i,f in enumerate(cds):
            q=f.qualifiers
            references.append(dict(accession=a,start=int(f.location.start)+1,end=int(f.location.end),
                strand='+' if f.location.strand==1 else '-',locus_id=q.get('protein_id',q.get('locus_tag',[f'CDS_{i+1}']))[0],
                product='; '.join(q.get('product',[])), location=str(f.location),
                compound_location=len(f.location.parts)>1, qualifiers=json.dumps(q,sort_keys=True)))
        pm = gff(root/'predictions/PhageMine'/a/'genes.gff3',a)
        ann = read(root/'predictions/PhageMine'/a/'annotation.tsv')
        by_id = {r['protein_id']:r for r in ann}
        assert len(by_id)==len(ann)==len(pm)
        states={r['protein_id']:r for r in read(root/'predictions/PhageMine'/a/'functional_classification.tsv')}
        for p in pm:
            row=by_id[p['locus_id']]
            assert locus(p)==(a,int(row['start']),int(row['end']),row['strand'])
            p.update(product=row['product'],confidence=row['confidence'],
                functional_state=states[p['locus_id']].get('functional_state',''),
                evidence=json.dumps({'annotation':row,'classification':states[p['locus_id']]},sort_keys=True))
        ph=gff(root/'predictions/Pharokka'/a/f'{a}.gff',a)
        assert {locus(r) for r in pm}=={locus(r) for r in ph}
        predictions['PhageMine'].extend(pm); predictions['Pharokka'].extend(ph)
    for tool,n in [('PhageMine',439),('Pharokka',398)]:
        assert len(predictions[tool])==804
        # Published yield has its own frozen lexical definition; truth scoring is stricter.
        import re
        assert sum(bool(p['product']) and not re.search('hypothetical|unknown function|uncharacterized|uncharacterised',p['product'],re.I) for p in predictions[tool])==n
    return references,predictions,provenance


def blinded_rows(rows, seed=12020260918):
    rng=random.Random(seed)
    review,key=[],[]
    for row in rows:
        if PENDING not in (row['PhageMine_category'],row['Pharokka_category']):
            continue
        tools=['PhageMine','Pharokka']; rng.shuffle(tools)
        ident=hashlib.sha256((row['matching_mode']+'|'+str(locus(row))).encode()).hexdigest()[:16]
        review.append(dict(adjudication_id=ident,matching_mode=row['matching_mode'],accession=row['accession'],
            start=row['start'],end=row['end'],strand=row['strand'],reference_start=row['reference_start'],
            reference_end=row['reference_end'],curated_reference_product=row['reference_product'],
            prediction_A=row[tools[0]+'_product'],prediction_B=row[tools[1]+'_product'],
            reviewer_category_A='',reviewer_category_B='',reviewer_notes=''))
        key.append(dict(adjudication_id=ident,prediction_A_tool=tools[0],prediction_B_tool=tools[1]))
    return review,key


def apply_reviews(rows, review, key, completed):
    if completed is None:
        return
    submitted=read(completed)
    expected={r['adjudication_id']:r for r in review}; mapping={r['adjudication_id']:r for r in key}
    if len({r['adjudication_id'] for r in submitted})!=len(submitted):
        raise ValueError('Duplicate adjudication ID')
    indexed={(r['matching_mode'],locus(r)):r for r in rows}
    for r in submitted:
        original=expected.get(r['adjudication_id'])
        if original is None:
            raise ValueError('Unknown adjudication ID')
        for k,v in original.items():
            if not k.startswith('reviewer_') and str(r[k])!=str(v):
                raise ValueError('Review identity/product changed: '+k)
        target=indexed[(r['matching_mode'],locus(r))]
        for side in ('A','B'):
            cat=r['reviewer_category_'+side]
            if not cat:
                continue
            tool=mapping[r['adjudication_id']]['prediction_'+side+'_tool']
            if target[tool+'_category']!=PENDING:
                raise ValueError('Cannot overwrite deterministic classification')
            if cat not in CATEGORIES or cat in {'NON_EVALUABLE_REFERENCE','ABSTENTION'} or not r['reviewer_notes'].strip():
                raise ValueError('Invalid category or missing reviewer justification')
            # Notes must document independent evidence for justified specificity.
            if cat=='COMPATIBLE_BUT_MORE_SPECIFIC' and 'evidence:' not in r['reviewer_notes'].lower():
                raise ValueError('More-specific judgments require an evidence: citation in notes')
            target[tool+'_category']=cat
            target[tool+'_review_notes']=r['reviewer_notes']


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--output',type=Path)
    parser.add_argument('--adjudications',type=Path)
    args=parser.parse_args()
    root=Path(__file__).resolve().parents[1]; out=args.output or root
    for source in read(root/'source_provenance.tsv'):
        assert digest(REPO/source['archived_path']) == source['sha256'], source['archived_path']
    refs,pred,provenance=load_inputs(root)
    rules=load_synonyms(root/'synonym_rules.tsv')
    write(out/'reference_provenance.tsv',provenance)
    truth=[{**r,'reference_category':'EVALUABLE_NAMED_REFERENCE' if informative(r['product']) else 'NON_EVALUABLE_REFERENCE'} for r in refs]
    write(out/'functional_truth_by_locus.tsv',truth)
    for tool,ps in pred.items():
        write(out/f'{tool.lower()}_predictions.tsv',ps)
    rows=[]; match_audit=[]
    # Historical structural scoring collapsed joined locations to a bounding span.
    # Do not use that span to assign a function to an unrelated overlapping CDS.
    matchable_refs=[r for r in refs if not r['compound_location']]
    ph_index={locus(p):p for p in pred['Pharokka']}
    for mode in ('exact','relaxed'):
        matches,ambiguous=match_loci(pred['PhageMine'],matchable_refs,mode)
        for i,p in enumerate(pred['PhageMine']):
            status='MATCHED' if i in matches else 'AMBIGUOUS_OVERLAP' if i in ambiguous else 'NO_REFERENCE_MATCH'
            match_audit.append(dict(matching_mode=mode,accession=p['accession'],start=p['start'],end=p['end'],strand=p['strand'],status=status,
                candidate_reference_ids=';'.join(matchable_refs[j]['locus_id'] for j in ambiguous.get(i,[]))))
            if i not in matches:
                continue
            r=matchable_refs[matches[i]]; ph=ph_index[locus(p)]
            row={k:p[k] for k in ('accession','start','end','strand')}
            row.update(matching_mode=mode,reference_id=r['locus_id'],reference_start=r['start'],reference_end=r['end'],reference_product=r['product'])
            for tool,call in [('PhageMine',p),('Pharokka',ph)]:
                row.update({tool+'_product':call['product'],tool+'_category':classify(call['product'],r['product'],rules),tool+'_review_notes':''})
            rows.append(row)
    review,key=blinded_rows(rows)
    write(out/'manual_adjudication_blinded.tsv',review,list(review[0]) if review else ['adjudication_id'])
    write(out/'manual_adjudication_key.tsv',key,list(key[0]) if key else ['adjudication_id','prediction_A_tool','prediction_B_tool'])
    apply_reviews(rows,review,key,args.adjudications)
    write(out/'functional_by_locus.tsv',rows)
    write(out/'matching_audit.tsv',match_audit)
    summaries=[]; paired={}
    for mode in ('exact','relaxed'):
        selected=[r for r in rows if r['matching_mode']==mode]
        for tool in pred:
            s=summarize([r[tool+'_category'] for r in selected])
            summaries.append(dict(tool=tool,matching_mode=mode,total_predicted_cds=len(pred[tool]),
                exact_reference_matched_cds=sum(r['matching_mode']=='exact' for r in rows),
                matched_reference_cds=len(selected),unmatched_or_ambiguous_cds=804-len(selected),**s))
        paired[mode]=paired_test([r['PhageMine_category'] for r in selected],[r['Pharokka_category'] for r in selected])
    fields=list(dict.fromkeys(k for s in summaries for k in s))
    write(out/'functional_summary.tsv',[{k:s.get(k) for k in fields} for s in summaries],fields)
    extra_summary={}
    for tool,other in [('PhageMine','Pharokka'),('Pharokka','PhageMine')]:
        extras=[]
        selected={locus(r):r for r in rows if r['matching_mode']=='exact'}
        other_index={locus(p):p for p in pred[other]}
        for p in pred[tool]:
            q=other_index[locus(p)]
            if not informative(p['product']) or informative(q['product']):
                continue
            r=selected.get(locus(p))
            extras.append({k:p[k] for k in ('accession','start','end','strand','locus_id','product')} | dict(
                other_product=q['product'],reference_product=r['reference_product'] if r else None,
                category=r[tool+'_category'] if r else 'NO_EXACT_REFERENCE_MATCH'))
        write(out/f'{tool.lower()}_only_functional_assignments.tsv',extras,list(extras[0]) if extras else ['accession','start','end','strand','category'])
        counts=dict(__import__('collections').Counter(r['category'] for r in extras))
        evaluable=[r for r in extras if r['category'] not in {'NON_EVALUABLE_REFERENCE','NO_EXACT_REFERENCE_MATCH'}]
        supported=sum(r['category'] in COMPATIBLE for r in evaluable); pending=sum(r['category']==PENDING for r in evaluable)
        n=len(evaluable)
        extra_summary[tool]=dict(total=len(extras),category_counts=counts,evaluable_exact_matched=n,
            supported_deterministically=supported,unresolved=pending,
            supported_proportion=supported/n if n and not pending else None,
            supported_ci95=wilson(supported,n) if not pending else [None,None],
            supported_identification_bounds=[supported/n,(supported+pending)/n] if n else [None,None])
    result=dict(status='PENDING_MANUAL_ADJUDICATION' if any(s['unresolved'] for s in summaries) else 'ADJUDICATED',
        excluded_compound_reference_cds=sum(r['compound_location'] for r in refs),
        implementation_sha256=digest(REPO/'src/phagemine/functional_benchmark.py'),
        runner_sha256=digest(Path(__file__)),
        source_inventory_sha256=digest(root/'source_provenance.tsv'),
        analysis_type='FROZEN_OUTPUT_REFERENCE_CONCORDANCE',summaries=summaries,paired_comparison=paired,
        exclusive_assignments=extra_summary,seed=12020260918,synonym_sha256=digest(root/'synonym_rules.tsv'),
        adjudications_sha256=digest(args.adjudications) if args.adjudications else None,
        confidence_interval_note='Wilson 95% intervals are descriptive locus-binomial intervals; within-genome dependence is not accounted for.',
        relaxed_note='Original overlap thresholds; exact matches locked first; ambiguous and many-to-one overlaps excluded from functional scoring.')
    (out/'functional_benchmark.json').write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    figure_rows=[]
    for s in summaries:
        for cat in (*CATEGORIES,PENDING):
            count=s['unresolved'] if cat==PENDING else s[cat]
            figure_rows.append(dict(tool=s['tool'],matching_mode=s['matching_mode'],category=cat,count=count))
    write(out/'figure_source.tsv',figure_rows)
    print(json.dumps({'status':result['status'],'extra':extra_summary,'summaries':[{k:s[k] for k in ['tool','matching_mode','matched_reference_cds','evaluable_named_reference_loci','named_assertions','EXACT_PRODUCT_AGREEMENT','EQUIVALENT_FUNCTION','ABSTENTION','unresolved']} for s in summaries]},indent=2))

if __name__=='__main__':
    main()
