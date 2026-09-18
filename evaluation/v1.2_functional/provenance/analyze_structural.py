#!/usr/bin/env python3
"""Read-only structural interpretation of frozen corrected benchmark outputs."""
from pathlib import Path
import csv, re, json, math, shutil
from collections import Counter, defaultdict
from Bio import SeqIO

ROOT=Path('/Users/emmanuelnnadi/Documents/PhageMine_v1.2_development')
BASE=ROOT/'evaluation/v1.2_tool_comparison'
CORR=BASE/'correction_validation_20260916'
OUT=BASE/'final_structural_interpretation_20260916'
TAB=OUT/'tables'; FIG=OUT/'figures'; VAL=OUT/'validation'; PROV=OUT/'provenance'
for d in (TAB,FIG,VAL,PROV): d.mkdir(parents=True,exist_ok=True)

PANEL=[('T4','NC_000866.4'),('Lambda','NC_001416.1'),('T7','NC_001604.1'),('T5','NC_005859.1'),('PhiX174','NC_001422.1'),('P22','NC_002371.2'),('Mu','NC_000929.1')]
ACC_NAME={a:g for g,a in PANEL}

def attrs(s):
    d={}
    for x in s.strip().split(';'):
        if '=' in x:
            k,v=x.split('=',1); d[k]=v
    return d

def parse_gff(path,tool,acc):
    out=[]
    with open(path) as f:
        for line in f:
            if not line.strip() or line.startswith('#'): continue
            p=line.rstrip('\n').split('\t')
            if len(p)<9 or p[2].lower()!='cds': continue
            st,en=sorted((int(p[3]),int(p[4]))); strand=p[6]; at=attrs(p[8])
            out.append({'tool':tool,'genome':acc,'id':at.get('ID',at.get('locus_tag','')),'start':st,'end':en,'strand':strand,'length_nt':en-st+1,'length_aa':(en-st+1)//3,'product':at.get('product',at.get('function','')),'attrs':at})
    return out

def parse_ref(path,acc):
    out=[]
    rec=next(SeqIO.parse(str(path),'genbank'))
    for i,feat in enumerate(rec.features):
        if feat.type!='CDS': continue
        st=int(feat.location.start)+1; en=int(feat.location.end)
        strand='+' if feat.location.strand==1 else '-' if feat.location.strand==-1 else '?'
        q=feat.qualifiers; product=';'.join(q.get('product',['']))
        out.append({'tool':'Reference','genome':acc,'id':f'REF_{i+1:04d}','start':min(st,en),'end':max(st,en),'strand':strand,'length_nt':abs(en-st)+1,'length_aa':(abs(en-st)+1)//3,'product':product,'attrs':q})
    return out

def overlap(a,b): return max(0,min(a['end'],b['end'])-max(a['start'],b['start'])+1)
def reciprocal(a,b):
    o=overlap(a,b); return (o/a['length_nt'] if a['length_nt'] else 0, o/b['length_nt'] if b['length_nt'] else 0)
def same_tuple(a,b): return a['start']==b['start'] and a['end']==b['end'] and a['strand']==b['strand']

def best_ref(pred,refs):
    # prioritize same-strand overlap, then antisense overlap, then nearest interval
    scored=[]
    for r in refs:
        o=overlap(pred,r); rs=reciprocal(pred,r)
        same=pred['strand']==r['strand']
        dist=0 if o else min(abs(pred['end']-r['start']),abs(r['end']-pred['start']))
        scored.append((1 if same and o else 0, o, min(rs), -dist, r))
    return max(scored,key=lambda x:x[:4])[-1] if scored else None

def category(pred,ref):
    if ref is None: return 'I_no_meaningful_overlap'
    o=overlap(pred,ref); rpred,rref=reciprocal(pred,ref)
    if same_tuple(pred,ref): return 'A_exact_start_stop'
    if pred['strand']==ref['strand']:
        if pred['end']==ref['end']: return 'B_same_stop_alternative_start'
        if pred['start']==ref['start']: return 'C_same_start_alternative_stop'
        if rpred>=0.5 and rref>=0.5: return 'D_substantial_same_strand_overlap'
        if o>0: return 'E_partial_same_strand_overlap'
        # no overlap: classify as intergenic only when within 1000 nt of a reference feature
        dist=min(abs(pred['end']-ref['start']),abs(ref['end']-pred['start']))
        return 'H_small_orf_absent' if pred['length_aa']<100 else ('G_intergenic_absent' if dist<=1000 else 'I_no_meaningful_overlap')
    if o>0: return 'F_antisense_overlap'
    dist=min(abs(pred['end']-ref['start']),abs(ref['end']-pred['start']))
    return 'H_small_orf_absent' if pred['length_aa']<100 else ('G_intergenic_absent' if dist<=1000 else 'I_no_meaningful_overlap')

def write_tsv(path,rows,fields=None):
    if fields is None: fields=list(rows[0]) if rows else []
    with open(path,'w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',extrasaction='ignore'); w.writeheader(); w.writerows(rows)
    with open(path.with_suffix('.csv'),'w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields,extrasaction='ignore'); w.writeheader(); w.writerows(rows)

allpred=defaultdict(list); refs={}; pm_ann={}
for genome,acc in PANEL:
    refs[acc]=parse_ref(BASE/'inputs'/f'{acc}.gb',acc)
    allpred[(acc,'PhageMine')]=parse_gff(CORR/'outputs/phagemine'/acc/'genes.gff3','PhageMine',acc)
    allpred[(acc,'Pharokka')]=parse_gff(BASE/'raw_outputs/pharokka'/acc/f'{acc}.gff','Pharokka',acc)
    allpred[(acc,'Prokka')]=parse_gff(CORR/'outputs/prokka'/acc/f'{acc}.gff','Prokka',acc)
    ap=CORR/'outputs/phagemine'/acc/'annotation.tsv'
    if ap.exists():
        with open(ap) as f:
            for row in csv.DictReader(f,delimiter='\t'): pm_ann[(acc,int(row['start']),int(row['end']),row['strand'])]=row

pred_rows=[]
for (acc,tool),preds in allpred.items():
    for p in preds:
        r=best_ref(p,refs[acc]); cat=category(p,r); o=overlap(p,r) if r else 0; rp,rr=reciprocal(p,r) if r else (0,0)
        ann=pm_ann.get((acc,p['start'],p['end'],p['strand']),{}) if tool=='PhageMine' else {}
        pred_rows.append({'genome':ACC_NAME[acc],'accession':acc,'tool':tool,'prediction_id':p['id'],'start':p['start'],'end':p['end'],'strand':p['strand'],'length_nt':p['length_nt'],'length_aa':p['length_aa'],'category':cat,'closest_reference':r['id'] if r else '','ref_start':r['start'] if r else '','ref_end':r['end'] if r else '','ref_strand':r['strand'] if r else '','overlap_nt':o,'pred_reciprocal_overlap':round(rp,6),'ref_reciprocal_overlap':round(rr,6),'start_difference':p['start']-r['start'] if r else '','stop_difference':p['end']-r['end'] if r else '','product':p['product'],'functional_evidence':ann.get('evidence_sources',''),'best_evidence':ann.get('best_evidence',''),'confidence':ann.get('confidence','')})
write_tsv(TAB/'all_prediction_structural_classifications.tsv',pred_rows)

# per-tool/per-genome strict and relaxed metrics
metric=[]
for genome,acc in PANEL:
  for tool in ['PhageMine','Pharokka','Prokka']:
    ps=allpred[(acc,tool)]; rs=refs[acc]; cats=[x for x in pred_rows if x['accession']==acc and x['tool']==tool]
    exact=sum(x['category']=='A_exact_start_stop' for x in cats)
    relaxed=sum(any(p['strand']==r['strand'] and overlap(p,r)/r['length_nt']>=0.5 and reciprocal(p,r)[0]>=0.2 for p in ps) for r in rs)
    relaxed_pred=sum(any(p['strand']==r['strand'] and overlap(p,r)/r['length_nt']>=0.5 and reciprocal(p,r)[0]>=0.2 for r in rs) for p in ps)
    # one-to-one relaxed matches for counts
    tp=relaxed; fp=len(ps)-relaxed_pred; fn=len(rs)-tp; prec=tp/(tp+fp) if tp+fp else 0; rec=tp/(tp+fn) if tp+fn else 0
    metric.append({'genome':genome,'accession':acc,'tool':tool,'predicted_cds':len(ps),'reference_cds':len(rs),'strict_tp':exact,'strict_fp':len(ps)-exact,'strict_fn':len(rs)-exact,'strict_precision':round(exact/len(ps),6) if ps else 0,'strict_recall':round(exact/len(rs),6) if rs else 0,'strict_f1':round(2*exact/(len(ps)+len(rs)),6) if ps and rs else 0,'relaxed_tp_reference':tp,'relaxed_fp_prediction':fp,'relaxed_fn_reference':fn,'relaxed_precision':round(prec,6),'relaxed_recall':round(rec,6),'relaxed_f1':round(2*prec*rec/(prec+rec),6) if prec+rec else 0})
write_tsv(TAB/'strict_relaxed_metrics.tsv',metric)

# category aggregate/per-genome
catrows=[]
for genome,acc in PANEL:
 for tool in ['PhageMine','Pharokka','Prokka']:
  c=Counter(x['category'] for x in pred_rows if x['accession']==acc and x['tool']==tool)
  catrows.append({'genome':genome,'accession':acc,'tool':tool,**{k:c[k] for k in ['A_exact_start_stop','B_same_stop_alternative_start','C_same_start_alternative_stop','D_substantial_same_strand_overlap','E_partial_same_strand_overlap','F_antisense_overlap','G_intergenic_absent','H_small_orf_absent','I_no_meaningful_overlap']}})
write_tsv(TAB/'prediction_category_counts.tsv',catrows)

# PhageMine non-exact audit with evidence and cross-tool presence
non=[]
for x in pred_rows:
 if x['tool']!='PhageMine' or x['category']=='A_exact_start_stop': continue
 acc=x['accession']; p=next(p for p in allpred[(acc,'PhageMine')] if p['id']==x['prediction_id'])
 others=[]
 for tool in ['Pharokka','Prokka']:
  q=next((q for q in allpred[(acc,tool)] if same_tuple(p,q)),None); others.append(q is not None)
 non.append({**x,'pharokka_exact_same_coordinate':others[0],'prokka_exact_same_coordinate':others[1],'is_small_lt30aa':p['length_aa']<30,'is_small_lt50aa':p['length_aa']<50,'is_small_lt100aa':p['length_aa']<100,'overlaps_other_prediction':any(overlap(p,q)>0 and not same_tuple(p,q) for q in allpred[(acc,'PhageMine')]),'functional_support':bool(x['functional_evidence']),'evidence_sources':x['functional_evidence']})
write_tsv(TAB/'phagemine_nonexact_prediction_audit.tsv',non)

# Reference-centric missed/alternative detections
fnrows=[]
for genome,acc in PANEL:
 for tool in ['PhageMine','Pharokka','Prokka']:
  for r in refs[acc]:
   exact=next((p for p in allpred[(acc,tool)] if same_tuple(p,r)),None)
   same=[p for p in allpred[(acc,tool)] if p['strand']==r['strand'] and overlap(p,r)>0]
   best=max(same,key=lambda p: reciprocal(p,r)[1],default=None)
   if exact: status='exact'
   elif best and best['end']==r['end']: status='alternative_start'
   elif best and best['start']==r['start']: status='alternative_stop'
   elif best and reciprocal(best,r)[1]>=0.5 and reciprocal(best,r)[0]>=0.2: status='substantial_overlap'
   elif best and overlap(best,r)>0: status='partial_overlap'
   else: status='missed'
   if status!='exact': fnrows.append({'genome':genome,'accession':acc,'tool':tool,'reference_id':r['id'],'ref_start':r['start'],'ref_end':r['end'],'ref_strand':r['strand'],'ref_length_nt':r['length_nt'],'ref_length_aa':r['length_aa'],'status':status,'prediction_id':best['id'] if best else '','pred_start':best['start'] if best else '','pred_end':best['end'] if best else '','overlap_nt':overlap(best,r) if best else 0,'pred_reciprocal_overlap':round(reciprocal(best,r)[0],6) if best else 0,'ref_reciprocal_overlap':round(reciprocal(best,r)[1],6) if best else 0,'small_ref_lt30aa':r['length_aa']<30,'small_ref_lt50aa':r['length_aa']<50,'small_ref_lt100aa':r['length_aa']<100})
write_tsv(TAB/'reference_nonexact_and_missed.tsv',fnrows)

# PhageMine versus Prokka additional predictions
add=[]
for genome,acc in PANEL:
 pro={(p['start'],p['end'],p['strand']) for p in allpred[(acc,'Prokka')]}
 for p in allpred[(acc,'PhageMine')]:
  if (p['start'],p['end'],p['strand']) in pro: continue
  r=best_ref(p,refs[acc]); x=next(x for x in pred_rows if x['accession']==acc and x['tool']=='PhageMine' and x['prediction_id']==p['id'])
  # evidence of a Prokka alternative-boundary/overlap model
  pp=[q for q in allpred[(acc,'Prokka')] if q['strand']==p['strand'] and overlap(p,q)>0]
  add.append({'genome':genome,'accession':acc,'phagemine_prediction_id':p['id'],'start':p['start'],'end':p['end'],'strand':p['strand'],'length_aa':p['length_aa'],'reference_category':x['category'],'closest_reference':r['id'] if r else '','reference_overlap_nt':overlap(p,r) if r else 0,'reference_overlap_fraction':round(reciprocal(p,r)[1],6) if r else 0,'prokka_overlapping_predictions':len(pp),'prokka_has_overlap':bool(pp),'functional_support':bool(x['functional_evidence']),'evidence_sources':x['functional_evidence'],'best_evidence':x['best_evidence']})
write_tsv(TAB/'phagemine_vs_prokka_additional_cds.tsv',add)

# Functional summary from existing output labels; avoid treating names as accuracy
func=[]
for genome,acc in PANEL:
 for tool in ['PhageMine','Pharokka','Prokka']:
  ps=allpred[(acc,tool)]
  if tool=='PhageMine': named=sum(bool(p['product']) and not re.search(r'hypothetical|uncharacterized|unknown',p['product'],re.I) for p in ps); hyp=len(ps)-named
  elif tool=='Prokka': named=sum(bool(p['product']) and not re.search(r'hypothetical|uncharacterized',p['product'],re.I) for p in ps); hyp=len(ps)-named
  else:
   named=sum(bool(p['product']) and not re.search(r'hypothetical|unknown function|uncharacterized',p['product'],re.I) for p in ps); hyp=len(ps)-named
  func.append({'genome':genome,'accession':acc,'tool':tool,'predicted_cds':len(ps),'named_product_calls':named,'hypothetical_or_uncharacterized':hyp,'named_fraction':round(named/len(ps),6) if ps else 0})
write_tsv(TAB/'functional_yield_by_genome.tsv',func)

# Reproducibility/provenance summary
prov={'reference_panel':str(BASE/'analysis/next_validation/REFERENCE_VALIDATION/REFERENCE_VALIDATION_PANEL.tsv'),'reference_source':'seven frozen inputs/*.gb files listed by the panel; no reference files modified','prediction_sources':{'PhageMine':str(CORR/'outputs/phagemine'),'Pharokka':str(BASE/'raw_outputs/pharokka'),'Prokka':str(CORR/'outputs/prokka')},'thresholds':{'strict':'same accession, start, end and strand','relaxed_primary':'same strand, >=50% of reference span and >=20% of prediction span','relaxed_sensitivity':'same strand, >=20% of reference span and >=20% reciprocal span','substantial_overlap':'same strand and reciprocal overlap >=0.5','partial_overlap':'same-strand overlap >0 below substantial threshold','small_orf':'predicted or reference length <100 aa','intergenic':'no overlap and nearest reference distance <=1000 nt'}}
json.dump(prov,open(PROV/'method_provenance.json','w'),indent=2)

# Figures
import matplotlib.pyplot as plt
import numpy as np
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
tools=['PhageMine','Pharokka','Prokka']; colors=['#2c7fb8','#7fcdbb','#d95f0e']
agg=[]
for t in tools:
 m=[x for x in metric if x['tool']==t]; agg.append((t,sum(x['strict_tp'] for x in m),sum(x['predicted_cds'] for x in m),sum(x['reference_cds'] for x in m),sum(x['relaxed_tp_reference'] for x in m),sum(x['relaxed_fp_prediction'] for x in m)))
fig,ax=plt.subplots(figsize=(6.5,4)); x=np.arange(3); strict=[2*a[1]/(a[2]+a[3]) for a in agg]; relaxed=[2*a[3]/(a[4]+a[5]+a[3]) for a in agg]; ax.bar(x-.18,strict,.36,label='Strict exact-coordinate F1',color=colors); ax.bar(x+.18,relaxed,.36,label='Relaxed overlap F1',color='#636363'); ax.set_xticks(x,tools); ax.set_ylim(0,1); ax.set_ylabel('F1'); ax.set_title('Strict versus relaxed structural performance'); ax.legend(frameon=False); fig.tight_layout()
for ext in ['png','svg','pdf']: fig.savefig(FIG/f'strict_vs_relaxed_f1.{ext}',dpi=300 if ext=='png' else None)

cats=['A_exact_start_stop','B_same_stop_alternative_start','C_same_start_alternative_stop','D_substantial_same_strand_overlap','E_partial_same_strand_overlap','F_antisense_overlap','G_intergenic_absent','H_small_orf_absent','I_no_meaningful_overlap']
fig,ax=plt.subplots(figsize=(8,4.5)); bot=np.zeros(3)
for cat in cats:
 vals=[sum(1 for x in pred_rows if x['tool']==t and x['category']==cat) for t in tools]; ax.bar(tools,vals,bottom=bot,label=cat.split('_',1)[1].replace('_',' ')); bot+=vals
ax.set_ylabel('Predictions'); ax.set_title('Exact, alternative-boundary, overlap and non-reference categories'); ax.legend(fontsize=7,ncol=2,frameon=False); fig.tight_layout()
for ext in ['png','svg','pdf']: fig.savefig(FIG/f'prediction_category_composition.{ext}',dpi=300 if ext=='png' else None)

fig,axs=plt.subplots(1,2,figsize=(9,4));
for i,t in enumerate(['PhageMine','Prokka']):
 vals=[x['length_aa'] for x in pred_rows if x['tool']==t and x['category']!='A_exact_start_stop']; axs[0 if i==0 else 1].hist(vals,bins=20,color=colors[i if i<3 else 0],alpha=.85); axs[0 if i==0 else 1].set_title(t+' non-exact CDSs'); axs[0 if i==0 else 1].set_xlabel('Length (aa)'); axs[0 if i==0 else 1].set_ylabel('Count')
fig.tight_layout();
for ext in ['png','svg','pdf']: fig.savefig(FIG/('nonexact_length_distributions.'+ext),dpi=300 if ext=='png' else None)

fig,ax=plt.subplots(figsize=(6,4)); vals=[]
for t in tools:
 named=sum(x['named_product_calls'] for x in func if x['tool']==t); total=sum(x['predicted_cds'] for x in func if x['tool']==t); vals.append(named/total)
ax.bar(tools,vals,color=colors); ax.set_ylim(0,1); ax.set_ylabel('Named-product fraction'); ax.set_title('Functional annotation yield (lexical category)'); fig.tight_layout()
for ext in ['png','svg','pdf']: fig.savefig(FIG/f'functional_yield.{ext}',dpi=300 if ext=='png' else None)

fig,ax=plt.subplots(figsize=(9,4.5)); width=.25; x=np.arange(7)
for j,t in enumerate(tools):
 vals=[next(z['strict_f1'] for z in metric if z['tool']==t and z['genome']==g) for g,a in PANEL]; ax.bar(x+(j-1)*width,vals,width,label=t,color=colors[j])
ax.set_xticks(x,[g for g,a in PANEL]); ax.set_ylim(0,1); ax.set_ylabel('Strict F1'); ax.set_title('Per-genome exact-coordinate structural performance'); ax.legend(frameon=False); fig.tight_layout()
for ext in ['png','svg','pdf']: fig.savefig(FIG/f'per_genome_strict_f1.{ext}',dpi=300 if ext=='png' else None)

print(json.dumps({'prediction_rows':len(pred_rows),'nonexact_phagemine':len(non),'additional_phagemine_vs_prokka':len(add),'reference_nonexact_rows':len(fnrows)},indent=2))
