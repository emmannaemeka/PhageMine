"""Evaluation-only reference comparison. No reference data enters runtime code."""
import argparse,csv,json
from pathlib import Path
from collections import Counter
from Bio import SeqIO
from Bio.Seq import Seq
STARTS={'ATG','GTG','TTG','CTG','ATT','ATC','ATA'}
STOPS={'TAA','TAG','TGA'}
def load_tsv(path):
    with Path(path).open() as f: return list(csv.DictReader(f,delimiter='\t'))
def write_tsv(path,rows):
    with Path(path).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
def model(identifier,low,high,strand):
    return dict(id=identifier,low=int(low),high=int(high),strand=strand,start=int(low if strand=='+' else high),stop=int(high if strand=='+' else low))
def cds(seq,m):
    x=seq[m['low']-1:m['high']];return x if m['strand']=='+' else str(Seq(x).reverse_complement())
def rbs(seq,m):
    s=m['start']; upstream=seq[max(0,s-31):s-1] if m['strand']=='+' else str(Seq(seq[s:s+30]).reverse_complement())
    hits=[]
    for motif in ('AGGAGG','GGAGG','AGGAG','GGAG','GAGG','AGGA'):
        for i in range(len(upstream)-len(motif)+1):
            spacer=len(upstream)-i-len(motif)
            if upstream[i:i+len(motif)]==motif and 3<=spacer<=15: hits.append(f'{motif}:spacer={spacer}')
    return upstream,';'.join(hits) or 'none_detected_not_evidence_of_absence'
def match(refs,preds):
    edges=[]
    for i,r in enumerate(refs):
        for j,p in enumerate(preds):
            overlap=max(0,min(r['high'],p['high'])-max(r['low'],p['low'])+1)
            if r['strand']!=p['strand']: continue
            exactstop=r['stop']==p['stop']
            reciprocal=overlap/max(r['high']-r['low']+1,p['high']-p['low']+1)
            if exactstop or reciprocal>=.5:
                edges.append((exactstop,r['start']==p['start'],reciprocal,overlap,-abs(r['start']-p['start']),i,j))
    used=set(); result={}
    for *_,i,j in sorted(edges,reverse=True):
        if i not in result and j not in used:result[i]=j;used.add(j)
    return result,used

def benchmark(reference,run,output,final_gbk=None,start_changes=None):
    changes = {e["protein_id"]: e for e in json.loads(Path(start_changes).read_text())} if start_changes else {}
    rec=SeqIO.read(reference,'genbank');seq=str(rec.seq).upper(); analysis=SeqIO.read(run/'analysis_genome.fasta','fasta')
    assert seq==str(analysis.seq).upper()
    refs=[model(f.qualifiers.get('protein_id',f.qualifiers.get('locus_tag',[f'CDS_{i}']))[0],int(f.location.start)+1,int(f.location.end),'+' if f.location.strand==1 else '-') for i,f in enumerate(rec.features) if f.type=='CDS']
    final=[model(x['final_protein_id'],x['start'],x['end'],x['strand']) for x in load_tsv(run/'final_gene_models.tsv')]
    if final_gbk:
        record=SeqIO.read(final_gbk,'genbank');final=[model(f.qualifiers['locus_tag'][0],int(f.location.start)+1,int(f.location.end),'+' if f.location.strand==1 else '-') for f in record.features if f.type=='CDS']
    ph=[model(x['raw_identifier'],x['start'],x['end'],x['strand']) for x in load_tsv(run/'gene_calls/raw/phanotate.tsv')]
    gv=[model(x['raw_identifier'],x['begin'],x['end'],'+' if x['strand']=='1' else '-') for x in load_tsv(run/'gene_calls/raw/prodigal_gv.tsv')]
    summary={};rows=[];candidates=[]
    for name,models in [('final',final),('PHANOTATE',ph),('Prodigal-GV',gv)]:
        matches,used=match(refs,models); exactstart=exactstop=0;causes=Counter(); strand_counts=Counter()
        for i,r in enumerate(refs):
            p=models[matches[i]] if i in matches else None
            offset=(p['start']-r['start'])*(1 if r['strand']=='+' else -1) if p else None
            startmatch=p is not None and p['start']==r['start'];stopmatch=p is not None and p['stop']==r['stop'];exactstart+=startmatch;exactstop+=stopmatch
            causes_here=[]
            rawph=next((x for x in ph if p and (x['low'],x['high'],x['strand'])==(p['low'],p['high'],p['strand'])),None)
            rawgv=next((x for x in gv if p and (x['low'],x['high'],x['strand'])==(p['low'],p['high'],p['strand'])),None)
            if not startmatch:
                strand_counts[r['strand']]+=1
                same_stop_gv=[x for x in gv if x['stop']==r['stop'] and x['strand']==r['strand']]
                if not p: causes_here=['missing_reference_gene']
                else:
                    if any(x['start']!=p['start'] for x in same_stop_gv): causes_here.append('caller_disagreement')
                    pc,rc=cds(seq,p)[:3],cds(seq,r)[:3]
                    if pc=='ATG' and rc!='ATG':causes_here.append('start_codon_preference_pattern_not_proven')
                    if offset<0:causes_here.append('longer_predicted_ORF_pattern_not_proven_longest_bias')
                    if pc not in STARTS:causes_here.append('noncanonical_boundary_start')
                    if not causes_here:causes_here=['unresolved_caller_start_choice']
                causes.update(causes_here)
                # Explicit alternatives from reference and both recorded callers.
                alternatives=[('reference',r)]+[(caller,c) for caller,arr in [('PHANOTATE',ph),('Prodigal-GV',gv)] for c in arr if c['strand']==r['strand'] and c['stop']==r['stop']]
                for caller,c in alternatives:
                    up,sd=rbs(seq,c); delta=(c['start']-r['start'])*(1 if c['strand']=='+' else -1)
                    candidates.append(dict(compared_set=name,reference_id=r['id'],candidate_source=caller,candidate_id=c['id'],start=c['start'],stop=c['stop'],strand=c['strand'],start_codon=cds(seq,c)[:3],length_nt=c['high']-c['low']+1,length_gain_loss_vs_reference_nt=-delta,upstream_30nt=up,sd_motif_screen=sd,rbs_caller_score='not_saved',internal_stop_codons=';'.join(f'{j+1}:{cds(seq,c)[j:j+3]}' for j in range(0,len(cds(seq,c))-3,3) if cds(seq,c)[j:j+3] in STOPS),recoding_evidence='not_established'))
            source=name if name!='final' else ('PHANOTATE+Prodigal-GV' if rawph and rawgv else 'PHANOTATE' if rawph else 'Prodigal-GV' if rawgv else 'export_repair_or_start_selection')
            path='phagemine.pipeline.run -> predictor.predict; PHANOTATE retained; run_automatic_observational_reconciliation does not replace final models'
            if name!='final':path='raw caller output; no reconciliation selection'
            if name=='final' and p and p['id'] in changes:
                source=changes[p['id']]['caller']
                path='phagemine.start_selection.select_supported_starts: '+changes[p['id']]['rule']
            if final_gbk and source=='export_repair_or_start_selection':path='phagemine.start_selection.select_supported_starts -> phagemine.genbank.repair_boundaries'
            rows.append(dict(reference_id=r['id'],reference_start=r['start'],reference_stop=r['stop'],reference_low=r['low'],reference_high=r['high'],reference_strand=r['strand'],prediction_set=name,predicted_id=p['id'] if p else '',predicted_start=p['start'] if p else '',predicted_stop=p['stop'] if p else '',predicted_strand=p['strand'] if p else '',start_offset_nt=offset if p else '',in_frame=(offset%3==0) if p else '',exact_start=startmatch,stop_match=stopmatch,start_supplier=source if p else '',selection_code_path=path if p else '',reference_start_codon=cds(seq,r)[:3],predicted_start_codon=cds(seq,p)[:3] if p else '',length_gain_loss_nt=(p['high']-p['low'])-(r['high']-r['low']) if p else '',mismatch_causes=';'.join(causes_here),overlap_trimming_evidence='not_recorded'))
        summary[name]=dict(reference_cds=len(refs),predicted_cds=len(models),matched=len(matches),missing_reference_genes=[r['id'] for i,r in enumerate(refs) if i not in matches],extra_predicted_genes=[p['id'] for j,p in enumerate(models) if j not in used],exact_starts=exactstart,exact_stops=exactstop,exact_start_pct=round(exactstart/len(refs)*100,4),exact_stop_pct=round(exactstop/len(refs)*100,4),mismatch_cause_labels=dict(causes),mismatch_strands=dict(strand_counts))
    output.parent.mkdir(parents=True,exist_ok=True);write_tsv(output,rows);write_tsv(output.with_name(output.stem+'_candidates.tsv'),candidates)
    summary['method']='Deterministic one-to-one matching: same strand, prioritize identical biological stop; otherwise reciprocal overlap >=50%; ties exact start, overlap, start distance. Percent denominator all reference CDSs. Cause labels describe observations, not proven caller mechanisms. RBS screen: literal SD-like motifs, 3-15 nt spacer; no calibrated RBS score available.'
    output.with_suffix('.json').write_text(json.dumps(summary,indent=2)+'\n');return summary
if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('reference',type=Path);a.add_argument('run',type=Path);a.add_argument('output',type=Path);a.add_argument('--final-gbk',type=Path);a.add_argument('--start-changes',type=Path);v=a.parse_args();print(json.dumps(benchmark(v.reference,v.run,v.output,v.final_gbk,v.start_changes),indent=2))
