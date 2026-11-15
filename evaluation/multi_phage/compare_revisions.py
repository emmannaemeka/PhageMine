"""Independent before/after regression and cross-export consistency report."""
import argparse
import csv
import json
from pathlib import Path
from Bio import SeqIO


def check_exports(root):
    errors=[]
    record=SeqIO.read(root/'genbank_submission/genome.gbk','genbank')
    features={f.qualifiers['locus_tag'][0]:f for f in record.features if f.type=='CDS'}
    with (root/'genbank_submission/orfs.tsv').open() as handle:
        rows=list(csv.DictReader(handle,delimiter='\t'))
    if set(features)!={r['locus_tag'] for r in rows}: errors.append('ORF table/GenBank locus sets differ')
    for r in rows:
        f=features.get(r['locus_tag'])
        if f is None:continue
        comparisons={'start':int(f.location.start)+1,'end':int(f.location.end),
                     'strand':'+' if f.location.strand==1 else '-',
                     'nucleotide_sequence':str(f.extract(record.seq)),
                     'protein_sequence':f.qualifiers.get('translation',[''])[0]}
        for k,v in comparisons.items():
            if str(v)!=r[k]:errors.append(f"{r['locus_tag']}: {k} differs between ORF table and GenBank")
    gff={}
    for line in (root/'genes.gff3').read_text().splitlines():
        if not line or line.startswith('#'):continue
        c=line.split('\t')
        if len(c)!=9 or c[2]!='CDS':continue
        attrs=dict(x.split('=',1) for x in c[8].split(';') if '=' in x)
        gff[attrs['ID']]=c
    for tag,f in features.items():
        c=gff.get(tag)
        if c is None:errors.append(f'{tag}: missing GFF3 CDS');continue
        want=(int(f.location.start)+1,int(f.location.end),'+' if f.location.strand==1 else '-',int(f.qualifiers['codon_start'][0])-1)
        got=(int(c[3]),int(c[4]),c[6],int(c[7]))
        if want!=got:errors.append(f'{tag}: GFF3 bounds/strand/phase {got} != GenBank {want}')
    return errors


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--before',type=Path,required=True);ap.add_argument('--after',type=Path,required=True);ap.add_argument('--runs',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);a=ap.parse_args()
    old={r['accession']:r for r in json.loads(a.before.read_text())['genomes']}
    results=[]
    for new in json.loads(a.after.read_text())['genomes']:
        acc=new['accession'];before=old[acc]
        if new['state']!='COMPLETE': results.append({'accession':acc,'error':'NO_EXPORT'});continue
        def models(r):
            return {p['id']:p for p in [m['prediction'] for group in ('exact','stop_only','structure_mismatch') for m in r[group]]+r['unmatched_predictions']}
        previous,current=models(before),models(new);changes=[]
        for tag,p in current.items():
            q=previous.get(tag)
            if q and (p['start'],p['end'],p['strand'])!=(q['start'],q['end'],q['strand']):
                changes.append({'id':tag,'before':[q['start'],q['end'],q['strand']],'after':[p['start'],p['end'],p['strand']],'before_start_codon':q['start_codon'],'after_start_codon':p['start_codon'],
                    'same_biological_stop':p['strand']==q['strand'] and (p['end'] if p['strand']=='+' else p['start'])==(q['end'] if q['strand']=='+' else q['start'])})
        was_exact={m['reference']['id'] for m in before['exact']};now_exact={m['reference']['id'] for m in new['exact']}
        consistency=check_exports(a.runs/acc)
        evidence=json.loads((a.runs/acc/'evidence.json').read_text())
        rna=json.loads((a.runs/acc/'rna_features.json').read_text())['features']
        rescued_rna=[{'id':p['protein_id'],'rna':f['feature_id']} for p in evidence if p['gene_call_source']=='short_orf_rescue' for f in rna if min(p['end'],f['end'])>=max(p['start'],f['start'])]
        results.append({'accession':acc,'before_exact':before['exact_count'],'after_exact':new['exact_count'],
            'reference_count':new['reference_count'],'strict_verifier':new['codon_verifier'],
            'lost_exact':sorted(was_exact-now_exact),'gained_exact':sorted(now_exact-was_exact),
            'coordinate_changes':changes,'cross_export_errors':consistency,'rescued_rna_overlaps':rescued_rna,
            'start_context':json.loads((a.runs/acc/'gene_calls/start_refinement/manifest.json').read_text()).get('start_context'),
            'overlap_promotions':json.loads((a.runs/acc/'gene_calls/start_refinement/overlap_rescue.json').read_text())['promoted']})
    a.output.mkdir(parents=True,exist_ok=True)
    (a.output/'revision_comparison.json').write_text(json.dumps(results,indent=2)+'\n')
    lines=['# Before/after verification','','| Genome | Exact before → after | Strict CDS checks | Cross-export consistency | Rescued RNA overlaps |','|---|---:|---:|---|---:|']
    for r in results:
        if 'error' in r:lines.append(f"| {r['accession']} | {r['error']} | — | — | — |");continue
        v=r['strict_verifier'];lines.append(f"| {r['accession']} | {r['before_exact']} → {r['after_exact']} / {r['reference_count']} | {v['passed']}/{v['total']} | {'FAIL' if r['cross_export_errors'] else 'PASS'} | {len(r['rescued_rna_overlaps'])} |")
        if v['returncode'] or r['cross_export_errors']:lines.extend(['','```json',json.dumps(r,indent=2),'```'])
    (a.output/'revision_comparison.md').write_text('\n'.join(lines)+'\n')
    print('\n'.join(lines))

if __name__=='__main__':main()
