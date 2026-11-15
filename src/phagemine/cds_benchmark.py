"""Reference-only evaluation; never imported by gene calling or curation."""
from pathlib import Path
import argparse
import json
from urllib.parse import unquote
from Bio import SeqIO
from Bio.Seq import Seq


def read_cds(path, fasta=None):
    path=Path(path); rows=[]
    if path.suffix.lower() in {'.gff','.gff3'}:
        if not fasta: raise ValueError('GFF3 requires --fasta')
        seqs={r.id:str(r.seq) for r in SeqIO.parse(fasta,'fasta')}
        for line in path.read_text().splitlines():
            if line.startswith('##FASTA'): break
            if not line or line.startswith('#'): continue
            f=line.split('\t')
            if len(f)!=9: raise ValueError('Invalid GFF3 line')
            if f[2]!='CDS': continue
            attrs=dict(part.split('=',1) for part in f[8].split(';') if '=' in part)
            a,b=int(f[3]),int(f[4]);nt=seqs[f[0]][a-1:b]
            if f[6]=='-':nt=str(Seq(nt).reverse_complement())
            rows.append(dict(id=unquote(attrs.get('locus_tag',attrs.get('ID','unknown'))),contig=f[0],start=a,end=b,strand=f[6],start_codon=nt[:3]))
    else:
        for record in SeqIO.parse(path,'genbank'):
            for f in record.features:
                if f.type!='CDS': continue
                nt=str(f.extract(record.seq))
                rows.append(dict(id=f.qualifiers.get('locus_tag',f.qualifiers.get('protein_id',['unknown']))[0],contig=record.id,start=int(f.location.start)+1,end=int(f.location.end),strand='+' if f.location.strand==1 else '-',start_codon=nt[:3],stop_codon=nt[-3:],length_nt=len(nt),parts=[[int(part.start)+1,int(part.end),part.strand] for part in f.location.parts],location=str(f.location),qualifiers=f.qualifiers))
    return rows


def compare(predictions, reference):
    used=set();exact=[];stop_only=[];missed=[];structure_mismatch=[]
    def parts(p):
        return p.get('parts', [[p['start'],p['end'],1 if p['strand']=='+' else -1]])
    for ref in reference:
        stop=ref['end'] if ref['strand']=='+' else ref['start']
        choices=[(i,p) for i,p in enumerate(predictions) if i not in used and p['contig']==ref['contig'] and p['strand']==ref['strand'] and (p['end'] if p['strand']=='+' else p['start'])==stop]
        choices.sort(key=lambda ip:(parts(ip[1])!=parts(ref),ip[1]['start']!=ref['start'] or ip[1]['end']!=ref['end'],abs(ip[1]['start']-ref['start'])+abs(ip[1]['end']-ref['end'])))
        if not choices: missed.append(ref);continue
        i,p=choices[0];used.add(i)
        row={'prediction':p,'reference':ref}
        (exact if parts(p)==parts(ref) else structure_mismatch if (p['start'],p['end'])==(ref['start'],ref['end']) else stop_only).append(row)
    n=len(reference)
    return {'reference_count':n,'prediction_count':len(predictions),'exact_count':len(exact),'exact_rate':len(exact)/n if n else None,'stop_only_count':len(stop_only),'stop_only_rate':len(stop_only)/n if n else None,'stop_only':stop_only,'structure_mismatch':structure_mismatch,'exact':exact,'missed':missed,'unmatched_predictions':[p for i,p in enumerate(predictions) if i not in used], 'note':'Unmatched predictions are benchmark false positives, not proof that the genes are spurious.'}


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--prediction',required=True);p.add_argument('--reference',required=True);p.add_argument('--fasta');p.add_argument('--output',required=True)
    a=p.parse_args();report=compare(read_cds(a.prediction,a.fasta),read_cds(a.reference));Path(a.output).write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))

if __name__=='__main__': main()
