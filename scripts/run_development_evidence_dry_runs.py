#!/usr/bin/env python3
"""Run the next-generation evidence engine on supplied development FASTA files.

This utility never reads benchmark/adjudication directories.  Genome IDs are
derived from the supplied FASTA headers or paths; no expected outcomes are
hard-coded.
"""
from __future__ import annotations
import argparse, json, shutil
from pathlib import Path
from phagemine.gene_callers import get_gene_model_provider, ProviderError
from phagemine.evidence_engine import adjudicate_structure, assess_architecture, aggregate_module_evidence

def records(path):
    ident=None; seq=[]
    for line in Path(path).read_text().splitlines():
        if line.startswith('>'):
            if ident is not None: yield ident, ''.join(seq)
            ident=line[1:].split()[0]; seq=[]
        elif line.strip(): seq.append(line.strip())
    if ident is not None: yield ident, ''.join(seq)

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('inputs',nargs='+'); ap.add_argument('--output',required=True); ap.add_argument('--phanotate'); a=ap.parse_args()
    root=Path(a.output); root.mkdir(parents=True,exist_ok=True)
    pg=get_gene_model_provider('prodigal_gv'); ph=get_gene_model_provider('phanotate', executable=a.phanotate)
    rows=[]
    for source in a.inputs:
        for path in sorted(Path(source).rglob('*') if Path(source).is_dir() else [Path(source)]):
            if path.suffix.lower() not in {'.fa','.fna','.fasta'}: continue
            for ident,seq in records(path):
                result={'genome_id':ident,'input':str(path),'phanotate_status':'NOT_RUN','prodigal_gv_status':'NOT_RUN','blockers':[]}
                try:
                    primary=ph.predict(ident,seq,path); result['phanotate_status']='SUCCESS'; result['phanotate_cds_count']=len(primary.models)
                except Exception as exc:
                    result['phanotate_status']='FAILED'; result['blockers'].append(f'PHANOTATE: {exc}'); primary=None
                try:
                    secondary=pg.predict(ident,seq,path); result['prodigal_gv_status']='SUCCESS'; result['prodigal_gv_cds_count']=len(secondary.models)
                except Exception as exc:
                    result['prodigal_gv_status']='FAILED'; result['blockers'].append(f'PYRODIGAL_GV: {exc}'); secondary=None
                if primary:
                    relationships=[]
                    for model in primary.models:
                        match=next((x for x in (secondary.models if secondary else []) if x.start==model.start and x.end==model.end and x.strand==model.strand),None)
                        row=adjudicate_structure(f'{ident}:{model.raw_identifier}',model,match,[])
                        relationships.append(row.relationship)
                        result.setdefault('structural_confidence',{}).setdefault(row.structural_cds_confidence,0); result['structural_confidence'][row.structural_cds_confidence]+=1
                    result['relationship_counts']={x:relationships.count(x) for x in sorted(set(relationships))}
                    result['exact_matches']=relationships.count('EXACT_MATCH')
                    result['same_stop_different_start']=relationships.count('SAME_STOP_DIFFERENT_START')
                    result['near_boundary_matches']=relationships.count('NEAR_BOUNDARY_MATCH')
                    result['phanotate_only']=relationships.count('PHANOTATE_ONLY')
                    result['conflicting_orfs']=relationships.count('CONFLICTING_ORF')+relationships.count('OVERLAPPING_ALTERNATIVE')
                    result['short_cds_reconciliation']={'short_phanotate_cds_count':sum(1 for m in primary.models if m.length_nt < 150), 'short_cds_supported':0, 'short_cds_review_required':sum(1 for m in primary.models if m.length_nt < 150)}
                result['pyrodigal_gv_only']=max(0, (len(secondary.models) if secondary else 0) - result.get('exact_matches',0))
                result['functional_class_distribution']={}
                modules=aggregate_module_evidence([])
                result['module_evidence']=modules
                result['architecture_assessment']=assess_architecture(modules)
                result['inphared']={'status':'NOT_RUN','reason':'No comparative database was supplied to this dry run.'}
                result['comparative_evidence']={'state':'INSUFFICIENT_DATA','significant':False}
                result['review_required_cases']=sum(1 for x in result.get('structural_confidence',{}) if x in {'LOW','UNRESOLVED'})
                rows.append(result)
    (root/'development_dry_run.json').write_text(json.dumps({'development_only':True,'results':rows},indent=2,sort_keys=True)+'\n')
    print(json.dumps({'development_only':True,'genomes':len(rows),'pyrodigal_gv_operational':pg.available(),'phanotate_operational':ph.available(),'output':str(root/'development_dry_run.json')},indent=2))

if __name__=='__main__': main()
