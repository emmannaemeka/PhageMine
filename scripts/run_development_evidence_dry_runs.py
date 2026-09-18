#!/usr/bin/env python3
"""Run the next-generation evidence engine on supplied development FASTA files.

This utility never reads benchmark/adjudication directories.  Genome IDs are
derived from the supplied FASTA headers or paths; no expected outcomes are
hard-coded.
"""
from __future__ import annotations
import argparse, csv, json, shutil, re
from pathlib import Path
from phagemine.gene_callers import get_gene_model_provider, ProviderError
from phagemine.evidence_engine import adjudicate_structure, assess_architecture, architecture_hallmarks, aggregate_module_evidence

def records(path):
    ident=None; seq=[]
    for line in Path(path).read_text().splitlines():
        if line.startswith('>'):
            if ident is not None: yield ident, ''.join(seq)
            ident=line[1:].split()[0]; seq=[]
        elif line.strip(): seq.append(line.strip())
    if ident is not None: yield ident, ''.join(seq)

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('inputs',nargs='+'); ap.add_argument('--output',required=True); ap.add_argument('--phanotate'); ap.add_argument('--results-root', help='Optional existing development result root containing evidence.json; never a benchmark directory'); a=ap.parse_args()
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
                        candidates=[x for x in (secondary.models if secondary else []) if x.strand==model.strand and max(0,min(model.end,x.end)-max(model.start,x.start)+1)>0]
                        match=max(candidates,key=lambda x:max(0,min(model.end,x.end)-max(model.start,x.start)+1),default=None)
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
                result['pyrodigal_gv_only']=sum(1 for m in (secondary.models if secondary else []) if not any(m.strand==p.strand and max(0,min(m.end,p.end)-max(m.start,p.start)+1)>0 for p in (primary.models if primary else [])))
                existing = Path(a.results_root)/ident if a.results_root else None
                modules=[]; functional_classes={}
                if existing and (existing/'evidence.json').is_file():
                    evidence_records=json.loads((existing/'evidence.json').read_text())
                    short_ids={r.get('protein_id') for r in evidence_records if any(str(e.get('evidence_strength','')).upper() in {'STRONG','HIGH','CURATED','EXPERIMENTAL'} and e.get('supports') is not False for e in r.get('evidence',[]))}
                    if primary:
                        short_models=[m for m in primary.models if m.length_nt < 150]
                        result['short_cds_reconciliation']['short_cds_supported']=sum(1 for m in short_models if m.raw_identifier in short_ids)
                    for record in evidence_records:
                        annotation=str(record.get('annotation') or '')
                        functional_classes['HYPOTHETICAL' if re.search(r'hypothetical|unknown|uncharacter',annotation,re.I) else 'SUPPORTED_SPECIFIC']=functional_classes.get('HYPOTHETICAL' if re.search(r'hypothetical|unknown|uncharacter',annotation,re.I) else 'SUPPORTED_SPECIFIC',0)+1
                        text=' '.join([annotation, str(record.get('gene') or '')]).lower()
                        annotation_modules=[]
                        if re.search(r'replication|rep\b|polymerase|helicase|primase',text): annotation_modules.append('REP_REPLICATION')
                        if re.search(r'zot|p/i|extrusion|morphogenesis',text): annotation_modules.append('ZOT_EXTRUSION')
                        if re.search(r'coat|capsid|head|virion|minor-head',text): annotation_modules.append('COAT_VIRION')
                        if re.search(r'membrane|structural protein',text): annotation_modules.append('MEMBRANE_STRUCTURAL')
                        if re.search(r'integrase|integration',text): annotation_modules.append('INTEGRATION')
                        if re.search(r'regulator|regulation|transcription',text): annotation_modules.append('REGULATION')
                        for e in record.get('evidence',[]):
                            mapped=[]
                            et=' '.join(str(e.get(k) or '') for k in ('statement','description','family_name','functional_category')).lower()
                            if re.search(r'replication|rep\b|polymerase|helicase|primase',et): mapped.append('REP_REPLICATION')
                            if re.search(r'zot|p/i|extrusion|morphogenesis',et): mapped.append('ZOT_EXTRUSION')
                            if re.search(r'coat|capsid|head|virion|minor-head',et): mapped.append('COAT_VIRION')
                            if re.search(r'membrane|structural protein',et): mapped.append('MEMBRANE_STRUCTURAL')
                            if re.search(r'integrase|integration',et): mapped.append('INTEGRATION')
                            if re.search(r'regulator|regulation|transcription',et): mapped.append('REGULATION')
                            if re.search(r'tail|tape measure|tape-measure',et): mapped.append('TAPE_MEASURE' if 'tape' in et else 'TAIL')
                            if re.search(r'holin|endolysin|lysin|spanin',et): mapped.append('LYSIS')
                            strength=e.get('evidence_strength') or 'WEAK'; role='SUPPORT' if e.get('supports') else 'CONFLICT'
                            if not mapped: mapped=list(annotation_modules); strength='WEAK'; role='CONFLICT' if not e.get('supports') else 'SUPPORT'
                            for mod in mapped: modules.append({'module':mod,'source_database':e.get('source'),'evidence_strength':strength,'role':role,'locus_id':record.get('protein_id')})
                result['functional_class_distribution']=functional_classes
                modules=aggregate_module_evidence([{'locus_id':x.get('locus_id',''),'evidence':[x]} for x in modules]) if modules else aggregate_module_evidence([])
                result['module_evidence']=modules
                result['architecture_assessment']=assess_architecture(modules)
                result['architecture_hallmarks']=architecture_hallmarks(result['architecture_assessment']['architecture_hypothesis'], modules)
                comparative_file=existing/'comparative'/'inphared_summary.tsv' if existing else None
                if comparative_file and comparative_file.is_file():
                    with comparative_file.open() as handle: comparative_rows=list(csv.DictReader(handle,delimiter='\t'))
                    best=comparative_rows[0] if comparative_rows else {}
                    hashes=str(best.get('matching_hashes','')).split('/')[0]
                    zero_hashes=hashes.isdigit() and int(hashes)==0
                    comparative_state='NO_SIGNIFICANT_REFERENCE' if zero_hashes else ('SIGNIFICANT_CLOSE_REFERENCE' if str(best.get('taxonomy_eligible','')).lower()=='true' and float(best.get('intergenomic_similarity_percent') or 0)>=95 else 'WEAK_REFERENCE')
                    result['inphared']={'status':'AVAILABLE','comparative_state':comparative_state,'nearest_candidate':best.get('reference_description'),'best_mash_distance':best.get('best_mash_distance'),'matching_hashes':best.get('matching_hashes'),'taxonomic_interpretation':best.get('taxonomic_interpretation'),'provenance':str(comparative_file)}
                    result['comparative_evidence']={'state':comparative_state,'significant':comparative_state.startswith('SIGNIFICANT'),'nearest_candidate':best.get('reference_description'),'interpretation':'Nearest candidate is not a biological relative unless comparative significance criteria are met.'}
                else:
                    result['inphared']={'status':'NOT_RUN','reason':'No comparative database was supplied to this dry run.'}
                    result['comparative_evidence']={'state':'INSUFFICIENT_DATA','significant':False}
                result['review_required_cases']=sum(v for k,v in result.get('structural_confidence',{}).items() if k in {'LOW','UNRESOLVED'})
                rows.append(result)
    (root/'development_dry_run.json').write_text(json.dumps({'development_only':True,'results':rows},indent=2,sort_keys=True)+'\n')
    print(json.dumps({'development_only':True,'genomes':len(rows),'pyrodigal_gv_operational':pg.available(),'phanotate_operational':ph.available(),'output':str(root/'development_dry_run.json')},indent=2))

if __name__=='__main__': main()
