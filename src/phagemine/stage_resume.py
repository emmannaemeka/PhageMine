"""Minimal mixed-checkpoint reconstruction for evidence stages."""
from __future__ import annotations
import json, shutil
import os
from pathlib import Path
from dataclasses import fields, asdict
from .resume import _evidence, checkpoint_reusable
from .models import Protein, EvidenceLevel
from .genome_representation import GenomeRepresentation, Topology, Orientation, Rotation
from .fusion import classify_proteins
from .context import build_context
from .mining import mine, ranked_candidates
from .quality import assess
from .reporting import write_outputs
from .genbank import write_package
from .io import read_fasta
from .alternative_evidence import alternative_models, acquire_alternative_evidence, write_alternative_evidence
from .adjudication import adjudicate, write_adjudication
from .reconciliation import reconcile_models, GeneModel
from .pfam import PfamHMMAdapter
from .vog import VOGHMMAdapter
from .phrogs import PHROGSMMseqsAdapter
from .resources import EvidenceResourceManager, ResourceType
from .swissprot import SwissProtEvidenceAdapter

def resume_stage(source: str|Path, swissprot_path: str|Path|None=None, diamond: str|None=None) -> dict:
    root=Path(source); cp=root/'checkpoints'; gene=cp/'gene_prediction'; manifest=json.loads((gene/'checkpoint_manifest.json').read_text())
    if not (gene/'original_input.fasta').is_file() or not (gene/'proteins.faa').is_file() or not (gene/'evidence.json').is_file(): raise ValueError('stage resume validation failed: incomplete gene_prediction checkpoint')
    raw=json.loads((gene/'evidence.json').read_text()); proteins=[]
    fasta={}; current=None
    for line in (gene/'proteins.faa').read_text().splitlines():
        if line.startswith('>'): current=line[1:].split()[0]; fasta[current]=''
        elif current: fasta[current]+=line.strip()
    if set(fasta) != {x.get('protein_id') for x in raw}: raise ValueError('stage resume validation failed: proteins.faa IDs do not match checkpoint evidence')
    for item in raw:
        item=dict(item); item['sequence']=fasta[item['protein_id']]; item['annotation_level']=item.get('annotation_level','UNKNOWN')
        try: item['annotation_level']=__import__('phagemine.models',fromlist=['EvidenceLevel']).EvidenceLevel(item['annotation_level'])
        except ValueError: item['annotation_level']=__import__('phagemine.models',fromlist=['EvidenceLevel']).EvidenceLevel.COMPUTATIONAL
        item['evidence']=[_evidence(e) for e in item.get('evidence',[])]
        proteins.append(Protein(**{k:item[k] for k in Protein.__dataclass_fields__ if k in item}))
    reused=[]; rerun=[]
    for name in ('pfam','vogdb','phrogs'):
        if (cp/name/'evidence.json').is_file(): reused.append(name.upper())
    swiss=cp/'swissprot'/'checkpoint_manifest.json'; old=json.loads(swiss.read_text()) if swiss.is_file() else {}
    registry_path=os.environ.get('PHAGEMINE_RESOURCE_REGISTRY',str(Path.home()/'.phagemine'/'resources.json'))
    manager=EvidenceResourceManager(registry_path)
    resources={t:manager.find(t) for t in (ResourceType.PFAM,ResourceType.VOGDB,ResourceType.SWISSPROT,ResourceType.PHROGS)}
    print('Resource resolution:')
    for t,r in resources.items(): print(f"{t.value}: {'AVAILABLE' if r else 'UNAVAILABLE'} + {r.get('path') if r else 'none'}")
    reg=resources[ResourceType.SWISSPROT]
    swissprot_path=swissprot_path or (reg or {}).get('path')
    if not swissprot_path: raise ValueError('stage resume validation failed: Swiss-Prot resource unavailable')
    adapter=SwissProtEvidenceAdapter(swissprot_path, diamond=diamond)
    if old.get('provenance',{}).get('status')=='UNAVAILABLE' or not checkpoint_reusable({'status':'REAL','provenance':old.get('provenance',{})},adapter.provenance()):
        result=adapter.analyze(proteins); rerun.append('SWISS-PROT')
        for p in proteins: p.evidence.extend(result.evidence)
    else: reused.append('SWISS-PROT')
    (root/'evidence.json').write_text(json.dumps([asdict(p) for p in proteins],indent=2,default=str))
    gid, gseq = read_fasta(gene/'analysis_genome.fasta')
    reconciliation_rows=[]
    recon_path=root/'orf_reconciliation.json'
    if recon_path.is_file(): reconciliation_rows=json.loads(recon_path.read_text()).get('records',[])
    alt_models=alternative_models(reconciliation_rows,gseq)
    stage_states={'alternative ORF evidence':'PLANNED','ORF adjudication':'PLANNED'}
    if alt_models:
        stage_states['alternative ORF evidence']='STARTED'
        pf=resources[ResourceType.PFAM]; vo=resources[ResourceType.VOGDB]; ph=resources[ResourceType.PHROGS]
        adapters=[PfamHMMAdapter(pf.get('path') if pf else None, hmmscan=shutil.which('hmmscan') if pf else None, evalue_threshold=1e-5, coverage_threshold=0.5),
                  VOGHMMAdapter(vo.get('path') if vo else None, (vo or {}).get('provenance',{}).get('annotations_path'), hmmscan=shutil.which('hmmscan') if vo else None, database_version=(vo or {}).get('version')),
                  SwissProtEvidenceAdapter(swissprot_path,diamond=diamond),
                  PHROGSMMseqsAdapter(ph.get('path') if ph else None,(ph or {}).get('provenance',{}).get('annotations_path'),shutil.which('mmseqs') if ph else None,(ph or {}).get('version'))]
        # Reconstruct from validated per-model caches; canonical adapters are
        # intentionally not rerun during interrupted-stage recovery.
        alt_models=acquire_alternative_evidence(alt_models,adapters,cp/'alternative_evidence',root/'alternative_orf_search_status.tsv',force_fresh=False)
        write_alternative_evidence(root,alt_models,{'recovery':True})
        if not (root/'alternative_orf_evidence.json').is_file() or not (root/'alternative_orf_evidence.tsv').is_file(): raise ValueError('alternative evidence outputs were not persisted')
        stage_states['alternative ORF evidence']='COMPLETED'; stage_states['ORF adjudication']='STARTED'
        evidence_map={p.protein_id:[asdict(e) for e in p.evidence] for p in proteins}
        for model in alt_models: evidence_map[model['caller_id']]=model.get('evidence',[])
        write_adjudication(root,adjudicate(reconciliation_rows,evidence_map),{'recovery':True})
        if not (root/'orf_adjudication.json').is_file() or not (root/'orf_adjudication.tsv').is_file(): raise ValueError('adjudication outputs were not persisted')
        stage_states['ORF adjudication']='COMPLETED'
    else:
        stage_states['alternative ORF evidence']='COMPLETED'; stage_states['ORF adjudication']='COMPLETED'
    classifications=classify_proteins(proteins); contexts,modules=build_context(proteins,classifications); mine(proteins); candidates=ranked_candidates(proteins); quality=assess('',proteins)
    representation=GenomeRepresentation.original(gid, gseq)
    manifest={'pipeline':'PhageMine','command':'resume-stage','stage_status':{'gene prediction':'REUSED','ORF reconciliation':'REUSED','Pfam':'REUSED','VOGDB':'REUSED','Swiss-Prot':'RERUN','PHROGs':'REUSED'},'recovery':True}
    write_outputs(root,representation, __import__('phagemine.sequencing_provenance',fromlist=['SequencingProvenance']).SequencingProvenance(),proteins,candidates,manifest,quality,gene/'original_input.fasta',classifications,contexts,modules)
    write_package(root,representation.analysis_sequence_id,representation.analysis_sequence,proteins,manifest,sequencing_provenance=__import__('phagemine.sequencing_provenance',fromlist=['SequencingProvenance']).SequencingProvenance())
    audit={'recovery_mode':'mixed_checkpoints','reused_stages':reused+['GENE PREDICTION','ORF RECONCILIATION'],'rerun_stages':rerun+['evidence integration','functional classification','context/modules','candidate ranking','QC/reporting','GenBank/export'],'stage_status':stage_states,'invalidation_reasons':{'SWISS-PROT':'previous checkpoint UNAVAILABLE' if 'SWISS-PROT' in rerun else None}}
    (root/'stage_recovery_manifest.json').write_text(json.dumps(audit,indent=2,sort_keys=True)); return audit
