"""Model-specific evidence acquisition for discordant alternative ORFs."""
from __future__ import annotations
import hashlib, json, csv, shutil
import time
from pathlib import Path
from dataclasses import asdict
from .gene_prediction import reverse_complement
from .genome import translate
from .models import Protein

RELEVANT={"START_DISCORDANCE","STOP_DISCORDANCE","START_AND_STOP_DISCORDANCE","STRAND_DISCORDANCE","PRODIGAL_ONLY","SPLIT_MODEL","MERGED_MODEL","COMPLEX_CONFLICT"}

def _canonical_source(value):
    """Normalize adapter/source spellings for accounting only.

    Evidence objects themselves are left untouched; this prevents a producer's
    historical source spelling from splitting the authoritative counts.
    """
    v = str(value or "")
    key = v.lower().replace("_", "-").replace(" ", "")
    return {
        "pfam": "Pfam", "pfamhmmadapter": "Pfam",
        "vogdb": "VOGDB", "voghmmdapter": "VOGDB", "voghhmmadapter": "VOGDB",
        "voghm": "VOGDB", "vog": "VOGDB",
        "swiss-prot": "Swiss-Prot", "swissprot": "Swiss-Prot",
        "swissprotevidenceadapter": "Swiss-Prot",
        "phrogs": "PHROGs", "phrog": "PHROGs",
        "phrogsmmsaadapter": "PHROGs", "phrogsmmsesadapter": "PHROGs",
        "phrogsmmsseqadapter": "PHROGs",
    }.get(key, v)

def _evidence_key(e):
    # Evidence records from one hit may share source/accession and even raw
    # rows while differing in coordinates, thresholds, or provenance.  Hash
    # the complete serialized record so genuine multi-hit evidence survives
    # cache-fragment merging; only byte-identical insertions are collapsed.
    return json.dumps(e, sort_keys=True, default=str, separators=(",", ":"))

def _merge_evidence(model, records, dedupe=False):
    if not dedupe:
        model.setdefault('evidence',[]).extend(records); return
    seen={_evidence_key(e) for e in model.get('evidence',[])}
    for e in records:
        if _evidence_key(e) not in seen:
            model.setdefault('evidence',[]).append(e); seen.add(_evidence_key(e))

def _cached_records_for_model(cache, model, source_name, exact_path=None):
    """Load all validated cache fragments for one model/source.

    Interrupted batched searches can leave multiple disjoint fragments for the
    same query.  Selecting only the current hash silently dropped valid hits.
    Non-empty fragments are attributable by their stable query id; an empty
    exact-key fragment is retained to represent an executed zero-hit search.
    """
    if not cache:
        return None
    records=[]; found=False
    for path in sorted(Path(cache).glob('*.json')):
        try: payload=json.loads(path.read_text())
        except (OSError, ValueError): continue
        if not isinstance(payload,list): continue
        matching=[e for e in payload
                  if _canonical_source(e.get('source')) == source_name and
                  (e.get('provenance',{}).get('protein_id') == model['alternative_model_id'] or
                   e.get('metrics',{}).get('query_protein_id') == model['alternative_model_id'])]
        if matching:
            found=True; records.extend(matching)
        elif exact_path is not None and path == exact_path and not payload:
            found=True
    if not found:
        return None
    unique=[]; seen=set()
    for record in records:
        key=_evidence_key(record)
        if key not in seen:
            seen.add(key); unique.append(record)
    return unique

def adapter_identity(adapter):
    """Return stable implementation/resource/configuration identity."""
    provenance = adapter.provenance() if callable(getattr(adapter, "provenance", None)) else getattr(adapter, "provenance", {})
    if not isinstance(provenance, dict): provenance = {}
    return {"class": adapter.__class__.__name__, "name": getattr(adapter, "name", None), "provenance": provenance}

def cache_key(model, adapter):
    payload={"translation_sha256":model["translation_sha256"],"adapter":adapter_identity(adapter)}
    return hashlib.sha256(json.dumps(payload,sort_keys=True,default=str,separators=(",",":")).encode()).hexdigest()

def extract_translation(sequence: str, start: int, end: int, strand: str) -> tuple[str,str]:
    genomic=sequence[start-1:end]; cds=reverse_complement(genomic) if strand == "-" else genomic
    return cds, translate(cds)

def alternative_models(rows, sequence: str):
    result=[]; index=0
    for row in rows:
        if row.get("conflict_type") not in RELEVANT or not row.get("prodigal_id"): continue
        index+=1; cds, protein=extract_translation(sequence,row["prodigal_start"],row["prodigal_end"],row["prodigal_strand"])
        result.append({"locus_id":row["locus_id"],"alternative_model_id":f"ALT_PRODIGAL_{index:05d}","caller":"Prodigal","caller_id":row["prodigal_id"],"start":row["prodigal_start"],"end":row["prodigal_end"],"strand":row["prodigal_strand"],"cds":cds,"protein_sequence":protein,"translation_sha256":hashlib.sha256(protein.encode()).hexdigest(),"evidence":[]})
    return result

def acquire_alternative_evidence(models, adapters=(), cache_dir=None, status_output=None, force_fresh=False):
    cache=Path(cache_dir) if cache_dir else None
    if cache:
        if force_fresh and cache.exists():
            for item in cache.glob('*.json'): item.unlink()
        cache.mkdir(parents=True,exist_ok=True)
    statuses=[]
    for adapter in adapters:
        pending=[]; cached={}; started=time.monotonic(); identity=adapter_identity(adapter)
        source_name=_canonical_source({'PfamHMMAdapter':'Pfam','VOGHMMAdapter':'VOGDB','SwissProtEvidenceAdapter':'Swiss-Prot','PHROGSMMseqsAdapter':'PHROGs'}.get(identity['class'], identity['name'] or identity['class']))
        for model in models:
            key=cache_key(model, adapter); path=cache/f"{key}.json" if cache else None
            if path and path.is_file() and not force_fresh:
                records=_cached_records_for_model(cache, model, source_name, path)
                if records is not None:
                    cached[model['alternative_model_id']]=records; continue
            pending.append(model)
        if pending:
            model_by_id={m["alternative_model_id"]:m for m in pending}
            proteins=[Protein("alternative",m["alternative_model_id"],m["start"],m["end"],m["strand"],m["cds"],m["protein_sequence"],"Prodigal") for m in pending]
            result=adapter.analyze(proteins); grouped={model_id:[] for model_id in model_by_id}
            for e in result.evidence:
                query_id=e.provenance.get('protein_id') or e.metrics.get('query_protein_id')
                if query_id in model_by_id: grouped[query_id].append(asdict(e))
            for model in pending:
                records=grouped.get(model['alternative_model_id'],[]); _merge_evidence(model, records)
                if cache: (cache/f"{cache_key(model,adapter)}.json").write_text(json.dumps(records,indent=2,sort_keys=True,default=str))
            status='RESOURCE_UNAVAILABLE' if result.status=='UNAVAILABLE' else ('SEARCH_EXECUTED_ZERO_HITS' if not result.evidence else 'SEARCH_EXECUTED')
            reason=result.message
        else: status='CACHE_REUSED'; reason='validated cache entries reused'
        for model in models:
            records=cached.get(model['alternative_model_id'])
            if records is not None: _merge_evidence(model, records, dedupe=True)
            # Use the records produced for this adapter, rather than the
            # adapter's result object (which is shared across the batch) or
            # the model's complete accumulated evidence.  The latter caused
            # source-name mismatches and reported zero for every row.
            if records is None:
                records = [e for e in model.get('evidence', []) if _canonical_source(e.get('source')) == source_name]
            row_status = status
            if status == 'SEARCH_EXECUTED' and not records:
                row_status = 'SEARCH_EXECUTED_ZERO_HITS'
            statuses.append({'alternative_model_id':model['alternative_model_id'],'source':source_name,'search_status':row_status,'records_returned':len(records),'accepted_records':sum(bool(e.get('supports')) for e in records),'cache_status':'REUSED' if model['alternative_model_id'] in cached else ('WRITTEN' if pending else 'NOT_USED'),'reason':reason,'runtime_seconds':time.monotonic()-started})
        # Guard the per-source bookkeeping against future regressions.
        source_records=sum(len([e for e in m.get('evidence',[]) if _canonical_source(e.get('source')) == source_name]) for m in models)
        source_status=sum(s['records_returned'] for s in statuses if s['source'] == source_name)
        source_accepted=sum(bool(e.get('supports')) for m in models for e in m.get('evidence',[]) if _canonical_source(e.get('source')) == source_name)
        source_status_accepted=sum(s['accepted_records'] for s in statuses if s['source'] == source_name)
        # Counts are finalized below from the merged model evidence.  Cache
        # rows can legitimately differ before that merge is reconciled.
    # Final reconciliation is authoritative, especially when a model already
    # contains evidence from an interrupted/cache-backed attempt.
    model_by_id={m['alternative_model_id']:m for m in models}
    for s in statuses:
        es=[e for e in model_by_id[s['alternative_model_id']].get('evidence',[]) if _canonical_source(e.get('source'))==s['source']]
        s['records_returned']=len(es); s['accepted_records']=sum(bool(e.get('supports')) for e in es)
    total_records = sum(len(m.get('evidence', [])) for m in models)
    total_accepted = sum(bool(e.get('supports')) for m in models for e in m.get('evidence', []))
    if (sum(s['records_returned'] for s in statuses), sum(s['accepted_records'] for s in statuses)) != (total_records, total_accepted):
        raise RuntimeError(f"alternative evidence status totals do not match persisted evidence: status={sum(s['records_returned'] for s in statuses)}/{sum(s['accepted_records'] for s in statuses)} evidence={total_records}/{total_accepted}")
    if status_output:
        Path(status_output).write_text('alternative_model_id\tsource\tsearch_status\trecords_returned\taccepted_records\tcache_status\treason\truntime_seconds\n'+'\n'.join(f"{s['alternative_model_id']}\t{s['source']}\t{s['search_status']}\t{s['records_returned']}\t{s['accepted_records']}\t{s['cache_status']}\t{s['reason'] or ''}\t{s['runtime_seconds']:.3f}" for s in statuses)+'\n')
    return models

def write_alternative_evidence(root, models, provenance=None):
    root=Path(root); (root/"alternative_orf_evidence.json").write_text(json.dumps({"provenance":provenance or {},"models":models},indent=2,sort_keys=True,default=str))
    cols=["locus_id","alternative_model_id","caller","caller_id","start","end","strand","protein_length","source","identifier","description","evidence_strength","supports"]
    with (root/"alternative_orf_evidence.tsv").open("w",newline="") as h:
        w=csv.DictWriter(h,fieldnames=cols,delimiter="\t"); w.writeheader()
        for m in models:
            for e in m.get("evidence",[]): w.writerow({"locus_id":m["locus_id"],"alternative_model_id":m["alternative_model_id"],"caller":m["caller"],"caller_id":m["caller_id"],"start":m["start"],"end":m["end"],"strand":m["strand"],"protein_length":len(m["protein_sequence"]),**{k:e.get(k) for k in cols if k in e}})
