"""Apply explicit human corrections to a completed submission package only."""
from __future__ import annotations
import copy, hashlib, json, shutil
from datetime import datetime, timezone
from pathlib import Path
from .genbank import write_package
from .io import read_fasta
from .models import Evidence, EvidenceLevel, Protein, SubmissionMetadata

ALLOWED_ANNOTATION_FIELDS = {"product", "note", "partial", "locus_tag"}

def _load_corrections(path: Path) -> dict:
    try:
        data = json.loads(path.read_text())
    except Exception as exc:
        raise ValueError(f"correction file is not valid JSON: {exc}") from exc
    if not isinstance(data, dict) or set(data) - {"metadata", "annotations"}:
        raise ValueError("correction file must contain only metadata and annotations sections")
    if not isinstance(data.get("metadata", {}), dict) or not isinstance(data.get("annotations", {}), dict):
        raise ValueError("metadata and annotations must be objects")
    return data

def _proteins(run: Path) -> list[Protein]:
    source = run / "checkpoints" / "gene_prediction" / "evidence.json"
    if not source.exists():
        raise ValueError("missing original gene-prediction provenance: checkpoints/gene_prediction/evidence.json")
    rows = json.loads(source.read_text()); result=[]
    for row in rows:
        allowed={k:v for k,v in row.items() if k in Protein.__dataclass_fields__}
        allowed["annotation_level"] = EvidenceLevel(allowed.get("annotation_level", EvidenceLevel.HYPOTHESIS.value))
        allowed["evidence"]=[Evidence(**{**e, "level": EvidenceLevel(e.get("level", EvidenceLevel.COMPUTATIONAL.value))}) for e in allowed.get("evidence", [])]
        result.append(Protein(**allowed))
    return result

def revise(run_dir: str | Path, corrections_path: str | Path, output: str | Path) -> dict:
    run, corr_path, out = Path(run_dir), Path(corrections_path), Path(output)
    if not run.exists(): raise ValueError("completed PhageMine result directory does not exist")
    corrections = _load_corrections(corr_path)
    fasta = run / "original_input.fasta"
    if not fasta.exists(): raise ValueError("original_input.fasta is required")
    genome_id, genome = read_fasta(fasta)
    digest = hashlib.sha256(fasta.read_bytes()).hexdigest()
    manifest_path = run / "run_manifest.json"
    if not manifest_path.exists(): raise ValueError("missing original run provenance: run_manifest.json")
    manifest=json.loads(manifest_path.read_text())
    expected=manifest.get("input_sha256")
    if expected and expected != digest: raise ValueError("Genome sequence has changed; a new PhageMine analysis is required.")
    proteins=_proteins(run); ids={p.protein_id for p in proteins}
    annotations=corrections.get("annotations", {})
    unknown=set(annotations)-ids
    if unknown: raise ValueError("unknown protein ID(s): " + ", ".join(sorted(unknown)))
    overrides={}; audit=[]; now=datetime.now(timezone.utc).isoformat()
    for pid, change in annotations.items():
        if not isinstance(change, dict): raise ValueError(f"annotation correction for {pid} must be an object")
        unsupported=set(change)-ALLOWED_ANNOTATION_FIELDS-{"reason","source"}
        if unsupported: raise ValueError(f"unsupported correction field(s) for {pid}: {', '.join(sorted(unsupported))}")
        if "product" in change and not isinstance(change["product"], str): raise ValueError(f"product correction for {pid} must be a string")
        overrides[pid]={k:change[k] for k in ALLOWED_ANNOTATION_FIELDS if k in change}
        for field,value in overrides[pid].items():
            original = getattr(next(p for p in proteins if p.protein_id==pid), "annotation", None) if field=="product" else None
            audit.append({"protein_id":pid,"field":field,"original_value":original,"corrected_value":value,"reason":change.get("reason"),"correction_source":change.get("source","manual_review"),"timestamp":now,"analysis_unchanged":True})
    metadata=SubmissionMetadata.from_dict(corrections.get("metadata", {}))
    out.mkdir(parents=True, exist_ok=True)
    version_dir=out / "submission_v2"; version_dir.mkdir(parents=True, exist_ok=True)
    provenance={"input_sha256":digest,"genome_id":genome_id,"analysis_run_id":manifest.get("created_at"),"submission_version":"v2","parent_submission_version":"v1","correction_file_sha256":hashlib.sha256(corr_path.read_bytes()).hexdigest(),"corrections_applied":len(audit),"analysis_unchanged":True}
    validation=write_package(version_dir, genome_id, genome, proteins, provenance, metadata=metadata, feature_overrides=overrides)
    (out/"submission_correction_audit.json").write_text(json.dumps(audit,indent=2,sort_keys=True)+"\n")
    if audit:
        fields=list(audit[0]); (out/"submission_correction_audit.tsv").write_text("\t".join(fields)+"\n"+"\n".join("\t".join(str(r.get(f,"")) for f in fields) for r in audit)+"\n")
    (out/"submission_revision_manifest.json").write_text(json.dumps(provenance,indent=2,sort_keys=True)+"\n")
    return validation
