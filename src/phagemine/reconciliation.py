"""Observational comparison of PHANOTATE and Prodigal gene models."""
from __future__ import annotations
import csv, json, shutil, subprocess, tempfile, os
from dataclasses import asdict, dataclass
from pathlib import Path
from .gene_prediction import reverse_complement
from .genome import translate

@dataclass(frozen=True)
class GeneModel:
    caller: str; identifier: str; start: int; end: int; strand: str
    sequence: str = ""; frame: int | None = None; caller_version: str | None = None
    command: list[str] | None = None; options: dict | None = None

def _overlap(a: GeneModel, b: GeneModel) -> int:
    return max(0, min(a.end,b.end)-max(a.start,b.start)+1)

def _classify(p: GeneModel|None, d: GeneModel|None, overlap: int) -> str:
    if p is None: return "PRODIGAL_ONLY"
    if d is None: return "PHANOTATE_ONLY"
    if p.strand != d.strand: return "STRAND_DISCORDANCE"
    if p.start == d.start and p.end == d.end: return "EXACT_CONCORDANCE"
    if p.end == d.end: return "START_DISCORDANCE"
    if p.start == d.start: return "STOP_DISCORDANCE"
    return "START_AND_STOP_DISCORDANCE"

def reconcile_models(phanotate: list[GeneModel], prodigal: list[GeneModel], genome_sha256: str | None = None) -> list[dict]:
    rows=[]; used_d=set(); locus=0
    for p in sorted(phanotate,key=lambda x:(x.start,x.end,x.identifier)):
        candidates=[(i,d,_overlap(p,d)) for i,d in enumerate(prodigal) if _overlap(p,d)>0]
        substantial=[(i,d,o) for i,d,o in candidates if o/max(p.end-p.start+1,d.end-d.start+1)>=0.5]
        if len(substantial)==1:
            i,d,o=substantial[0]; used_d.add(i); conflict=_classify(p,d,o)
        elif len(substantial)>1:
            locus+=1; rows.append(_row(locus,p,None,0,"SPLIT_MODEL" if len(substantial)>1 else "COMPLEX_CONFLICT",genome_sha256)); used_d.update(i for i,_,_ in substantial); continue
        else:
            d=None; o=0; conflict="PHANOTATE_ONLY"
        locus+=1; rows.append(_row(locus,p,d,o,conflict,genome_sha256))
    for i,d in enumerate(sorted(prodigal,key=lambda x:(x.start,x.end,x.identifier))):
        if i not in used_d:
            locus+=1; rows.append(_row(locus,None,d,0,"PRODIGAL_ONLY",genome_sha256))
    return rows

def _row(locus,p,d,o,conflict,sha):
    start_agreement=bool(p and d and p.start==d.start); stop_agreement=bool(p and d and p.end==d.end)
    return {"locus_id":f"LOCUS_{locus:06d}","conflict_type":conflict,"phanotate_id":p.identifier if p else None,"phanotate_start":p.start if p else None,"phanotate_end":p.end if p else None,"phanotate_strand":p.strand if p else None,"prodigal_id":d.identifier if d else None,"prodigal_start":d.start if d else None,"prodigal_end":d.end if d else None,"prodigal_strand":d.strand if d else None,"overlap_bp":o,"phanotate_overlap_fraction":o/(p.end-p.start+1) if p and o else 0.0,"prodigal_overlap_fraction":o/(d.end-d.start+1) if d and o else 0.0,"start_agreement":start_agreement,"stop_agreement":stop_agreement,"orf_existence_status":"CONCORDANT" if p and d else "CALLER_SPECIFIC","start_status":"START_CONFIRMED" if start_agreement else ("START_DISCORDANT" if p and d else "START_NOT_COMPARABLE"),"stop_status":"STOP_CONFIRMED" if stop_agreement else ("STOP_DISCORDANT" if p and d else "STOP_NOT_COMPARABLE"),"selected_model":"PHANOTATE","decision_status":"OBSERVATIONAL_NOT_ADJUDICATED" if conflict not in {"EXACT_CONCORDANCE"} else "OBSERVATIONAL_CONCORDANT","manual_review":conflict not in {"EXACT_CONCORDANCE"},"input_genome_sha256":sha}

def write_reconciliation(root: str|Path, rows: list[dict], provenance: dict) -> None:
    root=Path(root); (root/"orf_reconciliation.json").write_text(json.dumps({"provenance":provenance,"records":rows},indent=2,sort_keys=True))
    cols=list(rows[0]) if rows else ["locus_id","conflict_type"]
    with (root/"orf_reconciliation.tsv").open("w",newline="") as h:
        w=csv.DictWriter(h,fieldnames=cols,delimiter="\t"); w.writeheader(); w.writerows(rows)

class ProdigalPredictor:
    name="Prodigal"
    def __init__(self, executable: str|None=None): self.executable=executable or shutil.which("prodigal")
    def predict(self, fasta: str|Path, sequence: str) -> list[GeneModel]:
        if not self.executable: raise RuntimeError("Prodigal is required for --reconcile-orfs")
        handle=tempfile.NamedTemporaryFile(prefix="phagemine_prodigal_",suffix=".gff",delete=False); gff=Path(handle.name); handle.close(); cmd=[self.executable,"-i",str(fasta),"-o",str(gff),"-f","gff","-q"]
        result=subprocess.run(cmd,capture_output=True,text=True,check=False)
        if result.returncode:
            gff.unlink(missing_ok=True)
            raise RuntimeError(f"Prodigal failed (exit {result.returncode}): {(result.stderr or result.stdout).strip()}")
        models=[]
        for line in gff.read_text().splitlines():
            if line.startswith("#"): continue
            f=line.split("\t");
            if len(f)<9 or f[2]!="CDS": continue
            attrs=dict(x.split("=",1) for x in f[8].split(";") if "=" in x); s,e=int(f[3]),int(f[4]); strand=f[6]
            models.append(GeneModel("Prodigal",attrs.get("ID",f"PRODIGAL_{len(models)+1}"),s,e,strand,caller_version="reported",command=cmd,options={"format":"gff"}))
        gff.unlink(missing_ok=True)
        return models
