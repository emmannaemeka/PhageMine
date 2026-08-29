"""Observational comparison of PHANOTATE and Prodigal gene models."""
from __future__ import annotations
import csv, json, shutil, subprocess, tempfile, os
from pathlib import Path
from .gene_prediction import reverse_complement
from .genome import translate
from .gene_models import GeneModel

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
    prodigal=sorted(prodigal,key=lambda x:(x.start,x.end,x.identifier))
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
    for i,d in enumerate(prodigal):
        if i not in used_d:
            locus+=1; rows.append(_row(locus,None,d,0,"PRODIGAL_ONLY",genome_sha256))
    return rows

def _row(locus,p,d,o,conflict,sha):
    start_agreement=bool(p and d and p.start==d.start); stop_agreement=bool(p and d and p.end==d.end)
    return {"locus_id":f"LOCUS_{locus:06d}","conflict_type":conflict,"phanotate_id":p.identifier if p else None,"phanotate_start":p.start if p else None,"phanotate_end":p.end if p else None,"phanotate_strand":p.strand if p else None,"prodigal_id":d.identifier if d else None,"prodigal_start":d.start if d else None,"prodigal_end":d.end if d else None,"prodigal_strand":d.strand if d else None,"overlap_bp":o,"phanotate_overlap_fraction":o/(p.end-p.start+1) if p and o else 0.0,"prodigal_overlap_fraction":o/(d.end-d.start+1) if d and o else 0.0,"start_agreement":start_agreement,"stop_agreement":stop_agreement,"orf_existence_status":"CONCORDANT" if p and d else "CALLER_SPECIFIC","start_status":"START_CONFIRMED" if start_agreement else ("START_DISCORDANT" if p and d else "START_NOT_COMPARABLE"),"stop_status":"STOP_CONFIRMED" if stop_agreement else ("STOP_DISCORDANT" if p and d else "STOP_NOT_COMPARABLE"),"selected_model":"PHANOTATE","primary_model_policy":"PHANOTATE_RETAINED_WITHOUT_TRUTH_ADJUDICATION","model_resolution":"CALLERS_AGREE" if conflict == "EXACT_CONCORDANCE" else "NOT_RESOLVED","decision_status":"OBSERVATIONAL_NOT_ADJUDICATED" if conflict not in {"EXACT_CONCORDANCE"} else "OBSERVATIONAL_CONCORDANT","manual_review":conflict not in {"EXACT_CONCORDANCE"},"input_genome_sha256":sha}

def write_reconciliation(root: str|Path, rows: list[dict], provenance: dict) -> None:
    root=Path(root); (root/"orf_reconciliation.json").write_text(json.dumps({"provenance":provenance,"records":rows},indent=2,sort_keys=True))
    cols=list(rows[0]) if rows else ["locus_id","conflict_type"]
    with (root/"orf_reconciliation.tsv").open("w",newline="") as h:
        w=csv.DictWriter(h,fieldnames=cols,delimiter="\t"); w.writeheader(); w.writerows(rows)

def gene_call_review(rows: list[dict], evidence_by_candidate: dict[str, list[dict]], protein_lengths: dict[str, int]) -> list[dict]:
    """Summarise caller agreement without silently deleting predicted CDSs.

    ``POSSIBLE_FALSE_CALL`` is deliberately a review flag rather than an
    automatic rejection: short, caller-specific phage genes can be real.
    """
    reviewed=[]
    for row in rows:
        phanotate_id=row.get("phanotate_id"); prodigal_id=row.get("prodigal_id")
        protein_id=phanotate_id or prodigal_id
        evidence=evidence_by_candidate.get(protein_id, [])
        accepted=[item for item in evidence if item.get("supports")]
        strong=[item for item in accepted if item.get("evidence_strength") in {"STRONG", "EXPERIMENTAL"}]
        length=protein_lengths.get(protein_id)
        conflict=row["conflict_type"]
        if conflict == "EXACT_CONCORDANCE":
            confidence, flag, rationale = "HIGH", "NONE", "PHANOTATE and Prodigal agree on both CDS boundaries and strand."
        elif conflict in {"START_DISCORDANCE", "STOP_DISCORDANCE", "START_AND_STOP_DISCORDANCE"}:
            confidence, flag, rationale = "MODERATE", "REVIEW_BOUNDARIES", "Both callers detect an overlapping CDS but disagree on one or both boundaries."
        elif conflict == "STRAND_DISCORDANCE":
            confidence, flag, rationale = "LOW", "REVIEW_STRAND", "The callers place the overlapping CDS on opposite strands."
        elif strong:
            confidence, flag, rationale = "MODERATE", "REVIEW_CALLER_SPECIFIC", "Only one caller predicts this CDS, but independent strong functional evidence supports a translated product."
        elif conflict == "PHANOTATE_ONLY" and length is not None and length < 40:
            confidence, flag, rationale = "LOW", "POSSIBLE_FALSE_CALL", "Short PHANOTATE-only CDS has no accepted strong functional evidence; retain pending manual review."
        else:
            confidence, flag, rationale = "LOW", "REVIEW_CALLER_SPECIFIC", "Only one caller predicts this CDS and strong independent support is absent."
        best = strong[0] if strong else (accepted[0] if accepted else None)
        reviewed.append({
            "locus_id": row["locus_id"], "protein_id": protein_id,
            "phanotate_id": phanotate_id, "prodigal_id": prodigal_id,
            "conflict_type": conflict, "length_aa": length,
            "accepted_evidence_count": len(accepted),
            "strong_evidence_count": len(strong),
            "best_evidence": (f"{best.get('source')}:{best.get('identifier') or best.get('family_name') or 'match'}" if best else "No accepted evidence"),
            "gene_call_confidence": confidence, "review_flag": flag,
            "confidence_calibrated": False,
            "confidence_interpretation": "Rule-based caller agreement category; not an empirical probability that the ORF is real.",
            "rationale": rationale,
        })
    return reviewed

def write_gene_call_review(root: str|Path, records: list[dict]) -> None:
    root=Path(root)
    columns=["locus_id","protein_id","phanotate_id","prodigal_id","conflict_type","length_aa","accepted_evidence_count","strong_evidence_count","best_evidence","gene_call_confidence","confidence_calibrated","confidence_interpretation","review_flag","rationale"]
    with (root/"gene_call_confidence.tsv").open("w", newline="") as handle:
        writer=csv.DictWriter(handle, fieldnames=columns, delimiter="\t"); writer.writeheader(); writer.writerows(records)
    flagged=[record for record in records if record["review_flag"] != "NONE"]
    with (root/"gene_calls_for_review.tsv").open("w", newline="") as handle:
        writer=csv.DictWriter(handle, fieldnames=columns, delimiter="\t"); writer.writeheader(); writer.writerows(flagged)

class ProdigalPredictor:
    name="Prodigal"
    def __init__(self, executable: str|None=None):
        self.executable=executable or shutil.which("prodigal")
        self.last_command: list[str] | None = None
        self.last_stdout = ""
        self.last_stderr = ""
        self.last_raw_gff = ""
    def version(self) -> str:
        if not self.executable:
            return "unavailable"
        result = subprocess.run([self.executable, "-v"], capture_output=True, text=True, check=False)
        return ((result.stdout or result.stderr).strip().splitlines() or ["reported"])[0]
    def parameters(self) -> dict:
        return {"executable": self.executable, "format": "gff", "quiet": True}
    def predict(self, fasta: str|Path, sequence: str) -> list[GeneModel]:
        if not self.executable: raise RuntimeError("Prodigal is required for --reconcile-orfs")
        handle=tempfile.NamedTemporaryFile(prefix="phagemine_prodigal_",suffix=".gff",delete=False); gff=Path(handle.name); handle.close(); cmd=[self.executable,"-i",str(fasta),"-o",str(gff),"-f","gff","-q"]
        self.last_command = cmd
        result=subprocess.run(cmd,capture_output=True,text=True,check=False)
        self.last_stdout, self.last_stderr = result.stdout or "", result.stderr or ""
        if result.returncode:
            gff.unlink(missing_ok=True)
            raise RuntimeError(f"Prodigal failed (exit {result.returncode}): {(result.stderr or result.stdout).strip()}")
        self.last_raw_gff = gff.read_text()
        models=[]
        for line in self.last_raw_gff.splitlines():
            if line.startswith("#"): continue
            f=line.split("\t");
            if len(f)<9 or f[2]!="CDS": continue
            attrs=dict(x.split("=",1) for x in f[8].split(";") if "=" in x); s,e=int(f[3]),int(f[4]); strand=f[6]
            models.append(GeneModel("Prodigal",attrs.get("ID",f"PRODIGAL_{len(models)+1}"),s,e,strand,caller_version="reported",command=cmd,options={"format":"gff"}))
        gff.unlink(missing_ok=True)
        return models
