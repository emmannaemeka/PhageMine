"""Conservative, observational adjudication of competing ORF models."""
from __future__ import annotations
import csv, json
from pathlib import Path

def _support(items):
    return [x for x in items if x.get("supports") and x.get("evidence_strength") in {"STRONG", "EXPERIMENTAL"}]

def adjudicate(rows: list[dict], evidence_by_candidate: dict[str, list[dict]] | None = None) -> list[dict]:
    evidence_by_candidate = evidence_by_candidate or {}
    out=[]
    for row in rows:
        p=row.get("phanotate_id"); d=row.get("prodigal_id"); pe=evidence_by_candidate.get(p,[]); de=evidence_by_candidate.get(d,[])
        ps=_support(pe); ds=_support(de); conflict=row["conflict_type"]
        start_decision="START_NOT_COMPARABLE"; stop_decision="STOP_NOT_COMPARABLE"; strand_decision="STRAND_AMBIGUOUS"
        recommended="PHANOTATE"; status="INSUFFICIENT_EVIDENCE"; reason="No reliable model-specific evidence was available."
        comparable = bool(p and d and pe and de)
        if p and d and not comparable and conflict != "EXACT_CONCORDANCE":
            out.append({**row,"phanotate_orf_support":"SUPPORTED" if ps else "UNRESOLVED","prodigal_orf_support":"SUPPORTED" if ds else "UNRESOLVED","start_decision":"START_AMBIGUOUS" if conflict.startswith("START") else "START_NOT_COMPARABLE","stop_decision":"STOP_AMBIGUOUS" if conflict.startswith("STOP") else "STOP_NOT_COMPARABLE","strand_decision":"STRAND_AMBIGUOUS","evidence_agreement":"INSUFFICIENT_COMPARATIVE_EVIDENCE","recommended_model":"PHANOTATE","adjudication_status":"INSUFFICIENT_COMPARATIVE_EVIDENCE","reason":"Competing models did not receive comparable model-specific evidence access.","manual_review":True}); continue
        if conflict=="EXACT_CONCORDANCE":
            status="CONCORDANT_EVIDENCE"; reason="Both callers produced identical coordinates and strand."
            start_decision="START_CONFIRMED"; stop_decision="STOP_CONFIRMED"; strand_decision="STRAND_CONFIRMED"
            recommended="PHANOTATE"
        elif conflict=="START_DISCORDANCE":
            start_decision="START_RESOLVED_PHANOTATE" if ps and not ds else ("START_RESOLVED_PRODIGAL" if ds and not ps else "START_AMBIGUOUS")
            status="PARTIALLY_CONCORDANT_EVIDENCE" if start_decision!="START_AMBIGUOUS" else "INSUFFICIENT_EVIDENCE"; recommended="PHANOTATE" if start_decision!="START_RESOLVED_PRODIGAL" else "PRODIGAL"; reason="Strong model-specific evidence supports the selected start." if status!="INSUFFICIENT_EVIDENCE" else "Start boundary evidence is insufficient."
        elif conflict=="STOP_DISCORDANCE":
            stop_decision="STOP_RESOLVED_PHANOTATE" if ps and not ds else ("STOP_RESOLVED_PRODIGAL" if ds and not ps else "STOP_AMBIGUOUS")
            status="PARTIALLY_CONCORDANT_EVIDENCE" if stop_decision!="STOP_AMBIGUOUS" else "INSUFFICIENT_EVIDENCE"; recommended="PHANOTATE" if stop_decision!="STOP_RESOLVED_PRODIGAL" else "PRODIGAL"; reason="Strong model-specific evidence supports the selected stop." if status!="INSUFFICIENT_EVIDENCE" else "Stop boundary evidence is insufficient."
        elif conflict=="STRAND_DISCORDANCE":
            if ps and not ds: strand_decision="STRAND_SUPPORTS_PHANOTATE"; recommended="PHANOTATE"; status="CONCORDANT_EVIDENCE"; reason="Strong evidence supports the PHANOTATE orientation."
            elif ds and not ps: strand_decision="STRAND_SUPPORTS_PRODIGAL"; recommended="PRODIGAL"; status="CONCORDANT_EVIDENCE"; reason="Strong evidence supports the Prodigal orientation."
            else: reason="Opposite-strand models lack decisive independent evidence."
        elif conflict=="PRODIGAL_ONLY":
            if ds: status="CONCORDANT_EVIDENCE"; recommended="PRODIGAL"; reason="Independent strong evidence supports the Prodigal-only ORF."; row["prodigal_orf_support"]="PRODIGAL_ONLY_RESCUE_CANDIDATE"
            else: recommended="UNRESOLVED"; status="INSUFFICIENT_EVIDENCE"; reason="No independent evidence supports this Prodigal-only ORF."; row["prodigal_orf_support"]="PRODIGAL_ONLY_UNRESOLVED"
        elif conflict=="PHANOTATE_ONLY": row["phanotate_orf_support"]="PHANOTATE_ONLY_SUPPORTED" if ps else "PHANOTATE_ONLY_UNRESOLVED"
        coordinate_concordance = conflict == "EXACT_CONCORDANCE"
        evidence_agreement = ("COORDINATE_CONCORDANCE" if coordinate_concordance else
                              ("CONCORDANT_EVIDENCE" if ps and ds else
                               ("PARTIALLY_CONCORDANT_EVIDENCE" if ps or ds else "INSUFFICIENT_EVIDENCE")))
        out.append({**row,"phanotate_orf_support":row.get("phanotate_orf_support", "SUPPORTED" if ps else "UNRESOLVED"),"prodigal_orf_support":row.get("prodigal_orf_support", "SUPPORTED" if ds else "UNRESOLVED"),"start_decision":start_decision,"stop_decision":stop_decision,"strand_decision":strand_decision,"evidence_agreement":evidence_agreement,"phanotate_best_evidence":ps[0] if ps else None,"prodigal_best_evidence":ds[0] if ds else None,"supporting_evidence":{"PHANOTATE":pe,"PRODIGAL":de},"contradictory_evidence":{},"recommended_model":recommended,"adjudication_status":status,"reason":reason,"manual_review":False if coordinate_concordance else status=="INSUFFICIENT_EVIDENCE"})
    return out

def write_adjudication(root: str|Path, records: list[dict], provenance: dict|None=None):
    root=Path(root); (root/"orf_adjudication.json").write_text(json.dumps({"provenance":provenance or {},"records":records},indent=2,sort_keys=True,default=str))
    cols=["locus_id","conflict_type","phanotate_id","prodigal_id","phanotate_orf_support","prodigal_orf_support","start_decision","stop_decision","strand_decision","evidence_agreement","recommended_model","adjudication_status","reason","manual_review"]
    with (root/"orf_adjudication.tsv").open("w",newline="") as h:
        w=csv.DictWriter(h,fieldnames=cols,delimiter="\t"); w.writeheader(); w.writerows({k:r.get(k) for k in cols} for r in records)
