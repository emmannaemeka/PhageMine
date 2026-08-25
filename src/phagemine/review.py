"""Build a focused, auditable manual-annotation review queue."""
from __future__ import annotations

import csv
import json
from pathlib import Path


def build_annotation_review(proteins, classifications: list[dict], gene_calls: list[dict]|None=None, hallmarks: list[dict]|None=None) -> list[dict]:
    proteins_by_id={protein.protein_id:protein for protein in proteins}
    rows=[]
    for item in classifications:
        state=item.get("functional_state"); flags=item.get("ambiguity_flags") or []
        if state not in {"UNRESOLVED","CONFLICTING_EVIDENCE"} and not flags:
            continue
        protein=proteins_by_id[item["protein_id"]]
        issue="; ".join(flags) if flags else ("Accepted evidence supports incompatible functions." if state=="CONFLICTING_EVIDENCE" else "No reliable functional assignment was established.")
        rows.append({"priority":"HIGH" if state=="CONFLICTING_EVIDENCE" else "MEDIUM","review_type":"FUNCTION_ASSIGNMENT","protein_id":protein.protein_id,"coordinates":f"{protein.start}..{protein.end}","strand":protein.strand,"proposed_function":item.get("display_product"),"classification":item.get("display_classification"),"issue":issue,"recommended_action":"Inspect the protein evidence page and compare conserved neighbourhood and independent annotation outputs."})
    for item in gene_calls or []:
        if item.get("review_flag")=="NONE": continue
        protein=proteins_by_id.get(item.get("phanotate_id"))
        rows.append({"priority":"HIGH" if item.get("review_flag") in {"POSSIBLE_FALSE_CALL","REVIEW_STRAND"} else "MEDIUM","review_type":"GENE_CALL","protein_id":item.get("protein_id"),"coordinates":f"{protein.start}..{protein.end}" if protein else "alternative Prodigal model","strand":protein.strand if protein else "","proposed_function":protein.annotation if protein else "","classification":item.get("gene_call_confidence"),"issue":item.get("rationale"),"recommended_action":"Compare both translated models, caller support, overlaps and model-specific evidence; do not delete automatically."})
    for item in hallmarks or []:
        if item.get("status")=="DETECTED": continue
        rows.append({"priority":"HIGH" if item.get("hallmark") in {"major_capsid","portal","terminase"} else "LOW","review_type":"HALLMARK_NOT_ESTABLISHED","protein_id":"","coordinates":"","strand":"","proposed_function":item.get("hallmark"),"classification":item.get("status"),"issue":item.get("interpretation"),"recommended_action":"Inspect unresolved proteins and independent phage annotations; absence is not established."})
    priority={"HIGH":0,"MEDIUM":1,"LOW":2}
    return sorted(rows,key=lambda row:(priority[row["priority"]],row["review_type"],row["protein_id"],row["proposed_function"] or ""))


def write_annotation_review(root: str|Path, rows: list[dict]) -> None:
    root=Path(root)
    (root/"annotation_review.json").write_text(json.dumps(rows,indent=2,sort_keys=True)+"\n")
    columns=["priority","review_type","protein_id","coordinates","strand","proposed_function","classification","issue","recommended_action"]
    with (root/"annotation_review.tsv").open("w",newline="") as handle:
        writer=csv.DictWriter(handle,fieldnames=columns,delimiter="\t"); writer.writeheader(); writer.writerows(rows)
