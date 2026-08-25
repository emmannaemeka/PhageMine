"""Conservative bacteriophage hallmark-system checks.

Failure to detect a hallmark is reported as ``NOT_ESTABLISHED`` rather than
biological absence.  Phage architectures differ and database searches can miss
divergent proteins.
"""
from __future__ import annotations

import csv
import json
from pathlib import Path


HALLMARK_RULES_VERSION = "1.0"
HALLMARKS = {
    "major_capsid": ("major capsid", "major head protein"),
    "portal": ("portal protein", "portal vertex"),
    "terminase": ("terminase", "dna packaging atpase"),
    "tail_structure": ("tail protein", "tail fiber", "tail fibre", "tail spike", "tailspike", "baseplate"),
    "tape_measure": ("tape measure", "tape-measure"),
    "dna_replication": ("dna polymerase", "helicase", "primase", "replication protein"),
    "lysis": ("holin", "endolysin", "lysin", "spanin", "cell wall hydrolase", "peptidoglycan hydrolase", "amidase"),
}


def assess_hallmarks(classifications: list[dict]) -> list[dict]:
    rows=[]
    for hallmark, keywords in HALLMARKS.items():
        matches=[]
        for item in classifications:
            text=" ".join(str(item.get(key) or "") for key in ("display_product","proposed_function","domain_summary","functional_category")).lower()
            if any(keyword in text for keyword in keywords):
                matches.append(item)
        confident=[item for item in matches if item.get("functional_state") in {"KNOWN_FUNCTION","PROBABLE_FUNCTION"}]
        selected=confident or matches
        status="DETECTED" if confident else ("POSSIBLE_DOMAIN_OR_CATEGORY_SUPPORT" if matches else "NOT_ESTABLISHED")
        rows.append({
            "hallmark": hallmark,
            "status": status,
            "protein_ids": ";".join(item["protein_id"] for item in selected),
            "proposed_functions": "; ".join(str(item.get("display_product") or item.get("proposed_function") or "") for item in selected),
            "highest_confidence": next((level for level in ("HIGH","MODERATE","LOW") if any(item.get("confidence")==level for item in selected)), "NONE"),
            "interpretation": "Computational evidence supports this component." if status=="DETECTED" else ("Broad evidence is present but does not establish a specific product." if matches else "Not established by the current annotation; this is not evidence of biological absence."),
            "rules_version": HALLMARK_RULES_VERSION,
        })
    return rows


def write_hallmarks(root: str|Path, rows: list[dict]) -> None:
    root=Path(root)
    (root/"hallmark_completeness.json").write_text(json.dumps(rows,indent=2,sort_keys=True)+"\n")
    columns=["hallmark","status","protein_ids","proposed_functions","highest_confidence","interpretation","rules_version"]
    with (root/"hallmark_completeness.tsv").open("w",newline="") as handle:
        writer=csv.DictWriter(handle,fieldnames=columns,delimiter="\t"); writer.writeheader(); writer.writerows(rows)
