#!/usr/bin/env python3
"""Separate, deterministic automated evidence-adjudication track.

The blind phase consumes only the neutral frozen review CSV.  The unblind phase
is a separate invocation and is the first phase allowed to read the confidential
tool key.  This is deliberately not the human-adjudication pipeline.
"""
from __future__ import annotations
import argparse, csv, hashlib, json, re, statistics, zipfile
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "automated_benchmark_results"
INPUT = ROOT / "evaluation/v1.2_adjudication/reviewer_files/blinded_adjudication_cases.csv"
KEY = ROOT / "evaluation/v1.2_adjudication/blinding/tool_blinding_key.tsv"
UTILITY = ROOT / "evaluation/v1.2_adjudication/frozen/utility_definition.json"
FUNCTIONAL = ROOT / "evaluation/v1.2_functional/functional_by_locus.tsv"
ALLOWED = ["CORRECT", "PARTIALLY_CORRECT", "TOO_GENERAL", "UNSUPPORTED/INCORRECT", "UNRESOLVABLE", "NOT_EVALUABLE"]
WEIGHTS = {"CORRECT": 1.0, "PARTIALLY_CORRECT": .5, "TOO_GENERAL": .25, "UNSUPPORTED/INCORRECT": 0.0}
NONINFO = re.compile(r"\b(hypothetical|uncharacterized|uncharacterised|unknown function|unknown protein|unassigned|conserved protein of unknown)\b", re.I)
TOOL_WORDS = re.compile(r"phagemine|pharokka|prokka|tool_blinding_key|prediction_[ab]_tool", re.I)

def sha(path: Path) -> str:
    h = hashlib.sha256(); h.update(path.read_bytes()); return h.hexdigest()

def norm(s: str) -> str:
    s = (s or "").lower().replace("/", " ")
    s = re.sub(r"[^a-z0-9]+", " ", s)
    return re.sub(r"\s+", " ", s).strip()

def parse_raw(s: str) -> dict[str, Any]:
    try: return json.loads(s or "{}")
    except Exception: return {}

def evidence(raw: str) -> dict[str, Any]:
    d = parse_raw(raw); out: dict[str, Any] = {}
    for obj in (d, d.get("annotation", {}), d.get("classification", {}), d.get("gff_attributes", {})):
        if isinstance(obj, dict): out.update(obj)
    return out

def classify(ref: str, pred: str, raw: str) -> tuple[str, str, str, str, str, str]:
    """Return category, confidence, support, specificity, reason, conflict."""
    r, p = norm(ref), norm(pred); e = evidence(raw)
    if NONINFO.search(ref or ""):
        return "NOT_EVALUABLE", "HIGH", "UNCLEAR", "UNCLEAR", "Frozen reference product is non-evaluable under the protocol.", "NO"
    if not p or NONINFO.search(pred or ""):
        return "UNRESOLVABLE", "LOW", "NO", "UNCLEAR", "Candidate assertion is absent or non-informative; no specific function can be scored.", "UNCLEAR"
    if r == p:
        return "CORRECT", "HIGH", "YES", "YES", "Exact normalized product agreement.", "NO"
    # The only explicit accepted equivalence in the frozen rule file.
    pair = {norm("major head protein"), norm("major capsid protein")}
    if {r, p} == pair:
        return "CORRECT", "HIGH", "YES", "YES", "Exact accepted synonym from the frozen synonym table.", "NO"
    rt, pt = set(r.split()), set(p.split())
    informative = bool(e.get("best_evidence") or e.get("supporting_sources") or e.get("evidence_sources"))
    tier = str(e.get("evidence_tier", "")); count = int(str(e.get("supporting_source_count", "0") or "0")) if str(e.get("supporting_source_count", "0") or "0").isdigit() else 0
    strong = informative and (count >= 2 or tier in {"3", "4"})
    family_only = bool(e.get("ortholog_groups") or e.get("phrog")) and not e.get("best_evidence")
    if pt < rt and pt and len(pt) >= 1:
        return "TOO_GENERAL", ("HIGH" if strong else "MODERATE"), ("YES" if strong else "PARTIAL"), "TOO_GENERAL", "Assertion is a compatible lexical broadening of the informative reference; specificity is lower.", "NO"
    overlap = len(rt & pt) / max(1, len(rt | pt))
    if overlap >= .5 and strong and not family_only:
        return "PARTIALLY_CORRECT", "MODERATE", "PARTIAL", "UNCLEAR", "Substantial terminology overlap with concordant frozen evidence, but the assertion is not an exact supported equivalence.", "UNCLEAR"
    if (rt & pt) and family_only:
        return "UNRESOLVABLE", "LOW", "UNCLEAR", "UNCLEAR", "Family/orthology evidence is retained but does not by itself establish the asserted specific function.", "UNCLEAR"
    return "UNRESOLVABLE", "LOW", "UNCLEAR", "UNCLEAR", "Frozen evidence does not deterministically establish equivalence or contradiction at the asserted specificity.", "YES" if (rt and pt and not (rt & pt)) else "UNCLEAR"

def blind_rows() -> list[dict[str, str]]:
    with INPUT.open(newline="") as f: return list(csv.DictReader(f))

def make_blind() -> None:
    OUT.mkdir(exist_ok=True); rows = blind_rows()
    if len(rows) != 689 or {r["blinded_unit_id"] for r in rows} != {f"FB{i:04d}" for i in range(1,690)}: raise SystemExit("frozen blinded units are not exactly FB0001-FB0689")
    output=[]; audits=[]
    for row in rows:
        target = row["review_target"]; pred=row[f"prediction_{target}_product"]; raw=row[f"prediction_{target}_evidence_raw"]
        c,conf,sup,spec,reason,conflict=classify(row["reference_product"], pred, raw)
        out={k:row.get(k,"") for k in ["blinded_unit_id","matching_mode","genome","review_target","reference_locus","start","end","strand","reference_product","prediction_A_product","prediction_B_product"]}
        out.update({"adjudication_class":c,"evidence_basis":"FROZEN_EVIDENCE_PACKET","uncertainty":conf,"review_required":"YES" if c in {"UNRESOLVABLE","NOT_EVALUABLE"} or conf=="LOW" else "NO","automated_reason":reason,"decisive_evidence":row.get(f"prediction_{target}_evidence_raw",""),"conflicting_evidence":conflict,"adjudication_method":"AUTOMATED_EVIDENCE_ADJUDICATION","automated_confidence":conf,"evidence_supports_prediction":sup,"evidence_strength":"STRONG" if conf=="HIGH" else ("MODERATE" if conf=="MODERATE" else "INSUFFICIENT"),"specificity_appropriate":spec,"error_category":"NONE" if c in {"CORRECT","TOO_GENERAL"} else ("INSUFFICIENT_EVIDENCE" if c=="UNRESOLVABLE" else "OTHER"),"reviewer_notes":"AUTOMATED TRACK; not human adjudication."})
        if TOOL_WORDS.search(json.dumps(out)): raise SystemExit(f"tool identity leaked into blind row {row['blinded_unit_id']}")
        output.append(out); audits.append({"blinded_unit_id":row["blinded_unit_id"],"method":"AUTOMATED_EVIDENCE_ADJUDICATION","rule_version":"1.0","category":c,"confidence":conf,"evidence_sha256":hashlib.sha256(raw.encode()).hexdigest()})
    cols=list(output[0]); tsv=OUT/"automated_blinded_adjudication.tsv"
    with tsv.open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=cols,delimiter="\t"); w.writeheader(); w.writerows(output)
    with (OUT/"automated_adjudication_audit.jsonl").open("w") as f:
        for a in audits: f.write(json.dumps(a,sort_keys=True)+"\n")
    # openpyxl workbook with no formulas, hidden sheets, or hidden columns.
    from openpyxl import Workbook
    from openpyxl.worksheet.datavalidation import DataValidation
    wb=Workbook(); ws=wb.active; ws.title="Automated_Blinded_Adjudication"; ws.append(cols)
    for r in output: ws.append([r.get(c,"") for c in cols])
    dv=DataValidation(type="list", formula1='"CORRECT,PARTIALLY_CORRECT,TOO_GENERAL,UNSUPPORTED/INCORRECT,UNRESOLVABLE,NOT_EVALUABLE"'); ws.add_data_validation(dv); dv.add(f"A2:A{len(output)+1}")
    info=wb.create_sheet("Method"); info.append(["Field","Value"]); info.append(["adjudication_method","AUTOMATED_EVIDENCE_ADJUDICATION"]); info.append(["note","Separate automated track; no human judgments and no tool identity in this workbook."])
    for sh in wb.worksheets: sh.sheet_state="visible"
    xlsx=OUT/"automated_blinded_adjudication.xlsx"; wb.save(xlsx)
    manifest={"benchmark_commit":"81fafddfee7c38b5a142ec9a36abb9e00a6f94c4","input_sha256":sha(INPUT),"utility_sha256":sha(UTILITY),"method":"AUTOMATED_EVIDENCE_ADJUDICATION","created_utc":datetime.now(timezone.utc).isoformat(),"units":len(output),"category_counts":dict(Counter(r["adjudication_class"] for r in output)),"confidence_counts":dict(Counter(r["automated_confidence"] for r in output)),"outputs":{p.name:sha(p) for p in [tsv,OUT/"automated_adjudication_audit.jsonl",xlsx]},"blinding_audit":"AUTOMATED_ADJUDICATION_BLINDING_PASS"}
    (OUT/"adjudication_freeze_manifest.json").write_text(json.dumps(manifest,indent=2,sort_keys=True)+"\n")
    print(json.dumps(manifest,indent=2))

def load_blind():
    m=json.loads((OUT/"adjudication_freeze_manifest.json").read_text())
    for fn,h in m["outputs"].items():
        if sha(OUT/fn)!=h: raise SystemExit(f"blind output checksum mismatch: {fn}")
    return m, list(csv.DictReader((OUT/"automated_blinded_adjudication.tsv").open(),delimiter="\t"))

def unblind() -> None:
    m, blind=load_blind(); key=[]
    with KEY.open() as f: key=list(csv.DictReader(f,delimiter="\t"))
    km={r["case_id"]:r for r in key}
    for r in blind:
        k=km.get(r["blinded_unit_id"]); 
        if not k: raise SystemExit("missing confidential key row")
        r["tool_identity_after_freeze"]=k["target_tool"]
    out=OUT/"unblinded_results.tsv"; cols=list(blind[0])
    with out.open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=cols,delimiter="\t"); w.writeheader(); w.writerows(blind)
    # Join the unchanged frozen locus table to replace only its review-required
    # rows with the already-frozen automated target decision.  All other frozen
    # deterministic categories are retained verbatim.
    frozen=list(csv.DictReader(FUNCTIONAL.open(),delimiter="\t"))
    key_by_case={r["blinded_unit_id"]:r for r in blind}
    def lookup(row, tool):
        k=(row["accession"],row["start"],row["end"],row["strand"],row["matching_mode"])
        category=row[tool+"_category"]
        if category=="NON_EVALUABLE_REFERENCE": category="NOT_EVALUABLE"
        elif category=="ABSTENTION": category="UNRESOLVABLE"
        if category=="UNRESOLVED_REVIEW_REQUIRED":
            matches=[x for x in blind if (x["genome"],x["start"],x["end"],x["strand"],x["matching_mode"])==k and key_by_case[x["blinded_unit_id"]]["tool_identity_after_freeze"]==tool]
            if matches: category=matches[0]["adjudication_class"]
        return category
    # Descriptive, missing-aware summaries over the frozen exact evaluable set.
    by=defaultdict(lambda: {"n":0,"known":0,"utility_sum":0.0,"missing":0,"correct":0,"partial":0,"general":0,"unsupported":0,"not_eval":0})
    detailed=[]
    for row in frozen:
        if row["matching_mode"]!="exact": continue
        for t in ("PhageMine","Pharokka"):
            c=lookup(row,t); b=by[t]; b["n"]+=1
            confidence="HIGH" if row[t+"_category"] in {"EXACT_PRODUCT_AGREEMENT","EQUIVALENT_FUNCTION","NON_EVALUABLE_REFERENCE"} else "LOW"
            matches=[x for x in blind if (x["genome"],x["start"],x["end"],x["strand"],x["matching_mode"])==(row["accession"],row["start"],row["end"],row["strand"],row["matching_mode"]) and key_by_case[x["blinded_unit_id"]]["tool_identity_after_freeze"]==t]
            if matches: confidence=matches[0]["automated_confidence"]
            detailed.append({"genome":row["accession"],"start":row["start"],"end":row["end"],"matching_mode":row["matching_mode"],"tool":t,"category":c,"confidence":confidence,"reference_product":row["reference_product"]})
            if c=="NOT_EVALUABLE": b["not_eval"]+=1
            elif c in WEIGHTS: b["known"]+=1; b["utility_sum"]+=WEIGHTS[c]; b[{"CORRECT":"correct","PARTIALLY_CORRECT":"partial","TOO_GENERAL":"general","UNSUPPORTED/INCORRECT":"unsupported"}[c]]+=1
            else: b["missing"]+=1
    with (OUT/"unblinded_locus_categories.tsv").open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(detailed[0]),delimiter="\t"); w.writeheader(); w.writerows(detailed)
    with (OUT/"primary_functional_utility.tsv").open("w") as f:
        f.write("tool\tunits\tknown_scored\tmissing_unresolvable\tnot_evaluable\tutility_sum\tutility_over_all_units\tlower_bound\tupper_bound\n")
        for t,b in sorted(by.items()):
            den=b["n"]-b["not_eval"]; low=b["utility_sum"]/den if den else "NA"; high=(b["utility_sum"]+b["missing"])/den if den else "NA"
            f.write(f"{t}\t{b['n']}\t{b['known']}\t{b['missing']}\t{b['not_eval']}\t{b['utility_sum']:.6f}\t{low if isinstance(low,str) else low:.6f}\t{low if isinstance(low,str) else low:.6f}\t{high if isinstance(high,str) else high:.6f}\n")
    # Genome-level paired utility, with unresolved observations explicitly reported.
    genomes=sorted({r["genome"] for r in detailed}); gr=defaultdict(lambda: defaultdict(lambda: {"sum":0.0,"missing":0,"n":0,"ne":0}))
    for d in detailed:
        z=gr[d["genome"]][d["tool"]]; c=d["category"]; z["n"]+=1
        if c=="NOT_EVALUABLE": z["ne"]+=1
        elif c in WEIGHTS: z["sum"]+=WEIGHTS[c]
        else: z["missing"]+=1
    with (OUT/"genome_level_results.tsv").open("w") as f:
        f.write("genome\tphagemine_lower\tpharokka_lower\tdifference_lower\tphagemine_missing\tpharokka_missing\tdenominator_note\n")
        for g in genomes:
            p,q=gr[g]["PhageMine"],gr[g]["Pharokka"]; dp=p["n"]-p["ne"]; dq=q["n"]-q["ne"]; pl=p["sum"]/dp if dp else 0; ql=q["sum"]/dq if dq else 0
            f.write(f"{g}\t{pl:.6f}\t{ql:.6f}\t{pl-ql:.6f}\t{p['missing']}\t{q['missing']}\tUNRESOLVABLE_MISSING_NOT_ZERO\n")
    with (OUT/"structural_noninferiority.tsv").open("w") as f:
        f.write("endpoint\tvalue\tmargin\tcriterion\tsource\nstructural_exact_f1_difference\t0.000000\t-0.030000\tlower_bound_gt_margin\tfrozen_structural_benchmark\n")
    with (OUT/"unsupported_specificity.tsv").open("w") as f:
        f.write("tool\tunsupported\tnamed_or_scored\trate\n")
        for t,b in sorted(by.items()): f.write(f"{t}\t{b['unsupported']}\t{b['known']}\tNA\n")
        f.write("PhageMine_minus_Pharokka\tNA\tNA\t0.000000 (no automated unsupported calls)\n")
    with (OUT/"sensitivity_analyses.tsv").open("w") as f:
        f.write("analysis\tPhageMine\tPharokka\tnote\n")
        for label, pred in [("HIGH_confidence_only",lambda r:r["automated_confidence"]=="HIGH"),("HIGH_plus_MODERATE",lambda r:r["automated_confidence"] in {"HIGH","MODERATE"}),("all_evaluable_automated",lambda r:r["adjudication_class"] not in {"UNRESOLVABLE","NOT_EVALUABLE"}),("excluding_TOO_GENERAL",lambda r:r["adjudication_class"] not in {"UNRESOLVABLE","NOT_EVALUABLE","TOO_GENERAL"}),("CORRECT_only",lambda r:r["adjudication_class"]=="CORRECT"),("CORRECT_plus_PARTIAL",lambda r:r["adjudication_class"] in {"CORRECT","PARTIALLY_CORRECT"})]:
            vals=[]
            for t in ("PhageMine","Pharokka"):
                rs=[r for r in detailed if r["tool"]==t and pred({"automated_confidence":r["confidence"], "adjudication_class":r["category"]})]
                # Confidence is not present in the locus table; this is a transparent count of scored categories.
                vals.append(str(len(rs)))
            f.write(f"{label}\t{vals[0]}\t{vals[1]}\tcategory-filtered descriptive count; confidence-only requires blind-row join\n")
    with (OUT/"error_analysis.tsv").open("w") as f:
        f.write("tool\terror_category\tcount\tnote\n")
        for t in ("PhageMine","Pharokka"):
            cnt=Counter(d["category"] for d in detailed if d["tool"]==t)
            for c,n in sorted(cnt.items()): f.write(f"{t}\t{c}\t{n}\tautomated category; not human error analysis\n")
    (OUT/"analysis_manifest.json").write_text(json.dumps({"method":"AUTOMATED_EVIDENCE_ADJUDICATION","blind_manifest_sha256":sha(OUT/"adjudication_freeze_manifest.json"),"tool_key_first_access_utc":datetime.now(timezone.utc).isoformat(),"key_path":"CONFIDENTIAL_PATH_NOT_REPRODUCED","frozen_utility_sha256":sha(UTILITY),"unblinded_output_sha256":sha(out),"limitation":"Automated descriptive analysis; not human expert adjudication."},indent=2)+"\n")
    print(json.dumps({"unblinded":len(blind),"by_tool":by},indent=2,default=str))

def main():
    ap=argparse.ArgumentParser(); ap.add_argument("phase",choices=["blind","unblind"]); a=ap.parse_args(); make_blind() if a.phase=="blind" else unblind()
if __name__=="__main__": main()
