#!/usr/bin/env python3
"""Prepare neutral reviewer files from the immutable v1.2 evidence snapshot.

This script only reads frozen files.  It does not run an annotator or alter a
benchmark input.  Reviewer-facing evidence keeps raw values but replaces
tool-specific locus identifiers with neutral case-local identifiers.
"""
from __future__ import annotations

import csv
import hashlib
import json
import random
from datetime import datetime, timezone
from pathlib import Path

from openpyxl import Workbook
from openpyxl.worksheet.datavalidation import DataValidation
from openpyxl.styles import Font, PatternFill, Alignment

ROOT = Path(__file__).resolve().parents[1]
REPO = ROOT.parents[1]
FROZEN = REPO / "evaluation/v1.2_functional"
SEED = 12020260918
JUDGMENTS = ["CORRECT", "PARTIALLY_CORRECT", "TOO_GENERAL", "UNSUPPORTED_SPECIFIC", "INCORRECT", "UNRESOLVABLE", "NOT_EVALUABLE"]
YESNO = ["YES", "PARTIAL", "NO", "UNCLEAR"]
STRENGTH = ["STRONG", "MODERATE", "WEAK", "INSUFFICIENT"]
SPECIFICITY = ["YES", "TOO_SPECIFIC", "TOO_GENERAL", "UNCLEAR"]
CONFLICT = ["YES", "NO", "UNCLEAR"]
ERRORS = ["NONE", "OVER_SPECIFIC", "UNDER_SPECIFIC", "WRONG_FUNCTION", "MISLEADING_DATABASE_HIT", "EVIDENCE_CONFLICT", "TERMINOLOGY_OR_SYNONYM", "HYPOTHETICAL_DESPITE_EVIDENCE", "INSUFFICIENT_EVIDENCE", "OTHER"]


def rows(path):
    with path.open() as f:
        return list(csv.DictReader(f, delimiter="\t"))


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def neutralize(value: str, case_id: str) -> str:
    # Preserve all database/source values while preventing common tool/locus
    # identifiers from revealing which side is PhageMine or Pharokka.
    for prefix in ("PM_", "PHAR_"):
        value = value.replace(prefix, f"{case_id}_LOCUS_")
    return value

def archived_evidence(tool, accession, record, case_id):
    if tool == 'PhageMine':
        obj=json.loads(record['evidence'])
        ann=obj.get('annotation',{}); cls=obj.get('classification',{})
        return {
            'raw': neutralize(record['evidence'],case_id),
            'scores': json.dumps({'best_evidence':ann.get('best_evidence',''),'evidence_tier':cls.get('evidence_tier',''),'supporting_record_count':cls.get('supporting_record_count','')},sort_keys=True),
            'thresholds': json.dumps({'thresholds_recorded':False},sort_keys=True),
            'provenance': json.dumps({'evidence_sources':ann.get('evidence_sources',''),'ortholog_groups':cls.get('ortholog_groups',''),'annotation_review_flag':ann.get('functional_review_flag','')},sort_keys=True),
            'context': json.dumps({'domain_note':ann.get('domain_note',''),'domain_summary':cls.get('domain_summary',''),'functional_category':cls.get('functional_category','')},sort_keys=True),
        }
    attrs=json.loads(record['evidence']); score='NOT_RETAINED_IN_PREDICTIONS_TSV'; source='frozen Pharokka GFF'
    gff=FROZEN/'predictions/Pharokka'/accession/f'{accession}.gff'
    for line in gff.read_text().splitlines():
        fields=line.split('\t')
        if len(fields)>=9 and f"ID={attrs.get('ID','')}" in fields[8]: score=fields[5]; source=fields[1]; break
    return {'raw':neutralize(json.dumps({'gff_attributes':attrs},sort_keys=True),case_id),'scores':json.dumps({'gff_score':score},sort_keys=True),'thresholds':json.dumps({'thresholds_recorded':False},sort_keys=True),'provenance':json.dumps({'source':source,'phrog':attrs.get('phrog','')},sort_keys=True),'context':json.dumps({'function':attrs.get('function','')},sort_keys=True)}


def build_cases():
    loci = rows(FROZEN / "functional_by_locus.tsv")
    pm = {(r["accession"], r["start"], r["end"], r["strand"]): r for r in rows(FROZEN / "phagemine_predictions.tsv")}
    ph = {(r["accession"], r["start"], r["end"], r["strand"]): r for r in rows(FROZEN / "pharokka_predictions.tsv")}
    rng = random.Random(SEED)
    cases, key = [], []
    number = 0
    for row in loci:
        for tool in ("PhageMine", "Pharokka"):
            if row[f"{tool}_category"] != "UNRESOLVED_REVIEW_REQUIRED":
                continue
            number += 1
            case_id = f"FB{number:04d}"
            locus = (row["accession"], row["start"], row["end"], row["strand"])
            other = "Pharokka" if tool == "PhageMine" else "PhageMine"
            order = [tool, other]
            rng.shuffle(order)
            a = pm[locus] if order[0] == "PhageMine" else ph[locus]
            b = pm[locus] if order[1] == "PhageMine" else ph[locus]
            ev_a=archived_evidence(order[0],row["accession"],a,case_id); ev_b=archived_evidence(order[1],row["accession"],b,case_id)
            cases.append({
                "case_id": case_id, "blinded_unit_id": case_id, "matching_mode": row["matching_mode"], "genome": row["accession"],
                "review_target": "A" if order[0] == tool else "B",
                "reference_locus": row["reference_id"], "start": row["start"], "end": row["end"], "strand": row["strand"],
                "reference_product": row["reference_product"], "reference_start": row["reference_start"], "reference_end": row["reference_end"],
                "prediction_A_product": row[f"{order[0]}_product"], "prediction_B_product": row[f"{order[1]}_product"],
                "prediction_A_locus": f"{case_id}_A", "prediction_B_locus": f"{case_id}_B",
                "prediction_A_evidence_raw": ev_a['raw'], "prediction_B_evidence_raw": ev_b['raw'],
                "prediction_A_existing_scores": ev_a['scores'], "prediction_B_existing_scores": ev_b['scores'],
                "prediction_A_existing_thresholds": ev_a['thresholds'], "prediction_B_existing_thresholds": ev_b['thresholds'],
                "prediction_A_annotation_provenance": ev_a['provenance'], "prediction_B_annotation_provenance": ev_b['provenance'],
                "prediction_A_context_evidence": ev_a['context'], "prediction_B_context_evidence": ev_b['context'],
                "prediction_A_evidence_source": "archived functional evidence / GFF attributes",
                "prediction_B_evidence_source": "archived functional evidence / GFF attributes",
                "existing_evidence_note": "Raw value copied from frozen prediction archive; no new search or interpretation performed.",
                "final_judgment": "", "adjudication_class": "", "evidence_basis": "", "uncertainty": "", "review_required": "YES", "reviewer_notes": "",
                "evidence_supports_prediction": "", "evidence_strength": "", "specificity_appropriate": "",
                "conflicting_evidence": "", "preferred_product_name": "", "error_category": "", "reviewer_rationale": "",
                "reviewer_id": "", "review_timestamp": "", "reviewer_confidence": "", "reviewer_comment": "",
            })
            key.append({"case_id": case_id, "prediction_A_tool": order[0], "prediction_B_tool": order[1], "target_tool": tool, "key_classification": "CONFIDENTIAL_NOT_FOR_REVIEWER"})
    return cases, key


def write_csv(path, data):
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(data[0]), lineterminator="\n")
        writer.writeheader(); writer.writerows(data)


def workbook(cases, out):
    wb = Workbook(); instructions = wb.active; instructions.title = "Instructions"
    instructions.append(["Frozen v1.2 blinded adjudication", "Do not rename case IDs or alter evidence fields."])
    instructions.append(["Purpose", "Review unresolved functional assertions only. This workbook does not change the frozen benchmark."])
    instructions.append(["Blinding", "Prediction_A and Prediction_B are neutral labels. Do not seek or infer tool identity."])
    instructions.append(["Submission", "Fill reviewer fields in Adjudication, retain this workbook unchanged, and export the completed sheet as a separate raw submission."])
    instructions.append(["Categories", "See Decision_Definitions. Leave a field blank only when it is not applicable; do not convert missing decisions to INCORRECT."])
    sheet = wb.create_sheet("Adjudication"); headers = list(cases[0]); sheet.append(headers)
    for case in cases: sheet.append([case[h] for h in headers])
    definitions = wb.create_sheet("Decision_Definitions")
    definitions.append(["field", "allowed_value", "definition"])
    definitions_data = {
        "final_judgment": {"CORRECT":"Same biological function or accepted equivalent.", "PARTIALLY_CORRECT":"Meaningful function captured but incomplete or partly mismatched.", "TOO_GENERAL":"Compatible but materially less informative.", "UNSUPPORTED_SPECIFIC":"Specificity exceeds available support.", "INCORRECT":"Conflicts with supported function.", "UNRESOLVABLE":"Evidence cannot distinguish plausible interpretations.", "NOT_EVALUABLE":"Cannot contribute under frozen evaluability rules."},
        "evidence_supports_prediction": {"YES":"Evidence supports the assertion.", "PARTIAL":"Evidence supports only part of the assertion.", "NO":"Evidence does not support it.", "UNCLEAR":"Support cannot be determined."},
        "evidence_strength": {x:x.title() for x in STRENGTH},
        "specificity_appropriate": {"YES":"Appropriate.", "TOO_SPECIFIC":"More specific than evidence supports.", "TOO_GENERAL":"Less specific than supported evidence.", "UNCLEAR":"Cannot determine."},
        "conflicting_evidence": {x:x.title() for x in CONFLICT},
        "error_category": {x:x.replace("_", " ").title() for x in ERRORS},
    }
    for field, values in definitions_data.items():
        for value, definition in values.items(): definitions.append([field, value, definition])
    for ws in wb.worksheets:
        ws.freeze_panes = "A2" if ws.title == "Adjudication" else "A1"
        ws.auto_filter.ref = ws.dimensions
        for cell in ws[1]: cell.font = Font(bold=True, color="FFFFFF"); cell.fill = PatternFill("solid", fgColor="1F4E78")
        for column in ws.columns:
            letter = column[0].column_letter; ws.column_dimensions[letter].width = min(55, max(14, max(len(str(c.value or "")) for c in column[:30]) + 2))
    validations = {"final_judgment": JUDGMENTS, "evidence_supports_prediction": YESNO, "evidence_strength": STRENGTH, "specificity_appropriate": SPECIFICITY, "conflicting_evidence": CONFLICT, "error_category": ERRORS}
    for field, allowed in validations.items():
        col = headers.index(field) + 1; letter = sheet.cell(1, col).column_letter
        dv = DataValidation(type="list", formula1='"' + ",".join(allowed) + '"', allow_blank=True); sheet.add_data_validation(dv); dv.add(f"{letter}2:{letter}{len(cases)+1}")
    wb.save(out)


def main():
    out = ROOT / "reviewer_files"; out.mkdir(exist_ok=True)
    cases, key = build_cases()
    write_csv(out / "blinded_adjudication_cases.csv", cases)
    workbook(cases, out / "blinded_adjudication_cases.xlsx")
    write_csv(ROOT / "blinding" / "tool_blinding_key.tsv", key)
    (ROOT / "blinding" / "tool_blinding_key.tsv").write_text((ROOT / "blinding" / "tool_blinding_key.tsv").read_text().replace(",", "\t"))
    (ROOT / "frozen" / "review_dataset_manifest.json").write_text(json.dumps({"case_count":len(cases),"seed":SEED,"generated_at_utc":datetime.now(timezone.utc).isoformat(),"source_commit":"4933d8fa532a0157c24c3c01dbf1d8adfa787ed7","source_functional_by_locus_sha256":sha(FROZEN/"functional_by_locus.tsv"),"source_prediction_sha256":{"PhageMine":sha(FROZEN/"phagemine_predictions.tsv"),"Pharokka":sha(FROZEN/"pharokka_predictions.tsv")},"blinding":"Prediction_A/Prediction_B randomized independently per unresolved tool-locus unit; key is confidential and excluded from reviewer package."}, indent=2)+"\n")
    print(f"prepared {len(cases)} neutral unresolved judgment units")


if __name__ == "__main__": main()
