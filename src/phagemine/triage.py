"""Human-readable, provenance-only annotator triage report."""
from __future__ import annotations
import csv, json
from pathlib import Path

CONFLICTS = {"STRAND_CONFLICT_UNRESOLVED", "SPLIT_MERGE_UNRESOLVED", "BOUNDARY_CONFLICT_UNRESOLVED", "COMPLEX_CONFLICT", "STRAND_DISCORDANCE", "SPLIT_MODEL", "MERGED_MODEL"}

def _read(path):
    if not Path(path).exists(): return []
    with Path(path).open(newline="") as h: return list(csv.DictReader(h, delimiter="\t"))

def build_triage(root: str | Path) -> dict:
    root = Path(root); calls = root / "gene_calls"
    decisions = _read(calls / "model_decisions.tsv")
    rescue = _read(calls / "rescue_candidates.tsv")
    final = _read(calls / "final_gene_models.tsv")
    review = _read(calls / "gene_model_review.tsv")
    annotation = _read(root / "annotation.tsv")
    conflicts = [row for row in decisions + review if (row.get("decision_class") or row.get("reconciliation_class")) in CONFLICTS]
    changes = [row for row in final if str(row.get("changed_from_legacy", "")).lower() == "true"]
    unresolved = [row for row in annotation if (row.get("functional_state") or row.get("classification")) in {"UNRESOLVED", "UNKNOWN_OR_UNINFORMATIVE", "UNKNOWN"} or row.get("functional_review_flag") == "REVIEW_REQUIRED" or row.get("classification") == "No reliable function identified"]
    # Some legacy/single-RNA runs do not emit final_gene_models.tsv but do
    # emit the authoritative annotation table. Use it only as a count fallback.
    final_count = len(final) or len(annotation)
    return {"final_genes": final_count, "conflicts": conflicts, "rescue_candidates": rescue, "boundary_changes": changes, "unresolved_functions": unresolved,
            "counts": {"final_genes": final_count, "unresolved_conflicts": len(conflicts), "rescue_candidates": len(rescue), "boundary_changes": len(changes), "unresolved_function_proteins": len(unresolved)}}

def write_triage_report(root: str | Path) -> Path:
    root = Path(root); data = build_triage(root); c = data["counts"]
    lines = ["# Annotator's triage report", "", "## Summary", "", f"- Final genes: **{c['final_genes']}**", f"- Unresolved conflicts: **{c['unresolved_conflicts']}**", f"- Rescue candidates (not included by default): **{c['rescue_candidates']}**", f"- Boundary changes from legacy PHANOTATE: **{c['boundary_changes']}**", f"- Unresolved-function proteins: **{c['unresolved_function_proteins']}**", ""]
    lines += ["## Unresolved conflicts", ""]
    if data["conflicts"]:
        for row in data["conflicts"]:
            lines.append(f"- **{row.get('locus_id','unknown locus')}** ({row.get('start','?')}..{row.get('end','?')}): {row.get('decision_reason') or row.get('reason') or row.get('decision_class') or row.get('reconciliation_class')}")
    else: lines.append("No unresolved gene-model conflicts were recorded.")
    lines += ["", "## Rescue candidates", ""]
    if data["rescue_candidates"]:
        for row in data["rescue_candidates"]:
            lines.append(f"- **{row.get('candidate_id','unknown')}** ({row.get('provider','unknown')}, {row.get('start','?')}..{row.get('end','?')}): NOT included in final annotation by default. Evidence/status: {row.get('evidence_status') or row.get('adjudication_status') or 'not recorded'}.")
    else: lines.append("No rescue candidates were recorded.")
    lines += ["", "## Boundary changes", ""]
    if data["boundary_changes"]:
        for row in data["boundary_changes"]: lines.append(f"- {row.get('locus_id','unknown locus')}: selected {row.get('start','?')}..{row.get('end','?')} ({row.get('selection_reason') or row.get('selection_rule','reason not recorded')}).")
    else: lines.append("No boundary changes from the legacy PHANOTATE model were recorded.")
    lines += ["", "## Unresolved functions", ""]
    if data["unresolved_functions"]:
        for row in data["unresolved_functions"]: lines.append(f"- {row.get('protein_id','unknown protein')} ({row.get('coordinates','coordinates unavailable')}): functional state remains unresolved.")
    else: lines.append("No unresolved-function proteins were recorded (or functional classification data were unavailable).")
    out = root / "annotator_triage_report.md"; out.write_text("\n".join(lines) + "\n")
    (root / "annotator_triage_report.json").write_text(json.dumps(data, indent=2, sort_keys=True))
    return out
