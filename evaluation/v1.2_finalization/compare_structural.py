#!/usr/bin/env python3
"""Regenerate the PHANOTATE-primary structural comparison table.

This evaluation-only utility consumes completed structural run directories; it
does not invoke callers, modify PhageMine, or alter matching rules.  A run
directory must contain ``reconciliation.json`` and ``selection.json`` from a
frozen evaluation harness.
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path


def main(argv=None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("cohort_root", type=Path, help="cohort directory containing runs/")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)
    rows = []
    for rec_path in sorted((args.cohort_root / "runs").glob("*/reconciliation.json")):
        accession = rec_path.parent.name
        loci = json.loads(rec_path.read_text())
        phanotate = sum(1 for locus in loci for model in locus.get("candidate_models", []) if model.get("caller") == "phanotate")
        selections_path = rec_path.parent / "selection.json"
        selections = json.loads(selections_path.read_text()) if selections_path.exists() else []
        if isinstance(selections, dict):
            selections = selections.get("selections", selections.get("records", []))
        selected = sum(1 for item in selections if item.get("selected_candidate_id"))
        final_path = rec_path.parent / "final_gene_models.tsv"
        final_count = None
        if final_path.exists():
            with final_path.open(newline="") as handle:
                final_count = max(0, sum(1 for _ in handle) - 1)
        rows.append({"accession": accession, "phanotate_primary_cds": phanotate,
                     "phagemine_final_cds": final_count if final_count is not None else "NOT_AVAILABLE",
                     "selection_records": len(selections), "selected_records": selected,
                     "invariant_holds": (str(final_count == phanotate).lower() if final_count is not None else "NOT_AVAILABLE")})
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]) if rows else ["accession"], delimiter="\t")
        writer.writeheader(); writer.writerows(rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
