"""DeepOrigins: conservative audit of externally computed protein comparison evidence.

This module does not predict structures, infer function, or establish ancestry.
It consumes independently produced comparison results and preserves uncertainty.
"""
from __future__ import annotations

import argparse
import csv
import json
from collections import defaultdict
from pathlib import Path

REQUIRED = ("query", "target", "source_group", "sequence_identity",
            "alignment_coverage", "structural_similarity", "method", "reference")
GROUPS = {"phage", "bacteria", "eukaryote", "other"}


def read_evidence(path: str | Path) -> list[dict]:
    """Read a normalized, user-supplied TSV of comparison observations."""
    with open(path, newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = set(REQUIRED) - set(reader.fieldnames or ())
        if missing:
            raise ValueError("Missing columns: " + ", ".join(sorted(missing)))
        rows = []
        for line, row in enumerate(reader, 2):
            if not row["query"] or not row["target"] or not row["method"] or not row["reference"]:
                raise ValueError(f"Line {line}: missing identifiers, method, or reference")
            if row["source_group"] not in GROUPS:
                raise ValueError(f"Line {line}: invalid source_group")
            for field in ("sequence_identity", "alignment_coverage", "structural_similarity"):
                try:
                    value = float(row[field])
                except (ValueError, TypeError):
                    raise ValueError(f"Line {line}: invalid {field}") from None
                if not 0 <= value <= 1:
                    raise ValueError(f"Line {line}: {field} must be in [0,1]")
                row[field] = value
            rows.append({key: row[key] for key in REQUIRED})
    return rows


def summarize(rows: list[dict]) -> list[dict]:
    """Produce descriptive candidate records, not homology verdicts."""
    grouped = defaultdict(list)
    for row in rows:
        grouped[row["query"]].append(row)
    result = []
    for query, hits in sorted(grouped.items()):
        groups = sorted({h["source_group"] for h in hits})
        methods = sorted({h["method"] for h in hits})
        # This is a review queue only. Structural similarity alone is not ancestry.
        cross_domain = len(set(groups) & {"phage", "bacteria", "eukaryote"}) >= 2
        low_sequence = any(h["sequence_identity"] < 0.25 for h in hits)
        candidate = cross_domain and low_sequence
        result.append({
            "query": query,
            "hit_count": len(hits),
            "source_groups": groups,
            "methods": methods,
            "review_priority": "manual_review" if candidate else "background",
            "interpretation": "relationship_unresolved",
            "caveats": [
                "Structural resemblance does not establish common ancestry.",
                "Assess fold topology, coverage, confidence, sequence profiles, and genomic context.",
                "Check for database redundancy, generic folds, and convergent evolution.",
                "No biological function is inferred from this report.",
            ],
            "evidence": hits,
        })
    return result


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Audit DeepOrigins comparison evidence")
    parser.add_argument("evidence_tsv", help="Normalized comparison evidence TSV")
    parser.add_argument("--output", required=True, help="Output JSON report path")
    args = parser.parse_args(argv)
    records = summarize(read_evidence(args.evidence_tsv))
    destination = Path(args.output)
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(json.dumps({
        "schema_version": "deeporigins-evidence-v0",
        "status": "exploratory_not_validated",
        "candidates": records,
    }, indent=2) + "\n", encoding="utf-8")
    print(f"Wrote {len(records)} candidate summaries to {destination}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
