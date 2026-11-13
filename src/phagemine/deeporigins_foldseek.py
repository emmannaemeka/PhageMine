"""Import Foldseek easy-search TSV results without claiming evolutionary homology."""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

COLUMNS = ("query", "target", "fident", "alnlen", "qlen", "tlen",
           "qcov", "tcov", "evalue", "bits", "qtmscore", "ttmscore")
REQUIRED = set(COLUMNS)


def parse_foldseek(path: str | Path, *, target_groups: dict[str, str] | None = None) -> list[dict]:
    """Parse a headerless TSV produced with --format-output matching COLUMNS.

    Foldseek fractions are not normalized to the evidence-audit schema:
    original values and provenance remain intact.
    """
    results = []
    with open(path, newline="", encoding="utf-8") as handle:
        reader = csv.reader(handle, delimiter="\t")
        for line, fields in enumerate(reader, 1):
            if not fields:
                continue
            if len(fields) != len(COLUMNS):
                raise ValueError(f"Line {line}: expected {len(COLUMNS)} Foldseek fields, got {len(fields)}")
            row = dict(zip(COLUMNS, fields))
            if not row["query"] or not row["target"]:
                raise ValueError(f"Line {line}: empty protein identifier")
            for key in ("fident", "qcov", "tcov", "qtmscore", "ttmscore"):
                try:
                    val = float(row[key])
                except ValueError:
                    raise ValueError(f"Line {line}: invalid {key}") from None
                if not math.isfinite(val) or not 0 <= val <= 1:
                    raise ValueError(f"Line {line}: {key} outside [0,1]")
                row[key] = val
            for key in ("alnlen", "qlen", "tlen"):
                try:
                    val = int(row[key])
                except ValueError:
                    raise ValueError(f"Line {line}: invalid {key}") from None
                if val <= 0:
                    raise ValueError(f"Line {line}: {key} must be positive")
                row[key] = val
            for key in ("evalue", "bits"):
                try:
                    val = float(row[key])
                except ValueError:
                    raise ValueError(f"Line {line}: invalid {key}") from None
                if not math.isfinite(val) or val < 0:
                    raise ValueError(f"Line {line}: invalid {key}")
                row[key] = val
            if row["alnlen"] > max(row["qlen"], row["tlen"]) * 2:
                raise ValueError(f"Line {line}: implausible alignment length")
            row["target_group"] = (target_groups or {}).get(row["target"], "unknown")
            row["interpretation"] = "structural_match_only"
            results.append(row)
    return results


def read_groups(path: str | Path) -> dict[str, str]:
    with open(path, newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not {"target", "group"}.issubset(reader.fieldnames or []):
            raise ValueError("Target mapping must have target and group columns")
        groups = {}
        for row in reader:
            if not row["target"] or row["group"] not in {"phage", "bacteria", "eukaryote", "other"}:
                raise ValueError("Invalid target mapping")
            if row["target"] in groups and groups[row["target"]] != row["group"]:
                raise ValueError("Conflicting target groups")
            groups[row["target"]] = row["group"]
    return groups


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Import Foldseek structural matches into DeepOrigins")
    parser.add_argument("foldseek_tsv")
    parser.add_argument("--output", required=True)
    parser.add_argument("--target-groups", help="Optional TSV mapping target to group")
    parser.add_argument("--database", required=True, help="Exact target database name/version")
    parser.add_argument("--foldseek-version", required=True)
    args = parser.parse_args(argv)
    groups = read_groups(args.target_groups) if args.target_groups else {}
    rows = parse_foldseek(args.foldseek_tsv, target_groups=groups)
    report = {
        "schema_version": "deeporigins-foldseek-v0",
        "status": "unvalidated_structural_matches",
        "provenance": {"source_file": str(args.foldseek_tsv),
                       "database": args.database, "foldseek_version": args.foldseek_version,
                       "format_output": ",".join(COLUMNS)},
        "hits": rows,
        "warning": "Foldseek matches do not prove ancestry, function, or novelty.",
    }
    destination = Path(args.output)
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(f"Imported {len(rows)} structural matches")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
