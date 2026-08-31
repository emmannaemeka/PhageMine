"""Read-only sentinel-locus diagnostics for completed PhageMine runs."""
from __future__ import annotations

import csv
import json
from pathlib import Path


def _rows(path: Path):
    if not path.exists():
        return []
    if path.suffix.lower() == ".json":
        data = json.loads(path.read_text())
        return data if isinstance(data, list) else data.get("models", data.get("candidates", []))
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _num(row, *names):
    for name in names:
        if row.get(name) not in (None, ""):
            try:
                return int(row[name])
            except (TypeError, ValueError):
                pass
    return None


def _overlap(a, b):
    if a.get("segment_id") and b.get("segment_id") and a["segment_id"] != b["segment_id"]:
        return 0
    start = max(a["start"], b["start"]); end = min(a["end"], b["end"])
    return max(0, end - start + 1)


def load_sentinels(path: str | Path) -> list[dict]:
    rows = _rows(Path(path)); result = []
    for row in rows:
        start = _num(row, "start", "sentinel_start"); end = _num(row, "end", "sentinel_end")
        if start is None or end is None:
            continue
        result.append({"sentinel_id": row.get("sentinel_id") or row.get("id") or f"sentinel_{len(result)+1}",
                       "segment_id": row.get("segment_id") or row.get("seqid") or "",
                       "start": start, "end": end, "strand": row.get("strand", "")})
    return result


def diagnose(run_dir: str | Path, sentinels: list[dict]) -> list[dict]:
    root = Path(run_dir); raw = root / "gene_calls" / "raw"
    caller_paths = {"phanotate": [raw / "phanotate.tsv", raw / "phanotate.raw.txt"],
                    "pyrodigal": [raw / "pyrodigal.tsv"],
                    "prodigal_gv": [raw / "prodigal_gv.tsv", raw / "prodigal_gv.gff"],
                    "prodigal": [raw / "prodigal.tsv", raw / "prodigal.gff"]}
    calls = {}; availability = {}
    for caller, paths in caller_paths.items():
        rows = []
        for path in paths:
            rows = _rows(path)
            if rows: break
        availability[caller] = bool(rows)
        normalized = []
        for row in rows:
            start = _num(row, "start", "raw_start", "source_start"); end = _num(row, "end", "raw_end", "source_end")
            if start is not None and end is not None:
                normalized.append({"start": start, "end": end, "strand": row.get("strand", row.get("raw_strand", "")), "id": row.get("raw_identifier") or row.get("identifier") or row.get("protein_id", "")})
        calls[caller] = normalized
    # v1.1 commonly preserved only the final GFF, not raw caller streams.
    # Parse it as final-output evidence, never as fabricated raw provenance.
    final_gff = []
    gff_path = root / "genes.gff3"
    if gff_path.exists():
        for line in gff_path.read_text().splitlines():
            if not line or line.startswith("#"): continue
            fields = line.split("\t")
            if len(fields) >= 6:
                try: final_gff.append({"segment_id": fields[0], "start": int(fields[3]), "end": int(fields[4]), "strand": fields[6] if len(fields) > 6 else "", "id": fields[8] if len(fields) > 8 else ""})
                except ValueError: pass
    final_rows = _rows(root / "gene_calls" / "final_gene_models.tsv")
    reconciliation_rows = _rows(root / "gene_calls" / "reconciliation_v2.tsv")
    decisions = _rows(root / "gene_calls" / "model_decisions.tsv")
    output = []
    for sentinel in sentinels:
        caller_results = {}
        any_raw = False
        for caller, models in calls.items():
            overlaps = [m for m in models if _overlap(sentinel, m) > 0]
            any_raw |= bool(overlaps)
            nearest = min(models, key=lambda m: abs(m["start"] - sentinel["start"])) if models else None
            caller_results[caller] = {"status": "AVAILABLE" if availability[caller] else "PROVENANCE_UNAVAILABLE", "raw_overlap": bool(overlaps), "overlapping_calls": overlaps,
                                      "nearest_call": nearest}
        final_records = final_rows + final_gff
        final_overlap = any(_overlap(sentinel, {"start": _num(r, "start"), "end": _num(r, "end"), "segment_id": r.get("segment_id", "")}) > 0 for r in final_records if _num(r, "start") and _num(r, "end"))
        if final_overlap:
            stage = "FINAL_SELECTED"
        elif any_raw and decisions:
            stage = "DROPPED_AT_SELECTION_OR_ADJUDICATION"
        elif any_raw and reconciliation_rows:
            stage = "DROPPED_AT_RECONCILIATION_OR_LATER"
        elif any_raw:
            stage = "RAW_CALL_PRESENT_STAGE_UNKNOWN"
        elif final_gff and not any(availability.values()):
            stage = "STAGE_UNKNOWN_INSUFFICIENT_PROVENANCE"
        else:
            stage = "NEVER_CALLED_BY_AVAILABLE_V1.1_CALLER" if any(availability.values()) else "STAGE_UNKNOWN_INSUFFICIENT_PROVENANCE"
        output.append({"sentinel_id": sentinel["sentinel_id"], "segment_id": sentinel.get("segment_id", ""),
                       "start": sentinel["start"], "end": sentinel["end"], "stage": stage,
                       "callers": caller_results})
    return output


def write_report(results: list[dict], path: str | Path):
    Path(path).write_text(json.dumps(results, indent=2, sort_keys=True))
