"""Conservative presentation of existing doctor and checkpoint state."""
from __future__ import annotations

import json
from pathlib import Path


def doctor_rows(payload: dict) -> list[dict]:
    executable_labels = {"phanotate.py": "PHANOTATE", "hmmscan": "HMMER / hmmscan",
                         "mmseqs": "MMseqs2", "diamond": "DIAMOND", "mash": "Mash", "table2asn": "table2asn"}
    rows = []
    for item in payload.get("executables", []):
        if item.get("name") not in executable_labels:
            continue
        raw = item.get("status", "INVALID")
        state = "OPTIONAL" if item.get("name") == "table2asn" and raw != "READY" else ("READY" if raw == "READY" else "NOT READY")
        rows.append({"component": executable_labels[item["name"]], "state": state,
                     "detail": item.get("version") or item.get("path") or "Not found"})
    labels = {"PFAM": "Pfam", "VOGDB": "VOGDB", "SWISSPROT": "Swiss-Prot", "PHROGS": "PHROGs", "PMFDB": "PMFDB", "INPHARED_GENOMES": "INPHARED genomes"}
    resources = payload.get("resources", [])
    for kind, label in labels.items():
        matches = [r for r in resources if str(r.get("resource_type", "")).upper() == kind]
        ready = [r for r in matches if r.get("status") == "READY"]
        if ready:
            state, detail = "READY", ready[0].get("path", "Validated")
        elif matches:
            state = "ERROR" if any(r.get("status") == "INVALID" for r in matches) else "NOT READY"
            detail = "; ".join(error for r in matches for error in r.get("validation_errors", [])) or matches[0].get("status", "Not ready")
        else:
            state, detail = "NOT CONFIGURED", "No registered resource"
        rows.append({"component": label, "state": state, "detail": detail})
    return rows


def read_run_status(output: str | Path) -> dict:
    root = Path(output)
    if not root.is_dir():
        return {"state": "NOT FOUND", "output_directory": str(root), "completed_stages": [], "warnings": [], "errors": []}
    statuses = sorted(root.rglob("sample_status.json"))
    samples, warnings, errors, completed, current = [], [], [], [], None
    for path in statuses:
        try:
            item = json.loads(path.read_text()); samples.append(item)
            current = current or item.get("current_stage")
            for checkpoint in item.get("checkpoints", {}).values():
                if checkpoint.get("status") == "COMPLETE": completed.append(checkpoint.get("stage"))
            if item.get("error_message"): errors.append(item["error_message"])
        except (OSError, ValueError, TypeError) as exc:
            warnings.append(f"Could not read {path.name}: {exc}")
    manifest = next((p for p in (root / "batch_manifest.json", root / "run_manifest.json", root / "pooled_manifest.json") if p.is_file()), None)
    state = "COMPLETE" if manifest and not any(s.get("status") in {"RUNNING", "FAILED"} for s in samples) else ("FAILED" if errors else ("RUNNING" if samples else "INCOMPLETE"))
    return {"state": state, "output_directory": str(root), "current_stage": current,
            "completed_stages": sorted(set(filter(None, completed))), "warnings": warnings,
            "errors": errors, "samples": samples, "resume_available": bool(statuses)}
