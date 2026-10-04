"""Conservative presentation of existing doctor and checkpoint state."""
from __future__ import annotations

import json
from pathlib import Path


def doctor_rows(payload: dict) -> list[dict]:
    executable_labels = {"phanotate.py": "PHANOTATE", "hmmscan": "HMMER / hmmscan",
                         "mmseqs": "MMseqs2", "diamond": "DIAMOND", "mash": "Mash",
                         "blastn": "BLASTN", "prodigal": "Prodigal", "table2asn": "table2asn"}
    rows = []
    for item in payload.get("executables", []):
        if item.get("name") not in executable_labels:
            continue
        raw = item.get("status", "INVALID")
        state = "OPTIONAL" if item.get("name") == "table2asn" and raw != "READY" else ("READY" if raw == "READY" else "NOT READY")
        rows.append({"component": executable_labels[item["name"]], "state": state,
                     "detail": item.get("diagnostic") or item.get("version") or item.get("path") or "Not found"})
    for item in payload.get("python_modules", []):
        rows.append({"component": item["name"], "state": item.get("status", "NOT READY"),
                     "detail": item.get("diagnostic") or item.get("version") or "Not found"})
    for name, state in payload.get("capabilities", {}).items():
        rows.append({"component": name.replace("_", " ").title(), "state": state,
                     "detail": "Readiness for the installed workflow"})
    for item in payload.get("deep_checks", []):
        rows.append({"component": item["name"], "state": item["status"],
                     "detail": item.get("diagnostic") or "Operational format check"})
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
    root = Path(output).expanduser()
    if not root.is_dir():
        return {"state": "NOT FOUND", "output_directory": str(root), "completed_stages": [], "warnings": [], "errors": []}
    statuses = sorted(root.rglob("sample_status.json"))
    samples, warnings, errors, completed, current = [], [], [], [], None
    for path in statuses:
        try:
            item = json.loads(path.read_text())
            if not isinstance(item, dict) or not isinstance(item.get("checkpoints", {}), dict):
                raise ValueError("Expected a sample status object with checkpoints")
            samples.append(item)
            current = current or item.get("current_stage")
            for checkpoint in item.get("checkpoints", {}).values():
                if not isinstance(checkpoint, dict):
                    raise ValueError("Expected a checkpoint object")
                if checkpoint.get("status") == "COMPLETE": completed.append(checkpoint.get("stage"))
            if item.get("error_message"): errors.append(item["error_message"])
        except (OSError, ValueError, TypeError) as exc:
            warnings.append(f"Could not read {path.relative_to(root)}: {exc}")

    def manifest_state(folder: Path) -> str:
        for name in ("batch_manifest.json", "pooled_manifest.json", "run_manifest.json"):
            path = folder / name
            if not path.is_file():
                continue
            try:
                payload = json.loads(path.read_text())
                if not isinstance(payload, dict):
                    raise ValueError("Expected a manifest object")
                if name == "batch_manifest.json":
                    if payload.get("failures"):
                        errors.append(f"{path.relative_to(root)} records failed samples")
                        return "FAILED"
                    return "COMPLETE" if payload.get("ended_at") else "RUNNING"
                if name == "pooled_manifest.json":
                    return "COMPLETE" if "pmf_count" in payload and (folder / "discovery_report.html").is_file() else "INCOMPLETE"
                return "COMPLETE" if "stage_timings_seconds" in payload and (folder / "report.html").is_file() else "INCOMPLETE"
            except (OSError, ValueError, TypeError) as exc:
                warnings.append(f"Could not read {path.relative_to(root)}: {exc}")
                return "INCOMPLETE"
        return "INCOMPLETE"

    if (root / "annotation" / "batch_manifest.json").is_file():
        states = [manifest_state(root / "annotation"), manifest_state(root / "discovery")]
        state = "COMPLETE" if all(s == "COMPLETE" for s in states) else ("FAILED" if "FAILED" in states else "INCOMPLETE")
    else:
        state = manifest_state(root)
    if errors or any(s.get("status") == "FAILED" for s in samples):
        state = "FAILED"
    elif warnings:
        state = "INCOMPLETE"
    elif any(s.get("status") == "RUNNING" for s in samples):
        state = "RUNNING"
    return {"state": state, "output_directory": str(root), "current_stage": current,
            "completed_stages": sorted(set(filter(None, completed))), "warnings": warnings,
            "errors": errors, "samples": samples, "resume_available": bool(statuses)}
