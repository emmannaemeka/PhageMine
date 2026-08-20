"""Conservative presentation of existing doctor and checkpoint state."""
from __future__ import annotations

import json
import platform
from pathlib import Path


def doctor_rows(payload: dict) -> list[dict]:
    executable_labels = {"phanotate.py": "PHANOTATE", "hmmscan": "HMMER / hmmscan",
                         "mmseqs": "MMseqs2", "diamond": "DIAMOND", "table2asn": "table2asn"}
    rows = []
    for item in payload.get("executables", []):
        if item.get("name") not in executable_labels:
            continue
        raw = item.get("status", "INVALID")
        state = "OPTIONAL" if item.get("name") == "table2asn" and raw != "READY" else ("READY" if raw == "READY" else "NOT READY")
        rows.append({"component": executable_labels[item["name"]], "state": state,
                     "detail": item.get("version") or item.get("path") or "Not found"})
    labels = {"PFAM": "Pfam", "VOGDB": "VOGDB", "SWISSPROT": "Swiss-Prot", "PHROGS": "PHROGs", "PMFDB": "PMFDB"}
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


def capability_rows(payload: dict) -> list[dict]:
    capabilities = payload.get("capabilities", {})
    labels = {"CORE_ANALYSIS": "Core Analysis", "STANDARD_EVIDENCE": "Standard Evidence",
              "FULL_EVIDENCE": "Full Evidence", "GENBANK_PRE_SUBMISSION": "GenBank Pre-submission",
              "NCBI_TABLE2ASN_VALIDATION": "NCBI table2asn validation"}
    rows = []
    for key, label in labels.items():
        raw = capabilities.get(key, "UNAVAILABLE")
        state = "OPTIONAL" if key == "NCBI_TABLE2ASN_VALIDATION" and raw != "READY" else ("READY" if raw == "READY" else "NOT READY")
        rows.append({"capability": label, "state": state})
    return rows


def workflow_readiness(payload: dict, *, mode: str, evidence: str, cohort: bool) -> tuple[bool, str]:
    capabilities = payload.get("capabilities", {})
    if capabilities.get("CORE_ANALYSIS") != "READY":
        return False, "Core analysis needs a validated PHANOTATE installation."
    tools = {item.get("name"): item.get("status") for item in payload.get("executables", [])}
    if mode in {"discover", "both"} and cohort and tools.get("mmseqs") != "READY":
        return False, "Cohort discovery needs MMseqs2 for PMF clustering."
    required = {"standard": "STANDARD_EVIDENCE", "full": "FULL_EVIDENCE"}.get(evidence)
    if required and capabilities.get(required) != "READY":
        return False, f"The {evidence} evidence profile is not ready. Open Environment / Database Status for details."
    return True, "Ready to run with the selected workflow and evidence profile."


def platform_backend_note() -> str:
    if platform.system() == "Windows":
        return ("Windows application status is separate from scientific-backend readiness. "
                "HMMER is not natively supported upstream; full evidence may require WSL2 or a validated local port. "
                "PhageMine never substitutes another algorithm silently.")
    return ("Scientific executables and evidence databases are installed separately from the application "
            "and validated by PhageMine doctor.")


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
