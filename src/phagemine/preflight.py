"""Fail-closed public preflight and execution-environment diagnostics."""
from __future__ import annotations

import json
import os
import shutil
import subprocess
from pathlib import Path
from typing import Any

from . import __version__
from .resources import EvidenceResourceManager, ResourceStatus, validate_resource


def executable_status(name: str, explicit: str | None = None) -> dict[str, Any]:
    candidate = Path(explicit).expanduser() if explicit else None
    resolved = str(candidate) if candidate and candidate.exists() else (shutil.which(name) if not explicit else None)
    if not resolved:
        return {"name": name, "status": "MISSING", "path": explicit or None, "version": None}
    probes = {
        "prodigal": [[resolved, "-v"], [resolved, "--version"]],
        "hmmscan": [[resolved, "-h"]],
        "mmseqs": [[resolved, "version"]],
        "diamond": [[resolved, "version"]],
        "phanotate.py": [[resolved, "--version"], [resolved, "-h"]],
        "table2asn": [[resolved, "-version"], [resolved, "-h"]],
        "python": [[resolved, "--version"]],
    }
    version = None
    probe_succeeded = False
    for command in probes.get(name, [[resolved, "--version"], [resolved, "-h"]]):
        try:
            result = subprocess.run(command, capture_output=True, text=True, timeout=10, check=False, shell=False)
            text = (result.stdout or result.stderr).strip()
            if result.returncode == 0:
                probe_succeeded = True
            if text and result.returncode == 0:
                lines = [line.strip() for line in text.splitlines() if line.strip()]
                # HMMER and Prodigal print useful version banners on stderr.
                version = next((line for line in lines if any(token in line.lower() for token in ("version", "hmmer", "prodigal"))), lines[0])[:240]
                break
        except (OSError, subprocess.SubprocessError):
            continue
    status = "READY" if os.access(resolved, os.X_OK) else "INVALID"
    # The packaged backend is an explicit path under the application. It must
    # execute successfully; file presence alone is not scientific readiness.
    if explicit and name == "phanotate.py" and not probe_succeeded:
        status = "ERROR"
    return {"name": name, "status": status, "path": resolved, "version": version}


def preflight_resources(explicit: dict[str, dict[str, Any] | None], *, allow_degraded: bool = False) -> dict[str, Any]:
    """Validate explicitly requested evidence before gene prediction."""
    checked: list[dict[str, Any]] = []
    failures: list[str] = []
    for kind, spec in explicit.items():
        if not spec:
            continue
        resource = {"name": kind, "resource_type": kind.upper(), **spec}
        result = validate_resource(resource)
        checked.append(result)
        if result.get("status") != ResourceStatus.READY.value:
            failures.append(f"{kind}: {', '.join(result.get('validation_errors', [])) or result.get('status')}")
    payload = {"status": "READY" if not failures else ("DEGRADED" if allow_degraded else "FAILED"), "resources": checked, "failures": failures}
    if failures and not allow_degraded:
        raise RuntimeError("Preflight failed before gene prediction: " + "; ".join(failures))
    return payload


def preflight_profile(profile: str) -> dict[str, Any]:
    """Resolve a public batch evidence profile against registered resources."""
    profile = profile.lower()
    if profile not in {"core", "standard", "full"}:
        raise ValueError(f"unknown evidence profile: {profile}")
    if profile == "core":
        return {"profile": profile, "resources": [], "status": "READY"}
    manager = EvidenceResourceManager()
    required = ["PHROGS"] if profile == "standard" else ["PFAM", "VOGDB", "SWISSPROT", "PHROGS"]
    resources = []
    failures = []
    for kind in required:
        matches = [r for r in manager.validate_all() if r.get("resource_type") == kind]
        ready = next((r for r in matches if r.get("status") == ResourceStatus.READY.value), None)
        if ready is None:
            failures.append(f"{kind}: no registered READY resource")
        else:
            tool = {"PFAM": "hmmscan", "VOGDB": "hmmscan", "SWISSPROT": "diamond", "PHROGS": "mmseqs"}[kind]
            tool_state = executable_status(tool)
            if tool_state["status"] != "READY":
                failures.append(f"{kind}: required executable {tool} is {tool_state['status']}")
                continue
            resources.append(ready)
    if failures:
        raise RuntimeError("Evidence profile preflight failed: " + "; ".join(failures))
    return {"profile": profile, "resources": resources, "status": "READY"}


def doctor() -> dict[str, Any]:
    manager = EvidenceResourceManager()
    tools = ["phanotate.py", "prodigal", "hmmscan", "mmseqs", "diamond", "table2asn"]
    from .gene_prediction import bundled_phanotate_executable
    bundled = bundled_phanotate_executable()
    executables = [executable_status(tool, bundled if tool == "phanotate.py" else None) for tool in tools]
    resources = manager.validate_all(check_checksum=True)
    by_type = {kind: [r for r in resources if r.get("resource_type") == kind and r.get("status") == "READY"]
               for kind in ("PFAM", "VOGDB", "SWISSPROT", "PHROGS")}
    tool_by_name = {item["name"]: item for item in executables}
    table2asn_ready = tool_by_name["table2asn"]["status"] == "READY"
    capabilities = {
        "CORE_ANALYSIS": "READY" if tool_by_name["phanotate.py"]["status"] == "READY" else "UNAVAILABLE",
        "STANDARD_EVIDENCE": "READY" if by_type["PHROGS"] and tool_by_name["mmseqs"]["status"] == "READY" else "UNAVAILABLE",
        "FULL_EVIDENCE": "READY" if all(by_type[k] for k in by_type) and all(tool_by_name[t]["status"] == "READY" for t in ("hmmscan", "mmseqs", "diamond")) else "UNAVAILABLE",
        "GENBANK_PRE_SUBMISSION": "READY",
        "NCBI_TABLE2ASN_VALIDATION": "READY" if table2asn_ready else "UNAVAILABLE",
    }
    return {"phagemine_version": __version__, "python": executable_status("python"),
            "executables": executables, "resources": resources, "capabilities": capabilities}


def write_doctor_report(path: str | Path) -> dict[str, Any]:
    payload = doctor()
    Path(path).write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    return payload


def doctor_text(payload: dict[str, Any]) -> str:
    lines = [f"PhageMine {payload['phagemine_version']}", "", "Executables"]
    labels = {"phanotate.py": "PHANOTATE", "prodigal": "Prodigal", "hmmscan": "HMMER", "mmseqs": "MMseqs2", "diamond": "DIAMOND", "table2asn": "table2asn"}
    for item in payload["executables"]:
        status = item["status"] if item["name"] != "table2asn" or item["status"] == "READY" else "OPTIONAL/MISSING"
        lines.append(f"{labels[item['name']]:<14}{status:<17}{item.get('version') or ''}")
    lines += ["", "Evidence Resources"]
    labels = {"PFAM": "Pfam", "VOGDB": "VOGDB", "SWISSPROT": "Swiss-Prot", "PHROGS": "PHROGs"}
    for kind, label in labels.items():
        ready = [r for r in payload["resources"] if r.get("resource_type") == kind and r.get("status") == "READY"]
        lines.append(f"{label:<14}{'READY' if len(ready) == 1 else ('AMBIGUOUS' if len(ready) > 1 else 'UNAVAILABLE')}")
    lines += ["", "Capabilities"]
    for key, value in payload["capabilities"].items():
        lines.append(f"{key.replace('_', ' ').title():<28}{value}")
    return "\n".join(lines)
