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
        "mash": [[resolved, "--version"]],
        "blastn": [[resolved, "-version"]],
        "phanotate.py": [[resolved, "--version"], [resolved, "-h"]],
        "table2asn": [[resolved, "-version"], [resolved, "-h"]],
        "python": [[resolved, "--version"]],
    }
    version = None
    diagnostics: list[str] = []
    for command in probes.get(name, [[resolved, "--version"], [resolved, "-h"]]):
        try:
            result = subprocess.run(command, capture_output=True, text=True, timeout=10, check=False)
            output = "\n".join(part.strip() for part in (result.stdout, result.stderr) if part and part.strip())
            if result.returncode == 0:
                lines = [line.strip() for line in output.splitlines() if line.strip()]
                # HMMER and Prodigal print useful version banners on stderr.
                if lines:
                    version = next((line for line in lines if any(token in line.lower() for token in ("version", "hmmer", "prodigal"))), lines[0])[:240]
                break
            diagnostic = output or f"probe exited with status {result.returncode}"
            diagnostics.append(diagnostic[:500])
        except (OSError, subprocess.SubprocessError) as exc:
            diagnostics.append(str(exc)[:500])
    if not os.access(resolved, os.X_OK):
        status = "INVALID"
    elif version is not None:
        status = "READY"
    else:
        status = "BROKEN"
    return {"name": name, "status": status, "path": resolved, "version": version,
            "diagnostic": diagnostics[-1] if diagnostics and status != "READY" else None}


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
    tools = ["phanotate.py", "prodigal", "hmmscan", "mmseqs", "diamond", "mash", "blastn", "table2asn"]
    executables = [executable_status(tool) for tool in tools]
    resources = manager.validate_all(check_checksum=True)
    by_type = {kind: [r for r in resources if r.get("resource_type") == kind and r.get("status") == "READY"]
               for kind in ("PFAM", "VOGDB", "SWISSPROT", "PHROGS", "PMFDB", "INPHARED_GENOMES")}
    tool_by_name = {item["name"]: item for item in executables}
    table2asn_ready = tool_by_name["table2asn"]["status"] == "READY"
    capabilities = {
        "CORE_ANALYSIS": "READY" if tool_by_name["phanotate.py"]["status"] == "READY" else "UNAVAILABLE",
        "STANDARD_EVIDENCE": "READY" if by_type["PHROGS"] and tool_by_name["mmseqs"]["status"] == "READY" else "UNAVAILABLE",
        "FULL_EVIDENCE": "READY" if all(by_type[k] for k in ("PFAM", "VOGDB", "SWISSPROT", "PHROGS")) and all(tool_by_name[t]["status"] == "READY" for t in ("hmmscan", "mmseqs", "diamond")) else "UNAVAILABLE",
        "PMF_REFERENCE_COMPARISON": "READY" if by_type["PMFDB"] and tool_by_name["mmseqs"]["status"] == "READY" else "UNAVAILABLE",
        "WHOLE_GENOME_REFERENCE_COMPARISON": "READY" if by_type["INPHARED_GENOMES"] and all(tool_by_name[t]["status"] == "READY" for t in ("mash", "blastn")) else "UNAVAILABLE",
        "GENBANK_PRE_SUBMISSION": "READY",
        "NCBI_TABLE2ASN_VALIDATION": "READY" if table2asn_ready else "UNAVAILABLE",
    }
    install_names = {
        "PFAM": "pfam", "VOGDB": "vogdb", "SWISSPROT": "swissprot", "PHROGS": "phrogs",
        "PMFDB": "pmfdb", "INPHARED_GENOMES": "inphared",
    }
    registered_by_type = {kind: [r for r in resources if r.get("resource_type") == kind]
                          for kind in install_names}
    missing_resources = [install_names[kind] for kind in install_names if not registered_by_type[kind]]
    recommendations = []
    if missing_resources:
        recommendations.append({
            "action": "INSTALL_EVIDENCE_DATABASES",
            "missing": missing_resources,
            "commands": [f"phagemine databases install {name}" for name in missing_resources],
            "all_command": "phagemine databases install --all",
            "verify_command": "phagemine doctor",
        })
    repair_resources = []
    for kind, entries in registered_by_type.items():
        if by_type[kind] or not entries:
            continue
        errors = [str(error) for entry in entries for error in entry.get("validation_errors", [])]
        if any(not error.startswith("required executable unavailable:") for error in errors):
            repair_resources.append(install_names[kind])
    if repair_resources:
        recommendations.append({
            "action": "REPAIR_EVIDENCE_DATABASES",
            "resources": repair_resources,
            "commands": [f"phagemine databases install {name} --force" for name in repair_resources],
            "verify_command": "phagemine doctor",
        })
    broken_tools = [item for item in executables
                    if item["name"] != "table2asn" and item["status"] != "READY"]
    if broken_tools:
        conda_packages = {"phanotate.py": "phanotate", "prodigal": "prodigal", "hmmscan": "hmmer",
                          "mmseqs": "mmseqs2", "diamond": "diamond", "mash": "mash", "blastn": "blast"}
        packages = [conda_packages[item["name"]] for item in broken_tools if item["name"] in conda_packages]
        recommendations.append({
            "action": "INSTALL_OR_REPAIR_EXECUTABLES",
            "tools": [{"name": item["name"], "status": item["status"],
                       "diagnostic": item.get("diagnostic")} for item in broken_tools],
            "command": "conda install --channel conda-forge --channel bioconda --strict-channel-priority " + " ".join(packages),
            "verify_command": "phagemine doctor",
        })
    return {"phagemine_version": __version__, "python": executable_status("python"),
            "executables": executables, "resources": resources, "capabilities": capabilities,
            "recommendations": recommendations}


def write_doctor_report(path: str | Path) -> dict[str, Any]:
    payload = doctor()
    Path(path).write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    return payload


def doctor_text(payload: dict[str, Any]) -> str:
    lines = [f"PhageMine {payload['phagemine_version']}", "", "Executables"]
    labels = {"phanotate.py": "PHANOTATE", "prodigal": "Prodigal", "hmmscan": "HMMER", "mmseqs": "MMseqs2", "diamond": "DIAMOND", "mash": "Mash", "blastn": "BLASTN", "table2asn": "table2asn"}
    for item in payload["executables"]:
        status = item["status"] if item["name"] != "table2asn" or item["status"] == "READY" else "OPTIONAL/MISSING"
        lines.append(f"{labels[item['name']]:<14}{status:<17}{item.get('version') or ''}")
        if item.get("diagnostic") and item["name"] != "table2asn":
            lines.append(f"  Error: {item['diagnostic']}")
    lines += ["", "Evidence Resources"]
    labels = {"PFAM": "Pfam", "VOGDB": "VOGDB", "SWISSPROT": "Swiss-Prot", "PHROGS": "PHROGs", "PMFDB": "PMFDB", "INPHARED_GENOMES": "INPHARED genomes"}
    for kind, label in labels.items():
        registered = [r for r in payload["resources"] if r.get("resource_type") == kind]
        ready = [r for r in registered if r.get("status") == "READY"]
        state = "READY" if len(ready) == 1 else ("AMBIGUOUS" if len(ready) > 1 else ("BLOCKED/INVALID" if registered else "NOT INSTALLED"))
        lines.append(f"{label:<18}{state}")
    lines += ["", "Capabilities"]
    for key, value in payload["capabilities"].items():
        lines.append(f"{key.replace('_', ' ').title():<28}{value}")
    recommendations = payload.get("recommendations") or []
    database_recommendation = next((r for r in recommendations if r.get("action") == "INSTALL_EVIDENCE_DATABASES"), None)
    if database_recommendation:
        recommendation = database_recommendation
        lines += ["", "Database setup required"]
        lines.extend(f"  {command}" for command in recommendation["commands"])
        lines += ["", "Or install every evidence database:",
                  f"  {recommendation['all_command']}",
                  "Then verify:", f"  {recommendation['verify_command']}"]
    repair_recommendation = next((r for r in recommendations if r.get("action") == "REPAIR_EVIDENCE_DATABASES"), None)
    if repair_recommendation:
        lines += ["", "Database repair required"]
        lines.extend(f"  {command}" for command in repair_recommendation["commands"])
        lines += ["Then verify:", f"  {repair_recommendation['verify_command']}"]
    tool_recommendation = next((r for r in recommendations if r.get("action") == "INSTALL_OR_REPAIR_EXECUTABLES"), None)
    if tool_recommendation:
        lines += ["", "Executable setup or repair required",
                  f"  {tool_recommendation['command']}",
                  "Then verify:", f"  {tool_recommendation['verify_command']}"]
    return "\n".join(lines)
