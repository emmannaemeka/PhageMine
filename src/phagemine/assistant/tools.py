"""Structured tools the assistant may call.

Two kinds, kept apart:

* read tools answer from existing results and never write analysis outputs;
* job tools start a tracked worker process and return a job ID, never results.

There is no shell tool and no code-execution tool. Every handler receives a
``ToolContext`` bound to one project by the server; the model cannot choose
the project or reach paths outside it.
"""
from __future__ import annotations

import csv
import json
import os
import re
import shutil
import uuid
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable

from . import approvals, jobs
from .ncbi import CURATION_FIELDS, REVISABLE_FIELDS, check_metadata, check_product, classify
from .store import FASTA_SUFFIXES, Project, now, safe_name

ROW_LIMIT = 50
FUNCTIONAL_STATES = ["KNOWN_FUNCTION", "PROBABLE_FUNCTION", "FUNCTIONAL_CLASS_ONLY",
                     "CONSERVED_UNKNOWN", "CONFLICTING_EVIDENCE", "UNRESOLVED"]
UNKNOWN_STATES = {"CONSERVED_UNKNOWN", "UNRESOLVED"}


def demo_mode() -> bool:
    """Demonstration runs use PhageMine's fixture gene caller and mock evidence."""
    return os.environ.get("PHAGEMINE_ASSISTANT_DEMO") == "1"


@dataclass
class ToolContext:
    project: Project
    user_text: str = ""          # everything the researcher has typed in this project
    turn: int = 0                 # index of the current user turn
    grounding: Callable[[str], list[str]] | None = None


@dataclass
class Tool:
    name: str
    kind: str                     # "read" or "job"
    description: str
    parameters: dict
    handler: Callable[..., Any]
    required: list[str] = field(default_factory=list)

    def schema(self) -> dict:
        return {"name": self.name, "description": self.description,
                "input_schema": {"type": "object", "properties": self.parameters,
                                 "required": self.required, "additionalProperties": False}}


REGISTRY: dict[str, Tool] = {}


def tool(name: str, kind: str, description: str, parameters: dict | None = None, required: list[str] | None = None):
    def register(func):
        REGISTRY[name] = Tool(name, kind, description, parameters or {}, func, required or [])
        return func
    return register


def schemas() -> list[dict]:
    return [t.schema() for t in REGISTRY.values()]


def execute(ctx: ToolContext, name: str, arguments: dict | None) -> dict:
    """Run one tool. Errors come back as data so the model can explain them."""
    if name not in REGISTRY:
        return {"ok": False, "error": f"Unknown tool {name}"}
    spec = REGISTRY[name]
    arguments = dict(arguments or {})
    unexpected = set(arguments) - set(spec.parameters)
    missing = [key for key in spec.required if arguments.get(key) in (None, "", [])]
    problems = ([f"unexpected {', '.join(sorted(unexpected))}"] if unexpected else []) + \
               ([f"missing {', '.join(missing)}"] if missing else [])
    for key, value in arguments.items():
        if key in spec.parameters:
            problems += _check_value(key, value, spec.parameters[key])
    if problems:
        return {"ok": False, "error": "Invalid arguments: " + "; ".join(problems)}
    try:
        result = spec.handler(ctx, **arguments)
    except (ValueError, OSError, KeyError, json.JSONDecodeError) as exc:
        return {"ok": False, "error": str(exc)}
    return {"ok": True, **result}


_TYPES = {"string": str, "integer": int, "boolean": bool, "array": list, "object": dict}


def _check_value(name: str, value, schema: dict) -> list[str]:
    """Runtime check of the JSON-schema subset the tools use: type, enum, bounds, items."""
    expected = schema.get("type")
    if value is None:
        return [f"{name} must not be null"]
    kind = _TYPES.get(expected)
    if kind and (not isinstance(value, kind) or (kind is int and isinstance(value, bool))):
        return [f"{name} must be {expected}"]
    problems = []
    if "enum" in schema and value not in schema["enum"]:
        problems.append(f"{name} must be one of {', '.join(map(str, schema['enum']))}")
    if expected == "integer":
        if "minimum" in schema and value < schema["minimum"]:
            problems.append(f"{name} must be at least {schema['minimum']}")
        if "maximum" in schema and value > schema["maximum"]:
            problems.append(f"{name} must be at most {schema['maximum']}")
    if expected == "array":
        if len(value) < schema.get("minItems", 0):
            problems.append(f"{name} needs at least {schema['minItems']} item(s)")
        for index, item in enumerate(value):
            problems += _check_value(f"{name}[{index}]", item, schema.get("items", {}))
    if expected == "string" and len(value) > schema.get("maxLength", 20000):
        problems.append(f"{name} is too long")
    return problems


# ---------------------------------------------------------------- helpers ----

def _tsv(path: Path) -> list[dict]:
    if not path.is_file():
        return []
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _json(path: Path, default):
    try:
        return json.loads(path.read_text())
    except (OSError, ValueError):
        return default


def _records(ctx: ToolContext, sample: str, run: str | Path | None = None) -> list[dict]:
    """Annotation rows for a sample: from ``run`` when given (for example the exact run a
    comparison used), otherwise from the sample's latest completed annotation."""
    from phagemine.gui.services.results import annotation_records
    run = ctx.project.contained(run) if run else ctx.project.run_dir(sample)
    rows = annotation_records(run)
    for row in rows:
        row["sample"] = sample
    return rows


def _brief(row: dict) -> dict:
    keys = ("sample", "protein_id", "start", "end", "strand", "length_aa", "product",
            "functional_state", "confidence", "evidence_sources", "review_flag")
    return {k: row.get(k) for k in keys if row.get(k) not in (None, "")}


def _evidence_mode(run: Path) -> str:
    manifest = _json(run / "run_manifest.json", {})
    statuses = {str(a.get("status")) for a in manifest.get("evidence_adapters", [])}
    caller = (manifest.get("gene_caller") or {}).get("name", "")
    if "MOCK" in statuses or "demo" in caller.lower():
        return "DEMO (mock evidence or fixture gene caller; not biological evidence)"
    return "REAL"


def _job_brief(record: dict) -> dict:
    keys = ("id", "kind", "status", "sample", "samples", "created", "finished", "error", "depends_on")
    return {k: record.get(k) for k in keys if record.get(k) not in (None, [], "")}


# ------------------------------------------------------------- read tools ----

@tool("list_project", "read",
      "List uploaded files, analysis jobs with their status, completed samples, comparisons, figures and saved summaries in this project. Call this first when unsure what exists.")
def list_project(ctx: ToolContext) -> dict:
    project = ctx.project
    return {"project": project.name, "demo_mode": demo_mode(), "inputs": project.inputs(),
            "jobs": [_job_brief(jobs.status(project, j["id"])) for j in project.jobs()],
            "completed_samples": sorted(project.annotation_runs()),
            "figures": [p.name for p in sorted((project.root / "figures").glob("*.png"))],
            "summaries": [p.name for p in sorted((project.root / "summaries").glob("*.md"))]}


@tool("inspect_fasta", "read",
      "Check an uploaded FASTA file: record names, lengths, GC content, ambiguous bases and warnings. Use before annotating.",
      {"file": {"type": "string", "description": "Uploaded file name"}}, ["file"])
def inspect_fasta(ctx: ToolContext, file: str) -> dict:
    path = ctx.project.input_path(file)
    from phagemine.io import read_fasta  # noqa: F401  (validates that PhageMine can parse)
    records, current, seq = [], None, []
    for line in path.read_text(errors="replace").splitlines():
        if line.startswith(">"):
            if current is not None:
                records.append((current, "".join(seq)))
            current, seq = line[1:].split()[0] if line[1:].split() else "(unnamed)", []
        elif line.strip():
            seq.append(line.strip().upper())
    if current is not None:
        records.append((current, "".join(seq)))
    if not records:
        raise ValueError(f"{file} contains no FASTA records")
    summary, warnings = [], []
    for name, sequence in records:
        length = len(sequence)
        gc = (sequence.count("G") + sequence.count("C")) / length * 100 if length else 0
        ambiguous = sum(1 for base in sequence if base not in "ACGT")
        summary.append({"record": name, "length_bp": length, "gc_percent": round(gc, 2), "non_ACGT_bases": ambiguous})
        if length < 1000:
            warnings.append(f"{name} is only {length} bp; most phage genomes are 3-500 kb")
        if ambiguous:
            warnings.append(f"{name} has {ambiguous} ambiguous bases")
    if len(records) > 1:
        warnings.append("More than one record: the single-genome workflow expects one genome per file (segmented RNA genomes are a separate workflow)")
    return {"file": file, "records": summary, "warnings": warnings}


@tool("check_environment", "read",
      "Report which external tools (PHANOTATE, MMseqs2, HMMER, DIAMOND, table2asn) and evidence databases are installed and ready.")
def check_environment(ctx: ToolContext) -> dict:
    from phagemine.preflight import doctor
    payload = doctor()
    tools = {item["name"]: item.get("status") for item in payload.get("executables", [])}
    databases = {}
    for item in payload.get("resources", []):
        databases[str(item.get("resource_type"))] = {"status": item.get("status"), "version": item.get("version")}
    return {"phagemine_version": payload.get("phagemine_version"), "executables": tools,
            "databases": databases or "none registered", "capabilities": payload.get("capabilities", {}),
            "demo_mode": demo_mode(),
            "notes": ["Missing databases reduce evidence; they do not stop annotation.",
                      "Without MMseqs2, comparisons find homologs only through shared PHROG assignments and identical sequences."]}


@tool("get_job_status", "read", "Status, progress stage, error and last log lines of a job.",
      {"job_id": {"type": "string"}}, ["job_id"])
def get_job_status(ctx: ToolContext, job_id: str) -> dict:
    record = jobs.status(ctx.project, job_id)
    keys = ("id", "kind", "status", "sample", "samples", "progress", "error", "exit_codes",
            "depends_on", "created", "started", "finished", "log_tail")
    return {k: record.get(k) for k in keys if record.get(k) not in (None, [])}


@tool("get_run_summary", "read",
      "Summary of one completed annotation: protein counts by functional state, evidence databases used and their versions, gene caller, review priorities, submission readiness.",
      {"sample": {"type": "string", "description": "Sample name from completed_samples"}}, ["sample"])
def get_run_summary(ctx: ToolContext, sample: str) -> dict:
    run = ctx.project.run_dir(sample)
    manifest = _json(run / "run_manifest.json", {})
    classes = _json(run / "functional_classification.json", [])
    states = Counter(row.get("functional_state") or "UNRESOLVED" for row in classes)
    adapters = [{"database": a.get("adapter"), "status": a.get("status"),
                 "version": (a.get("provenance") or {}).get("pfam_version") or (a.get("provenance") or {}).get("database_version")
                 or (a.get("provenance") or {}).get("adapter_version")}
                for a in manifest.get("evidence_adapters", [])]
    triage = _json(run / "annotator_triage_report.json", {})
    validation = _json(run / "genbank_submission" / "validation.json", {})
    return {"sample": sample, "evidence_mode": _evidence_mode(run), "proteins": len(classes),
            "functional_states": dict(states), "evidence_databases": adapters,
            "gene_caller": manifest.get("gene_caller"), "pipeline_version": manifest.get("pipeline_version"),
            "input_sha256": manifest.get("input_sha256"),
            "hallmarks": manifest.get("hallmark_summary"), "annotation_review": manifest.get("annotation_review"),
            "triage_summary": triage.get("summary", triage if len(json.dumps(triage)) < 1500 else None),
            "submission_readiness": validation.get("submission_readiness"),
            "caution": "UNRESOLVED or CONSERVED_UNKNOWN means no reliable function was assigned; it does not mean novel."}


@tool("query_proteins", "read",
      "Find proteins in completed annotations. Filter by sample, functional state, text in product/function, or review flag. Returns up to 50 rows and the total count.",
      {"sample": {"type": "string", "description": "Omit to search all completed samples"},
       "functional_states": {"type": "array", "items": {"type": "string", "enum": FUNCTIONAL_STATES}},
       "unknown_function_only": {"type": "boolean", "description": "Shortcut for CONSERVED_UNKNOWN and UNRESOLVED"},
       "text": {"type": "string", "description": "Words that must all appear in product, proposed function or protein ID"},
       "needs_review": {"type": "boolean"},
       "limit": {"type": "integer", "minimum": 1, "maximum": ROW_LIMIT}})
def query_proteins(ctx: ToolContext, sample: str | None = None, functional_states: list[str] | None = None,
                   unknown_function_only: bool = False, text: str | None = None,
                   needs_review: bool = False, limit: int = ROW_LIMIT) -> dict:
    samples = [sample] if sample else sorted(ctx.project.annotation_runs())
    if not samples:
        raise ValueError("No completed annotations in this project yet")
    wanted = set(functional_states or []) | (UNKNOWN_STATES if unknown_function_only else set())
    terms = re.findall(r"[\w.-]+", (text or "").lower())
    matches = []
    for name in samples:
        for row in _records(ctx, name):
            if wanted and row.get("functional_state") not in wanted:
                continue
            haystack = " ".join(str(row.get(k) or "") for k in ("protein_id", "product", "proposed_function", "functional_category")).lower()
            if terms and not all(t in haystack for t in terms):
                continue
            if needs_review and str(row.get("review_flag") or "NONE") in {"NONE", ""}:
                continue
            matches.append(_brief(row))
    limit = max(1, min(int(limit or ROW_LIMIT), ROW_LIMIT))
    return {"total_matches": len(matches), "returned": min(len(matches), limit), "proteins": matches[:limit],
            "filters": {"samples": samples, "functional_states": sorted(wanted), "text": text, "needs_review": needs_review}}


@tool("get_protein_evidence", "read",
      "All recorded evidence for one protein: database hits with versions and scores, functional state, conflicts, genomic neighbours, and homologs in compared genomes.",
      {"sample": {"type": "string"}, "protein_id": {"type": "string"},
       "include_sequence": {"type": "boolean", "description": "Only when the researcher asks for the sequence"}},
      ["sample", "protein_id"])
def get_protein_evidence(ctx: ToolContext, sample: str, protein_id: str, include_sequence: bool = False) -> dict:
    row = next((r for r in _records(ctx, sample) if r.get("protein_id") == protein_id), None)
    if row is None:
        raise ValueError(f"{protein_id} is not in {sample}")
    evidence = [{k: item.get(k) for k in ("source", "source_version", "status", "supports", "level", "statement",
                                          "identifier", "description", "family_name", "metrics")
                 if item.get(k) not in (None, "", {})} for item in row.get("evidence", [])]
    context = row.get("genomic_context") or {}
    result = {"sample": sample, "evidence_mode": _evidence_mode(ctx.project.run_dir(sample)),
              "annotation": _brief(row), "proposed_function": row.get("proposed_function"),
              "classification": row.get("classification"), "evidence": evidence,
              "conflicts": row.get("conflict_review"),
              "genomic_context": {k: v for k, v in context.items() if k not in {"sequence"}} if context else None,
              "homologs_in_compared_genomes": _homologs_for(ctx, sample, protein_id)}
    if include_sequence:
        result["protein_sequence"] = row.get("sequence")
    return result


def _comparison(ctx: ToolContext, comparison_job_id: str | None) -> dict:
    done = [j for j in ctx.project.jobs() if j.get("kind") == "comparison" and j.get("status") == "COMPLETE"]
    if comparison_job_id:
        record = ctx.project.job(comparison_job_id)
        if record.get("kind") != "comparison":
            raise ValueError(f"{comparison_job_id} is not a comparison job")
        if record.get("status") != "COMPLETE":
            raise ValueError(f"Comparison {comparison_job_id} is {record.get('status')}, not COMPLETE")
        return record
    if not done:
        raise ValueError("No completed comparison in this project yet")
    return done[-1]


def _links(ctx: ToolContext, record: dict) -> list[dict]:
    """Ortholog links from every direction, with run paths mapped back to sample names."""
    by_path = {str(ctx.project.contained(path)): name for name, path in record["run_dirs"].items()}
    links = []
    for direction in record["directions"]:
        for row in _tsv(ctx.project.contained(direction["output"]) / "ortholog_links.tsv"):
            target = by_path.get(str(Path(row["reference_genome_id"]).resolve()))
            if target:
                links.append({"sample_a": direction["query"], "protein_a": row["query_protein_id"],
                              "sample_b": target, "protein_b": row["reference_protein_id"],
                              "method": row["orthology_method"], "identity": row["sequence_identity"],
                              "query_coverage": row["query_coverage"], "target_coverage": row["target_coverage"]})
    return links


def _homologs_for(ctx: ToolContext, sample: str, protein_id: str) -> list[dict] | dict | None:
    try:
        record = _comparison(ctx, None)
    except ValueError:
        return None
    if sample not in record.get("run_dirs", {}):
        return None
    if ctx.project.contained(record["run_dirs"][sample]) != ctx.project.run_dir(sample):
        return {"note": f"The latest comparison used an earlier annotation of {sample}; re-run the comparison to see homologs for the current one."}
    return [p for p in _canonical_pairs(ctx, record, *_pair_of(record, sample)) if p["protein_a"] == protein_id][:20]


def _pair_of(record: dict, sample: str) -> tuple[str, str]:
    other = next(s for s in record["samples"] if s != sample)
    return sample, other


def _canonical_pairs(ctx: ToolContext, record: dict, a: str, b: str) -> list[dict]:
    """Deduplicated A-B protein pairs from both directions, with which direction(s) found each."""
    pairs: dict[tuple[str, str], dict] = {}
    for link in _links(ctx, record):
        if (link["sample_a"], link["sample_b"]) == (a, b):
            key, direction = (link["protein_a"], link["protein_b"]), f"{a}->{b}"
        elif (link["sample_a"], link["sample_b"]) == (b, a):
            key, direction = (link["protein_b"], link["protein_a"]), f"{b}->{a}"
        else:
            continue
        entry = pairs.setdefault(key, {"protein_a": key[0], "protein_b": key[1], "directions": set(), "methods": set(),
                                       "identity": link["identity"]})
        entry["directions"].add(direction)
        entry["methods"].add(link["method"])
    result = []
    for entry in pairs.values():
        directions = entry.pop("directions")
        entry["support"] = "both" if len(directions) == 2 else next(iter(directions))
        entry["methods"] = sorted(entry["methods"])
        result.append(entry)
    return sorted(result, key=lambda e: (e["protein_a"], e["protein_b"]))


def _method_status(ctx: ToolContext, record: dict) -> dict:
    executed = bool(record.get("mmseqs2_executed"))
    status = {"mmseqs2_used": executed, "mmseqs_searches": record.get("mmseqs_searches")}
    if not executed:
        status["warning"] = ("MMseqs2 searches did not run successfully for this comparison, so homologs come only from "
                             "shared PHROG assignments and identical sequences. Shared counts are a lower bound.")
    return status


@tool("get_shared_genes", "read",
      "Which proteins two compared genomes share: deduplicated protein pairs found in either direction (with which direction found each), counts of matched proteins in each genome, whether MMseqs2 actually ran, and how many shared proteins have unknown function. Uses the exact annotation runs the comparison was made from.",
      {"comparison_job_id": {"type": "string", "description": "Omit to use the latest completed comparison"},
       "sample_a": {"type": "string"}, "sample_b": {"type": "string"},
       "unknown_function_only": {"type": "boolean"},
       "limit": {"type": "integer", "minimum": 1, "maximum": ROW_LIMIT}})
def get_shared_genes(ctx: ToolContext, comparison_job_id: str | None = None, sample_a: str | None = None,
                     sample_b: str | None = None, unknown_function_only: bool = False, limit: int = ROW_LIMIT) -> dict:
    record = _comparison(ctx, comparison_job_id)
    samples = record["samples"]
    a = sample_a or samples[0]
    b = sample_b or next(s for s in samples if s != a)
    for name in (a, b):
        if name not in samples:
            raise ValueError(f"{name} was not part of comparison {record['id']} ({', '.join(samples)})")
    runs = {name: record["run_dirs"][name] for name in (a, b)}
    annotations = {name: {r["protein_id"]: r for r in _records(ctx, name, runs[name])} for name in (a, b)}
    pairs = _canonical_pairs(ctx, record, a, b)
    matched = {a: {p["protein_a"] for p in pairs}, b: {p["protein_b"] for p in pairs}}
    rows = []
    for pair in pairs:
        pa, pb = annotations[a].get(pair["protein_a"], {}), annotations[b].get(pair["protein_b"], {})
        rows.append({f"{a}_protein": pair["protein_a"], f"{a}_product": pa.get("product"),
                     f"{a}_functional_state": pa.get("functional_state"),
                     f"{b}_protein": pair["protein_b"], f"{b}_product": pb.get("product"),
                     "support": pair["support"], "methods": pair["methods"], "identity": pair["identity"]})
    if unknown_function_only:
        rows = [r for r in rows if r[f"{a}_functional_state"] in UNKNOWN_STATES]
    unknown = sum(1 for pid in matched[a] if annotations[a].get(pid, {}).get("functional_state") in UNKNOWN_STATES)
    limit = max(1, min(int(limit or ROW_LIMIT), ROW_LIMIT))
    stale = [n for n in (a, b) if ctx.project.annotation_runs().get(n, {}).get("output") not in (None, runs[n])]
    return {"comparison_job_id": record["id"], "sample_a": a, "sample_b": b, "source_runs": runs,
            "proteins_total": {a: len(annotations[a]), b: len(annotations[b])},
            "proteins_with_homolog_in_other": {a: len(matched[a]), b: len(matched[b])},
            "distinct_pairs": len(pairs), "shared_with_unknown_function_in_" + a: unknown,
            **_method_status(ctx, record),
            "thresholds": "MMseqs2: E-value <= 1e-3, query and target coverage >= 0.5 (PhageMine compare defaults)",
            "pairs_returned": min(len(rows), limit), "pairs_total": len(rows), "pairs": rows[:limit],
            **({"stale_note": f"{', '.join(stale)} has been re-annotated since this comparison; results describe the earlier run."} if stale else {})}


@tool("get_submission_status", "read",
      "NCBI pre-submission package status for a sample: readiness, missing metadata, PhageMine validation errors, table2asn result, and corrections applied so far.",
      {"sample": {"type": "string"}}, ["sample"])
def get_submission_status(ctx: ToolContext, sample: str) -> dict:
    run = ctx.project.run_dir(sample)
    revision = _latest_revision(ctx, sample)
    package = (ctx.project.contained(revision["output"]) / "submission_v2" / "genbank_submission") if revision else run / "genbank_submission"
    validation = _json(package / "validation.json", None)
    if validation is None:
        raise ValueError(f"No submission package found for {sample}")
    pre = validation.get("phagemine_pre_submission_validation", {})
    t2a = validation.get("ncbi_table2asn_validation", {})
    val_lines = []
    for name in t2a.get("output_files", []) or []:
        if name.endswith(".val"):
            val_lines += (package / name).read_text(errors="replace").splitlines()[:40]
    return {"sample": sample, "package_version": "revision " + revision["id"] if revision else "original",
            "package_files": sorted(p.name for p in package.iterdir() if p.is_file()),
            "submission_readiness": validation.get("submission_readiness"),
            "missing_required_metadata": validation.get("missing_required_metadata"),
            "phagemine_errors": pre.get("errors", [])[:20], "phagemine_warnings": (pre.get("warnings") or [])[:20],
            "table2asn": {"state": t2a.get("state"), "message": t2a.get("message"), "val_messages": val_lines},
            "corrections_applied": _json(ctx.project.contained(revision["corrections"]), {}) if revision else {},
            "note": "PhageMine never submits to NCBI. The researcher uploads the package through the NCBI Submission Portal."}


@tool("classify_ncbi_errors", "read",
      "Classify pasted NCBI validator output or GenBank staff messages into: FIXABLE_LABEL (product/note/partial/locus tag), FIXABLE_METADATA (source details the researcher must supply), NEEDS_REANALYSIS (coordinates, frames, stops), or UNKNOWN.",
      {"error_text": {"type": "string"}}, ["error_text"])
def classify_ncbi_errors(ctx: ToolContext, error_text: str) -> dict:
    return classify(error_text)


# -------------------------------------------------------------- job tools ----

@tool("start_annotation", "job",
      "Start PhageMine annotation of uploaded FASTA genomes, one job per file. Returns job IDs immediately; it does not return results. Samples already annotated are skipped unless rerun is true.",
      {"files": {"type": "array", "items": {"type": "string"}, "minItems": 1},
       "threads": {"type": "integer", "minimum": 1, "maximum": 64},
       "rerun": {"type": "boolean"}}, ["files"])
def start_annotation(ctx: ToolContext, files: list[str], threads: int | None = None, rerun: bool = False) -> dict:
    project = ctx.project
    threads = int(threads or min(4, os.cpu_count() or 1))
    completed = project.annotation_runs()
    active = {j.get("sample"): j for j in project.jobs() if j.get("kind") == "annotation" and j.get("status") not in jobs.FINAL}
    started, skipped = [], []
    for file in files:
        path = project.input_path(file)
        if path.suffix.lower() not in FASTA_SUFFIXES:
            raise ValueError(f"{file} is not a FASTA file ({', '.join(sorted(FASTA_SUFFIXES))})")
        sample = path.stem
        if sample in active:
            skipped.append({"file": file, "sample": sample, "reason": "already running", "job_id": active[sample]["id"]})
            continue
        if sample in completed and not rerun:
            skipped.append({"file": file, "sample": sample, "reason": "already annotated", "job_id": completed[sample]["id"]})
            continue
        step = ["phagemine", "annotate", f"inputs/{path.name}", "--output", "runs/{job_id}", "--threads", str(threads)]
        if demo_mode():
            step += ["--gene-predictor", "demo", "--mock-evidence"]
        record = jobs.create_job(project, "annotation", [step], output="runs/{job_id}",
                                 meta={"sample": sample, "input": path.name, "demo": demo_mode()})
        started.append({"file": file, "sample": sample, "job_id": record["id"], "status": record["status"]})
    return {"started": started, "skipped": skipped,
            "note": "Jobs run in the background. Results are available only after status is COMPLETE."}


@tool("start_comparison", "job",
      "Compare 2-10 annotated genomes in both directions. Accepts completed sample names and/or IDs of annotation jobs still running (it waits for them). Returns a job ID.",
      {"samples": {"type": "array", "items": {"type": "string"}},
       "annotation_job_ids": {"type": "array", "items": {"type": "string"}},
       "allow_exact_only": {"type": "boolean", "description": "Proceed without MMseqs2 (homologs only from shared PHROG assignments and identical sequences). Works only after the researcher presses Approve on the request this tool creates."}})
def start_comparison(ctx: ToolContext, samples: list[str] | None = None, annotation_job_ids: list[str] | None = None,
                     allow_exact_only: bool = False) -> dict:
    project = ctx.project
    run_dirs, depends = {}, []
    for name in samples or []:
        run_dirs[name] = str(project.run_dir(name).relative_to(project.root))
    for job_id in annotation_job_ids or []:
        record = project.job(job_id)
        if record.get("kind") != "annotation":
            raise ValueError(f"{job_id} is not an annotation job")
        if record.get("status") in {"FAILED", "CANCELLED"}:
            raise ValueError(f"Annotation job {job_id} {record['status'].lower()}; it cannot be compared")
        run_dirs[record["sample"]] = record["output"]
        if record.get("status") != "COMPLETE":
            depends.append(job_id)
    if not 2 <= len(run_dirs) <= 10:
        raise ValueError("A comparison needs 2 to 10 genomes")
    mmseqs = shutil.which("mmseqs")
    consent = None
    if not mmseqs:
        consent = approvals.unused_approved(project, "exact_only_comparison") if allow_exact_only else None
        if consent is None:
            if not any(r["kind"] == "exact_only_comparison" for r in approvals.pending(project)):
                approvals.request(project, "exact_only_comparison", {"samples": sorted(run_dirs)},
                                  "Compare without MMseqs2? Homologs will come only from shared PHROG assignments "
                                  "and identical sequences, so shared-gene counts will be a lower bound.")
            raise ValueError("MMseqs2 is not installed. This needs the researcher's approval: an Approve button is now "
                             "shown in the chat. Explain the limitation and ask them to approve or install MMseqs2; "
                             "after they approve, call start_comparison again with allow_exact_only=true.")
        approvals.consume(project, consent["id"])
    names = list(run_dirs)
    directions, steps = [], []
    for query in names:
        others = [run_dirs[n] for n in names if n != query]
        output = f"runs/{{job_id}}/{safe_name(query)}"
        directions.append({"query": query, "output": output})
        step = ["phagemine", "compare", run_dirs[query], *others, "--output", output]
        if mmseqs:
            step += ["--mmseqs", mmseqs]
        steps.append(step)
    record = jobs.create_job(project, "comparison", steps, output="runs/{job_id}", depends_on=depends,
                             meta={"samples": names, "run_dirs": run_dirs, "directions": directions,
                                   "mmseqs2_available": bool(mmseqs),
                                   "exact_only_approval": consent["id"] if consent else None})
    return {"job_id": record["id"], "samples": names, "waiting_for": depends, "status": record["status"],
            "mmseqs2_available": bool(mmseqs)}


@tool("cancel_job", "job", "Cancel a queued or running job.", {"job_id": {"type": "string"}}, ["job_id"])
def cancel_job(ctx: ToolContext, job_id: str) -> dict:
    return {"job": _job_brief(jobs.cancel(ctx.project, job_id))}


@tool("make_comparison_figure", "job",
      "Draw a gene map of compared genomes, colouring proteins with a homolog in another genome and linking homolog pairs. Saves PNG, SVG and the source table.",
      {"comparison_job_id": {"type": "string", "description": "Omit to use the latest completed comparison"}})
def make_comparison_figure(ctx: ToolContext, comparison_job_id: str | None = None) -> dict:
    from .figures import comparison_map
    record = _comparison(ctx, comparison_job_id)
    links = _links(ctx, record)
    shared = {(l["sample_a"], l["protein_a"]) for l in links} | {(l["sample_b"], l["protein_b"]) for l in links}
    genomes = []
    for sample in record["samples"]:
        run = ctx.project.contained(record["run_dirs"][sample])
        rows = _json(run / "evidence.json", [])
        genome_length = sum(len(line.strip()) for line in (run / "analysis_genome.fasta").read_text().splitlines()
                            if line and not line.startswith(">"))
        genomes.append({"sample": sample, "length": genome_length,
                        "genes": [{"protein_id": r["protein_id"], "start": r["start"], "end": r["end"],
                                   "strand": r.get("strand", "+"), "shared": (sample, r["protein_id"]) in shared} for r in rows]})
    adjacent = [l for l in links if abs(record["samples"].index(l["sample_a"]) - record["samples"].index(l["sample_b"])) == 1]
    base = ctx.project.root / "figures" / f"comparison_{record['id']}"
    title = "Gene map: " + " vs ".join(record["samples"])
    if not record.get("mmseqs2_executed"):
        title += " (without MMseqs2: PHROG anchors and identical sequences only)"
    paths = comparison_map(genomes, adjacent, base, title)
    return {"files": [p.name for p in paths], "comparison_job_id": record["id"],
            "shared_counts": {g["sample"]: sum(x["shared"] for x in g["genes"]) for g in genomes},
            "gene_counts": {g["sample"]: len(g["genes"]) for g in genomes}}


@tool("write_summary", "job",
      "Save a plain-language summary for the researcher as a Markdown file. Every protein ID and figure in it must come from tool results in this conversation; it is checked before saving.",
      {"title": {"type": "string"}, "markdown": {"type": "string"}}, ["title", "markdown"])
def write_summary(ctx: ToolContext, title: str, markdown: str) -> dict:
    problems = ctx.grounding(markdown) if ctx.grounding else []
    if problems:
        raise ValueError("Summary not saved; it is not grounded in tool results: " + "; ".join(problems))
    runs = ctx.project.annotation_runs()
    footer = ["", "---", "", f"Generated by the PhageMine assistant on {now()}.",
              "Sources: " + (", ".join(f"{s} (job {j['id']}, PhageMine {j.get('pipeline_version', 'unknown')})" for s, j in runs.items()) or "none"),
              "Computational predictions are not experimental confirmation."]
    if demo_mode():
        footer.append("DEMO MODE: mock evidence; not biological results.")
    name = f"{safe_name(title.lower().replace(' ', '-'))[:50]}-{uuid.uuid4().hex[:4]}.md"
    path = ctx.project.root / "summaries" / name
    path.write_text(f"# {title}\n\n{markdown.strip()}\n" + "\n".join(footer) + "\n")
    return {"file": name}


def _latest_revision(ctx: ToolContext, sample: str) -> dict | None:
    done = [j for j in ctx.project.jobs() if j.get("kind") == "ncbi_revision"
            and j.get("sample") == sample and j.get("status") == "COMPLETE"]
    return done[-1] if done else None


def _original_metadata(run: Path) -> dict:
    data = _json(run / "genbank_submission" / "submission_metadata.json", {})
    return {k: v for k, v in data.items() if v not in (None, "", [])}


def _approval_payload(preview: dict) -> dict:
    """Exactly what the researcher approves; any change afterwards voids the approval."""
    return {key: preview[key] for key in ("id", "sample", "source_run", "input_sha256", "corrections")}


@tool("preview_ncbi_corrections", "job",
      "Prepare (not apply) corrections to a sample's NCBI submission package: submission metadata the researcher typed, and product/note/locus_tag changes. Validates them, shows the researcher an Approve/Decline request, and returns a preview ID. Coordinates and partial flags cannot be changed here.",
      {"sample": {"type": "string"},
       "metadata": {"type": "object", "description": "Submission fields such as organism, host, collection_date, geographic_location; values must be the researcher's own words"},
       "annotations": {"type": "object", "description": "{protein_id: {product?, note?, locus_tag?, reason}}"},
       "addresses_errors": {"type": "string", "description": "The NCBI messages these corrections respond to"}},
      ["sample"])
def preview_ncbi_corrections(ctx: ToolContext, sample: str, metadata: dict | None = None,
                             annotations: dict | None = None, addresses_errors: str | None = None) -> dict:
    run = ctx.project.run_dir(sample)
    previous = _latest_revision(ctx, sample)
    base = _json(ctx.project.contained(previous["corrections"]), {}) if previous else {}
    metadata, annotations = dict(metadata or {}), dict(annotations or {})
    if not metadata and not annotations:
        raise ValueError("No corrections were given")
    problems = check_metadata(metadata, ctx.user_text)
    records = {r["protein_id"]: r for r in _records(ctx, sample)}
    changes = []
    for pid, change in annotations.items():
        if pid not in records:
            problems.append(f"{pid} is not a protein in {sample}")
            continue
        if not isinstance(change, dict):
            problems.append(f"correction for {pid} must be an object")
            continue
        curation = set(change) & CURATION_FIELDS
        if curation:
            problems.append(f"{pid}: {', '.join(sorted(curation))} is a curation decision about the gene model, "
                            "not an NCBI wording fix; review the gene model instead")
        unsupported = set(change) - REVISABLE_FIELDS - CURATION_FIELDS - {"reason", "source"}
        if unsupported:
            problems.append(f"{pid}: {', '.join(sorted(unsupported))} cannot be revised")
        current = records[pid].get("product")
        if "product" in change:
            reason = check_product(str(change["product"]), current)
            if reason:
                problems.append(f"{pid}: product '{change['product']}' refused, {reason}")
        for key in REVISABLE_FIELDS & set(change):
            changes.append({"protein_id": pid, "field": key, "before": current if key == "product" else None, "after": change[key]})
    if problems:
        raise ValueError("Corrections refused: " + " | ".join(problems))
    # Merge onto the original package metadata and earlier corrections, field by field.
    merged_annotations = {pid: dict(fields) for pid, fields in base.get("annotations", {}).items()}
    for pid, change in annotations.items():
        entry = merged_annotations.setdefault(pid, {})
        entry.update({k: v for k, v in change.items() if k in REVISABLE_FIELDS})
        entry["reason"] = change.get("reason") or entry.get("reason") or "NCBI validation response"
        entry["source"] = "researcher-approved assistant proposal"
    merged = {"metadata": {**_original_metadata(run), **base.get("metadata", {}), **metadata},
              "annotations": merged_annotations}
    preview_id = f"preview-{uuid.uuid4().hex[:8]}"
    manifest = _json(run / "run_manifest.json", {})
    preview = {"id": preview_id, "sample": sample, "created": now(), "corrections": merged,
               "source_run": str(run.relative_to(ctx.project.root)), "input_sha256": manifest.get("input_sha256"),
               "addresses_errors": addresses_errors, "applied": False}
    request = approvals.request(ctx.project, "ncbi_corrections", _approval_payload(preview),
                                f"Apply these corrections to the {sample} NCBI package?")
    preview["approval_id"] = request["id"]
    ctx.project.write_state(f"ncbi/{preview_id}.json", preview)
    return {"preview_id": preview_id, "approval_id": request["id"], "sample": sample,
            "metadata_changes": [{"field": k, "after": v} for k, v in metadata.items()],
            "annotation_changes": changes, "builds_on": previous["id"] if previous else "original package",
            "next_step": "Explain the changes. The researcher approves or declines with the buttons shown in the chat; "
                         "only then can apply_ncbi_corrections run."}


@tool("apply_ncbi_corrections", "job",
      "Apply a previewed set of corrections that the researcher approved with the Approve button. Runs 'phagemine revise' to write a new submission package and re-validate it. Returns a job ID. A typed message is not approval.",
      {"preview_id": {"type": "string"}}, ["preview_id"])
def apply_ncbi_corrections(ctx: ToolContext, preview_id: str) -> dict:
    if not re.fullmatch(r"preview-[0-9a-f]{8}", preview_id or ""):
        raise ValueError("Unknown preview ID")
    preview = ctx.project.read_state(f"ncbi/{preview_id}.json", None)
    if not preview:
        raise ValueError("Unknown preview ID")
    if preview.get("applied"):
        raise ValueError("This preview was already applied")
    approvals.check(ctx.project, preview["approval_id"], "ncbi_corrections", _approval_payload(preview))
    run = ctx.project.contained(preview["source_run"])
    if _json(run / "run_manifest.json", {}).get("input_sha256") != preview["input_sha256"]:
        raise ValueError("The annotation this preview was built from has changed; prepare a new preview")
    corrections_file = f"ncbi/{preview_id}.corrections.json"
    ctx.project.write_state(corrections_file, preview["corrections"])
    step = ["phagemine", "revise", preview["source_run"], "--corrections", corrections_file, "--output", "runs/{job_id}"]
    approvals.consume(ctx.project, preview["approval_id"])
    record = jobs.create_job(ctx.project, "ncbi_revision", [step], output="runs/{job_id}",
                             meta={"sample": preview["sample"], "corrections": corrections_file, "preview_id": preview_id})
    preview["applied"] = record["id"]
    ctx.project.write_state(f"ncbi/{preview_id}.json", preview)
    return {"job_id": record["id"], "status": record["status"],
            "note": "This job writes and re-validates a package. Its completion does not mean the submission passes NCBI review."}
