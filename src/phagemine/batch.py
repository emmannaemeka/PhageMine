"""Sequential, checkpointed orchestration of the single-genome workflow."""
from __future__ import annotations

import csv
import html
import copy
import hashlib
import json
import os
import re
import traceback
import shutil
import pickle
from concurrent.futures import ThreadPoolExecutor, as_completed
import time
from datetime import datetime, timezone
from pathlib import Path

from .gene_prediction import create_predictor
from .pipeline import run
from .progress import ProgressReporter
from .resume import recover_evidence_complete
from .reporting import prefix_checkpoint_artifacts
from .io import read_fasta, checksum
from .models import Protein
from .genome_representation import GenomeRepresentation
from .sequencing_provenance import SequencingProvenance
from .fusion import classify_proteins
from .context import build_context
from .mining import mine, ranked_candidates
from .quality import assess
from .reporting import write_outputs
from .genbank import write_package
from .pooled import deduplicate_proteins, sequence_sha256
from .pfam import PfamHMMAdapter
from .vog import VOGHMMAdapter
from .swissprot import SwissProtEvidenceAdapter
from .phrogs import PHROGSMMseqsAdapter
from .resources import EvidenceResourceManager, ResourceType
from .hmmer_chunking import run_chunked
from .discovery import build_discovery_outputs

FASTA_EXTENSIONS = {".fasta", ".fa", ".fna"}
REQUIRED_RUN = ("run_manifest.json", "evidence.json", "functional_classification.json",
                "genomic_context.json", "modules.json")
SUMMARY_FIELDS = ("sample_id", "input_file", "status", "error_message", "genome_length",
                  "predicted_proteins", "KNOWN_FUNCTION", "PROBABLE_FUNCTION",
                  "FUNCTIONAL_CLASS_ONLY", "CONSERVED_UNKNOWN", "CONFLICTING_EVIDENCE",
                  "UNRESOLVED", "module_count", "output_directory")


def _comparative_resources() -> tuple[str | None, dict | None]:
    """Resolve only fully validated comparative resources for automatic use."""
    resources = EvidenceResourceManager().validate_all(check_checksum=True)
    pmfdb = next((item for item in resources if item.get("resource_type") == "PMFDB" and item.get("status") == "READY"), None)
    inphared = next((item for item in resources if item.get("resource_type") == "INPHARED_GENOMES" and item.get("status") == "READY"), None)
    return (pmfdb.get("path") if pmfdb else None), inphared


def pooled_batch(input_paths, output, *, profile="full", threads=1, gene_predictor="phanotate", phanotate=None, progress=None, resume_existing=False):
    """Run gene prediction per genome, then each requested adapter once on exact representatives."""
    from .preflight import preflight_profile
    resolution = preflight_profile(profile)
    progress = progress or ProgressReporter(quiet=False)
    overall_started = time.monotonic()
    def stage_start(name):
        progress._write(f"[Discovery] RUNNING: {name}")
        return time.monotonic()
    def stage_done(name, started, detail=""):
        suffix = f"; {detail}" if detail else ""
        progress._write(f"[Discovery] DONE: {name} ({time.monotonic()-started:.1f}s){suffix}")
    project = Path(output); project.mkdir(parents=True, exist_ok=True)
    checkpoint_path = project / "pooled_checkpoints.json"
    checkpoints = json.loads(checkpoint_path.read_text()) if resume_existing and checkpoint_path.is_file() else {}
    def mark(stage, state="COMPLETE", detail=None):
        checkpoints[stage] = {"state": state, "detail": detail, "time": datetime.now(timezone.utc).isoformat()}
        checkpoint_path.write_text(json.dumps(checkpoints, indent=2, sort_keys=True))
    state_path = project / "pooled_state.pkl"
    fingerprint = {str(Path(p)): checksum(p) for p in input_paths}
    reusable = resume_existing and state_path.is_file() and checkpoints.get("GENE_PREDICTION", {}).get("state") == "COMPLETE"
    state = pickle.loads(state_path.read_bytes()) if reusable else {}
    predictor_factory = lambda: create_predictor(gene_predictor, phanotate)
    def predict(path):
        gid, seq = read_fasta(path)
        proteins = predictor_factory().predict(gid, seq, path)
        if not proteins: raise ValueError(f"No proteins predicted for {path}")
        return path, gid, seq, proteins
    pred_started = stage_start("Gene prediction")
    if reusable and state.get("input_fingerprint") == fingerprint:
        predicted = state["predicted"]
    else:
        with ThreadPoolExecutor(max_workers=max(1, int(threads))) as pool:
            futures = {pool.submit(predict, path): path for path in input_paths}
            predicted = []
            for index, future in enumerate(as_completed(futures), 1):
                predicted.append(future.result())
                progress._write(f"[Discovery] Gene prediction {index}/{len(input_paths)} complete: {futures[future].name}")
            predicted.sort(key=lambda item: str(item[0]))
        state = {"input_fingerprint": fingerprint, "predicted": predicted}
        mark("GENE_PREDICTION")
    stage_done("Gene prediction", pred_started, f"{len(predicted)}/{len(input_paths)} genomes")
    pool_started = stage_start("Protein pooling")
    all_proteins = [p for _, _, _, ps in predicted for p in ps]
    stage_done("Protein pooling", pool_started, f"{len(all_proteins)} protein occurrences")
    dedup_started = stage_start("Exact protein deduplication")
    if reusable and state.get("input_fingerprint") == fingerprint and state.get("representatives"):
        representatives, occurrences = state["representatives"], state["occurrences"]
    else:
        representatives, occurrences = deduplicate_proteins(all_proteins)
        state.update(representatives=representatives, occurrences=occurrences)
    state_path.write_bytes(pickle.dumps(state))
    mark("POOL_DEDUP")
    reduction = 100.0 * (len(all_proteins)-len(representatives)) / len(all_proteins) if all_proteins else 0.0
    stage_done("Exact protein deduplication", dedup_started, f"{len(representatives)} unique; {reduction:.1f}% reduction")
    by_digest = {sequence_sha256(p.sequence): p for p in representatives}
    resource_map = {r["resource_type"]: r for r in resolution.get("resources", [])}
    evidence_adapters = []
    adapters = []
    if profile == "standard":
        kinds = ["PHROGS"]
    elif profile == "full":
        kinds = ["PFAM", "VOGDB", "SWISSPROT", "PHROGS"]
    else:
        kinds = []
    for kind in kinds:
        evidence_started = stage_start(kind)
        r = resource_map[kind]; prov = r.get("provenance") or {}
        # CLI-resolved explicit resources use adapter-default identity unless a
        # user supplies a version; keep pooled and single-genome semantics equal.
        if kind == "PFAM": adapter = PfamHMMAdapter(r["path"], evalue_threshold=None, coverage_threshold=None, threads=threads)
        elif kind == "VOGDB": adapter = VOGHMMAdapter(r["path"], prov.get("annotations_path"), database_version="unknown", threads=threads)
        elif kind == "SWISSPROT": adapter = SwissProtEvidenceAdapter(r["path"], prov.get("metadata_path"), database_version="unknown", threads=threads)
        else: adapter = PHROGSMMseqsAdapter(r["path"], prov.get("annotations_path"), database_version="unknown", threads=threads)
        adapters.append((kind, adapter))
        if reusable and checkpoints.get(kind, {}).get("state") == "COMPLETE" and state.get("evidence", {}).get(kind) is not None:
            result = state["evidence"][kind]
        else:
            result = adapter.analyze(representatives)
            state.setdefault("evidence", {})[kind] = result
            state_path.write_bytes(pickle.dumps(state))
            mark(kind)
        evidence_adapters.append({"adapter": result.adapter, "status": result.status, "provenance": result.provenance, "pooled_execution": True})
        stage_done(kind, evidence_started)
        rep_by_id = {p.protein_id: p for p in representatives}
        for evidence in result.evidence:
            pid = evidence.provenance.get("protein_id") or evidence.metrics.get("query_protein_id")
            if pid in rep_by_id:
                evidence.provenance["pooled_execution"] = True
                rep_by_id[pid].evidence.append(evidence)
    evidence_by_digest = {sequence_sha256(p.sequence): list(p.evidence) for p in representatives}
    for p in all_proteins:
        p.evidence = [type(e)(**{**e.__dict__, "provenance": {**e.provenance, "pooled_execution": True}}) for e in evidence_by_digest[sequence_sha256(p.sequence)]]
    mark("EVIDENCE_REMAP")
    map_started = stage_start("Evidence mapping/classification")
    rows=[]
    for path, gid, seq, proteins in predicted:
        destination = project / _sample_id(path); destination.mkdir(parents=True, exist_ok=True)
        local = [p for p in all_proteins if p.genome_id == gid]
        classifications = classify_proteins(local); contexts, modules = build_context(local, classifications)
        mine(local); candidates = ranked_candidates(local); quality = assess(seq, local)
        manifest = {"pipeline":"PhageMine", "command":"batch", "pooled_execution":True,
                    "threads":threads, "evidence_profile":profile, "input":str(path), "input_sha256":checksum(path),
                    "gene_caller":{"name":predictor_factory().name,"version":predictor_factory().version(),"parameters":predictor_factory().parameters()},
                    "evidence_adapters":evidence_adapters, "pooled_proteins":{"total_occurrences":len(all_proteins),"unique_sequences":len(representatives),"occurrence_map":occurrences},
                    "discovery_ranking":{"status":"RANKED"}}
        rep = GenomeRepresentation.original(gid, seq)
        write_outputs(destination, rep, SequencingProvenance(), local, candidates, manifest, quality, path, classifications, contexts, modules)
        write_package(destination, gid, seq, local, manifest, sequencing_provenance=SequencingProvenance())
        rows.append({"sample_id":_sample_id(path),"status":"SUCCESS","predicted_proteins":len(local),"output_directory":str(destination)})
    stage_done("Evidence mapping/classification", map_started)
    progress._write(f"[Discovery] DONE: Discovery analysis ({time.monotonic()-overall_started:.1f}s); genomes={len(predicted)} proteins={len(all_proteins)} unique={len(representatives)} reduction={reduction:.1f}%")
    (project/"pooled_manifest.json").write_text(json.dumps({"pooled_execution":True,"threads":threads,"profile":profile,"total_proteins":len(all_proteins),"unique_proteins":len(representatives),"occurrences":occurrences,"samples":rows},indent=2,sort_keys=True))
    return rows


def discovery_from_annotation(annotation_root, output):
    """Build the discovery view from completed independent annotations.

    No gene caller or evidence adapter is invoked; annotation artifacts are
    copied verbatim and provenance records the reuse.
    """
    source = Path(annotation_root); destination = Path(output)
    destination.mkdir(parents=True, exist_ok=True)
    samples = []
    for sample in sorted(p for p in source.iterdir() if p.is_dir() and (p / "proteins.faa").is_file()):
        target = destination / sample.name
        if target.exists(): shutil.rmtree(target)
        shutil.copytree(sample, target)
        samples.append({"sample_id": sample.name, "status": "REUSED_ANNOTATION",
                        "output_directory": str(target)})
    (destination / "pooled_manifest.json").write_text(json.dumps({
        "pooled_execution": False, "discovery_from_annotation": True,
        "evidence_reused": True, "samples": samples}, indent=2, sort_keys=True))
    return samples


def _sample_id(path: Path) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", path.stem).strip("._-") or "sample"


def discover_inputs(input_dir: str | Path, recursive: bool = False) -> list[Path]:
    root = Path(input_dir)
    if not root.is_dir():
        raise ValueError(f"Batch input directory does not exist: {root}")
    paths = root.rglob("*") if recursive else root.iterdir()
    return sorted((p for p in paths if p.is_file() and p.suffix.lower() in FASTA_EXTENSIONS),
                  key=lambda p: str(p.relative_to(root)).lower())


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _valid_output(path: Path, input_path: Path | None = None) -> bool:
    if not all((path / name).is_file() for name in REQUIRED_RUN):
        return False
    if input_path is None:
        return True
    status_path = path / "sample_status.json"
    if not status_path.is_file():
        return False
    try:
        status = json.loads(status_path.read_text())
        return status.get("input_sha256") == _sha256(input_path) and status.get("status") in {"SUCCESS", "REUSED"}
    except (OSError, ValueError, TypeError):
        return False


def _atomic_json(path: Path, value) -> None:
    temp = path.with_name(path.name + ".tmp")
    temp.write_text(json.dumps(value, indent=2, sort_keys=True))
    os.replace(temp, path)


def _write_batch_presentation(project: Path, rows: list[dict]) -> None:
    """Create compact cohort annotation presentation from completed sample artifacts."""
    completed = [r for r in rows if r.get("status") in {"SUCCESS", "REUSED"}]
    states = ["KNOWN_FUNCTION", "PROBABLE_FUNCTION", "FUNCTIONAL_CLASS_ONLY", "CONSERVED_UNKNOWN", "UNRESOLVED", "CONFLICTING_EVIDENCE"]
    totals = {s: sum(int(r.get(s) or 0) for r in completed) for s in states}; total = sum(totals.values())
    figures = project / "figures"; figures.mkdir(exist_ok=True)
    data = project / "figure_data"; data.mkdir(exist_ok=True)
    with (data / "genome_functional_states.tsv").open("w") as h:
        h.write("genome\tfunctional_state\tcount\n")
        for r in completed:
            for s in states: h.write(f"{r['sample_id']}\t{s}\t{int(r.get(s) or 0)}\n")
    try:
        from .figures import _plotting, _save
        plt = _plotting(figures)
        labels = ["characterized/probable", "conserved unknown", "unresolved", "other/ambiguous"]
        values = [totals["KNOWN_FUNCTION"] + totals["PROBABLE_FUNCTION"] + totals["FUNCTIONAL_CLASS_ONLY"], totals["CONSERVED_UNKNOWN"], totals["UNRESOLVED"], totals["CONFLICTING_EVIDENCE"]]
        fig, ax = plt.subplots(figsize=(7.5,4.5)); ax.bar(labels, values, color=["#4c78a8","#f58518","#777777","#b2182b"]); ax.set_title("Functional landscape across annotated genomes"); ax.set_ylabel("Proteins"); ax.tick_params(axis="x", rotation=25); fig.tight_layout(); _save(fig, figures / "genome_by_functional_state"); plt.close(fig)
    except Exception:
        pass
    lines = ["<html><head><meta charset='utf-8'><style>body{font-family:Arial;max-width:1400px;margin:auto;padding:2em}table{border-collapse:collapse;width:100%}td,th{border:1px solid #bbb;padding:.3em}th{background:#eee}</style></head><body><h1>PHAGEMINE BATCH ANNOTATION</h1>", f"<p>Genomes analysed: <b>{len(rows)}</b>; completed successfully: <b>{len(completed)}</b>; GenBank-ready genomes: <b>{sum((Path(r['output_directory'])/'genbank_submission').is_dir() for r in completed)}</b></p>", "<h2>Functional-state totals</h2><table><tr><th>State</th><th>Count</th><th>Percent</th></tr>"]
    for s, value in totals.items(): lines.append(f"<tr><td>{s}</td><td>{value}</td><td>{value/total:.1%}</td></tr>" if total else f"<tr><td>{s}</td><td>0</td><td>0%</td></tr>")
    lines.append("</table><h2>Per-genome annotations</h2>")
    for r in completed:
        sample = Path(r["output_directory"]); lines.append(f"<h3>{html.escape(r['sample_id'])}</h3><table><tr><th>Protein</th><th>Coordinates</th><th>Strand</th><th>Classification</th><th>Proposed function</th><th>Confidence</th></tr>")
        try: records=json.loads((sample/'functional_classification.json').read_text())
        except Exception: records=[]
        for c in records:
            pid=c.get('protein_id')
            if not pid: continue
            lines.append(f"<tr><td><a href='{html.escape(str(sample.relative_to(project) / 'protein_details' / (pid+'.html')))}'>{html.escape(pid)}</a></td><td>{c.get('start')}-{c.get('end')}</td><td>{html.escape(str(c.get('strand') or ''))}</td><td>{html.escape(str(c.get('functional_state') or 'UNRESOLVED'))}</td><td>{html.escape(str(c.get('proposed_function') or 'Function unresolved'))}</td><td>{html.escape(str(c.get('confidence') or 'NONE'))}</td></tr>")
        lines.append("</table>")
    lines.append("<h2>Downloads</h2><p><a href='batch_summary.tsv'>Batch summary TSV</a> · <a href='figure_data/genome_functional_states.tsv'>Figure source data</a></p></body></html>")
    (project / "batch_annotation_report.html").write_text("".join(lines))
    (project / "batch_annotation_report.md").write_text(f"# PHAGEMINE BATCH ANNOTATION\n\nGenomes analysed: {len(rows)}\n\nCompleted: {len(completed)}\n\nFunctional states: {json.dumps(totals, sort_keys=True)}\n")


def batch(input_dir: str | Path, output: str | Path, recursive=False, resume_existing=False,
          fail_fast=False, gene_predictor="phanotate", phanotate=None, progress=None,
          reconcile_orfs=False, prodigal=None, threads: int = 1, evidence_profile: str = "core",
          mode: str = "annotate") -> list[dict]:
    root = Path(input_dir).resolve()
    from .preflight import preflight_profile
    profile_resolution = preflight_profile(evidence_profile)
    if mode not in {"annotate", "discover", "both"}:
        raise ValueError("mode must be annotate, discover, or both")
    project = Path(output).resolve()
    project.mkdir(parents=True, exist_ok=True)
    inputs = discover_inputs(root, recursive)
    if mode == "discover":
        rows = pooled_batch(inputs, project, profile=evidence_profile, threads=threads,
                            gene_predictor=gene_predictor, phanotate=phanotate, progress=progress,
                            resume_existing=resume_existing)
        # Pooled annotation is an input to the cohort-level biological workflow;
        # always materialize PMFs, recurrence, context, PMFDB validation, ranking,
        # reports, and figures from the completed per-genome artifacts.
        sample_dirs = [project / _sample_id(path) for path in inputs]
        pmfdb, inphared = _comparative_resources()
        build_discovery_outputs(sample_dirs, project, progress=progress, pmfdb=pmfdb, inphared=inphared,
                                resume_existing=resume_existing, source_mode="discover")
        return rows
    if mode == "both":
        batch(input_dir, project / "annotation", recursive, resume_existing, fail_fast,
              gene_predictor, phanotate, progress, reconcile_orfs, prodigal, threads,
              evidence_profile, mode="annotate")
        samples = discovery_from_annotation(project / "annotation", project / "discovery")
        sample_dirs = [project / "discovery" / _sample_id(path) for path in inputs]
        pmfdb, inphared = _comparative_resources()
        build_discovery_outputs(sample_dirs, project / "discovery", progress=progress, pmfdb=pmfdb, inphared=inphared,
                                resume_existing=resume_existing, evidence_reused=True,
                                source_mode="both")
        return samples
    progress = progress or ProgressReporter(quiet=True)
    started = datetime.now(timezone.utc).isoformat()
    rows: list[dict] = []
    manifest = {"batch_schema_version": "1.2", "pipeline": "PhageMine", "started_at": started,
                "ended_at": None, "input_directory": str(root), "recursive": recursive,
                "inputs": [str(p) for p in inputs], "samples": [], "failures": [],
                "configuration": {"gene_predictor": gene_predictor, "resume_existing": resume_existing,
                                   "fail_fast": fail_fast, "threads": threads,
                                   "evidence_profile": evidence_profile,
                                   "mode": mode,
                                   "evidence_profile_resolution": profile_resolution,
                                   "pooled_execution": False,
                                   "pooled_execution_note": "Reserved until pooled adapter equivalence is validated."}}

    def write_summary() -> None:
        _atomic_json(project / "batch_summary.json", rows)
        temp = project / "batch_summary.tsv.tmp"
        with temp.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=SUMMARY_FIELDS, delimiter="\t")
            writer.writeheader()
            writer.writerows({key: row.get(key) for key in SUMMARY_FIELDS} for row in rows)
        os.replace(temp, project / "batch_summary.tsv")

    def persist() -> None:
        manifest["samples"] = rows
        manifest["failures"] = [row for row in rows if row["status"] == "FAILED"]
        _atomic_json(project / "batch_manifest.json", manifest)
        write_summary()

    _atomic_json(project / "batch_manifest.json", manifest)
    write_summary()
    used_ids: set[str] = set()

    for index, path in enumerate(inputs, 1):
        sid = _sample_id(path)
        base_sid = sid
        suffix = 2
        while sid in used_ids:
            sid = f"{base_sid}_{suffix}"
            suffix += 1
        used_ids.add(sid)
        destination = project / sid
        preexisting = destination.exists()
        can_reuse = resume_existing and _valid_output(destination, path)
        previous_status = None
        if preexisting and (destination / "sample_status.json").is_file():
            try:
                previous_status = json.loads((destination / "sample_status.json").read_text())
            except (OSError, ValueError):
                previous_status = None
        destination.mkdir(parents=True, exist_ok=True)
        for folder in ("checkpoints", "logs", ".working"):
            (destination / folder).mkdir(exist_ok=True)
        status_path = destination / "sample_status.json"
        input_hash = _sha256(path)
        status = {"sample_id": sid, "input_file": str(path), "input_sha256": input_hash,
                  "status": "RUNNING", "failed_stage": None, "current_stage": None,
                  "checkpoints": {}, "started_at": datetime.now(timezone.utc).isoformat()}
        if previous_status and resume_existing and previous_status.get("input_sha256") == input_hash:
            status["checkpoints"] = previous_status.get("checkpoints", {})
            if previous_status.get("status") == "FAILED" and status["checkpoints"].get("evidence_integration", {}).get("status") == "COMPLETE":
                status["audit_history"] = previous_status.get("audit_history", []) + [{"status": previous_status.get("status"), "failed_stage": previous_status.get("failed_stage"), "error_message": previous_status.get("error_message")}]
        _atomic_json(status_path, status)
        row = {"sample_id": sid, "input_file": str(path), "status": "FAILED", "error_message": None,
               "genome_length": None, "predicted_proteins": None, "module_count": None,
               "output_directory": str(destination)}
        for state in ("KNOWN_FUNCTION", "PROBABLE_FUNCTION", "FUNCTIONAL_CLASS_ONLY",
                      "CONSERVED_UNKNOWN", "CONFLICTING_EVIDENCE", "UNRESOLVED"):
            row[state] = None
        rows.append(row)
        log = destination / "logs" / "run.log"
        log.write_text(f"START {sid} {path}\n")

        def checkpoint(state_name, stage, detail):
            key = stage.replace("/", "_").replace(" ", "_")
            now = datetime.now(timezone.utc).isoformat()
            status["checkpoints"][key] = {"stage": stage, "status": state_name, "time": now, "detail": detail}
            status["current_stage"] = stage if state_name == "RUNNING" else None
            if state_name == "COMPLETE":
                prefix_checkpoint_artifacts(destination, sid)
            _atomic_json(status_path, status)
            with log.open("a") as handle:
                handle.write(f"{state_name} {stage}\n")

        sample_progress = ProgressReporter(stream=progress.stream, quiet=progress.quiet,
                                            no_progress=progress.no_progress, callback=checkpoint)
        progress._write(f"[Batch {index}/{len(inputs)}] RUNNING: {sid}")
        recovering = False
        try:
            if can_reuse:
                row["status"] = "REUSED"
                status["status"] = "REUSED"
                status["reused"] = True
                _atomic_json(status_path, status)
            elif resume_existing and previous_status and previous_status.get("input_sha256") == input_hash and status.get("audit_history") and status["checkpoints"].get("evidence_integration", {}).get("status") == "COMPLETE":
                recovering = True
                recover_evidence_complete(destination, sample_progress)
                row["status"] = "SUCCESS"
                status = json.loads(status_path.read_text())
            elif resume_existing and preexisting and status.get("checkpoints"):
                raise ValueError("checkpointed partial run cannot be resumed by the monolithic single-genome workflow")
            elif preexisting and any(destination.iterdir()):
                raise ValueError(f"output directory exists but is not a valid completed run: {destination}")
            else:
                run(path, destination, command="run", predictor=create_predictor(gene_predictor, phanotate),
                    reconcile_orfs=reconcile_orfs, prodigal=prodigal,
                    progress=sample_progress)
                row["status"] = "SUCCESS"
                status["status"] = "SUCCESS"
                _atomic_json(status_path, status)
            classification = json.loads((destination / "functional_classification.json").read_text())
            modules = json.loads((destination / "modules.json").read_text())
            fasta_lines = [line for line in (destination / "analysis_genome.fasta").read_text().splitlines() if not line.startswith(">")]
            row["genome_length"] = len("".join(fasta_lines))
            row["predicted_proteins"] = len(classification)
            row["module_count"] = len(modules)
            for state_name in ("KNOWN_FUNCTION", "PROBABLE_FUNCTION", "FUNCTIONAL_CLASS_ONLY",
                               "CONSERVED_UNKNOWN", "CONFLICTING_EVIDENCE", "UNRESOLVED"):
                row[state_name] = sum(item.get("functional_state") == state_name for item in classification)
            with log.open("a") as handle:
                handle.write("Completed single-genome workflow\n")
        except Exception as exc:
            row["error_message"] = str(exc)
            traceback_text = traceback.format_exc()
            failed_stage = status.get("current_stage")
            if recovering:
                failed_stage = "resume validation"
            if "fusion.py" in traceback_text or "classify_protein" in traceback_text:
                failed_stage = "functional classification / evidence fusion"
            status.update({"status": "FAILED", "failed_stage": failed_stage,
                           "error_type": type(exc).__name__, "error_message": str(exc),
                           "finished_at": datetime.now(timezone.utc).isoformat()})
            traceback_path = destination / "logs" / "traceback.log"
            traceback_path.write_text(traceback_text)
            status["traceback_path"] = str(traceback_path)
            _atomic_json(status_path, status)
            with log.open("a") as handle:
                handle.write(f"FAILED {status.get('failed_stage')}: {exc}\n")
            progress._write(f"[Batch {index}/{len(inputs)}] FAILED: {sid}")
            progress._write(f"Stage: {status.get('failed_stage') or 'unknown'}")
            progress._write(f"Reason: {exc}")
            progress._write(f"Resume available from: {destination}")
            persist()
            if fail_fast:
                break
            continue
        progress._write(f"[Batch {index}/{len(inputs)}] DONE: {sid}")
        persist()

    manifest["ended_at"] = datetime.now(timezone.utc).isoformat()
    persist()
    _write_batch_presentation(project, rows)
    return rows
