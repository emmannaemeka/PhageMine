"""Fail-closed continuation of a completed PhageMine run."""
from __future__ import annotations

import hashlib
import csv
import json
import os
import shutil
import tempfile
from datetime import datetime, timezone
from dataclasses import fields
from pathlib import Path
from typing import Any

from .genome_representation import GenomeRepresentation, Orientation, Rotation, Topology, TransformEvent
from .genbank import write_package
from .io import checksum, read_fasta
from .mining import mine, ranked_candidates
from .models import Evidence, EvidenceLevel, Protein
from .phrogs import PHROGSMMseqsAdapter, PHROGSPyHMMERAdapter, merge_phrogs_evidence
from .progress import ProgressReporter
from .quality import assess
from .reporting import update_comparative_report, write_outputs
from .resources import EvidenceResourceManager, ResourceType
from .sequencing_provenance import SequencingProvenance
from .fusion import attach_gene_call_assessments, classify_proteins
from .context import build_context


REUSED = ("input/genome validation", "genome representation", "gene prediction", "Pfam", "VOGDB", "Swiss-Prot")
RESUME_STAGES = REUSED + ("PHROGs", "evidence integration", "candidate ranking/mining", "QC/report generation", "GenBank package")


def _load_gene_call_assessments(root: Path) -> list[dict[str, Any]]:
    path = root / "gene_call_confidence.tsv"
    if not path.is_file():
        return []
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _copy_gene_call_artifacts(source: Path, destination: Path) -> None:
    """Preserve caller-review evidence in regenerated result directories."""
    for name in (
        "orf_reconciliation.tsv", "orf_reconciliation.json",
        "alternative_orf_evidence.tsv", "alternative_orf_evidence.json",
        "orf_adjudication.tsv", "orf_adjudication.json",
        "gene_call_confidence.tsv", "gene_calls_for_review.tsv",
    ):
        source_path = source / name
        if source_path.is_file():
            shutil.copy2(source_path, destination / name)


def reclassify(source: str | Path, output: str | Path,
               progress: ProgressReporter | None = None) -> int:
    """Regenerate annotation products from persisted evidence without searches.

    This is intentionally separate from ``resume``: resume acquires missing or
    explicitly refreshed evidence, whereas reclassify never instantiates or
    invokes an evidence adapter.
    """
    source, output = Path(source).resolve(), Path(output).resolve()
    if source == output or output.exists():
        raise ValueError("Reclassify output must be a new, non-existent directory")
    manifest, representation, proteins, _, sequencing = _load_source(source)
    progress = progress or ProgressReporter(quiet=True)
    progress.STAGES = ("persisted evidence validation", "functional classification",
                       "genomic context", "candidate ranking", "report generation",
                       "GenBank package")
    progress.start("persisted evidence validation"); progress.finish(f"{len(proteins)} proteins validated")
    progress.start("functional classification")
    classifications = attach_gene_call_assessments(
        classify_proteins(proteins), _load_gene_call_assessments(source))
    progress.finish("regenerated without database searches")
    progress.start("genomic context")
    context_records, modules = build_context(proteins, classifications)
    progress.finish("regenerated")
    progress.start("candidate ranking")
    mine(proteins); candidates = ranked_candidates(proteins)
    progress.finish("regenerated")
    quality = assess(representation.analysis_sequence, proteins)
    ranking_status = "INSUFFICIENT_EVIDENCE" if candidates and not any(
        e.supports and e.evidence_strength in {"STRONG", "EXPERIMENTAL"}
        for protein in candidates for e in protein.evidence
    ) else "RANKED"
    new_manifest = {
        **manifest,
        "command": "reclassify",
        "input": str(source / "original_input.fasta"),
        "input_sha256": checksum(source / "original_input.fasta"),
        "stage_status": {
            "input/genome validation": "REUSED", "gene prediction": "REUSED",
            "Pfam": "REUSED", "VOGDB": "REUSED", "Swiss-Prot": "REUSED",
            "PHROGs": "REUSED", "functional classification": "RUN",
            "genomic context/modules": "RUN", "candidate ranking/mining": "RUN",
            "QC/report generation": "RUN", "GenBank package": "RUN",
        },
        "reclassification": {
            "source": str(source), "database_searches_run": [],
            "evidence_reused": True, "source_manifest_preserved": True,
        },
        "quality_control": quality,
        "discovery_ranking": {"status": ranking_status,
                              "message": "Candidates re-ranked from persisted evidence."},
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    temp = Path(tempfile.mkdtemp(prefix=f".{output.name}.", dir=str(output.parent)))
    try:
        progress.start("report generation")
        write_outputs(temp, representation, sequencing, proteins, candidates,
                      new_manifest, quality, source / "original_input.fasta",
                      classifications, context_records, modules)
        _copy_gene_call_artifacts(source, temp)
        comparative_source = source / "comparative"
        if comparative_source.is_dir():
            shutil.copytree(comparative_source, temp / "comparative")
            inphared_json = temp / "comparative" / "inphared_nearest_phages.json"
            if inphared_json.is_file():
                update_comparative_report(temp, {"inphared": json.loads(inphared_json.read_text())})
        progress.finish("regenerated")
        progress.start("GenBank package")
        write_package(temp, representation.analysis_sequence_id,
                      representation.analysis_sequence, proteins, new_manifest,
                      sequencing_provenance=sequencing)
        progress.finish("regenerated")
        os.replace(temp, output)
    except Exception:
        shutil.rmtree(temp, ignore_errors=True)
        raise
    return len(proteins)

def checkpoint_reusable(checkpoint: dict, current_provenance: dict | None) -> bool:
    """A previously unavailable resource is never reusable once requested now."""
    if not checkpoint or checkpoint.get("status") not in {"REAL", "COMPLETE"}:
        return False
    previous = checkpoint.get("provenance", {})
    if previous.get("status") == "UNAVAILABLE":
        return False
    if current_provenance:
        for key in ("adapter", "adapter_version", "swissprot_path", "swissprot_version", "diamond", "diamond_version", "thresholds"):
            if key in previous and key in current_provenance and previous[key] != current_provenance[key]:
                return False
    return True


def _representation(data: dict[str, Any]) -> GenomeRepresentation:
    events = tuple(TransformEvent(**event) for event in data.get("transform_history", []))
    return GenomeRepresentation(data["original_sequence_id"], data["analysis_sequence_id"], "", "",
                                Topology(data["topology"]), Orientation(data["orientation"]), Rotation(data["rotation"]), events,
                                tuple(data.get("evidence", [])), data.get("reference"))


def _evidence(data: dict[str, Any]) -> Evidence:
    data = dict(data)
    data["level"] = EvidenceLevel(data["level"])
    allowed = {field.name for field in fields(Evidence)}
    return Evidence(**{key: value for key, value in data.items() if key in allowed})


def _load_source(source: Path) -> tuple[dict[str, Any], GenomeRepresentation, list[Protein], str, SequencingProvenance]:
    manifest_path = source / "run_manifest.json"
    if not manifest_path.is_file() and (source / "checkpoint_manifest.json").is_file():
        manifest_path = source / "checkpoint_manifest.json"
    manifest = json.loads(manifest_path.read_text()) if manifest_path.is_file() else {}
    sample_id = manifest.get("sample_id")
    def artifact(generic: str, suffix: str) -> Path:
        direct = source / generic
        if direct.is_file():
            return direct
        matches = sorted(source.glob(f"*_{suffix}"))
        return matches[0] if len(matches) == 1 else direct
    original_path = artifact("original_input.fasta", "original.fasta")
    analysis_path = artifact("analysis_genome.fasta", "analysis.fasta")
    proteins_path = artifact("proteins.faa", "proteins.faa")
    evidence_path = artifact("evidence.json", "evidence.json")
    genome_path = artifact("genome_representation.json", "genome.json")
    required = [manifest_path, original_path, analysis_path, proteins_path, evidence_path, genome_path]
    missing = [str(path) for path in required if not path.is_file()]
    if missing:
        raise ValueError("Resume validation failed; missing source artifacts: " + ", ".join(missing))
    input_id, input_sequence = read_fasta(original_path)
    if manifest.get("input_sha256") != checksum(original_path):
        raise ValueError("Resume validation failed; original input genome SHA-256 does not match run_manifest.json")
    analysis_id, analysis_sequence = read_fasta(analysis_path)
    rep_data = json.loads(genome_path.read_text())
    if rep_data.get("original_sequence_id") != input_id or rep_data.get("analysis_sequence_id") != analysis_id:
        raise ValueError("Resume validation failed; genome representation identifiers do not match FASTA artifacts")
    rep = _representation(rep_data)
    rep = GenomeRepresentation(rep.original_sequence_id, rep.analysis_sequence_id, input_sequence, analysis_sequence, rep.topology, rep.orientation, rep.rotation, rep.transform_history, rep.evidence, rep.reference)
    if rep.original_sequence != input_sequence or rep.analysis_sequence != analysis_sequence:
        raise ValueError("Resume validation failed; genome representation sequence content is inconsistent")
    raw = json.loads(evidence_path.read_text())
    if not isinstance(raw, list) or not raw:
        raise ValueError("Resume validation failed; evidence.json does not contain protein records")
    fasta_records = []
    header = None
    seq = []
    for line in proteins_path.read_text().splitlines():
        if line.startswith(">"):
            if header is not None: fasta_records.append((header, "".join(seq)))
            header, seq = line[1:].split()[0], []
        else: seq.append(line.strip())
    if header is not None: fasta_records.append((header, "".join(seq)))
    if len(fasta_records) != len(raw):
        raise ValueError("Resume validation failed; persisted protein FASTA count differs from evidence records")
    proteins: list[Protein] = []
    for item, (protein_id, sequence) in zip(raw, fasta_records):
        if item.get("protein_id") != protein_id or item.get("sequence") != sequence or len(item.get("sequence", "")) != len(sequence):
            raise ValueError(f"Resume validation failed; protein sequence/length mismatch for {protein_id}")
        item = dict(item)
        item["annotation_level"] = EvidenceLevel(item["annotation_level"])
        item["evidence"] = [_evidence(e) for e in item.get("evidence", [])]
        proteins.append(Protein(**{key: item[key] for key in Protein.__dataclass_fields__ if key in item}))
    if manifest.get("gene_caller", {}).get("name") != "PHANOTATE":
        raise ValueError("Resume validation failed; original gene-caller provenance is missing or not PHANOTATE")
    sequencing = SequencingProvenance.from_dict(json.loads((source / "sequencing_provenance.json").read_text())) if (source / "sequencing_provenance.json").is_file() else SequencingProvenance()
    return manifest, rep, proteins, input_sequence, sequencing


def recover_evidence_complete(source: str | Path, progress: ProgressReporter | None = None) -> int:
    """Regenerate only fusion and downstream outputs from validated evidence."""
    source = Path(source).resolve()
    status_path = source / "sample_status.json"
    if not status_path.is_file():
        raise ValueError("Evidence-complete recovery requires sample_status.json")
    status = json.loads(status_path.read_text())
    checkpoint = status.get("checkpoints", {}).get("evidence_integration")
    if not checkpoint or checkpoint.get("status") != "COMPLETE":
        raise ValueError("Evidence-complete recovery requires a COMPLETE evidence integration checkpoint")
    manifest, representation, proteins, _, sequencing = _load_source(source)
    checkpoint_source = source / "checkpoints" / "evidence_complete"
    if not checkpoint_source.is_dir():
        raise ValueError("Evidence-complete recovery failed; persisted checkpoint snapshot is missing")
    manifest, representation, proteins, _, sequencing = _load_source(checkpoint_source)
    if not manifest.get("evidence_adapters"):
        raise ValueError("Evidence-complete recovery failed; adapter provenance is missing")
    progress = progress or ProgressReporter(quiet=True)
    progress.STAGES = RESUME_STAGES
    for stage in REUSED:
        progress.start(f"REUSED: {stage}"); progress.finish("validated")
    for stage in ("PHROGs", "evidence integration"):
        progress.start(f"REUSED: {stage}"); progress.finish("validated")
    progress.start("functional classification")
    classifications = attach_gene_call_assessments(
        classify_proteins(proteins), _load_gene_call_assessments(source))
    context_records, modules = build_context(proteins, classifications)
    progress.finish("regenerated")
    progress.start("candidate ranking/mining")
    mine(proteins)
    candidates = ranked_candidates(proteins)
    progress.finish("regenerated")
    ranking_status = "INSUFFICIENT_EVIDENCE" if candidates and not any(e.supports and e.evidence_strength in {"STRONG", "EXPERIMENTAL"} for p in candidates for e in p.evidence) else "RANKED"
    progress.start("QC/report generation")
    quality = assess(representation.analysis_sequence, proteins)
    recovered = {**manifest, "command": "batch-recovery", "stage_status": {stage: "REUSED" for stage in REUSED + ("PHROGs", "evidence integration")} | {"functional classification": "RUN", "genomic context/modules": "RUN", "candidate ranking/mining": "RUN", "QC/report generation": "RUN", "GenBank package": "RUN"}, "recovery": {"mode": "evidence_complete", "source": str(source)}, "quality_control": quality, "discovery_ranking": {"status": ranking_status, "message": "Candidates ranked by available evidence."}}
    write_outputs(source, representation, sequencing, proteins, candidates, recovered, quality, source / "original_input.fasta", classifications, context_records, modules)
    progress.finish("regenerated")
    progress.start("GenBank package")
    write_package(source, representation.analysis_sequence_id, representation.analysis_sequence, proteins, recovered, sequencing_provenance=sequencing)
    progress.finish("regenerated")
    history = status.setdefault("audit_history", [])
    history.append({"status": status.get("status"), "failed_stage": status.get("failed_stage"), "error_message": status.get("error_message"), "recovered_at": datetime.now(timezone.utc).isoformat()})
    status.update({"status": "SUCCESS", "failed_stage": None, "current_stage": None, "recovery": "evidence_complete"})
    _atomic_json(status_path, status)
    (source / "run_manifest.json").write_text(json.dumps(recovered, indent=2, sort_keys=True))
    return len(proteins)


def resume(source: str | Path, output: str | Path, run_missing_evidence: bool = False,
           refresh_evidence: str | None = None, mmseqs: str | None = None,
           phrogs_path: str | Path | None = None, phrogs_annotations: str | Path | None = None,
           phrogs_hmm_path: str | Path | None = None,
           phrogs_evalue: float | None = 1e-5, phrogs_coverage: float | None = 0.5,
           phrogs_score: float | None = None, progress: ProgressReporter | None = None,
           threads: int = 1) -> int:
    source, output = Path(source).resolve(), Path(output).resolve()
    if source == output or output.exists():
        raise ValueError("Resume output must be a new, non-existent directory")
    manifest, representation, proteins, input_sequence, sequencing = _load_source(source)
    adapters = manifest.get("evidence_adapters", [])
    old_by_name = {item.get("adapter"): item for item in adapters}
    should_run = run_missing_evidence or refresh_evidence == "PHROGS"
    old_phrogs = old_by_name.get("PHROGSMMseqsAdapter", {})
    old_hmm = old_by_name.get("PHROGSPyHMMERAdapter", {})
    if not should_run and old_phrogs.get("status") != "UNAVAILABLE" and old_hmm.get("status") != "UNAVAILABLE":
        raise ValueError("No missing PHROGs adapter to run; use --refresh-evidence PHROGS to force refresh")
    manager = EvidenceResourceManager()
    registered = manager.find(ResourceType.PHROGS)
    explicit_hmm = phrogs_hmm_path is not None
    if phrogs_path is None and registered:
        provenance = registered.get("provenance", {})
        phrogs_path = registered["path"]
        phrogs_annotations = phrogs_annotations or provenance.get("annotations_path")
        phrogs_hmm_path = phrogs_hmm_path or provenance.get("hmm_profiles_path")
    if phrogs_path is None:
        raise ValueError("Resume validation failed; PHROGs database is unavailable")
    phrogs_version = registered.get("version") if registered else None
    adapter = PHROGSMMseqsAdapter(
        phrogs_path, phrogs_annotations, mmseqs, phrogs_version,
        phrogs_evalue, phrogs_coverage, phrogs_score, threads=threads)
    hmm_adapter = PHROGSPyHMMERAdapter(
        phrogs_hmm_path, phrogs_annotations, phrogs_version,
        phrogs_evalue, phrogs_coverage, phrogs_score, threads=threads)
    progress = progress or ProgressReporter(quiet=True)
    # Resume has its own eleven-stage ledger; normal run progress is unchanged.
    progress.STAGES = RESUME_STAGES
    for stage in REUSED:
        progress.start(f"REUSED: {stage}"); progress.finish("validated")
    progress.start("PHROGs")
    phrogs_result = adapter.analyze(proteins)
    hmm_result = hmm_adapter.analyze(proteins)
    if phrogs_result.status != "REAL" and hmm_result.status != "REAL":
        raise ValueError(
            f"PHROGs resume failed: MMseqs2: {phrogs_result.message}; PyHMMER: {hmm_result.message}")
    if explicit_hmm and hmm_result.status != "REAL":
        raise ValueError(f"Explicit PHROGs HMM refresh failed: {hmm_result.message}")
    merged_phrogs = merge_phrogs_evidence(phrogs_result.evidence, hmm_result.evidence)
    by_id = {protein.protein_id: protein for protein in proteins}
    # A refresh replaces prior PHROGs calls. Retaining them would silently
    # duplicate evidence and could overweight the same biological source.
    for protein in proteins:
        protein.evidence = [evidence for evidence in protein.evidence if evidence.source != "PHROGs"]
    for evidence in merged_phrogs:
        protein_id = evidence.provenance.get("protein_id") or evidence.metrics.get("query_protein_id")
        if protein_id in by_id: by_id[protein_id].evidence.append(evidence)
    progress.finish(
        f"{sum(e.supports for e in merged_phrogs)} accepted deduplicated hits; "
        f"MMseqs2={phrogs_result.state.value}; PyHMMER={hmm_result.state.value}")
    stages = ["evidence integration", "candidate ranking/mining", "QC/report generation", "GenBank package"]
    progress.start("evidence integration")
    progress.finish("integrated")
    classifications = attach_gene_call_assessments(
        classify_proteins(proteins), _load_gene_call_assessments(source))
    context_records, modules = build_context(proteins, classifications)
    progress.start("candidate ranking/mining")
    mine(proteins)
    candidates = ranked_candidates(proteins)
    progress.finish("ranked")
    progress.start("QC/report generation")
    ranking_status = "INSUFFICIENT_EVIDENCE" if candidates and not any(e.supports and e.evidence_strength in {"STRONG", "EXPERIMENTAL"} for p in candidates for e in p.evidence) else "RANKED"
    quality = assess(representation.analysis_sequence, proteins)
    progress.finish("regenerated")
    stage_status = {stage: "REUSED" for stage in REUSED} | {"PHROGs": "RUN", **{stage: "RUN" for stage in stages}}
    new_adapters = [item for item in adapters if item.get("adapter") not in {
        "PHROGSMMseqsAdapter", "PHROGSPyHMMERAdapter"}]
    new_adapters.append({"adapter": phrogs_result.adapter, "status": phrogs_result.status, "provenance": phrogs_result.provenance, "message": phrogs_result.message})
    new_adapters.append({"adapter": hmm_result.adapter, "status": hmm_result.status, "provenance": hmm_result.provenance, "message": hmm_result.message})
    new_manifest = {**manifest, "command": "resume", "input": str(source / "original_input.fasta"), "input_sha256": checksum(source / "original_input.fasta"), "evidence_adapters": new_adapters, "stage_status": stage_status, "resume": {"source": str(source), "reused_stages": [*REUSED], "run_stages": ["PHROGs", *stages]}, "discovery_ranking": {"status": ranking_status, "message": "Candidates ranked by available evidence."}}
    new_manifest["quality_control"] = quality
    temp = Path(tempfile.mkdtemp(prefix=f".{output.name}.", dir=str(output.parent)))
    try:
        write_outputs(temp, representation, sequencing, proteins, candidates, new_manifest, quality, source / "original_input.fasta", classifications, context_records, modules)
        _copy_gene_call_artifacts(source, temp)
        progress.start("GenBank package")
        write_package(temp, representation.analysis_sequence_id, representation.analysis_sequence, proteins, new_manifest, sequencing_provenance=sequencing)
        progress.finish("regenerated")
        os.replace(temp, output)
    except Exception:
        shutil.rmtree(temp, ignore_errors=True)
        raise
    return len(proteins)
