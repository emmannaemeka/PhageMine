"""Fail-closed continuation of a completed PhageMine run."""
from __future__ import annotations

import hashlib
import json
import os
import shutil
import tempfile
from dataclasses import fields
from pathlib import Path
from typing import Any

from .genome_representation import GenomeRepresentation, Orientation, Rotation, Topology, TransformEvent
from .genbank import write_package
from .io import checksum, read_fasta
from .mining import mine, ranked_candidates
from .models import Evidence, EvidenceLevel, Protein
from .phrogs import PHROGSMMseqsAdapter
from .progress import ProgressReporter
from .quality import assess
from .reporting import write_outputs
from .resources import EvidenceResourceManager, ResourceType
from .sequencing_provenance import SequencingProvenance
from .fusion import classify_proteins


REUSED = ("input/genome validation", "genome representation", "gene prediction", "Pfam", "VOGDB", "Swiss-Prot")
RESUME_STAGES = REUSED + ("PHROGs", "evidence integration", "candidate ranking/mining", "QC/report generation", "GenBank package")


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
    required = [manifest_path, source / "original_input.fasta", source / "analysis_genome.fasta", source / "proteins.faa", source / "evidence.json", source / "genome_representation.json"]
    missing = [str(path) for path in required if not path.is_file()]
    if missing:
        raise ValueError("Resume validation failed; missing source artifacts: " + ", ".join(missing))
    manifest = json.loads(manifest_path.read_text())
    input_id, input_sequence = read_fasta(source / "original_input.fasta")
    if manifest.get("input_sha256") != checksum(source / "original_input.fasta"):
        raise ValueError("Resume validation failed; original input genome SHA-256 does not match run_manifest.json")
    analysis_id, analysis_sequence = read_fasta(source / "analysis_genome.fasta")
    rep_data = json.loads((source / "genome_representation.json").read_text())
    if rep_data.get("original_sequence_id") != input_id or rep_data.get("analysis_sequence_id") != analysis_id:
        raise ValueError("Resume validation failed; genome representation identifiers do not match FASTA artifacts")
    rep = _representation(rep_data)
    rep = GenomeRepresentation(rep.original_sequence_id, rep.analysis_sequence_id, input_sequence, analysis_sequence, rep.topology, rep.orientation, rep.rotation, rep.transform_history, rep.evidence, rep.reference)
    if rep.original_sequence != input_sequence or rep.analysis_sequence != analysis_sequence:
        raise ValueError("Resume validation failed; genome representation sequence content is inconsistent")
    raw = json.loads((source / "evidence.json").read_text())
    if not isinstance(raw, list) or len(raw) != 99:
        raise ValueError("Resume validation failed; evidence.json does not contain the expected 99 protein records")
    fasta_records = []
    header = None
    seq = []
    for line in (source / "proteins.faa").read_text().splitlines():
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


def resume(source: str | Path, output: str | Path, run_missing_evidence: bool = False,
           refresh_evidence: str | None = None, mmseqs: str | None = None,
           phrogs_path: str | Path | None = None, phrogs_annotations: str | Path | None = None,
           phrogs_evalue: float | None = 1e-5, phrogs_coverage: float | None = 0.5,
           phrogs_score: float | None = None, progress: ProgressReporter | None = None) -> int:
    source, output = Path(source).resolve(), Path(output).resolve()
    if source == output or output.exists():
        raise ValueError("Resume output must be a new, non-existent directory")
    manifest, representation, proteins, input_sequence, sequencing = _load_source(source)
    adapters = manifest.get("evidence_adapters", [])
    old_by_name = {item.get("adapter"): item for item in adapters}
    should_run = run_missing_evidence or refresh_evidence == "PHROGS"
    old_phrogs = old_by_name.get("PHROGSMMseqsAdapter", {})
    if not should_run and old_phrogs.get("status") != "UNAVAILABLE":
        raise ValueError("No missing PHROGs adapter to run; use --refresh-evidence PHROGS to force refresh")
    manager = EvidenceResourceManager()
    registered = manager.find(ResourceType.PHROGS)
    if phrogs_path is None and registered:
        phrogs_path, phrogs_annotations = registered["path"], phrogs_annotations or registered.get("provenance", {}).get("annotations_path")
    if phrogs_path is None:
        raise ValueError("Resume validation failed; PHROGs database is unavailable")
    phrogs_version = registered.get("version") if registered else None
    adapter = PHROGSMMseqsAdapter(phrogs_path, phrogs_annotations, mmseqs, phrogs_version, phrogs_evalue, phrogs_coverage, phrogs_score)
    if not adapter.available():
        raise ValueError("Resume validation failed; PHROGs MMseqs2 database or executable is unavailable")
    progress = progress or ProgressReporter(quiet=True)
    # Resume has its own eleven-stage ledger; normal run progress is unchanged.
    progress.STAGES = RESUME_STAGES
    for stage in REUSED:
        progress.start(f"REUSED: {stage}"); progress.finish("validated")
    progress.start("PHROGs")
    phrogs_result = adapter.analyze(proteins)
    if phrogs_result.status != "REAL":
        raise ValueError(f"PHROGs resume failed: {phrogs_result.message}")
    by_id = {protein.protein_id: protein for protein in proteins}
    for evidence in phrogs_result.evidence:
        protein_id = evidence.provenance.get("protein_id") or evidence.metrics.get("query_protein_id")
        if protein_id in by_id: by_id[protein_id].evidence.append(evidence)
    progress.finish(f"{sum(e.supports for e in phrogs_result.evidence)} accepted hits")
    stages = ["evidence integration", "candidate ranking/mining", "QC/report generation", "GenBank package"]
    progress.start("evidence integration")
    progress.finish("integrated")
    classifications = classify_proteins(proteins)
    progress.start("candidate ranking/mining")
    mine(proteins)
    candidates = ranked_candidates(proteins)
    progress.finish("ranked")
    progress.start("QC/report generation")
    ranking_status = "INSUFFICIENT_EVIDENCE" if candidates and not any(e.supports and e.evidence_strength in {"STRONG", "EXPERIMENTAL"} for p in candidates for e in p.evidence) else "RANKED"
    quality = assess(representation.analysis_sequence, proteins)
    progress.finish("regenerated")
    stage_status = {stage: "REUSED" for stage in REUSED} | {"PHROGs": "RUN", **{stage: "RUN" for stage in stages}}
    new_adapters = [item for item in adapters if item.get("adapter") != "PHROGSMMseqsAdapter"]
    new_adapters.append({"adapter": phrogs_result.adapter, "status": phrogs_result.status, "provenance": phrogs_result.provenance, "message": phrogs_result.message})
    new_manifest = {**manifest, "command": "resume", "input": str(source / "original_input.fasta"), "input_sha256": checksum(source / "original_input.fasta"), "evidence_adapters": new_adapters, "stage_status": stage_status, "resume": {"source": str(source), "reused_stages": [*REUSED], "run_stages": ["PHROGs", *stages]}, "discovery_ranking": {"status": ranking_status, "message": "Candidates ranked by available evidence."}}
    new_manifest["quality_control"] = quality
    temp = Path(tempfile.mkdtemp(prefix=f".{output.name}.", dir=str(output.parent)))
    try:
        write_outputs(temp, representation, sequencing, proteins, candidates, new_manifest, quality, source / "original_input.fasta", classifications)
        progress.start("GenBank package")
        write_package(temp, representation.analysis_sequence_id, representation.analysis_sequence, proteins, new_manifest, sequencing_provenance=sequencing)
        progress.finish("regenerated")
        os.replace(temp, output)
    except Exception:
        shutil.rmtree(temp, ignore_errors=True)
        raise
    return len(proteins)
