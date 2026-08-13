from __future__ import annotations

import csv
import html
import json
from dataclasses import asdict
from datetime import datetime, timezone
from pathlib import Path

from .models import Protein
from .genome_representation import GenomeRepresentation
from .sequencing_provenance import SequencingProvenance


def write_outputs(output: str | Path, representation: GenomeRepresentation, sequencing_provenance: SequencingProvenance, proteins: list[Protein], candidates: list[Protein], manifest: dict, quality_control: dict | None = None, original_fasta: str | Path | None = None) -> None:
    root = Path(output)
    root.mkdir(parents=True, exist_ok=True)
    if original_fasta is not None:
        (root / "original_input.fasta").write_bytes(Path(original_fasta).read_bytes())
    (root / "analysis_genome.fasta").write_text(f">{representation.analysis_sequence_id}\n{representation.analysis_sequence}\n")
    (root / "genome_representation.json").write_text(json.dumps(representation.manifest(), indent=2, sort_keys=True))
    (root / "sequencing_provenance.json").write_text(json.dumps(sequencing_provenance.manifest(), indent=2, sort_keys=True))
    (root / "proteins.faa").write_text("".join(f">{p.protein_id} genome={p.genome_id} start={p.start} end={p.end}\n{p.sequence}\n" for p in proteins))
    (root / "cds.fna").write_text("".join(f">{p.protein_id}\n{p.cds}\n" for p in proteins))
    (root / "genes.gff3").write_text("##gff-version 3\n" + "".join(f"{representation.analysis_sequence_id}\tPhageMine\tCDS\t{p.start}\t{p.end}\t.\t{p.strand}\t0\tID={p.protein_id};Name={p.protein_id};calling_source={p.gene_call_source};coordinate_representation={representation.analysis_sequence_id};coordinate_system=1-based-inclusive;raw_start={p.gene_call_parameters.get("raw_start", p.start)};raw_end={p.gene_call_parameters.get("raw_end", p.end)}\n" for p in proteins))
    columns = ["analysis_sequence_id", "protein_id", "start", "end", "strand", "length", "annotation", "annotation_level", "functional_confidence", "biological_interest", "evidence_diversity"]
    with (root / "annotation.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t")
        writer.writeheader()
        for p in proteins:
            writer.writerow({key: p.genome_id if key == "analysis_sequence_id" else (p.annotation_level.value if key == "annotation_level" else getattr(p, key)) for key in columns})
    with (root / "candidate_ranking.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["rank", "protein_id", "annotation", "biological_interest", "functional_confidence", "evidence_diversity", "score_components"])
        insufficient = manifest.get("discovery_ranking", {}).get("status") == "INSUFFICIENT_EVIDENCE"
        for rank, p in enumerate(candidates, 1):
            writer.writerow(["NA" if insufficient else rank, p.protein_id, p.annotation, p.biological_interest, p.functional_confidence, p.evidence_diversity, json.dumps(p.score_components, sort_keys=True)])
    (root / "evidence.json").write_text(json.dumps([asdict(p) for p in proteins], indent=2, default=str))
    manifest["created_at"] = datetime.now(timezone.utc).isoformat()
    (root / "run_manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True))
    if quality_control is not None:
        (root / "quality_control.json").write_text(json.dumps(quality_control, indent=2, sort_keys=True))
    markdown = _markdown(representation, sequencing_provenance, proteins, candidates, manifest)
    (root / "report.md").write_text(markdown)
    (root / "report.html").write_text("<html><body><pre>" + html.escape(markdown) + "</pre></body></html>")


def _markdown(representation: GenomeRepresentation, sequencing_provenance: SequencingProvenance, proteins: list[Protein], candidates: list[Protein], manifest: dict) -> str:
    transforms = [f"- `{event.operation}`: {event.rationale}" for event in representation.transform_history] or ["- No transformations; analysis uses the original input representation."]
    statuses = [item.get("status") for item in manifest.get("evidence_adapters", [])]
    warning = "> **Evidence status:** MOCK evidence is demonstration-only and not biological evidence." if "MOCK" in statuses else "> **Evidence status:** Functional evidence is limited to explicitly available adapters; unavailable evidence is not fabricated."
    ranking = manifest.get("discovery_ranking", {})
    ranking_status = ranking.get("status", "RANKED")
    ranking_message = ranking.get("message", "Candidates ranked by available evidence.")
    lines = [f"# PhageMine report: {representation.analysis_sequence_id}", "", warning, "", "## Genome representation", "", f"- Original input sequence: `{representation.original_sequence_id}`", f"- Analysis sequence: `{representation.analysis_sequence_id}`", f"- Topology: `{representation.topology.value}`", f"- Orientation: `{representation.orientation.value}`", f"- Rotation: `{representation.rotation.value}`", "- All reported gene, protein, neighborhood, and annotation coordinates are relative to the analysis sequence.", "", "### Explicit transformation history", *transforms, "", "## Sequencing provenance", "", f"- Sequencing platform: `{sequencing_provenance.sequencing_platform.value}`", f"- Assembler: `{sequencing_provenance.assembler or 'UNKNOWN'}`", f"- Polishing method: `{sequencing_provenance.polishing_method or 'UNKNOWN'}`", "- Sequencing provenance does not determine genome topology, orientation, or rotation.", "", f"Predicted proteins: **{len(proteins)}**  ", f"Poorly characterised proteins: **{len(candidates)}**", "", "## Discovery ranking", "", f"- Status: **{ranking_status}**", f"- {ranking_message}", "", "## Top candidates worth investigating", "", "| Rank | Protein | Current annotation | Biological interest | Functional confidence | Evidence diversity |", "|---:|---|---|---:|---|---|"]
    for rank, p in enumerate(candidates, 1):
        display_rank = "NA" if ranking_status == "INSUFFICIENT_EVIDENCE" else rank
        lines.append(f"| {display_rank} | {p.protein_id} | {p.annotation} | {p.biological_interest} | {p.functional_confidence} | {p.evidence_diversity} |")
    for p in candidates:
        supporting = [e.statement for e in p.evidence if e.supports]
        contradicting = [e.statement for e in p.evidence if not e.supports]
        lines += ["", f"## Candidate: {p.protein_id}", "", f"**Current annotation:** {p.annotation}  ", f"**Biological interest:** {p.biological_interest}/100  ", f"**Functional confidence:** {p.functional_confidence}  ", f"**Evidence diversity:** {p.evidence_diversity}", "", "### Why prioritised", *[f"- {item}" for item in supporting], "", "### Contradicting evidence", *([f"- {item}" for item in contradicting] or ["- No explicit contradictory evidence was generated; absence is not confirmation."]), "", "### Missing evidence", *[f"- {item}" for item in p.missing_evidence], "", "### Alternative hypotheses", *[f"- {item}" for item in p.alternatives], "", "**Confidence statement:** Computational hypothesis only. Experimental validation is required before assigning biological function."]
    return "\n".join(lines) + "\n"
