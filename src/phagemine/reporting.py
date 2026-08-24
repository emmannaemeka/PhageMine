from __future__ import annotations

import csv
import html
import json
import os
import shutil
import hashlib
from dataclasses import asdict
from datetime import datetime, timezone
from pathlib import Path

from .models import Protein
from .genome_representation import GenomeRepresentation
from .sequencing_provenance import SequencingProvenance
from .fusion import write_classification
from .context import write_context
from .figures import generate_annotation_figures


def _evidence_lines(protein: Protein) -> list[str]:
    by_source = {"Pfam": [], "VOGDB": [], "PHROGs": [], "Swiss-Prot": []}
    for evidence in protein.evidence:
        source = evidence.source
        if source in by_source:
            label = evidence.description or evidence.identifier or evidence.family_name or "no informative hit"
            by_source[source].append((label, bool(evidence.supports)))
    return [f"{source}: " + ("; ".join(label + (" [support]" if ok else " [not informative]") for label, ok in values) if values else "no informative hit") for source, values in by_source.items()]


def _write_protein_details(root: Path, proteins: list[Protein], classifications: list[dict], contexts: list[dict] | None) -> list[str]:
    """Publish deterministic, directly navigable records without changing scientific data."""
    detail_dir = root / "protein_details"; fasta_dir = detail_dir / "fasta"
    detail_dir.mkdir(parents=True, exist_ok=True); fasta_dir.mkdir(exist_ok=True)
    cls = {r["protein_id"]: r for r in classifications}; ctx = {r.get("protein_id"): r for r in (contexts or [])}
    links = []
    for protein in proteins:
        record = cls.get(protein.protein_id, {})
        aa = f">{protein.protein_id} genome={protein.genome_id}\n{protein.sequence}\n"
        cds = f">{protein.protein_id}\n{protein.cds}\n"
        (fasta_dir / f"{protein.protein_id}.faa").write_text(aa)
        (fasta_dir / f"{protein.protein_id}.fna").write_text(cds)
        evidence = [{"source": e.source, "identifier": e.identifier, "description": e.description, "supports": e.supports, "strength": e.evidence_strength, "metrics": e.metrics} for e in protein.evidence]
        payload = {"protein_id": protein.protein_id, "genome_id": protein.genome_id, "start": protein.start, "end": protein.end, "strand": protein.strand, "length_aa": protein.length, "classification": record, "evidence": evidence, "genomic_context": ctx.get(protein.protein_id)}
        (detail_dir / f"{protein.protein_id}.json").write_text(json.dumps(payload, indent=2, sort_keys=True, default=str))
        reason = record.get("reasoning_summary") or "No deterministic evidence-based explanation was generated."
        evidence_html = "".join(f"<li>{html.escape(line)}</li>" for line in _evidence_lines(protein))
        context_text = json.dumps(ctx.get(protein.protein_id), sort_keys=True) if ctx.get(protein.protein_id) else "No genomic-context record available."
        page = f"<!doctype html><html><head><meta charset='utf-8'><title>{html.escape(protein.protein_id)}</title></head><body><h1>{html.escape(protein.protein_id)}</h1><p><b>Coordinates:</b> {protein.start}-{protein.end} &nbsp; <b>Strand:</b> {html.escape(protein.strand)} &nbsp; <b>Length:</b> {protein.length} aa</p><h2>Classification</h2><p>{html.escape(str(record.get('functional_state') or 'UNRESOLVED'))}</p><p><b>Proposed product:</b> {html.escape(record.get('display_product') or record.get('proposed_function') or 'hypothetical protein')}</p><p><b>Confidence:</b> {html.escape(str(record.get('confidence') or 'NONE'))}</p><p><b>Reason for annotation:</b> {html.escape(reason)}</p><h2>Evidence</h2><ul>{evidence_html}</ul><h2>Genomic context</h2><pre>{html.escape(context_text)}</pre><p><a download href='fasta/{html.escape(protein.protein_id)}.faa'>Protein FASTA</a> | <a download href='fasta/{html.escape(protein.protein_id)}.fna'>CDS FASTA</a> | <a href='{html.escape(protein.protein_id)}.json'>Evidence JSON</a></p></body></html>"
        (detail_dir / f"{protein.protein_id}.html").write_text(page)
        links.append(str(detail_dir / f"{protein.protein_id}.html"))
    return links


def extract_protein_record(run: str | Path, protein_id: str) -> dict:
    root = Path(run)
    candidates = [root] + [p for p in root.iterdir() if p.is_dir()] if root.is_dir() else []
    for sample in candidates:
        evidence_path = sample / "evidence.json"
        if not evidence_path.is_file(): continue
        records = json.loads(evidence_path.read_text())
        for record in records:
            if record.get("protein_id") == protein_id:
                return {"protein_id": protein_id, "protein_fasta": f">{protein_id}\n{record.get('sequence','')}\n", "cds_fasta": f">{protein_id}\n{record.get('cds','')}\n", "evidence": record, "detail_html": str(sample / "protein_details" / f"{protein_id}.html")}
    raise KeyError(f"protein ID not found in run: {protein_id}")


def write_checkpoint_snapshot(root: str | Path, checkpoint_dir: str | Path, representation: GenomeRepresentation,
                              sequencing_provenance: SequencingProvenance, proteins: list[Protein],
                              manifest: dict, original_fasta: str | Path) -> None:
    """Atomically persist the minimum scientific state needed for recovery."""
    destination = Path(checkpoint_dir); destination.mkdir(parents=True, exist_ok=True)
    files = {
        "original_input.fasta": Path(original_fasta).read_bytes(),
        "analysis_genome.fasta": f">{representation.analysis_sequence_id}\n{representation.analysis_sequence}\n".encode(),
        "genome_representation.json": json.dumps(representation.manifest(), indent=2, sort_keys=True).encode(),
        "proteins.faa": "".join(f">{p.protein_id} genome={p.genome_id} start={p.start} end={p.end}\n{p.sequence}\n" for p in proteins).encode(),
        "cds.fna": "".join(f">{p.protein_id}\n{p.cds}\n" for p in proteins).encode(),
        "evidence.json": json.dumps([asdict(p) for p in proteins], indent=2, default=str).encode(),
        "checkpoint_manifest.json": json.dumps(manifest, indent=2, sort_keys=True).encode(),
    }
    for name, content in files.items():
        temporary = destination / (name + ".tmp")
        temporary.write_bytes(content)
        os.replace(temporary, destination / name)


def write_stage_checkpoint(checkpoint_dir: str | Path, stage: str, payload: dict, provenance: dict | None = None) -> None:
    """Atomically persist a stage payload and its provenance."""
    destination = Path(checkpoint_dir) / stage
    destination.mkdir(parents=True, exist_ok=True)
    encoded = json.dumps(payload, indent=2, sort_keys=True, default=str).encode()
    manifest = {"stage": stage, "provenance": provenance or {}, "schema_version": "1.1",
                "output_checksums": {"evidence.json": hashlib.sha256(encoded).hexdigest()},
                "dependencies": ["GENE_PREDICTION"]}
    for name, value in (("evidence.json", payload), ("checkpoint_manifest.json", manifest)):
        temporary = destination / (name + ".tmp")
        temporary.write_text(json.dumps(value, indent=2, sort_keys=True, default=str))
        os.replace(temporary, destination / name)


def prefix_checkpoint_artifacts(sample_root: str | Path, sample_id: str) -> None:
    """Publish sample-identifying aliases while retaining legacy artifact names."""
    root = Path(sample_root)
    mappings = {
        "checkpoints/gene_prediction/original_input.fasta": f"checkpoints/gene_prediction/{sample_id}_original.fasta",
        "checkpoints/gene_prediction/analysis_genome.fasta": f"checkpoints/gene_prediction/{sample_id}_analysis.fasta",
        "checkpoints/gene_prediction/proteins.faa": f"checkpoints/gene_prediction/{sample_id}_proteins.faa",
        "checkpoints/gene_prediction/cds.fna": f"checkpoints/gene_prediction/{sample_id}_cds.fna",
        "checkpoints/gene_prediction/genome_representation.json": f"checkpoints/gene_prediction/{sample_id}_genome.json",
        "checkpoints/evidence_complete/evidence.json": f"checkpoints/evidence_complete/{sample_id}_evidence.json",
        "checkpoints/evidence_complete/checkpoint_manifest.json": f"checkpoints/evidence_complete/{sample_id}_checkpoint.json",
        "checkpoints/pfam/evidence.json": f"checkpoints/pfam/{sample_id}_pfam_evidence.json",
        "checkpoints/vogdb/evidence.json": f"checkpoints/vogdb/{sample_id}_vogdb_evidence.json",
        "checkpoints/swissprot/evidence.json": f"checkpoints/swissprot/{sample_id}_swissprot_evidence.json",
        "checkpoints/phrogs/evidence.json": f"checkpoints/phrogs/{sample_id}_phrogs_evidence.json",
    }
    for source, target in mappings.items():
        source_path = root / source
        if source_path.is_file():
            target_path = root / target
            temporary = target_path.with_name(target_path.name + ".tmp")
            if target.endswith(".json"):
                try:
                    value = json.loads(source_path.read_text())
                    if isinstance(value, list):
                        value = [{**item, "sample_id": sample_id} if isinstance(item, dict) else item for item in value]
                    elif isinstance(value, dict):
                        value = {**value, "sample_id": sample_id}
                    temporary.write_text(json.dumps(value, indent=2, sort_keys=True, default=str))
                except (OSError, ValueError, TypeError):
                    shutil.copyfile(source_path, temporary)
            else:
                shutil.copyfile(source_path, temporary)
            os.replace(temporary, target_path)


def write_outputs(output: str | Path, representation: GenomeRepresentation, sequencing_provenance: SequencingProvenance, proteins: list[Protein], candidates: list[Protein], manifest: dict, quality_control: dict | None = None, original_fasta: str | Path | None = None, classifications: list[dict] | None = None, context_records: list[dict] | None = None, modules: list[dict] | None = None) -> None:
    root = Path(output)
    root.mkdir(parents=True, exist_ok=True)
    if original_fasta is not None:
        (root / "original_input.fasta").write_bytes(Path(original_fasta).read_bytes())
    (root / "analysis_genome.fasta").write_text(f">{representation.analysis_sequence_id}\n{representation.analysis_sequence}\n")
    (root / "genome_representation.json").write_text(json.dumps(representation.manifest(), indent=2, sort_keys=True))
    (root / "sequencing_provenance.json").write_text(json.dumps(sequencing_provenance.manifest(), indent=2, sort_keys=True))
    (root / "proteins.faa").write_text("".join(f">{p.protein_id} genome={p.genome_id} start={p.start} end={p.end}\n{p.sequence}\n" for p in proteins))
    (root / "cds.fna").write_text("".join(f">{p.protein_id}\n{p.cds}\n" for p in proteins))
    gff_records = []
    for p in proteins:
        raw_start = p.gene_call_parameters.get("raw_start", p.start)
        raw_end = p.gene_call_parameters.get("raw_end", p.end)
        gff_records.append(
            f"{representation.analysis_sequence_id}\tPhageMine\tCDS\t{p.start}\t{p.end}\t.\t{p.strand}\t0\t"
            f"ID={p.protein_id};Name={p.protein_id};calling_source={p.gene_call_source};"
            f"coordinate_representation={representation.analysis_sequence_id};coordinate_system=1-based-inclusive;"
            f"raw_start={raw_start};raw_end={raw_end}\n"
        )
    (root / "genes.gff3").write_text("##gff-version 3\n" + "".join(gff_records))
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
    write_classification(root, proteins, classifications)
    if context_records is not None and modules is not None:
        write_context(root, context_records, modules)
    if quality_control is not None:
        (root / "quality_control.json").write_text(json.dumps(quality_control, indent=2, sort_keys=True))
    manifest["figures"] = generate_annotation_figures(root, proteins, classifications or [], context_records, modules)
    detail_links = _write_protein_details(root, proteins, classifications or [], context_records)
    manifest["protein_detail_records"] = detail_links
    manifest["created_at"] = datetime.now(timezone.utc).isoformat()
    (root / "run_manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True))
    markdown = _markdown(representation, sequencing_provenance, proteins, candidates, manifest)
    (root / "report.md").write_text(markdown)
    rows = []
    cls_by_id = {r.get("protein_id"): r for r in (classifications or [])}
    for protein in proteins:
        c = cls_by_id.get(protein.protein_id, {})
        rows.append(f"<tr><td><a href='protein_details/{html.escape(protein.protein_id)}.html'>{html.escape(protein.protein_id)}</a></td><td>{protein.start}-{protein.end}</td><td>{html.escape(protein.strand)}</td><td>{html.escape(str(c.get('functional_state') or 'UNRESOLVED'))}</td><td>{html.escape(str(c.get('display_product') or c.get('proposed_function') or 'hypothetical protein'))}</td><td>{html.escape(str(c.get('confidence') or 'NONE'))}</td></tr>")
    report_html = f"<html><head><meta charset='utf-8'><style>body{{font-family:Arial;max-width:1400px;margin:auto;padding:2em}}table{{border-collapse:collapse;width:100%}}td,th{{border:1px solid #bbb;padding:.3em}}th{{background:#eee}}</style></head><body><h1>PhageMine annotation report</h1><h2>Summary</h2><p>Predicted proteins: {len(proteins)}</p><h2>Protein annotation table</h2><table><tr><th>Protein</th><th>Coordinates</th><th>Strand</th><th>Classification</th><th>Proposed product</th><th>Confidence</th></tr>{''.join(rows)}</table><h2>Methods and provenance</h2><pre>{html.escape(markdown)}</pre></body></html>"
    (root / "report.html").write_text(report_html)


def _markdown(representation: GenomeRepresentation, sequencing_provenance: SequencingProvenance, proteins: list[Protein], candidates: list[Protein], manifest: dict) -> str:
    transforms = [f"- `{event.operation}`: {event.rationale}" for event in representation.transform_history] or ["- No transformations; analysis uses the original input representation."]
    statuses = [item.get("status") for item in manifest.get("evidence_adapters", [])]
    warning = "> **Evidence status:** MOCK evidence is demonstration-only and not biological evidence." if "MOCK" in statuses else "> **Evidence status:** Functional evidence is limited to explicitly available adapters; unavailable evidence is not fabricated."
    ranking = manifest.get("discovery_ranking", {})
    ranking_status = ranking.get("status", "RANKED")
    ranking_message = ranking.get("message", "Candidates ranked by available evidence.")
    lines = [f"# PhageMine report: {representation.analysis_sequence_id}", "", warning, "", "## Genome representation", "", f"- Original input sequence: `{representation.original_sequence_id}`", f"- Analysis sequence: `{representation.analysis_sequence_id}`", f"- Topology: `{representation.topology.value}`", f"- Orientation: `{representation.orientation.value}`", f"- Rotation: `{representation.rotation.value}`", "- All reported gene, protein, neighborhood, and annotation coordinates are relative to the analysis sequence.", "", "### Explicit transformation history", *transforms, "", "## Sequencing provenance", "", f"- Sequencing platform: `{sequencing_provenance.sequencing_platform.value}`", f"- Assembler: `{sequencing_provenance.assembler or 'UNKNOWN'}`", f"- Polishing method: `{sequencing_provenance.polishing_method or 'UNKNOWN'}`", "- Sequencing provenance does not determine genome topology, orientation, or rotation.", "", f"Predicted proteins: **{len(proteins)}**  ", f"Poorly characterised proteins: **{len(candidates)}**", "", "## Figures", ""]
    for figure in manifest.get("figures", {}).get("created", []):
        try: lines.append(f"- [{Path(figure).name}]({Path(figure).relative_to(Path(manifest.get('output', '.')) if manifest.get('output') else Path(figure).parent.parent)})")
        except ValueError: lines.append(f"- `{figure}`")
    for item in manifest.get("figures", {}).get("skipped", []): lines.append(f"- Skipped `{item['figure']}`: {item['reason']}")
    lines += ["", "## Discovery ranking", "", f"- Status: **{ranking_status}**", f"- {ranking_message}", "", "## Top candidates worth investigating", "", "| Rank | Protein | Current annotation | Biological interest | Functional confidence | Evidence diversity |", "|---:|---|---|---:|---|---|"]
    for rank, p in enumerate(candidates, 1):
        display_rank = "NA" if ranking_status == "INSUFFICIENT_EVIDENCE" else rank
        lines.append(f"| {display_rank} | {p.protein_id} | {p.annotation} | {p.biological_interest} | {p.functional_confidence} | {p.evidence_diversity} |")
    for p in candidates:
        supporting = [e.statement for e in p.evidence if e.supports]
        contradicting = [e.statement for e in p.evidence if not e.supports]
        lines += ["", f"## Candidate: {p.protein_id}", "", f"**Current annotation:** {p.annotation}  ", f"**Biological interest:** {p.biological_interest}/100  ", f"**Functional confidence:** {p.functional_confidence}  ", f"**Evidence diversity:** {p.evidence_diversity}", "", "### Why prioritised", *[f"- {item}" for item in supporting], "", "### Contradicting evidence", *([f"- {item}" for item in contradicting] or ["- No explicit contradictory evidence was generated; absence is not confirmation."]), "", "### Missing evidence", *[f"- {item}" for item in p.missing_evidence], "", "### Alternative hypotheses", *[f"- {item}" for item in p.alternatives], "", "**Confidence statement:** Computational hypothesis only. Experimental validation is required before assigning biological function."]
    return "\n".join(lines) + "\n"
