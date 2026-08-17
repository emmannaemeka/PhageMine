from __future__ import annotations

import json
from pathlib import Path
from dataclasses import asdict

from .gene_prediction import GenePredictor, create_predictor
from .io import checksum, read_fasta
from .genbank import write_package
from .models import SubmissionMetadata
from .mining import mine, ranked_candidates
from .reporting import write_outputs, write_checkpoint_snapshot, write_stage_checkpoint
from .quality import assess
from .genome_representation import GenomeRepresentation
from .sequencing_provenance import SequencingProvenance
from .annotation import MockEvidenceBackend
from .pfam import PfamHMMAdapter
from .vog import VOGHMMAdapter
from .swissprot import SwissProtEvidenceAdapter
from .phrogs import PHROGSMMseqsAdapter
from .resources import EvidenceResourceManager, ResourceType
from .progress import ProgressReporter
from .fusion import classify_proteins
from .context import build_context
from .reconciliation import GeneModel, ProdigalPredictor, reconcile_models, write_reconciliation
from .adjudication import adjudicate, write_adjudication
from .alternative_evidence import alternative_models, acquire_alternative_evidence, write_alternative_evidence


def run(fasta: str | Path, output: str | Path, command: str = "run", metadata: SubmissionMetadata | None = None, table2asn_executable: str | None = None, predictor: GenePredictor | None = None, representation: GenomeRepresentation | None = None, sequencing_provenance: SequencingProvenance | None = None, pfam_path: str | Path | None = None, pfam_hmmscan: str | None = None, pfam_evalue: float | None = None, pfam_coverage: float | None = None, pfam_trusted_cutoff: bool = False, use_mock_evidence: bool = False, pfam_threshold_mode: str | None = None, vog_path: str | Path | None = None, vog_annotations: str | Path | None = None, vog_hmmscan: str | None = None, vog_evalue: float | None = 1e-5, vog_coverage: float | None = 0.5, swissprot_path: str | Path | None = None, swissprot_metadata: str | Path | None = None, diamond: str | None = None, swissprot_evalue: float = 1e-5, phrogs_path: str | Path | None = None, phrogs_annotations: str | Path | None = None, mmseqs: str | None = None, phrogs_evalue: float | None = 1e-5, phrogs_coverage: float | None = 0.5, phrogs_score: float | None = None, phrogs_identity: float | None = None, phrogs_alignment_length: int | None = None, reconcile_orfs: bool = False, prodigal: str | None = None, progress: ProgressReporter | None = None) -> int:
    progress = progress or ProgressReporter(quiet=True)
    progress.start("input/genome validation")
    genome_id, genome = read_fasta(fasta)
    representation = representation or GenomeRepresentation.original(genome_id, genome)
    sequencing_provenance = sequencing_provenance or SequencingProvenance()
    if representation.original_sequence_id != genome_id or representation.original_sequence != genome:
        raise ValueError("GenomeRepresentation must be derived from the supplied authoritative input FASTA.")
    progress.finish(f"genome {genome_id}; {len(genome):,} bp")
    predictor = predictor or create_predictor("phanotate")
    predictor_input = Path(fasta)
    if representation.analysis_sequence != genome or representation.analysis_sequence_id != genome_id:
        Path(output).mkdir(parents=True, exist_ok=True)
        predictor_input = Path(output) / "analysis_predictor_input.fasta"
        predictor_input.write_text(f">{representation.analysis_sequence_id}\n{representation.analysis_sequence}\n")
    progress.start("gene prediction")
    proteins = predictor.predict(representation.analysis_sequence_id, representation.analysis_sequence, predictor_input)
    if not proteins:
        raise ValueError("No ORFs met the MVP minimum length; use a genome with coding sequences or lower the configured threshold in a future adapter.")
    gene_manifest = {"stage": "gene_prediction", "gene_caller": {"name": predictor.name, "version": predictor.version(), "parameters": predictor.parameters()}, "input_sha256": checksum(fasta)}
    write_checkpoint_snapshot(output, Path(output) / "checkpoints" / "gene_prediction", representation, sequencing_provenance, proteins, gene_manifest, fasta)
    reconciliation_rows = None
    prodigal_models = None
    if reconcile_orfs:
        progress.start("Prodigal secondary gene prediction")
        prodigal_models = ProdigalPredictor(prodigal).predict(fasta, representation.analysis_sequence)
        progress.finish(f"{len(prodigal_models)} proteins")
        progress.start("ORF reconciliation")
        phanotate_models = [GeneModel("PHANOTATE", p.protein_id, p.start, p.end, p.strand, p.cds, caller_version=predictor.version(), options=p.gene_call_parameters) for p in proteins]
        reconciliation_rows = reconcile_models(phanotate_models, prodigal_models, checksum(fasta))
        write_reconciliation(output, reconciliation_rows, {"input_sha256": checksum(fasta), "phanotate": predictor.parameters(), "prodigal": {"executable": prodigal or "PATH"}})
        (Path(output) / "checkpoints" / "orf_reconciliation").mkdir(parents=True, exist_ok=True)
        (Path(output) / "checkpoints" / "orf_reconciliation" / "predictions.json").write_text(json.dumps([m.__dict__ for m in prodigal_models], indent=2, sort_keys=True))
        progress.finish("reconciliation persisted")
    progress.finish(f"{len(proteins)} proteins")
    evidence_adapters = []
    if use_mock_evidence:
        mock_result = MockEvidenceBackend().analyze(proteins)
        evidence_adapters.append({"adapter": mock_result.adapter, "status": mock_result.status, "provenance": mock_result.provenance, "message": mock_result.message})
    pfam_origin = "explicit_cli" if pfam_path else "unavailable"
    if pfam_path is None:
        registered = EvidenceResourceManager().find(ResourceType.PFAM)
        if registered:
            pfam_path = registered["path"]
            pfam_origin = "registered_resource"
    pfam_adapter = PfamHMMAdapter(pfam_path, pfam_hmmscan, pfam_evalue, pfam_coverage, pfam_trusted_cutoff, threshold_mode=pfam_threshold_mode)
    progress.start("Pfam")
    pfam_result = pfam_adapter.analyze(proteins)
    pfam_result.provenance["resource_origin"] = pfam_origin
    proteins_by_id = {protein.protein_id: protein for protein in proteins}
    for evidence in pfam_result.evidence:
        # The query identifier is authoritative in parsed HMMER metrics. Keep
        # provenance-only adapters compatible, but do not drop a valid hit if
        # an adapter has not duplicated that identifier into provenance yet.
        protein_id = evidence.provenance.get("protein_id") or evidence.metrics.get("query_protein_id")
        if protein_id in proteins_by_id:
            evidence.provenance.setdefault("protein_id", protein_id)
            proteins_by_id[protein_id].evidence.append(evidence)
    evidence_adapters.append({"adapter": pfam_result.adapter, "status": pfam_result.status, "provenance": pfam_result.provenance, "message": pfam_result.message})
    write_stage_checkpoint(Path(output) / "checkpoints", "pfam", [asdict(e) for e in pfam_result.evidence], pfam_result.provenance)
    progress.finish(f"{sum(e.supports for e in pfam_result.evidence)} accepted hits / {len({e.provenance.get('protein_id') for e in pfam_result.evidence if e.supports})} proteins")
    vog_origin = "explicit_cli" if vog_path else "unavailable"
    vog_version = None
    if vog_path is None:
        registered_vog = EvidenceResourceManager().find(ResourceType.VOGDB)
        if registered_vog:
            vog_path = registered_vog["path"]
            vog_origin = "registered_resource"
            vog_version = registered_vog.get("version")
            vog_annotations = vog_annotations or registered_vog.get("provenance", {}).get("annotations_path")
    vog_adapter = VOGHMMAdapter(vog_path, vog_annotations, vog_hmmscan, vog_evalue, vog_coverage, database_version=vog_version)
    progress.start("VOGDB")
    vog_result = vog_adapter.analyze(proteins)
    vog_result.provenance["resource_origin"] = vog_origin
    for evidence in vog_result.evidence:
        protein_id = evidence.provenance.get("protein_id") or evidence.metrics.get("query_protein_id")
        if protein_id in proteins_by_id:
            evidence.provenance.setdefault("protein_id", protein_id)
            proteins_by_id[protein_id].evidence.append(evidence)
    evidence_adapters.append({"adapter": vog_result.adapter, "status": vog_result.status, "provenance": vog_result.provenance, "message": vog_result.message})
    write_stage_checkpoint(Path(output) / "checkpoints", "vogdb", [asdict(e) for e in vog_result.evidence], vog_result.provenance)
    progress.finish(f"{sum(e.supports for e in vog_result.evidence)} accepted hits / {len({e.provenance.get('protein_id') for e in vog_result.evidence if e.supports})} proteins")
    swiss_origin = "explicit_cli" if swissprot_path else "unavailable"
    swiss_version = None
    if swissprot_path is None:
        registered_swiss = EvidenceResourceManager().find(ResourceType.SWISSPROT)
        if registered_swiss:
            swissprot_path = registered_swiss["path"]
            swiss_origin = "registered_resource"
            swiss_version = registered_swiss.get("version")
            swissprot_metadata = swissprot_metadata or registered_swiss.get("provenance", {}).get("metadata_path")
    swiss_adapter = SwissProtEvidenceAdapter(swissprot_path, swissprot_metadata, diamond, database_version=swiss_version, evalue_threshold=swissprot_evalue)
    progress.start("Swiss-Prot")
    swiss_result = swiss_adapter.analyze(proteins)
    swiss_result.provenance["resource_origin"] = swiss_origin
    for evidence in swiss_result.evidence:
        protein_id = evidence.provenance.get("protein_id") or evidence.metrics.get("query_protein_id")
        if protein_id in proteins_by_id:
            evidence.provenance.setdefault("protein_id", protein_id)
            proteins_by_id[protein_id].evidence.append(evidence)
    evidence_adapters.append({"adapter": swiss_result.adapter, "status": swiss_result.status, "provenance": swiss_result.provenance, "message": swiss_result.message})
    write_stage_checkpoint(Path(output) / "checkpoints", "swissprot", [asdict(e) for e in swiss_result.evidence], swiss_result.provenance)
    if swiss_result.status == "UNAVAILABLE":
        progress.finish(f"UNAVAILABLE — {swiss_result.message or 'resource or DIAMOND unavailable'}")
    else:
        progress.finish(f"{sum(e.supports for e in swiss_result.evidence)} accepted hits / {len({e.provenance.get('protein_id') for e in swiss_result.evidence if e.supports})} proteins")
    phrogs_origin = "explicit_cli" if phrogs_path else "unavailable"
    phrogs_version = None
    if phrogs_path is None:
        registered_phrogs = EvidenceResourceManager().find(ResourceType.PHROGS)
        if registered_phrogs:
            phrogs_path = registered_phrogs["path"]
            phrogs_origin = "registered_resource"
            phrogs_version = registered_phrogs.get("version")
            phrogs_annotations = phrogs_annotations or registered_phrogs.get("provenance", {}).get("annotations_path")
    phrogs_adapter = PHROGSMMseqsAdapter(phrogs_path, phrogs_annotations, mmseqs, phrogs_version, phrogs_evalue, phrogs_coverage, phrogs_score, phrogs_identity, phrogs_alignment_length)
    progress.start("PHROGs")
    phrogs_result = phrogs_adapter.analyze(proteins)
    phrogs_result.provenance["resource_origin"] = phrogs_origin
    for evidence in phrogs_result.evidence:
        protein_id = evidence.provenance.get("protein_id") or evidence.metrics.get("query_protein_id")
        if protein_id in proteins_by_id:
            evidence.provenance.setdefault("protein_id", protein_id)
            proteins_by_id[protein_id].evidence.append(evidence)
    evidence_adapters.append({"adapter": phrogs_result.adapter, "status": phrogs_result.status, "provenance": phrogs_result.provenance, "message": phrogs_result.message})
    write_stage_checkpoint(Path(output) / "checkpoints", "phrogs", [asdict(e) for e in phrogs_result.evidence], phrogs_result.provenance)
    if reconcile_orfs:
        alt = alternative_models(reconciliation_rows, representation.analysis_sequence)
        progress.start("alternative ORF evidence")
        # Adapter objects are reused with isolated one-protein inputs; canonical evidence is untouched.
        alt = acquire_alternative_evidence(alt, (pfam_adapter, vog_adapter, swiss_adapter, phrogs_adapter), Path(output)/"checkpoints"/"alternative_evidence")
        write_alternative_evidence(output, alt, {"input_sha256": checksum(fasta), "cached": True})
        progress.finish("alternative evidence persisted")
        evidence_map = {p.protein_id: [asdict(e) for e in p.evidence] for p in proteins}
        evidence_map.update({m["alternative_model_id"]: m.get("evidence",[]) for m in alt})
        evidence_map.update({m["caller_id"]: m.get("evidence",[]) for m in alt})
        progress.start("ORF adjudication")
        write_adjudication(output, adjudicate(reconciliation_rows, evidence_map), {"input_sha256": checksum(fasta), "observational": True})
        progress.finish("adjudication persisted")
    progress.finish(f"{sum(e.supports for e in phrogs_result.evidence)} accepted hits / {len({e.provenance.get('protein_id') for e in phrogs_result.evidence if e.supports})} proteins")
    progress.start("evidence integration")
    checkpoint_manifest = {"pipeline": "PhageMine", "pipeline_version": "0.1.0", "command": command,
                           "input": str(fasta), "input_sha256": checksum(fasta),
                           "gene_caller": {"name": predictor.name, "version": predictor.version(), "parameters": predictor.parameters()},
                           "evidence_adapters": evidence_adapters}
    write_checkpoint_snapshot(output, Path(output) / "checkpoints" / "evidence_complete",
                              representation, sequencing_provenance, proteins, checkpoint_manifest, fasta)
    progress.finish(f"{sum(len(p.evidence) for p in proteins)} evidence records")
    classifications = classify_proteins(proteins)
    context_records, modules = build_context(proteins, classifications)
    progress.start("candidate ranking/mining")
    mine(proteins, mock=use_mock_evidence)
    candidates = ranked_candidates(proteins)
    ranking_status = "INSUFFICIENT_EVIDENCE" if candidates and not any(e.supports and e.evidence_strength in {"STRONG", "EXPERIMENTAL"} for protein in candidates for e in protein.evidence) else "RANKED"
    progress.finish(f"{len(candidates)} candidates; {ranking_status}")
    manifest = {"pipeline": "PhageMine", "pipeline_version": "0.1.0", "command": command, "input": str(fasta), "input_sha256": checksum(fasta), "genome_representation": representation.manifest(), "sequencing_provenance": sequencing_provenance.manifest(), "gene_caller": {"name": predictor.name, "version": predictor.version(), "parameters": predictor.parameters()}, "evidence_adapters": evidence_adapters, "discovery_ranking": {"status": ranking_status, "message": "Candidate prioritization was not performed because sufficient evidence was unavailable." if ranking_status == "INSUFFICIENT_EVIDENCE" else "Candidates ranked by available evidence."}}
    progress.start("QC/report generation")
    quality_control = assess(representation.analysis_sequence, proteins)
    manifest["quality_control"] = quality_control
    write_outputs(output, representation, sequencing_provenance, proteins, candidates, manifest, quality_control, fasta, classifications, context_records, modules)
    progress.finish(f"outputs written to {output}")
    # Local package generation follows annotation, mining, ranking, and QC evidence collection.
    progress.start("GenBank pre-submission package")
    write_package(output, representation.analysis_sequence_id, representation.analysis_sequence, proteins, manifest, metadata, table2asn_executable, sequencing_provenance)
    progress.finish("package generated")
    return len(proteins)
