from __future__ import annotations

import json
import time
import csv
import hashlib
from pathlib import Path
from dataclasses import asdict

from .gene_prediction import GenePredictor, create_predictor
from .io import checksum, read_fasta
from .genbank import write_package
from .models import EvidenceLevel, SubmissionMetadata
from .mining import mine, ranked_candidates
from .reporting import write_outputs, write_checkpoint_snapshot, write_stage_checkpoint, update_comparative_report
from .quality import assess
from .genome_representation import GenomeRepresentation
from .sequencing_provenance import SequencingProvenance
from .annotation import MockEvidenceBackend
from .pfam import PfamHMMAdapter
from .vog import VOGHMMAdapter
from .swissprot import SwissProtEvidenceAdapter
from .phrogs import PHROGSMMseqsAdapter, PHROGSPyHMMERAdapter, merge_phrogs_evidence
from .resources import EvidenceResourceManager, ResourceType, resolve_validated_inphared
from .progress import ProgressReporter
from .fusion import attach_gene_call_assessments, classify_proteins
from .context import build_context
from .reconciliation import GeneModel, ProdigalPredictor, gene_call_review, reconcile_models, write_gene_call_review, write_reconciliation
from .adjudication import adjudicate, write_adjudication
from .alternative_evidence import alternative_models, acquire_alternative_evidence, write_alternative_evidence
from . import __version__
from .hallmarks import assess_hallmarks, write_hallmarks
from .review import build_annotation_review, write_annotation_review


def _protein_from_gene_model(model, protein_id: str, *, candidate_id: str | None = None, locus_id: str | None = None):
    """Materialize a normalized GeneModel as the downstream Protein record."""
    from .models import Protein
    params = dict(model.options or {})
    params.update({"candidate_id": candidate_id, "locus_id": locus_id, "selection_reason": "consensus-selected model"})
    return Protein(model.genome_id or "", protein_id, model.start, model.end, model.strand,
                   model.cds_sequence or "", model.protein_sequence or model.sequence or "",
                   model.caller, gene_call_parameters=params,
                   start_codon=model.start_codon, stop_codon=model.stop_codon)


def _consensus_provider_ids(profile: str = "standard") -> tuple[str, ...]:
    """Return the one shared v1.2 DNA consensus provider profile."""
    if profile not in {"standard", "extended"}:
        raise ValueError(f"unknown gene-model profile: {profile}")
    return ("phanotate", "pyrodigal", "prodigal_gv")


def _run_consensus_gene_models(fasta, output, representation, predictor, *, profile, progress):
    """Run providers, reconciliation, observational adjudication and selection."""
    from .gene_callers import PHANOTATEProvider, PyrodigalProvider, ProdigalGVProvider, ProviderUnavailable
    from .reconciliation_engine import run_provider_reconciliation
    from .model_adjudication import CandidateModel, candidates_from_locus, decide_locus, write_model_adjudication
    from .gene_model_selection import select_final_gene_models, write_selection_outputs, CONSENSUS

    if profile not in {"standard", "extended"}:
        raise ValueError(f"unknown gene-model profile: {profile}")
    # v1.2's single DNA consensus profile intentionally uses all three
    # caller implementations.  ``profile`` is retained as an explicit
    # provenance field; both supported profiles currently resolve to this
    # provider set (extended is reserved for future additions, not a hidden
    # fourth Prodigal vote).
    providers = [PHANOTATEProvider(getattr(predictor, "executable", None)), PyrodigalProvider(), ProdigalGVProvider()]
    unavailable = [provider.provider_id for provider in providers if not provider.available()]
    if unavailable:
        raise ProviderUnavailable("Consensus %s profile requires %s; unavailable: %s" % (profile, ", ".join(p.provider_id for p in providers), ", ".join(unavailable)))
    progress.start("consensus gene callers")
    _, loci = run_provider_reconciliation(providers, representation.analysis_sequence_id,
                                           representation.analysis_sequence, fasta,
                                           Path(output) / "gene_calls", molecule_type="dna")
    candidates = [candidate for locus in loci for candidate in candidates_from_locus(locus)]
    decisions = [decide_locus(locus, candidates_from_locus(locus)) for locus in loci]
    write_model_adjudication(Path(output) / "gene_calls", candidates, decisions)
    selections = select_final_gene_models(loci, decisions, policy=CONSENSUS, profile=profile)
    write_selection_outputs(Path(output) / "gene_calls", loci, selections, decisions)
    by_locus = {locus.locus_id: locus for locus in loci}
    root = Path(output) / "gene_calls"
    selection_counts = {
        "number_loci": len(selections),
        "number_boundary_changes": sum(item.changed_from_legacy for item in selections),
        "number_fallback_phanotate": sum(item.fallback_used for item in selections),
        "number_rescue_candidates": sum(item.selection_rule == "RESCUE_CANDIDATE_NOT_AUTOMATICALLY_SELECTED" for item in selections),
        "number_review_required": sum(item.review_required for item in selections),
    }
    (root / "gene_call_manifest.json").write_text(json.dumps({
        "schema_version": "1.2-step7b", "input_fasta": str(Path(fasta).resolve()),
        "input_sequence_sha256": hashlib.sha256(representation.analysis_sequence.encode()).hexdigest(),
        "genome_id": representation.analysis_sequence_id, "declared_molecule_type": "dna",
        "gene_model_policy": "consensus", "gene_model_profile": profile,
        "selection_policy_version": "1.0", "active_providers": [provider.provider_id for provider in providers],
        "provider_versions": {provider.provider_id: provider.version() for provider in providers},
        "reconciliation_performed": True, "model_specific_evidence_performed": False,
        "adjudication_performed": True, "selection_performed": True,
        "final_cds_source": "CONSENSUS_SELECTED_MODELS", "selection_summary": selection_counts,
    }, indent=2, sort_keys=True))
    # Consensus-specific review/rescue ledgers and a complete final trace.
    with (root / "rescue_candidates.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["locus_id", "candidate_id", "provider", "start", "end", "evidence_status", "adjudication_status", "reason"])
        for locus, decision in zip(loci, decisions):
            cmap = {candidate.candidate_id: candidate for candidate in candidates_from_locus(locus)}
            for candidate in cmap.values():
                if candidate.provider_id != "phanotate" and candidate.candidate_id != (next((s.selected_candidate_id for s in selections if s.locus_id == locus.locus_id), None)):
                    writer.writerow([locus.locus_id, candidate.candidate_id, candidate.provider_id, candidate.start, candidate.end, candidate.evidence_status, decision.decision_status, "Non-PHANOTATE caller-specific model is not automatically rescued."])
    with (root / "gene_model_review.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["locus_id", "reconciliation_class", "decision_class", "review_required", "reason"])
        for locus, decision, selection in zip(loci, decisions, selections):
            if selection.review_required:
                writer.writerow([locus.locus_id, locus.reconciliation_class, decision.decision_class, "true", selection.selection_reason])
    with (root / "final_gene_model_trace.tsv").open("w", newline="") as handle:
        columns = ["protein_id", "locus_id", "candidate_id", "start", "end", "strand", "selected_source", "changed_from_legacy", "selection_policy", "selection_rule", "selection_reason", "adjudication_class", "review_required", "input_sha256"]
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t"); writer.writeheader()
        for index, selection in enumerate(selections, 1):
            selected_candidate = next((candidate for candidate in candidates_from_locus(by_locus[selection.locus_id]) if candidate.candidate_id == selection.selected_candidate_id), None) if selection.selected_candidate_id else None
            writer.writerow({"protein_id": f"PM_{index:06d}", "locus_id": selection.locus_id, "candidate_id": selection.selected_candidate_id or "", "start": selected_candidate.start if selected_candidate else "", "end": selected_candidate.end if selected_candidate else "", "strand": selected_candidate.strand if selected_candidate else "", "selected_source": selected_candidate.provider_id if selected_candidate else "UNRESOLVED", "changed_from_legacy": str(selection.changed_from_legacy).lower(), "selection_policy": selection.selection_policy, "selection_rule": selection.selection_rule, "selection_reason": selection.selection_reason, "adjudication_class": selection.adjudication_decision_class or "", "review_required": str(selection.review_required).lower(), "input_sha256": selected_candidate.input_sequence_sha256 if selected_candidate else ""})
    selected = []
    for index, selection in enumerate(selections, 1):
        locus = by_locus[selection.locus_id]
        if not selection.selected_candidate_id:
            continue
        candidate = next((item for item in candidates_from_locus(locus) if item.candidate_id == selection.selected_candidate_id), None)
        if candidate is None:
            raise ValueError(f"Selected candidate {selection.selected_candidate_id} is not present in locus {selection.locus_id}")
        model = next(model for model in locus.candidate_models if CandidateModel.from_model(locus.locus_id, model).candidate_id == candidate.candidate_id)
        selected.append(_protein_from_gene_model(model, f"PM_{len(selected)+1:06d}", candidate_id=candidate.candidate_id, locus_id=locus.locus_id))
    if not selected:
        raise ValueError("Consensus selection produced no final models")
    progress.finish(f"{len(selected)} consensus-selected proteins")
    return selected, providers, loci, decisions, selections


def _write_gene_call_provenance(output, fasta, representation, predictor, proteins, *,
                                prodigal_predictor=None, prodigal_models=None,
                                reconciliation_enabled=False):
    """Persist raw caller streams and a machine-readable final-model trace.

    This is intentionally observational in v1.2 Step 1: the existing PHANOTATE
    proteins remain the final CDS set and no model-selection rule is changed.
    """
    root = Path(output) / "gene_calls"
    raw = root / "raw"
    raw.mkdir(parents=True, exist_ok=True)
    phanotate_raw = getattr(predictor, "last_raw_output", "")
    (raw / "phanotate.raw.txt").write_text(phanotate_raw)
    with (raw / "phanotate.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["caller", "raw_identifier", "start", "end", "strand", "length_nt", "length_aa"])
        for protein in proteins:
            writer.writerow(["PHANOTATE", protein.protein_id, protein.start, protein.end,
                             protein.strand, len(protein.cds), len(protein.sequence)])

    if reconciliation_enabled and prodigal_predictor is not None:
        (raw / "prodigal.gff").write_text(getattr(prodigal_predictor, "last_raw_gff", ""))
        with (raw / "prodigal.tsv").open("w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(["caller", "raw_identifier", "start", "end", "strand", "length_nt"])
            for model in prodigal_models or []:
                writer.writerow([model.caller, model.raw_identifier, model.start, model.end,
                                 model.strand, model.length_nt])

    input_sha = checksum(fasta)
    manifest = {
        "schema_version": "1.2-step1",
        "input_fasta": str(Path(fasta).resolve()),
        "input_sequence_sha256": input_sha,
        "analysis_sequence_sha256": hashlib.sha256(representation.analysis_sequence.encode()).hexdigest(),
        "genome_id": representation.analysis_sequence_id,
        "segment_id": None,
        "genome_length": len(representation.analysis_sequence),
        "declared_molecule_type": None,
        "selected_gene_caller_policy": "phanotate-only-legacy-compatible",
        "callers_invoked": ["PHANOTATE"] + (["Prodigal"] if reconciliation_enabled else []),
        "caller_versions": {"PHANOTATE": predictor.version(), **({"Prodigal": prodigal_predictor.version()} if reconciliation_enabled else {})},
        "exact_commands": {"PHANOTATE": getattr(predictor, "last_command", None), **({"Prodigal": getattr(prodigal_predictor, "last_command", None)} if reconciliation_enabled else {})},
        "parameters": {"PHANOTATE": predictor.parameters(), **({"Prodigal": prodigal_predictor.parameters()} if reconciliation_enabled else {})},
        "raw_output_paths": {"PHANOTATE": [str(raw / "phanotate.raw.txt"), str(raw / "phanotate.tsv")], **({"Prodigal": [str(raw / "prodigal.gff"), str(raw / "prodigal.tsv")]} if reconciliation_enabled else {})},
        "coordinate_conventions": {"PHANOTATE": "1-based-inclusive", "Prodigal": "1-based-inclusive"},
        "phagemine_version": __version__,
        "reconciliation_policy": "observational-only; PHANOTATE remains final",
        "decision_policy": "v1.1 final CDS behavior preserved",
    }
    manifest["providers"] = [{
        "provider_id": "phanotate", "name": "PHANOTATE", "version": predictor.version(),
        "role": "PRIMARY", "status": "SUCCESS", "command": getattr(predictor, "last_command", None),
        "parameters": predictor.parameters(),
    }]
    if reconciliation_enabled and prodigal_predictor is not None:
        manifest["providers"].append({
            "provider_id": "prodigal", "name": "Prodigal", "version": prodigal_predictor.version(),
            "role": "SECONDARY_OBSERVATIONAL", "status": "SUCCESS",
            "command": getattr(prodigal_predictor, "last_command", None),
            "parameters": prodigal_predictor.parameters(),
        })
    (root / "gene_call_manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True))
    with (root / "final_gene_model_trace.tsv").open("w", newline="") as handle:
        columns = ["protein_id", "locus_id", "start", "end", "strand", "selected_source",
                   "phanotate_id", "prodigal_id", "caller_agreement", "decision_class",
                   "decision_reason", "review_required", "input_sha256"]
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t")
        writer.writeheader()
        for index, protein in enumerate(proteins, 1):
            writer.writerow({"protein_id": protein.protein_id, "locus_id": f"LOCUS_{index:06d}",
                             "start": protein.start, "end": protein.end, "strand": protein.strand,
                             "selected_source": "PHANOTATE", "phanotate_id": protein.protein_id,
                             "prodigal_id": "", "caller_agreement": "PHANOTATE_FINAL",
                             "decision_class": "LEGACY_PHANOTATE_FINAL",
                             "decision_reason": "Step 1 preserves v1.1 final CDS selection",
                             "review_required": "false", "input_sha256": input_sha})


def _evidence_progress_summary(result, unavailable_fallback: str) -> str:
    """Keep resource absence distinct from a completed zero-hit search."""
    if result.status == "UNAVAILABLE":
        return f"UNAVAILABLE — {result.message or unavailable_fallback}"
    accepted = [evidence for evidence in result.evidence if evidence.supports]
    proteins = {evidence.provenance.get("protein_id") for evidence in accepted}
    return f"{len(accepted)} accepted hits / {len(proteins)} proteins"


def _run_validated_inphared(
    output: str | Path,
    manifest: dict,
    progress: ProgressReporter,
    validated_resolution: tuple[dict | None, str] | None = None,
) -> dict:
    """Run the automatic genome comparison, or record an explicit safe skip."""
    from .inphared import compare_genomes

    progress.start("INPHARED genome comparison")
    resource, reason = validated_resolution or resolve_validated_inphared()
    comparative_directory = Path(output) / "comparative"
    if resource is None:
        result = {"status": "SKIPPED", "reason": reason, "matches": []}
        progress.skip(reason)
    else:
        paths = resource["runtime_paths"]
        analysis_fasta = Path(output) / "analysis_genome.fasta"
        result = compare_genomes(
            {Path(output).name: analysis_fasta},
            mash_index=paths["mash_index"],
            metadata=paths["metadata"],
            output=comparative_directory,
            reference_fasta=paths["reference_fasta"],
        )
        result["resource_version"] = resource.get("version")
        result["resource_manifest"] = paths["manifest"]
        progress.finish(f"{len(result.get('matches') or [])} nearest-reference rows")

    manifest.setdefault("comparative_analysis", {})["inphared"] = result
    manifest["comparative_analysis"]["output_directory"] = str(comparative_directory)
    update_comparative_report(output, {"inphared": result})
    return result


def run(fasta: str | Path, output: str | Path, command: str = "run", metadata: SubmissionMetadata | None = None, table2asn_executable: str | None = None, predictor: GenePredictor | None = None, representation: GenomeRepresentation | None = None, sequencing_provenance: SequencingProvenance | None = None, pfam_path: str | Path | None = None, pfam_hmmscan: str | None = None, pfam_evalue: float | None = None, pfam_coverage: float | None = None, pfam_trusted_cutoff: bool = False, use_mock_evidence: bool = False, pfam_threshold_mode: str | None = None, vog_path: str | Path | None = None, vog_annotations: str | Path | None = None, vog_hmmscan: str | None = None, vog_evalue: float | None = 1e-5, vog_coverage: float | None = 0.5, swissprot_path: str | Path | None = None, swissprot_metadata: str | Path | None = None, diamond: str | None = None, swissprot_evalue: float = 1e-5, phrogs_path: str | Path | None = None, phrogs_annotations: str | Path | None = None, phrogs_hmm_path: str | Path | None = None, mmseqs: str | None = None, phrogs_evalue: float | None = 1e-5, phrogs_coverage: float | None = 0.5, phrogs_score: float | None = None, phrogs_identity: float | None = None, phrogs_alignment_length: int | None = None, reconcile_orfs: bool = False, prodigal: str | None = None, progress: ProgressReporter | None = None, threads: int = 1, inphared_resolution: tuple[dict | None, str] | None = None, gene_model_policy: str = "phanotate-only", gene_model_profile: str = "standard") -> int:
    progress = progress or ProgressReporter(quiet=True)
    consensus_active = gene_model_policy == "consensus"
    if reconcile_orfs and not consensus_active:
        stages = [
            "input/genome validation", "gene prediction", "Prodigal secondary gene prediction",
            "ORF reconciliation", "Pfam", "VOGDB", "Swiss-Prot", "PHROGs", "alternative ORF evidence",
            "ORF adjudication", "evidence integration", "candidate ranking/mining", "QC/report generation",
        ]
        if command in {"annotate", "run"}:
            stages.append("INPHARED genome comparison")
        stages.append("GenBank pre-submission package")
        progress.STAGES = tuple(stages)
    elif command in {"annotate", "run"}:
        stages = list(progress.STAGES)
        if "INPHARED genome comparison" not in stages:
            stages.insert(-1, "INPHARED genome comparison")
        progress.STAGES = tuple(stages)
    timings = {}
    def timed_start(name): timings[name] = {"start": time.time()}
    def timed_end(name): timings[name]["end"] = time.time(); timings[name]["seconds"] = timings[name]["end"] - timings[name]["start"]
    timed_start("genome_validation")
    progress.start("input/genome validation")
    genome_id, genome = read_fasta(fasta)
    representation = representation or GenomeRepresentation.original(genome_id, genome)
    sequencing_provenance = sequencing_provenance or SequencingProvenance()
    if representation.original_sequence_id != genome_id or representation.original_sequence != genome:
        raise ValueError("GenomeRepresentation must be derived from the supplied authoritative input FASTA.")
    progress.finish(f"genome {genome_id}; {len(genome):,} bp")
    timed_end("genome_validation")
    predictor = predictor or create_predictor("phanotate")
    predictor_input = Path(fasta)
    if representation.analysis_sequence != genome or representation.analysis_sequence_id != genome_id:
        Path(output).mkdir(parents=True, exist_ok=True)
        predictor_input = Path(output) / "analysis_predictor_input.fasta"
        predictor_input.write_text(f">{representation.analysis_sequence_id}\n{representation.analysis_sequence}\n")
    consensus_context = None
    if consensus_active:
        # Generic providers own prediction/provenance in this branch.  Legacy
        # --reconcile-orfs is intentionally not duplicated here.
        proteins, active_providers, consensus_loci, consensus_decisions, consensus_selections = _run_consensus_gene_models(
            fasta, output, representation, predictor, profile=gene_model_profile, progress=progress)
        consensus_context = (active_providers, consensus_loci, consensus_decisions, consensus_selections)
    else:
        progress.start("gene prediction")
        timed_start("phanotate")
        proteins = predictor.predict(representation.analysis_sequence_id, representation.analysis_sequence, predictor_input)
        timed_end("phanotate")
    if not proteins:
        raise ValueError("No ORFs met the MVP minimum length; use a genome with coding sequences or lower the configured threshold in a future adapter.")
    gene_manifest = {"stage": "gene_prediction", "gene_caller": {"name": predictor.name, "version": predictor.version(), "parameters": predictor.parameters()}, "input_sha256": checksum(fasta), "gene_model_policy": gene_model_policy, "gene_model_profile": gene_model_profile, "selection_policy_version": "1.0"}
    write_checkpoint_snapshot(output, Path(output) / "checkpoints" / "gene_prediction", representation, sequencing_provenance, proteins, gene_manifest, fasta)
    reconciliation_rows = None
    prodigal_models = None
    gene_review_records = []
    progress.finish(f"{len(proteins)} proteins")
    if reconcile_orfs and not consensus_active:
        progress.start("Prodigal secondary gene prediction")
        prodigal_predictor = ProdigalPredictor(prodigal)
        prodigal_models = prodigal_predictor.predict(fasta, representation.analysis_sequence)
        progress.finish(f"{len(prodigal_models)} proteins")
        progress.start("ORF reconciliation")
        phanotate_models = [GeneModel(
            caller="PHANOTATE", identifier=p.protein_id, start=p.start, end=p.end,
            strand=p.strand, sequence=p.sequence, frame=None,
            caller_version=predictor.version(), command=getattr(predictor, "last_command", None),
            options=p.gene_call_parameters, genome_id=p.genome_id,
            start_codon=p.start_codon, stop_codon=p.stop_codon,
            cds_sequence=p.cds, protein_sequence=p.sequence,
            input_sequence_sha256=checksum(fasta), raw_start=p.gene_call_parameters.get("raw_start"),
            raw_end=p.gene_call_parameters.get("raw_end"), raw_strand=p.gene_call_parameters.get("reported_strand"),
            source_file=str(getattr(predictor, "last_input_fasta", fasta)),
        ) for p in proteins]
        reconciliation_rows = reconcile_models(phanotate_models, prodigal_models, checksum(fasta))
        write_reconciliation(output, reconciliation_rows, {"input_sha256": checksum(fasta), "phanotate": predictor.parameters(), "prodigal": {"executable": prodigal or "PATH"}})
        (Path(output) / "checkpoints" / "orf_reconciliation").mkdir(parents=True, exist_ok=True)
        (Path(output) / "checkpoints" / "orf_reconciliation" / "predictions.json").write_text(json.dumps([m.__dict__ for m in prodigal_models], indent=2, sort_keys=True))
        progress.finish("reconciliation persisted")
    if not consensus_active:
        _write_gene_call_provenance(
            output, fasta, representation, predictor, proteins,
            prodigal_predictor=locals().get("prodigal_predictor"),
            prodigal_models=prodigal_models,
            reconciliation_enabled=reconcile_orfs,
        )
    # Selection is an additive policy record at this stage.  The legacy
    # PHANOTATE final path remains unchanged; consensus selection is exposed
    # through the dedicated selection layer and is opt-in for future wiring.
    manifest_path = Path(output) / "gene_calls" / "gene_call_manifest.json"
    if manifest_path.is_file():
        manifest = json.loads(manifest_path.read_text())
        manifest.update({"gene_model_policy": gene_model_policy, "gene_model_profile": gene_model_profile,
                         "selection_policy_version": "1.0"})
        manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True))
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
    pfam_adapter = PfamHMMAdapter(pfam_path, pfam_hmmscan, pfam_evalue, pfam_coverage, pfam_trusted_cutoff, threshold_mode=pfam_threshold_mode, threads=threads)
    progress.start("Pfam")
    timed_start("pfam")
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
    progress.finish(_evidence_progress_summary(pfam_result, "Pfam/HMMER unavailable"))
    timed_end("pfam")
    vog_origin = "explicit_cli" if vog_path else "unavailable"
    vog_version = None
    if vog_path is None:
        registered_vog = EvidenceResourceManager().find(ResourceType.VOGDB)
        if registered_vog:
            vog_path = registered_vog["path"]
            vog_origin = "registered_resource"
            vog_version = registered_vog.get("version")
            vog_annotations = vog_annotations or registered_vog.get("provenance", {}).get("annotations_path")
    vog_adapter = VOGHMMAdapter(vog_path, vog_annotations, vog_hmmscan, vog_evalue, vog_coverage, database_version=vog_version, threads=threads)
    progress.start("VOGDB")
    timed_start("vogdb")
    vog_result = vog_adapter.analyze(proteins)
    vog_result.provenance["resource_origin"] = vog_origin
    for evidence in vog_result.evidence:
        protein_id = evidence.provenance.get("protein_id") or evidence.metrics.get("query_protein_id")
        if protein_id in proteins_by_id:
            evidence.provenance.setdefault("protein_id", protein_id)
            proteins_by_id[protein_id].evidence.append(evidence)
    evidence_adapters.append({"adapter": vog_result.adapter, "status": vog_result.status, "provenance": vog_result.provenance, "message": vog_result.message})
    write_stage_checkpoint(Path(output) / "checkpoints", "vogdb", [asdict(e) for e in vog_result.evidence], vog_result.provenance)
    progress.finish(_evidence_progress_summary(vog_result, "VOGDB/HMMER unavailable"))
    timed_end("vogdb")
    swiss_origin = "explicit_cli" if swissprot_path else "unavailable"
    swiss_version = None
    if swissprot_path is None:
        registered_swiss = EvidenceResourceManager().find(ResourceType.SWISSPROT)
        if registered_swiss:
            swissprot_path = registered_swiss["path"]
            swiss_origin = "registered_resource"
            swiss_version = registered_swiss.get("version")
            swissprot_metadata = swissprot_metadata or registered_swiss.get("provenance", {}).get("metadata_path")
    swiss_adapter = SwissProtEvidenceAdapter(swissprot_path, swissprot_metadata, diamond, database_version=swiss_version, evalue_threshold=swissprot_evalue, threads=threads)
    progress.start("Swiss-Prot")
    timed_start("swissprot")
    swiss_result = swiss_adapter.analyze(proteins)
    swiss_result.provenance["resource_origin"] = swiss_origin
    for evidence in swiss_result.evidence:
        protein_id = evidence.provenance.get("protein_id") or evidence.metrics.get("query_protein_id")
        if protein_id in proteins_by_id:
            evidence.provenance.setdefault("protein_id", protein_id)
            proteins_by_id[protein_id].evidence.append(evidence)
    evidence_adapters.append({"adapter": swiss_result.adapter, "status": swiss_result.status, "provenance": swiss_result.provenance, "message": swiss_result.message})
    write_stage_checkpoint(Path(output) / "checkpoints", "swissprot", [asdict(e) for e in swiss_result.evidence], swiss_result.provenance)
    progress.finish(_evidence_progress_summary(swiss_result, "Swiss-Prot/DIAMOND unavailable"))
    timed_end("swissprot")
    phrogs_origin = "explicit_cli" if phrogs_path else "unavailable"
    phrogs_hmm_origin = "explicit_cli" if phrogs_hmm_path else "unavailable"
    phrogs_version = None
    if phrogs_path is None:
        registered_phrogs = EvidenceResourceManager().find(ResourceType.PHROGS)
        if registered_phrogs:
            phrogs_path = registered_phrogs["path"]
            phrogs_origin = "registered_resource"
            phrogs_version = registered_phrogs.get("version")
            registered_provenance = registered_phrogs.get("provenance", {})
            phrogs_annotations = phrogs_annotations or registered_provenance.get("annotations_path")
            if phrogs_hmm_path is None and registered_provenance.get("hmm_profiles_path"):
                phrogs_hmm_path = registered_provenance["hmm_profiles_path"]
                phrogs_hmm_origin = "registered_resource"
    phrogs_adapter = PHROGSMMseqsAdapter(phrogs_path, phrogs_annotations, mmseqs, phrogs_version, phrogs_evalue, phrogs_coverage, phrogs_score, phrogs_identity, phrogs_alignment_length, threads=threads)
    phrogs_hmm_adapter = PHROGSPyHMMERAdapter(
        phrogs_hmm_path, phrogs_annotations, phrogs_version,
        phrogs_evalue, phrogs_coverage, phrogs_score, threads=threads)
    progress.start("PHROGs")
    timed_start("phrogs")
    phrogs_result = phrogs_adapter.analyze(proteins)
    phrogs_result.provenance["resource_origin"] = phrogs_origin
    phrogs_hmm_result = phrogs_hmm_adapter.analyze(proteins)
    phrogs_hmm_result.provenance["resource_origin"] = phrogs_hmm_origin
    merged_phrogs = merge_phrogs_evidence(phrogs_result.evidence, phrogs_hmm_result.evidence)
    for evidence in merged_phrogs:
        protein_id = evidence.provenance.get("protein_id") or evidence.metrics.get("query_protein_id")
        if protein_id in proteins_by_id:
            evidence.provenance.setdefault("protein_id", protein_id)
            proteins_by_id[protein_id].evidence.append(evidence)
    evidence_adapters.append({"adapter": phrogs_result.adapter, "status": phrogs_result.status, "provenance": phrogs_result.provenance, "message": phrogs_result.message})
    evidence_adapters.append({"adapter": phrogs_hmm_result.adapter, "status": phrogs_hmm_result.status, "provenance": phrogs_hmm_result.provenance, "message": phrogs_hmm_result.message})
    write_stage_checkpoint(
        Path(output) / "checkpoints", "phrogs", [asdict(e) for e in merged_phrogs],
        {"backends": [phrogs_result.provenance, phrogs_hmm_result.provenance],
         "deduplicated_evidence_count": len(merged_phrogs)})
    accepted = sum(e.supports for e in merged_phrogs)
    backend_states = f"MMseqs2={phrogs_result.state.value}; PyHMMER={phrogs_hmm_result.state.value}"
    progress.finish(f"{accepted} accepted deduplicated hits; {backend_states}")
    timed_end("phrogs")
    if reconcile_orfs and not consensus_active:
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
        lengths = {protein.protein_id: protein.length for protein in proteins}
        lengths.update({model["caller_id"]: len(model.get("protein_sequence", "")) for model in alt})
        gene_review_records = gene_call_review(reconciliation_rows, evidence_map, lengths)
        write_gene_call_review(output, gene_review_records)
        progress.finish("adjudication persisted")
    progress.start("evidence integration")
    timed_start("evidence_fusion")
    checkpoint_manifest = {"pipeline": "PhageMine", "pipeline_version": __version__, "command": command,
                           "input": str(fasta), "input_sha256": checksum(fasta),
                           "gene_caller": {"name": predictor.name, "version": predictor.version(), "parameters": predictor.parameters()},
                           "evidence_adapters": evidence_adapters}
    write_checkpoint_snapshot(output, Path(output) / "checkpoints" / "evidence_complete",
                              representation, sequencing_provenance, proteins, checkpoint_manifest, fasta)
    progress.finish(f"{sum(len(p.evidence) for p in proteins)} annotation evidence records integrated")
    timed_end("evidence_fusion")
    classifications = attach_gene_call_assessments(
        classify_proteins(proteins), gene_review_records)
    classification_by_id = {item["protein_id"]: item for item in classifications}
    # Keep every public output on the same final product vocabulary.  The
    # previous pipeline left Protein.annotation at its constructor default,
    # causing ranking and annotation.tsv to say "Hypothetical protein" even
    # when evidence fusion had selected a defensible product.
    for protein in proteins:
        classification = classification_by_id[protein.protein_id]
        protein.annotation = classification["display_product"]
        protein.functional_confidence = classification["confidence"].title()
        if classification["functional_state"] == "KNOWN_FUNCTION":
            protein.annotation_level = EvidenceLevel.CURATED
        elif classification["functional_state"] in {"PROBABLE_FUNCTION", "FUNCTIONAL_CLASS_ONLY", "CONSERVED_UNKNOWN"}:
            protein.annotation_level = EvidenceLevel.COMPUTATIONAL
        else:
            protein.annotation_level = EvidenceLevel.HYPOTHESIS
    context_records, modules = build_context(proteins, classifications)
    hallmarks = assess_hallmarks(classifications)
    write_hallmarks(output, hallmarks)
    annotation_review = build_annotation_review(proteins, classifications, gene_review_records, hallmarks)
    write_annotation_review(output, annotation_review)
    progress.start("candidate ranking/mining")
    timed_start("ranking_mining")
    mine(proteins, mock=use_mock_evidence)
    candidates = ranked_candidates(proteins)
    ranking_status = "INSUFFICIENT_EVIDENCE" if candidates and not any(e.supports and e.evidence_strength in {"STRONG", "EXPERIMENTAL"} for protein in candidates for e in protein.evidence) else "RANKED"
    progress.finish(f"{len(candidates)} candidates; {ranking_status}")
    timed_end("ranking_mining")
    manifest = {"pipeline": "PhageMine", "pipeline_version": __version__, "command": command, "input": str(fasta), "input_sha256": checksum(fasta), "genome_representation": representation.manifest(), "sequencing_provenance": sequencing_provenance.manifest(), "gene_caller": {"name": predictor.name, "version": predictor.version(), "parameters": predictor.parameters()}, "evidence_adapters": evidence_adapters, "discovery_ranking": {"status": ranking_status, "message": "Candidate prioritization was not performed because sufficient evidence was unavailable." if ranking_status == "INSUFFICIENT_EVIDENCE" else "Candidates ranked by available evidence."}}
    manifest["hallmark_summary"] = {row["hallmark"]: row["status"] for row in hallmarks}
    manifest["annotation_review"] = {"records": len(annotation_review), "high_priority": sum(row["priority"]=="HIGH" for row in annotation_review), "path": str(Path(output)/"annotation_review.tsv")}
    progress.start("QC/report generation")
    timed_start("reporting")
    quality_control = assess(representation.analysis_sequence, proteins)
    manifest["quality_control"] = quality_control
    write_outputs(output, representation, sequencing_provenance, proteins, candidates, manifest, quality_control, fasta, classifications, context_records, modules)
    progress.finish(f"outputs written to {output}")
    timed_end("reporting")
    if command == "annotate":
        timed_start("inphared")
        _run_validated_inphared(output, manifest, progress, inphared_resolution)
        timed_end("inphared")

    # ``run`` is the complete single-genome workflow.  Installed comparative
    # resources must be consumed, not merely reported as READY by ``doctor``.
    if command == "run":
        from .discovery import build_discovery_outputs
        manager = EvidenceResourceManager()
        pmfdb_resource = manager.find(ResourceType.PMFDB)
        comparative = {"pmfdb": {"status": "PMFDB_UNAVAILABLE"},
                       "inphared": {"status": "SKIPPED", "matches": []}}
        if pmfdb_resource:
            comparative = build_discovery_outputs(
                [Path(output)], Path(output) / "comparative",
                mmseqs=mmseqs or "mmseqs",
                pmfdb=pmfdb_resource.get("path") if pmfdb_resource else None,
                inphared=None,
                progress=progress,
                source_mode="single-run",
            )
        manifest["comparative_analysis"] = {
            "output_directory": str(Path(output) / "comparative"),
            "pmfdb": comparative.get("pmfdb"),
            "inphared": comparative.get("inphared"),
        }
        timed_start("inphared")
        inphared_result = _run_validated_inphared(output, manifest, progress, inphared_resolution)
        timed_end("inphared")
        comparative["inphared"] = inphared_result
        update_comparative_report(output, comparative)
    # Local package generation follows annotation, mining, ranking, and QC evidence collection.
    progress.start("GenBank pre-submission package")
    timed_start("genbank")
    write_package(output, representation.analysis_sequence_id, representation.analysis_sequence, proteins, manifest, metadata, table2asn_executable, sequencing_provenance)
    progress.finish("package generated")
    timed_end("genbank")
    manifest["stage_timings_seconds"] = timings
    (Path(output) / "run_manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True))
    return len(proteins)
