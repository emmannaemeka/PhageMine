from __future__ import annotations

from pathlib import Path

from .gene_prediction import GenePredictor, create_predictor
from .io import checksum, read_fasta
from .genbank import write_package
from .models import SubmissionMetadata
from .mining import mine, ranked_candidates
from .reporting import write_outputs
from .quality import assess
from .genome_representation import GenomeRepresentation
from .sequencing_provenance import SequencingProvenance
from .annotation import MockEvidenceBackend
from .pfam import PfamHMMAdapter


def run(fasta: str | Path, output: str | Path, command: str = "run", metadata: SubmissionMetadata | None = None, table2asn_executable: str | None = None, predictor: GenePredictor | None = None, representation: GenomeRepresentation | None = None, sequencing_provenance: SequencingProvenance | None = None, pfam_path: str | Path | None = None, pfam_hmmscan: str | None = None, pfam_evalue: float | None = None, pfam_coverage: float | None = None, pfam_trusted_cutoff: bool = False, use_mock_evidence: bool = False) -> int:
    genome_id, genome = read_fasta(fasta)
    representation = representation or GenomeRepresentation.original(genome_id, genome)
    sequencing_provenance = sequencing_provenance or SequencingProvenance()
    if representation.original_sequence_id != genome_id or representation.original_sequence != genome:
        raise ValueError("GenomeRepresentation must be derived from the supplied authoritative input FASTA.")
    predictor = predictor or create_predictor("phanotate")
    predictor_input = Path(fasta)
    if representation.analysis_sequence != genome or representation.analysis_sequence_id != genome_id:
        Path(output).mkdir(parents=True, exist_ok=True)
        predictor_input = Path(output) / "analysis_predictor_input.fasta"
        predictor_input.write_text(f">{representation.analysis_sequence_id}\n{representation.analysis_sequence}\n")
    proteins = predictor.predict(representation.analysis_sequence_id, representation.analysis_sequence, predictor_input)
    if not proteins:
        raise ValueError("No ORFs met the MVP minimum length; use a genome with coding sequences or lower the configured threshold in a future adapter.")
    evidence_adapters = []
    if use_mock_evidence:
        mock_result = MockEvidenceBackend().analyze(proteins)
        evidence_adapters.append({"adapter": mock_result.adapter, "status": mock_result.status, "provenance": mock_result.provenance, "message": mock_result.message})
    pfam_result = PfamHMMAdapter(pfam_path, pfam_hmmscan, pfam_evalue, pfam_coverage, pfam_trusted_cutoff).analyze(proteins)
    proteins_by_id = {protein.protein_id: protein for protein in proteins}
    for evidence in pfam_result.evidence:
        protein_id = evidence.provenance.get("protein_id")
        if protein_id in proteins_by_id:
            proteins_by_id[protein_id].evidence.append(evidence)
    evidence_adapters.append({"adapter": pfam_result.adapter, "status": pfam_result.status, "provenance": pfam_result.provenance, "message": pfam_result.message})
    mine(proteins, mock=use_mock_evidence)
    candidates = ranked_candidates(proteins)
    ranking_status = "INSUFFICIENT_EVIDENCE" if candidates and all(protein.biological_interest == 0 and not any(e.supports for e in protein.evidence) for protein in candidates) else "RANKED"
    manifest = {"pipeline": "PhageMine", "pipeline_version": "0.1.0", "command": command, "input": str(fasta), "input_sha256": checksum(fasta), "genome_representation": representation.manifest(), "sequencing_provenance": sequencing_provenance.manifest(), "gene_caller": {"name": predictor.name, "version": predictor.version(), "parameters": predictor.parameters()}, "evidence_adapters": evidence_adapters, "discovery_ranking": {"status": ranking_status, "message": "Candidate prioritization was not performed because sufficient evidence was unavailable." if ranking_status == "INSUFFICIENT_EVIDENCE" else "Candidates ranked by available evidence."}}
    quality_control = assess(representation.analysis_sequence, proteins)
    manifest["quality_control"] = quality_control
    write_outputs(output, representation, sequencing_provenance, proteins, candidates, manifest, quality_control, fasta)
    # Local package generation follows annotation, mining, ranking, and QC evidence collection.
    write_package(output, representation.analysis_sequence_id, representation.analysis_sequence, proteins, manifest, metadata, table2asn_executable, sequencing_provenance)
    return len(proteins)
