from __future__ import annotations

from pathlib import Path

from .annotation import MockEvidenceBackend
from .gene_prediction import GenePredictor, create_predictor
from .io import checksum, read_fasta
from .genbank import write_package
from .models import SubmissionMetadata
from .mining import mine, ranked_candidates
from .reporting import write_outputs
from .quality import assess
from .genome_representation import GenomeRepresentation
from .sequencing_provenance import SequencingProvenance


def run(fasta: str | Path, output: str | Path, command: str = "run", metadata: SubmissionMetadata | None = None, table2asn_executable: str | None = None, predictor: GenePredictor | None = None, representation: GenomeRepresentation | None = None, sequencing_provenance: SequencingProvenance | None = None) -> int:
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
    MockEvidenceBackend().annotate(proteins)
    mine(proteins)
    candidates = ranked_candidates(proteins)
    manifest = {"pipeline": "PhageMine", "pipeline_version": "0.1.0", "command": command, "input": str(fasta), "input_sha256": checksum(fasta), "genome_representation": representation.manifest(), "sequencing_provenance": sequencing_provenance.manifest(), "gene_caller": {"name": predictor.name, "version": predictor.version(), "parameters": predictor.parameters()}, "evidence_backend": {"name": "mock-phage-evidence", "version": "demo-1", "status": "MOCK; not biological evidence"}}
    quality_control = assess(representation.analysis_sequence, proteins)
    manifest["quality_control"] = quality_control
    write_outputs(output, representation, sequencing_provenance, proteins, candidates, manifest, quality_control, fasta)
    # Local package generation follows annotation, mining, ranking, and QC evidence collection.
    write_package(output, representation.analysis_sequence_id, representation.analysis_sequence, proteins, manifest, metadata, table2asn_executable, sequencing_provenance)
    return len(proteins)
