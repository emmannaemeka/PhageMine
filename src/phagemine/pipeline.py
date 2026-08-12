from __future__ import annotations

from pathlib import Path

from .annotation import MockEvidenceBackend
from .genome import predict_orfs
from .io import checksum, read_fasta
from .genbank import write_package
from .models import SubmissionMetadata
from .mining import mine, ranked_candidates
from .reporting import write_outputs
from .quality import assess


def run(fasta: str | Path, output: str | Path, command: str = "run", metadata: SubmissionMetadata | None = None, table2asn_executable: str | None = None) -> int:
    genome_id, genome = read_fasta(fasta)
    proteins = predict_orfs(genome_id, genome)
    if not proteins:
        raise ValueError("No ORFs met the MVP minimum length; use a genome with coding sequences or lower the configured threshold in a future adapter.")
    MockEvidenceBackend().annotate(proteins)
    mine(proteins)
    candidates = ranked_candidates(proteins)
    manifest = {"pipeline": "PhageMine", "pipeline_version": "0.1.0", "command": command, "input": str(fasta), "input_sha256": checksum(fasta), "gene_caller": {"name": "simple_orf_demo", "version": "0.1.0", "status": "demonstration only"}, "evidence_backend": {"name": "mock-phage-evidence", "version": "demo-1", "status": "MOCK; not biological evidence"}, "parameters": {"min_orf_nt": 90}}
    quality_control = assess(genome, proteins)
    manifest["quality_control"] = quality_control
    write_outputs(output, genome_id, genome, proteins, candidates, manifest, quality_control)
    # Local package generation follows annotation, mining, ranking, and QC evidence collection.
    write_package(output, genome_id, genome, proteins, manifest, metadata, table2asn_executable)
    return len(proteins)
