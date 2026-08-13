from __future__ import annotations

import argparse
from pathlib import Path

from .pipeline import run
from .resume import resume
from .models import SubmissionMetadata
from .resources import EvidenceResourceManager, ResourceType
from .progress import ProgressReporter
import json


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="PhageMine: annotation followed by cautious discovery mining")
    subcommands = parser.add_subparsers(dest="command", required=True)
    databases = subcommands.add_parser("databases", help="Manage registered local evidence resources")
    database_commands = databases.add_subparsers(dest="database_command", required=True)
    register = database_commands.add_parser("register", help="Register a local evidence resource")
    register.add_argument("resource_type", choices=[item.value.lower() for item in ResourceType])
    register.add_argument("path")
    register.add_argument("--name")
    register.add_argument("--version")
    register.add_argument("--checksum")
    register.add_argument("--notes")
    register.add_argument("--annotations", help="Optional annotation table path (used by VOGDB resources)")
    register.add_argument("--metadata", help="Optional metadata path (used by Swiss-Prot resources)")
    database_commands.add_parser("status", help="List registered evidence resources")
    remove = database_commands.add_parser("remove", help="Remove a resource registration")
    remove.add_argument("name")
    resume_command = subcommands.add_parser("resume", help="Resume from a completed run and run missing evidence")
    resume_command.add_argument("source")
    resume_command.add_argument("--output", required=True)
    resume_command.add_argument("--run-missing-evidence", action="store_true")
    resume_command.add_argument("--refresh-evidence", choices=("PHROGS",))
    resume_command.add_argument("--mmseqs")
    resume_command.add_argument("--phrogs")
    resume_command.add_argument("--phrogs-annotations")
    resume_command.add_argument("--phrogs-evalue", type=float, default=1e-5)
    resume_command.add_argument("--phrogs-coverage", type=float, default=0.5)
    resume_command.add_argument("--phrogs-score", type=float)
    for name in ("annotate", "mine", "run", "genbank"):
        command = subcommands.add_parser(name, help=f"Run the MVP {name} workflow")
        command.add_argument("fasta", help="Single-record phage genome FASTA")
        command.add_argument("--output", default=None, help="Output directory (default: results/<genome stem>)")
        command.add_argument("--metadata", help="Optional JSON submission metadata; missing fields are not inferred")
        command.add_argument("--sequencing-provenance", help="Optional JSON sequencing/assembly provenance; absent values remain UNKNOWN")
        command.add_argument("--table2asn", help="Optional path to official NCBI table2asn executable")
        command.add_argument("--gene-predictor", choices=("phanotate", "demo"), default="phanotate", help="Gene caller; demo is for fixtures/tests only")
        command.add_argument("--phanotate", help="Path to PHANOTATE executable")
        command.add_argument("--pfam", help="Optional local Pfam HMM database path")
        command.add_argument("--pfam-hmmscan", help="Optional hmmscan executable path")
        command.add_argument("--pfam-evalue", type=float, help="Optional Pfam domain i-Evalue threshold")
        command.add_argument("--pfam-coverage", type=float, help="Optional Pfam query coverage threshold (0-1)")
        command.add_argument("--pfam-trusted-cutoff", action="store_true", help="Use configured Pfam trusted-cutoff policy when available")
        command.add_argument("--pfam-threshold-mode", choices=("ga", "manual", "none"), help="Pfam threshold mode")
        command.add_argument("--vogdb", help="Optional prepared combined VOGDB HMM database path")
        command.add_argument("--vog-annotations", help="Optional VOGDB annotation TSV/TSV.GZ path")
        command.add_argument("--vog-hmmscan", help="Optional hmmscan executable path for VOGDB")
        command.add_argument("--vog-evalue", type=float, default=1e-5, help="VOGDB independent E-value threshold")
        command.add_argument("--vog-coverage", type=float, default=0.5, help="VOGDB query coverage threshold")
        command.add_argument("--swissprot", help="Optional Swiss-Prot DIAMOND database path")
        command.add_argument("--swissprot-metadata", help="Optional Swiss-Prot DAT metadata path")
        command.add_argument("--diamond", help="Optional DIAMOND executable path")
        command.add_argument("--swissprot-evalue", type=float, default=1e-5, help="Swiss-Prot E-value threshold")
        command.add_argument("--phrogs", help="Optional prepared PHROGs MMseqs2 profile database prefix")
        command.add_argument("--phrogs-annotations", help="Optional PHROGs annotation table path")
        command.add_argument("--mmseqs", help="Optional MMseqs2 executable path for PHROGs")
        command.add_argument("--phrogs-evalue", type=float, default=1e-5, help="Manual PHROGs E-value threshold")
        command.add_argument("--phrogs-coverage", type=float, default=0.5, help="Manual PHROGs query coverage threshold (0-1)")
        command.add_argument("--phrogs-score", type=float, help="Optional manual PHROGs MMseqs2 bit-score threshold")
        command.add_argument("--phrogs-identity", type=float, help="Optional manual PHROGs identity threshold (MMseqs2 fident percent)")
        command.add_argument("--phrogs-alignment-length", type=int, help="Optional manual PHROGs alignment-length threshold")
        command.add_argument("--quiet", action="store_true", help="Suppress progress display")
        command.add_argument("--no-progress", action="store_true", help="Disable dynamic progress rendering")
        command.add_argument("--mock-evidence", action="store_true", help="Use demonstration evidence; fixture/testing only")
    args = parser.parse_args(argv)
    if args.command == "resume":
        try:
            count = resume(args.source, args.output, args.run_missing_evidence, args.refresh_evidence, args.mmseqs, args.phrogs, args.phrogs_annotations, args.phrogs_evalue, args.phrogs_coverage, args.phrogs_score, ProgressReporter(quiet=False))
        except (OSError, ValueError, RuntimeError) as exc:
            parser.error(str(exc))
        print(f"PhageMine resume complete: {count} predicted proteins. Outputs: {args.output}")
        return 0
    if args.command == "databases":
        manager = EvidenceResourceManager()
        if args.database_command == "register":
            name = args.name or args.resource_type.upper()
            provenance = {key: value for key, value in (("annotations_path", args.annotations), ("metadata_path", args.metadata)) if value}
            resource = manager.register(name, args.resource_type, args.path, version=args.version, checksum=args.checksum, notes=args.notes, provenance=provenance)
            print(json.dumps(resource.metadata(), indent=2, sort_keys=True))
            return 0
        if args.database_command == "status":
            print(json.dumps(manager.list(), indent=2, sort_keys=True))
            return 0
        removed = manager.unregister(args.name)
        print(f"Removed {args.name}" if removed else f"No registration found for {args.name}")
        return 0 if removed else 1
    output = args.output or str(Path("results") / Path(args.fasta).stem)
    try:
        metadata = SubmissionMetadata.from_dict(json.loads(Path(args.metadata).read_text())) if args.metadata else None
        from .sequencing_provenance import SequencingProvenance
        sequencing_provenance = SequencingProvenance.from_dict(json.loads(Path(args.sequencing_provenance).read_text())) if args.sequencing_provenance else SequencingProvenance()
        from .gene_prediction import create_predictor
        count = run(args.fasta, output, args.command, metadata, args.table2asn, create_predictor(args.gene_predictor, args.phanotate), sequencing_provenance=sequencing_provenance, pfam_path=args.pfam, pfam_hmmscan=args.pfam_hmmscan, pfam_evalue=args.pfam_evalue, pfam_coverage=args.pfam_coverage, pfam_trusted_cutoff=args.pfam_trusted_cutoff, use_mock_evidence=args.mock_evidence, pfam_threshold_mode=args.pfam_threshold_mode, vog_path=args.vogdb, vog_annotations=args.vog_annotations, vog_hmmscan=args.vog_hmmscan, vog_evalue=args.vog_evalue, vog_coverage=args.vog_coverage, swissprot_path=args.swissprot, swissprot_metadata=args.swissprot_metadata, diamond=args.diamond, swissprot_evalue=args.swissprot_evalue, phrogs_path=args.phrogs, phrogs_annotations=args.phrogs_annotations, mmseqs=args.mmseqs, phrogs_evalue=args.phrogs_evalue, phrogs_coverage=args.phrogs_coverage, phrogs_score=args.phrogs_score, phrogs_identity=args.phrogs_identity, phrogs_alignment_length=args.phrogs_alignment_length, progress=ProgressReporter(quiet=args.quiet, no_progress=args.no_progress))
    except (OSError, ValueError, RuntimeError) as exc:
        parser.error(str(exc))
    print(f"PhageMine complete: {count} predicted proteins. Outputs: {output}")
    print(f"GenBank pre-submission package: {Path(output) / 'genbank_submission'}")
    if args.mock_evidence:
        print("WARNING: mock evidence was enabled; it is demonstration-only and not biological evidence.")
    else:
        print("Functional evidence is limited to explicitly available adapters; unavailable evidence was not fabricated.")
    return 0
