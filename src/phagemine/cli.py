from __future__ import annotations

import argparse
from pathlib import Path

from .pipeline import run
from .models import SubmissionMetadata
import json


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="PhageMine: annotation followed by cautious discovery mining")
    subcommands = parser.add_subparsers(dest="command", required=True)
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
        command.add_argument("--mock-evidence", action="store_true", help="Use demonstration evidence; fixture/testing only")
    args = parser.parse_args(argv)
    output = args.output or str(Path("results") / Path(args.fasta).stem)
    try:
        metadata = SubmissionMetadata.from_dict(json.loads(Path(args.metadata).read_text())) if args.metadata else None
        from .sequencing_provenance import SequencingProvenance
        sequencing_provenance = SequencingProvenance.from_dict(json.loads(Path(args.sequencing_provenance).read_text())) if args.sequencing_provenance else SequencingProvenance()
        from .gene_prediction import create_predictor
        count = run(args.fasta, output, args.command, metadata, args.table2asn, create_predictor(args.gene_predictor, args.phanotate), sequencing_provenance=sequencing_provenance, pfam_path=args.pfam, pfam_hmmscan=args.pfam_hmmscan, pfam_evalue=args.pfam_evalue, pfam_coverage=args.pfam_coverage, pfam_trusted_cutoff=args.pfam_trusted_cutoff, use_mock_evidence=args.mock_evidence)
    except (OSError, ValueError, RuntimeError) as exc:
        parser.error(str(exc))
    print(f"PhageMine complete: {count} predicted proteins. Outputs: {output}")
    print(f"GenBank pre-submission package: {Path(output) / 'genbank_submission'}")
    if args.mock_evidence:
        print("WARNING: mock evidence was enabled; it is demonstration-only and not biological evidence.")
    else:
        print("Functional evidence is limited to explicitly available adapters; unavailable evidence was not fabricated.")
    return 0
