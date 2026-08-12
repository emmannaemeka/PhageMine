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
        command.add_argument("--table2asn", help="Optional path to official NCBI table2asn executable")
    args = parser.parse_args(argv)
    output = args.output or str(Path("results") / Path(args.fasta).stem)
    try:
        metadata = SubmissionMetadata.from_dict(json.loads(Path(args.metadata).read_text())) if args.metadata else None
        count = run(args.fasta, output, args.command, metadata, args.table2asn)
    except (OSError, ValueError) as exc:
        parser.error(str(exc))
    print(f"PhageMine complete: {count} predicted proteins. Outputs: {output}")
    print(f"GenBank pre-submission package: {Path(output) / 'genbank_submission'}")
    print("WARNING: this MVP run uses mock evidence and does not establish biological function.")
    return 0
