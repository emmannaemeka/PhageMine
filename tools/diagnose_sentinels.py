#!/usr/bin/env python
"""Diagnose sentinel loci from existing caller/provenance outputs (read-only)."""
import argparse
import json
from phagemine.sentinel_diagnostic import diagnose, load_sentinels, write_report

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("run_dir", help="completed PhageMine output directory")
parser.add_argument("sentinels", help="TSV/JSON sentinel table with start/end columns")
parser.add_argument("--output", help="write JSON report; defaults to stdout")
args = parser.parse_args()
result = diagnose(args.run_dir, load_sentinels(args.sentinels))
if args.output:
    write_report(result, args.output)
else:
    print(json.dumps(result, indent=2, sort_keys=True))
