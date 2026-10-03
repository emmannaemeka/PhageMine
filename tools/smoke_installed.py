"""Exercise an installed distribution from outside the source checkout.

The default synthetic fixture checks software execution, not biological accuracy.
Use --external in an environment with PHANOTATE for a real-caller Core check.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--external", action="store_true")
    args = parser.parse_args()
    source = args.input.resolve()
    output = args.output.resolve()
    if output.exists():
        parser.error("Use a new output directory so stale artifacts cannot pass this check")
    environment = os.environ.copy()
    environment.pop("PYTHONPATH", None)
    environment["PATH"] = str(Path(sys.executable).parent) + os.pathsep + environment.get("PATH", "")
    with tempfile.TemporaryDirectory(prefix="phagemine-installed-smoke-") as directory:
        environment["PHAGEMINE_REGISTRY_PATH"] = str(Path(directory) / "resources.json")
        command = [sys.executable, "-m", "phagemine", "run", str(source),
                   "--output", str(output), "--quiet"]
        if not args.external:
            command += ["--gene-predictor", "demo"]
        subprocess.run(command, cwd=directory, env=environment, check=True)
        required = ["annotation.tsv", "functional_classification.json", "evidence.json",
                    "proteins.faa", "cds.fna", "genes.gff3", "quality_control.json",
                    "report.html", "report.md", "run_manifest.json",
                    "scientific_validation_status.json", "annotation_review.tsv",
                    "genbank_submission/genome.fsa", "genbank_submission/features.tbl",
                    "genbank_submission/validation.json"]
        for name in required:
            if not (output / name).is_file() or not (output / name).stat().st_size:
                raise RuntimeError(f"Missing or empty installed-package output: {name}")
        with (output / "annotation.tsv").open(newline="") as handle:
            annotation = list(csv.DictReader(handle, delimiter="\t"))
        if not annotation:
            raise RuntimeError("The smoke input produced no annotation records")
        ids = [row["protein_id"] for row in annotation]
        evidence = json.loads((output / "evidence.json").read_text())
        classifications = json.loads((output / "functional_classification.json").read_text())
        if len(ids) != len(set(ids)) or set(ids) != {r["protein_id"] for r in evidence} or set(ids) != {r["protein_id"] for r in classifications}:
            raise RuntimeError("Annotation, evidence and classification IDs are inconsistent")
        manifest = json.loads((output / "run_manifest.json").read_text())
        if "stage_timings_seconds" not in manifest:
            raise RuntimeError("Final run manifest was not completed")
        report = {"check": "INSTALLED_PACKAGE_EXECUTION", "status": "PASS",
                  "scope": "REAL_CALLER_CORE" if args.external else "SYNTHETIC_FIXTURE",
                  "scientific_accuracy_validated": False,
                  "package_version": manifest["pipeline_version"], "command": command,
                  "input_sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
                  "protein_count": len(ids), "artifacts_checked": required}
        (output / "installed_smoke_check.json").write_text(json.dumps(report, indent=2) + "\n")
        print(json.dumps(report, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
