"""Run a real core annotation through the frozen PhageMine CLI entry point."""
from __future__ import annotations

import argparse
import csv
import json
import subprocess
import tempfile
from pathlib import Path


REQUIRED_OUTPUTS = (
    "annotation.tsv",
    "evidence.json",
    "proteins.faa",
    "cds.fna",
    "genes.gff3",
    "run_manifest.json",
)


def smoke(executable: Path, fasta: Path, *, demo: bool = False) -> dict:
    with tempfile.TemporaryDirectory(prefix="phagemine-packaged-smoke-") as temporary:
        output = Path(temporary) / "core annotation"
        command = [
            str(executable), "--phagemine-cli", "run", str(fasta),
            "--output", str(output), "--threads", "1", "--no-progress",
        ]
        if demo:
            command.extend(["--gene-predictor", "demo"])
        completed = subprocess.run(
            command, capture_output=True, text=True, timeout=180,
            check=False, shell=False,
        )
        if completed.returncode:
            raise RuntimeError(
                f"Frozen core analysis failed ({completed.returncode}).\n"
                f"stdout:\n{completed.stdout}\nstderr:\n{completed.stderr}"
            )
        missing = [name for name in REQUIRED_OUTPUTS if not (output / name).is_file()]
        empty = [name for name in REQUIRED_OUTPUTS if (output / name).is_file() and not (output / name).stat().st_size]
        if missing or empty:
            raise RuntimeError(f"Frozen analysis outputs invalid; missing={missing}, empty={empty}")
        with (output / "annotation.tsv").open(newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        if not rows or not all(row.get("protein_id") for row in rows):
            raise RuntimeError("Frozen analysis produced no stable protein annotations")
        manifest = json.loads((output / "run_manifest.json").read_text())
        caller = manifest.get("gene_caller", {})
        if demo and caller.get("name") != "simple_orf_demo":
            raise RuntimeError(f"Frozen technical-preview smoke used an unexpected caller: {caller}")
        if not demo and (caller.get("name") != "PHANOTATE" or "1.6.7" not in str(caller.get("version"))):
            raise RuntimeError(f"Frozen analysis did not use bundled PHANOTATE 1.6.7: {caller}")
        return {
            "command": command,
            "protein_count": len(rows),
            "gene_caller": caller,
            "verified_outputs": list(REQUIRED_OUTPUTS),
        }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("executable", type=Path)
    parser.add_argument("fasta", type=Path)
    parser.add_argument("--demo", action="store_true", help="Packaging-only smoke for platforms without PHANOTATE")
    args = parser.parse_args()
    result = smoke(args.executable.resolve(), args.fasta.resolve(), demo=args.demo)
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
