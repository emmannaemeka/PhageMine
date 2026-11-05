"""Optional local Foldseek runner for DeepOrigins; no Ask PhageMine dependency."""
from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import subprocess
from pathlib import Path

from .deeporigins_foldseek import COLUMNS, parse_foldseek


def run_foldseek(query: str, database: str, output: str, *,
                 executable: str = "foldseek", threads: int = 1,
                 timeout: int = 3600) -> dict:
    """Run Foldseek easy-search against a user-provided local structure database.

    Input query is a structure file or directory of structures, not amino-acid FASTA.
    Database must already be prepared by Foldseek.
    """
    if threads < 1 or timeout < 1:
        raise ValueError("threads and timeout must be positive")
    binary = shutil.which(executable)
    if binary is None:
        raise FileNotFoundError("Foldseek not found. Install Foldseek separately.")
    query_path = Path(query).expanduser().resolve()
    db_path = Path(database).expanduser().resolve()
    out = Path(output).expanduser().resolve()
    if not query_path.exists():
        raise FileNotFoundError(f"Query structures not found: {query_path}")
    if not db_path.exists():
        raise FileNotFoundError(f"Foldseek database not found: {db_path}")
    if out.exists():
        raise FileExistsError(f"Output already exists: {out}")
    out.parent.mkdir(parents=True, exist_ok=True)
    tmp = out.parent / (out.stem + "_foldseek_tmp")
    if tmp.exists():
        raise FileExistsError(f"Temporary directory already exists: {tmp}")
    cmd = [binary, "easy-search", str(query_path), str(db_path), str(out),
           str(tmp), "--threads", str(threads),
           "--format-output", ",".join(COLUMNS)]
    version = subprocess.run([binary, "version"], capture_output=True,
                             text=True, timeout=30, check=True).stdout.strip()
    try:
        subprocess.run(cmd, capture_output=True, text=True,
                       timeout=timeout, check=True)
        hits = parse_foldseek(out)
        digest = hashlib.sha256(out.read_bytes()).hexdigest()
        report = {
            "schema_version": "deeporigins-foldseek-run-v0",
            "status": "unvalidated_structural_matches",
            "foldseek_version": version,
            "database": str(db_path),
            "query": str(query_path),
            "command": cmd,
            "results_sha256": digest,
            "hit_count": len(hits),
            "warning": "Structural similarity alone does not establish ancestry or function.",
        }
        manifest = out.with_suffix(out.suffix + ".manifest.json")
        manifest.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
        return report
    finally:
        if tmp.is_dir():
            shutil.rmtree(tmp)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Run optional local Foldseek structure search")
    parser.add_argument("structures", help="PDB/mmCIF file or directory; not protein FASTA")
    parser.add_argument("database", help="Existing Foldseek structure database")
    parser.add_argument("--output", required=True, help="New TSV output file")
    parser.add_argument("--foldseek", default="foldseek")
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--timeout", type=int, default=3600)
    args = parser.parse_args(argv)
    report = run_foldseek(args.structures, args.database, args.output,
                          executable=args.foldseek, threads=args.threads,
                          timeout=args.timeout)
    print(f"Foldseek completed; {report['hit_count']} structural matches")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
