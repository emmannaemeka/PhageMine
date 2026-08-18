"""Chunked HMMER execution with deterministic domtblout merging."""
from __future__ import annotations

import subprocess
import tempfile
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path


def _write_fasta(path: Path, records):
    path.write_text("".join(f">{pid}\n{seq}\n" for pid, seq in records))


def run_chunked(executable: str, database: str, records: list[tuple[str, str]], *, processes: int, cpus: int, parser, adapter_provenance: dict | None = None, extra_args: list[str] | None = None):
    """Run complete HMM database searches per chunk and parse through adapter logic."""
    if processes * cpus != 8:
        raise ValueError("chunk benchmark CPU budget must equal 8")
    if not records:
        return {"seconds": 0.0, "chunk_seconds": [], "commands": [], "evidence": [], "domtblout": ""}
    chunks = [records[i::processes] for i in range(min(processes, len(records)))]
    extra_args = list(extra_args or [])
    with tempfile.TemporaryDirectory(prefix="phagemine-hmmer-chunks-") as td:
        root = Path(td)
        def worker(index_records):
            index, chunk = index_records
            fasta, domtbl = root / f"chunk_{index}.faa", root / f"chunk_{index}.domtblout"
            _write_fasta(fasta, chunk)
            command = [executable, *extra_args, "--cpu", str(cpus), "--domtblout", str(domtbl), database, str(fasta)]
            started = time.perf_counter()
            result = subprocess.run(command, capture_output=True, text=True, check=False)
            elapsed = time.perf_counter() - started
            if result.returncode:
                raise RuntimeError(f"HMMER chunk {index} failed: {' '.join(command)}\n{result.stderr}")
            return index, domtbl.read_text(), elapsed, command
        started = time.perf_counter()
        with ThreadPoolExecutor(max_workers=processes) as pool:
            outputs = list(pool.map(worker, enumerate(chunks)))
        merged = "\n".join(text for _, text, _, _ in sorted(outputs))
        evidence = parser(merged, adapter_provenance or {})
        return {"seconds": time.perf_counter() - started, "chunk_seconds": [x[2] for x in outputs], "commands": [x[3] for x in sorted(outputs)], "evidence": evidence, "domtblout": merged}
