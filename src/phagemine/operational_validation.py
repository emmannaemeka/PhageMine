"""Artifact consistency and resource measurements; no accuracy claims."""
from __future__ import annotations
import csv
import json
from pathlib import Path
import platform
import subprocess
import time
from . import __version__
from .artifact_ops import checksum, completed_run, read_json
from .io import read_fasta_records


def validate_run(run: str | Path) -> dict:
    root = completed_run(run)
    proteins = read_json(root / "evidence.json", list)
    ids = [row["protein_id"] for row in proteins]
    if len(ids) != len(set(ids)): raise ValueError("Duplicate evidence protein IDs")
    expected = set(ids)
    for name in ("annotation.tsv", "functional_classification.tsv"):
        with (root / name).open(newline="") as handle: rows = list(csv.DictReader(handle, delimiter="\t"))
        observed = [row["protein_id"] for row in rows]
        if len(observed) != len(set(observed)) or set(observed) != expected: raise ValueError(f"Protein IDs do not agree in {name}")
    for name, key in (("proteins.faa", "sequence"), ("cds.fna", "cds")):
        sequences = {}; ident = None; parts = []
        for line in (root / name).read_text().splitlines():
            if line.startswith(">"):
                if ident is not None:
                    if ident in sequences: raise ValueError(f"Duplicate FASTA ID in {name}")
                    sequences[ident] = "".join(parts)
                ident = line[1:].split()[0]; parts = []
            elif ident is not None: parts.append(line.strip())
        if ident is not None:
            if ident in sequences: raise ValueError(f"Duplicate FASTA ID in {name}")
            sequences[ident] = "".join(parts)
        if set(sequences) != expected or any(sequences[row["protein_id"]] != row[key] for row in proteins): raise ValueError(f"Sequences do not agree in {name}")
    genome = dict(read_fasta_records(root / "analysis_genome.fasta"))
    from .gene_prediction import reverse_complement
    for protein in proteins:
        seq = genome.get(protein["genome_id"])
        if seq is None and len(genome) == 1: seq = next(iter(genome.values()))
        if seq is None: raise ValueError("Protein genome identifier has no analysis sequence")
        if not 1 <= protein["start"] <= protein["end"] <= len(seq): raise ValueError("Protein coordinates outside analysis sequence")
        cds = seq[protein["start"] - 1:protein["end"]]
        if protein["strand"] == "-": cds = reverse_complement(cds)
        elif protein["strand"] != "+": raise ValueError("Invalid protein strand")
        if cds != protein["cds"]: raise ValueError("CDS does not match the analysis sequence")
    gff_ids = []
    by_id = {row["protein_id"]: row for row in proteins}
    from urllib.parse import unquote
    for line in (root / "genes.gff3").read_text().splitlines():
        if not line or line.startswith("#"): continue
        fields = line.split("\t")
        if len(fields) != 9: raise ValueError("Malformed GFF3")
        if fields[2] != "CDS": continue
        attributes = {unquote(k): unquote(v) for k, v in (item.split("=", 1) for item in fields[8].split(";") if "=" in item)}
        ident = attributes.get("ID"); gff_ids.append(ident); record = by_id.get(ident)
        if record is None or (int(fields[3]), int(fields[4]), fields[6]) != (record["start"], record["end"], record["strand"]): raise ValueError("GFF3 coordinates do not match evidence records")
    if len(gff_ids) != len(set(gff_ids)) or set(gff_ids) != expected: raise ValueError("GFF3 protein IDs do not agree")
    if (root / "curation_audit.jsonl").is_file():
        from .curation import verify_audit
        verify_audit(root / "curation_audit.jsonl")
    if (root / "rna_features.json").is_file():
        from .rna_features import validate_features
        rna = read_json(root / "rna_features.json")
        if rna.get("analysis_fasta_sha256") != checksum(root / "analysis_genome.fasta"): raise ValueError("RNA feature fingerprint is stale")
        validate_features(rna["features"], root / "analysis_genome.fasta")
    return {"status": "ARTIFACTS_CONSISTENT", "protein_count": len(proteins), "gene_coordinates_checked": True,
            "curation_audit_checked": (root / "curation_audit.jsonl").is_file(), "scientific_accuracy_validated": False, "official_ncbi_validation": False}


def profile(command: list[str], output: str | Path) -> dict:
    import psutil
    if not command: raise ValueError("Supply a command after --")
    destination = Path(output).expanduser()
    if destination.exists(): raise ValueError("Resource measurement output already exists")
    started = time.monotonic(); samples = 0; peak = 0
    process = subprocess.Popen(command, shell=False); monitored = psutil.Process(process.pid)
    while process.poll() is None:
        rss = 0
        try:
            for child in [monitored, *monitored.children(recursive=True)]:
                try: rss += child.memory_info().rss
                except (psutil.NoSuchProcess, psutil.AccessDenied): pass
            samples += 1; peak = max(peak, rss)
        except (psutil.NoSuchProcess, psutil.AccessDenied): pass
        time.sleep(0.1)
    result = {"schema_version": 1, "command": command, "exit_code": process.returncode,
              "status": "COMPLETED" if process.returncode == 0 else "FAILED", "wall_seconds": time.monotonic() - started,
              "sampled_peak_process_tree_rss_bytes": peak if samples else None, "sampling_interval_seconds": 0.1, "samples": samples,
              "memory_measurement": "Sampled RSS sum: short-lived peaks can be missed; shared pages can be double-counted",
              "phagemine_version": __version__, "platform": platform.platform(), "scientific_accuracy_validated": False}
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n"); return result
