"""Prepare and query reproducible INPHARED-derived comparative resources."""
from __future__ import annotations

import csv
import gzip
import hashlib
import html
import json
import re
import shutil
import sqlite3
import subprocess
from datetime import datetime, timezone
from pathlib import Path
from typing import Callable, Iterable, Iterator


PMFDB_SCHEMA_VERSION = "1.0"
INPHARED_RELEASE = "2026-04-07"
UNKNOWN_ANNOTATION = re.compile(
    r"\b(hypothetical|unknown|uncharacteri[sz]ed|unannotated|predicted protein)\b",
    re.IGNORECASE,
)
VALID_AA = set("ABCDEFGHIKLMNPQRSTVWXYZJUO*")


class INPHAREDPreparationError(RuntimeError):
    """Raised when an INPHARED-derived resource cannot be prepared safely."""


def _sha256(path: str | Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _accession(value: str) -> str:
    text = html.unescape(value or "").strip()
    match = re.search(r">([^<>]+)</a>", text, flags=re.IGNORECASE)
    if match:
        return match.group(1).strip()
    return re.sub(r"<[^>]+>", "", text).strip()


def _annotation_state(annotation: str) -> str:
    if not annotation or UNKNOWN_ANNOTATION.search(annotation):
        return "PREDICTED_UNCHARACTERIZED"
    return "PREDICTED_FUNCTION"


def _fasta_records(handle: Iterable[str]) -> Iterator[tuple[str, str, str]]:
    identifier = description = None
    sequence: list[str] = []
    for raw in handle:
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            if identifier is not None:
                yield identifier, description or "", "".join(sequence).upper()
            header = line[1:].strip()
            fields = header.split(maxsplit=1)
            identifier = fields[0] if fields else ""
            description = fields[1] if len(fields) > 1 else ""
            sequence = []
        elif identifier is None:
            raise INPHAREDPreparationError("protein FASTA contains sequence before its first header")
        else:
            sequence.append(line)
    if identifier is not None:
        yield identifier, description or "", "".join(sequence).upper()


def _write_qc(path: Path, rows: list[tuple[str, object, str]]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["metric", "value", "status"])
        writer.writerows(rows)


def _source_manifest(artifacts: list[dict]) -> list[dict]:
    return [{key: value for key, value in item.items() if value is not None} for item in artifacts]


def build_pmfdb(
    proteins_gz: str | Path,
    mapping_gz: str | Path,
    genomes_gz: str | Path,
    output: str | Path,
    *,
    mmseqs: str,
    run: Callable[[list[str]], None],
    threads: int = 1,
    source_artifacts: list[dict] | None = None,
    release: str = INPHARED_RELEASE,
) -> dict:
    """Convert the pinned INPHARED proteins and metadata into PMFDB schema 1.0."""
    root = Path(output)
    root.mkdir(parents=True, exist_ok=True)
    database = root / ".pmfdb-build.sqlite"
    connection = sqlite3.connect(database)
    connection.execute("PRAGMA journal_mode=OFF")
    connection.execute("PRAGMA synchronous=OFF")
    connection.executescript(
        """
        CREATE TABLE mapping (
            protein_id TEXT PRIMARY KEY,
            source_genome_id TEXT,
            keywords TEXT
        );
        CREATE TABLE genomes (
            accession TEXT PRIMARY KEY,
            description TEXT,
            genome_length_kb TEXT,
            gc_percent TEXT,
            genus TEXT,
            subfamily TEXT,
            family TEXT,
            host TEXT
        );
        CREATE TABLE proteins (
            protein_id TEXT PRIMARY KEY,
            annotation TEXT,
            sequence_length INTEGER,
            sequence_sha256 TEXT,
            sequence TEXT
        );
        """
    )
    duplicate_mappings = duplicate_proteins = invalid_proteins = 0
    mapping_count = genome_count = protein_count = 0
    try:
        with gzip.open(mapping_gz, "rt", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle)
            batch = []
            for row in reader:
                protein_id = (row.get("protein_id") or "").strip()
                genome_id = (row.get("contig_id") or row.get("genome_id") or "").strip()
                if not protein_id or not genome_id:
                    continue
                batch.append((protein_id, genome_id, (row.get("keywords") or "").strip()))
                if len(batch) >= 10000:
                    before = connection.total_changes
                    connection.executemany("INSERT OR IGNORE INTO mapping VALUES (?,?,?)", batch)
                    duplicate_mappings += len(batch) - (connection.total_changes - before)
                    mapping_count += connection.total_changes - before
                    batch.clear()
            if batch:
                before = connection.total_changes
                connection.executemany("INSERT OR IGNORE INTO mapping VALUES (?,?,?)", batch)
                duplicate_mappings += len(batch) - (connection.total_changes - before)
                mapping_count += connection.total_changes - before

        with gzip.open(genomes_gz, "rt", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            batch = []
            for row in reader:
                accession = _accession(row.get("Accession", ""))
                if not accession:
                    continue
                batch.append((
                    accession,
                    (row.get("Description") or "").strip(),
                    (row.get("Genome Length (KB)") or "").strip(),
                    (row.get("molGC (%)") or "").strip(),
                    (row.get("Genus") or "").strip(),
                    (row.get("Sub-family") or "").strip(),
                    (row.get("Family") or "").strip(),
                    (row.get("Host") or "").strip(),
                ))
            connection.executemany("INSERT OR REPLACE INTO genomes VALUES (?,?,?,?,?,?,?,?)", batch)
            genome_count = len(batch)

        fasta = root / "reference_phage_proteins.faa"
        with gzip.open(proteins_gz, "rt", encoding="utf-8") as source:
            batch = []
            for protein_id, annotation, sequence in _fasta_records(source):
                normalized = sequence[:-1] if sequence.endswith("*") else sequence
                if not protein_id or not normalized or "*" in normalized or set(normalized) - VALID_AA:
                    invalid_proteins += 1
                    continue
                batch.append((protein_id, annotation, len(normalized), hashlib.sha256(normalized.encode()).hexdigest(), normalized))
                if len(batch) >= 10000:
                    before = connection.total_changes
                    connection.executemany("INSERT OR IGNORE INTO proteins VALUES (?,?,?,?,?)", batch)
                    duplicate_proteins += len(batch) - (connection.total_changes - before)
                    protein_count += connection.total_changes - before
                    batch.clear()
            if batch:
                before = connection.total_changes
                connection.executemany("INSERT OR IGNORE INTO proteins VALUES (?,?,?,?,?)", batch)
                duplicate_proteins += len(batch) - (connection.total_changes - before)
                protein_count += connection.total_changes - before
        connection.commit()

        with fasta.open("w") as destination:
            for protein_id, annotation, sequence in connection.execute(
                "SELECT protein_id,annotation,sequence FROM proteins ORDER BY rowid"
            ):
                destination.write(f">{protein_id} {annotation}\n")
                for start in range(0, len(sequence), 80):
                    destination.write(sequence[start:start + 80] + "\n")

        missing_mapping = connection.execute(
            "SELECT COUNT(*) FROM proteins p LEFT JOIN mapping m USING(protein_id) WHERE m.protein_id IS NULL"
        ).fetchone()[0]
        missing_genome_metadata = connection.execute(
            """SELECT COUNT(*) FROM proteins p
               JOIN mapping m USING(protein_id)
               LEFT JOIN genomes g ON g.accession=m.source_genome_id
               WHERE g.accession IS NULL"""
        ).fetchone()[0]

        metadata = root / "reference_metadata.tsv"
        columns = [
            "external_protein_id", "source_genome_id", "phage_description", "host_genus",
            "phage_taxonomy", "phage_genus", "phage_subfamily", "phage_family",
            "annotation", "annotation_status", "characterized", "evidence_type",
            "sequence_length", "sequence_sha256", "source_database", "source_release",
        ]
        with metadata.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t")
            writer.writeheader()
            cursor = connection.execute(
                """SELECT p.protein_id,m.source_genome_id,g.description,g.host,g.genus,g.subfamily,g.family,
                          p.annotation,p.sequence_length,p.sequence_sha256
                   FROM proteins p LEFT JOIN mapping m USING(protein_id)
                   LEFT JOIN genomes g ON g.accession=m.source_genome_id ORDER BY p.rowid"""
            )
            for protein_id, genome_id, description, host, genus, subfamily, family, annotation, length, digest in cursor:
                taxonomy = ";".join(
                    value for value in (
                        f"family:{family}" if family else "",
                        f"subfamily:{subfamily}" if subfamily else "",
                        f"genus:{genus}" if genus else "",
                    ) if value
                )
                writer.writerow({
                    "external_protein_id": protein_id,
                    "source_genome_id": genome_id or "",
                    "phage_description": description or "",
                    "host_genus": host or "",
                    "phage_taxonomy": taxonomy,
                    "phage_genus": genus or "",
                    "phage_subfamily": subfamily or "",
                    "phage_family": family or "",
                    "annotation": annotation or "",
                    "annotation_status": _annotation_state(annotation or ""),
                    # INPHARED annotations are computational; do not promote them to experimental evidence.
                    "characterized": "false",
                    "evidence_type": "COMPUTATIONAL_PREDICTION",
                    "sequence_length": length,
                    "sequence_sha256": digest,
                    "source_database": "INPHARED",
                    "source_release": release,
                })

        qc = root / "reference_qc.tsv"
        qc_rows = [
            ("protein_records", protein_count, "PASS" if protein_count else "FAIL"),
            ("protein_to_genome_mappings", mapping_count, "PASS" if mapping_count else "FAIL"),
            ("reference_genomes_with_metadata", genome_count, "PASS" if genome_count else "FAIL"),
            ("invalid_proteins_excluded", invalid_proteins, "RECORDED"),
            ("duplicate_protein_ids_excluded", duplicate_proteins, "RECORDED"),
            ("duplicate_mapping_ids_excluded", duplicate_mappings, "RECORDED"),
            ("proteins_missing_genome_mapping", missing_mapping, "PASS" if not missing_mapping else "WARNING"),
            ("proteins_missing_genome_metadata", missing_genome_metadata, "PASS" if not missing_genome_metadata else "WARNING"),
            ("annotation_evidence", "computational predictions only", "CAUTION"),
        ]
        _write_qc(qc, qc_rows)
        if not protein_count or not mapping_count:
            raise INPHAREDPreparationError("INPHARED PMFDB conversion produced no usable proteins or mappings")

        mmseqs_dir = root / "mmseqs"
        mmseqs_dir.mkdir()
        target = mmseqs_dir / "target_db"
        temporary = mmseqs_dir / "tmp"
        run([mmseqs, "createdb", str(fasta), str(target)])
        if not Path(str(target) + ".dbtype").is_file():
            raise INPHAREDPreparationError("MMseqs2 did not create the PMFDB target database")
        run([mmseqs, "createindex", str(target), str(temporary), "--threads", str(max(1, int(threads)))])
        shutil.rmtree(temporary, ignore_errors=True)
        index_manifest = {
            "index_status": "SUCCESS",
            "pmfdb_version": f"INPHARED-{release}",
            "pmfdb_reference_fasta_sha256": _sha256(fasta),
            "created_at": datetime.now(timezone.utc).isoformat(),
            "commands": [
                ["mmseqs", "createdb", fasta.name, "mmseqs/target_db"],
                ["mmseqs", "createindex", "mmseqs/target_db", "mmseqs/tmp", "--threads", str(max(1, int(threads)))],
            ],
        }
        (mmseqs_dir / "mmseqs_index_manifest.json").write_text(json.dumps(index_manifest, indent=2, sort_keys=True) + "\n")

        manifest = {
            "schema_version": PMFDB_SCHEMA_VERSION,
            "pmfdb_version": f"INPHARED-{release}",
            "creation_date": datetime.now(timezone.utc).isoformat(),
            "source_database": "INPHARED",
            "source_release": release,
            "retrieval_date": datetime.now(timezone.utc).date().isoformat(),
            "filters": {
                "invalid_or_empty_protein_sequences": "excluded",
                "duplicate_protein_identifiers": "first retained and counted",
                "annotation_interpretation": "INPHARED product labels remain computational predictions",
            },
            "genome_count": genome_count,
            "protein_count": protein_count,
            "checksums": {path.name: _sha256(path) for path in (fasta, metadata, qc)},
            "source_artifacts": _source_manifest(source_artifacts or []),
        }
        (root / "reference_manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
        return manifest
    finally:
        connection.close()
        database.unlink(missing_ok=True)


def build_genome_reference(
    genomes_gz: str | Path,
    metadata_gz: str | Path,
    output: str | Path,
    *,
    mash: str,
    run: Callable[[list[str]], None],
    source_artifacts: list[dict] | None = None,
    release: str = INPHARED_RELEASE,
) -> dict:
    """Prepare the INPHARED whole-genome FASTA, metadata, and Mash sketch."""
    root = Path(output)
    root.mkdir(parents=True, exist_ok=True)
    fasta = root / "reference_phage_genomes.fna"
    genome_count = total_bases = invalid_genomes = 0
    with gzip.open(genomes_gz, "rt", encoding="utf-8") as source, fasta.open("w") as destination:
        for accession, description, sequence in _fasta_records(source):
            sequence = sequence.upper()
            if not accession or not sequence or set(sequence) - set("ACGTNURYKMSWBDHV-"):
                invalid_genomes += 1
                continue
            genome_count += 1
            total_bases += len(sequence.replace("-", ""))
            destination.write(f">{accession} {description}\n")
            for start in range(0, len(sequence), 80):
                destination.write(sequence[start:start + 80] + "\n")
    if not genome_count:
        raise INPHAREDPreparationError("INPHARED genome conversion produced no valid genomes")

    metadata = root / "genome_metadata.tsv"
    metadata_count = 0
    with gzip.open(metadata_gz, "rt", encoding="utf-8", newline="") as source, metadata.open("w", newline="") as destination:
        reader = csv.DictReader(source, delimiter="\t")
        columns = ["accession", "description", "genome_length_kb", "gc_percent", "phage_genus", "phage_subfamily", "phage_family", "host_genus", "source_database", "source_release"]
        writer = csv.DictWriter(destination, fieldnames=columns, delimiter="\t")
        writer.writeheader()
        for row in reader:
            accession = _accession(row.get("Accession", ""))
            if not accession:
                continue
            metadata_count += 1
            writer.writerow({
                "accession": accession,
                "description": (row.get("Description") or "").strip(),
                "genome_length_kb": (row.get("Genome Length (KB)") or "").strip(),
                "gc_percent": (row.get("molGC (%)") or "").strip(),
                "phage_genus": (row.get("Genus") or "").strip(),
                "phage_subfamily": (row.get("Sub-family") or "").strip(),
                "phage_family": (row.get("Family") or "").strip(),
                "host_genus": (row.get("Host") or "").strip(),
                "source_database": "INPHARED",
                "source_release": release,
            })

    prefix = root / "inphared"
    # INPHARED is a multi-FASTA reference. ``-i`` creates one sketch per
    # sequence; otherwise Mash would treat the complete release as one genome.
    run([mash, "sketch", "-i", "-o", str(prefix), str(fasta)])
    sketch = root / "inphared.msh"
    if not sketch.is_file():
        raise INPHAREDPreparationError("Mash did not create inphared.msh")
    qc = root / "genome_qc.tsv"
    _write_qc(qc, [
        ("reference_genomes", genome_count, "PASS"),
        ("reference_bases", total_bases, "PASS"),
        ("genome_metadata_records", metadata_count, "PASS" if metadata_count else "FAIL"),
        ("invalid_genomes_excluded", invalid_genomes, "RECORDED"),
        ("interpretation", "Mash distance is screening evidence, not formal taxonomy or ANI", "CAUTION"),
    ])
    manifest = {
        "schema_version": "1.0",
        "database_version": f"INPHARED-{release}",
        "creation_date": datetime.now(timezone.utc).isoformat(),
        "source_database": "INPHARED",
        "source_release": release,
        "genome_count": genome_count,
        "total_bases": total_bases,
        "checksums": {path.name: _sha256(path) for path in (fasta, metadata, qc, sketch)},
        "source_artifacts": _source_manifest(source_artifacts or []),
        "scientific_interpretation": "Mash nearest-neighbour results require confirmatory alignment/ANI and do not establish taxonomy.",
    }
    (root / "genome_manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    return manifest


def compare_genomes(
    sample_fastas: dict[str, str | Path],
    *,
    mash_index: str | Path,
    metadata: str | Path,
    output: str | Path,
    mash: str = "mash",
    top_n: int = 10,
) -> dict:
    """Find nearest INPHARED references and retain conservative provenance."""
    executable = str(Path(mash)) if Path(mash).is_file() else shutil.which(mash)
    if not executable:
        raise RuntimeError("Mash unavailable; install mash or provide its executable path")
    index = Path(mash_index)
    if not index.is_file():
        raise ValueError(f"INPHARED Mash index does not exist: {index}")
    meta = {}
    with Path(metadata).open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row.get("accession"):
                meta[row["accession"]] = row
    rows = []
    commands = []
    for sample_id, fasta in sorted(sample_fastas.items()):
        command = [executable, "dist", str(index), str(fasta)]
        commands.append(command)
        result = subprocess.run(command, capture_output=True, text=True, check=False)
        if result.returncode:
            raise RuntimeError(result.stderr.strip() or f"Mash comparison failed for {sample_id}")
        hits = []
        for line in result.stdout.splitlines():
            fields = line.split("\t")
            if len(fields) < 5:
                continue
            try:
                distance = float(fields[2])
                pvalue = float(fields[3])
            except ValueError:
                continue
            reference = fields[0].split()[0]
            hits.append((distance, pvalue, reference, fields[4]))
        for rank, (distance, pvalue, reference, shared) in enumerate(sorted(hits)[:max(1, int(top_n))], 1):
            details = meta.get(reference, {})
            rows.append({
                "sample_id": sample_id,
                "rank": rank,
                "reference_accession": reference,
                "mash_distance": distance,
                "mash_similarity_screen": 1.0 - distance,
                "p_value": pvalue,
                "matching_hashes": shared,
                "reference_description": details.get("description"),
                "host_genus": details.get("host_genus"),
                "phage_genus": details.get("phage_genus"),
                "phage_subfamily": details.get("phage_subfamily"),
                "phage_family": details.get("phage_family"),
                "interpretation": "SCREENING_ONLY_REQUIRES_CONFIRMATORY_ALIGNMENT_OR_ANI",
            })
    destination = Path(output)
    destination.mkdir(parents=True, exist_ok=True)
    columns = [
        "sample_id", "rank", "reference_accession", "mash_distance", "mash_similarity_screen",
        "p_value", "matching_hashes", "reference_description", "host_genus", "phage_genus",
        "phage_subfamily", "phage_family", "interpretation",
    ]
    with (destination / "inphared_nearest_phages.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    payload = {
        "status": "COMPLETE",
        "database": "INPHARED",
        "commands": commands,
        "scientific_interpretation": "Mash distance is nearest-neighbour screening evidence, not ANI or a taxonomic assignment.",
        "matches": rows,
    }
    (destination / "inphared_nearest_phages.json").write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    return payload
