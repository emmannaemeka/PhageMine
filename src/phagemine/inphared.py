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
import tempfile
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
ICTV_SPECIES_THRESHOLD = 95.0
# Genus demarcation is family/proposal-specific and changes over time.  No
# boundary is embedded until a versioned, cited policy file is supplied.
ICTV_FAMILY_GENUS_THRESHOLDS: dict[str, float] = {}
ICTV_METHOD = "PhageMine bidirectional BLASTN length-normalized similarity"


def _covered_length(intervals: list[tuple[int, int]]) -> int:
    merged: list[list[int]] = []
    for start, end in sorted((min(a, b), max(a, b)) for a, b in intervals):
        if not merged or start > merged[-1][1] + 1:
            merged.append([start, end])
        else:
            merged[-1][1] = max(merged[-1][1], end)
    return sum(end - start + 1 for start, end in merged)


def _directional_blast_identity(blastn: str, query: Path, subject: Path) -> tuple[float, int, list[str]]:
    """Return non-overlapping identical bases and query coverage."""
    command = [blastn, "-query", str(query), "-subject", str(subject),
               "-word_size", "7", "-reward", "2", "-penalty", "-3",
               "-gapopen", "5", "-gapextend", "2",
               "-outfmt", "6 qstart qend length nident bitscore"]
    result = subprocess.run(command, capture_output=True, text=True, check=False)
    if result.returncode:
        raise RuntimeError(result.stderr.strip() or "BLASTN intergenomic comparison failed")
    hsps = []
    for line in result.stdout.splitlines():
        fields = line.split("\t")
        if len(fields) != 5:
            continue
        try:
            start, end, length, nident = map(int, fields[:4]); score = float(fields[4])
        except ValueError:
            continue
        hsps.append((score, min(start, end), max(start, end), length, nident))
    occupied: set[int] = set(); identical = 0.0; intervals = []
    for _score, start, end, length, nident in sorted(hsps, reverse=True):
        uncovered = sum(position not in occupied for position in range(start, end + 1))
        if not uncovered or not length:
            continue
        identical += nident * (uncovered / length)
        occupied.update(range(start, end + 1)); intervals.append((start, end))
    return identical, _covered_length(intervals), command


def _ictv_interpretation(similarity: float | None, family: str | None,
                         taxonomy_eligible: bool = True) -> tuple[float, float | None, str, str]:
    """Describe numerical boundaries without making a taxonomic assignment.

    Genus criteria are not universal across bacteriophages.  PhageMine applies
    no generic 70% genus rule and abstains for partial/poorly aligned queries.
    """
    genus = ICTV_FAMILY_GENUS_THRESHOLDS.get(family or "")
    source = ("configured family-specific working boundary; verify against the current ICTV proposal"
              if genus is not None else "no family-specific genus boundary configured")
    if similarity is None:
        label = "NOT_CALCULATED"
    elif not taxonomy_eligible:
        label = "NUMERICAL_TAXONOMY_NOT_APPLICABLE_TO_PARTIAL_OR_POORLY_ALIGNED_QUERY"
    elif similarity >= ICTV_SPECIES_THRESHOLD:
        label = "ABOVE_95_PERCENT_WORKING_SPECIES_BOUNDARY_REQUIRES_ICTV_REVIEW"
    elif genus is None:
        label = "NO_FAMILY_SPECIFIC_GENUS_BOUNDARY_CONFIGURED"
    elif similarity >= genus:
        label = "ABOVE_CONFIGURED_FAMILY_GENUS_BOUNDARY_REQUIRES_ICTV_REVIEW"
    else:
        label = "BELOW_CONFIGURED_FAMILY_GENUS_BOUNDARY"
    return ICTV_SPECIES_THRESHOLD, genus, source, label


def _calculate_intergenomic_similarity(query: str, reference: str, blastn: str) -> dict:
    with tempfile.TemporaryDirectory(prefix="phagemine-ictv-") as temporary:
        root = Path(temporary); query_path = root / "query.fna"; reference_path = root / "reference.fna"
        query_path.write_text(f">query\n{query}\n"); reference_path.write_text(f">reference\n{reference}\n")
        id_ab, aligned_query, command_ab = _directional_blast_identity(blastn, query_path, reference_path)
        id_ba, aligned_reference, command_ba = _directional_blast_identity(blastn, reference_path, query_path)
    denominator = len(query) + len(reference)
    return {
        "intergenomic_similarity_percent": ((id_ab + id_ba) * 100.0 / denominator) if denominator else None,
        "query_aligned_percent": aligned_query * 100.0 / len(query) if query else None,
        "reference_aligned_percent": aligned_reference * 100.0 / len(reference) if reference else None,
        "genome_length_ratio": min(len(query), len(reference)) / max(len(query), len(reference)) if query and reference else None,
        "similarity_method": ICTV_METHOD,
        "similarity_commands": [command_ab, command_ba],
    }


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


def _single_fasta_sequence(path: str|Path) -> str:
    sequence=[]
    for line in Path(path).read_text().splitlines():
        line=line.strip()
        if line and not line.startswith(">"):
            sequence.append(line.upper().replace("-",""))
    return "".join(sequence)


def _reference_sequences(path: str|Path, accessions: set[str]) -> dict[str,str]:
    selected={}
    if not path or not Path(path).is_file() or not accessions:
        return selected
    with Path(path).open() as handle:
        for accession, _description, sequence in _fasta_records(handle):
            if accession in accessions:
                selected[accession]=sequence.upper().replace("-","")
                if len(selected)==len(accessions): break
    return selected


def _confirm_sequence(query: str, reference: str) -> tuple[str, float|None]:
    if not query or not reference or len(query)!=len(reference):
        return "NOT_IDENTICAL_LENGTH_OR_SEQUENCE", None
    complement=str.maketrans("ACGTRYSWKMBDHVN", "TGCAYRSWMKVHDBN")
    reverse=reference.translate(complement)[::-1]
    if query==reference: return "CONFIRMED_IDENTICAL_SEQUENCE", 1.0
    if query==reverse: return "CONFIRMED_REVERSE_COMPLEMENT_SEQUENCE", 1.0
    if query in reference+reference: return "CONFIRMED_ROTATION_EQUIVALENT_SEQUENCE", 1.0
    if query in reverse+reverse: return "CONFIRMED_REVERSE_COMPLEMENT_ROTATION_EQUIVALENT", 1.0
    return "NOT_IDENTICAL_LENGTH_OR_SEQUENCE", None


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
    reference_fasta: str | Path | None = None,
    blastn: str = "blastn",
) -> dict:
    """Screen INPHARED with Mash, then calculate ICTV-compatible similarity."""
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
                "p_value": pvalue,
                "matching_hashes": shared,
                "reference_description": details.get("description"),
                "host_genus": details.get("host_genus"),
                "phage_genus": details.get("phage_genus"),
                "phage_subfamily": details.get("phage_subfamily"),
                "phage_family": details.get("phage_family"),
                "interpretation": "MASH_SCREENING_ONLY",
            })
    # Confirm zero-distance Mash results by direct nucleotide comparison when
    # the prepared reference FASTA is available.  Rotation equivalence is a
    # sequence observation and does not itself establish circular topology.
    accessions={row["reference_accession"] for row in rows}
    reference_sequences=_reference_sequences(reference_fasta,accessions)
    query_sequences={sample_id:_single_fasta_sequence(path) for sample_id,path in sample_fastas.items()}
    blastn_executable = str(Path(blastn)) if Path(blastn).is_file() else shutil.which(blastn)
    for row in rows:
        if row["mash_distance"]!=0.0:
            confirmation, identity="NOT_RUN_MASH_DISTANCE_NONZERO", None
        elif row["reference_accession"] not in reference_sequences:
            confirmation, identity="NOT_RUN_REFERENCE_SEQUENCE_UNAVAILABLE", None
        else:
            confirmation, identity=_confirm_sequence(query_sequences.get(row["sample_id"],""),reference_sequences[row["reference_accession"]])
        row["sequence_confirmation"]=confirmation
        row["confirmed_identity"]=identity
        row["confirmation_method"]="built-in exact nucleotide comparison; no ANI inferred"
        comparison = None
        if blastn_executable and row["reference_accession"] in reference_sequences:
            comparison = _calculate_intergenomic_similarity(query_sequences.get(row["sample_id"], ""), reference_sequences[row["reference_accession"]], blastn_executable)
        row.update(comparison or {"intergenomic_similarity_percent": None, "query_aligned_percent": None,
            "reference_aligned_percent": None, "genome_length_ratio": None,
            "similarity_method": "NOT_CALCULATED_BLASTN_OR_REFERENCE_UNAVAILABLE", "similarity_commands": []})
        length_ratio = row.get("genome_length_ratio")
        query_aligned = row.get("query_aligned_percent")
        reference_aligned = row.get("reference_aligned_percent")
        taxonomy_eligible = bool(
            comparison and length_ratio is not None and length_ratio >= 0.90
            and query_aligned is not None and query_aligned >= 70.0
            and reference_aligned is not None and reference_aligned >= 70.0
        )
        comparison_scope = ("APPROXIMATELY_COMPLETE_WHOLE_GENOME_COMPARISON"
                            if taxonomy_eligible else "PARTIAL_OR_POORLY_ALIGNED_QUERY")
        species_cutoff, genus_cutoff, threshold_source, taxonomic_interpretation = _ictv_interpretation(
            row["intergenomic_similarity_percent"], row.get("phage_family"), taxonomy_eligible)
        row.update({"species_threshold_percent": species_cutoff, "genus_threshold_percent": genus_cutoff,
                    "threshold_source": threshold_source, "taxonomic_interpretation": taxonomic_interpretation,
                    "taxonomy_eligible": taxonomy_eligible, "comparison_scope": comparison_scope})
        row["interpretation"] = taxonomic_interpretation if comparison else "MASH_SCREENING_ONLY; ICTV_SIMILARITY_NOT_CALCULATED"
    destination = Path(output)
    destination.mkdir(parents=True, exist_ok=True)
    columns = [
        "sample_id", "rank", "reference_accession", "mash_distance",
        "p_value", "matching_hashes", "reference_description", "host_genus", "phage_genus",
        "phage_subfamily", "phage_family", "intergenomic_similarity_percent", "query_aligned_percent",
        "reference_aligned_percent", "genome_length_ratio", "similarity_method", "species_threshold_percent",
        "genus_threshold_percent", "threshold_source", "taxonomy_eligible", "comparison_scope", "taxonomic_interpretation", "sequence_confirmation",
        "confirmed_identity", "confirmation_method", "interpretation", "similarity_commands",
    ]
    with (destination / "inphared_nearest_phages.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    # The detailed table intentionally retains accession-level records.  A
    # biological genome can, however, be represented by both GenBank and
    # RefSeq accessions.  This second table groups those duplicates so that a
    # researcher does not mistake them for independent nearest neighbours.
    unique = []
    grouped = {}
    for row in rows:
        key = (
            row["sample_id"], (row.get("reference_description") or "").strip().lower(),
            row.get("host_genus") or "", row.get("phage_genus") or "",
            row.get("phage_subfamily") or "", row.get("phage_family") or "",
        )
        grouped.setdefault(key, []).append(row)
    for group_rows in grouped.values():
        best = min(group_rows, key=lambda item: (item["mash_distance"], item["rank"]))
        shared = str(best.get("matching_hashes") or "")
        parts = shared.split("/", 1)
        exact_hashes = len(parts) == 2 and parts[0] == parts[1] and parts[0] not in {"", "0"}
        sequence_confirmed=str(best.get("sequence_confirmation") or "").startswith("CONFIRMED_")
        unique.append({
            "sample_id": best["sample_id"],
            "relationship": "CONFIRMED_SEQUENCE_EQUIVALENT_REFERENCE" if sequence_confirmed else ("EXACT_MASH_SKETCH_MATCH" if best["mash_distance"] == 0.0 and exact_hashes else "NEAREST_NEIGHBOUR_SCREEN"),
            "reference_description": best.get("reference_description"),
            "reference_accessions": ";".join(sorted({item["reference_accession"] for item in group_rows})),
            "best_mash_distance": best["mash_distance"],
            "matching_hashes": best["matching_hashes"],
            "intergenomic_similarity_percent": best.get("intergenomic_similarity_percent"),
            "query_aligned_percent": best.get("query_aligned_percent"),
            "reference_aligned_percent": best.get("reference_aligned_percent"),
            "genome_length_ratio": best.get("genome_length_ratio"),
            "similarity_method": best.get("similarity_method"),
            "species_threshold_percent": best.get("species_threshold_percent"),
            "genus_threshold_percent": best.get("genus_threshold_percent"),
            "threshold_source": best.get("threshold_source"),
            "taxonomy_eligible": best.get("taxonomy_eligible"),
            "comparison_scope": best.get("comparison_scope"),
            "taxonomic_interpretation": best.get("taxonomic_interpretation"),
            "sequence_confirmation": best.get("sequence_confirmation"),
            "confirmed_identity": best.get("confirmed_identity"),
            "host_genus": best.get("host_genus"), "phage_genus": best.get("phage_genus"),
            "phage_subfamily": best.get("phage_subfamily"), "phage_family": best.get("phage_family"),
            "interpretation": best.get("interpretation"),
        })
    unique.sort(key=lambda item: (item["sample_id"], item["best_mash_distance"], item["reference_description"] or ""))
    summary_columns = ["sample_id", "relationship", "reference_description", "reference_accessions", "best_mash_distance", "matching_hashes", "intergenomic_similarity_percent", "query_aligned_percent", "reference_aligned_percent", "genome_length_ratio", "similarity_method", "species_threshold_percent", "genus_threshold_percent", "threshold_source", "taxonomy_eligible", "comparison_scope", "taxonomic_interpretation", "sequence_confirmation", "confirmed_identity", "host_genus", "phage_genus", "phage_subfamily", "phage_family", "interpretation"]
    with (destination / "inphared_summary.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=summary_columns, delimiter="\t")
        writer.writeheader(); writer.writerows(unique)
    payload = {
        "status": "COMPLETE",
        "database": "INPHARED",
        "commands": commands,
        "scientific_interpretation": "Mash selects candidate references only. PhageMine reports its own bidirectional BLASTN length-normalized similarity and does not claim VIRIDIC equivalence or make an ICTV assignment. Taxonomic boundary interpretation is withheld for partial or poorly aligned queries.",
        "matches": rows,
        "unique_reference_summaries": unique,
    }
    (destination / "inphared_nearest_phages.json").write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    return payload
