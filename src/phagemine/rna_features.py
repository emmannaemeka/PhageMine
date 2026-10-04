"""Optional tRNAscan-SE execution and validated noncoding-RNA imports."""
from __future__ import annotations

import csv
import hashlib
import html
import json
import math
from pathlib import Path
import shutil
import subprocess
from urllib.parse import quote, unquote

from .artifact_ops import checksum, completed_run, derived_run, read_json, record_sidecar
from .io import read_fasta_records

RNA_TYPES = {"tRNA", "tmRNA", "rRNA", "ncRNA"}


def validate_features(features: list[dict], fasta: str | Path) -> list[dict]:
    lengths = {name: len(sequence) for name, sequence in read_fasta_records(fasta)}
    seen = set()
    for feature in features:
        if feature.get("type") not in RNA_TYPES:
            raise ValueError("Only tRNA, tmRNA, rRNA and ncRNA features are supported")
        seqid = feature.get("seqid")
        if seqid not in lengths:
            raise ValueError(f"RNA sequence identifier {seqid!r} does not match the analysis FASTA")
        if not (1 <= feature["start"] <= feature["end"] <= lengths[seqid]):
            raise ValueError("RNA coordinates fall outside the analysis sequence")
        if feature.get("strand") not in {"+", "-"}:
            raise ValueError("RNA strand must be + or -")
        identity = (seqid, feature["type"], feature["start"], feature["end"], feature["strand"])
        if identity in seen:
            raise ValueError("Duplicate RNA feature")
        seen.add(identity)
        if not feature.get("source") or not feature.get("product"):
            raise ValueError("RNA features require source and product provenance")
        if any(ord(char) < 32 for char in str(feature["product"])):
            raise ValueError("RNA products must be single-line strings")
        feature["feature_id"] = "PM_RNA_" + hashlib.sha256(json.dumps(identity).encode()).hexdigest()[:12]
    return features


def parse_gff(path: str | Path, source: str, version: str) -> list[dict]:
    features = []
    for line in Path(path).read_text().splitlines():
        if line.startswith("##FASTA"): break
        if not line or line.startswith("#"): continue
        fields = line.split("\t")
        if len(fields) != 9:
            raise ValueError("RNA input is not valid GFF3")
        if fields[2] not in RNA_TYPES: continue
        attributes = {unquote(k): unquote(v) for k, v in (part.split("=", 1) for part in fields[8].split(";") if "=" in part)}
        features.append({"seqid": fields[0], "type": fields[2], "start": int(fields[3]),
                         "end": int(fields[4]), "strand": fields[6], "score": fields[5],
                         "product": attributes.get("product") or attributes.get("Name") or fields[2],
                         "source": source, "source_version": version,
                         "original_id": attributes.get("ID"), "attributes": attributes})
    if not features:
        raise ValueError("No supported RNA features in GFF3; do not treat an unsupported file as a zero-hit search")
    return features


def parse_trnascan(path: str | Path, version: str, *, confirmed_zero_hits: bool = False) -> tuple[list[dict], list[dict]]:
    features, rejected = [], []
    if confirmed_zero_hits and not Path(path).read_text().strip():
        return features, rejected
    header_seen = False
    for line in Path(path).read_text().splitlines():
        if "Sequence" in line and "tRNA" in line:
            header_seen = True
        fields = line.split()
        if len(fields) < 9 or not fields[1].isdigit(): continue
        seqid, number, begin, end, amino, anticodon, intron_begin, intron_end, score = fields[:9]
        first, last = int(begin), int(end)
        value = float(score)
        if not math.isfinite(value): raise ValueError("Invalid tRNAscan-SE score")
        if int(intron_begin) or int(intron_end) or amino.lower() in {"undet", "unknown", "pseudo"} or any("pseudo" in field.lower() for field in fields[9:]):
            rejected.append({"seqid": seqid, "number": number, "reason": "Intron-containing or uncertain tRNA requires specialist review", "raw": line})
            continue
        features.append({"seqid": seqid, "type": "tRNA", "start": min(first, last), "end": max(first, last),
                         "strand": "+" if first <= last else "-", "score": value,
                         "product": "tRNA-" + amino, "source": "tRNAscan-SE", "source_version": version,
                         "original_id": number, "anticodon": anticodon, "attributes": {}})
    if not header_seen:
        raise ValueError("Unrecognized tRNAscan-SE table format")
    return features, rejected


def append_submission_features(root: Path, features: list[dict]) -> None:
    table = root / "genbank_submission/features.tbl"
    if not table.is_file(): return
    with table.open("a") as handle:
        for feature in features:
            start, end = (feature["start"], feature["end"]) if feature["strand"] == "+" else (feature["end"], feature["start"])
            note = f"Computational RNA prediction; {feature['source']} {feature['source_version']}; record {feature['feature_id']}"
            handle.write(f"{start}\t{end}\t{feature['type']}\n\t\t\tproduct\t{feature['product']}\n\t\t\tnote\t{note}\n")
    validation = root / "genbank_submission/validation.json"
    if validation.is_file():
        payload = read_json(validation)
        payload["ncbi_table2asn_validation"] = {"state": "not_run", "message": "RNA features were added; run official validation on the current package"}
        payload["submission_readiness"] = "READY_WITH_WARNINGS" if payload.get("submission_readiness") == "READY" else payload.get("submission_readiness")
        payload["rna_feature_review"] = {"count": len(features), "status": "COMPUTATIONAL_REQUIRES_REVIEW"}
        validation.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    for pattern in ("*.sqn", "*.val", "*.gbf", "*.discr"):
        for stale in table.parent.glob(pattern): stale.unlink()


def attach(root: Path, features: list[dict], provenance: dict, rejected: list[dict] | None = None) -> dict:
    features = validate_features(features, root / "analysis_genome.fasta")
    if (root / "rna_features.json").is_file():
        raise ValueError("RNA features already attached; start from the original run")
    payload = {"status": "COMPLETED", "feature_count": len(features), "features": features,
               "rejected": rejected or [], "provenance": provenance,
               "analysis_fasta_sha256": checksum(root / "analysis_genome.fasta"),
               "interpretation": "Computational RNA predictions; zero accepted features is not proof of absence"}
    (root / "rna_features.json").write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    columns = ["feature_id", "seqid", "type", "start", "end", "strand", "product", "source", "source_version", "score"]
    with (root / "rna_features.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", extrasaction="ignore")
        writer.writeheader(); writer.writerows(features)
    gff_lines = []
    for feature in features:
        attributes = {"ID": feature["feature_id"], "product": feature["product"], "source_version": feature["source_version"]}
        attrs = ";".join(f"{key}={quote(str(value), safe='')}" for key, value in attributes.items())
        gff_lines.append(f"{feature['seqid']}\t{feature['source']}\t{feature['type']}\t{feature['start']}\t{feature['end']}\t{feature['score']}\t{feature['strand']}\t.\t{attrs}\n")
    (root / "rna_features.gff3").write_text("##gff-version 3\n" + "".join(gff_lines))
    with (root / "genes.gff3").open("a") as handle: handle.write("".join(gff_lines))
    append_submission_features(root, features)
    status_path = root / "scientific_validation_status.json"
    scientific = read_json(status_path)
    searched = provenance.get("searched_types", sorted({feature["type"] for feature in features}))
    for kind in searched: scientific["feature_scope"][kind] = "COMPUTATIONALLY_ASSESSED_REQUIRES_REVIEW"
    status_path.write_text(json.dumps(scientific, indent=2, sort_keys=True) + "\n")
    manifest = read_json(root / "run_manifest.json")
    manifest["scientific_validation"] = scientific
    manifest["rna_annotation"] = {"status": "COMPLETED", "feature_count": len(features), "provenance": provenance}
    (root / "run_manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    section = f"<section><h2>Noncoding RNA annotations</h2><p>{len(features)} computational features. {len(rejected or [])} records require separate review.</p><a href='rna_features.tsv'>RNA table</a> · <a href='rna_features.json'>Source provenance</a></section>"
    report = root / "report.html"
    report.write_text(report.read_text().replace("</body>", section + "</body>"))
    return payload


def scan(root: Path, executable: str = "tRNAscan-SE", threads: int = 1) -> dict:
    if threads < 1: raise ValueError("threads must be positive")
    program = shutil.which(executable)
    if not program: raise ValueError("tRNAscan-SE is unavailable; install trnascan-se or supply its executable path")
    directory = root / "rna_search"; directory.mkdir()
    version_probe = subprocess.run([program, "--help"], capture_output=True, text=True, timeout=30)
    if version_probe.returncode != 0: raise RuntimeError("tRNAscan-SE version probe failed")
    version_lines = (version_probe.stderr + "\n" + version_probe.stdout).strip().splitlines()
    if not version_lines: raise RuntimeError("tRNAscan-SE returned no version information")
    version = next((line.strip() for line in version_lines if line.strip().startswith("tRNAscan-SE")), version_lines[0])
    output = directory / "trnascan.tsv"
    stats = directory / "trnascan.stats"
    command = [program, "-B", "-m", str(stats), "--thread", str(threads), "-o", str(output), str(root / "analysis_genome.fasta")]
    result = subprocess.run(command, capture_output=True, text=True, timeout=3600)
    (directory / "stdout.txt").write_text(result.stdout)
    (directory / "stderr.txt").write_text(result.stderr)
    if result.returncode != 0 or not output.is_file():
        raise RuntimeError(f"tRNAscan-SE failed with exit {result.returncode}: {result.stderr[-1000:]}")
    import re
    stats_text = stats.read_text() if stats.is_file() else ""
    confirmed_zero = bool(re.search(r"(?m)^Total tRNAs:\s+0\s*$", stats_text) and
                          re.search(r"Sequences read:\s+[1-9][0-9]*", stats_text))
    features, rejected = parse_trnascan(output, version, confirmed_zero_hits=confirmed_zero)
    return attach(root, features, {"provider": "tRNAscan-SE", "version": version, "command": command,
                    "raw_table_sha256": checksum(output), "zero_hits_confirmed_by_statistics": confirmed_zero, "stats_sha256": checksum(stats) if stats.is_file() else None, "searched_types": ["tRNA"], "mode": "bacterial"}, rejected)


def scan_run(run: str | Path, output: str | Path, executable: str = "tRNAscan-SE", threads: int = 1) -> dict:
    with derived_run(run, output, "trnascan_rna_annotation") as root:
        payload = scan(root, executable, threads)
    return payload


def import_run(run: str | Path, gff: str | Path, output: str | Path, source: str, version: str) -> dict:
    if not source.strip() or not version.strip() or any(char in source + version for char in "\t\n\r"):
        raise ValueError("Supply single-line RNA source and version")
    features = parse_gff(gff, source, version)
    with derived_run(run, output, "external_rna_annotation") as root:
        payload = attach(root, features, {"provider": source, "version": version,
                         "source_gff_sha256": checksum(gff), "search_execution": "EXTERNAL_IMPORT"})
    return payload
