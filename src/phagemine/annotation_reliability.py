"""Observational annotation checks; never edits accepted CDS models."""
from __future__ import annotations

import csv
import hashlib
import html
import json
from pathlib import Path
import shutil
import subprocess
import tempfile


def coding_capacity(length, intervals):
    """Union coverage in 1-based inclusive coordinates, counting overlaps once."""
    merged = []
    for start, end in sorted(intervals):
        if not 1 <= start <= end <= length:
            continue
        if merged and start <= merged[-1][1] + 1:
            merged[-1][1] = max(end, merged[-1][1])
        else:
            merged.append([start, end])
    covered = sum(end - start + 1 for start, end in merged)
    return {"coding_bases": covered, "coding_density_percent": 100 * covered / length if length else None,
            "method": "union of valid CDS intervals; overlapping bases counted once"}


def _tsv(path, rows, columns):
    with Path(path).open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def caller_review(genome, proteins, secondary):
    primary = coding_capacity(len(genome), [(p.start, p.end) for p in proteins])
    records = secondary.get("records", [])
    secondary_capacity = coding_capacity(len(genome), [(r["begin"], r["end"]) for r in records]) if records else None
    tables = sorted({int(r["translation_table"]) for r in records})
    alternate = [t for t in tables if t != 11]
    short = []
    for p in proteins:
        if len(p.cds) < 100:
            agreement = any((r["begin"], r["end"], r["strand"]) == (p.start, p.end, 1 if p.strand == "+" else -1) for r in records)
            short.append({"id": p.protein_id, "source": "PRIMARY", "start": p.start, "end": p.end,
                          "strand": p.strand, "length_nt": len(p.cds),
                          "assessment": "EXACT_SECONDARY_SUPPORT" if agreement else "REVIEW_SHORT_PRIMARY_CDS",
                          "functional_evidence_sources": ",".join(sorted({e.source for e in p.evidence if e.supports}))})
    for r in records:
        if r["end"] - r["begin"] + 1 >= 100:
            continue
        strand = "+" if r["strand"] == 1 else "-"
        exact = any((p.start, p.end, p.strand) == (r["begin"], r["end"], strand) for p in proteins)
        if not exact:
            overlap = any(p.strand == strand and p.start <= r["end"] and p.end >= r["begin"] for p in proteins)
            short.append({"id": r["raw_identifier"], "source": "SECONDARY_ONLY", "start": r["begin"], "end": r["end"],
                          "strand": strand, "length_nt": r["end"] - r["begin"] + 1,
                          "assessment": "DISPUTED_BOUNDARY" if overlap else "POSSIBLE_MISSED_SHORT_CDS",
                          "functional_evidence_sources": ""})
    return {"primary": primary, "secondary": secondary_capacity,
            "alternative_code_status": "REVIEW_REQUIRED" if alternate else "NO_ALTERNATIVE_IN_NATIVE_CALLS" if records else "NOT_ASSESSED",
            "alternative_code_evidence": [{"id": r["raw_identifier"], "start": r["begin"], "end": r["end"], "strand": r["strand"], "translation_table": r["translation_table"]} for r in records if int(r["translation_table"]) != 11],
            "observed_translation_tables": tables, "alternative_translation_tables": alternate,
            "secondary_provider": secondary.get("provider_id"), "secondary_version": secondary.get("version"),
            "interpretation": "Low coding density is a review signal, not proof of an alternative code. Native secondary calls are diagnostic; final translations remain unchanged.",
            "short_cds_limit_nt": 100, "short_cds": short,
            "short_cds_limitation": "Only caller-predicted candidates are reviewed; genes missed by both callers cannot be recovered by agreement checks."}


def parse_alignments(text, query_length, reference_length):
    """Parse BLAST qseq/sseq; keep best HSP per query base for screening."""
    hsps = []
    for line in text.splitlines():
        if not line.strip():
            continue
        fields = line.split("\t")
        if len(fields) != 7:
            raise ValueError("Expected seven BLAST alignment columns")
        qs, qe, ss, se = map(int, fields[:4])
        score = float(fields[4]); qseq, sseq = fields[5:]
        if len(qseq) != len(sseq) or not (1 <= qs <= query_length and 1 <= qe <= query_length and 1 <= ss <= reference_length and 1 <= se <= reference_length):
            raise ValueError("Invalid BLAST alignment coordinates")
        qstep = 1 if qe >= qs else -1
        sstep = 1 if se >= ss else -1
        q, s = qs, ss; mapping = {}; bases = {}
        for a, b in zip(qseq.upper(), sseq.upper()):
            if a != "-":
                bases[q] = int(a == b and a in "ACGT")
                if b != "-":
                    mapping[s] = q
                q += qstep
            if b != "-":
                s += sstep
        if q - qstep != qe or s - sstep != se:
            raise ValueError("BLAST sequence/coordinate mismatch")
        hsps.append({"score": score, "mapping": mapping, "bases": bases, "orientation": qstep * sstep})
    bases = {}
    for hsp in sorted(hsps, key=lambda h: h["score"], reverse=True):
        for position, identical in hsp["bases"].items():
            bases.setdefault(position, identical)
    coverage = 100 * len(bases) / query_length if query_length else 0
    identity = 100 * sum(bases.values()) / len(bases) if bases else 0
    return coverage, identity, hsps


def reference_gene_review(record, hsps, proteins):
    """Project simple, exact-boundary CDSs through one nucleotide HSP."""
    from Bio.SeqFeature import ExactPosition
    rows = []
    for index, feature in enumerate(record.features, 1):
        if feature.type != "CDS":
            continue
        loc = feature.location
        row = {"reference_id": record.id, "reference_gene": feature.qualifiers.get("locus_tag", feature.qualifiers.get("protein_id", [str(index)]))[0],
               "reference_start": int(loc.start) + 1 if loc else "", "reference_end": int(loc.end) if loc else "",
               "query_start": "", "query_end": "", "query_strand": "", "primary_ids": "", "status": "UNMAPPABLE_REFERENCE_CDS"}
        if loc and len(loc.parts) == 1 and isinstance(loc.start, ExactPosition) and isinstance(loc.end, ExactPosition) and loc.strand in {-1, 1}:
            start, end = int(loc.start) + 1, int(loc.end)
            eligible = [h for h in hsps if all(pos in h["mapping"] for pos in range(start, end + 1))]
            projections = {(min(h["mapping"][start], h["mapping"][end]), max(h["mapping"][start], h["mapping"][end]), loc.strand * h["orientation"]) for h in eligible}
            if len(projections) > 1:
                row["status"] = "AMBIGUOUS_REFERENCE_MAPPING"
            elif projections:
                qs, qe, strand = next(iter(projections)); strand = "+" if strand == 1 else "-"
                row.update(query_start=qs, query_end=qe, query_strand=strand)
                hits = [p for p in proteins if p.strand == strand and p.start <= qe and p.end >= qs]
                row["primary_ids"] = ",".join(p.protein_id for p in hits)
                row["status"] = "EXACT_BOUNDARY_SUPPORT" if any(p.start == qs and p.end == qe for p in hits) else "BOUNDARY_DISAGREEMENT" if hits else "POSSIBLE_MISSED_CDS"
                # Insertions can disrupt the frame even when endpoints map.
                if qe - qs != end - start:
                    row["status"] = "INDEL_REQUIRES_REVIEW"
        rows.append(row)
    return rows


def reference_screen(genome, proteins, paths, blastn, min_coverage=70, min_identity=80):
    from Bio import SeqIO
    result = {"status": "NOT_REQUESTED", "candidates": [], "gene_reviews": [], "sources": [], "commands": [],
              "minimum_query_coverage": min_coverage, "minimum_identity": min_identity,
              "identity_definition": "identical A/C/G/T bases / aligned query bases, including query bases opposite deletions; best-scoring HSP per query base",
              "scope": "User-supplied annotated references. Screening cutoffs are heuristics, not taxonomy or statistical significance. No annotation transfer."}
    if not paths:
        return result
    executable = shutil.which(blastn)
    if not executable:
        result.update(status="UNAVAILABLE", reason="BLASTN executable unavailable")
        return result
    seen = set()
    with tempfile.TemporaryDirectory(prefix="phagemine-reference-review-") as tmp:
        query = Path(tmp) / "query.fna"; subject = Path(tmp) / "reference.fna"
        query.write_text(">query\n" + genome + "\n")
        for path in paths:
            path = Path(path)
            result["sources"].append({"path": str(path.resolve()), "sha256": hashlib.sha256(path.read_bytes()).hexdigest()})
            with path.open() as handle:
                reference_records = list(SeqIO.parse(handle, "genbank"))
            for record in reference_records:
                if record.id in seen:
                    raise ValueError(f"Duplicate reference ID: {record.id}")
                seen.add(record.id)
                if not len(record):
                    raise ValueError(f"Empty reference: {record.id}")
                subject.write_text(">reference\n" + str(record.seq) + "\n")
                command = [executable, "-query", str(query), "-subject", str(subject), "-outfmt", "6 qstart qend sstart send bitscore qseq sseq"]
                completed = subprocess.run(command, capture_output=True, text=True, timeout=300)
                result["commands"].append({"reference_id": record.id, "argv": command})
                if completed.returncode:
                    raise RuntimeError(completed.stderr.strip() or "Reference BLASTN failed")
                coverage, identity, hsps = parse_alignments(completed.stdout, len(genome), len(record))
                selected = coverage >= min_coverage and identity >= min_identity
                result["candidates"].append({"reference_id": record.id, "description": record.description,
                    "query_coverage_percent": coverage, "aligned_query_identity_percent": identity, "selected": selected,
                    "reference_provenance": "REFSEQ_ACCESSION" if record.id.startswith(("NC_", "AC_")) else "USER_SUPPLIED",
                    "ictv_reference_status": "NOT_VERIFIED", "sequence_sha256": hashlib.sha256(str(record.seq).encode()).hexdigest()})
                if selected:
                    result["gene_reviews"].extend(reference_gene_review(record, hsps, proteins))
    result["status"] = "COMPLETE" if seen else "NO_REFERENCE_RECORDS"
    return result


def scatter_svg(candidates, coverage, identity):
    parts = ['<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 660 470" role="img" aria-label="Reference query coverage versus identity"><rect width="660" height="470" fill="white"/>']
    for tick in range(0, 101, 20):
        x, y = 70 + tick * 5, 400 - tick * 3.5
        parts += [f'<path d="M{x} 50V400 M70 {y}H570" stroke="#ddd"/>', f'<text x="{x}" y="425" text-anchor="middle">{tick}</text>', f'<text x="60" y="{y}" text-anchor="end">{tick}</text>']
    parts += [f'<path d="M{70 + coverage * 5} 50V400 M70 {400 - identity * 3.5}H570" stroke="#999" stroke-dasharray="5 4"/>', '<text x="320" y="455" text-anchor="middle">Query coverage (%)</text>', '<text transform="translate(20 230) rotate(-90)" text-anchor="middle">Aligned-query identity (%)</text>']
    for row in candidates:
        color = "#087e8b" if row["selected"] else "#888"
        parts.append(f'<circle cx="{70 + row["query_coverage_percent"] * 5}" cy="{400 - row["aligned_query_identity_percent"] * 3.5}" r="5" fill="{color}"><title>{html.escape(row["reference_id"])}</title></circle>')
    parts.append('</svg>')
    return "".join(parts)


def write_reliability_review(output, genome, proteins, *, reference_genbank=None, blastn="blastn", min_coverage=70, min_identity=80):
    root = Path(output); target = root / "annotation_reliability"; target.mkdir(exist_ok=True)
    secondary = {}; secondary_status = "NOT_ASSESSED"
    manifest_path = root / "gene_calls" / "gene_call_manifest.json"
    if manifest_path.exists():
        manifest = json.loads(manifest_path.read_text())
        providers = [p for p in manifest.get("providers", []) if p.get("provider_id") == "prodigal_gv" and p.get("status") == "SUCCESS"]
        raw = root / "gene_calls" / "raw" / "prodigal_gv.json"
        consensus_success = manifest.get("provider_outcomes", {}).get("prodigal_gv", {}).get("status") == "SUCCESS"
        if (providers or consensus_success) and raw.exists():
            secondary = json.loads(raw.read_text()); secondary_status = "SUCCESS"
    review = caller_review(genome, proteins, secondary)
    review["secondary_status"] = secondary_status
    review["input_sequence_sha256"] = hashlib.sha256(genome.encode()).hexdigest()
    try:
        references = reference_screen(genome, proteins, reference_genbank, blastn, min_coverage, min_identity)
    except (OSError, ValueError, RuntimeError, subprocess.TimeoutExpired) as exc:
        references = {"status": "FAILED", "reason": str(exc), "candidates": [], "gene_reviews": []}
    review["references"] = references
    _tsv(target / "alternative_code_evidence.tsv", review["alternative_code_evidence"], ["id", "start", "end", "strand", "translation_table"])
    _tsv(target / "short_cds.tsv", review["short_cds"], ["id", "source", "start", "end", "strand", "length_nt", "assessment", "functional_evidence_sources"])
    _tsv(target / "reference_candidates.tsv", references["candidates"], ["reference_id", "description", "query_coverage_percent", "aligned_query_identity_percent", "selected", "reference_provenance", "ictv_reference_status", "sequence_sha256"])
    _tsv(target / "reference_gene_review.tsv", references["gene_reviews"], ["reference_id", "reference_gene", "reference_start", "reference_end", "query_start", "query_end", "query_strand", "primary_ids", "status"])
    (target / "reference_distribution.svg").write_text(scatter_svg(references["candidates"], min_coverage, min_identity))
    (target / "review.json").write_text(json.dumps(review, indent=2, sort_keys=True))
    (target / "index.html").write_text('<!doctype html><html><meta charset="utf-8"><title>Annotation reliability</title><body><h1>Annotation reliability review</h1><p>Observational evidence only; accepted CDSs and translations are unchanged. Low density does not establish an alternative code. Short CDSs missed by both callers remain undetected. Reference support is not experimental validation.</p><p><a href="short_cds.tsv">Short CDS review</a> | <a href="reference_candidates.tsv">Reference candidates</a> | <a href="reference_gene_review.tsv">Reference gene review</a></p><img src="reference_distribution.svg" alt="Reference coverage versus identity; teal points pass screening" width="660"><pre>' + html.escape(json.dumps(review, indent=2)) + '</pre></body></html>')
    return {"report": "annotation_reliability/index.html", "coding_density_percent": review["primary"]["coding_density_percent"],
            "alternative_code_status": review["alternative_code_status"], "short_cds_review_count": len(review["short_cds"]), "reference_status": references["status"]}
