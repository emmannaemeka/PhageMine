"""Import external annotation evidence with exact sequence/model identity checks."""
from __future__ import annotations

import csv
import json
from pathlib import Path
import shutil

from .artifact_ops import checksum, completed_run, derived_run, read_json, record_sidecar
from .io import read_fasta

UNKNOWN = {"", "hypothetical protein", "unknown function", "uncharacterized protein", "uncharacterised protein"}


def import_genbank(run: str | Path, genbank: str | Path, output: str | Path, *,
                   source: str = "Phold", version: str, database_version: str) -> dict:
    """Retain proposals without converting imported predictions into truth."""
    from Bio import SeqIO
    root = completed_run(run)
    for field, value in (("source", source), ("version", version), ("database_version", database_version)):
        if not isinstance(value, str) or not value.strip() or any(ord(char) < 32 for char in value):
            raise ValueError(f"{field} must be a nonempty single-line string")
    genome_id, genome = read_fasta(root / "analysis_genome.fasta")
    proteins = read_json(root / "evidence.json", list)
    with (root / "annotation.tsv").open(newline="") as handle:
        current_annotations = {row["protein_id"]: row for row in csv.DictReader(handle, delimiter="\t")}
    by_model = {}
    for protein in proteins:
        key = (protein["start"], protein["end"], protein["strand"])
        by_model.setdefault(key, []).append(protein)
    accepted, rejected = [], []
    external_records = list(SeqIO.parse(str(genbank), "genbank"))
    if not external_records: raise ValueError("External annotation contains no GenBank records")
    for record in external_records:
        try:
            same_genome = str(record.seq).upper() == genome.upper()
        except Exception as exc:
            raise ValueError("External GenBank must contain the nucleotide sequence for identity validation") from exc
        if not same_genome:
            raise ValueError(f"External sequence {record.id} differs from the analysis sequence")
        for feature in record.features:
            if feature.type != "CDS": continue
            qualifier = feature.qualifiers
            ident = (qualifier.get("locus_tag") or qualifier.get("protein_id") or ["unknown"])[0]
            if not feature.location or len(feature.location.parts) != 1:
                rejected.append({"external_id": ident, "reason": "Compound or missing location requires manual alignment"}); continue
            strand = "+" if feature.location.strand == 1 else "-" if feature.location.strand == -1 else "?"
            key = (int(feature.location.start) + 1, int(feature.location.end), strand)
            candidates = by_model.get(key, [])
            sequence = str((qualifier.get("translation") or [""])[0]).replace(" ", "").rstrip("*")
            if len(candidates) != 1 or not sequence or sequence != candidates[0]["sequence"].rstrip("*"):
                rejected.append({"external_id": ident, "start": key[0], "end": key[1],
                                 "reason": "No unique exact coordinate and protein-sequence match"}); continue
            protein = candidates[0]
            product = str((qualifier.get("product") or [""])[0]).strip()
            if product.lower() in UNKNOWN: continue
            if any(ord(char) < 32 for char in product):
                raise ValueError("External product contains control characters")
            current = current_annotations.get(protein["protein_id"], {}).get("product") or "hypothetical protein"
            from .benchmark import classify_product_relation
            relation = classify_product_relation(current, product)
            accepted.append({"protein_id": protein["protein_id"], "genome_id": genome_id,
                             "start": key[0], "end": key[1], "strand": key[2],
                             "current_product": current, "proposed_product": product,
                             "source": source, "source_version": version,
                             "database_version": database_version, "external_id": ident,
                             "relation": relation, "evidence_type": "EXTERNAL_COMPUTATIONAL_ANNOTATION",
                             "review_status": "PENDING_REVIEW", "exact_sequence_match": True,
                             "source_qualifiers": qualifier})
    payload = {"source": source, "source_version": version, "database_version": database_version,
               "source_genbank_sha256": checksum(genbank), "proposal_count": len(accepted),
               "rejected": rejected, "records": accepted,
               "automatic_product_changes": False, "experimentally_validated": False}
    with derived_run(root, output, "external_annotation_import") as destination:
        if (destination / "external_evidence.json").exists():
            raise ValueError("External evidence already present; start from the un-enriched run")
        directory = destination / "external_evidence_sources"; directory.mkdir()
        shutil.copyfile(genbank, directory / "annotations.gbk")
        (destination / "external_evidence.json").write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
        columns = ["protein_id", "current_product", "proposed_product", "source", "source_version",
                   "database_version", "external_id", "relation", "review_status"]
        with (destination / "annotation_proposals.tsv").open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", extrasaction="ignore")
            writer.writeheader(); writer.writerows(accepted)
        from .conflict_review import write_conflicts
        write_conflicts(destination)
        record_sidecar(destination, "external_evidence", {key: value for key, value in payload.items() if key != "records"})
    return {"status": "IMPORTED_FOR_REVIEW", "proposal_count": len(accepted), "rejected_count": len(rejected), "output": str(output)}
