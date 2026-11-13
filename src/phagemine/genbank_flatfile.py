"""Write a local GenBank view of the submission package, without table2asn."""
from pathlib import Path
import json
import hashlib

from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, SimpleLocation, BeforePosition, AfterPosition


def write_flatfile(package: Path) -> Path:
    """Keep sequence, products, overrides and RNA features aligned with features.tbl."""
    package = Path(package)
    record = SeqIO.read(package / "genome.fsa", "fasta")
    metadata = json.loads((package / "submission_metadata.json").read_text())
    # Keep the full accession in ACCESSION/VERSION; bound only the LOCUS label.
    if len(record.name) > 16:
        record.name = record.name[:7] + "_" + hashlib.sha256(record.id.encode()).hexdigest()[:8]
    record.description = "PhageMine annotated genome"
    record.annotations = {"molecule_type": "DNA", "data_file_division": "PHG",
                          "comment": "Computational annotations; this file does not establish NCBI submission readiness."}
    topology = metadata.get("topology")
    representation = package.parent / "genome_representation.json"
    if topology is None and representation.is_file():
        topology = json.loads(representation.read_text()).get("topology", "unknown").lower()
    if topology in {"linear", "circular"}:
        record.annotations["topology"] = topology
    elif topology not in {None, "unknown"}:
        raise ValueError("topology must be linear, circular or unknown")
    else:
        record.annotations["comment"] += " Genome topology is unknown; no topology was inferred."
    source = {}
    for field, qualifier in (("organism", "organism"), ("host", "host"),
                             ("isolation_source", "isolation_source"), ("collection_date", "collection_date"),
                             ("geographic_location", "geo_loc_name"), ("isolate", "isolate"), ("strain", "strain")):
        if metadata.get(field):
            source[qualifier] = [metadata[field]]
    if metadata.get("organism"):
        record.annotations["organism"] = metadata["organism"]
    record.features = [SeqFeature(SimpleLocation(0, len(record)), type="source", qualifiers=source)]
    with (package / "proteins.faa").open() as handle:
        translations = {p.id: str(p.seq) for p in SeqIO.parse(handle, "fasta")}
    feature = None
    for line in (package / "features.tbl").read_text().splitlines():
        if not line or line.startswith(">"):
            continue
        fields = line.split("\t")
        if fields[0]:
            first, last, kind = fields[:3]
            a, b = int(first.lstrip("<>")), int(last.lstrip("<>"))
            strand = 1 if a <= b else -1
            low, high = (first, last) if strand == 1 else (last, first)
            start = BeforePosition(min(a, b) - 1) if low[0] in "<>" else min(a, b) - 1
            end = AfterPosition(max(a, b)) if high[0] in "<>" else max(a, b)
            feature = SeqFeature(SimpleLocation(start, end, strand=strand), type=kind)
            record.features.append(feature)
        elif feature is not None and len(fields) >= 5:
            feature.qualifiers.setdefault(fields[3], []).append(fields[4])
    for feature in record.features:
        if feature.type == "CDS":
            protein_id = feature.qualifiers.get("protein_id", [""])[0]
            if protein_id in translations:
                feature.qualifiers["translation"] = [translations[protein_id]]
            feature.qualifiers.setdefault("codon_start", ["1"])
            feature.qualifiers["transl_table"] = ["11"]
    destination = package / "genome.gbk"
    SeqIO.write(record, destination, "genbank")
    return destination
