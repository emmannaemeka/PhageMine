"""Local GenBank pre-submission packaging with optional official table2asn validation.

No function in this module submits data to NCBI.
"""
from __future__ import annotations

import json
import shutil
import subprocess
from dataclasses import asdict
from pathlib import Path

from .models import EvidenceLevel, Protein, SubmissionMetadata
from .sequencing_provenance import SequencingProvenance


def _reverse_complement(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


def _safe(value: str | None) -> str:
    return (value or "").replace('"', "'").replace("\n", " ")


def validate(genome_id: str, genome: str, proteins: list[Protein], provenance: dict) -> dict:
    """PhageMine pre-submission checks, explicitly not official NCBI validation."""
    errors: list[dict] = []
    warnings: list[dict] = []
    seen_ids: set[str] = set()
    for protein in proteins:
        prefix = {"protein_id": protein.protein_id, "start": protein.start, "end": protein.end}
        if protein.protein_id in seen_ids:
            errors.append({**prefix, "code": "duplicate_protein_id", "message": "Protein identifiers must be unique."})
        seen_ids.add(protein.protein_id)
        if not (1 <= protein.start <= protein.end <= len(genome)):
            errors.append({**prefix, "code": "coordinate_out_of_range", "message": "CDS coordinates are outside the submitted FASTA sequence."})
            continue
        observed_cds = genome[protein.start - 1:protein.end]
        if protein.strand not in {"+", "-"}:
            errors.append({**prefix, "code": "invalid_strand", "message": "Strand must be + or -."})
        elif protein.strand == "-":
            observed_cds = _reverse_complement(observed_cds)
        if observed_cds != protein.cds:
            errors.append({**prefix, "code": "cds_fasta_mismatch", "message": "CDS sequence does not match its stated coordinates in submitted FASTA."})
        if len(protein.cds) % 3 != 0 and not (protein.partial_5prime or protein.partial_3prime):
            errors.append({**prefix, "code": "NCBI_CDS_FRAME", "severity": "ERROR", "message": "Complete CDS length is not divisible by three."})
        if not protein.sequence or "*" in protein.sequence:
            errors.append({**prefix, "code": "invalid_translation", "message": "Protein translation is empty or contains an internal stop."})
        expected_translation = __import__("phagemine.genome", fromlist=["translate"]).translate(protein.cds)
        if protein.sequence != expected_translation:
            errors.append({**prefix, "code": "translation_mismatch", "message": "Protein translation does not match the submitted CDS sequence."})
        if protein.strand in {"+", "-"} and len(protein.cds) >= 6:
            if protein.cds[:3] not in {"ATG", "GTG", "TTG"}:
                warnings.append({**prefix, "code": "noncanonical_start_codon", "message": "CDS does not begin with a standard PHANOTATE start codon."})
            if protein.cds[-3:] not in {"TAA", "TAG", "TGA"}:
                diagnostic = {**prefix, "code": "NCBI_CDS_NOSTOP", "message": "CDS does not end with a valid terminal stop codon.", "evidence": {"terminal_codon": protein.cds[-3:], "strand": protein.strand}}
                if protein.partial_3prime:
                    warnings.append({**diagnostic, "severity": "WARNING", "recommended_action": "Confirm that 3-prime partialness is biologically justified."})
                elif (protein.strand == "+" and protein.end == len(genome)) or (protein.strand == "-" and protein.start == 1):
                    errors.append({**diagnostic, "severity": "ERROR", "recommended_action": "Boundary CDS is marked complete. Review terminal-repeat/origin context before marking partial or changing coordinates."})
                    warnings.append({**prefix, "code": "POSSIBLE_ORIGIN_CROSSING_CDS", "severity": "REQUIRES_REVIEW", "message": "CDS touches a genome boundary; inspect DTR/circularization/origin context before changing the annotation."})
                else:
                    errors.append({**diagnostic, "severity": "ERROR", "recommended_action": "Correct the CDS coordinates/partial state; do not add a stop codon to the genome merely to satisfy validation."})
        if "mock" in protein.annotation.lower() or any(e.status == "mock" for e in protein.evidence):
            warnings.append({**prefix, "code": "mock_or_unsupported_function", "message": "Mock/computational annotation cannot support a GenBank functional claim; feature product will be emitted as hypothetical protein."})
        if protein.annotation_level.value in {"weak inference", "hypothesis requiring experimental validation"}:
            warnings.append({**prefix, "code": "uncharacterized_product", "message": "No supported functional product name is available; feature product will be hypothetical protein."})
    if provenance.get("input_sha256") is None:
        errors.append({"code": "missing_provenance", "message": "Input checksum is required for provenance linkage."})
    return {"validator": "PhageMine NCBI pre-submission validator", "official_ncbi_validation": False, "genome_id": genome_id, "sequence_length": len(genome), "coordinate_system": "1-based-inclusive", "valid": not errors, "ncbi_submission_ready": not errors, "fatal_error_count": len(errors), "warning_count": len(warnings), "errors": errors, "warnings": warnings, "provenance": provenance}


def product_name(protein: Protein) -> str:
    """Only experimentally/curated, non-mock claims may become a submitted product name."""
    if protein.annotation_level in {EvidenceLevel.EXPERIMENTAL, EvidenceLevel.CURATED} and "mock" not in protein.annotation.lower() and not any(e.status == "mock" for e in protein.evidence):
        return protein.annotation
    return "hypothetical protein"


def locus_tag(prefix: str | None, protein_id: str) -> str | None:
    """NCBI form PREFIX_uniqueID; only with a registered prefix (no prefix, no locus_tag)."""
    if not prefix:
        return None
    return f"{prefix}_{protein_id.rsplit('_', 1)[-1]}"


def protein_fasta_id(prefix: str | None, protein_id: str) -> str:
    tag = locus_tag(prefix, protein_id)
    return f"gnl|{prefix}|{tag}" if tag else f"gnl|PhageMine|{protein_id}"


def submitted_product(protein: Protein, overrides: dict | None = None) -> str:
    return (overrides or {}).get(protein.protein_id, {}).get("product") or product_name(protein)


def feature_table(genome_id: str, proteins: list[Protein], overrides: dict | None = None, locus_tag_prefix: str | None = None) -> str:
    """Five-column table: a gene feature (carrying the locus_tag) and a CDS for every protein.

    NCBI expects locus_tags on gene features only (DISC: MISSING_GENES, NON_GENE_LOCUS_TAG).
    Notes are written only when a researcher supplied one; local file references never are.
    """
    lines = [f">Feature {genome_id}"]
    for protein in proteins:
        start, end = (protein.start, protein.end) if protein.strand == "+" else (protein.end, protein.start)
        # NCBI 5-column feature tables use < and > on the biological 5-prime/3-prime ends.
        if protein.partial_5prime:
            start = f"<{start}" if protein.strand == "+" else f">{start}"
        if protein.partial_3prime:
            end = f">{end}" if protein.strand == "+" else f"<{end}"
        change = (overrides or {}).get(protein.protein_id, {})
        tag = change.get("locus_tag") or locus_tag(locus_tag_prefix, protein.protein_id)
        lines.append(f"{start}\t{end}\tgene")
        if tag:
            lines.append(f"\t\t\tlocus_tag\t{tag}")
        lines.extend([f"{start}\t{end}\tCDS", f"\t\t\tproduct\t{_safe(submitted_product(protein, overrides))}",
                      f"\t\t\tprotein_id\t{protein_fasta_id(locus_tag_prefix, protein.protein_id)}"])
        if change.get("note"):
            lines.append(f"\t\t\tnote\t{_safe(change['note'])}")
    return "\n".join(lines) + "\n"


SOURCE_MODIFIERS = (("organism", "organism"), ("host", "host"), ("isolation_source", "isolation_source"),
                    ("collection_date", "collection_date"), ("geographic_location", "geo_loc_name"),
                    ("isolate", "isolate"), ("strain", "strain"))


def fasta_defline(genome_id: str, metadata: SubmissionMetadata) -> str:
    """Source information goes on the FASTA definition line, where table2asn reads it.
    Phages use translation table 11."""
    parts = [f">{genome_id}"]
    for field, modifier in SOURCE_MODIFIERS:
        value = getattr(metadata, field, None)
        if value:
            parts.append(f"[{modifier}={_safe(str(value))}]")
    parts += ["[mol_type=genomic DNA]", "[gcode=11]"]
    return " ".join(parts)


def submitted_translation(protein: Protein) -> str:
    """NCBI translates an alternative start codon (GTG, TTG...) of a complete CDS as M."""
    sequence = protein.sequence or ""
    if sequence and not protein.partial_5prime and sequence[0] != "M":
        return "M" + sequence[1:]
    return sequence


def table2asn_status(package_root: Path, executable: str | None = None) -> dict:
    executable = executable or shutil.which("table2asn")
    if not executable or (executable != "table2asn" and not Path(executable).exists()):
        return {"performed": False, "state": "unavailable", "official_ncbi_validation": False, "message": "NCBI table2asn was not found; official NCBI validation has not been performed.", "installation": "Install NCBI table2asn from the NCBI Submission Portal tools documentation, add it to PATH, then rerun phagemine genbank."}
    output = package_root / "table2asn_output"
    output.mkdir(exist_ok=True)
    command = [executable, "-indir", str(package_root), "-outdir", str(output), "-t", str(package_root / "submission.sbt"), "-M", "n", "-Z", str(output / "table2asn.val")]
    completed = subprocess.run(command, capture_output=True, text=True, check=False)
    files = sorted(str(path.relative_to(package_root)) for pattern in ("*.val", "*.stats", "*.gbf", "*.sqn") for path in output.glob(pattern))
    return {"performed": True, "state": "passed" if completed.returncode == 0 else "failed", "official_ncbi_validation": completed.returncode == 0, "command": command, "returncode": completed.returncode, "stdout": completed.stdout, "stderr": completed.stderr, "output_files": files, "message": "table2asn executed locally; inspect captured outputs before final NCBI submission."}


def _submission_template(metadata: SubmissionMetadata) -> str:
    # A minimal, inspectable template input. Users still need to review it with NCBI tooling.
    return "Submit-block ::= {\n  contact { contact { name name { last \"%s\" }, email \"%s\" } },\n  cit { authors { names std { %s } } }\n}\n" % (_safe(metadata.submitter_name), _safe(metadata.submitter_email), ", ".join('{ name name { last \"%s\" } }' % _safe(author) for author in metadata.authors))


def readiness(pre_validation: dict, metadata: SubmissionMetadata, table2asn: dict) -> tuple[str, list[str]]:
    missing = metadata.missing_required()
    if not pre_validation["valid"]:
        return "VALIDATION_FAILED", missing
    if missing:
        return "INCOMPLETE_METADATA", missing
    if table2asn.get("state") == "failed":
        return "VALIDATION_FAILED", missing
    if pre_validation["warnings"] or table2asn.get("state") != "passed":
        return "READY_WITH_WARNINGS", missing
    return "READY", missing


def write_package(output: str | Path, genome_id: str, genome: str, proteins: list[Protein], provenance: dict, metadata: SubmissionMetadata | None = None, table2asn_executable: str | None = None, sequencing_provenance: SequencingProvenance | None = None, feature_overrides: dict | None = None) -> dict:
    root = Path(output) / "genbank_submission"
    root.mkdir(parents=True, exist_ok=True)
    metadata = metadata or SubmissionMetadata()
    sequencing_provenance = sequencing_provenance or SequencingProvenance()
    pre_validation = validate(genome_id, genome, proteins, provenance)
    if metadata.sequence_id and metadata.sequence_id != genome_id:
        pre_validation["errors"].append({"code": "metadata_sequence_id_mismatch", "message": "Metadata sequence_id does not match the FASTA identifier."})
        pre_validation["valid"] = False
    prefix = getattr(metadata, "locus_tag_prefix", None)
    from .submission_review import review
    products = {p.protein_id: submitted_product(p, feature_overrides) for p in proteins}
    ncbi_review = review(genome, proteins, products, metadata)
    pre_validation["errors"].extend(ncbi_review["errors"])
    pre_validation["warnings"].extend(ncbi_review["warnings"])
    pre_validation["valid"] = pre_validation["ncbi_submission_ready"] = not pre_validation["errors"]
    pre_validation["fatal_error_count"], pre_validation["warning_count"] = len(pre_validation["errors"]), len(pre_validation["warnings"])
    pre_validation["ncbi_requirements_review"] = {"rules_source": "phagemine.submission_review (NCBI/INSDC guidance, October 2026)",
                                                  "errors": len(ncbi_review["errors"]), "warnings": len(ncbi_review["warnings"])}
    (root / "genome.fsa").write_text(f"{fasta_defline(genome_id, metadata)}\n{genome}\n")
    (root / "features.tbl").write_text(feature_table(genome_id, proteins, feature_overrides, prefix))
    (root / "proteins.faa").write_text("".join(f">{protein_fasta_id(prefix, p.protein_id)} {p.protein_id}\n{submitted_translation(p)}\n" for p in proteins))
    (root / "submission_metadata.json").write_text(json.dumps(asdict(metadata), indent=2, sort_keys=True))
    (root / "sequencing_provenance.json").write_text(json.dumps(sequencing_provenance.manifest(), indent=2, sort_keys=True))
    official = Path(metadata.submission_template).expanduser() if metadata.submission_template else None
    if official and official.is_file():
        shutil.copyfile(official, root / "submission.sbt")
    else:
        (root / "submission.sbt").write_text(_submission_template(metadata))
    cds_provenance = {p.protein_id: {"coordinates": {"start": p.start, "end": p.end, "strand": p.strand}, "coordinate_system": "1-based-inclusive", "gene_call_parameters": p.gene_call_parameters, "protein_fasta_id": protein_fasta_id(prefix, p.protein_id), "evidence_record": f"../evidence.json#{p.protein_id}", "evidence": [asdict(e) for e in p.evidence]} for p in proteins}
    (root / "cds_provenance.json").write_text(json.dumps(cds_provenance, indent=2, default=str, sort_keys=True))
    table2asn = table2asn_status(root, table2asn_executable)
    state, missing = readiness(pre_validation, metadata, table2asn)
    validation = {"submission_readiness": state, "missing_required_metadata": missing, "phagemine_pre_submission_validation": pre_validation, "ncbi_table2asn_validation": table2asn, "final_ncbi_submission": {"performed": False, "message": "No submission to NCBI was attempted. Final NCBI review occurs only after a user submits through NCBI."}}
    (root / "validation.json").write_text(json.dumps(validation, indent=2, sort_keys=True))
    report_lines = [
        "PhageMine NCBI Pre-submission Validation",
        "========================================",
        f"Genome: {genome_id}",
        f"Length: {len(genome):,} bp",
        f"CDSs: {len(proteins)}",
        "",
        f"OVERALL: {'READY FOR NCBI SUBMISSION' if state == 'READY' else 'NOT READY FOR NCBI SUBMISSION'}",
        f"Fatal errors: {len(pre_validation['errors'])}",
        f"Warnings/review items: {len(pre_validation['warnings'])}",
    ]
    for item in pre_validation["errors"] + pre_validation["warnings"]:
        report_lines.extend(["", f"{item.get('severity', 'ERROR' if item in pre_validation['errors'] else 'WARNING')} {item['code']}"
                             + (f" (NCBI: {item['ncbi_check']})" if item.get("ncbi_check") else ""),
                             f"{item.get('protein_id', '')} {item.get('message', '')}".strip()]
                            + ([f"Fix: {item['recommended_action']}"] if item.get("recommended_action") else [])
                            + ([f"Source: {item['ncbi_reference']}"] if item.get("ncbi_reference") else []))
    (root / "NCBI_VALIDATION_REPORT.txt").write_text("\n".join(report_lines) + "\n")
    (root / "provenance.json").write_text(json.dumps(provenance, indent=2, sort_keys=True))
    (root / "README.txt").write_text("This is a local PhageMine pre-submission package; it has not been submitted to NCBI. PhageMine QC, PhageMine pre-submission validation, table2asn validation, and final NCBI review are distinct stages. If table2asn is unavailable, install the official NCBI table2asn distribution, add it to PATH, and rerun this command. Never treat mock, weak, or hypothesis-level annotations as asserted product names.\n")
    return validation
