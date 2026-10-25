"""NCBI requirement review: one test per rule, plus the package it produces."""
from __future__ import annotations

import json
from datetime import date
from pathlib import Path

import pytest

from phagemine import submission_review as sr
from phagemine.genbank import fasta_defline, feature_table, submitted_translation, write_package
from phagemine.genome import predict_orfs
from phagemine.io import read_fasta
from phagemine.models import SubmissionMetadata

ROOT = Path(__file__).resolve().parents[1]
TODAY = date(2026, 10, 25)


def codes(findings):
    return {f["code"] for f in findings}


@pytest.mark.parametrize("value", ["14-Mar-2025", "Mar-2025", "2025", "2025-03-14", "2025-03", "missing: not collected", "not provided"])
def test_accepted_collection_dates(value):
    assert sr.check_collection_date(value, TODAY) == []


@pytest.mark.parametrize("value", ["14/03/2025", "March 2025", "2025-13-01", "31-Feb-2025", "14-Mar-2027", "14-mar-2025"])
def test_rejected_collection_dates(value):
    assert codes(sr.check_collection_date(value, TODAY)) == {"BAD_COLLECTION_DATE"}


def test_geo_loc_name_uses_insdc_list_and_colon():
    assert sr.check_geo_loc_name("Nigeria: Plateau State, Jos") == []
    assert sr.check_geo_loc_name("Viet Nam") == []
    assert codes(sr.check_geo_loc_name("Nigeria, Jos")) == {"BAD_GEO_LOC_NAME"}
    bad = sr.check_geo_loc_name("nigeria: Jos")
    assert codes(bad) == {"BAD_GEO_LOC_NAME"} and "Nigeria" in bad[0]["message"]
    assert codes(sr.check_geo_loc_name("Vietnam")) == {"BAD_GEO_LOC_NAME"}


def test_metadata_rules():
    good = SubmissionMetadata(organism="Xanthomonas phage vB_XcaM_JOS1", submitter_email="a@b.org", locus_tag_prefix="JOS1")
    assert codes(sr.check_metadata(good)) == {"LOCAL_SUBMISSION_TEMPLATE"}
    bad = SubmissionMetadata(organism="my new phage", submitter_email="not-an-email", locus_tag_prefix="1X",
                             host="E. coli [K12]")
    assert {"PHAGE_ORGANISM_NAME", "BAD_SUBMITTER_EMAIL", "BAD_LOCUS_TAG_PREFIX", "BRACKET_IN_SOURCE_VALUE"} <= codes(sr.check_metadata(bad))
    assert "NO_LOCUS_TAG_PREFIX" in codes(sr.check_metadata(SubmissionMetadata()))


def test_sequence_rules():
    good = "ACGT" * 100
    assert sr.check_sequence(good) == []
    assert "SHORT_SEQUENCE" in codes(sr.check_sequence("ACGT" * 10))
    assert "TERMINAL_NS" in codes(sr.check_sequence("NN" + good))
    assert "INTERNAL_N_RUN" in codes(sr.check_sequence(good + "N" * 120 + good))
    assert "HIGH_PERCENT_N" in codes(sr.check_sequence(good + "N" * 60 + good))
    assert "INVALID_CHARACTERS" in codes(sr.check_sequence(good + "X" + good))
    assert "ZERO_BASECOUNT" in codes(sr.check_sequence("ACA" * 100))


@pytest.mark.parametrize("name,expected", [
    ("hypothetical protein", set()),
    ("terminase large subunit", set()),
    ("putative tail fiber protein", set()),
    ("DUF1234 domain-containing protein", set()),
    ("conserved hypothetical protein", {"HYPOTHETICAL_VARIANT"}),
    ("hypothetical protein, conserved", {"HYPOTHETICAL_VARIANT", "PUNCTUATION"}),
    ("probable portal protein", {"QUALIFIER_WORDING"}),
    ("terminase homolog", {"EVOLUTIONARY_IMPLICATION"}),
    ("phrog_123 protein", {"DATABASE_ID_IN_PRODUCT"}),
    ("DUF1234", {"DATABASE_ID_IN_PRODUCT"}),
    ("PM_000001", {"DATABASE_ID_IN_PRODUCT"}),
    ("endolysin EC 3.2.1.17", {"EC_NUMBER_IN_PRODUCT"}),
    ("partial tail protein", {"SUSPICIOUS_PHRASE"}),
    ("major capsid protein.", {"PUNCTUATION"}),
    ("holin 50% identity", {"PUNCTUATION"}),
    ("Major capsid protein", {"CAPITALIZATION"}),
    ("RecA protein", set()),
    ("tail fibre protein", {"SPELLING"}),
    ("lysozyme gene", {"GENE_IN_PRODUCT"}),
    ("nuclease protein", {"REDUNDANT_PROTEIN"}),
    ("", {"EMPTY_PRODUCT"}),
])
def test_product_name_rules(name, expected):
    assert {code for code, _, _ in sr.check_product(name)} == expected


def demo():
    genome_id, genome = read_fasta(ROOT / "examples/demo_phage.fasta")
    return genome_id, genome, predict_orfs(genome_id, genome)


def test_feature_rules():
    genome_id, genome, proteins = demo()
    products = {p.protein_id: "hypothetical protein" for p in proteins}
    baseline = codes(sr.check_features(genome, proteins, products))
    assert baseline <= {"SHORT_PROTEIN", "CONTAINED_CDS"}
    inner = proteins[1]
    inner.partial_5prime = True  # not at a sequence end
    assert "INTERNAL_PARTIAL_CDS" in codes(sr.check_features(genome, proteins, products))
    inner.partial_5prime = False
    inner.cds = "CCC" + inner.cds[3:]
    assert "NO_START_CODON" in codes(sr.check_features(genome, proteins, products))


def test_feature_table_puts_locus_tags_on_genes_only():
    genome_id, genome, proteins = demo()
    table = feature_table(genome_id, proteins[:1], locus_tag_prefix="JOS1")
    lines = table.splitlines()
    assert lines[1].endswith("\tgene") and lines[2] == "\t\t\tlocus_tag\tJOS1_000001"
    assert lines[3].endswith("\tCDS") and "locus_tag" not in "\n".join(lines[3:])
    assert "\t\t\tprotein_id\tgnl|JOS1|JOS1_000001" in lines
    assert "cds_provenance" not in table and "note" not in table
    untagged = feature_table(genome_id, proteins[:1])
    assert "locus_tag" not in untagged and "gene" in untagged


def test_defline_carries_source_modifiers_and_table_11():
    metadata = SubmissionMetadata(organism="Xanthomonas phage vB_XcaM_JOS1", host="Xanthomonas campestris",
                                  collection_date="14-Mar-2025", geographic_location="Nigeria: Plateau State, Jos",
                                  isolation_source="tomato leaf lesion")
    line = fasta_defline("JOS1", metadata)
    for part in ("[organism=Xanthomonas phage vB_XcaM_JOS1]", "[host=Xanthomonas campestris]", "[collection_date=14-Mar-2025]",
                 "[geo_loc_name=Nigeria: Plateau State, Jos]", "[isolation_source=tomato leaf lesion]", "[gcode=11]"):
        assert part in line


def test_alternative_start_is_translated_as_methionine():
    _, _, proteins = demo()
    protein = proteins[0]
    protein.sequence = "V" + protein.sequence[1:]
    assert submitted_translation(protein).startswith("M")
    protein.partial_5prime = True
    assert submitted_translation(protein).startswith("V")


def test_package_reports_ncbi_findings_with_sources(tmp_path):
    genome_id, genome, proteins = demo()
    metadata = SubmissionMetadata(organism="my phage", collection_date="2025/03/14", geographic_location="Nigeria, Jos",
                                  locus_tag_prefix="JOS1")
    validation = write_package(tmp_path, genome_id, genome, proteins, {"input_sha256": "x"}, metadata=metadata,
                               feature_overrides={proteins[0].protein_id: {"product": "probable holin."}})
    pre = validation["phagemine_pre_submission_validation"]
    error_codes = codes(pre["errors"])
    assert {"BAD_COLLECTION_DATE", "BAD_GEO_LOC_NAME", "QUALIFIER_WORDING", "PUNCTUATION"} <= error_codes
    assert validation["submission_readiness"] == "VALIDATION_FAILED"
    assert all(f["ncbi_reference"].startswith("https://") for f in pre["errors"] if "ncbi_check" in f)
    report = (tmp_path / "genbank_submission" / "NCBI_VALIDATION_REPORT.txt").read_text()
    assert "SEQ_DESCR_BadCollectionDate" in report and "Fix:" in report


def test_official_template_replaces_local_one(tmp_path):
    official = tmp_path / "official.sbt"
    official.write_text("Submit-block ::= { contact { } }\n")
    metadata = SubmissionMetadata(organism="Escherichia phage JOS2", locus_tag_prefix="JOS2", submission_template=str(official))
    assert "LOCAL_SUBMISSION_TEMPLATE" not in codes(sr.check_metadata(metadata))
    genome_id, genome, proteins = demo()
    write_package(tmp_path / "out", genome_id, genome, proteins, {"input_sha256": "x"}, metadata=metadata)
    assert (tmp_path / "out" / "genbank_submission" / "submission.sbt").read_text() == official.read_text()
    missing = SubmissionMetadata(submission_template=str(tmp_path / "nope.sbt"))
    assert "SUBMISSION_TEMPLATE_NOT_FOUND" in codes(sr.check_metadata(missing))
