import json
import os
import stat
import tempfile
import unittest
from pathlib import Path

from phagemine.genome import predict_orfs, translate
from phagemine.io import read_fasta
from phagemine.models import Evidence, EvidenceLevel
from phagemine.pipeline import run
from phagemine.genbank import feature_table, product_name, table2asn_status, validate, write_package
from phagemine.annotation import MockEvidenceBackend
from phagemine.cli import main
from phagemine.models import SubmissionMetadata


ROOT = Path(__file__).resolve().parents[1]


class PhageMineTests(unittest.TestCase):
    def test_translation(self):
        self.assertEqual(translate("ATGGCTTAA"), "MA")

    def test_fasta_validation(self):
        genome_id, sequence = read_fasta(ROOT / "examples/demo_phage.fasta")
        self.assertEqual(genome_id, "demo_phage")
        self.assertTrue(sequence.startswith("ATG"))

    def test_orfs_receive_stable_ids(self):
        _, sequence = read_fasta(ROOT / "examples/demo_phage.fasta")
        proteins = predict_orfs("demo", sequence)
        self.assertEqual(len(proteins), 5)
        self.assertEqual(proteins[0].protein_id, "PM_000001")

    def test_evidence_has_provenance(self):
        evidence = Evidence("domain", "signal", EvidenceLevel.WEAK, "source", "v1")
        self.assertEqual(evidence.source_version, "v1")
        self.assertTrue(evidence.supports)

    def test_full_pipeline_generates_explainable_outputs(self):
        with tempfile.TemporaryDirectory() as temp:
            output = Path(temp) / "run"
            count = run(ROOT / "examples/demo_phage.fasta", output)
            self.assertEqual(count, 5)
            for filename in ("genes.gff3", "proteins.faa", "cds.fna", "annotation.tsv", "evidence.json", "candidate_ranking.tsv", "report.md", "report.html", "run_manifest.json", "quality_control.json"):
                self.assertTrue((output / filename).exists(), filename)
            report = (output / "report.md").read_text()
            self.assertIn("Computational hypothesis only", report)
            evidence = json.loads((output / "evidence.json").read_text())
            self.assertTrue(any(item["missing_evidence"] for item in evidence))
            for filename in ("genome.fsa", "features.tbl", "proteins.faa", "validation.json", "provenance.json", "README.txt"):
                self.assertTrue((output / "genbank_submission" / filename).exists(), filename)

    def test_genbank_feature_table_and_validation(self):
        genome_id, genome = read_fasta(ROOT / "examples/demo_phage.fasta")
        proteins = predict_orfs(genome_id, genome)
        validation = validate(genome_id, genome, proteins, {"input_sha256": "fixture"})
        self.assertTrue(validation["valid"])
        table = feature_table(genome_id, proteins)
        self.assertIn(f"{proteins[0].start}\t{proteins[0].end}\tCDS", table)
        self.assertIn("gnl|PhageMine|PM_000001", table)
        self.assertIn("product\thypothetical protein", table)
        self.assertEqual(validation["provenance"]["input_sha256"], "fixture")

    def test_genbank_reverse_strand_coordinates(self):
        genome_id, genome = read_fasta(ROOT / "examples/demo_phage.fasta")
        protein = predict_orfs(genome_id, genome)[0]
        protein.strand = "-"
        protein.cds = protein.cds.translate(str.maketrans("ACGT", "TGCA"))[::-1]
        table = feature_table(genome_id, [protein])
        self.assertIn(f"{protein.end}\t{protein.start}\tCDS", table)
        self.assertTrue(validate(genome_id, genome, [protein], {"input_sha256": "fixture"})["valid"])

    def test_genbank_reports_fasta_and_claim_problems(self):
        genome_id, genome = read_fasta(ROOT / "examples/demo_phage.fasta")
        proteins = predict_orfs(genome_id, genome)
        proteins[0].cds = "ATGAAA"
        invalid = validate(genome_id, genome, proteins, {"input_sha256": "fixture"})
        self.assertFalse(invalid["valid"])
        self.assertIn("cds_fasta_mismatch", {item["code"] for item in invalid["errors"]})

    def test_genbank_warns_on_unsupported_mock_claims(self):
        genome_id, genome = read_fasta(ROOT / "examples/demo_phage.fasta")
        proteins = predict_orfs(genome_id, genome)
        MockEvidenceBackend().annotate(proteins)
        warnings = validate(genome_id, genome, proteins, {"input_sha256": "fixture"})["warnings"]
        self.assertIn("mock_or_unsupported_function", {item["code"] for item in warnings})

    def test_genbank_cli_is_independently_callable(self):
        with tempfile.TemporaryDirectory() as temp:
            output = Path(temp) / "genbank"
            self.assertEqual(main(["genbank", str(ROOT / "examples/demo_phage.fasta"), "--output", str(output)]), 0)
            self.assertTrue((output / "genbank_submission" / "validation.json").exists())

    def _complete_metadata(self):
        return SubmissionMetadata.from_dict(json.loads((ROOT / "examples/submission_metadata.json").read_text()))

    def _fake_table2asn(self, directory, exit_code=0):
        script = Path(directory) / f"fake_table2asn_{exit_code}.sh"
        script.write_text(f"#!/bin/sh\nmkdir -p \"$4\"\necho validation > \"$4/fake.val\"\necho stats > \"$4/fake.stats\"\nexit {exit_code}\n")
        script.chmod(script.stat().st_mode | stat.S_IXUSR)
        return str(script)

    def _submission_acceptable_proteins(self, genome_id, genome):
        proteins = predict_orfs(genome_id, genome)
        for protein in proteins:
            protein.annotation = "capsid protein"
            protein.annotation_level = EvidenceLevel.CURATED
        return proteins

    def test_missing_metadata_is_not_submission_ready(self):
        genome_id, genome = read_fasta(ROOT / "examples/demo_phage.fasta")
        with tempfile.TemporaryDirectory() as temp:
            result = write_package(temp, genome_id, genome, predict_orfs(genome_id, genome), {"input_sha256": "fixture"}, table2asn_executable="not-installed-table2asn")
            self.assertEqual(result["submission_readiness"], "INCOMPLETE_METADATA")
            self.assertFalse(result["ncbi_table2asn_validation"]["official_ncbi_validation"])

    def test_complete_metadata_and_table2asn_success_are_ready(self):
        genome_id, genome = read_fasta(ROOT / "examples/demo_phage.fasta")
        with tempfile.TemporaryDirectory() as temp:
            result = write_package(temp, genome_id, genome, self._submission_acceptable_proteins(genome_id, genome), {"input_sha256": "fixture"}, self._complete_metadata(), self._fake_table2asn(temp))
            self.assertEqual(result["submission_readiness"], "READY")
            self.assertEqual(result["ncbi_table2asn_validation"]["state"], "passed")
            self.assertIn("table2asn_output/fake.val", result["ncbi_table2asn_validation"]["output_files"])

    def test_table2asn_failure_fails_readiness(self):
        genome_id, genome = read_fasta(ROOT / "examples/demo_phage.fasta")
        with tempfile.TemporaryDirectory() as temp:
            result = write_package(temp, genome_id, genome, predict_orfs(genome_id, genome), {"input_sha256": "fixture"}, self._complete_metadata(), self._fake_table2asn(temp, 1))
            self.assertEqual(result["submission_readiness"], "VALIDATION_FAILED")
            self.assertEqual(result["ncbi_table2asn_validation"]["state"], "failed")

    def test_table2asn_unavailable_is_explicit(self):
        result = table2asn_status(Path(tempfile.mkdtemp()), "not-installed-table2asn")
        self.assertEqual(result["state"], "unavailable")
        self.assertFalse(result["official_ncbi_validation"])

    def test_conservative_product_and_cds_provenance(self):
        genome_id, genome = read_fasta(ROOT / "examples/demo_phage.fasta")
        proteins = predict_orfs(genome_id, genome)
        MockEvidenceBackend().annotate(proteins)
        with tempfile.TemporaryDirectory() as temp:
            write_package(temp, genome_id, genome, proteins, {"input_sha256": "fixture"})
            package = Path(temp) / "genbank_submission"
            self.assertEqual(product_name(proteins[0]), "hypothetical protein")
            self.assertNotIn("Putative structural", (package / "features.tbl").read_text())
            provenance = json.loads((package / "cds_provenance.json").read_text())
            self.assertTrue(provenance["PM_000001"]["evidence"])
