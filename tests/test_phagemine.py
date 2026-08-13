import json
import os
import stat
import tempfile
import unittest
from pathlib import Path

from phagemine.genome import predict_orfs, translate
from phagemine.gene_prediction import DemoORFPredictor, PHANOTATEPredictor
from phagemine.io import read_fasta
from phagemine.models import Evidence, EvidenceLevel
from phagemine.pipeline import run
from phagemine.genbank import feature_table, product_name, table2asn_status, validate, write_package
from phagemine.annotation import MockEvidenceBackend
from phagemine.cli import main
from phagemine.models import SubmissionMetadata
from phagemine.genome_representation import GenomeRepresentation, Orientation, Rotation, Topology
from phagemine.sequencing_provenance import SequencingPlatform, SequencingProvenance


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
            count = run(ROOT / "examples/demo_phage.fasta", output, predictor=DemoORFPredictor())
            self.assertEqual(count, 5)
            for filename in ("original_input.fasta", "analysis_genome.fasta", "genome_representation.json", "sequencing_provenance.json", "genes.gff3", "proteins.faa", "cds.fna", "annotation.tsv", "evidence.json", "candidate_ranking.tsv", "report.md", "report.html", "run_manifest.json", "quality_control.json"):
                self.assertTrue((output / filename).exists(), filename)
            report = (output / "report.md").read_text()
            self.assertIn("Computational hypothesis only", report)
            evidence = json.loads((output / "evidence.json").read_text())
            self.assertTrue(any(item["missing_evidence"] for item in evidence))
            for filename in ("genome.fsa", "features.tbl", "proteins.faa", "sequencing_provenance.json", "validation.json", "provenance.json", "README.txt"):
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
            self.assertEqual(main(["genbank", str(ROOT / "examples/demo_phage.fasta"), "--gene-predictor", "demo", "--output", str(output)]), 0)
            self.assertTrue((output / "genbank_submission" / "validation.json").exists())

    def test_phanotate_adapter_parses_coordinates_strands_and_translation(self):
        genome_id, genome = read_fasta(ROOT / "examples/demo_phage.fasta")
        proteins = PHANOTATEPredictor.parse_output(genome_id, genome, "# START STOP FRAME CONTIG\n1 150 + demo_phage\n311 475 + demo_phage\n")
        self.assertEqual([protein.protein_id for protein in proteins], ["PM_000001", "PM_000002"])
        self.assertEqual(proteins[0].gene_call_source, "PHANOTATE")
        self.assertEqual(proteins[0].locus_tag, "PM_000001")
        self.assertEqual(proteins[0].start_codon, "ATG")
        self.assertEqual(proteins[0].stop_codon, "TAA")
        self.assertEqual(proteins[0].sequence, translate(genome[:150]))
        self.assertEqual(proteins[1].gene_call_parameters["reported_strand"], "+")

    def test_phanotate_adapter_handles_reverse_strand(self):
        sequence = "TTATTTCAT"
        proteins = PHANOTATEPredictor.parse_output("reverse", sequence, "9 1 - reverse\n")
        self.assertEqual(proteins[0].strand, "-")
        self.assertEqual(proteins[0].cds, "ATGAAATAA")
        self.assertEqual(proteins[0].sequence, "MK")

    def test_phanotate_unavailable_gives_actionable_error(self):
        predictor = PHANOTATEPredictor("/definitely/not/a/phanotate")
        with self.assertRaisesRegex(RuntimeError, "PHANOTATE is required"):
            predictor.predict("demo", "ATGAAATAA", ROOT / "examples/demo_phage.fasta")

    def test_original_representation_preserves_authoritative_sequence(self):
        representation = GenomeRepresentation.original("assembly", "ATGCCCTAA")
        self.assertEqual(representation.original_sequence, "ATGCCCTAA")
        self.assertEqual(representation.analysis_sequence, "ATGCCCTAA")
        self.assertEqual(representation.topology, Topology.UNKNOWN)
        self.assertEqual(representation.orientation, Orientation.ORIGINAL)
        self.assertEqual(representation.rotation, Rotation.NONE)
        self.assertEqual(representation.transform_history, ())

    def test_reverse_complement_representation_has_explicit_history(self):
        representation = GenomeRepresentation.original("assembly", "ATGCCCTAA").with_reverse_complement("User requested comparison to reference orientation", [{"type": "reference_alignment", "status": "supporting"}], {"accession": "REF_1"})
        self.assertEqual(representation.original_sequence, "ATGCCCTAA")
        self.assertEqual(representation.analysis_sequence, "TTAGGGCAT")
        self.assertEqual(representation.orientation, Orientation.REVERSE_COMPLEMENT)
        self.assertEqual(representation.transform_history[0].operation, "reverse_complement")
        self.assertEqual(representation.manifest()["reference"]["accession"], "REF_1")

    def test_circular_rotation_representation_records_coordinate_scope(self):
        representation = GenomeRepresentation.original("assembly", "AAACCCGGG", topology=Topology.CIRCULAR).with_rotation(4, "User selected coordinate 4 as analysis origin")
        self.assertEqual(representation.analysis_sequence, "CCCGGGAAA")
        self.assertEqual(representation.rotation, Rotation.ROTATED)
        self.assertEqual(representation.transform_history[0].parameters["analysis_origin"], 4)
        with self.assertRaises(ValueError):
            GenomeRepresentation.original("linear", "AAACCC").with_rotation(2, "not allowed")

    def test_combined_transform_history_and_analysis_coordinate_scope(self):
        representation = GenomeRepresentation.original("assembly", "ATGAAATAA", topology=Topology.CIRCULAR).with_reverse_complement("explicit test transform").with_rotation(2, "explicit test origin")
        self.assertEqual([event.operation for event in representation.transform_history], ["reverse_complement", "rotate"])
        self.assertEqual(representation.analysis_sequence, "TATTTCATT")
        self.assertIn("analysis_sequence_id", representation.manifest()["coordinate_scope"])

    def test_pipeline_preserves_original_fasta_and_coordinates_use_analysis_representation(self):
        genome_id, genome = read_fasta(ROOT / "examples/demo_phage.fasta")
        representation = GenomeRepresentation.original(genome_id, genome, topology=Topology.CIRCULAR).with_rotation(151, "Fixture-only explicit rotation")
        with tempfile.TemporaryDirectory() as temp:
            output = Path(temp) / "represented"
            run(ROOT / "examples/demo_phage.fasta", output, predictor=DemoORFPredictor(), representation=representation)
            self.assertEqual((output / "original_input.fasta").read_bytes(), (ROOT / "examples/demo_phage.fasta").read_bytes())
            manifest = json.loads((output / "run_manifest.json").read_text())
            self.assertEqual(manifest["genome_representation"]["analysis_sequence_id"], representation.analysis_sequence_id)
            gff = (output / "genes.gff3").read_text()
            self.assertIn(representation.analysis_sequence_id + "\tPhageMine\tCDS", gff)
            self.assertIn("Analysis sequence", (output / "report.md").read_text())
            annotations = (output / "annotation.tsv").read_text()
            self.assertIn("analysis_sequence_id", annotations)
            self.assertIn(representation.analysis_sequence_id, annotations)

    def test_illumina_provenance(self):
        provenance = SequencingProvenance.from_dict({"sequencing_platform": "ILLUMINA", "assembler": "SPAdes"})
        self.assertEqual(provenance.sequencing_platform, SequencingPlatform.ILLUMINA)
        self.assertEqual(provenance.assembler, "SPAdes")

    def test_oxford_nanopore_provenance(self):
        provenance = SequencingProvenance.from_dict({"sequencing_platform": "OXFORD_NANOPORE", "polishing_method": "Medaka"})
        self.assertEqual(provenance.sequencing_platform, SequencingPlatform.OXFORD_NANOPORE)
        self.assertEqual(provenance.polishing_method, "Medaka")

    def test_pacbio_provenance(self):
        provenance = SequencingProvenance.from_dict({"sequencing_platform": "PACBIO", "read_type": "HiFi"})
        self.assertEqual(provenance.sequencing_platform, SequencingPlatform.PACBIO)
        self.assertEqual(provenance.read_type, "HiFi")

    def test_hybrid_provenance(self):
        provenance = SequencingProvenance.from_dict({"sequencing_platform": "HYBRID", "assembly_method": "hybrid de novo"})
        self.assertEqual(provenance.sequencing_platform, SequencingPlatform.HYBRID)
        self.assertEqual(provenance.assembly_method, "hybrid de novo")

    def test_unknown_provenance_does_not_invent_metadata(self):
        provenance = SequencingProvenance()
        self.assertEqual(provenance.sequencing_platform, SequencingPlatform.UNKNOWN)
        self.assertIsNone(provenance.assembler)
        self.assertIsNone(provenance.raw_reads_available)

    def test_fasta_analysis_without_sequencing_metadata_uses_unknown(self):
        with tempfile.TemporaryDirectory() as temp:
            output = Path(temp) / "unknown-provenance"
            run(ROOT / "examples/demo_phage.fasta", output, predictor=DemoORFPredictor())
            manifest = json.loads((output / "run_manifest.json").read_text())
            self.assertEqual(manifest["sequencing_provenance"]["sequencing_platform"], "UNKNOWN")

    def test_sequencing_provenance_is_preserved_and_separate_from_representation(self):
        genome_id, genome = read_fasta(ROOT / "examples/demo_phage.fasta")
        representation = GenomeRepresentation.original(genome_id, genome)
        provenance = SequencingProvenance.from_dict({"sequencing_platform": "OXFORD_NANOPORE", "assembler": "Flye", "raw_reads_available": True, "metadata_source": "lab notebook"})
        with tempfile.TemporaryDirectory() as temp:
            output = Path(temp) / "ont"
            run(ROOT / "examples/demo_phage.fasta", output, predictor=DemoORFPredictor(), representation=representation, sequencing_provenance=provenance)
            manifest = json.loads((output / "run_manifest.json").read_text())
            self.assertEqual(manifest["sequencing_provenance"]["sequencing_platform"], "OXFORD_NANOPORE")
            self.assertEqual(manifest["genome_representation"]["topology"], "UNKNOWN")
            self.assertEqual(manifest["genome_representation"]["orientation"], "ORIGINAL")
            package = json.loads((output / "genbank_submission" / "sequencing_provenance.json").read_text())
            self.assertEqual(package["assembler"], "Flye")

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
