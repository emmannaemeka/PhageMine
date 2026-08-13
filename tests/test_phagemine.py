import json
import os
import stat
import shutil
import tempfile
import unittest
from io import StringIO
from unittest.mock import patch
from pathlib import Path

from phagemine.genome import predict_orfs, translate
from phagemine.gene_prediction import DemoORFPredictor, PHANOTATEPredictor
from phagemine.io import read_fasta
from phagemine.models import Evidence, EvidenceLevel, Protein
from phagemine.pipeline import run
from phagemine.genbank import feature_table, product_name, table2asn_status, validate, write_package
from phagemine.annotation import MockEvidenceBackend
from phagemine.cli import main
from phagemine.models import SubmissionMetadata
from phagemine.genome_representation import GenomeRepresentation, Orientation, Rotation, Topology
from phagemine.sequencing_provenance import SequencingPlatform, SequencingProvenance
from phagemine.pfam import PfamHMMAdapter
from phagemine.vog import VOGHMMAdapter
from phagemine.phrogs import MMSEQS_FORMAT, PHROGSMMseqsAdapter
from phagemine.swissprot import SwissProtEvidenceAdapter
from phagemine.progress import ProgressReporter
from phagemine.resume import _load_source, RESUME_STAGES
from phagemine.fusion import classify_protein, classify_proteins, normalize_function, write_classification
from phagemine.evidence import EvidenceAdapterResult
from phagemine.mining import mine
from phagemine.resources import EvidenceResourceManager, ResourceStatus, ResourceType, default_registry_path


ROOT = Path(__file__).resolve().parents[1]


class PhageMineTests(unittest.TestCase):
    def _fusion_protein(self, evidence=()):
        return Protein("g", "P", 1, 30, "+", "ATG" * 10, "M" * 10, "test", evidence=list(evidence))

    def _fusion_e(self, source, desc, strength="STRONG", supports=True, modality="test", category=None, identifier="x"):
        metrics = {"functional_category": category} if category else {}
        return Evidence(modality, "statement", EvidenceLevel.COMPUTATIONAL, source, "1", status="REAL", supports=supports, identifier=identifier, evidence_strength=strength, description=desc, metrics=metrics)

    def test_fusion_classification_rules_and_determinism(self):
        unknown = classify_protein(self._fusion_protein([self._fusion_e("PHROGs", "unknown function", category="unknown")]))
        self.assertEqual(unknown["functional_state"], "CONSERVED_UNKNOWN")
        probable = classify_protein(self._fusion_protein([self._fusion_e("PHROGs", "RNA polymerase")]))
        self.assertEqual(probable["functional_state"], "PROBABLE_FUNCTION")
        orthology_pair = classify_protein(self._fusion_protein([self._fusion_e("PHROGs", "toxin", identifier="p"), self._fusion_e("VOGDB", "toxin", identifier="v")]))
        self.assertEqual(orthology_pair["functional_state"], "PROBABLE_FUNCTION")
        self.assertEqual(orthology_pair["supporting_source_count"], 2)
        self.assertEqual(orthology_pair["supporting_modality_count"], 2)
        known = classify_protein(self._fusion_protein([self._fusion_e("Swiss-Prot", "major capsid protein"), self._fusion_e("PHROGs", "Major capsid protein {ECO:0001}", identifier="y")]))
        self.assertEqual(known["functional_state"], "KNOWN_FUNCTION")
        broad = classify_protein(self._fusion_protein([self._fusion_e("PHROGs", "unknown function", category="head and packaging")]))
        self.assertEqual(broad["functional_state"], "FUNCTIONAL_CLASS_ONLY")
        unresolved = classify_protein(self._fusion_protein())
        self.assertEqual(unresolved["functional_state"], "UNRESOLVED")
        rejected = classify_protein(self._fusion_protein([self._fusion_e("Swiss-Prot", "DNA polymerase", supports=False)]))
        self.assertEqual(rejected["functional_state"], "UNRESOLVED")
        conflict = classify_protein(self._fusion_protein([self._fusion_e("Swiss-Prot", "integrase"), self._fusion_e("PHROGs", "major capsid protein", identifier="y")]))
        self.assertEqual(conflict["functional_state"], "CONFLICTING_EVIDENCE")
        same_source = classify_protein(self._fusion_protein([self._fusion_e("PHROGs", "integrase", identifier="a"), self._fusion_e("PHROGs", "integrase", identifier="b")]))
        self.assertEqual(same_source["supporting_source_count"], 1)
        self.assertEqual(same_source["supporting_modality_count"], 1)
        rejected_alt = classify_protein(self._fusion_protein([self._fusion_e("Swiss-Prot", "integrase"), self._fusion_e("PHROGs", "major capsid protein", supports=False, identifier="rejected")]))
        self.assertNotEqual(rejected_alt["functional_state"], "CONFLICTING_EVIDENCE")
        hypothetical = classify_protein(self._fusion_protein([self._fusion_e("VOGDB", "REFSEQ hypothetical protein")]))
        self.assertEqual(hypothetical["functional_state"], "CONSERVED_UNKNOWN")
        unavailable = classify_protein(self._fusion_protein())
        self.assertEqual(unavailable["functional_state"], "UNRESOLVED")
        original = self._fusion_e("PHROGs", "RNA polymerase")
        before = original.to_dict() if hasattr(original, "to_dict") else original.__dict__.copy()
        classify_protein(self._fusion_protein([original]))
        self.assertEqual(before, original.__dict__)
        self.assertEqual(normalize_function("Major capsid protein {ECO:0001}"), "major capsid protein")
        self.assertEqual(unknown, classify_protein(self._fusion_protein([self._fusion_e("PHROGs", "unknown function", category="unknown")])))

    def test_fusion_outputs_cover_integration_fixture(self):
        proteins = [self._fusion_protein() for _ in range(99)]
        for i, protein in enumerate(proteins, 1): protein.protein_id = f"P{i:03d}"
        with tempfile.TemporaryDirectory() as temp:
            write_classification(temp, proteins)
            self.assertEqual(len(json.loads((Path(temp) / "functional_classification.json").read_text())), 99)
            self.assertEqual(len((Path(temp) / "functional_classification.tsv").read_text().splitlines()), 100)
            self.assertEqual([x["protein_id"] for x in classify_proteins(proteins)], [x["protein_id"] for x in classify_proteins(proteins)])

    def test_fusion_conserved_unknown_from_phrogs_and_vogdb(self):
        result = classify_protein(self._fusion_protein([self._fusion_e("PHROGs", "unknown function"), self._fusion_e("VOGDB", "hypothetical protein")]))
        self.assertEqual(result["functional_state"], "CONSERVED_UNKNOWN")

    def test_fusion_phrogs_and_vogdb_same_function_is_probable(self):
        result = classify_protein(self._fusion_protein([self._fusion_e("PHROGs", "portal protein"), self._fusion_e("VOGDB", "portal protein")]))
        self.assertEqual(result["functional_state"], "PROBABLE_FUNCTION")
        self.assertNotEqual(result["functional_state"], "KNOWN_FUNCTION")

    def test_fusion_category_only(self):
        result = classify_protein(self._fusion_protein([self._fusion_e("PHROGs", "unknown function", category="head and packaging")]))
        self.assertEqual(result["functional_state"], "FUNCTIONAL_CLASS_ONLY")

    def test_fusion_supporting_source_and_modality_counts(self):
        result = classify_protein(self._fusion_protein([self._fusion_e("PHROGs", "portal protein"), self._fusion_e("PHROGs", "portal protein", identifier="2"), self._fusion_e("VOGDB", "portal protein")]))
        self.assertEqual(result["supporting_source_count"], 2)
        self.assertEqual(result["supporting_modality_count"], 2)

    def test_fusion_rejected_alternative_cannot_conflict(self):
        result = classify_protein(self._fusion_protein([self._fusion_e("Swiss-Prot", "integrase"), self._fusion_e("PHROGs", "major capsid protein", supports=False)]))
        self.assertNotEqual(result["functional_state"], "CONFLICTING_EVIDENCE")

    def test_fusion_unavailable_is_not_negative_evidence(self):
        result = classify_protein(self._fusion_protein([self._fusion_e("PHROGs", "portal protein")]))
        self.assertEqual(result["functional_state"], "PROBABLE_FUNCTION")

    def test_fusion_ambiguous_labels_are_not_conflict(self):
        result = classify_protein(self._fusion_protein([self._fusion_e("PHROGs", "DNA-dependent RNA polymerase"), self._fusion_e("VOGDB", "T7 RNA polymerase")]))
        self.assertNotEqual(result["functional_state"], "CONFLICTING_EVIDENCE")
        self.assertTrue(result["ambiguity_flags"])

    def test_fusion_resume_output_layer_regenerates_both_files(self):
        proteins = [self._fusion_protein()]
        with tempfile.TemporaryDirectory() as temp:
            write_classification(temp, proteins)
            self.assertTrue((Path(temp) / "functional_classification.json").is_file())
            self.assertTrue((Path(temp) / "functional_classification.tsv").is_file())
    def test_resume_progress_uses_resume_stage_ledger_without_duplicate_labels(self):
        stream = StringIO()
        progress = ProgressReporter(stream=stream)
        progress.STAGES = RESUME_STAGES
        progress.start("PHROGs")
        progress.finish("accepted hits")
        output = stream.getvalue()
        self.assertIn("[PhageMine 0/11] RUNNING: PHROGs", output)
        self.assertIn("[PhageMine 1/11] DONE: PHROGs", output)
        self.assertNotIn("RUNNING: RUNNING:", output)

    def test_resume_source_loader_reuses_proteins_and_evidence_verbatim(self):
        source = ROOT / "results" / "phage_c6_full_evidence"
        manifest, representation, proteins, sequence, _ = _load_source(source)
        self.assertEqual(len(proteins), 99)
        self.assertEqual(proteins[0].protein_id, "PM_000001")
        self.assertEqual(proteins[0].evidence[0].provenance["adapter"], "VOGHMMAdapter")
        self.assertEqual(manifest["evidence_adapters"][0]["provenance"]["threshold_mode"], "GA")
        self.assertEqual(representation.analysis_sequence, sequence)

    def test_resume_source_loader_rejects_protein_sequence_mismatch(self):
        with tempfile.TemporaryDirectory() as temp:
            source = Path(temp) / "source"
            shutil.copytree(ROOT / "results" / "phage_c6_full_evidence", source)
            fasta = source / "proteins.faa"
            text = fasta.read_text()
            fasta.write_text(text.replace("MISQDKFEYEISAMK", "MSSQDKFEYEISAMK", 1))
            with self.assertRaisesRegex(ValueError, "protein sequence/length mismatch"):
                _load_source(source)

    def test_resume_source_loader_rejects_missing_artifact(self):
        with tempfile.TemporaryDirectory() as temp:
            source = Path(temp) / "source"
            shutil.copytree(ROOT / "results" / "phage_c6_full_evidence", source)
            (source / "evidence.json").unlink()
            with self.assertRaisesRegex(ValueError, "missing source artifacts"):
                _load_source(source)
    def test_progress_stage_transitions_and_non_tty_output(self):
        stream = StringIO()
        progress = ProgressReporter(stream=stream)
        progress.start("Pfam")
        progress.finish("2 accepted hits / 1 proteins")
        text = stream.getvalue()
        self.assertIn("RUNNING: Pfam", text)
        self.assertIn("DONE: Pfam", text)
        self.assertIn("2 accepted hits", text)

    def test_progress_quiet_mode(self):
        stream = StringIO()
        progress = ProgressReporter(stream=stream, quiet=True)
        progress.start("Pfam")
        progress.finish("summary")
        self.assertEqual(stream.getvalue(), "")

    def test_progress_does_not_enter_scientific_outputs(self):
        with tempfile.TemporaryDirectory() as temp:
            output = Path(temp) / "progress-output"
            run(ROOT / "examples/demo_phage.fasta", output, predictor=DemoORFPredictor(), progress=ProgressReporter(stream=StringIO()))
            self.assertNotIn("RUNNING", (output / "run_manifest.json").read_text())
            self.assertNotIn("DONE", (output / "annotation.tsv").read_text())
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
        protein.sequence = translate(protein.cds)
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

    def test_pfam_domtblout_parser_preserves_domain_evidence(self):
        _, genome = read_fasta(ROOT / "examples/demo_phage.fasta")
        proteins = predict_orfs("demo", genome)
        proteins[0].protein_id = "PM_000001"
        proteins[0].sequence = "M" * 104
        adapter = PfamHMMAdapter(database_version="Pfam-test-1")
        evidence = adapter.parse_domtblout((ROOT / "tests/fixtures/pfam_domtblout.txt").read_text(), proteins, {"status": "REAL", "fixture": True})
        self.assertEqual(len(evidence), 3)
        record = next(item for item in evidence if item.provenance["protein_id"] == "PM_000001")
        self.assertEqual(record.identifier, "PF01813.21")
        self.assertEqual(record.family_name, "ATP-synt_D")
        self.assertEqual(record.source, "Pfam")
        self.assertEqual(record.status, "REAL")
        self.assertEqual(record.evidence_strength, "WEAK")
        self.assertEqual(record.description, "ATP synthase subunit D")
        self.assertEqual(record.coordinates, {"start": 12, "end": 62})
        self.assertAlmostEqual(record.metrics["independent_domain_e_value"], 6.7e-7)
        self.assertAlmostEqual(record.metrics["conditional_domain_e_value"], 0.00082)
        self.assertEqual(record.metrics["query_protein_id"], "PM_000001")
        self.assertEqual(record.metrics["envelope_coordinates"], {"start": 10, "end": 102})
        self.assertEqual(record.provenance["protein_id"], "PM_000001")
        weak = next(item for item in evidence if item.provenance["protein_id"] == "PM_000002")
        self.assertEqual(weak.evidence_strength, "WEAK")
        self.assertFalse("Putative" in weak.description or "protein" == weak.description)

    def test_pfam_thresholds_retain_rejected_raw_hits(self):
        _, genome = read_fasta(ROOT / "examples/demo_phage.fasta")
        proteins = predict_orfs("demo", genome)
        adapter = PfamHMMAdapter(evalue_threshold=1e-5, coverage_threshold=0.4)
        evidence = adapter.parse_domtblout((ROOT / "tests/fixtures/pfam_domtblout.txt").read_text(), proteins, {"status": "REAL"})
        rejected = next(item for item in evidence if item.provenance["protein_id"] == "PM_000002")
        self.assertEqual(rejected.evidence_strength, "REJECTED")
        self.assertFalse(rejected.supports)
        self.assertIn("raw_hmmscan_row", rejected.metrics)

    def test_pfam_ga_mode_classifies_hits_strong_and_records_mode(self):
        adapter = PfamHMMAdapter(threshold_mode="GA", hmmscan="hmmscan")
        self.assertEqual(adapter.threshold_mode, "GA")
        self.assertTrue(adapter.trusted_cutoff)
        self.assertIn("threshold_mode", adapter.provenance())
        self.assertEqual(adapter.provenance()["threshold_mode"], "GA")
        self.assertTrue(adapter.provenance()["trusted_cutoff"])
        command = adapter._search_command(Path("domtblout"), Path("proteins.faa"))
        self.assertIn("--cut_ga", command)

    def test_pfam_ga_parser_marks_reported_hits_strong_and_mining_can_use_them(self):
        _, genome = read_fasta(ROOT / "examples/demo_phage.fasta")
        proteins = predict_orfs("demo", genome)
        proteins[0].sequence = "M" * 104
        adapter = PfamHMMAdapter(threshold_mode="GA", database_version="Pfam-test-1")
        evidence = adapter.parse_domtblout((ROOT / "tests/fixtures/pfam_domtblout.txt").read_text(), proteins)
        ga = next(item for item in evidence if item.provenance["protein_id"] == "PM_000001")
        self.assertEqual(ga.evidence_strength, "STRONG")
        self.assertTrue(ga.supports)
        proteins[0].evidence.append(ga)
        mine(proteins)
        self.assertGreater(proteins[0].biological_interest, 0)

    def test_weak_only_pfam_does_not_enable_ranking(self):
        _, genome = read_fasta(ROOT / "examples/demo_phage.fasta")
        proteins = predict_orfs("demo", genome)
        adapter = PfamHMMAdapter(threshold_mode="NONE")
        evidence = adapter.parse_domtblout((ROOT / "tests/fixtures/pfam_domtblout.txt").read_text(), proteins)
        proteins[0].evidence.append(next(item for item in evidence if item.provenance["protein_id"] == "PM_000001"))
        mine(proteins)
        self.assertEqual(proteins[0].biological_interest, 0)
        self.assertEqual(proteins[0].evidence_diversity, "Low")

    def test_pipeline_persists_accepted_ga_pfam_evidence(self):
        with tempfile.TemporaryDirectory() as temp:
            output = Path(temp) / "ga-persistence"
            _, genome = read_fasta(ROOT / "examples/demo_phage.fasta")
            proteins = predict_orfs("demo", genome)
            accepted = Evidence(
                "domain", "Pfam domain evidence", EvidenceLevel.COMPUTATIONAL,
                "Pfam", "Pfam-test", status="REAL", supports=True,
                identifier="PF00001", evidence_strength="STRONG",
                family_name="Test_family", description="Test domain",
                metrics={"query_protein_id": proteins[0].protein_id,
                         "independent_domain_e_value": 1e-8},
                provenance={"threshold_mode": "GA", "trusted_cutoff": True},
            )
            result = EvidenceAdapterResult(
                "PfamHMMAdapter", "REAL", evidence=[accepted],
                provenance={"threshold_mode": "GA", "trusted_cutoff": True},
            )
            with patch("phagemine.pipeline.PfamHMMAdapter") as adapter:
                adapter.return_value.analyze.return_value = result
                run(ROOT / "examples/demo_phage.fasta", output,
                    predictor=DemoORFPredictor(), pfam_threshold_mode="GA")
            records = json.loads((output / "evidence.json").read_text())
            persisted = [e for p in records for e in p["evidence"] if e["source"] == "Pfam"]
            self.assertEqual(len(persisted), 1)
            self.assertEqual(persisted[0]["status"], "REAL")
            self.assertEqual(persisted[0]["provenance"]["threshold_mode"], "GA")
            self.assertTrue(persisted[0]["provenance"]["trusted_cutoff"])

    def test_pfam_unavailable_has_no_fabricated_evidence(self):
        result = PfamHMMAdapter("/definitely/missing/Pfam-A.hmm", "/definitely/missing/hmmscan").analyze([])
        self.assertEqual(result.status, "UNAVAILABLE")
        self.assertEqual(result.evidence, [])
        self.assertIn("no domain evidence was fabricated", result.message)

    def test_vog_parser_maps_orthology_and_annotations(self):
        _, genome = read_fasta(ROOT / "examples/demo_phage.fasta")
        proteins = predict_orfs("demo", genome)
        with tempfile.TemporaryDirectory() as temp:
            annotations = Path(temp) / "vog.annotations.tsv"
            annotations.write_text("#GroupName\tProteinCount\tSpeciesCount\tFunctionalCategory\tConsensusFunctionalDescription\nVOG00001\t1\t1\tS\tstructural viral protein\n")
            adapter = VOGHMMAdapter(annotations_path=annotations, evalue_threshold=1e-5, coverage_threshold=0.5)
            evidence = adapter.parse_domtblout((ROOT / "tests/fixtures/vog_domtblout.txt").read_text(), proteins, {"status": "REAL"})
        strong = next(item for item in evidence if item.identifier == "VOG00001")
        self.assertEqual(strong.source, "VOGDB")
        self.assertEqual(strong.evidence_strength, "STRONG")
        self.assertEqual(strong.metrics["query_coordinates"], {"start": 8, "end": 96})
        self.assertEqual(strong.metrics["vog_id"], "VOG00001")
        self.assertEqual(strong.metrics["functional_category"], "S")
        self.assertEqual(strong.description, "structural viral protein")
        self.assertEqual(strong.metrics["consensus_functional_description"], "structural viral protein")
        self.assertEqual(strong.provenance["threshold_mode"], "MANUAL")

    def test_vog_unavailable_does_not_fabricate_evidence(self):
        result = VOGHMMAdapter("/definitely/missing/vog.hmm").analyze([])
        self.assertEqual(result.status, "UNAVAILABLE")
        self.assertEqual(result.evidence, [])
        self.assertIn("no VOG evidence was fabricated", result.message)

    def test_vog_rejected_rows_are_not_supportive(self):
        _, genome = read_fasta(ROOT / "examples/demo_phage.fasta")
        proteins = predict_orfs("demo", genome)
        evidence = VOGHMMAdapter(evalue_threshold=1e-5, coverage_threshold=0.5).parse_domtblout(
            (ROOT / "tests/fixtures/vog_domtblout.txt").read_text(), proteins)
        accepted = next(item for item in evidence if item.identifier == "VOG00001")
        rejected = next(item for item in evidence if item.identifier == "VOG00002")
        self.assertTrue(accepted.supports)
        self.assertEqual(accepted.evidence_strength, "STRONG")
        self.assertFalse(rejected.supports)
        self.assertEqual(rejected.evidence_strength, "REJECTED")

    def test_registered_vog_version_propagates(self):
        with tempfile.TemporaryDirectory() as temp:
            registry = Path(temp) / "resources.json"
            hmm = Path(temp) / "VOGDB.hmm"
            hmm.write_text("synthetic")
            manager = EvidenceResourceManager(registry)
            manager.register("VOGDB-236", ResourceType.VOGDB, hmm, version="236")
            registered = manager.find(ResourceType.VOGDB)
            adapter = VOGHMMAdapter(registered["path"], database_version=registered["version"], hmmscan="/missing/hmmscan")
            self.assertEqual(adapter.provenance()["vogdb_version"], "236")
            _, genome = read_fasta(ROOT / "examples/demo_phage.fasta")
            proteins = predict_orfs("demo", genome)
            row = adapter.parse_domtblout((ROOT / "tests/fixtures/vog_domtblout.txt").read_text(), proteins)[0]
            self.assertEqual(row.source_version, "236")
            self.assertEqual(row.provenance["vogdb_version"], "236")

    def test_vog_is_distinct_evidence_source(self):
        evidence = Evidence("viral_orthology", "VOG evidence", EvidenceLevel.COMPUTATIONAL, "VOGDB", "test", status="REAL", evidence_strength="STRONG", provenance={"protein_id": "PM_000001"})
        self.assertNotEqual(evidence.source, "Pfam")

    def test_phrogs_parser_preserves_phrogs_semantics(self):
        _, genome = read_fasta(ROOT / "examples/demo_phage.fasta")
        proteins = predict_orfs("demo", genome)
        with tempfile.TemporaryDirectory() as temp:
            annotations = Path(temp) / "phrogs.tsv"
            annotations.write_text("phrog\tcolor\tannot\tcategory\n1\t#fff\tunknown function\tunknown\n2\t#000\tDNA-associated protein\tDNA, RNA and nucleotide metabolism\n")
            adapter = PHROGSMMseqsAdapter(annotations_path=annotations, evalue_threshold=1e-5, coverage_threshold=0.5)
            evidence = adapter.parse_tabular((ROOT / "tests/fixtures/phrogs_mmseqs.tsv").read_text(), proteins)
        strong = next(item for item in evidence if item.identifier == "1")
        rejected = next(item for item in evidence if item.identifier == "2")
        self.assertEqual(strong.modality, "phage_orthology")
        self.assertEqual(strong.metrics["functional_category"], "unknown")
        self.assertIsNone(strong.description)
        self.assertTrue(strong.supports)
        self.assertEqual(strong.metrics["mmseqs_score"], 75.0)
        self.assertEqual(strong.metrics["sequence_identity"], 32.5)
        self.assertEqual(strong.metrics["bit_score"], 75.0)
        self.assertEqual(strong.metrics["percent_identity"], 32.5)
        self.assertEqual(strong.metrics["raw_mmseqs_row"].split("\t")[0], "PM_000001")
        self.assertEqual(adapter.provenance()["output_format"], MMSEQS_FORMAT)
        self.assertEqual(strong.provenance["threshold_mode"], "MANUAL")
        self.assertEqual(rejected.evidence_strength, "REJECTED")
        self.assertFalse(rejected.supports)

    def test_phrogs_unavailable_does_not_fabricate_evidence(self):
        result = PHROGSMMseqsAdapter("/missing/phrogs", mmseqs="/missing/mmseqs").analyze([])
        self.assertEqual(result.status, "UNAVAILABLE")
        self.assertEqual(result.evidence, [])

    def test_phrogs_manual_identity_and_alignment_thresholds_reject(self):
        _, genome = read_fasta(ROOT / "examples/demo_phage.fasta")
        proteins = predict_orfs("demo", genome)
        adapter = PHROGSMMseqsAdapter(
            evalue_threshold=None, coverage_threshold=None,
            identity_threshold=40.0, alignment_length_threshold=80)
        evidence = adapter.parse_tabular(
            (ROOT / "tests/fixtures/phrogs_mmseqs.tsv").read_text(), proteins)
        self.assertTrue(evidence)
        self.assertTrue(all(not item.supports for item in evidence))
        self.assertTrue(all(item.evidence_strength == "REJECTED" for item in evidence))

    def test_phrogs_conflicting_strong_annotations_are_preserved(self):
        _, genome = read_fasta(ROOT / "examples/demo_phage.fasta")
        proteins = predict_orfs("demo", genome)
        with tempfile.TemporaryDirectory() as temp:
            annotations = Path(temp) / "phrogs.tsv"
            annotations.write_text(
                "phrog\tcolor\tannot\tcategory\n"
                "1\t#fff\tintegrase\tintegration and excision\n"
                "2\t#000\tportal protein\thead and packaging\n")
            adapter = PHROGSMMseqsAdapter(annotations_path=annotations)
            text = (
                "PM_000001\t1\t50\t90\t1\t90\t104\t1\t90\t100\t1e-20\t100\t0.86\t0.9\n"
                "PM_000001\t2\t45\t85\t2\t86\t104\t3\t87\t120\t1e-15\t90\t0.82\t0.71\n")
            evidence = adapter.parse_tabular(text, proteins)
        adapter._mark_conflicts(evidence)
        self.assertEqual({item.description for item in evidence}, {"integrase", "portal protein"})
        self.assertTrue(all(item.metrics["conflict"] for item in evidence))

    def test_swissprot_parser_metadata_and_strengths(self):
        _, genome = read_fasta(ROOT / "examples/demo_phage.fasta")
        proteins = predict_orfs("demo", genome)
        with tempfile.TemporaryDirectory() as temp:
            dat = Path(temp) / "uniprot.dat"
            dat.write_text("ID   001R_FRG3G              Reviewed;         100 AA.\nAC   Q6GZX4;\nDE   RecName: Full=Curated viral protein;\nOS   Test virus.\nOX   NCBI_TaxID=12345;\nGN   Name=geneA;\nPE   1: Evidence at protein level;\n//\nID   002L_FRG3G              Reviewed;         100 AA.\nAC   Q6GZX3;\nDE   RecName: Full=Second curated protein;\n//\n")
            adapter = SwissProtEvidenceAdapter(metadata_path=dat)
            evidence = adapter.parse_tabular((ROOT / "tests/fixtures/swissprot_diamond.tsv").read_text(), proteins)
        strong = next(item for item in evidence if item.identifier == "Q6GZX4")
        moderate = next(item for item in evidence if item.identifier == "Q6GZX3")
        rejected = next(item for item in evidence if item.provenance["protein_id"] == "PM_000003")
        self.assertEqual(strong.evidence_strength, "STRONG")
        self.assertTrue(strong.supports)
        self.assertEqual(strong.metrics["entry_name"], "001R_FRG3G")
        self.assertEqual(strong.metrics["protein_name"], "Curated viral protein")
        self.assertEqual(strong.metrics["taxonomy_id"], "12345")
        self.assertAlmostEqual(strong.metrics["query_coverage"], 0.8)
        self.assertAlmostEqual(strong.metrics["subject_coverage"], 0.8)
        self.assertAlmostEqual(strong.metrics["evalue"], 1e-30)
        self.assertLessEqual(strong.metrics["query_coverage"], 1.0)
        self.assertLessEqual(strong.metrics["subject_coverage"], 1.0)
        self.assertEqual(moderate.evidence_strength, "MODERATE")
        self.assertTrue(moderate.supports)
        self.assertEqual(rejected.evidence_strength, "REJECTED")
        self.assertFalse(rejected.supports)

    def test_swissprot_full_length_coverage_is_one(self):
        _, genome = read_fasta(ROOT / "examples/demo_phage.fasta")
        proteins = predict_orfs("demo", genome)
        row = "PM_000001\tsp|Q6GZX4|001R_FRG3G\t80.0\t120\t120\t120\t1\t120\t1\t120\t1e-30\t120.0\n"
        record = SwissProtEvidenceAdapter().parse_tabular(row, proteins)[0]
        self.assertEqual(record.metrics["query_coverage"], 1.0)
        self.assertEqual(record.metrics["subject_coverage"], 1.0)

    def test_swissprot_subject_id_and_unavailable(self):
        self.assertEqual(SwissProtEvidenceAdapter.parse_subject_id("sp|Q6GZX4|001R_FRG3G"), ("Q6GZX4", "001R_FRG3G"))
        result = SwissProtEvidenceAdapter("/missing/db.dmnd", diamond="/missing/diamond").analyze([])
        self.assertEqual(result.status, "UNAVAILABLE")
        self.assertEqual(result.evidence, [])

    def test_pipeline_manifest_marks_pfam_unavailable_and_no_mock_by_default(self):
        with tempfile.TemporaryDirectory() as temp:
            output = Path(temp) / "real-no-db"
            run(ROOT / "examples/demo_phage.fasta", output, predictor=DemoORFPredictor(), pfam_path="/definitely/missing/Pfam-A.hmm")
            manifest = json.loads((output / "run_manifest.json").read_text())
            self.assertEqual(manifest["evidence_adapters"][0]["status"], "UNAVAILABLE")
            self.assertNotIn("MOCK", {adapter["status"] for adapter in manifest["evidence_adapters"]})
            evidence = json.loads((output / "evidence.json").read_text())
            self.assertFalse(any(item["source"] == "mock-phage-evidence" for protein in evidence for item in protein["evidence"]))

    def test_no_evidence_does_not_assign_misleading_ordinal_ranks(self):
        with tempfile.TemporaryDirectory() as temp:
            output = Path(temp) / "insufficient"
            run(ROOT / "examples/demo_phage.fasta", output, predictor=DemoORFPredictor())
            manifest = json.loads((output / "run_manifest.json").read_text())
            self.assertEqual(manifest["discovery_ranking"]["status"], "INSUFFICIENT_EVIDENCE")
            ranking = (output / "candidate_ranking.tsv").read_text().splitlines()
            self.assertTrue(all(line.startswith("NA\t") for line in ranking[1:]))
            self.assertIn("INSUFFICIENT_EVIDENCE", (output / "report.md").read_text())

    def test_resource_manager_registers_pfam_and_persists_metadata(self):
        with tempfile.TemporaryDirectory() as temp:
            registry = Path(temp) / "resources.json"
            database = Path(temp) / "Pfam-A.hmm"
            database.write_text("synthetic HMM")
            manager = EvidenceResourceManager(registry)
            manager.register("PFAM", ResourceType.PFAM, database, version="Pfam-test-1", required_tools=["hmmscan"])
            loaded = EvidenceResourceManager(registry).get("PFAM")
            self.assertEqual(loaded["resource_type"], "PFAM")
            self.assertEqual(loaded["status"], ResourceStatus.READY.value)
            self.assertEqual(loaded["version"], "Pfam-test-1")
            self.assertIsNone(loaded["checksum"])

    def test_resource_manager_supports_types_missing_paths_and_removal(self):
        with tempfile.TemporaryDirectory() as temp:
            registry = Path(temp) / "resources.json"
            manager = EvidenceResourceManager(registry)
            manager.register("SWISSPROT", ResourceType.SWISSPROT, Path(temp) / "missing.fasta")
            manager.register("REFSEQ", "REFSEQ", Path(temp) / "missing-refseq")
            self.assertEqual({item["resource_type"] for item in manager.list()}, {"SWISSPROT", "REFSEQ"})
            self.assertTrue(all(item["status"] == "UNAVAILABLE" for item in manager.list()))
            self.assertTrue(manager.unregister("SWISSPROT"))
            self.assertFalse(manager.unregister("SWISSPROT"))

    def test_resource_manager_does_not_download_or_mutate_database(self):
        with tempfile.TemporaryDirectory() as temp:
            registry = Path(temp) / "resources.json"
            database = Path(temp) / "dummy.hmm"
            database.write_text("unchanged")
            manager = EvidenceResourceManager(registry)
            manager.register("PFAM", "PFAM", database)
            self.assertEqual(database.read_text(), "unchanged")
            self.assertFalse((Path(temp) / "download").exists())

    def test_pfam_registered_fallback_and_explicit_precedence(self):
        with tempfile.TemporaryDirectory() as temp:
            registry = Path(temp) / "resources.json"
            registered = Path(temp) / "registered.hmm"
            explicit = Path(temp) / "explicit.hmm"
            registered.write_text("registered")
            explicit.write_text("explicit")
            manager = EvidenceResourceManager(registry)
            manager.register("PFAM", "PFAM", registered)
            found = manager.find(ResourceType.PFAM)
            self.assertEqual(found["path"], str(registered))
            # The pipeline's precedence decision is explicit path first; this mirrors the adapter input contract.
            explicit_adapter = PfamHMMAdapter(explicit, "/missing/hmmscan")
            self.assertEqual(explicit_adapter.pfam_path, explicit)

    def test_default_registry_is_user_level_not_repository(self):
        registry = default_registry_path()
        self.assertNotIn(str(ROOT), str(registry))

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
