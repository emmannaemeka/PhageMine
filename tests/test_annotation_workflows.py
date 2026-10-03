from __future__ import annotations

import csv
import json
from pathlib import Path
import shutil
import sys
from unittest.mock import patch
import zipfile

import pytest


@pytest.fixture(scope="module")
def source_run(tmp_path_factory):
    from phagemine.pipeline import run
    from phagemine.gene_prediction import DemoORFPredictor
    root = tmp_path_factory.mktemp("annotation-source")
    source = Path(__file__).resolve().parents[1] / "examples/demo_phage.fasta"
    with patch.dict("os.environ", {"PHAGEMINE_REGISTRY_PATH": str(root / "registry.json")}):
        run(source, root / "result", predictor=DemoORFPredictor(), rna_annotation="none")
    return root / "result"


def test_gff_exports_products_and_does_not_compare_different_genomes(source_run):
    from phagemine.benchmark import import_gff, metrics
    records = import_gff(source_run / "genes.gff3", "PHAGEMINE")
    assert all(record["product"] == "hypothetical protein" for record in records)
    other = [{**record, "genome_id": "different_genome"} for record in records]
    assert metrics(records, other)["exact_matches"] == 0


def test_negative_strand_start_and_stop_labels():
    from phagemine.benchmark import compare_models
    left = [{"start": 10, "end": 100, "strand": "-"}]
    right = [{"start": 10, "end": 110, "strand": "-"}]
    assert compare_models(left, right)[0]["relationship"] == "START_DIFFERENCE"


def test_curation_updates_exports_preserves_evidence_and_chains_audit(source_run, tmp_path):
    from phagemine.curation import apply, template, verify_audit
    from phagemine.operational_validation import validate_run
    original = (source_run / "evidence.json").read_bytes()
    changes = template(source_run)
    changes["reviewer"] = "Test reviewer"
    changes["annotations"] = {"PM_000001": {"product": "reviewed synthetic product; alpha=beta", "reason": "Fixture review", "evidence_reference": "synthetic:test"}}
    out = tmp_path / "reviewed"
    apply(source_run, changes, out)
    assert (out / "evidence.json").read_bytes() == original
    assert "product=reviewed%20synthetic%20product%3B%20alpha%3Dbeta" in (out / "genes.gff3").read_text()
    assert "reviewed synthetic product; alpha=beta" in (out / "genbank_submission/features.tbl").read_text()
    assert verify_audit(out / "curation_audit.jsonl")
    assert validate_run(out)["status"] == "ARTIFACTS_CONSISTENT"
    changes = template(out); changes["reviewer"] = "Second reviewer"
    changes["annotations"] = {"PM_000001": {"product": "hypothetical protein", "reason": "Withdraw assertion"}}
    second = tmp_path / "second"; apply(out, changes, second)
    assert len((second / "curation_audit.jsonl").read_text().splitlines()) == 2
    event = json.loads((second / "curation_audit.jsonl").read_text().splitlines()[-1])
    assert event["before"] == "reviewed synthetic product; alpha=beta"


def test_curation_rejects_stale_templates_coordinates_and_missing_evidence(source_run, tmp_path):
    from phagemine.curation import apply, template
    changes = template(source_run); changes["reviewer"] = "R"
    changes["annotations"] = {"PM_000001": {"start": 2, "reason": "test"}}
    with pytest.raises(ValueError, match="Unsupported edits"): apply(source_run, changes, tmp_path / "bad")
    changes["annotations"] = {"PM_000001": {"product": "specific product", "reason": "test"}}
    with pytest.raises(ValueError, match="evidence reference"): apply(source_run, changes, tmp_path / "bad")
    changes["source_annotation_sha256"] = "wrong"
    with pytest.raises(ValueError, match="changed since review"): apply(source_run, changes, tmp_path / "bad")
    assert not (tmp_path / "bad").exists()


def test_rna_import_exports_features_and_curation_preserves_them(source_run, tmp_path):
    from phagemine.rna_features import import_run
    from phagemine.curation import apply, template
    gff = tmp_path / "rna.gff"
    gff.write_text("##gff-version 3\ndemo_phage\ttest\ttRNA\t10\t60\t42\t-\t.\tID=synthetic;product=tRNA-Ala\n")
    root = tmp_path / "rna"
    result = import_run(source_run, gff, root, "synthetic-provider", "test-1")
    assert result["feature_count"] == 1
    assert "60\t10\ttRNA" in (root / "genbank_submission/features.tbl").read_text()
    scope = json.loads((root / "scientific_validation_status.json").read_text())["feature_scope"]
    assert scope["tRNA"] == "COMPUTATIONALLY_ASSESSED_REQUIRES_REVIEW"
    changes = template(root); changes["reviewer"] = "R"
    changes["annotations"] = {"PM_000001": {"note": "reviewed", "reason": "test"}}
    apply(root, changes, tmp_path / "curated-rna")
    assert "60\t10\ttRNA" in (tmp_path / "curated-rna/genbank_submission/features.tbl").read_text()


def test_rna_rejects_wrong_sequence_and_coordinates(source_run, tmp_path):
    from phagemine.rna_features import import_run
    gff = tmp_path / "rna.gff"
    gff.write_text("wrong\ttest\ttRNA\t10\t60\t.\t+\t.\tID=x\n")
    with pytest.raises(ValueError, match="identifier"): import_run(source_run, gff, tmp_path / "bad", "test", "1")
    assert not (tmp_path / "bad").exists()


def test_trnascan_parser_preserves_zero_hits_and_flags_unsupported_features(tmp_path):
    from phagemine.rna_features import parse_trnascan
    table = tmp_path / "trna.tsv"
    table.write_text("Sequence tRNA Bounds Type Codon Intron Score\nseq 1 70 10 Ala TGC 0 0 61.5\nseq 2 90 110 Undet NNN 0 0 10\n")
    features, rejected = parse_trnascan(table, "test")
    assert features[0]["strand"] == "-"
    assert len(rejected) == 1
    table.write_text("Sequence tRNA Bounds Type Codon Intron Score\n")
    assert parse_trnascan(table, "test") == ([], [])
    table.write_text("not a table")
    with pytest.raises(ValueError, match="Unrecognized"): parse_trnascan(table, "test")


def test_external_genbank_import_requires_exact_sequences(source_run, tmp_path):
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio.SeqFeature import SeqFeature, FeatureLocation
    from phagemine.external_annotation import import_genbank
    from phagemine.io import read_fasta
    seqid, genome = read_fasta(source_run / "analysis_genome.fasta")
    protein = json.loads((source_run / "evidence.json").read_text())[0]
    record = SeqRecord(Seq(genome), id=seqid, annotations={"molecule_type": "DNA"})
    record.features = [SeqFeature(FeatureLocation(protein["start"] - 1, protein["end"], strand=1 if protein["strand"] == "+" else -1), type="CDS", qualifiers={"translation": [protein["sequence"]], "product": ["synthetic external proposal"], "locus_tag": ["external_1"]})]
    path = tmp_path / "external.gbk"; SeqIO.write(record, path, "genbank")
    result = import_genbank(source_run, path, tmp_path / "external", version="test", database_version="synthetic")
    assert result["proposal_count"] == 1
    assert (tmp_path / "external/annotation.tsv").read_bytes() == (source_run / "annotation.tsv").read_bytes()
    record.seq = Seq("A" + genome[1:]) if genome[0] != "A" else Seq("C" + genome[1:])
    SeqIO.write(record, path, "genbank")
    with pytest.raises(ValueError, match="differs"): import_genbank(source_run, path, tmp_path / "bad", version="test", database_version="synthetic")


def test_bundle_round_trip_and_tamper_detection(source_run, tmp_path):
    from phagemine.reproducibility import bundle, verify
    archive = tmp_path / "run.zip"
    assert bundle(source_run, archive)["verification"]["status"] == "VERIFIED"
    assert verify(archive)["scientific_accuracy_validated"] is False
    broken = tmp_path / "broken.zip"
    with zipfile.ZipFile(archive) as source, zipfile.ZipFile(broken, "w") as destination:
        for name in source.namelist():
            destination.writestr(name, b"modified" if name == "annotation.tsv" else source.read(name))
    with pytest.raises(ValueError, match="Checksum mismatch"): verify(broken)


def test_artifact_validator_rejects_corrupt_translation(source_run, tmp_path):
    from phagemine.operational_validation import validate_run
    copy = tmp_path / "copy"; shutil.copytree(source_run, copy)
    (copy / "proteins.faa").write_text(">PM_000001\nINVALID\n")
    with pytest.raises(ValueError, match="Sequences do not agree"): validate_run(copy)


def test_profile_records_failed_command(tmp_path):
    from phagemine.operational_validation import profile
    result = profile([sys.executable, "-c", "raise SystemExit(7)"], tmp_path / "metrics.json")
    assert result["exit_code"] == 7
    assert result["status"] == "FAILED"


def test_benchmark_does_not_double_count_genes_or_capture_rna_qualifiers(tmp_path):
    from phagemine.benchmark import import_gff, import_genbank
    gff = tmp_path / 'genes.gff3'
    gff.write_text('g\ttest\tgene\t1\t9\t.\t+\t.\tID=g1\ng\ttest\tCDS\t1\t9\t.\t+\t0\tID=c1;Parent=g1\n')
    assert len(import_gff(gff, 'test')) == 1
    gb = tmp_path / 'annotations.gbk'
    gb.write_text('LOCUS       genome_one\n     CDS             1..9\n                     /product="CDS product"\n     tRNA            20..30\n                     /product="tRNA-Ala"\nLOCUS       genome_two\n     CDS             1..9\n                     /product="second product"\n//\n')
    rows = import_genbank(gb, 'test')
    assert [(r['genome_id'], r['product']) for r in rows] == [('genome_one', 'CDS product'), ('genome_two', 'second product')]


def test_rna_executable_zero_hit_workflow_and_empty_version(source_run, tmp_path):
    from phagemine.rna_features import scan_run
    executable = tmp_path / 'trnascan'
    executable.write_text('#!/usr/bin/env python3\nimport sys\nfrom pathlib import Path\nif "--version" in sys.argv: print("tRNAscan-SE test"); sys.exit(0)\nPath(sys.argv[sys.argv.index("-o")+1]).write_text("Sequence tRNA Begin End Type Codon Intron Begin End Score\\n")\n')
    executable.chmod(0o755)
    result = scan_run(source_run, tmp_path / 'scanned', str(executable), 1)
    assert result['feature_count'] == 0
    executable.write_text('#!/usr/bin/env python3\n')
    with pytest.raises(RuntimeError, match='no version'): scan_run(source_run, tmp_path / 'failed', str(executable), 1)
    assert not (tmp_path / 'failed').exists()
