from __future__ import annotations

import csv
from pathlib import Path
from types import SimpleNamespace

from phagemine.batch import _comparative_resources
from phagemine.family_external import validate_external
from phagemine.inphared import compare_genomes


def test_compare_genomes_writes_ranked_conservative_results(tmp_path, monkeypatch):
    mash = tmp_path / "mash"
    mash.write_text("executable fixture")
    index = tmp_path / "inphared.msh"
    index.write_bytes(b"index")
    metadata = tmp_path / "genome_metadata.tsv"
    metadata.write_text(
        "accession\tdescription\tgenome_length_kb\tgc_percent\tphage_genus\tphage_subfamily\tphage_family\thost_genus\tsource_database\tsource_release\n"
        "REF1\tReference one\t40\t50\tTestvirus\tTestvirinae\tTestviridae\tPseudomonas\tINPHARED\t2026-04-07\n"
    )
    query = tmp_path / "query.fna"
    query.write_text(">query\nACGT\n")

    def fake_run(command, **_kwargs):
        return SimpleNamespace(returncode=0, stdout="REF1\tquery\t0.05\t1e-20\t900/1000\n", stderr="")

    monkeypatch.setattr("phagemine.inphared.subprocess.run", fake_run)
    payload = compare_genomes(
        {"sample": query}, mash_index=index, metadata=metadata,
        output=tmp_path / "output", mash=str(mash),
    )
    assert payload["status"] == "COMPLETE"
    assert payload["matches"][0]["reference_accession"] == "REF1"
    assert payload["matches"][0]["interpretation"] == "SCREENING_ONLY_REQUIRES_CONFIRMATORY_ALIGNMENT_OR_ANI"
    rows = list(csv.DictReader((tmp_path / "output" / "inphared_nearest_phages.tsv").open(), delimiter="\t"))
    assert rows[0]["host_genus"] == "Pseudomonas"


def test_batch_resolves_only_validated_comparative_resources(monkeypatch):
    resources = [
        {"resource_type": "PMFDB", "status": "READY", "path": "/db/pmfdb"},
        {"resource_type": "INPHARED_GENOMES", "status": "READY", "path": "/db/inphared", "provenance": {}},
        {"resource_type": "PMFDB", "status": "INVALID", "path": "/db/broken"},
    ]
    monkeypatch.setattr("phagemine.batch.EvidenceResourceManager.validate_all", lambda self, check_checksum=True: resources)
    pmfdb, inphared = _comparative_resources()
    assert pmfdb == "/db/pmfdb"
    assert inphared["path"] == "/db/inphared"


def test_inphared_product_label_remains_predicted_evidence(tmp_path, monkeypatch):
    executable = tmp_path / "mmseqs"
    executable.write_text("executable fixture")
    family_db = tmp_path / "families"
    family_db.mkdir()
    (family_db / "families.json").write_text(
        '[{"family_id":"PMF1","representative_protein_id":"P1",'
        '"representative_sequence":"MPEPTIDE"}]'
    )
    reference = tmp_path / "reference.faa"
    reference.write_text(">REF1 terminase large subunit\nMPEPTIDE\n")
    metadata = tmp_path / "metadata.tsv"
    metadata.write_text(
        "external_protein_id\tannotation\tannotation_status\tcharacterized\n"
        "REF1\tterminase large subunit\tPREDICTED_FUNCTION\tfalse\n"
    )

    def fake_run(command, **_kwargs):
        Path(command[4]).write_text("PMF1\tREF1\t100\t8\t8\t8\t1e-20\t80\t1\t1\n")
        return SimpleNamespace(returncode=0, stdout="", stderr="")

    monkeypatch.setattr("phagemine.family_external.subprocess.run", fake_run)
    rows = validate_external(
        family_db, reference=reference, metadata=metadata,
        output=tmp_path / "results", mmseqs=str(executable),
    )
    assert rows[0]["characterized_homolog_count"] == 0
    assert rows[0]["predicted_function_homolog_count"] == 1
    assert rows[0]["external_validation_status"] == "EXTERNAL_PREDICTED_FUNCTION_HOMOLOG_FOUND"
