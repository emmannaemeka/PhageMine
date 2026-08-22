from __future__ import annotations

import json
from pathlib import Path

import pytest


def test_gui_services_import_without_streamlit():
    import phagemine.gui
    from phagemine.gui.services import execution, results, status
    assert phagemine.gui.GUI_VERSION == "0.1"
    assert execution.MODES == ("annotate", "discover", "both")


def test_cli_translation_single_and_batch(tmp_path):
    from phagemine.gui.services.execution import build_cli_args, display_command
    fasta = tmp_path / "genome.fasta"; fasta.write_text(">g\nATG\n")
    assert build_cli_args(fasta, tmp_path / "one", mode="both", threads=3) == [
        "phagemine", "run", str(fasta), "--output", str(tmp_path / "one"), "--threads", "3"]
    args = build_cli_args(tmp_path, tmp_path / "batch result", mode="discover", evidence="full", threads=8)
    assert args == ["phagemine", "batch", str(tmp_path), "--output", str(tmp_path / "batch result"),
                    "--mode", "discover", "--evidence", "full", "--threads", "8"]
    assert "'" in display_command(args)
    with pytest.raises(ValueError, match="require cohort"):
        build_cli_args(fasta, tmp_path / "one", evidence="full")


def test_command_metacharacters_remain_literal_arguments(tmp_path):
    from phagemine.gui.services.execution import build_cli_args
    source = tmp_path / "genomes; touch SHOULD_NOT_EXIST"
    output = tmp_path / "results $(touch ALSO_NOT_EXECUTED)"
    args = build_cli_args(source, output, cohort=True)
    assert args[2] == str(source)
    assert args[4] == str(output)
    assert "touch" not in args


def test_read_annotation_outputs(tmp_path):
    from phagemine.gui.services.results import annotation_records
    (tmp_path / "annotation.tsv").write_text("protein_id\tstart\tend\tstrand\tannotation\nP1\t1\t30\t+\thypothetical protein\n")
    (tmp_path / "functional_classification.json").write_text(json.dumps([
        {"protein_id": "P1", "functional_state": "UNRESOLVED", "proposed_function": "Function unresolved", "confidence": "NONE"}]))
    (tmp_path / "evidence.json").write_text(json.dumps([
        {"protein_id": "P1", "sequence": "MKK", "cds": "ATGAAAAAA", "evidence": []}]))
    rows = annotation_records(tmp_path)
    assert rows[0]["protein_id"] == "P1"
    assert rows[0]["functional_state"] == "UNRESOLVED"
    assert rows[0]["sequence"] == "MKK"


def test_read_discovery_outputs(tmp_path):
    from phagemine.gui.services.results import discovery_tables
    (tmp_path / "pmf_families.tsv").write_text("pmf_id\tmember_count\tgenome_count\nPMF1\t3\t2\n")
    (tmp_path / "discovery_ranking.tsv").write_text("rank\tpmf_id\tdiscovery_priority\n1\tPMF1\tPRIORITIZED_RECURRENT_UNKNOWN\n")
    tables = discovery_tables(tmp_path)
    assert tables["families"][0]["pmf_id"] == "PMF1"
    assert tables["ranking"][0]["discovery_priority"] == "PRIORITIZED_RECURRENT_UNKNOWN"
    assert tables["validation"] == []


def test_doctor_status_is_conservative():
    from phagemine.gui.services.status import doctor_rows
    payload = {"executables": [
        {"name": "phanotate.py", "status": "READY", "path": "/bin/phanotate"},
        {"name": "table2asn", "status": "MISSING", "path": None}],
        "resources": [{"resource_type": "PFAM", "status": "INVALID", "validation_errors": ["indexes incomplete"]}]}
    rows = {row["component"]: row for row in doctor_rows(payload)}
    assert rows["PHANOTATE"]["state"] == "READY"
    assert rows["table2asn"]["state"] == "OPTIONAL"
    assert rows["Pfam"]["state"] == "ERROR"
    assert rows["PHROGs"]["state"] == "NOT CONFIGURED"


def test_missing_and_invalid_outputs(tmp_path):
    from phagemine.gui.services.results import annotation_records, discovery_tables
    from phagemine.gui.services.status import read_run_status
    assert read_run_status(tmp_path / "missing")["state"] == "NOT FOUND"
    with pytest.raises(ValueError, match="does not exist"):
        discovery_tables(tmp_path / "missing")
    (tmp_path / "annotation.tsv").write_text("protein_id\nP1\n")
    (tmp_path / "functional_classification.json").write_text("not json")
    with pytest.raises(ValueError, match="Invalid PhageMine output"):
        annotation_records(tmp_path)


def test_optional_dependency_message(monkeypatch, capsys):
    from phagemine.gui import launcher
    monkeypatch.setattr(launcher.importlib.util, "find_spec", lambda name: None)
    assert launcher.main() == 2
    assert 'pip install "phagemine[gui]"' in capsys.readouterr().err
