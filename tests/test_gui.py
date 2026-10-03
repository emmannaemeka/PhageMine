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


def test_launcher_binds_to_localhost(monkeypatch):
    from phagemine.gui import launcher
    import sys
    from types import ModuleType
    cli = ModuleType("streamlit.web.cli")
    commands = []
    cli.main = lambda: commands.append(list(sys.argv)) or 0
    streamlit = ModuleType("streamlit")
    streamlit.__path__ = []
    web = ModuleType("streamlit.web")
    web.__path__ = []
    web.cli = cli
    streamlit.web = web
    monkeypatch.setitem(sys.modules, "streamlit", streamlit)
    monkeypatch.setitem(sys.modules, "streamlit.web", web)
    monkeypatch.setitem(sys.modules, "streamlit.web.cli", cli)
    monkeypatch.setattr(launcher.importlib.util, "find_spec", lambda name: True)
    monkeypatch.setattr(sys, "argv", ["phagemine-gui"])
    assert launcher.main([]) == 0
    assert "--server.address=127.0.0.1" in commands[0]


def test_batch_status_requires_end_and_respects_failures(tmp_path):
    from phagemine.gui.services.status import read_run_status
    path = tmp_path / "batch_manifest.json"
    path.write_text(json.dumps({"ended_at": None, "failures": []}))
    assert read_run_status(tmp_path)["state"] == "RUNNING"
    path.write_text(json.dumps({"ended_at": "2026-10-03", "failures": [{"sample_id": "a"}]}))
    assert read_run_status(tmp_path)["state"] == "FAILED"
    path.write_text("not json")
    status = read_run_status(tmp_path)
    assert status["state"] == "INCOMPLETE"
    assert status["warnings"]


def test_failed_sample_without_error_message_is_failed(tmp_path):
    from phagemine.gui.services.status import read_run_status
    (tmp_path / "sample_status.json").write_text(json.dumps({"status": "FAILED"}))
    assert read_run_status(tmp_path)["state"] == "FAILED"


def test_both_layout_does_not_duplicate_annotations(tmp_path):
    from phagemine.gui.services.results import annotation_records, discovery_tables
    from phagemine.gui.services.status import read_run_status
    annotation = tmp_path / "annotation"
    discovery = tmp_path / "discovery"
    annotation.mkdir(); discovery.mkdir()
    (annotation / "batch_manifest.json").write_text(json.dumps({"ended_at": "2026-10-03"}))
    for sample in (annotation / "a", annotation / "b", discovery / "a", discovery / "b"):
        sample.mkdir()
        (sample / "annotation.tsv").write_text("protein_id\tproduct\nP1\thypothetical protein\n")
    assert read_run_status(tmp_path)["state"] == "INCOMPLETE"
    (discovery / "pmf_families.tsv").write_text("pmf_id\nPMF1\n")
    (discovery / "pooled_manifest.json").write_text(json.dumps({"pmf_count": 1}))
    (discovery / "discovery_report.html").write_text("<html></html>")
    assert read_run_status(tmp_path)["state"] == "COMPLETE"
    records = annotation_records(tmp_path)
    assert len(records) == 2
    assert {row["record_id"] for row in records} == {"a / P1", "b / P1"}
    assert discovery_tables(tmp_path)["families"][0]["pmf_id"] == "PMF1"


def test_exports_include_submission_fasta_and_exclude_symlinks(tmp_path):
    from phagemine.gui.services.results import downloadable_files, export_files
    root = tmp_path / "run"; root.mkdir()
    (root / "genome.fsa").write_text(">g\nATG\n")
    private = tmp_path / "outside.txt"; private.write_text("outside")
    (root / "linked.txt").symlink_to(private)
    assert [p.name for p in export_files(root)] == ["genome.fsa"]
    assert [p.name for p in downloadable_files(root)] == ["genome.fsa"]


def test_upload_collision_does_not_overwrite_input(tmp_path):
    from phagemine.gui.services.execution import save_uploads
    from types import SimpleNamespace
    uploads = [SimpleNamespace(name="genome.fasta", getvalue=lambda: b">g\nATG\n") for _ in range(2)]
    with pytest.raises(ValueError, match="unique"):
        save_uploads(uploads, tmp_path)
    assert not list(tmp_path.iterdir())


def test_input_validation_rejects_wrong_type_and_empty_directory(tmp_path):
    from phagemine.gui.services.execution import validate_input, build_cli_args
    with pytest.raises(ValueError, match="FASTA file"):
        validate_input(tmp_path, cohort=False)
    with pytest.raises(ValueError, match="no FASTA"):
        validate_input(tmp_path, cohort=True)
    with pytest.raises(ValueError, match="empty"):
        build_cli_args("", tmp_path)


def test_malformed_record_schema_is_actionable(tmp_path):
    from phagemine.gui.services.results import annotation_records
    (tmp_path / "annotation.tsv").write_text("protein_id\nP1\n")
    (tmp_path / "evidence.json").write_text('{"unexpected": true}')
    with pytest.raises(ValueError, match="Expected an array"):
        annotation_records(tmp_path)
