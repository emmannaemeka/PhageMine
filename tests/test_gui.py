from __future__ import annotations

import json
import inspect
import zipfile
from io import BytesIO
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


@pytest.mark.parametrize("suffix", [".fa", ".fasta", ".fna", ".fas", ".FASTA"])
def test_supported_upload_extensions_and_filename_sanitizing(suffix):
    from phagemine.gui.services.inputs import sanitize_filename
    assert sanitize_filename(f"../../unsafe genome;name{suffix}") == "unsafe_genome_name" + suffix.lower()


def test_one_uploaded_fasta_is_single_genome(tmp_path):
    from phagemine.gui.services.inputs import UploadedGenome, stage_uploads
    staged = stage_uploads([UploadedGenome("one genome.fasta", b">one\nATGAAATAA\n")])
    try:
        assert staged.cohort is False
        assert staged.path.is_file()
        assert staged.path.name == "one_genome.fasta"
    finally:
        import shutil
        shutil.rmtree(staged.root)


def test_multiple_uploaded_fastas_are_cohort():
    from phagemine.gui.services.inputs import UploadedGenome, stage_uploads
    staged = stage_uploads([UploadedGenome("one.fa", b">one\nATG\n"), UploadedGenome("two.fna", b">two\nATG\n")])
    try:
        assert staged.cohort is True
        assert staged.path == staged.root
        assert len(staged.files) == 2
    finally:
        import shutil
        shutil.rmtree(staged.root)


def test_missing_or_malformed_upload_is_rejected():
    from phagemine.gui.services.inputs import UploadedGenome, stage_uploads, validate_fasta_bytes
    with pytest.raises(ValueError, match="at least one"):
        stage_uploads([])
    with pytest.raises(ValueError, match="exactly one FASTA record"):
        validate_fasta_bytes("broken.fa", b"not fasta\n")
    with pytest.raises(ValueError, match="standard IUPAC"):
        validate_fasta_bytes("broken.fasta", b">g\nATGX\n")


def test_duplicate_uploaded_names_remain_distinct():
    from phagemine.gui.services.inputs import UploadedGenome, stage_uploads
    staged = stage_uploads([UploadedGenome("same.fa", b">a\nATG\n"), UploadedGenome("same.fa", b">b\nATG\n")])
    try:
        assert [path.name for path in staged.files] == ["same.fa", "same_2.fa"]
    finally:
        import shutil
        shutil.rmtree(staged.root)


def test_launcher_is_loopback_only_and_browser_behavior_is_testable(tmp_path):
    from phagemine.gui.launcher import LOOPBACK, streamlit_arguments
    args = streamlit_arguments(tmp_path / "app.py", [], 8765)
    assert LOOPBACK == "127.0.0.1"
    assert "--server.address=127.0.0.1" in args
    assert "--browser.serverAddress=127.0.0.1" in args
    assert "--server.port=8765" in args
    assert "0.0.0.0" not in " ".join(args)
    assert "::" not in " ".join(args)


def test_launcher_opens_browser_only_after_local_health_check(monkeypatch):
    from phagemine.gui import launcher

    class HealthyResponse:
        status = 200
        def __enter__(self): return self
        def __exit__(self, *args): return False

    opened = []
    monkeypatch.setattr(launcher.urllib.request, "urlopen", lambda url, timeout: HealthyResponse())
    monkeypatch.setattr(launcher.webbrowser, "open", lambda url, new: opened.append((url, new)))
    launcher._open_when_ready("http://127.0.0.1:8765", timeout=0.1)
    assert opened == [("http://127.0.0.1:8765", 1)]


def test_workflow_readiness_uses_doctor_capabilities():
    from phagemine.gui.services.status import workflow_readiness
    missing = {"capabilities": {"CORE_ANALYSIS": "UNAVAILABLE"}, "executables": []}
    assert workflow_readiness(missing, mode="annotate", evidence="core", cohort=False)[0] is False
    ready = {"capabilities": {"CORE_ANALYSIS": "READY", "STANDARD_EVIDENCE": "READY", "FULL_EVIDENCE": "READY"},
             "executables": [{"name": "mmseqs", "status": "READY"}]}
    assert workflow_readiness(ready, mode="annotate", evidence="core", cohort=False)[0] is True
    assert workflow_readiness(ready, mode="discover", evidence="full", cohort=True)[0] is True
    ready["executables"][0]["status"] = "MISSING"
    assert workflow_readiness(ready, mode="discover", evidence="core", cohort=True)[0] is False


def test_result_export_zip_preserves_native_files(tmp_path):
    from phagemine.gui.services.results import result_zip
    (tmp_path / "annotation.tsv").write_text("protein_id\nP1\n")
    (tmp_path / "figures").mkdir(); (tmp_path / "figures" / "map.svg").write_text("<svg/>")
    with zipfile.ZipFile(BytesIO(result_zip(tmp_path))) as archive:
        assert sorted(archive.namelist()) == ["annotation.tsv", "figures/map.svg"]
        assert archive.read("annotation.tsv") == b"protein_id\nP1\n"


def test_gui_subprocesses_explicitly_disable_shell_execution():
    from phagemine.gui.services import execution, runs
    assert "shell=False" in inspect.getsource(execution.execute)
    assert "shell=False" in inspect.getsource(runs.start_run)
    assert "shell=True" not in inspect.getsource(execution)
    assert "shell=True" not in inspect.getsource(runs)


def test_frozen_execution_reenters_same_phagemine_cli(monkeypatch):
    from phagemine.gui.services import execution
    monkeypatch.setattr(execution.sys, "frozen", True, raising=False)
    monkeypatch.setattr(execution.sys, "executable", r"C:\Program Files\PhageMine\PhageMine.exe")
    assert execution.execution_argv(["phagemine", "run", "input with spaces.fa"]) == [
        r"C:\Program Files\PhageMine\PhageMine.exe", "--phagemine-cli", "run", "input with spaces.fa"]


def test_desktop_entry_uses_package_launcher(monkeypatch):
    from phagemine.gui import desktop, launcher
    received = []
    monkeypatch.setattr(desktop.sys, "argv", ["PhageMine", "--no-browser"])
    monkeypatch.setattr(launcher, "main", lambda argv: received.append(argv) or 0)
    assert desktop.main() == 0
    assert received == [["--no-browser"]]
    assert desktop.os.environ["PHAGEMINE_DESKTOP"] == "1"
