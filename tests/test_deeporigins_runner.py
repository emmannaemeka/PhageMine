import json
from pathlib import Path
from unittest.mock import patch
import pytest
from phagemine.deeporigins_runner import run_foldseek


def test_missing_executable():
    with patch("phagemine.deeporigins_runner.shutil.which", return_value=None):
        with pytest.raises(FileNotFoundError, match="Foldseek not found"):
            run_foldseek("query", "db", "out")


def test_requires_structure_input(tmp_path):
    with patch("phagemine.deeporigins_runner.shutil.which", return_value="/bin/foldseek"):
        with pytest.raises(FileNotFoundError, match="Query structures"):
            run_foldseek(str(tmp_path / "missing"), str(tmp_path / "db"), str(tmp_path / "out"))


def test_runner_builds_command_and_manifest(tmp_path):
    query = tmp_path / "query.pdb"
    db = tmp_path / "referenceDB"
    output = tmp_path / "hits.tsv"
    query.write_text("MODEL\nENDMDL\n")
    db.write_text("placeholder")
    class Completed:
        stdout = "Foldseek Version 1\n"
    def fake_run(cmd, **kwargs):
        if "easy-search" in cmd:
            Path(cmd[4]).write_text("q1\tt1\t0.12\t80\t100\t110\t0.8\t0.73\t1e-8\t100\t0.65\t0.61\n")
        return Completed()
    with patch("phagemine.deeporigins_runner.shutil.which", return_value="/bin/foldseek"), patch(
        "phagemine.deeporigins_runner.subprocess.run", side_effect=fake_run
    ) as mocked:
        report = run_foldseek(str(query), str(db), str(output), threads=2)
    assert report["hit_count"] == 1
    assert report["foldseek_version"] == "Foldseek Version 1"
    assert json.loads((tmp_path / "hits.tsv.manifest.json").read_text())["results_sha256"]
    assert mocked.call_count == 2


def test_no_overwrite(tmp_path):
    query = tmp_path / "q.pdb"
    db = tmp_path / "db"
    out = tmp_path / "out.tsv"
    for file in (query, db, out):
        file.write_text("test")
    with patch("phagemine.deeporigins_runner.shutil.which", return_value="/bin/foldseek"):
        with pytest.raises(FileExistsError, match="Output already exists"):
            run_foldseek(str(query), str(db), str(out))
