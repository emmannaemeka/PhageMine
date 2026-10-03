from pathlib import Path
from unittest.mock import patch

import pytest

from phagemine.batch import batch
from phagemine.cli import main


def test_batch_cli_returns_failure_when_any_sample_fails(capsys):
    with patch("phagemine.cli.batch", return_value=[{"status": "SUCCESS"}, {"status": "FAILED"}]):
        assert main(["batch", "input", "--output", "output"]) == 1
    assert "1 failed sample" in capsys.readouterr().err


def test_batch_cli_success_remains_zero():
    with patch("phagemine.cli.batch", return_value=[{"status": "SUCCESS"}]):
        assert main(["batch", "input", "--output", "output"]) == 0


def test_empty_batch_is_an_error(tmp_path):
    source = tmp_path / "input"; source.mkdir()
    with pytest.raises(ValueError, match="No FASTA"):
        batch(source, tmp_path / "output", gene_predictor="demo")


def test_both_does_not_build_discovery_from_failed_annotations(tmp_path):
    source = tmp_path / "input"; source.mkdir()
    (source / "a.fasta").write_text(">a\nATGAAATAG\n")
    def failed_run(*args, **kwargs):
        raise RuntimeError("caller failed")
    with patch("phagemine.batch.run", side_effect=failed_run), patch("phagemine.batch.discovery_from_annotation") as discover:
        rows = batch(source, tmp_path / "output", mode="both", gene_predictor="demo")
    assert rows[0]["status"] == "FAILED"
    assert not discover.called
    assert (tmp_path / "output/annotation/a/sample_status.json").is_file()
