import csv
import json

import pytest

from phagemine.deeporigins import main, read_evidence, summarize

HEADER = ["query", "target", "source_group", "sequence_identity",
          "alignment_coverage", "structural_similarity", "method", "reference"]


def test_cross_group_low_identity_requires_review_not_ancestry():
    rows = [
        dict(zip(HEADER, ["p1", "b1", "bacteria", 0.12, 0.85, 0.8, "example", "ref1"])),
        dict(zip(HEADER, ["p1", "v1", "phage", 0.18, 0.9, 0.9, "example", "ref2"])),
    ]
    report = summarize(rows)[0]
    assert report["review_priority"] == "manual_review"
    assert report["interpretation"] == "relationship_unresolved"


def test_single_group_is_background():
    rows = [dict(zip(HEADER, ["p1", "v1", "phage", 0.1, 0.9, 0.9, "example", "ref1"]))]
    assert summarize(rows)[0]["review_priority"] == "background"


def test_invalid_fraction_rejected(tmp_path):
    file = tmp_path / "bad.tsv"
    file.write_text("\t".join(HEADER) + "\n" +
                    "\t".join(["p", "t", "phage", "110", "0.9", "0.8", "test", "ref"]) + "\n")
    with pytest.raises(ValueError, match="sequence_identity"):
        read_evidence(file)


def test_cli_output(tmp_path):
    file = tmp_path / "evidence.tsv"
    output = tmp_path / "report.json"
    file.write_text("\t".join(HEADER) + "\n" +
                    "\t".join(["p", "t", "bacteria", "0.1", "0.9", "0.8", "test", "ref"]) + "\n")
    assert main([str(file), "--output", str(output)]) == 0
    assert json.loads(output.read_text())["status"] == "exploratory_not_validated"
