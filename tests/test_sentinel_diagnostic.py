import csv
from pathlib import Path

from phagemine.sentinel_diagnostic import diagnose, load_sentinels


def test_sentinel_diagnostic_reports_raw_overlap_and_final_selection(tmp_path):
    raw = tmp_path / "gene_calls" / "raw"; raw.mkdir(parents=True)
    with (raw / "phanotate.tsv").open("w", newline="") as h:
        w = csv.writer(h, delimiter="\t"); w.writerow(["caller", "raw_identifier", "start", "end", "strand"]); w.writerow(["PHANOTATE", "p1", 100, 200, "+"])
    final = tmp_path / "gene_calls" / "final_gene_models.tsv"
    final.write_text("start\tend\tsegment_id\n100\t200\t\n")
    sent = tmp_path / "sentinels.tsv"; sent.write_text("sentinel_id\tstart\tend\nS1\t150\t160\n")
    result = diagnose(tmp_path, load_sentinels(sent))[0]
    assert result["callers"]["phanotate"]["raw_overlap"]
    assert result["stage"] == "FINAL_SELECTED"


def test_sentinel_diagnostic_is_explicit_when_outputs_are_missing(tmp_path):
    sent = tmp_path / "sentinels.tsv"; sent.write_text("sentinel_id\tstart\tend\nS1\t1\t10\n")
    result = diagnose(tmp_path, load_sentinels(sent))[0]
    assert result["stage"] == "STAGE_UNKNOWN_INSUFFICIENT_PROVENANCE"
    assert all(not info["raw_overlap"] for info in result["callers"].values())
