import json
import pytest
from phagemine.deeporigins_foldseek import COLUMNS, main, parse_foldseek, read_groups


def example():
    return "q1\tt1\t0.12\t80\t100\t110\t0.8\t0.73\t1e-8\t100\t0.65\t0.61\n"


def test_import_preserves_distinct_tm_scores(tmp_path):
    source = tmp_path / "foldseek.tsv"
    source.write_text(example())
    hits = parse_foldseek(source, target_groups={"t1": "bacteria"})
    assert hits[0]["qtmscore"] == 0.65
    assert hits[0]["ttmscore"] == 0.61
    assert hits[0]["target_group"] == "bacteria"
    assert hits[0]["interpretation"] == "structural_match_only"


def test_rejects_invalid_score(tmp_path):
    source = tmp_path / "foldseek.tsv"
    source.write_text(example().replace("0.65", "1.5"))
    with pytest.raises(ValueError, match="qtmscore"):
        parse_foldseek(source)


def test_cli_with_metadata(tmp_path):
    source = tmp_path / "foldseek.tsv"
    source.write_text(example())
    mapping = tmp_path / "groups.tsv"
    mapping.write_text("target\tgroup\nt1\tbacteria\n")
    out = tmp_path / "hits.json"
    assert main([str(source), "--output", str(out), "--target-groups", str(mapping),
                 "--database", "example-db-2026", "--foldseek-version", "example"]) == 0
    report = json.loads(out.read_text())
    assert report["hits"][0]["target_group"] == "bacteria"
    assert report["status"] == "unvalidated_structural_matches"


def test_conflicting_group_mapping(tmp_path):
    mapping = tmp_path / "groups.tsv"
    mapping.write_text("target\tgroup\nt1\tbacteria\nt1\tphage\n")
    with pytest.raises(ValueError, match="Conflicting"):
        read_groups(mapping)
