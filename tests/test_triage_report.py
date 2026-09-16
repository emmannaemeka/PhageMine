import csv
from phagemine.triage import build_triage, write_triage_report

def test_empty_triage_report_is_readable(tmp_path):
    path = write_triage_report(tmp_path)
    text = path.read_text()
    assert "Final genes: **0**" in text
    assert "No unresolved gene-model conflicts" in text

def test_triage_surfaces_conflict_rescue_and_boundary_change(tmp_path):
    calls = tmp_path / "gene_calls"; calls.mkdir()
    (calls / "model_decisions.tsv").write_text("locus_id\tdecision_class\tdecision_reason\nL1\tSTRAND_CONFLICT_UNRESOLVED\tOpposite strands\n")
    (calls / "rescue_candidates.tsv").write_text("candidate_id\tprovider\tstart\tend\tevidence_status\nC1\tpyrodigal\t20\t90\tWEAK\n")
    (calls / "final_gene_models.tsv").write_text("locus_id\tstart\tend\tchanged_from_legacy\tselection_reason\nL1\t10\t100\ttrue\tevidence\n")
    write_triage_report(tmp_path)
    data = build_triage(tmp_path)
    assert data["counts"] == {"final_genes": 1, "unresolved_conflicts": 1, "rescue_candidates": 1, "boundary_changes": 1, "unresolved_function_proteins": 0}
