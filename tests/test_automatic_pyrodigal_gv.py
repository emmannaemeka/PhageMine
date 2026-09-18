import json
from pathlib import Path

from phagemine.gene_callers import GenePredictionResult
from phagemine.gene_models import GeneModel
from phagemine.reconciliation_engine import run_automatic_observational_reconciliation


def _model(caller, identifier, start, end):
    return GeneModel(caller=caller, identifier=identifier, start=start, end=end,
                     strand="+", sequence="M", protein_sequence="M",
                     cds_sequence="ATG", genome_id="g")


def test_automatic_observation_preserves_phanotate_and_records_relationship(tmp_path, monkeypatch):
    class Provider:
        provider_id = "prodigal_gv"
        name = "Prodigal-gv"

        def version(self):
            return "test"

        def parameters(self):
            return {"mode": "meta", "viral_only": False}

        def predict(self, *args, **kwargs):
            return GenePredictionResult(
                provider_id=self.provider_id, provider_name=self.name,
                provider_version=self.version(), models=[_model("prodigal_gv", "G1", 1, 3)],
                parameters=self.parameters(), status="SUCCESS",
            )

        def persist_raw_output(self, result, output_dir):
            Path(output_dir).mkdir(parents=True, exist_ok=True)
            p = Path(output_dir) / "prodigal_gv.tsv"
            p.write_text("raw\n")
            result.raw_output_paths = [str(p)]
            return result.raw_output_paths

    monkeypatch.setattr("phagemine.gene_callers.ProdigalGVProvider", Provider)
    primary = [_model("phanotate", "P1", 1, 3)]
    results, loci, records = run_automatic_observational_reconciliation(
        primary, "g", "ATG", tmp_path / "g.fna", tmp_path / "gene_calls")
    assert [(m.start, m.end, m.strand) for m in primary] == [(1, 3, "+")]
    assert results["prodigal_gv"].status == "SUCCESS"
    assert loci[0].reconciliation_class == "EXACT_CONCORDANCE"
    row = (tmp_path / "gene_calls" / "structural_reconciliation.tsv").read_text()
    assert "EXACT_MATCH" in row
    assert records[0]["gene_call_confidence"] == "MODERATE"


def test_pyrodigal_failure_is_nonfatal_and_unresolved(tmp_path, monkeypatch):
    class Provider:
        provider_id = "prodigal_gv"
        name = "Prodigal-gv"
        def version(self): return "unavailable"
        def parameters(self): return {"mode": "meta"}
        def predict(self, *args, **kwargs): raise RuntimeError("synthetic failure")

    monkeypatch.setattr("phagemine.gene_callers.ProdigalGVProvider", Provider)
    results, loci, records = run_automatic_observational_reconciliation(
        [_model("phanotate", "P1", 1, 3)], "g", "ATG", tmp_path / "g.fna", tmp_path / "gene_calls")
    assert results["prodigal_gv"].status == "FAILED_PROVIDER"
    assert records[0]["structural_CDS_confidence"] == "UNRESOLVED"
    manifest = json.loads((tmp_path / "gene_calls" / "structural_reconciliation.json").read_text())
    assert manifest["secondary_status"] == "FAILED_PROVIDER"
