import json
from pathlib import Path

import pytest

from phagemine.gene_callers import (
    GenePredictionResult,
    MoleculeType,
    get_gene_model_provider,
    registered_provider_ids,
    compare_gene_model_sets,
)
from phagemine.gene_models import GeneModel


def test_registry_resolves_stable_provider_ids():
    assert {"phanotate", "prodigal"}.issubset(registered_provider_ids())
    assert get_gene_model_provider("phanotate").provider_id == "phanotate"
    assert get_gene_model_provider("prodigal").provider_id == "prodigal"


def test_unknown_provider_fails_clearly():
    with pytest.raises(ValueError, match="Unknown gene-model provider"):
        get_gene_model_provider("not-a-provider")


def test_capabilities_are_explicit_and_do_not_claim_rna():
    phanotate = get_gene_model_provider("phanotate", executable="/missing/phanotate.py")
    prodigal = get_gene_model_provider("prodigal", executable="/missing/prodigal")
    assert phanotate.supports_molecule_type(MoleculeType.DNA)
    assert not phanotate.supports_molecule_type(MoleculeType.RNA)
    assert prodigal.supports_molecule_type("dna")
    assert not prodigal.supports_molecule_type("rna")
    assert prodigal.capabilities.external_executable_required
    assert not prodigal.capabilities.native_python_provider


def test_prediction_result_serialization_is_deterministic():
    result = GenePredictionResult("phanotate", "PHANOTATE", "1.6.7", input_sequence_sha256="abc")
    assert result.to_dict()["provider_id"] == "phanotate"
    assert json.dumps(result.to_dict(), sort_keys=True) == json.dumps(result.to_dict(), sort_keys=True)


def test_legacy_predictor_api_remains_available():
    from phagemine.gene_prediction import create_predictor
    predictor = create_predictor("phanotate", "/missing/phanotate.py")
    assert predictor.name == "PHANOTATE"
    assert predictor.parameters()["executable"] == "/missing/phanotate.py"


def test_registry_resolves_pyrodigal_with_dna_only_native_capabilities():
    provider = get_gene_model_provider("pyrodigal")
    assert provider.provider_id == "pyrodigal"
    assert provider.available()
    assert provider.supports_molecule_type("dna")
    assert not provider.supports_molecule_type("rna")
    assert provider.capabilities.native_python_provider
    assert not provider.capabilities.external_executable_required


def test_registry_resolves_prodigal_gv_with_lineage_metadata():
    provider = get_gene_model_provider("prodigal_gv")
    assert provider.provider_id == "prodigal_gv"
    assert provider.name == "Prodigal-gv"
    assert provider.method_family == "prodigal_gv"
    assert provider.method_lineage == "prodigal"
    assert provider.available()
    assert provider.supports_molecule_type("dna")
    assert not provider.supports_molecule_type("rna")
    assert provider.capabilities.native_python_provider
    assert not provider.capabilities.external_executable_required


def test_pyrodigal_returns_normalized_models_and_native_raw_records(tmp_path):
    provider = get_gene_model_provider("pyrodigal")
    sequence = "C" * 1500 + "ATG" + "AAA" * 45 + "TAA" + "C" * 1500
    result = provider.predict("synthetic", sequence, tmp_path / "synthetic.fasta")
    assert isinstance(result, GenePredictionResult)
    assert result.models
    model = result.models[0]
    assert model.caller == "pyrodigal"
    assert model.coordinate_system.startswith("pyrodigal-0-based-inclusive")
    assert model.start >= 1 and model.end <= len(sequence)
    assert len(model.cds_sequence) == model.end - model.start + 1
    assert len(model.cds_sequence) % 3 == 0
    assert model.protein_sequence
    assert model.input_sequence_sha256
    paths = provider.persist_raw_output(result, tmp_path / "raw")
    assert [Path(p).name for p in paths] == ["pyrodigal.tsv", "pyrodigal.json"]
    assert "begin" in (tmp_path / "raw" / "pyrodigal.tsv").read_text()


def test_pyrodigal_rejects_rna_and_invalid_sequence():
    provider = get_gene_model_provider("pyrodigal")
    with pytest.raises(Exception, match="does not support"):
        provider.predict("rna", "A" * 120, molecule_type=MoleculeType.RNA)
    with pytest.raises(Exception, match="only A/C/G/T/N"):
        provider.predict("bad", "A" * 100 + "U")


def test_prodigal_gv_normalizes_and_preserves_native_raw_provenance(tmp_path):
    provider = get_gene_model_provider("prodigal_gv")
    sequence = "C" * 30 + "ATG" + "AAA" * 45 + "TAA" + "C" * 30
    result = provider.predict("synthetic", sequence, tmp_path / "synthetic.fasta")
    assert result.models
    assert all(model.caller == "prodigal_gv" for model in result.models)
    assert all(model.method_lineage == "prodigal" for model in result.models)
    paths = provider.persist_raw_output(result, tmp_path / "raw")
    assert [Path(p).name for p in paths] == ["prodigal_gv.tsv", "prodigal_gv.json"]
    assert "prodigal_gv" in (tmp_path / "raw" / "prodigal_gv.json").read_text()


def test_pairwise_provider_comparison_is_deterministic_and_descriptive():
    left = [GeneModel("phanotate", "p1", 100, 400, "+"), GeneModel("phanotate", "p2", 800, 900, "+")]
    right = [GeneModel("pyrodigal", "y1", 120, 400, "+"), GeneModel("pyrodigal", "y2", 1000, 1100, "+")]
    rows = compare_gene_model_sets(left, right, primary_id="phanotate", secondary_id="pyrodigal")
    assert rows[0]["comparison_class"] == "ALTERNATE_START"
    assert rows[1]["comparison_class"] == "PHANOTATE_ONLY"
    assert rows[2]["comparison_class"] == "PYRODIGAL_ONLY"
    assert rows == compare_gene_model_sets(left, right, primary_id="phanotate", secondary_id="pyrodigal")
