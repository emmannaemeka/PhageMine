import json

import pytest

from phagemine.gene_callers import (
    GenePredictionResult,
    MoleculeType,
    get_gene_model_provider,
    registered_provider_ids,
)


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
