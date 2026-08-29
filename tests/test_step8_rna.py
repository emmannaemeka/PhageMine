import pytest

from phagemine.gene_callers import get_gene_model_provider, MoleculeType, UnsupportedMoleculeType
from phagemine.genome_representation import GenomeRecord, SegmentRecord


def test_pyrodigal_rv_is_registered_and_rna_only():
    provider = get_gene_model_provider("pyrodigal_rv")
    assert provider.provider_id == "pyrodigal_rv"
    assert provider.supports_molecule_type(MoleculeType.RNA)
    assert not provider.supports_molecule_type(MoleculeType.DNA)
    assert provider.capabilities.native_python_provider
    assert provider.capabilities.external_executable_required is False


def test_pyrodigal_rv_unavailable_is_explicit_when_dependency_missing():
    provider = get_gene_model_provider("pyrodigal_rv")
    if not provider.available():
        with pytest.raises(UnsupportedMoleculeType):
            provider.predict("rna", "ATG" * 40, molecule_type="dna")


def test_dna_providers_do_not_advertise_rna():
    for name in ("phanotate", "pyrodigal", "prodigal_gv"):
        assert not get_gene_model_provider(name).supports_molecule_type(MoleculeType.RNA)


def test_segmented_genome_preserves_segment_local_sequences():
    genome = GenomeRecord("phi6", "rna", (SegmentRecord("L", "A" * 12), SegmentRecord("M", "C" * 9), SegmentRecord("S", "G" * 6)))
    assert [segment.segment_id for segment in genome.segments] == ["L", "M", "S"]
    assert [segment.length for segment in genome.segments] == [12, 9, 6]
    assert len({segment.sha256 for segment in genome.segments}) == 3
