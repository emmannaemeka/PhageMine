from __future__ import annotations

from types import SimpleNamespace

from phagemine.models import Evidence, EvidenceLevel, Protein
from phagemine.phrogs import PHROGSPyHMMERAdapter, merge_phrogs_evidence


class _Context:
    def __init__(self, value=None):
        self.value = value or self

    def __enter__(self):
        return self.value

    def __exit__(self, *_args):
        return False


class _Hits(list):
    pass


def _protein() -> Protein:
    return Protein(
        genome_id="phage", protein_id="PM_000001", start=1, end=300,
        strand="+", cds="ATG", sequence="M" * 100, gene_call_source="test",
    )


def test_pyhmmer_adapter_parses_profile_hit_without_product_fabrication(tmp_path):
    database = tmp_path / "all_phrogs.h3m"
    database.write_bytes(b"profile fixture")
    annotations = tmp_path / "phrog_annot_v4.tsv"
    annotations.write_text(
        "phrog\tannot\tcategory\n247\tmajor head protein\thead and packaging\n")

    alignment = SimpleNamespace(
        target_from=10, target_to=90, hmm_from=3, hmm_to=83, hmm_length=100)
    domain = SimpleNamespace(alignment=alignment, i_evalue=1e-31)
    hit = SimpleNamespace(
        name=b"phrog_247", score=145.0, evalue=1e-35,
        reported=True, best_domain=domain)
    hits = _Hits([hit])
    hits.query = SimpleNamespace(name=b"PM_000001")
    module = SimpleNamespace(
        __version__="0.12.1",
        easel=SimpleNamespace(
            Alphabet=SimpleNamespace(amino=lambda: object()),
            SequenceFile=lambda *_args, **_kwargs: _Context()),
        plan7=SimpleNamespace(HMMFile=lambda *_args, **_kwargs: _Context()),
        hmmer=SimpleNamespace(hmmscan=lambda *_args, **_kwargs: [hits]),
    )
    adapter = PHROGSPyHMMERAdapter(
        database, annotations, "v4", evalue_threshold=1e-5,
        coverage_threshold=0.5, threads=4)
    adapter._module = lambda: module

    result = adapter.analyze([_protein()])

    assert result.status == "REAL"
    assert len(result.evidence) == 1
    evidence = result.evidence[0]
    assert evidence.supports is True
    assert evidence.identifier == "phrog_247"
    assert evidence.description == "major head protein"
    assert evidence.metrics["search_backend"] == "PyHMMER"
    assert evidence.metrics["query_coverage"] == 0.81
    assert evidence.provenance["protein_id"] == "PM_000001"


def test_pyhmmer_unavailable_never_fabricates_evidence(tmp_path):
    database = tmp_path / "all_phrogs.h3m"
    database.write_bytes(b"profile fixture")
    adapter = PHROGSPyHMMERAdapter(database)
    adapter._module = lambda: None

    result = adapter.analyze([_protein()])

    assert result.status == "UNAVAILABLE"
    assert result.evidence == []
    assert "not installed" in result.message


def _evidence(backend: str, supports: bool = True) -> Evidence:
    return Evidence(
        modality="phage_orthology", statement="test", level=EvidenceLevel.COMPUTATIONAL,
        source="PHROGs", source_version="v4", supports=supports,
        identifier="phrog_247", family_name="phrog_247", description="major head protein",
        evidence_strength="STRONG" if supports else "REJECTED",
        metrics={"query_protein_id": "PM_000001", "search_backend": backend,
                 "search_backends": [backend], "bit_score": 100.0, "evalue": 1e-20},
        provenance={"protein_id": "PM_000001", "search_backend": backend},
    )


def test_backend_duplicate_is_counted_once_with_corroboration():
    merged = merge_phrogs_evidence([_evidence("MMseqs2")], [_evidence("PyHMMER")])

    assert len(merged) == 1
    assert merged[0].metrics["search_backends"] == ["MMseqs2", "PyHMMER"]
    assert len(merged[0].metrics["corroborating_backend_hits"]) == 1
