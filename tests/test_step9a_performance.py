from phagemine.models import Evidence, EvidenceLevel, Protein
from phagemine.evidence import EvidenceAdapterResult
from phagemine.pipeline import _SegmentEvidenceCache


class CountingAdapter:
    name = "counting"
    calls = 0

    def provenance(self):
        return {"database": "synthetic", "version": "1"}

    def analyze(self, proteins):
        type(self).calls += 1
        return EvidenceAdapterResult(
            adapter=self.name, status="SUCCESS_WITH_HITS",
            evidence=[Evidence("domain", "support", EvidenceLevel.COMPUTATIONAL,
                               "synthetic", "1", provenance={"protein_id": p.protein_id})
                      for p in proteins],
            provenance=self.provenance())


def _protein(pid, seq="MKK"):
    return Protein("g", pid, 1, len(seq) * 3, "+", "ATG" + "AAA" * (len(seq)-1), seq, "pyrodigal_rv")


def test_segment_cache_deduplicates_identical_sequences_and_preserves_identity():
    CountingAdapter.calls = 0
    adapter = CountingAdapter()
    cache = _SegmentEvidenceCache()
    first = cache.analyze(adapter, [_protein("A")])
    second = cache.analyze(adapter, [_protein("B")])
    assert CountingAdapter.calls == 1
    assert first.evidence[0].provenance["protein_id"] == "A"
    assert second.evidence[0].provenance["protein_id"] == "B"


def test_segment_cache_keeps_distinct_biological_proteins():
    adapter = CountingAdapter(); cache = _SegmentEvidenceCache(); CountingAdapter.calls = 0
    result = cache.analyze(adapter, [_protein("A", "MKK"), _protein("B", "MNN")])
    assert CountingAdapter.calls == 1
    assert {e.provenance["protein_id"] for e in result.evidence} == {"A", "B"}
