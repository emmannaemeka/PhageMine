from phagemine.genbank import feature_table, validate
from phagemine.models import Protein


def protein(cds, *, start=1, end=None, strand="+", p5=False, p3=False):
    end = end or len(cds)
    from phagemine.genome import translate
    return Protein("g", "PM_000001", start, end, strand, cds, translate(cds), "test", partial_5prime=p5, partial_3prime=p3)


def test_complete_cds_without_stop_is_fatal():
    genome = "ATG" + "AAA" * 8 + "AAA"
    p = protein(genome)
    result = validate("g", genome, [p], {"input_sha256": "x"})
    assert not result["ncbi_submission_ready"]
    assert any(x["code"] == "NCBI_CDS_NOSTOP" for x in result["errors"])


def test_k141_reverse_boundary_nostop_regression():
    # Mirrors the rejected c1050-1 geometry: reverse-strand CDS reaches coordinate 1.
    cds = "ATG" + "AAA" * 348 + "AAA"
    genome = cds.translate(str.maketrans("ACGT", "TGCA"))[::-1]
    p = protein(cds, start=1, end=1050, strand="-")
    result = validate("k141_371497", genome, [p], {"input_sha256": "x"})
    assert any(x["code"] == "NCBI_CDS_NOSTOP" for x in result["errors"])
    assert any(x["code"] == "POSSIBLE_ORIGIN_CROSSING_CDS" for x in result["warnings"])


def test_minus_strand_three_prime_partial_feature_table():
    p = protein("ATGAAATAA", start=1, end=9, strand="-", p3=True)
    table = feature_table("g", [p])
    assert "9\t<1\tCDS" in table


def test_plus_strand_three_prime_partial_feature_table():
    p = protein("ATGAAATAA", start=1, end=9, strand="+", p3=True)
    table = feature_table("g", [p])
    assert "1\t>9\tCDS" in table
