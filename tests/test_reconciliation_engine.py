from random import Random

from phagemine.gene_models import GeneModel, ReconciledLocus
from phagemine.reconciliation_engine import reconcile_gene_models, write_reconciliation_v2


def model(caller, ident, start, end, strand="+"):
    return GeneModel(caller, ident, start, end, strand, sequence="M" * ((end - start + 1) // 3))


def test_three_caller_exact_agreement_and_provider_model_counts():
    providers = {name: [model(name, name[0], 100, 400)] for name in ("a", "b", "c")}
    locus = reconcile_gene_models(providers)[0]
    assert isinstance(locus, ReconciledLocus)
    assert locus.reconciliation_class == "EXACT_CONCORDANCE"
    assert locus.provider_count == 3 and locus.candidate_count == 3
    assert len(locus.exact_coordinate_groups) == 1


def test_common_stop_and_common_start_classes():
    assert reconcile_gene_models({"a": [model("a", "1", 100, 400)], "b": [model("b", "1", 120, 400)]})[0].reconciliation_class == "COMMON_STOP_ALTERNATE_START"
    assert reconcile_gene_models({"a": [model("a", "1", 100, 400)], "b": [model("b", "1", 100, 420)]})[0].reconciliation_class == "COMMON_START_ALTERNATE_STOP"


def test_caller_specific_strand_and_adjacent_loci():
    loci = reconcile_gene_models({"a": [model("a", "1", 100, 200), model("a", "2", 300, 400)], "b": [model("b", "3", 100, 200, "-")]})
    assert loci[0].reconciliation_class == "STRAND_DISCORDANCE"
    assert len(loci) == 2


def test_split_and_merge_are_explicit():
    split = reconcile_gene_models({"a": [model("a", "1", 100, 500)], "b": [model("b", "1", 100, 320), model("b", "2", 280, 500)]})
    assert split[0].reconciliation_class == "SPLIT_MODEL"
    merge = reconcile_gene_models({"a": [model("a", "1", 100, 320), model("a", "2", 280, 500)], "b": [model("b", "1", 100, 500)]})
    assert merge[0].reconciliation_class == "MERGED_MODEL"


def test_order_invariance_for_providers_and_models():
    base = {"a": [model("a", "1", 100, 400), model("a", "2", 600, 700)], "b": [model("b", "1", 100, 400)], "c": [model("c", "1", 600, 700)]}
    shuffled = {"c": list(reversed(base["c"])), "a": list(reversed(base["a"])), "b": list(reversed(base["b"]))}
    def signature(rows):
        return [(r.locus_id, r.reconciliation_class, [m.raw_identifier for m in r.candidate_models], r.provider_count, r.candidate_count) for r in rows]
    assert signature(reconcile_gene_models(base)) == signature(reconcile_gene_models(shuffled))


def test_nested_short_and_one_base_overlap_are_not_automatically_merged():
    loci = reconcile_gene_models({"a": [model("a", "1", 100, 200)], "b": [model("b", "2", 200, 300)]})
    assert len(loci) == 2


def test_fourth_provider_plugs_in_without_core_changes_and_serializes(tmp_path):
    providers = {name: [model(name, name, 100, 400)] for name in ("phanotate", "pyrodigal", "prodigal", "future")}
    loci = reconcile_gene_models(providers)
    assert loci[0].provider_count == 4
    write_reconciliation_v2(tmp_path, providers, loci)
    assert (tmp_path / "candidate_gene_models.tsv").is_file()
    assert (tmp_path / "reconciliation_v2.tsv").is_file()
    assert (tmp_path / "reconciliation_v2.json").is_file()
