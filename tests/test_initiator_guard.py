"""Unit tests for phagemine.initiator_guard (synthetic sequences; no reference data)."""
from Bio.Seq import Seq

from phagemine.models import Protein
from phagemine.initiator_guard import enforce_canonical_initiators, RULE

PAD = 'C' * 30  # SD-free filler


def _plus_genome():
    # ... ATT (rare, upstream) ... 3 codons ... ATG (canonical, Prodigal) ... body ... TAA
    body = 'GCT' * 40
    orf_rare = 'ATT' + 'GCA' * 3 + 'ATG' + body + 'TAA'
    genome = PAD + 'AGGAGGTCACCC' + orf_rare + PAD
    rare_start = len(PAD) + 12 + 1
    canon_start = rare_start + 12
    stop_end = rare_start + len(orf_rare) - 1
    return genome, rare_start, canon_start, stop_end


def _cds(genome, start, end, strand='+'):
    s = genome[start - 1:end]
    return s if strand == '+' else str(Seq(s).reverse_complement())


def _protein(genome, pid, start, end, strand='+', **kw):
    cds = _cds(genome, start, end, strand)
    return Protein('g', pid, start, end, strand, cds, str(Seq(cds).translate(table=11)).rstrip('*'),
                   'PHANOTATE', **kw)


class Model:
    """Minimal stand-in for gene_models.GeneModel (caller, identifier, start, end, strand)."""
    def __init__(self, caller, identifier, start, end, strand):
        self.caller, self.identifier, self.start, self.end, self.strand = caller, identifier, start, end, strand


def test_rare_start_replaced_by_unique_canonical_caller_start():
    genome, rare, canon, stop = _plus_genome()
    p = _protein(genome, 'PM_000001', rare, stop)
    secondary = [Model('prodigal_gv', 'PG1', canon, stop, '+')]
    out, changes, reviews = enforce_canonical_initiators(genome, [p], secondary)
    assert out[0].start == canon and out[0].end == stop
    assert out[0].cds.startswith('ATG') and out[0].sequence.startswith('M')
    assert changes[0]['rule'] == RULE and changes[0]['old_start_codon'] == 'ATT'
    assert not reviews


def test_canonical_start_never_touched():
    genome, rare, canon, stop = _plus_genome()
    p = _protein(genome, 'PM_000001', canon, stop)
    secondary = [Model('prodigal_gv', 'PG1', rare, stop, '+')]
    out, changes, _ = enforce_canonical_initiators(genome, [p], secondary)
    assert out[0] is p and not changes


def test_homology_anchored_rare_start_kept():
    genome, rare, canon, stop = _plus_genome()
    params = {'start_refinement': {'code': 'START_REFINED_NTERMINAL_HOMOLOGY'}}
    p = _protein(genome, 'PM_000001', rare, stop, gene_call_parameters=params)
    secondary = [Model('prodigal_gv', 'PG1', canon, stop, '+')]
    out, changes, reviews = enforce_canonical_initiators(genome, [p], secondary)
    assert out[0] is p and not changes
    assert reviews[0]['decision'] == 'KEPT_INDEPENDENT_EVIDENCE'


def test_overlapping_architecture_flagged_not_changed():
    genome, rare, canon, stop = _plus_genome()
    p = _protein(genome, 'PM_000001', rare, stop)
    # same-strand partner overlapping >= 100 nt (e.g. Rz/Rz1-like arrangement)
    partner = Protein('g', 'PM_000002', rare + 3, stop - 2, '+', '', '', 'PHANOTATE')
    secondary = [Model('prodigal_gv', 'PG1', canon, stop, '+')]
    out, changes, reviews = enforce_canonical_initiators(genome, [p, partner], secondary)
    assert out[0] is p and not changes
    assert reviews[0]['decision'] == 'KEPT_OVERLAPPING_ARCHITECTURE_REVIEW'


def test_no_caller_start_kept_for_review():
    genome, rare, canon, stop = _plus_genome()
    p = _protein(genome, 'PM_000001', rare, stop)
    out, changes, reviews = enforce_canonical_initiators(genome, [p], [])
    assert out[0] is p and not changes
    assert reviews[0]['decision'] == 'KEPT_NO_UNIQUE_CANONICAL_CALLER_START'


def test_different_stop_or_other_caller_ignored():
    genome, rare, canon, stop = _plus_genome()
    p = _protein(genome, 'PM_000001', rare, stop)
    secondary = [Model('prodigal_gv', 'PG1', canon, stop - 3, '+'),
                 Model('phanotate', 'PH1', canon, stop, '+')]
    out, changes, _ = enforce_canonical_initiators(genome, [p], secondary)
    assert out[0] is p and not changes


def test_minus_strand_replacement():
    genome, rare, canon, stop = _plus_genome()
    rc = str(Seq(genome).reverse_complement())
    L = len(genome)
    # map plus-strand coordinates onto the reverse-complement genome
    r_start, r_end = L - stop + 1, L - rare + 1
    c_end = L - canon + 1
    p = _protein(rc, 'PM_000001', r_start, r_end, '-')
    secondary = [Model('prodigal_gv', 'PG1', r_start, c_end, '-')]
    out, changes, _ = enforce_canonical_initiators(rc, [p], secondary)
    assert (out[0].start, out[0].end) == (r_start, c_end)
    assert out[0].cds.startswith('ATG') and out[0].stop_codon == 'TAA'


def test_partial_cds_untouched():
    genome, rare, canon, stop = _plus_genome()
    p = _protein(genome, 'PM_000001', rare, stop, partial_5prime=True)
    secondary = [Model('prodigal_gv', 'PG1', canon, stop, '+')]
    out, changes, reviews = enforce_canonical_initiators(genome, [p], secondary)
    assert out[0] is p and not changes and not reviews


def test_pipeline_invokes_guard_after_start_refinement():
    from pathlib import Path
    import phagemine.pipeline as pipeline
    text = Path(pipeline.__file__).read_text()
    refine = text.index('proteins = run_refinement(')
    guard = text.index('enforce_canonical_initiators(')
    pfam = text.index('PfamHMMAdapter(pfam_path')
    assert refine < guard < pfam  # evidence is searched on guarded models
    assert 'initiator_guard.json' in text
