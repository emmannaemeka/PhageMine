"""Canonical-initiator guard for final CDS starts (reference-independent).

Rationale: start refinement scores a rare table-11 initiator (CTG/ATT/ATC/ATA,
codon prior 0.5) against canonical ATG/GTG/TTG mainly through SD-like motif and
spacing terms worth up to 12 points, so a chance purine motif can make a rare
initiator outrank a canonical, caller-supported start sharing the same stop.

This guard runs after all start-changing stages. A CDS that begins with a rare
initiator is switched to the unique canonical start proposed by the secondary
ab initio caller (Prodigal-GV) for the same strand and stop, unless:
  * the rare start is anchored by N-terminal homology or sequence-confirmed
    translational coupling (independent initiation evidence), or
  * the CDS overlaps a same-strand CDS by >= ``max_same_strand_overlap`` nt
    (overlapping/embedded architectures such as Rz/Rz1 spanins), which is
    flagged for review instead of being changed automatically.
Stops, frames, IDs and products are never changed. Reference annotations are
never consulted.
"""
from __future__ import annotations

from dataclasses import replace
from Bio.Seq import Seq

from .start_selection import STOPS

CANONICAL_INITIATORS = frozenset({'ATG', 'GTG', 'TTG'})
RARE_INITIATORS = frozenset({'CTG', 'ATT', 'ATC', 'ATA'})
PROTECTED_REFINEMENT_REASONS = frozenset({
    'START_REFINED_NTERMINAL_HOMOLOGY',
    'START_REFINED_TRANSLATIONAL_COUPLING',
})
SECONDARY_CALLERS = frozenset({'prodigal_gv', 'pyrodigal_gv'})
RULE = 'RARE_INITIATOR_REPLACED_BY_CANONICAL_CALLER_START'


def _oriented(genome, start, end, strand):
    seq = genome[start - 1:end]
    return seq if strand == '+' else str(Seq(seq).reverse_complement())


def _stop_coordinate(model):
    return model.end if model.strand == '+' else model.start


def _start_coordinate(model):
    return model.start if model.strand == '+' else model.end


def _protected(p):
    params = getattr(p, 'gene_call_parameters', {}) or {}
    event = params.get('start_refinement') or {}
    if event.get('code') in PROTECTED_REFINEMENT_REASONS:
        return event['code']
    components = event.get('score_components') or {}
    if components.get('homology_score'):
        return 'START_REFINED_NTERMINAL_HOMOLOGY'
    return None


def _max_same_strand_overlap(p, proteins):
    best = 0
    for q in proteins:
        if q is p or q.protein_id == p.protein_id or q.strand != p.strand:
            continue
        best = max(best, min(p.end, q.end) - max(p.start, q.start) + 1)
    return best


def _is_secondary_caller(model):
    caller = getattr(model, 'caller', None)
    if not caller:
        return True  # plain coordinate models supplied by the caller stage
    return caller.lower().replace('-', '_') in SECONDARY_CALLERS


def _valid_orf(cds):
    if len(cds) % 3 or cds[:3] not in CANONICAL_INITIATORS or cds[-3:] not in STOPS:
        return False
    return not any(cds[i:i + 3] in STOPS for i in range(0, len(cds) - 3, 3))


def enforce_canonical_initiators(genome, proteins, secondary=(), max_same_strand_overlap=100):
    """Return (proteins, changes, reviews). Never alters partial CDSs or stops."""
    genome = genome.upper()
    output, changes, reviews = [], [], []
    for p in proteins:
        output.append(p)
        if p.strand not in {'+', '-'} or p.partial_5prime or p.partial_3prime:
            continue
        if not 1 <= p.start < p.end <= len(genome):
            continue
        cds = _oriented(genome, p.start, p.end, p.strand)
        codon = cds[:3]
        if codon not in RARE_INITIATORS:
            continue
        reason = _protected(p)
        if reason:
            reviews.append({'protein_id': p.protein_id, 'start_codon': codon,
                            'decision': 'KEPT_INDEPENDENT_EVIDENCE', 'evidence': reason})
            continue
        overlap = _max_same_strand_overlap(p, proteins)
        if overlap >= max_same_strand_overlap:
            reviews.append({'protein_id': p.protein_id, 'start_codon': codon,
                            'decision': 'KEPT_OVERLAPPING_ARCHITECTURE_REVIEW',
                            'same_strand_overlap_nt': overlap})
            continue
        stop = _stop_coordinate(p)
        options = {}
        for c in secondary:
            if c.strand != p.strand or _stop_coordinate(c) != stop or not _is_secondary_caller(c):
                continue
            new_cds = _oriented(genome, c.start, c.end, c.strand)
            if _valid_orf(new_cds) and (len(new_cds) - len(cds)) % 3 == 0:
                options[_start_coordinate(c)] = (c, new_cds)
        if len(options) != 1:
            reviews.append({'protein_id': p.protein_id, 'start_codon': codon,
                            'decision': 'KEPT_NO_UNIQUE_CANONICAL_CALLER_START',
                            'options': sorted(options)})
            continue
        c, new_cds = next(iter(options.values()))
        translation = str(Seq(new_cds).translate(table=11, cds=True))
        event = {'rule': RULE, 'protein_id': p.protein_id, 'old_start': p.start, 'old_end': p.end,
                 'new_start': c.start, 'new_end': c.end, 'strand': p.strand,
                 'old_start_codon': codon, 'new_start_codon': new_cds[:3],
                 'distance_nt': (c.start - p.start) if p.strand == '+' else (p.end - c.end),
                 'same_strand_overlap_nt': overlap}
        params = dict(getattr(p, 'gene_call_parameters', {}) or {})
        params['initiator_guard'] = event
        params['start_selection_review'] = 'rare initiator replaced by canonical caller start - review'
        output[-1] = replace(p, start=c.start, end=c.end, cds=new_cds, sequence=translation,
                             start_codon=new_cds[:3], stop_codon=new_cds[-3:],
                             gene_call_parameters=params)
        changes.append(event)
    return output, changes, reviews
