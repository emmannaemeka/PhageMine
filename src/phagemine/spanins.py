"""Context-supported spanin product normalization; never changes gene boundaries.

The family allowlist is an explicit annotation policy, not a claim that every
member of a broad spanin family has experimentally established membrane topology.
"""
from __future__ import annotations
import re

RZ1_FAMILIES = {('phrogs', 'phrog_5672')}
LIPOBOX = re.compile(r'[LVI][ASTVI][GAS]C')


def canonical_spanin(product):
    text = ' '.join((product or '').split())
    if re.fullmatch(r'rz(?:-like)? spanin', text, re.I):
        return 'Rz-like spanin'
    if re.fullmatch(r'rz1(?:-like)? spanin', text, re.I):
        return 'Rz1-like spanin'
    return text


def outer_support(protein):
    motif = LIPOBOX.search(protein.sequence[:30].upper())
    support = []
    if motif:
        support.append({'kind': 'lipobox', 'motif': motif.group(),
                        'start_aa': motif.start()+1, 'end_aa': motif.end()})
    for evidence in protein.evidence:
        if evidence.status.upper() != 'REAL' or not evidence.supports:
            continue
        metrics = evidence.metrics
        value = metrics.get('evalue', metrics.get('e_value', metrics.get('full_sequence_e_value')))
        coverage = metrics.get('query_coverage')
        try:
            if value is None or float(value) > 1e-5 or coverage is None or float(coverage) < .5:
                continue
        except (TypeError, ValueError):
            continue
        text = ' '.join(str(x or '') for x in (evidence.description, evidence.family_name, evidence.statement))
        explicit = bool(re.search(r'\b(?:Rz1|o-spanin)\b|outer[- ]membrane spanin|spanin.*outer[- ]membrane', text, re.I))
        family = (evidence.source.lower(), (evidence.identifier or '').lower()) in RZ1_FAMILIES
        if explicit or family:
            support.append({'kind': 'family', 'source': evidence.source, 'identifier': evidence.identifier,
                            'evalue': float(value), 'query_coverage': float(coverage),
                            'basis': 'explicit Rz1 annotation' if explicit else 'configured Rz1 family allowlist'})
    return support


def normalize_spanins(proteins, classifications=None):
    """Update named spanins and matching classification rows, with an audit.

    Overlapping/embedded pairs must be in the +1 relative frame in transcription
    direction. Separated pairs may have any frame and gaps 0..15 nt. Candidate
    pairs must be on the same contig and strand. Ambiguous partners stay generic.
    """
    rows = {r['protein_id']: r for r in classifications or []}
    names = {p.protein_id: canonical_spanin(rows.get(p.protein_id, {}).get('display_product', p.annotation)) for p in proteins}
    events = []
    def update(p, product, event):
        old = rows.get(p.protein_id, {}).get('display_product', p.annotation)
        p.annotation = product
        if p.protein_id in rows:
            row = rows[p.protein_id]
            for key in ('display_product', 'proposed_function', 'normalized_function'):
                row[key] = product
            row['spanin_normalization'] = event
        if old != product:
            event = {**event, 'protein_id': p.protein_id, 'old_product': old, 'new_product': product}
            p.gene_call_parameters = {**p.gene_call_parameters, 'spanin_product_normalization': event}
            events.append(event)
    for p in proteins:
        name = names[p.protein_id]
        if name in {'Rz-like spanin', 'Rz1-like spanin'}:
            update(p, name, {'code': 'SPANIN_CAPITALIZATION', 'rule_version': '1.0'})
    pairs = []
    for up in proteins:
        if names[up.protein_id] not in {'Rz-like spanin', 'spanin', 'Spanin'}: continue
        for down in proteins:
            if down is up or down.genome_id != up.genome_id or down.strand != up.strand: continue
            if names[down.protein_id] not in {'Rz-like spanin', 'Rz1-like spanin', 'spanin', 'Spanin'}: continue
            u, d = (up.start, down.start) if up.strand == '+' else (-up.end, -down.end)
            if d <= u: continue
            gap = down.start-up.end-1 if up.strand == '+' else up.start-down.end-1
            frame = (d-u) % 3
            if not ((gap < 0 and frame == 1) or 0 <= gap <= 15): continue
            support = outer_support(down)
            if support:
                pairs.append((up, down, {'code': 'SPANIN_TWO_COMPONENT_SUBTYPE', 'rule_version': '1.0',
                    'upstream_id': up.protein_id, 'downstream_id': down.protein_id,
                    'gap_nt': gap, 'relative_frame': frame, 'support': support,
                    'interpretation': 'computational subtype assignment; not experimental confirmation'}))
    for up, down, event in pairs:
        # Do not resolve competing cassettes or chains by iteration order.
        if sum(up in (a,b) or down in (a,b) for a,b,_ in pairs) != 1: continue
        update(up, 'Rz-like spanin', event)
        update(down, 'Rz1-like spanin', event)
    return events
