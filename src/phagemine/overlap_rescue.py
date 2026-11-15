"""Evidence-gated overlapping ORF recovery after functional annotation.

RNA masks are absolute. A motif alone never overrides a protein overlap, and
reference coordinates are never consumed by this module.
"""
from dataclasses import replace
import json
from Bio.Seq import Seq
from .models import Protein, EvidenceLevel
from .spanins import outer_support


def family_support(protein):
    result = []
    for e in protein.evidence:
        if e.source.lower() not in {'phrogs', 'vogdb'} or e.status.upper() != 'REAL' or not e.supports:
            continue
        try:
            value = float(e.metrics.get('evalue', e.metrics.get('e_value', e.metrics.get('full_sequence_e_value'))))
            coverage = float(e.metrics['query_coverage'])
        except (TypeError, ValueError, KeyError):
            continue
        if value <= 1e-5 and .70 <= coverage <= 1 and e.identifier:
            result.append({'source': e.source.lower(), 'identifier': e.identifier,
                           'evalue': value, 'query_coverage': coverage})
    return result


def overlap_decision(candidate, row, proteins, rna_features=()):
    if any(min(candidate.end, f['end']) >= max(candidate.start, f['start']) for f in rna_features):
        return None
    overlaps = [p for p in proteins if p.genome_id == candidate.genome_id
                and min(p.end, candidate.end) >= max(p.start, candidate.start)]
    if not overlaps or len(candidate.cds) < 75:
        return None
    # Narrow overlap exception: same-strand, downstream +1-frame initiation.
    def related(p):
        delta = candidate.start-p.start if p.strand == '+' else p.end-candidate.end
        return p.strand == candidate.strand and delta > 0 and delta % 3 in (1, 2)
    if not all(related(p) for p in overlaps):
        return None
    families = family_support(candidate)
    outer_families = {(e.source.lower(), e.identifier) for p in overlaps for e in p.evidence
                      if e.status.upper() == 'REAL' and e.supports}
    distinct = [e for e in families if (e['source'], e['identifier']) not in outer_families]
    spanin = [p for p in overlaps if 'spanin' in p.annotation.lower()
              and not any(word in p.annotation.lower() for word in ('rz1', 'outer', 'unimolecular'))]
    support = outer_support(candidate)
    if len(overlaps) == 1 and spanin and row['rbs_score'] >= 8 and support:
        return {'code': 'SHORT_ORF_RZ1_OVERLAP_PROMOTED', 'product': 'Rz1-like spanin',
                'outer_ids': [p.protein_id for p in overlaps], 'support': support}
    if distinct:
        return {'code': 'SHORT_ORF_DISTINCT_FAMILY_OVERLAP_PROMOTED',
                'outer_ids': [p.protein_id for p in overlaps], 'support': distinct}
    return None


def promote_overlaps(genome, proteins, root, analyze, rna_features=()):
    """Analyze review candidates, select one start per stop, and retain a full audit."""
    from .fusion import classify_proteins
    from .start_refinement import cds_problems
    path = root / 'short_orf_candidates.json'
    if not path.exists(): return proteins
    names = {r['protein_id']: r['display_product'] for r in classify_proteins(proteins)}
    for p in proteins: p.annotation = names[p.protein_id]
    rows = json.loads(path.read_text()); candidates = []; by_id = {}; seen = set()
    for row in rows:
        key = (row['start'], row['end'], row['strand'])
        if row.get('accepted') or row.get('gap_eligible') or key in seen: continue
        if any(min(row['end'], f['end']) >= max(row['start'], f['start']) for f in rna_features): continue
        overlaps = [p for p in proteins if min(p.end, row['end']) >= max(p.start, row['start'])]
        if not overlaps or not all(p.strand == row['strand'] and
            (row['start']-p.start if p.strand == '+' else p.end-row['end']) > 0 and
            (row['start']-p.start if p.strand == '+' else p.end-row['end']) % 3 in (1, 2) for p in overlaps): continue
        row['overlap_policy_priority'] = 0 if all((row['start']-p.start if p.strand=='+' else p.end-row['end'])%3==1 for p in overlaps) else 1
        seen.add(key)
        cds = genome[key[0]-1:key[1]]
        if key[2] == '-': cds = str(Seq(cds).reverse_complement())
        p = Protein(proteins[0].genome_id, f'OVERLAP_{len(candidates)+1}', *key, cds,
                    str(Seq(cds).translate(table=11, cds=True)), 'short_orf_rescue',
                    start_codon=cds[:3], stop_codon=cds[-3:])
        if cds_problems(p, genome): continue
        candidates.append(p); by_id[p.protein_id] = (p, row)
    import time
    started = time.monotonic()
    statuses = []
    def search(batch, stage):
        if not batch: return
        for offset in range(0, len(batch), 128):
            for result in analyze(batch[offset:offset+128]):
                statuses.append({'adapter': result.adapter, 'status': result.status,
                                 'stage': stage, 'queries': len(batch[offset:offset+128]), 'provenance': result.provenance})
                for e in result.evidence:
                    pid = e.provenance.get('protein_id') or e.metrics.get('query_protein_id')
                    if pid in by_id: by_id[pid][0].evidence.append(e)
    # Group nested starts by stop. Search the longest peptide first so a
    # family-supported, SD-free ORF can qualify without reference coordinates.
    groups = {}
    direct = set()
    for p in candidates:
        row = by_id[p.protein_id][1]
        overlaps = [q for q in proteins if min(q.end,p.end)>=max(q.start,p.start)]
        context = any('spanin' in q.annotation.lower() for q in overlaps)
        if row['rbs_score'] > 0 or context: direct.add(p.protein_id)
        key = (p.strand, p.end if p.strand=='+' else p.start)
        groups.setdefault(key, []).append(p)
    seeds = [max(group, key=lambda p:len(p.sequence)) for group in groups.values()
             if any(p.protein_id not in direct for p in group)]
    search(seeds, 'family_prefilter')
    seeded = {p.protein_id for p in seeds}
    family_stops = set()
    for key, group in groups.items():
        for p in group:
            for e in p.evidence:
                try:
                    value=float(e.metrics.get('evalue',e.metrics.get('e_value',e.metrics.get('full_sequence_e_value'))))
                except (ValueError,TypeError): continue
                if e.status.upper()=='REAL' and e.source.lower() in {'phrogs','vogdb'} and value<=1e-3:
                    family_stops.add(key)
    eligible=[]
    for p in candidates:
        key=(p.strand,p.end if p.strand=='+' else p.start)
        row=by_id[p.protein_id][1]
        row['overlap_prefilter']='RBS_OR_SPANIN' if p.protein_id in direct else ('FAMILY_EVIDENCE_EXCEPTION' if key in family_stops else 'NO_INITIATION_OR_FAMILY_SUPPORT')
        if p.protein_id in direct or key in family_stops: eligible.append(p)
    confirmation_queries = sum(p.protein_id not in seeded for p in eligible)
    search([p for p in eligible if p.protein_id not in seeded], 'candidate_confirmation')
    allowed={p.protein_id for p in eligible}
    output = list(proteins); events = []; considered = []
    used = {(p.strand, p.end if p.strand == '+' else p.start) for p in proteins}
    serial = max(int(p.protein_id[3:]) for p in proteins if p.protein_id.startswith('PM_'))
    for p, row in sorted(by_id.values(), key=lambda pair: (pair[1].get('overlap_policy_priority',0), -pair[1]['score'])):
        decision = overlap_decision(p, row, output, rna_features) if p.protein_id in allowed else None
        considered.append({'start': p.start, 'end': p.end, 'strand': p.strand,
                           'rbs_score': row['rbs_score'], 'families': family_support(p), 'decision': decision})
        stop = (p.strand, p.end if p.strand == '+' else p.start)
        if not decision or stop in used: continue
        serial += 1; p.protein_id = p.locus_tag = f'PM_{serial:06d}'
        p.gene_call_parameters['overlap_rescue'] = decision
        if decision.get('product'):
            p.annotation = decision['product']; p.annotation_level = EvidenceLevel.COMPUTATIONAL
        output.append(p); used.add(stop)
        row.update(accepted=True, review_reason=decision['code'])
        events.append({'protein_id': p.protein_id, 'start': p.start, 'end': p.end, **decision})
    path.write_text(json.dumps(rows, indent=2)+'\n')
    (root/'overlap_rescue.json').write_text(json.dumps({'elapsed_seconds': time.monotonic()-started, 'prefilter': {'total_candidates': len(candidates), 'family_seed_queries': len(seeds), 'eligible_candidates': len(eligible), 'confirmation_queries': confirmation_queries}, 'adapters': statuses, 'candidates': considered, 'promoted': events}, indent=2)+'\n')
    return output
