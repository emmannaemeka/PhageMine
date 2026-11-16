"""Reference-independent CDS start scoring and conservative short-ORF rescue.

Scores are transparent heuristics, not probabilities. All coordinates in public
records are 1-based inclusive. Work internally in coding-strand orientation.
"""
from __future__ import annotations

import csv
import json
import subprocess
from dataclasses import dataclass, asdict, replace
from pathlib import Path
from Bio.Seq import Seq
from .models import Protein
from .start_selection import STARTS_11, STOPS


@dataclass(frozen=True)
class RefinementConfig:
    min_cds_nt: int = 90  # includes the terminal stop (90 nt encodes 29 residues)
    rescue_min_nt: int = 75
    rescue_max_nt: int = 300
    change_margin: float = 2.0
    rbs_weight: float = 1.0
    gap_penalty: float = 0.5
    large_overlap_penalty: float = 0.15
    caller_bonus: float = 1.5
    rescue_min_score: float = 8.0
    min_uncovered_nt: int = 30


CODON_PRIOR = {'ATG': 3.0, 'GTG': 2.0, 'TTG': 1.5,
               'ATT': 0.5, 'ATA': 0.5, 'ATC': 0.5, 'CTG': 0.5}


def oriented_bounds(p, length):
    return (p.start - 1, p.end) if p.strand == '+' else (length - p.end, length - p.start + 1)


def genomic_bounds(left, right, strand, length):
    return (left + 1, right) if strand == '+' else (length - right + 1, length - left)


def rbs_score(sequence, left):
    """Best SD-like motif wholly within -20..-4, spacer measured from motif end."""
    best = {'score': 0.0, 'motif': '', 'spacer_nt': None, 'matches': 0}
    consensus = 'AGGAGG'
    for width in (4, 5, 6):
        patterns = {consensus[i:i + width] for i in range(7 - width)}
        for i in range(max(0, left - 20), left - 4 - width + 1):
            word = sequence[i:i + width]
            if any(b not in 'ACGT' for b in word):
                continue
            matches = max(sum(a == b for a, b in zip(word, motif)) for motif in patterns)
            if word in patterns:
                score = 2.0 * width
            elif width == 6 and matches == 5 and sum(b in 'AG' for b in word) >= 5:
                score = 6.0
            else:
                continue
            spacer = left - (i + width)
            if not 5 <= spacer <= 12: score *= 0.6
            if score > best['score']:
                best = {'score': score, 'motif': word, 'spacer_nt': spacer, 'matches': matches}
    return best


def enumerate_starts(sequence, right, min_cds_nt=90):
    """Enumerate every in-frame start back to the preceding stop/gap; no wrapping."""
    if not 3 <= right <= len(sequence) or sequence[right-3:right] not in STOPS:
        return []
    result = []
    for left in range(right - 6, -1, -3):
        codon = sequence[left:left+3]
        if codon in STOPS or any(b not in 'ACGT' for b in codon): break
        if right - left >= min_cds_nt and codon in STARTS_11: result.append(left)
    return sorted(result)


def cds_problems(p, genome):
    """Strict sequence-level checks, including stops hidden by a supplied translation."""
    if p.strand not in {'+', '-'} or not 1 <= p.start <= p.end <= len(genome):
        return ['CDS_COORDINATE_INVALID']
    seq = genome[p.start-1:p.end].upper()
    if p.strand == '-': seq = str(Seq(seq).reverse_complement())
    errors = []
    if seq != p.cds.upper(): errors.append('CDS_SEQUENCE_MISMATCH')
    offset = p.codon_start - 1 if p.partial_5prime else 0
    if offset not in (0, 1, 2): return errors + ['CDS_CODON_START_INVALID']
    coding = seq[offset:]
    if not p.partial_5prime and coding[:3] not in STARTS_11: errors.append('CDS_START_INVALID')
    if not p.partial_3prime and coding[-3:] not in STOPS: errors.append('CDS_STOP_MISSING')
    if len(coding) % 3 and not p.partial_3prime: errors.append('CDS_FRAME_INVALID')
    usable = coding[:len(coding)//3*3]
    translation = str(Seq(usable).translate(table=11))
    body = translation[:-1] if not p.partial_3prime and translation.endswith('*') else translation
    if '*' in body: errors.append('CDS_INTERNAL_STOP_REVIEW_RECODING')
    if not p.partial_5prime and coding[:3] in STARTS_11 and body: body = 'M' + body[1:]
    if body != p.sequence: errors.append('CDS_TRANSLATION_MISMATCH')
    return errors


def _spacing(left, right, strand, proteins, length, exclude):
    preceding = []
    for p in proteins:
        if p.protein_id == exclude or p.strand != strand: continue
        a, b = oriented_bounds(p, length)
        if a < left and b < right: preceding.append((a, b))
    if not preceding: return None
    # Most recently initiated CDS, including overlaps; do not use the candidate itself.
    _, end = max(preceding)
    return left - end


def genome_start_context(genome, proteins, secondary=()):
    """Calibrate only from caller-agreement anchors, never reference annotations.

    SD absence is not evidence of leaderless transcription. Sparse anchor sets
    leave SD weighting unchanged. Moderate-GC genomes retain the original rule.
    """
    genome = genome.upper()
    gc = sum(b in 'GC' for b in genome) / max(1, sum(b in 'ACGT' for b in genome))
    agreement = {(p.start, p.end, p.strand) for p in secondary}
    anchors = [p for p in proteins if (p.start, p.end, p.strand) in agreement
               and not p.partial_5prime and not p.partial_3prime
               and 'start_selection' not in p.gene_call_parameters]
    hits = 0
    for p in anchors:
        seq = genome if p.strand == '+' else str(Seq(genome).reverse_complement())
        hits += rbs_score(seq, oriented_bounds(p, len(genome))[0])['score'] >= 8
    prevalence = hits / len(anchors) if anchors else None
    t = max(0.0, min(1.0, (gc - .58) / .12))
    adaptation = t*t*(3-2*t)  # smoothstep, no jump at 58% GC
    sparse_sd = len(anchors) >= 10 and prevalence < .35
    return {'gc_fraction': gc, 'anchor_count': len(anchors), 'sd_prevalence': prevalence,
            'adaptation': adaptation, 'sd_weight': 0.0 if gc > .58 and sparse_sd else 1.0,
            'sd_sparse': bool(gc > .58 and sparse_sd)}


def attenuate_rbs(rbs, context):
    """Attenuate short/background motifs only in high-GC genomes."""
    result = dict(rbs)
    result['raw_score'] = rbs['score']
    if context and context['gc_fraction'] > .58 and len(rbs['motif']) < 5:
        result['score'] *= .25
    return result


def coupling_score(sequence, left, gap):
    """Require an actual shared stop/start sequence, not just close coordinates."""
    if gap == -1 and sequence[max(0,left-2):left+3] in {'TAATG','TGATG','TAGTG','TATTG'}:
        return 4.0
    if gap == -4 and sequence[left:left+4] in {'ATGA','GTGA','TTGA'}:
        return 4.0
    return 0.0


def candidate_scores(genome, proteins, p, secondary=(), anchors=None, config=None, context=None):
    config = config or RefinementConfig()
    length = len(genome)
    sequence = genome if p.strand == '+' else str(Seq(genome).reverse_complement())
    original_left, right = oriented_bounds(p, length)
    caller_starts = {original_left}
    for c in secondary:
        if c.strand == p.strand:
            a, b = oriented_bounds(c, length)
            if b == right: caller_starts.add(a)
    rows = []
    for left in enumerate_starts(sequence, right, config.min_cds_nt):
        # Attenuation targets alternative background motifs, not caller-supported starts.
        rbs = attenuate_rbs(rbs_score(sequence, left), context if left not in caller_starts else None)
        gap = _spacing(left, right, p.strand, proteins, length, p.protein_id)
        packing = 0.0
        coupling = coupling_score(sequence, left, gap)
        if gap in (-1, -4):
            packing = coupling
            if coupling and rbs['raw_score'] > rbs['score']:
                packing += 2.0  # independent sequence-confirmed coupling evidence
        elif gap is not None and 0 <= gap <= 15: packing = 3.0
        elif gap is not None and gap > 15: packing = -min(20.0, config.gap_penalty*(gap-15))
        elif gap is not None and gap < -20:
            packing = -min(20.0, config.large_overlap_penalty*(-gap-20))
            if rbs['score'] >= 8: packing *= 0.25
        start, end = genomic_bounds(left, right, p.strand, length)
        key = (start, end, p.strand)
        homology = (anchors or {}).get(key)
        homology_score = 16.0 if homology else 0.0
        codon = sequence[left:left+3]
        prior = CODON_PRIOR[codon]
        if context and codon in {'GTG', 'TTG'}:
            prior += context['adaptation'] * (3.0 - prior) * .8
        sd_weight = (context or {}).get('sd_weight', 1.0)
        support = config.caller_bonus if left in caller_starts else 0.0
        score = prior + config.rbs_weight*sd_weight*rbs['score'] + packing + support + homology_score
        rows.append({'protein_id': p.protein_id, 'start': start, 'end': end, 'strand': p.strand,
                     'left': left, 'right': right, 'start_codon': sequence[left:left+3],
                     'length_nt': right-left, 'score': round(score, 4), 'codon_prior': prior,
                     'coupling_score': coupling, 'rbs_raw_score': rbs['raw_score'], 'rbs_weight': sd_weight, 'rbs_score': rbs['score'], 'rbs_motif': rbs['motif'], 'spacer_nt': rbs['spacer_nt'],
                     'upstream_gap_nt': gap, 'packing_score': packing, 'caller_score': support,
                     'homology_score': homology_score, 'homology': homology,
                     'original': left == original_left})
    return sorted(rows, key=lambda r: (-r['score'], not r['original'], r['left']))


def _changed_model(p, row, genome, reason):
    cds = genome[row['start']-1:row['end']]
    if p.strand == '-': cds = str(Seq(cds).reverse_complement())
    translation = str(Seq(cds).translate(table=11, cds=True))
    event = {'code': reason, 'protein_id': p.protein_id, 'old_start': p.start, 'old_end': p.end,
             'new_start': row['start'], 'new_end': row['end'], 'strand': p.strand,
             'distance_nt': (row['start']-p.start if p.strand=='+' else p.end-row['end']),
             'score_components': row}
    params = dict(p.gene_call_parameters)
    params['start_refinement'] = event
    q = replace(p, start=row['start'], end=row['end'], cds=cds, sequence=translation,
                start_codon=cds[:3], stop_codon=cds[-3:], partial_5prime=False,
                gene_call_parameters=params)
    return q, event


def refine_starts(genome, proteins, secondary=(), anchors=None, config=None):
    """Use an immutable neighbor snapshot, retaining ambiguous/unsupported changes."""
    genome = genome.upper(); config = config or RefinementConfig()
    output, audit, candidates = [], [], []
    context = genome_start_context(genome, proteins, secondary)
    for p in proteins:
        output.append(p)
        if p.strand not in {'+', '-'} or p.partial_5prime or p.partial_3prime: continue
        if not 1 <= p.start < p.end <= len(genome): continue
        rows = candidate_scores(genome, proteins, p, secondary, anchors, config, context)
        candidates.extend(rows)
        if not rows: continue
        best = rows[0]; original = next((r for r in rows if r['original']), None)
        if original and best['original']: continue
        # No long truncation solely from a chance purine motif: caller or N-terminal
        # alignment support is required for loss of >20% of the existing CDS.
        if best['length_nt'] < 0.8 * len(p.cds) and not (best['caller_score'] or best['homology_score']): continue
        margin = best['score'] - (original['score'] if original else -100)
        if margin < config.change_margin: continue
        # Protect a caller-supported translational coupling boundary from a
        # chance internal RBS unless an alignment independently anchors it.
        if original and original['upstream_gap_nt'] in (-1, -4) and best['upstream_gap_nt'] not in (-1, -4) and not best['homology_score']:
            continue
        compact_extension = bool(original and 0 <= (best['upstream_gap_nt'] if best['upstream_gap_nt'] is not None else -1) <= 30
                                 and 0 < original['left']-best['left'] <= 30
                                 and original['rbs_score'] == 0
                                 and original['upstream_gap_nt'] is not None
                                 and original['upstream_gap_nt']-best['upstream_gap_nt'] >= 15)
        if not ((not context['sd_sparse'] and best['rbs_score'] >= 8) or best['homology_score'] or (best.get('coupling_score',0) and best.get('rbs_raw_score',0)>=8) or compact_extension): continue
        reason = 'START_REFINED_NTERMINAL_HOMOLOGY' if best['homology_score'] else ('START_REFINED_COMPACT_SPACING_REVIEW' if compact_extension and best['rbs_score'] < 8 else 'START_REFINED_RBS_SPACING')
        if best.get('coupling_score',0) and best['rbs_score'] < 8 and not best['homology_score']:
            reason = 'START_REFINED_TRANSLATIONAL_COUPLING'
        output[-1], event = _changed_model(p, best, genome, reason)
        output[-1].gene_call_parameters['start_selection_review'] = 'heuristic start refinement - review'
        audit.append(event)
    return output, audit, candidates


def intergenic_candidate(start, end, proteins):
    """Only empty gaps, with at most four bases overlapping either boundary.

    No internal or antisense CDS is automatically added. Both strands occupy
    the same genomic space; homology does not override this safety constraint.
    """
    for p in proteins:
        overlap = min(end, p.end) - max(start, p.start) + 1
        if overlap <= 0:
            continue
        if overlap > 4 or (start <= p.start and end >= p.end) or (p.start <= start and p.end >= end):
            return False
    return True


def rescue_short_orfs(genome, proteins, config=None, anchors=None, rna_features=()):
    """Recover intergenic 75–300 nt stop-delimited ORFs with initiation evidence.

    Existing IDs never change. Unsupported candidates remain in the audit table
    rather than being asserted as additional genes merely because they are open.
    """
    config = config or RefinementConfig(); config = replace(config, min_cds_nt=config.rescue_min_nt)
    genome = genome.upper(); length = len(genome)
    occupied = bytearray(length)
    for p in proteins:
        occupied[p.start-1:p.end] = b'\1' * (p.end-p.start+1)
    used_stops = {(p.strand, oriented_bounds(p,length)[1]) for p in proteins}
    output = list(proteins); audit = []; candidates = []
    serial = max([int(p.protein_id[3:]) for p in proteins if p.protein_id.startswith('PM_') and p.protein_id[3:].isdigit()] or [0])
    for strand in ('+', '-'):
        sequence = genome if strand == '+' else str(Seq(genome).reverse_complement())
        for stop_left in range(length-2):
            right = stop_left+3
            if sequence[stop_left:right] not in STOPS or (strand,right) in used_stops: continue
            starts = [a for a in enumerate_starts(sequence,right,config.min_cds_nt) if right-a<=config.rescue_max_nt]
            if not starts: continue
            a,b=genomic_bounds(starts[0],right,strand,length)
            dummy=Protein(proteins[0].genome_id if proteins else 'genome','candidate',a,b,strand,'','', 'short_orf_rescue')
            rows=candidate_scores(genome, proteins, dummy, anchors=anchors, config=config)
            # A synthetic placeholder must not acquire a caller-support bonus.
            for row in rows:
                row['score'] -= row['caller_score']; row['caller_score']=0.0
            rows=[r for r in rows if r['length_nt']<=config.rescue_max_nt]
            rows.sort(key=lambda r:(-r['score'],r['start']))
            for row in rows:
                uncovered = sum(v==0 for v in occupied[row['start']-1:row['end']])
                row['uncovered_nt']=uncovered
                row['gap_eligible']=intergenic_candidate(row['start'], row['end'], output)
                row['rna_overlap'] = any(min(row['end'], f['end']) >= max(row['start'], f['start']) for f in rna_features)
                row['accepted']=False
                row['review_reason']='SHORT_ORF_INSUFFICIENT_INITIATION_EVIDENCE' if row['gap_eligible'] else 'SHORT_ORF_OVERLAPS_EXISTING_CDS_REVIEW'
            for row in rows:
                if row['rna_overlap']: row['review_reason'] = 'SHORT_ORF_OVERLAPS_RNA_REJECTED'
            candidates.extend(rows)
            eligible=[r for r in rows if not r['rna_overlap'] and r['gap_eligible'] and r['uncovered_nt']>=config.min_uncovered_nt
                      and (r['homology_score'] or
                           (r['start_codon'] in STARTS_11 and
                            ((r['score']>=config.rescue_min_score and r['rbs_score']>=8) or
                             (r['upstream_gap_nt'] in (-1,-4) and r['score']>=7))))]
            if not eligible: continue
            best=eligible[0];serial+=1
            while any(p.protein_id==f'PM_{serial:06d}' for p in output):serial+=1
            dummy=replace(dummy,protein_id=f'PM_{serial:06d}',locus_tag=f'PM_{serial:06d}',start=best['start'],end=best['end'])
            q,event=_changed_model(dummy,best,genome,'SHORT_ORF_RESCUED_REVIEW')
            q.gene_call_parameters['start_selection_review']='short ORF rescued - review'
            output.append(q);audit.append(event);best['accepted']=True;best['review_reason']='SHORT_ORF_RESCUED_REVIEW'
            occupied[q.start-1:q.end]=b'\1'*(q.end-q.start+1)
    return output,audit,candidates


def diamond_anchors(genome, candidate_rows, database, executable, output, threads=1):
    """Optional full-length protein N-terminal support; domain hits are insufficient."""
    root=Path(output);root.mkdir(parents=True,exist_ok=True)
    query=root/'candidate_starts.faa'; keys={}
    with query.open('w') as f:
        seen=set()
        for row in candidate_rows:
            key=(row['start'],row['end'],row['strand'])
            if key in seen:continue
            seen.add(key);name=f'C{len(keys)+1}';keys[name]=key
            seq=genome[key[0]-1:key[1]]
            if key[2]=='-':seq=str(Seq(seq).reverse_complement())
            f.write(f'>{name}\n{Seq(seq).translate(table=11,cds=True)}\n')
    hits=root/'nterminal_hits.tsv'
    command=[str(executable),'blastp','--db',str(database),'--query',str(query),'--out',str(hits),
             '--outfmt','6','qseqid','sseqid','pident','length','qlen','slen','qstart','qend','sstart','send','evalue','bitscore',
             '--max-target-seqs','5','--threads',str(threads)]
    if not keys:return {},{'state':'NO_CANDIDATES'}
    try:
        completed=subprocess.run(command,capture_output=True,text=True,check=False,timeout=300)
    except (OSError, subprocess.TimeoutExpired) as exc:
        status={'command':command,'state':'FAILED','error':str(exc)}
        (root/'status.json').write_text(json.dumps(status,indent=2))
        return {},status
    status={'command':command,'returncode':completed.returncode,'stderr':completed.stderr,'state':'COMPLETE' if completed.returncode==0 else 'FAILED'}
    (root/'status.json').write_text(json.dumps(status,indent=2))
    if completed.returncode or not hits.exists():return {},status
    anchors={}
    for line in hits.read_text().splitlines():
        f=line.split('\t')
        if len(f)!=12 or f[0] not in keys:continue
        identity,aln,qlen,slen,qs,qe,ss,se,evalue,bits=map(float,f[2:])
        if identity>=40 and evalue<=1e-5 and qs<=2 and ss<=2 and aln/qlen>=0.8 and aln/slen>=0.8:
            hit={'target':f[1],'identity':identity,'evalue':evalue,'bitscore':bits,'query_start':qs,'target_start':ss}
            if keys[f[0]] not in anchors or bits>anchors[keys[f[0]]]['bitscore']:anchors[keys[f[0]]]=hit
    return anchors,status


def run_refinement(genome, proteins, output, secondary=(), database=None, diamond=None, threads=1, config=None, rna_features=(), phrogs_start_db=None):
    config=config or RefinementConfig();root=Path(output);root.mkdir(parents=True,exist_ok=True)
    _,_,start_rows=refine_starts(genome,proteins,secondary,config=config)
    _,_,rescue_rows=rescue_short_orfs(genome,proteins,config=config,rna_features=rna_features)
    anchors={};status={'state':'NOT_CONFIGURED'}
    if database and diamond:
        anchors,status=diamond_anchors(genome,start_rows+rescue_rows,database,diamond,root/'homology',threads)
    if phrogs_start_db and diamond:
        extra, phrogs_status = diamond_anchors(genome, start_rows, phrogs_start_db, diamond, root/'phrogs_nterminal', threads)
        status = {'swissprot': status, 'phrogs': phrogs_status}
        for key, hit in extra.items():
            hit['source'] = 'PHROGs'
            anchors.setdefault(key, hit)
    refined,changes,candidates=refine_starts(genome,proteins,secondary,anchors,config)
    final,rescues,rescue_candidates=rescue_short_orfs(genome,refined,config,anchors,rna_features)
    with (root/'upstream_100nt.fna').open('w') as handle:
        for p in final:
            sequence=genome.upper() if p.strand=='+' else str(Seq(genome.upper()).reverse_complement())
            left,_=oriented_bounds(p,len(genome))
            handle.write(f'>{p.protein_id} upstream_nt={min(100,left)} includes_start_codon=true\n{sequence[max(0,left-100):left+3]}\n')
    for filename,data in [('start_candidates.json',candidates),('short_orf_candidates.json',rescue_candidates),('changes.json',changes+rescues)]:
        (root/filename).write_text(json.dumps(data,indent=2,sort_keys=True)+'\n')
    (root/'manifest.json').write_text(json.dumps({'config':asdict(config),'start_context':genome_start_context(genome, proteins, secondary),'homology':status,'changes':len(changes),'rescued':len(rescues),'reference_used':False},indent=2)+'\n')
    return final
