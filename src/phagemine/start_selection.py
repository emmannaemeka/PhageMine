"""Conservative, reference-independent selection among caller-supplied starts."""
from dataclasses import replace
from Bio.Seq import Seq
from Bio.Data.CodonTable import TranslationError

STARTS_11 = frozenset({'ATG', 'GTG', 'TTG', 'CTG', 'ATT', 'ATC', 'ATA'})
STOPS = frozenset({'TAA', 'TAG', 'TGA'})


def sd_motifs(genome, start, strand):
    """Literal SD-like motifs with a 3–15 nt spacer, not calibrated RBS scores."""
    upstream = (genome[max(0, start - 31):start - 1] if strand == '+'
                else str(Seq(genome[start:start + 30]).reverse_complement()))
    return [f'{motif}:{len(upstream) - i - len(motif)}'
            for motif in ('AGGAGG', 'GGAGG', 'AGGAG', 'GGAG', 'GAGG', 'AGGA')
            for i in range(len(upstream) - len(motif) + 1)
            if upstream[i:i + len(motif)] == motif
            and 3 <= len(upstream) - i - len(motif) <= 15]


def select_supported_starts(genome, proteins, secondary_models):
    """Keep gene identity/stop/product; change only uniquely supported starts.

    Eligible secondary models share the biological stop, strand and frame and
    encode a complete table-11 ORF. Prefer the unique Prodigal-GV start, retaining
    PHANOTATE when no eligible secondary model exists. Agreed starts are unchanged.
    Never cross stops or N gaps. Multiple alternatives remain unresolved.
    """
    genome = genome.upper()
    output, changes = [], []
    for p in proteins:
        output.append(p)
        if p.strand not in {'+', '-'} or not 1 <= p.start < p.end <= len(genome):
            continue
        if p.partial_5prime or p.partial_3prime or p.codon_start != 1:
            continue
        old_start, stop = (p.start, p.end) if p.strand == '+' else (p.end, p.start)
        original = genome[p.start - 1:p.end]
        if p.strand == '-': original = str(Seq(original).reverse_complement())
        if original != p.cds.upper() or len(original) % 3:
            continue
        # Existing stops/recoding candidates are not repaired by start selection.
        if original[-3:] not in STOPS or any(original[i:i+3] in STOPS for i in range(0, len(original)-3, 3)):
            continue
        eligible = {}
        for candidate in secondary_models:
            if candidate.caller.lower().replace("-", "_") not in {"prodigal_gv", "pyrodigal_gv"}:
                continue
            if candidate.strand != p.strand or not 1 <= candidate.start < candidate.end <= len(genome):
                continue
            new_start, new_stop = ((candidate.start, candidate.end) if p.strand == '+'
                                   else (candidate.end, candidate.start))
            distance = (new_start - old_start) * (1 if p.strand == '+' else -1)
            if new_stop != stop or distance % 3:
                continue
            sequence = genome[candidate.start - 1:candidate.end]
            if p.strand == '-': sequence = str(Seq(sequence).reverse_complement())
            if sequence[:3] not in STARTS_11 or 'N' in sequence:
                continue
            try:
                translation = str(Seq(sequence).translate(table=11, cds=True))
            except TranslationError:
                continue
            old_sd = sd_motifs(genome, old_start, p.strand)
            new_sd = sd_motifs(genome, new_start, p.strand)
            rule = 'SAME_STRAND_SAME_STOP_PRODIGAL_GV_START'
            eligible[new_start] = (candidate, sequence, translation, rule, distance, old_sd, new_sd)
        if len(eligible) != 1:
            params = dict(p.gene_call_parameters)
            params["start_selection_review"] = ("ambiguous secondary starts - review" if eligible
                                                else "no matching valid Prodigal-GV start - review")
            output[-1] = replace(p, gene_call_parameters=params)
            continue
        candidate, sequence, translation, rule, distance, old_sd, new_sd = next(iter(eligible.values()))
        if distance == 0:
            continue
        event = {'protein_id': p.protein_id, 'old_start': p.start, 'old_end': p.end,
                 'new_start': candidate.start, 'new_end': candidate.end, 'strand': p.strand,
                 'rule': rule, 'distance_nt': distance, 'caller': candidate.caller,
                 'candidate_id': candidate.identifier, 'old_sd_motifs': old_sd, 'new_sd_motifs': new_sd}
        params = dict(p.gene_call_parameters)
        params['start_selection'] = event
        output[-1] = replace(p, start=candidate.start, end=candidate.end, cds=sequence,
                             sequence=translation, start_codon=sequence[:3],
                             gene_call_parameters=params)
        changes.append(event)
    return output, changes
