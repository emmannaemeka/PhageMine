"""Conservative reference-concordance scoring; semantic cases require human review."""
from __future__ import annotations
import csv
import math
import re
from collections import Counter

CATEGORIES = (
    'EXACT_PRODUCT_AGREEMENT', 'EQUIVALENT_FUNCTION', 'COMPATIBLE_BUT_BROADER',
    'COMPATIBLE_BUT_MORE_SPECIFIC', 'ABSTENTION', 'UNSUPPORTED_SPECIFICITY',
    'GENUINE_FUNCTIONAL_DISAGREEMENT', 'NON_EVALUABLE_REFERENCE',
)
CORRECT = set(CATEGORIES[:2])
COMPATIBLE = CORRECT | set(CATEGORIES[2:4])
PENDING = 'UNRESOLVED_REVIEW_REQUIRED'


def normalize(product):
    # Preserve putative/probable, specificity, numbers and subunit designations.
    return re.sub(r'[^a-z0-9]+', ' ', str(product or '').lower()).strip()


def informative(product):
    text = normalize(product)
    return bool(text) and not re.search(
        r'\b(hypothetical|uncharacterized|uncharacterised|unknown|unassigned)\b', text
    ) and text not in {'protein', 'conserved protein', 'phage protein', 'conserved phage protein'}


def load_synonyms(path):
    rules = {}
    with open(path) as handle:
        for row in csv.DictReader(handle, delimiter='\t'):
            a, b = normalize(row['term_a']), normalize(row['term_b'])
            if not informative(a) or not informative(b) or not row['rationale'].strip():
                raise ValueError('Synonyms require informative terms and explicit rationale')
            if a == b or (a, b) in rules or (b, a) in rules:
                raise ValueError('Duplicate or normalized-identical synonym rule')
            rules[a, b] = row['rule_id']
            rules[b, a] = row['rule_id']
    return rules


def classify(prediction, reference, synonyms):
    if not informative(reference):
        return 'NON_EVALUABLE_REFERENCE'
    if not informative(prediction):
        return 'ABSTENTION'
    if normalize(prediction) == normalize(reference):
        return 'EXACT_PRODUCT_AGREEMENT'
    if (normalize(prediction), normalize(reference)) in synonyms:
        return 'EQUIVALENT_FUNCTION'
    return PENDING


def locus(record):
    return record['accession'], int(record['start']), int(record['end']), record['strand']


def match_loci(predictions, references, mode='exact'):
    """Exact first. Relaxed thresholds match structural protocol, but ambiguous
    overlaps are withheld, never forced onto an unrelated CDS. Return indices.
    """
    if mode not in {'exact', 'relaxed'}:
        raise ValueError(mode)
    for records in (predictions, references):
        if len({locus(r) for r in records}) != len(records):
            raise ValueError('Duplicate locus')
    index = {locus(r): j for j, r in enumerate(references)}
    matches = {i: index[locus(p)] for i, p in enumerate(predictions) if locus(p) in index}
    if mode == 'exact':
        return matches, {}
    used = set(matches.values())
    candidates = {}
    for i, p in enumerate(predictions):
        if i in matches:
            continue
        hits = []
        for j, r in enumerate(references):
            if j in used or (p['accession'], p['strand']) != (r['accession'], r['strand']):
                continue
            overlap = max(0, min(p['end'], r['end']) - max(p['start'], r['start']) + 1)
            if overlap / (r['end']-r['start']+1) >= .5 and overlap / (p['end']-p['start']+1) >= .2:
                hits.append(j)
        if hits:
            candidates[i] = hits
    frequency = Counter(j for hits in candidates.values() for j in hits)
    ambiguous = {}
    for i, hits in candidates.items():
        if len(hits) == 1 and frequency[hits[0]] == 1:
            matches[i] = hits[0]
        else:
            ambiguous[i] = hits
    return matches, ambiguous


def wilson(k, n):
    if not n:
        return None, None
    z = 1.959963984540054
    p = k/n
    center = (p + z*z/(2*n))/(1+z*z/n)
    half = z*math.sqrt(p*(1-p)/n+z*z/(4*n*n))/(1+z*z/n)
    return max(0, center-half), min(1, center+half)


def summarize(categories):
    counts = Counter(categories)
    n = len(categories)-counts['NON_EVALUABLE_REFERENCE']
    named = n-counts['ABSTENTION']
    correct = sum(counts[c] for c in CORRECT)
    compatible = sum(counts[c] for c in COMPATIBLE)
    pending = counts[PENDING]
    result = {'evaluable_named_reference_loci': n, 'named_assertions': named,
              'unresolved': pending, **{c: counts[c] for c in CATEGORIES}}
    for category in CATEGORIES:
        denominator = len(categories) if category == 'NON_EVALUABLE_REFERENCE' else n
        result[category+'_percent'] = 100*counts[category]/denominator if denominator else None
    ratios = {'functional_coverage': (named, n), 'strict_functional_precision': (correct, named),
              'expanded_compatible_precision': (compatible, named), 'functional_recall': (correct, n),
              'abstention_rate': (counts['ABSTENTION'], n),
              'unsupported_specificity_rate': (counts['UNSUPPORTED_SPECIFICITY'], named),
              'genuine_disagreement_rate': (counts['GENUINE_FUNCTIONAL_DISAGREEMENT'], named)}
    for metric, (k, d) in ratios.items():
        unresolved = pending and metric not in {'functional_coverage', 'abstention_rate'}
        value = k/d if d and not unresolved else None
        low, high = wilson(k, d) if not unresolved else (None, None)
        result.update({metric: value, metric+'_numerator': k, metric+'_denominator': d,
                       metric+'_percent': 100*value if value is not None else None,
                       metric+'_ci95_low': low, metric+'_ci95_high': high})
        if unresolved and d:
            result[metric+'_identification_low'] = k/d
            result[metric+'_identification_high'] = (k+pending)/d
    p, r = result['strict_functional_precision'], result['functional_recall']
    result['strict_functional_f1'] = None if p is None or r is None else 2*p*r/(p+r) if p+r else 0.0
    return result


def paired_test(left, right):
    if len(left) != len(right):
        raise ValueError('Paired outcomes must have equal lengths')
    if any(c == PENDING for c in left+right):
        return {'status': 'PENDING_MANUAL_ADJUDICATION', 'p_value': None}
    pairs = [(a in CORRECT, b in CORRECT) for a, b in zip(left, right)
             if 'NON_EVALUABLE_REFERENCE' not in (a, b)]
    b = sum(a and not c for a, c in pairs)
    c = sum(c and not a for a, c in pairs)
    n = b+c
    p = min(1.0, 2*sum(math.comb(n, k) for k in range(min(b,c)+1))/2**n) if n else 1.0
    return {'status': 'DESCRIPTIVE_EXACT_MCNEMAR_LOCUS_INDEPENDENCE_NOT_ESTABLISHED',
            'paired_loci': len(pairs), 'phagemine_correct_only': b, 'pharokka_correct_only': c,
            'p_value': p}
