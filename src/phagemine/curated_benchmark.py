"""Auditable benchmark reporting against explicitly curated reference annotations."""
from __future__ import annotations
import csv
import json
from pathlib import Path
from .artifact_ops import checksum
from .functional_benchmark import CATEGORIES, PENDING, classify, informative, load_synonyms, locus, match_loci, summarize


def read_records(path, curated=False):
    with Path(path).open(newline='') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        required = {'accession', 'start', 'end', 'strand', 'product'}
        if curated: required |= {'curated_by', 'reference_evidence'}
        if not required.issubset(reader.fieldnames or []): raise ValueError('Missing TSV columns: ' + ', '.join(sorted(required)))
        records = list(reader)
    if not records: raise ValueError('Benchmark table is empty')
    for row in records:
        row['start'], row['end'] = int(row['start']), int(row['end'])
        if not row['accession'].strip() or row['start'] < 1 or row['end'] < row['start'] or row['strand'] not in {'+', '-'}:
            raise ValueError('Invalid benchmark locus')
        if curated and (not row['curated_by'].strip() or not row['reference_evidence'].strip()):
            raise ValueError('Each reference requires curator and evidence provenance')
    if len({locus(row) for row in records}) != len(records): raise ValueError('Duplicate locus')
    return records


def evaluate(predictions, references, output, *, synonyms=None, adjudications=None):
    output = Path(output)
    if output.exists(): raise ValueError('Benchmark output must be a new directory')
    if checksum(predictions) == checksum(references): raise ValueError('Prediction and reference inputs cannot be identical')
    predicted = read_records(predictions); reference = read_records(references, curated=True)
    rules = load_synonyms(synonyms) if synonyms else {}
    matches, _ = match_loci(predicted, reference)
    by_reference = {j: i for i, j in matches.items()}
    decisions = {}
    if adjudications:
        with Path(adjudications).open(newline='') as handle:
            for row in csv.DictReader(handle, delimiter='\t'):
                key = locus(row)
                if key in decisions: raise ValueError('Duplicate adjudication')
                if row.get('category') not in CATEGORIES or not row.get('reviewer', '').strip() or not row.get('rationale', '').strip():
                    raise ValueError('Adjudications require a supported category, reviewer and rationale')
                decisions[key] = row
    rows = []; used = set()
    for j, ref in enumerate(reference):
        pred = predicted[by_reference[j]] if j in by_reference else None
        category = classify(pred['product'] if pred else None, ref['product'], rules)
        row = {**ref, 'predicted_product': pred['product'] if pred else '',
               'model_match': 'EXACT' if pred else 'MISSING_EXACT_MODEL', 'category': category,
               'adjudication_reviewer': '', 'adjudication_rationale': ''}
        decision = decisions.get(locus(ref))
        if decision:
            if category != PENDING: raise ValueError('Only unresolved matched product differences can be adjudicated')
            if decision['category'] in {'ABSTENTION', 'NON_EVALUABLE_REFERENCE'}:
                raise ValueError('Adjudication cannot change the evaluability or named-assertion denominator')
            row.update(category=decision['category'], adjudication_reviewer=decision['reviewer'], adjudication_rationale=decision['rationale'])
            used.add(locus(ref))
        rows.append(row)
    if set(decisions) != used: raise ValueError('Adjudication includes a locus absent from unresolved comparisons')
    summary = summarize([row['category'] for row in rows])
    summary.update(exact_model_matches=len(matches), predicted_models=len(predicted), reference_models=len(reference),
                   unmatched_prediction_models=len(predicted)-len(matches),
                   exact_model_precision=len(matches)/len(predicted), exact_model_recall=len(matches)/len(reference))
    manifest = {'schema_version': 1, 'status': 'PENDING_ADJUDICATION' if summary['unresolved'] else 'SCORED',
                'summary': summary, 'source_checksums': {'predictions': checksum(predictions), 'references': checksum(references)},
                'matching': 'EXACT_ACCESSION_COORDINATES_STRAND', 'accuracy_superiority_established': False,
                'interpretation': 'Reference quality and independence require external assessment. Wilson intervals describe loci; loci within genomes are not necessarily independent. Missing exact models are counted as abstentions at evaluable reference loci.'}
    for name, path in (('synonyms', synonyms), ('adjudications', adjudications)):
        if path: manifest['source_checksums'][name] = checksum(path)
    output.mkdir(parents=True)
    (output / 'curated_benchmark.json').write_text(json.dumps(manifest, indent=2, sort_keys=True) + '\n')
    with (output / 'locus_results.tsv').open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter='\t'); writer.writeheader(); writer.writerows(rows)
    with (output / 'review_queue.tsv').open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=['accession', 'start', 'end', 'strand', 'category', 'reviewer', 'rationale'], delimiter='\t', extrasaction='ignore'); writer.writeheader()
        for row in rows:
            if row['category'] == PENDING: writer.writerow({**row, 'category': '', 'reviewer': '', 'rationale': ''})
    return manifest
