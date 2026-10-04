#!/usr/bin/env python3
"""Evaluate immutable archived predictions; never claim a fresh 1.3 annotation run."""
from __future__ import annotations
import argparse
import csv
import json
from collections import Counter
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / 'src'))
from phagemine.artifact_ops import checksum
from phagemine.benchmark import import_gff
from phagemine.curated_benchmark import evaluate
from phagemine.functional_benchmark import informative, locus, match_loci


def read(path):
    with path.open(newline='') as handle: return list(csv.DictReader(handle, delimiter='\t'))


def write(path, rows, columns):
    with path.open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter='\t', extrasaction='ignore', lineterminator='\n')
        writer.writeheader(); writer.writerows(rows)


def error_cases(tool, predicted, references):
    matches, _ = match_loci(predicted, references)
    by_reference = {j: i for i, j in matches.items()}
    rows = []
    for j, reference in enumerate(references):
        if j in by_reference:
            prediction = predicted[by_reference[j]]
            if informative(reference['product']) and not informative(prediction['product']):
                kind = 'EXACT_MODEL_FUNCTION_ABSTENTION'
            elif informative(reference['product']) and prediction['product'].lower() != reference['product'].lower():
                kind = 'LABEL_DIFFERENCE_REVIEW_REQUIRED'
            else: continue
            count = 1
        else:
            overlap = [p for p in predicted if p['accession'] == reference['accession'] and
                       min(p['end'], reference['end']) >= max(p['start'], reference['start'])]
            same = [p for p in overlap if p['strand'] == reference['strand']]
            count = len(same)
            if len(same) == 1:
                prediction = same[0]
                left, right = prediction['start'] == reference['start'], prediction['end'] == reference['end']
                if left or right:
                    start_equal = left if reference['strand'] == '+' else right
                    kind = 'STOP_BOUNDARY_DIFFERENCE' if start_equal else 'START_BOUNDARY_DIFFERENCE'
                else: kind = 'BOTH_BOUNDARIES_DIFFER'
            elif same:
                prediction = None; kind = 'MULTIPLE_OVERLAPPING_MODELS_REVIEW_REQUIRED'
            else:
                prediction = None; kind = 'OPPOSITE_STRAND_OVERLAP' if overlap else 'NO_OVERLAPPING_MODEL'
        rows.append({'tool': tool, 'accession': reference['accession'], 'reference_start': reference['start'],
                     'reference_end': reference['end'], 'strand': reference['strand'], 'kind': kind,
                     'candidate_count': count, 'prediction_start': prediction['start'] if prediction else '',
                     'prediction_end': prediction['end'] if prediction else '',
                     'reference_product': reference['product'], 'prediction_product': prediction['product'] if prediction else '',
                     'interpretation': 'Review candidate; the reference and prediction have not been independently adjudicated'})
    return rows


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True)
    args = parser.parse_args(argv)
    output = Path(args.output)
    if output.exists(): parser.error('Output must be a new directory')
    source = ROOT / 'evaluation/v1.2_functional'
    # Protect historical inputs before using them in a separate new report.
    for item in read(source / 'source_provenance.tsv'):
        if checksum(ROOT / item['archived_path']) != item['sha256']:
            raise ValueError('Archived source checksum mismatch: ' + item['archived_path'])
    reference_raw = read(source / 'functional_truth_by_locus.tsv')
    references = [{key: row[key] for key in ('accession', 'start', 'end', 'strand', 'product')}
                  | {'curated_by': 'Archived RefSeq annotation; independent review not established',
                     'reference_evidence': 'evaluation/v1.2_functional/reference/' + row['accession'] + '.gb#' + row['locus_id']}
                  for row in reference_raw if row['compound_location'] != 'True']
    for row in references: row['start'], row['end'] = int(row['start']), int(row['end'])
    predictions = {'PhageMine': read(source / 'phagemine_predictions.tsv'), 'Pharokka': read(source / 'pharokka_predictions.tsv'), 'Prokka': []}
    for accession in sorted({r['accession'] for r in references}):
        for row in import_gff(source / 'predictions/Prokka' / accession / (accession + '.gff'), 'Prokka', accession):
            predictions['Prokka'].append({'accession': accession, 'start': row['start'], 'end': row['end'], 'strand': row['strand'], 'product': row.get('product') or ''})
    output.mkdir(parents=True)
    inputs = output / 'normalized_inputs'; inputs.mkdir()
    columns = ['accession', 'start', 'end', 'strand', 'product']
    write(inputs / 'references.tsv', references, columns + ['curated_by', 'reference_evidence'])
    summary, per_genome, errors = [], [], []
    sources = []
    for tool, records in predictions.items():
        for row in records: row['start'], row['end'] = int(row['start']), int(row['end'])
        path = inputs / (tool + '.tsv'); write(path, records, columns)
        result = evaluate(path, inputs / 'references.tsv', output / tool, synonyms=source / 'synonym_rules.tsv')
        summary.append({'tool': tool, **result['summary']})
        errors.extend(error_cases(tool, records, references))
        for accession in sorted({r['accession'] for r in references}):
            pp = [p for p in records if p['accession'] == accession]; rr = [r for r in references if r['accession'] == accession]
            matches, _ = match_loci(pp, rr)
            per_genome.append({'tool': tool, 'accession': accession, 'predicted_cds': len(pp), 'simple_reference_cds': len(rr),
                               'exact_matches': len(matches), 'missing_exact_models': len(rr)-len(matches),
                               'unmatched_prediction_models': len(pp)-len(matches),
                               'exact_model_precision': len(matches)/len(pp) if pp else None,
                               'exact_model_recall': len(matches)/len(rr)})
    write(output / 'summary.tsv', summary, list(summary[0]))
    write(output / 'per_genome.tsv', per_genome, list(per_genome[0]))
    write(output / 'error_candidates.tsv', errors, list(errors[0]))
    manifest = {'schema_version': 1, 'analysis_type': 'ARCHIVED_REFERENCE_CONCORDANCE', 'fresh_annotation_run': False,
                'prediction_release': 'Archived v1.2 outputs; not a v1.3 performance measurement',
                'reference_genomes': len({r['accession'] for r in references}), 'simple_reference_cds': len(references),
                'excluded_compound_reference_cds': len(reference_raw)-len(references),
                'held_out': False, 'independence_established': False, 'accuracy_superiority_established': False,
                'tool_status': {**{tool: 'ARCHIVED_OUTPUT_EVALUATED' for tool in predictions}, 'Phold': 'NOT_RUN_NO_OUTPUTS_AVAILABLE'},
                'error_kind_counts': {tool: dict(Counter(r['kind'] for r in errors if r['tool'] == tool)) for tool in predictions},
                'reference_input_sha256': checksum(source / 'functional_truth_by_locus.tsv'),
                'runner_sha256': checksum(Path(__file__)), 'scorer_sha256': checksum(ROOT / 'src/phagemine/curated_benchmark.py'),
                'missing_evidence': ['Independent completed reference review', 'Untouched held-out panel with sequence leakage checks', 'Phold predictions on the same genomes', 'Fresh v1.3 predictions with frozen database versions'],
                'interpretation': 'Coordinates measure agreement with archived references. Named products and label differences do not establish correctness. Unresolved product precision remains unset.'}
    (output / 'manifest.json').write_text(json.dumps(manifest, indent=2, sort_keys=True) + '\n')
    checks = [{'path':str(p.relative_to(output)), 'sha256':checksum(p)} for p in sorted(output.rglob('*')) if p.is_file()]
    (output / 'SHA256SUMS').write_text(''.join(row['sha256'] + '  ' + row['path'] + '\n' for row in checks))
    print(json.dumps({'output': str(output), 'tools': list(predictions), 'simple_reference_cds': len(references), 'status':'REFERENCE_CONCORDANCE_ONLY'}, sort_keys=True))
    return 0


if __name__ == '__main__': raise SystemExit(main())
