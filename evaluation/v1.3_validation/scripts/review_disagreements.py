#!/usr/bin/env python3
"""Audit paired annotation differences without inferring biological correctness."""
from __future__ import annotations
import argparse
from collections import Counter
import csv
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / 'src'))
from phagemine.artifact_ops import checksum
from phagemine.curated_benchmark import read_records
from phagemine.functional_benchmark import locus, normalize
from phagemine.annotation_labels import informative


def audit(left, right, references, output):
    output = Path(output)
    if output.exists():
        raise ValueError('Output must be a new directory')
    tables = [read_records(path, curated=i == 2)
              for i, path in enumerate((left, right, references))]
    indexes = [{locus(row): row for row in rows} for rows in tables]
    rows = []
    for key in sorted(set(indexes[0]) | set(indexes[1]) | set(indexes[2])):
        a, b, ref = (index.get(key) for index in indexes)
        if a is None or b is None:
            category = 'MODEL_MISSING_BOTH' if a is b is None else 'MODEL_MISSING_A' if a is None else 'MODEL_MISSING_B'
        elif normalize(a['product']) == normalize(b['product']):
            category = 'SAME_NORMALIZED_LABEL'
        elif not informative(a['product']) and not informative(b['product']):
            category = 'BOTH_UNNAMED'
        elif not informative(a['product']):
            category = 'A_UNNAMED_B_NAMED'
        elif not informative(b['product']):
            category = 'A_NAMED_B_UNNAMED'
        else:
            category = 'DIFFERENT_NAMED_LABELS'
        rows.append(dict(zip(('accession', 'start', 'end', 'strand'), key),
                         product_a=a['product'] if a else '', product_b=b['product'] if b else '',
                         reference_product=ref['product'] if ref else '', reference_model_present=bool(ref),
                         category=category, review_status='NOT_ADJUDICATED',
                         independent_evidence='', reviewer='', conclusion='', rationale=''))
    report = {'analysis_type': 'PAIRED_DISAGREEMENT_TRIAGE', 'held_out': False,
              'accuracy_superiority_established': False, 'functional_accuracy': None,
              'counts': dict(Counter(row['category'] for row in rows)),
              'sources': {name: {'path': str(path), 'sha256': checksum(path)}
                          for name, path in zip(('A', 'B', 'reference'), (left, right, references))},
              'interpretation': 'Labels and reference agreement are not independent correctness evidence. Missing models are separate from unnamed functions. No automatic biological adjudication.'}
    output.mkdir(parents=True)
    with (output / 'paired_review.tsv').open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter='\t')
        writer.writeheader(); writer.writerows(rows)
    (output / 'manifest.json').write_text(json.dumps(report, indent=2, sort_keys=True) + '\n')
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('left', 'right', 'references', 'output'):
        parser.add_argument('--' + name, required=True, type=Path)
    args = parser.parse_args()
    print(json.dumps(audit(args.left, args.right, args.references, args.output), indent=2))
