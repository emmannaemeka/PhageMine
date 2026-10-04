#!/usr/bin/env python3
"""Validate explicit evidence-review decisions; never infer biological truth."""
from __future__ import annotations
import argparse
from collections import Counter
import csv
import hashlib
import json
from pathlib import Path
import sys

from Bio import SeqIO

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / 'src'))
from phagemine.artifact_ops import checksum
from phagemine.functional_benchmark import locus

VERDICTS = {'CONSISTENT_WITH_RECORD', 'COMPATIBLE_BROAD_ROLE',
            'CONFLICTS_WITH_RECORD', 'SPECIFICITY_NOT_ESTABLISHED',
            'NOT_ESTABLISHED_BY_RECORD', 'ABSTENTION', 'UNRESOLVED_NO_EXACT_MATCH'}
DIFFERENCES = {'DIFFERENT_NAMED_LABELS', 'A_UNNAMED_B_NAMED', 'A_NAMED_B_UNNAMED'}


def read(path):
    with Path(path).open(newline='') as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def digest_sequence(sequence):
    return hashlib.sha256(sequence.encode()).hexdigest()


def validate_decision(case, decision, sequence, records):
    if any(decision.get(k) != case.get(k) for k in ('product_a', 'product_b')):
        raise ValueError('Reviewed prediction changed')
    if decision.get('protein_sequence_sha256') != digest_sequence(sequence):
        raise ValueError('Reviewed sequence changed')
    if any(decision.get(k) not in VERDICTS for k in ('verdict_a', 'verdict_b')):
        raise ValueError('Unsupported review verdict')
    if not decision.get('reviewer', '').strip() or not decision.get('rationale', '').strip():
        raise ValueError('Reviewer and rationale are required')
    accession = decision.get('reviewed_accession', '')
    if not accession:
        if any(decision[k] != 'UNRESOLVED_NO_EXACT_MATCH' for k in ('verdict_a', 'verdict_b')):
            raise ValueError('A supported verdict requires an exact reviewed sequence')
        return None
    record = records.get(accession)
    if record is None or record['sequence']['value'] != sequence:
        raise ValueError('Curated record is not an exact sequence match')
    record_hash = hashlib.sha256(json.dumps(record, sort_keys=True, separators=(',', ':')).encode()).hexdigest()
    if decision.get('reviewed_record_sha256') != record_hash:
        raise ValueError('Reviewed record changed')
    return record


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('paired', 'runs-root', 'snapshot', 'decisions', 'output'):
        parser.add_argument('--' + name, required=True, type=Path)
    args = parser.parse_args(argv)
    if args.output.exists():
        raise ValueError('Output must be a new directory')
    snapshot = json.loads(args.snapshot.read_text())
    records = {r['primaryAccession']: r for r in snapshot['results']}
    if len(records) != len(snapshot['results']) or any('reviewed' not in r['entryType'].lower() or
                                                     'unreviewed' in r['entryType'].lower()
                                                     for r in records.values()):
        raise ValueError('Snapshot must contain unique reviewed records')
    cases = [r for r in read(args.paired) if r['category'] in DIFFERENCES]
    decisions = read(args.decisions)
    by_locus = {locus(r): r for r in decisions}
    if len(by_locus) != len(decisions) or set(by_locus) != {locus(r) for r in cases}:
        raise ValueError('Review must cover exactly the paired differences, once each')
    sequences = {}
    for directory in sorted(args.runs_root.glob('NC_*/runs/PhageMine/*')):
        classes = {r['protein_id']: r for r in read(directory / 'functional_classification.tsv')}
        for protein in SeqIO.parse(directory / 'proteins.faa', 'fasta'):
            row = classes[protein.id]
            key = (directory.name, int(row['start']), int(row['end']), row['strand'])
            if key in sequences: raise ValueError('Duplicate predicted locus')
            sequences[key] = str(protein.seq)
    reviewed = []
    for case in cases:
        decision = by_locus[locus(case)]
        if locus(case) not in sequences: raise ValueError('Prediction sequence is missing')
        record = validate_decision(case, decision, sequences[locus(case)], records)
        reviewed.append({**case, **decision, 'review_status': 'EVIDENCE_REVIEW_COMPLETE',
                         'evidence_match': 'EXACT_SEQUENCE' if record else 'NO_EXACT_REVIEWED_SEQUENCE'})
    summary = {'analysis_type': 'CURATED_RECORD_REVIEW_OF_SELECTED_DISAGREEMENTS',
               'reviewed_cases': len(reviewed),
               'exact_reviewed_sequence_matches': sum(r['evidence_match'] == 'EXACT_SEQUENCE' for r in reviewed),
               'verdicts': {tool: dict(Counter(r[field] for r in reviewed))
                           for tool, field in (('PhageMine', 'verdict_a'), ('Pharokka', 'verdict_b'))},
               'functional_accuracy': None, 'accuracy_superiority_established': False,
               'held_out': False, 'independent_expert_adjudication': False,
               'sources': {name: checksum(path) for name, path in
                           (('paired', args.paired), ('snapshot', args.snapshot), ('decisions', args.decisions))},
               'limitations': ['Selected disagreements are not an unbiased accuracy sample.',
                               'Reviewed records can contain inferred or uncertain functions.',
                               'A missing function in a reviewed record does not prove a prediction false.',
                               'Database/reference overlap is not excluded.',
                               'Review is by Codex, not an independent domain expert.']}
    args.output.mkdir(parents=True)
    with (args.output / 'reviewed_cases.tsv').open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(reviewed[0]), delimiter='\t')
        writer.writeheader(); writer.writerows(reviewed)
    (args.output / 'review_summary.json').write_text(json.dumps(summary, indent=2, sort_keys=True) + '\n')
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
