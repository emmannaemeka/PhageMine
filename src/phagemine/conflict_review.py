"""Expose existing evidence conflicts and external label differences for review."""
from __future__ import annotations
import csv
import json
from pathlib import Path
from .artifact_ops import read_json, record_sidecar
from .functional_benchmark import normalize, informative


def write_conflicts(root):
    root = Path(root)
    rows = []
    classes = read_json(root / 'functional_classification.json', list)
    with (root / 'annotation.tsv').open(newline='') as handle:
        current = {r['protein_id']: r for r in csv.DictReader(handle, delimiter='\t')}
    evidence = {r['protein_id']: r['evidence'] for r in read_json(root / 'evidence.json', list)}
    for record in classes:
        if not record.get('conflicting_evidence_ids') and record.get('functional_state') != 'CONFLICTING_EVIDENCE': continue
        pid = record['protein_id']; annotation = current.get(pid, {})
        rows.append({'protein_id': pid, 'kind': 'ENGINE_EVIDENCE_CONFLICT', 'current_product': annotation.get('product'),
                     'proposed_product': None, 'sources': record.get('conflicting_sources') or [],
                     'evidence_ids': record.get('conflicting_evidence_ids') or [],
                     'reason': '; '.join(record.get('conflict_descriptions') or []) or 'The existing fusion engine flagged conflicting evidence',
                     'review_status': annotation.get('curation_state') or 'PENDING_REVIEW',
                     'evidence': evidence.get(pid, [])})
    external = read_json(root / 'external_evidence.json') if (root / 'external_evidence.json').is_file() else {}
    for record in external.get('records', []):
        pid = record['protein_id']; annotation = current.get(pid, {})
        product = annotation.get('product', record.get('current_product'))
        proposal = record['proposed_product']
        if not informative(product) or normalize(product) == normalize(proposal): continue
        rows.append({'protein_id': pid, 'kind': 'EXTERNAL_LABEL_DIFFERENCE', 'current_product': product,
                     'proposed_product': proposal, 'sources': [record['source']], 'evidence_ids': [record['external_id']],
                     'reason': 'Different product labels require review; wording alone does not establish biological incompatibility',
                     'review_status': 'PENDING_REVIEW', 'evidence': [record]})
    payload = {'schema_version': 1, 'record_count': len(rows), 'records': rows,
               'automatic_product_changes': False, 'confidence_changed': False}
    (root / 'annotation_conflicts.json').write_text(json.dumps(payload, indent=2, sort_keys=True) + '\n')
    columns = ['protein_id', 'kind', 'current_product', 'proposed_product', 'sources', 'evidence_ids', 'reason', 'review_status']
    with (root / 'annotation_conflicts.tsv').open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter='\t', extrasaction='ignore')
        writer.writeheader()
        for row in rows: writer.writerow({**row, 'sources': ';'.join(row['sources']), 'evidence_ids': ';'.join(row['evidence_ids'])})
    record_sidecar(root, 'conflict_review', {'records': len(rows), 'path': 'annotation_conflicts.tsv', 'automatic_resolution': False})
    return payload
