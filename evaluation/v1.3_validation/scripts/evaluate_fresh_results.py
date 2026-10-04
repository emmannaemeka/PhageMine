#!/usr/bin/env python3
"""Combine completed fresh genome runs and audit current fusion decisions.

Reclassification reuses newly generated raw evidence, not archived predictions.
It changes no gene coordinates and performs no additional database searches.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / 'src'))
from phagemine.artifact_ops import checksum
from phagemine.benchmark import import_phagemine
from phagemine.curated_benchmark import evaluate
from phagemine.fusion import FUSION_RULES_VERSION
from phagemine.operational_validation import validate_run
from phagemine.progress import ProgressReporter
from phagemine.resume import reclassify
from run_archived import read, write


def validate_sources(sources):
    manifests = [json.loads((source / 'manifest.json').read_text()) for source in sources]
    if not manifests or any(m.get('fresh_annotation_run') is not True or
                            m.get('analysis_type') != 'FRESH_REFERENCE_CONCORDANCE' for m in manifests):
        raise ValueError('Only completed fresh benchmark results are accepted')
    if any(any(m['tool_status'].get(tool) != 'FRESH_RUN_COMPLETED'
               for tool in ('PhageMine', 'Pharokka', 'Prokka')) for m in manifests):
        raise ValueError('Incomplete tool comparison')
    if len({m['source_commit'] for m in manifests}) != 1:
        raise ValueError('Fresh runs use different source commits')
    if len({checksum(source / 'database_checksums.tsv') for source in sources}) != 1:
        raise ValueError('Fresh runs use different databases')
    genomes = [g['accession'] for m in manifests for g in m['genomes']]
    if len(set(genomes)) != len(genomes): raise ValueError('Duplicate genome results')
    for source, manifest in zip(sources, manifests):
        manifest['omitted_rendering_caches'] = []
        for item in read(source / 'SHA256SUMS.tsv'):
            relative = Path(item['path'])
            if relative.is_absolute() or '..' in relative.parts:
                raise ValueError('Unsafe fresh artifact path')
            # upload-artifact excludes hidden files by default. Only the
            # disposable font cache may be absent; scientific data must match.
            if (not (source / relative).exists() and len(relative.parts) >= 3 and
                    relative.parts[-3:-1] == ('figures', '.matplotlib') and
                    relative.name.startswith('fontlist-') and relative.suffix == '.json'):
                manifest['omitted_rendering_caches'].append(str(relative))
                continue
            if checksum(source / relative) != item['sha256']:
                raise ValueError('Fresh artifact checksum mismatch')
        for genome in manifest['genomes']:
            if checksum(source / 'inputs' / (genome['accession'] + '.fasta')) != genome['fasta_sha256']:
                raise ValueError('Fresh input genome checksum mismatch')
    return manifests


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--results', nargs='+', required=True, help='Extracted complete per-genome artifact directories')
    parser.add_argument('--output', required=True)
    args = parser.parse_args(argv)
    sources = [Path(p).resolve() for p in args.results]
    manifests = validate_sources(sources)
    output = Path(args.output).resolve()
    if output.exists(): parser.error('Output must be a new directory')
    output.mkdir(parents=True)
    inputs = output / 'inputs'; inputs.mkdir()
    references = []
    predictions = {tool: [] for tool in ('PhageMine', 'Pharokka', 'Prokka', 'PhageMine-current-rules')}
    changes, provenance = [], []
    for source, manifest in zip(sources, manifests):
        references.extend(read(source / 'inputs/references.tsv'))
        for tool in ('PhageMine', 'Pharokka', 'Prokka'):
            predictions[tool].extend(read(source / 'inputs' / (tool + '.tsv')))
        for genome in manifest['genomes']:
            accession = genome['accession']
            before = source / 'runs/PhageMine' / accession
            after = output / 'reclassified' / accession
            reclassify(before, after, ProgressReporter(quiet=True))
            validation = validate_run(after)
            if validation['status'] != 'ARTIFACTS_CONSISTENT': raise ValueError('Reclassified artifacts are inconsistent')
            if checksum(before / 'evidence.json') != checksum(after / 'evidence.json'):
                raise ValueError('Reclassification changed raw evidence')
            old = {r['protein_id']: r for r in read(before / 'functional_classification.tsv')}
            new = read(after / 'functional_classification.tsv')
            for row in new:
                previous = old[row['protein_id']]
                if any(row[key] != previous[key] for key in ('start', 'end', 'strand')):
                    raise ValueError('Reclassification changed gene coordinates')
                if any(row[key] != previous[key] for key in ('display_product', 'gene', 'ec_number')):
                    changes.append({'accession': accession, 'protein_id': row['protein_id'],
                                    'start': row['start'], 'end': row['end'], 'strand': row['strand'],
                                    'before_product': previous['display_product'], 'after_product': row['display_product'],
                                    'before_gene': previous['gene'], 'after_gene': row['gene'],
                                    'before_ec': previous['ec_number'], 'after_ec': row['ec_number'],
                                    'review_status': 'DECISION_CHANGE_NOT_INDEPENDENTLY_ADJUDICATED'})
            predictions['PhageMine-current-rules'].extend(
                {'accession': accession, 'start': r['start'], 'end': r['end'], 'strand': r['strand'],
                 'product': r.get('product') or 'hypothetical protein', 'gene': r.get('gene') or ''}
                for r in import_phagemine(after))
            provenance.append({'accession': accession, 'fresh_manifest_sha256': checksum(source / 'manifest.json'),
                               'fresh_input_sha256': genome['fasta_sha256'],
                               'raw_evidence_sha256': checksum(before / 'evidence.json'),
                               'corrected_evidence_sha256': checksum(after / 'evidence.json')})
    reference_path = inputs / 'references.tsv'
    write(reference_path, references, ['accession', 'start', 'end', 'strand', 'product', 'gene', 'curated_by', 'reference_evidence'])
    summary = []
    for tool, rows in predictions.items():
        path = inputs / (tool + '.tsv'); write(path, rows, ['accession', 'start', 'end', 'strand', 'product', 'gene'])
        result = evaluate(path, reference_path, output / 'evaluation' / tool)
        summary.append({'tool': tool, **result['summary']})
    write(output / 'summary.tsv', summary, list(summary[0]))
    write(output / 'decision_changes.tsv', changes,
          ['accession', 'protein_id', 'start', 'end', 'strand', 'before_product', 'after_product', 'before_gene', 'after_gene', 'before_ec', 'after_ec', 'review_status'])
    report = {'analysis_type': 'FRESH_EVIDENCE_DECISION_AUDIT', 'raw_searches_generated_fresh': True,
              'source_commit': manifests[0]['source_commit'], 'current_fusion_rules_version': FUSION_RULES_VERSION,
              'current_fusion_sha256': checksum(ROOT / 'src/phagemine/fusion.py'),
              'runner_sha256': checksum(Path(__file__)), 'genomes': provenance,
              'database_checksums_table_sha256': checksum(sources[0] / 'database_checksums.tsv'),
              'omitted_non_scientific_rendering_caches': [p for m in manifests for p in m['omitted_rendering_caches']],
              'reclassification_method': 'Current rules applied to newly generated raw evidence; no repeat searches',
              'changed_assignments': len(changes), 'held_out': False, 'independence_established': False,
              'accuracy_superiority_established': False, 'scope': manifests[0]['scope'],
              'limitations': manifests[0]['comparison_limits'] + ['Decision changes require independent review; partial-profile filtering can increase abstention']}
    (output / 'manifest.json').write_text(json.dumps(report, indent=2) + '\n')


if __name__ == '__main__': main()
