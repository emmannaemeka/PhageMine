#!/usr/bin/env python3
"""Generate fresh annotations with every installed PhageMine resource.

Reference products remain evaluation targets, never annotation inputs.
Missing resources or failed searches invalidate the full-system benchmark.
"""
from __future__ import annotations
import argparse
from collections import Counter
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / 'src'))
from phagemine.artifact_ops import checksum
from phagemine.benchmark import import_phagemine
from phagemine.curated_benchmark import evaluate
from phagemine.operational_validation import validate_run
from phagemine.resources import EvidenceResourceManager
from run_fresh import prepare, run
from run_archived import write

REQUIRED_RESOURCES = {'PFAM', 'VOGDB', 'SWISSPROT', 'PHROGS', 'PMFDB', 'INPHARED_GENOMES'}
REQUIRED_ADAPTERS = {'PfamHMMAdapter', 'VOGHMMAdapter', 'SwissProtEvidenceAdapter',
                     'PHROGSMMseqsAdapter', 'PHROGSPyHMMERAdapter'}
SUCCESS = {'REAL', 'SUCCESS_WITH_HITS', 'SUCCESS_NO_HIT'}


def verify_full_run(destination):
    manifest = json.loads((destination / 'run_manifest.json').read_text())
    adapters = manifest.get('evidence_adapters', [])
    counts = Counter(a['adapter'] for a in adapters)
    if any(counts[name] != 1 for name in REQUIRED_ADAPTERS):
        raise ValueError('Full-system run must contain every required evidence adapter exactly once')
    if any(a['status'] not in SUCCESS for a in adapters if a['adapter'] in REQUIRED_ADAPTERS):
        raise ValueError('A required annotation database search did not succeed')
    comparative = manifest.get('comparative_analysis', {})
    if comparative.get('inphared', {}).get('status') != 'COMPLETE':
        raise ValueError('INPHARED genome comparison did not complete')
    pmfdb = comparative.get('pmfdb') or {}
    if pmfdb.get('reference_type') != 'PMFDB' or not pmfdb.get('command'):
        raise ValueError('PMFDB comparison did not execute MMseqs2')
    integrity = validate_run(destination)
    proteins = json.loads((destination / 'evidence.json').read_text())
    classifications = json.loads((destination / 'functional_classification.json').read_text())
    return {'adapters': adapters, 'comparative_analysis': comparative, 'artifact_integrity': integrity,
            'accepted_evidence_by_source': dict(Counter(e['source'] for p in proteins
                for e in p.get('evidence', []) if e.get('supports'))),
            'functional_states': dict(Counter(c['functional_state'] for c in classifications)),
            'selected_products_by_supporting_source': dict(Counter(source for c in classifications
                for source in c.get('product_supporting_sources', []))),
            'curated_anchor_selections': sum(bool(c.get('selected_by_curated_phage_anchor')) for c in classifications),
            'fusion_rules_versions': sorted({c['fusion_rules_version'] for c in classifications})}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True, type=Path)
    parser.add_argument('--registry', required=True, type=Path)
    parser.add_argument('--threads', type=int, default=2)
    args = parser.parse_args(argv)
    output = args.output.resolve()
    if output.exists(): parser.error('Output must be a new directory')
    os.environ['PHAGEMINE_REGISTRY_PATH'] = str(args.registry.resolve())
    manager = EvidenceResourceManager(args.registry.resolve())
    resources = manager.validate_all()
    ready = {r['resource_type'] for r in resources if r['status'] == 'READY'}
    if not REQUIRED_RESOURCES <= ready:
        raise ValueError('Missing READY resources: ' + ', '.join(sorted(REQUIRED_RESOURCES - ready)))
    output.mkdir(parents=True)
    started = datetime.now(timezone.utc).isoformat()
    genomes, references, excluded = [], [], 0
    for label, source in [('original', ROOT / 'evaluation/v1.2_functional/reference'),
                          ('external', ROOT / 'evaluation/v1.3_validation/external_panel')]:
        panel = output / label; panel.mkdir()
        g, r, e = prepare(source, panel)
        genomes.extend({**row, 'panel': label} for row in g)
        references.extend(r); excluded += e
    write(output / 'references.tsv', references, list(references[0]))
    commands, predictions, audits = [], [], []
    run([sys.executable, '-m', 'pip', 'freeze'], output, 'packages', commands)
    # Keep versions, provider checksums and preparation details, not bulky databases.
    (output / 'resources.json').write_text(json.dumps(resources, indent=2) + '\n')
    installs = output / 'database-manifests'; installs.mkdir()
    for resource in resources:
        if resource['resource_type'] not in REQUIRED_RESOURCES: continue
        source = Path(resource['provenance']['install_manifest'])
        (installs / (resource['resource_type'] + '.json')).write_bytes(source.read_bytes())
    for genome in genomes:
        accession = genome['accession']
        destination = output / 'runs' / accession
        destination.parent.mkdir(exist_ok=True)
        run(['phagemine', 'run', genome['fasta'], '--output', str(destination),
             '--threads', str(args.threads), '--rna-features', 'none', '--no-progress'],
            output, accession, commands)
        audit = verify_full_run(destination)
        audits.append({'accession': accession, **audit})
        rows = import_phagemine(destination)
        predictions.extend({'accession': accession, 'start': r['start'], 'end': r['end'],
                            'strand': r['strand'], 'product': r.get('product') or 'hypothetical protein',
                            'gene': r.get('gene') or ''} for r in rows)
        (output / 'source_audits.json').write_text(json.dumps(audits, indent=2) + '\n')
        write(output / 'PhageMine-full.tsv', predictions,
              ['accession', 'start', 'end', 'strand', 'product', 'gene'])
    for panel in ('original', 'external', 'combined'):
        accessions = {g['accession'] for g in genomes if panel == 'combined' or g['panel'] == panel}
        directory = output / 'evaluation' / panel; directory.mkdir(parents=True)
        write(directory / 'predictions.tsv', [r for r in predictions if r['accession'] in accessions], list(predictions[0]))
        write(directory / 'references.tsv', [r for r in references if r['accession'] in accessions], list(references[0]))
        evaluate(directory / 'predictions.tsv', directory / 'references.tsv', directory / 'scores',
                 synonyms=ROOT / 'evaluation/v1.3_validation/naming_equivalences.tsv')
    manifest = {'scope': 'All six PhageMine resources; fresh full-system runs; no manual product overrides',
                'started_at': started, 'completed_at': datetime.now(timezone.utc).isoformat(),
                'source_commit': os.environ.get('GITHUB_SHA', 'LOCAL_NOT_RECORDED'),
                'fresh_annotation_run': True, 'required_resources': sorted(REQUIRED_RESOURCES),
                'genomes': genomes, 'excluded_compound_reference_cds': excluded,
                'accuracy_superiority_established': False, 'independence_established': False,
                'limitations': ['Database overlap with reference genomes is not excluded',
                    'Reference-name concordance is not independent functional accuracy',
                    'PMFDB and INPHARED provide comparative context; their labels are not automatically transferred as functions']}
    (output / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
    write(output / 'SHA256SUMS.tsv', [{'path': str(p.relative_to(output)), 'sha256': checksum(p)}
        for p in sorted(output.rglob('*')) if p.is_file() and '.matplotlib' not in p.parts], ['path', 'sha256'])


if __name__ == '__main__': main()
