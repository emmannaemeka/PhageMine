#!/usr/bin/env python3
"""Generate new annotations from nucleotide-only inputs and compare names.

The familiar reference panel is diagnostic, not an independent held-out test.
Run in the PhageMine environment with a separate comparator conda environment.
"""
from __future__ import annotations
import argparse
import csv
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import subprocess
import sys
import time

from Bio import SeqIO

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / 'src'))
from phagemine.artifact_ops import checksum
from phagemine.benchmark import import_gff, import_phagemine
from phagemine.curated_benchmark import evaluate
from phagemine.resources import EvidenceResourceManager, ResourceType

from run_archived import error_cases, write

DATABASE_URL = 'https://zenodo.org/record/17110353/files/pharokka_v1.8.0_databases.tar.gz'
DATABASE_MD5 = 'a63c485241b900a11989bd1821bfbb09'


def prepare(source, output, accession=None):
    inputs = output / 'inputs'; inputs.mkdir()
    references, genomes, excluded = [], [], 0
    for path in sorted(source.glob('*.gb')):
        records = list(SeqIO.parse(path, 'genbank'))
        if len(records) != 1: raise ValueError('Expected one genome per reference file')
        record = records[0]
        if accession and record.id != accession: continue
        fasta = inputs / (record.id + '.fasta')
        # Neither features nor reference product names enter annotation inputs.
        fasta.write_text('>' + record.id + '\n' + str(record.seq) + '\n')
        genomes.append({'accession': record.id, 'fasta': str(fasta.resolve()),
                        'fasta_sha256': checksum(fasta), 'reference_sha256': checksum(path)})
        for feature in record.features:
            if feature.type != 'CDS': continue
            if len(feature.location.parts) != 1:
                excluded += 1; continue
            references.append({'accession': record.id, 'start': int(feature.location.start)+1,
                               'end': int(feature.location.end),
                               'strand': '+' if feature.location.strand == 1 else '-',
                               'product': feature.qualifiers.get('product', ['hypothetical protein'])[0],
                               'gene': feature.qualifiers.get('gene', [''])[0],
                               'curated_by': 'RefSeq annotation; independent review not established',
                               'reference_evidence': str(path)})
    if not genomes: raise ValueError('No reference genomes found')
    write(inputs / 'references.tsv', references,
          ['accession', 'start', 'end', 'strand', 'product', 'gene', 'curated_by', 'reference_evidence'])
    return genomes, references, excluded


def run(command, root, label, records, env=None):
    log = root / (label + '.log')
    print('START ' + label, flush=True)
    started = time.monotonic()
    with log.open('w') as handle:
        result = subprocess.run(command, stdout=handle, stderr=subprocess.STDOUT, env=env, check=False)
    record = {'command': command, 'returncode': result.returncode,
              'elapsed_seconds': round(time.monotonic()-started, 3), 'log': log.name,
              'log_sha256': checksum(log)}
    records.append(record)
    print(f'END {label}: exit={result.returncode}, seconds={record["elapsed_seconds"]}', flush=True)
    (root / 'commands.json').write_text(json.dumps(records, indent=2) + '\n')
    if result.returncode: raise RuntimeError(f'{label} failed; see {log}')


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True)
    parser.add_argument('--database', required=True)
    parser.add_argument('--comparator-environment', default='comparators')
    parser.add_argument('--threads', type=int, default=2)
    parser.add_argument('--accession', help='Run one panel genome for parallel evaluation')
    args = parser.parse_args(argv)
    output = Path(args.output).resolve()
    if output.exists(): parser.error('Output must be a new directory')
    output.mkdir(parents=True)
    started = datetime.now(timezone.utc).isoformat()
    genomes, references, excluded = prepare(ROOT / 'evaluation/v1.2_functional/reference', output, args.accession)
    database = Path(args.database).resolve()
    db_files = sorted(p for p in database.rglob('*') if p.is_file())
    if not db_files or not (database / 'VERSION_1_8_0').is_file():
        raise ValueError('Pharokka database 1.8.0 is required')
    write(output / 'database_checksums.tsv',
          [{'path': str(p.relative_to(database)), 'sha256': checksum(p)} for p in db_files], ['path', 'sha256'])
    registry = output / 'registry.json'
    EvidenceResourceManager(registry).register(
        'PHROGs', ResourceType.PHROGS, database / 'phrogs_profile_db',
        version='v4-pharokka-1.8.0', preparation_status='prepared', required_tools=['mmseqs'],
        provenance={'annotations_path': str(database / 'phrog_annot_v4.tsv'),
                    'hmm_profiles_path': str(database / 'all_phrogs.h3m'),
                    'source_url': DATABASE_URL, 'distribution_md5': DATABASE_MD5})
    env = {**os.environ, 'PHAGEMINE_REGISTRY_PATH': str(registry)}
    prefix = ['micromamba', 'run', '-n', args.comparator_environment]
    commands = []
    run([sys.executable, '-m', 'pip', 'freeze'], output, 'phagemine-packages', commands, env)
    run(['micromamba', 'list', '-n', args.comparator_environment, '--json'], output, 'comparator-packages', commands)
    predictions = {'PhageMine': [], 'Pharokka': [], 'Prokka': []}
    for genome in genomes:
        accession = genome['accession']; fasta = genome['fasta']
        for tool in predictions:
            destination = output / 'runs' / tool / accession
            destination.parent.mkdir(parents=True, exist_ok=True)
            if tool == 'PhageMine':
                command = ['phagemine', 'annotate', fasta, '--output', str(destination),
                           '--threads', str(args.threads), '--rna-features', 'none', '--no-progress']
            elif tool == 'Pharokka':
                command = prefix + ['pharokka', 'run', '-i', fasta, '-o', str(destination),
                                    '-d', str(database), '-t', str(args.threads), '-p', accession,
                                    '--skip_extra_annotations', '--skip_mash']
            else:
                command = prefix + ['prokka', '--outdir', str(destination), '--prefix', accession,
                                    '--kingdom', 'Viruses', '--gcode', '11', '--cpus', str(args.threads),
                                    '--norrna', '--notrna', fasta]
            run(command, output, tool + '-' + accession, commands, env)
            if tool == 'PhageMine':
                manifest = json.loads((destination / 'run_manifest.json').read_text())
                adapters = [a for a in manifest['evidence_adapters'] if 'PHROGS' in a['adapter'].upper()]
                if len(adapters) != 2 or any(a['status'] not in {'REAL', 'SUCCESS_WITH_HITS', 'SUCCESS_NO_HIT'} for a in adapters):
                    raise RuntimeError('Required PHROGs searches did not both succeed')
                rows = import_phagemine(destination)
            else:
                rows = import_gff(destination / (accession + '.gff'), tool, accession)
            if not rows: raise ValueError('Annotation returned no CDS: ' + tool + ' ' + accession)
            predictions[tool].extend({'accession': accession, 'start': row['start'], 'end': row['end'],
                                      'strand': row['strand'], 'product': row.get('product') or 'hypothetical protein',
                                      'gene': row.get('gene') or ''} for row in rows)
    summary, errors, names = [], [], []
    reference_by_locus = {(r['accession'], r['start'], r['end'], r['strand']): r for r in references}
    for tool, rows in predictions.items():
        path = output / 'inputs' / (tool + '.tsv')
        write(path, rows, ['accession', 'start', 'end', 'strand', 'product', 'gene'])
        result = evaluate(path, output / 'inputs/references.tsv', output / 'evaluation' / tool)
        summary.append({'tool': tool, **result['summary']})
        errors.extend(error_cases(tool, rows, references))
        for row in rows:
            reference = reference_by_locus.get((row['accession'], row['start'], row['end'], row['strand']))
            if reference and row['gene'] != reference['gene']:
                names.append({'tool': tool, **row, 'reference_gene': reference['gene'],
                              'reference_product': reference['product'], 'status': 'DIFFERENCE_REQUIRES_REVIEW'})
    write(output / 'summary.tsv', summary, list(summary[0]))
    write(output / 'error_candidates.tsv', errors, ['tool', 'accession', 'reference_start', 'reference_end', 'strand', 'kind', 'candidate_count', 'prediction_start', 'prediction_end', 'reference_product', 'prediction_product', 'interpretation'])
    write(output / 'gene_name_review.tsv', names, ['tool', 'accession', 'start', 'end', 'strand', 'product', 'gene', 'reference_gene', 'reference_product', 'status'])
    manifest = {'analysis_type': 'FRESH_REFERENCE_CONCORDANCE', 'fresh_annotation_run': True,
                'started_at': started, 'completed_at': datetime.now(timezone.utc).isoformat(),
                'source_commit': os.environ.get('GITHUB_SHA', 'LOCAL_NOT_RECORDED'), 'genomes': genomes,
                'excluded_compound_reference_cds': excluded, 'held_out': False,
                'independence_established': False, 'accuracy_superiority_established': False,
                'tool_status': {**{tool: 'FRESH_RUN_COMPLETED' for tool in predictions}, 'Phold': 'NOT_RUN'},
                'scope': 'PHROGs-only PhageMine vs standard Pharokka vs viral Prokka; CDS comparison',
                'comparison_limits': ['Search thresholds differ; commands and adapter manifests retain them',
                                      'Prokka uses its packaged databases; PHROGs is shared by PhageMine and Pharokka',
                                      'Reference gene symbols are locus-specific; differences are not proven errors',
                                      'This familiar panel is not a held-out accuracy benchmark']}
    (output / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
    write(output / 'SHA256SUMS.tsv', [{'path': str(p.relative_to(output)), 'sha256': checksum(p)}
                                   for p in sorted(output.rglob('*')) if p.is_file()], ['path', 'sha256'])


if __name__ == '__main__': main()
