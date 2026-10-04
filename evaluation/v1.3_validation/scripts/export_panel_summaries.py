#!/usr/bin/env python3
"""Report original and external diagnostic panels without changing predictions."""
from __future__ import annotations
import argparse
import csv
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / 'src'))
from phagemine.artifact_ops import checksum
from phagemine.curated_benchmark import evaluate, read_records


def write(path, rows):
    with path.open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter='\t')
        writer.writeheader(); writer.writerows(rows)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('audit', 'selection', 'equivalences', 'output'):
        parser.add_argument('--' + name, required=True, type=Path)
    args = parser.parse_args(argv)
    if args.output.exists(): raise ValueError('Output must be a new directory')
    manifest = json.loads((args.audit / 'manifest.json').read_text())
    if manifest.get('raw_searches_generated_fresh') is not True:
        raise ValueError('A verified fresh-evidence audit is required')
    references = read_records(args.audit / 'inputs/references.tsv', curated=True)
    accessions = {r['accession'] for r in references}
    selected = json.loads(args.selection.read_text())
    external = {r['accession'] for r in selected['genomes']}
    if not external or not external < accessions:
        raise ValueError('Both original and external reference panels are required')
    predictions = {tool: read_records(args.audit / 'inputs' / (source + '.tsv'))
                   for tool, source in (('PhageMine', 'PhageMine-current-rules'),
                                        ('Pharokka', 'Pharokka'), ('Prokka', 'Prokka'))}
    if any({r['accession'] for r in rows} - accessions for rows in predictions.values()):
        raise ValueError('Prediction contains a genome outside the reference panel')
    args.output.mkdir(parents=True)
    summary = []
    for panel, genomes in (('original', accessions - external), ('external', external), ('combined', accessions)):
        directory = args.output / panel; directory.mkdir()
        ref_path = directory / 'references.tsv'
        write(ref_path, [r for r in references if r['accession'] in genomes])
        for tool, rows in predictions.items():
            path = directory / (tool + '.tsv')
            write(path, [r for r in rows if r['accession'] in genomes])
            result = evaluate(path, ref_path, directory / ('evaluation-' + tool), synonyms=args.equivalences)
            s = result['summary']
            summary.append({'panel': panel, 'tool': tool, 'genomes': len(genomes),
                            'predicted_cds': s['predicted_models'], 'reference_cds': s['reference_models'],
                            'exact_reference_models': s['exact_model_matches'],
                            'named_reference_loci': s['evaluable_named_reference_loci'],
                            'named_assertions': s['named_assertions'],
                            'literal_name_agreements': s['EXACT_PRODUCT_AGREEMENT'],
                            'equivalent_name_agreements': s['EQUIVALENT_FUNCTION'],
                            'exact_or_equivalent_agreements': s['EXACT_PRODUCT_AGREEMENT'] + s['EQUIVALENT_FUNCTION'],
                            'unresolved_named_differences': s['unresolved'],
                            'missing_models_at_named_loci': s['missing_models_at_named_reference_loci'],
                            'abstentions_at_exact_named_loci': s['abstentions_on_exact_named_reference_loci'],
                            'functional_accuracy': None})
    write(args.output / 'panel_summary.tsv', summary)
    report = {'analysis_type': 'FRESH_REFERENCE_CONCORDANCE_BY_PANEL',
              'source_audit_sha256': checksum(args.audit / 'manifest.json'),
              'equivalences_sha256': checksum(args.equivalences),
              'selection_sha256': checksum(args.selection),
              'external_accessions': sorted(external), 'original_accessions': sorted(accessions - external),
              'accuracy_superiority_established': False, 'independent_ground_truth': False,
              'scope': manifest['scope'], 'current_fusion_rules_version': manifest['current_fusion_rules_version'],
              'interpretation': 'External genomes were excluded from rule tuning; database overlap and reference independence are not excluded. Named assertions and exact/equivalent labels are concordance, not biological accuracy.'}
    (args.output / 'manifest.json').write_text(json.dumps(report, indent=2, sort_keys=True) + '\n')
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
