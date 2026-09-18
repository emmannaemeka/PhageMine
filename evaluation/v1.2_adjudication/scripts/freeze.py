"""Read-only audit of the PR #14 snapshot. Never refresh frozen checksums."""
import csv
import hashlib
import json
import subprocess
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

REPO = Path(__file__).resolve().parents[3]
ROOT = Path(__file__).resolve().parents[1]
SOURCE = REPO / 'evaluation/v1.2_functional'
BASE_COMMIT = '4933d8fa532a0157c24c3c01dbf1d8adfa787ed7'
TOOLS = ('PhageMine', 'Pharokka')
PENDING = 'UNRESOLVED_REVIEW_REQUIRED'


def read(path, delimiter='\t'):
    with Path(path).open(newline='') as f:
        return list(csv.DictReader(f, delimiter=delimiter))


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def frozen_paths():
    paths = set()
    for manifest in ('docs/benchmark_v1.2/SHA256SUMS', 'evaluation/v1.2_functional/SHA256SUMS'):
        paths.add(manifest)
        for line in (REPO/manifest).read_text().splitlines():
            checksum, name = line.split('  ', 1)
            if digest(REPO/name) != checksum:
                raise ValueError('FROZEN CHECKSUM MISMATCH: '+name)
            paths.add(name)
    paths.add('src/phagemine/functional_benchmark.py')
    return sorted(paths)


def audit():
    paths = frozen_paths()
    # Compare against Git's immutable merge snapshot, not a newly refreshed list.
    for path in paths:
        blob = subprocess.check_output(['git', 'show', f'{BASE_COMMIT}:{path}'], cwd=REPO)
        if hashlib.sha256(blob).hexdigest() != digest(REPO/path):
            raise ValueError('FROZEN SNAPSHOT MODIFIED: '+path)
    rows = read(SOURCE/'functional_by_locus.tsv')
    exact = [r for r in rows if r['matching_mode']=='exact']
    counts = {}
    for tool, named, unresolved, correct in [('PhageMine',358,292,66),('Pharokka',341,38,303)]:
        c = Counter(r[tool+'_category'] for r in exact)
        actual = {'evaluable_loci':len(exact)-c['NON_EVALUABLE_REFERENCE'],
                  'named_assertions':len(exact)-c['NON_EVALUABLE_REFERENCE']-c['ABSTENTION'],
                  'unresolved':c[PENDING],
                  'automatic_agreements':c['EXACT_PRODUCT_AGREEMENT']+c['EQUIVALENT_FUNCTION']}
        if actual != dict(evaluable_loci=390,named_assertions=named,unresolved=unresolved,automatic_agreements=correct):
            raise ValueError('FROZEN COUNT MISMATCH: '+tool)
        counts[tool] = actual
    yields=read(REPO/'docs/benchmark_v1.2/tables/functional_yield_by_genome.tsv')
    recovered=[int(r['named_product_calls']) for r in yields if r['tool']=='PhageMine']
    if recovered != [144,56,50,95,8,48,38]:
        raise ValueError('Frozen named yield changed')
    if len(read(SOURCE/'phagemine_only_functional_assignments.tsv'))!=75 or len(read(SOURCE/'pharokka_only_functional_assignments.tsv'))!=34:
        raise ValueError('Frozen exclusive assignments changed')
    unique={(r['accession'],r['start'],r['end'],r['strand'],r['reference_id'],tool)
            for r in rows for tool in TOOLS if r[tool+'_category']==PENDING}
    if len(unique)!=359:
        raise ValueError('Frozen unique review-unit count changed')
    return dict(verified_file_count=len(paths), counts=counts, named_yield=recovered,
                unique_unresolved_judgments=len(unique), checksums={p:digest(REPO/p) for p in paths})


def create_manifest():
    data=audit()
    target=ROOT/'frozen/benchmark_manifest.json'
    manifest=dict(benchmark_version='v1.2-functional-PR14-pre-adjudication',
        repository_commit=BASE_COMMIT,
        freeze_time_utc='2026-09-18T03:20:47Z',
        freeze_time_basis='PR14 merge timestamp: identifies frozen source snapshot, not a retrospective statistical preregistration',
        audit_time_utc=datetime.now(timezone.utc).isoformat(),
        reference_genomes=read(SOURCE/'reference_provenance.tsv'),
        scoring_script_checksums={p:sha for p,sha in data['checksums'].items() if p.endswith('.py')},
        prediction_file_checksums={p:sha for p,sha in data['checksums'].items() if '/predictions/' in p},
        exclusions=dict(compound_reference_locations=11,exact_non_evaluable_reference_loci=198,
            unmatched_prediction_loci=216,phagemine_only_non_evaluable=21,phagemine_only_no_exact_match=18),
        predicted_cds_per_tool=804,named_product_yield_per_tool={'PhageMine':439,'Pharokka':398},
        **data)
    with target.open('x') as f:
        json.dump(manifest,f,indent=2,sort_keys=True);f.write('\n')
    return manifest


if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser();parser.add_argument('--create-manifest',action='store_true')
    args=parser.parse_args()
    result=create_manifest() if args.create_manifest else audit()
    print(json.dumps({k:v for k,v in result.items() if k not in {'checksums','prediction_file_checksums','scoring_script_checksums'}},indent=2))
