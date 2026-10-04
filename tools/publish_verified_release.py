#!/usr/bin/env python3
"""Publish immutable GitHub artifacts only after all checks on current main pass."""
from __future__ import annotations
import argparse
import json
import os
from pathlib import Path
import re
import subprocess
import hashlib
import zipfile
import tarfile

REQUIRED = ('PhageMine CI', 'External-tool integration', 'Container annotation smoke')


def api(path, optional=False):
    result = subprocess.run(['gh', 'api', path], capture_output=True, text=True)
    if result.returncode:
        if optional and '404' in result.stderr: return None
        raise RuntimeError('GitHub API check failed: ' + result.stderr)
    return json.loads(result.stdout)


def ready(repo, sha, get=api):
    if not re.fullmatch(r'[A-Za-z0-9_.-]+/[A-Za-z0-9_.-]+', repo) or not re.fullmatch(r'[a-f0-9]{40}', sha):
        raise ValueError('Invalid repository or commit identity')
    if get(f'repos/{repo}/git/ref/heads/main')['object']['sha'] != sha: return False
    runs = get(f'repos/{repo}/actions/runs?head_sha={sha}&per_page=100')['workflow_runs']
    for name in REQUIRED:
        candidates = [run for run in runs if run['name'] == name and run['head_branch'] == 'main' and
                      run['event'] == 'push' and run['head_sha'] == sha]
        if not candidates: return False
        latest = max(candidates, key=lambda run: run['id'])
        if latest['status'] != 'completed' or latest['conclusion'] != 'success': return False
    return True


def tag_matches(repo, tag, sha, get=api):
    ref = get(f'repos/{repo}/git/ref/tags/{tag}', optional=True)
    if ref is None: return
    obj = ref['object']
    if obj['type'] == 'tag': obj = get(f'repos/{repo}/git/tags/{obj["sha"]}')['object']
    if obj['type'] != 'commit' or obj['sha'] != sha:
        raise ValueError('Existing release tag points to a different commit; published versions cannot be reused')


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repo', required=True); parser.add_argument('--sha', required=True)
    parser.add_argument('--check', action='store_true')
    parser.add_argument('--evidence', action='store_true')
    parser.add_argument('--output', default='release-validation.json')
    parser.add_argument('--tag'); parser.add_argument('--notes'); parser.add_argument('--artifacts', default='dist')
    args=parser.parse_args(argv)
    checked = ready(args.repo, args.sha)
    if args.check:
        if os.environ.get('GITHUB_OUTPUT'):
            with open(os.environ['GITHUB_OUTPUT'], 'a') as handle: handle.write('ready=' + ('true' if checked else 'false') + '\n')
        print('READY' if checked else 'SKIPPED_CHECKS_PENDING_FAILED_OR_MAIN_MOVED'); return 0
    if not checked: raise RuntimeError('Release gates changed; publication stopped')
    if args.evidence:
        runs = api(f'repos/{args.repo}/actions/runs?head_sha={args.sha}&per_page=100')['workflow_runs']
        selected = [max((r for r in runs if r['name'] == name and r['head_branch'] == 'main'
                         and r['event'] == 'push' and r['head_sha'] == args.sha), key=lambda r: r['id'])
                    for name in REQUIRED]
        payload = {'release_commit': args.sha, 'status': 'PASS',
                   'checks': [{key: r[key] for key in ('name', 'id', 'html_url', 'conclusion')} for r in selected],
                   'full_system_benchmark': 'https://github.com/emmannaemeka/PhageMine/actions/runs/37202973780',
                   'scientific_accuracy_superiority_established': False}
        Path(args.output).write_text(json.dumps(payload, indent=2) + '\n')
        return 0
    if not args.tag or not re.fullmatch(r'v[0-9]+\.[0-9]+\.[0-9]+', args.tag): raise ValueError('Stable release tag required')
    tag_matches(args.repo, args.tag, args.sha)
    existing = api(f'repos/{args.repo}/releases/tags/{args.tag}', optional=True)
    if existing is not None and not existing['draft']:
        print('Existing release retained unchanged'); return 0
    notes=Path(args.notes or '')
    if not notes.is_file(): raise ValueError('Release notes are missing')
    artifacts=sorted(p for p in Path(args.artifacts).iterdir() if p.is_file())
    if not any(p.suffix == '.whl' for p in artifacts) or not any(p.name.endswith('.tar.gz') for p in artifacts) or not any(p.name == 'SHA256SUMS' for p in artifacts):
        raise ValueError('Wheel, source distribution and checksums are required')
    version = args.tag[1:]
    wheel = Path(args.artifacts) / f'phagemine-{version}-py3-none-any.whl'
    source = Path(args.artifacts) / f'phagemine-{version}.tar.gz'
    if not wheel.is_file() or not source.is_file():
        raise ValueError('Artifacts must match the stable release version')
    with zipfile.ZipFile(wheel) as archive:
        metadata = archive.read(f'phagemine-{version}.dist-info/METADATA').decode()
    with tarfile.open(source) as archive:
        member = archive.extractfile(f'phagemine-{version}/PKG-INFO')
        if member is None: raise ValueError('Source metadata is missing')
        source_metadata = member.read().decode()
    if not all(f'\nVersion: {version}\n' in text for text in (metadata, source_metadata)):
        raise ValueError('Artifact metadata does not match the release')
    checksums = (Path(args.artifacts) / 'SHA256SUMS').read_text().splitlines()
    expected = {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in artifacts if p.name != 'SHA256SUMS'}
    actual = {}
    for line in checksums:
        digest, name = line.split('  ', 1)
        if name in actual: raise ValueError('Duplicate checksum entry')
        actual[name] = digest
    if actual != expected: raise ValueError('Artifact checksum validation failed')
    # Recheck after inspecting artifacts, before the external write.
    if not ready(args.repo, args.sha): raise RuntimeError('Main moved or checks changed before publication')
    if existing is None:
        subprocess.run(['gh', 'release', 'create', args.tag, '--repo', args.repo, '--target', args.sha,
                        '--title', 'PhageMine ' + args.tag, '--notes-file', str(notes), '--draft'], check=True)
    subprocess.run(['gh', 'release', 'upload', args.tag, '--repo', args.repo, '--clobber', *map(str, artifacts)], check=True)
    if not ready(args.repo, args.sha): raise RuntimeError('Release gates changed before publishing the draft')
    subprocess.run(['gh', 'release', 'edit', args.tag, '--repo', args.repo, '--draft=false', '--latest'], check=True)
    return 0


if __name__ == '__main__': raise SystemExit(main())
