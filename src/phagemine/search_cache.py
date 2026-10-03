"""Opt-in persistent cache for sequence-only evidence adapters."""
from __future__ import annotations
import copy
from dataclasses import asdict
import hashlib
import json
import os
from pathlib import Path
import shutil
import tempfile

from . import __version__
from .artifact_ops import checksum
from .evidence import EvidenceAdapterResult
from .models import Evidence, EvidenceLevel

SUPPORTED = {'PfamHMMAdapter', 'VOGHMMAdapter', 'SwissProtEvidenceAdapter', 'PHROGSMMseqsAdapter', 'PHROGSPyHMMERAdapter'}
SUCCESS = {'REAL', 'SUCCESS_WITH_HITS', 'SUCCESS_NO_HIT'}


def _digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':'), default=str).encode()).hexdigest()


def fingerprint(adapter):
    """Hash implementation, parameters, executable, database and index bytes."""
    dependencies = {}
    for name, value in vars(adapter).items():
        if value is None: continue
        if isinstance(value, Path) or name.endswith('_path') or name in {'hmmscan', 'diamond', 'mmseqs'}:
            path = Path(str(value)).expanduser()
            if name in {'hmmscan', 'diamond', 'mmseqs'}:
                path = Path(shutil.which(str(value)) or str(path))
            if not path.is_file(): raise ValueError(f'Cache dependency is not a file: {name}')
            files = {path}
            if name.endswith('_path'):
                # MMseqs prefixes and HMMER pressed companions are inputs too.
                files.update(p for p in path.parent.iterdir() if p.name.startswith(path.name + '.') or p.name.startswith(path.name + '_'))
            for candidate in sorted(files):
                if candidate.is_file(): dependencies[str(candidate.resolve())] = checksum(candidate)
    implementation = {str(p.relative_to(Path(__file__).parent)): checksum(p) for p in sorted(Path(__file__).parent.rglob('*.py'))}
    return _digest({'schema': 1, 'phagemine': __version__, 'adapter': adapter.__class__.__module__ + '.' + adapter.__class__.__qualname__,
                    'configuration': vars(adapter), 'provenance': adapter.provenance(), 'files': dependencies, 'implementation': implementation})


def _bind(evidence, protein_id):
    evidence = copy.deepcopy(evidence)
    evidence.provenance['protein_id'] = protein_id
    if 'query_protein_id' in evidence.metrics: evidence.metrics['query_protein_id'] = protein_id
    return evidence


class PersistentEvidenceCache:
    def __init__(self, directory):
        self.directory = Path(directory).expanduser()
        if self.directory.is_symlink(): raise ValueError('Cache directory cannot be a symlink')
        self.directory.mkdir(parents=True, exist_ok=True)
        self.stats = []

    def analyze(self, adapter, proteins):
        if adapter.name not in SUPPORTED or not proteins or not adapter.available():
            return adapter.analyze(proteins)
        key = fingerprint(adapter)
        unique = {}
        for protein in proteins:
            unique.setdefault(hashlib.sha256(protein.sequence.encode()).hexdigest(), protein)
        entries = {}; missing = []; corrupt = 0
        for seqhash, protein in unique.items():
            path = self.directory / (key + '-' + seqhash + '.json')
            try:
                if path.is_symlink(): raise ValueError('Symlink cache entry')
                payload = json.loads(path.read_text())
                expected = payload.pop('sha256')
                if expected != _digest(payload) or payload['adapter_key'] != key or payload['sequence_sha256'] != seqhash or payload['status'] not in SUCCESS:
                    raise ValueError('Cache integrity or identity mismatch')
                evidence = [Evidence(**{**row, 'level': EvidenceLevel(row['level'])}) for row in payload['evidence']]
                entries[seqhash] = (payload, evidence)
            except FileNotFoundError: missing.append(protein)
            except (OSError, ValueError, KeyError, TypeError):
                corrupt += 1; missing.append(protein)
        fresh = adapter.analyze(missing) if missing else None
        # A failed new search cannot be presented as a successful mixed result.
        if fresh is not None and fresh.status not in SUCCESS:
            self.stats.append({'adapter': adapter.name, 'reused': len(entries), 'searched': len(missing), 'status': 'SEARCH_FAILED_NOT_CACHED'})
            return fresh
        if fresh is not None:
            if fingerprint(adapter) != key: raise RuntimeError('Search inputs changed during execution; evidence was not cached')
            known = {p.protein_id for p in missing}
            if any((e.provenance.get('protein_id') or e.metrics.get('query_protein_id')) not in known for e in fresh.evidence):
                raise ValueError('Search evidence contains an unknown query ID')
            for protein in missing:
                seqhash = hashlib.sha256(protein.sequence.encode()).hexdigest()
                hits = [e for e in fresh.evidence if (e.provenance.get('protein_id') or e.metrics.get('query_protein_id')) == protein.protein_id]
                payload = {'schema_version': 1, 'adapter_key': key, 'sequence_sha256': seqhash,
                           'status': fresh.status, 'provenance': fresh.provenance, 'message': fresh.message,
                           'evidence': [asdict(e) for e in hits]}
                payload['sha256'] = _digest(payload)
                target = self.directory / (key + '-' + seqhash + '.json')
                with tempfile.NamedTemporaryFile(mode='w', dir=self.directory, delete=False) as handle:
                    temp = Path(handle.name); json.dump(payload, handle, sort_keys=True)
                try: os.replace(temp, target)
                finally: temp.unlink(missing_ok=True)
                entries[seqhash] = (payload, hits)
        evidence = []
        for protein in proteins:
            seqhash = hashlib.sha256(protein.sequence.encode()).hexdigest()
            evidence.extend(_bind(e, protein.protein_id) for e in entries[seqhash][1])
        template = fresh.provenance if fresh else next(iter(entries.values()))[0]['provenance']
        provenance = {**template, 'cache': {'adapter_key': key, 'reused_sequences': len(unique)-len(missing),
                       'searched_sequences': len(missing), 'corrupt_entries_recomputed': corrupt,
                       'cached_commands_are_historical': True}}
        self.stats.append({'adapter': adapter.name, **provenance['cache'], 'status': 'SUCCESS'})
        return EvidenceAdapterResult(adapter.name, 'REAL', evidence, provenance)
