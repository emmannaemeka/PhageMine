"""Offline database diagnostics; freshness is never inferred from a file timestamp."""
from __future__ import annotations
from datetime import datetime, timezone
from pathlib import Path

from .resources import EvidenceResourceManager, validate_resource

ANNOTATION_TYPES = ('PFAM', 'VOGDB', 'SWISSPROT', 'PHROGS')


def health(*, manager=None, expected_versions=None, max_age_days=365, check_checksums=False):
    if max_age_days < 1:
        raise ValueError('max_age_days must be positive')
    expected_versions = expected_versions or {}
    if not isinstance(expected_versions, dict) or any(not isinstance(k, str) or not isinstance(v, str) or not v for k, v in expected_versions.items()):
        raise ValueError('Expected versions must map resource names or types to nonempty version strings')
    manager = manager or EvidenceResourceManager()
    now = datetime.now(timezone.utc)
    items = []
    for resource in manager.list():
        row = validate_resource(resource, check_checksum=check_checksums)
        warnings = []
        version = row.get('version')
        if not version or str(version).lower() in {'unknown', 'unavailable'}:
            warnings.append('VERSION_UNKNOWN: register the actual database release')
        expected = expected_versions.get(row['name'], expected_versions.get(row['resource_type']))
        if expected is not None and str(version) != expected:
            row['status'] = 'INVALID'
            row['validation_errors'].append(f'Expected database version {expected}; registered version is {version}')
        provenance = row.get('provenance') or {}
        date = provenance.get('release_date') or provenance.get('date_downloaded') or provenance.get('date_installed') or provenance.get('installed_at')
        if not date:
            primary = Path(row['path']).expanduser()
            manifest_path = (primary if primary.is_dir() else primary.parent) / 'install_manifest.json'
            if manifest_path.is_file():
                import json
                try: date = json.loads(manifest_path.read_text()).get('installed_at')
                except (OSError, ValueError, AttributeError): warnings.append('INVALID_INSTALL_MANIFEST: cannot assess installation date')
        date_basis = 'release_date' if provenance.get('release_date') else 'installation_date'
        if date:
            try:
                stamp = datetime.fromisoformat(str(date).replace('Z', '+00:00'))
                stamp = stamp.replace(tzinfo=timezone.utc) if stamp.tzinfo is None else stamp
                age = (now-stamp).total_seconds()/86400
                if age < 0: warnings.append('DATE_IN_FUTURE: check recorded database date')
                elif age > max_age_days: warnings.append(f'AGE_WARNING: recorded {date_basis} exceeds {max_age_days} days; check provider releases')
                row['recorded_age_days'] = round(age, 1)
            except (ValueError, TypeError): warnings.append('INVALID_DATE: use an ISO database date')
        else:
            warnings.append('FRESHNESS_UNKNOWN: no release/download/install date recorded')
        row.update(warnings=warnings, freshness_basis=date_basis if date else None,
                   checksum_verification='REQUESTED' if check_checksums else 'NOT_REQUESTED')
        items.append(row)
    missing = [kind for kind in ANNOTATION_TYPES if not any(r['resource_type'] == kind for r in items)]
    ambiguous = [kind for kind in ANNOTATION_TYPES if sum(r['resource_type'] == kind and r['status'] == 'READY' for r in items) > 1]
    unmet_expectations = [key for key in expected_versions if not any(key in {r['name'], r['resource_type']} for r in items)]
    problems = bool(unmet_expectations) or any(row['status'] != 'READY' for row in items) or bool(ambiguous)
    return {'schema_version': 1, 'status': 'NEEDS_ATTENTION' if problems else 'READY_WITH_LIMITED_EVIDENCE' if missing else 'READY',
            'resources': items, 'missing_annotation_resources': missing, 'ambiguous_resource_types': ambiguous,
            'unmet_version_expectations': unmet_expectations, 'freshness_checked_online': False, 'max_age_days': max_age_days,
            'interpretation': 'Age warnings are maintenance prompts, not proof of an outdated release. Missing resources reduce available annotation evidence.'}
