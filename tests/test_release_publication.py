import importlib.util
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('release_gate', ROOT / 'tools/publish_verified_release.py')
gate = importlib.util.module_from_spec(spec)
spec.loader.exec_module(gate)
SHA = 'a' * 40
REPO = 'owner/project'


def fixture_api(runs, main=SHA):
    def get(path, optional=False):
        if path.endswith('git/ref/heads/main'):
            return {'object': {'sha': main}}
        return {'workflow_runs': runs}
    return get


def successful_runs():
    return [{'id': i, 'name': name, 'head_branch': 'main', 'event': 'push',
             'head_sha': SHA, 'status': 'completed', 'conclusion': 'success'}
            for i, name in enumerate(gate.REQUIRED, 1)]


def test_release_requires_every_check_on_current_main():
    runs = successful_runs()
    assert gate.ready(REPO, SHA, fixture_api(runs))
    assert not gate.ready(REPO, SHA, fixture_api(runs[:-1]))
    assert not gate.ready(REPO, SHA, fixture_api(runs, main='b' * 40))


@pytest.mark.parametrize('field,value', [('event', 'pull_request'), ('head_branch', 'feature'),
                                        ('head_sha', 'b' * 40), ('status', 'in_progress'),
                                        ('conclusion', 'failure')])
def test_release_rejects_wrong_or_failed_check(field, value):
    runs = successful_runs()
    runs[0][field] = value
    assert not gate.ready(REPO, SHA, fixture_api(runs))


def test_new_failed_run_supersedes_old_success():
    runs = successful_runs()
    runs.append(dict(runs[0], id=100, conclusion='failure'))
    assert not gate.ready(REPO, SHA, fixture_api(runs))


def test_existing_release_tag_cannot_be_retargeted():
    def get(path, optional=False):
        return {'object': {'type': 'commit', 'sha': 'b' * 40}}
    with pytest.raises(ValueError, match='cannot be reused'):
        gate.tag_matches(REPO, 'v1.3.0', SHA, get)
