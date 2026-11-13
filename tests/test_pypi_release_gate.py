import importlib.util
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'tools'))
from pypi_release_gate import published_and_verified
from publish_verified_release import REQUIRED

SHA='a'*40


def api_fixture(*, draft=False, tag_sha=SHA, checks=True, release=True):
    def get(path, optional=False):
        if path.endswith('git/ref/heads/main'): return {'object': {'sha': SHA}}
        if 'actions/runs?' in path:
            return {'workflow_runs': [{'id':i, 'name':name, 'head_branch':'main', 'event':'push', 'head_sha':SHA, 'status':'completed', 'conclusion':'success' if checks else 'failure'} for i,name in enumerate(REQUIRED)]}
        if '/releases/tags/' in path:
            return {'tag_name':'v1.3.3', 'draft':draft, 'prerelease':False} if release else None
        if '/git/ref/tags/' in path: return {'object':{'type':'commit', 'sha':tag_sha}}
        raise AssertionError(path)
    return get


def test_publication_requires_matching_verified_public_release():
    assert published_and_verified('owner/repo', SHA, '1.3.3', api_fixture())
    for settings in ({'draft':True}, {'tag_sha':'b'*40}, {'checks':False}, {'release':False}):
        assert not published_and_verified('owner/repo', SHA, '1.3.3', api_fixture(**settings))
    assert not published_and_verified('owner/repo', SHA, '1.3.3rc1', api_fixture())
