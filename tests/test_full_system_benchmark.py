import importlib.util
import json
from pathlib import Path
import sys

import pytest


def runner():
    directory = Path(__file__).resolve().parents[1] / 'evaluation/v1.3_validation/scripts'
    spec = importlib.util.spec_from_file_location('full_system_runner', directory / 'run_full_system.py')
    module = importlib.util.module_from_spec(spec)
    sys.path.insert(0, str(directory))
    try:
        spec.loader.exec_module(module)
    finally:
        sys.path.remove(str(directory))
    return module


def test_full_benchmark_rejects_missing_database_before_creating_output(tmp_path):
    module = runner()
    with pytest.raises(ValueError, match='Missing READY resources'):
        module.main(['--registry', str(tmp_path / 'missing.json'), '--output', str(tmp_path / 'result')])
    assert not (tmp_path / 'result').exists()


@pytest.mark.parametrize('failure', ['absent_adapter', 'unavailable_adapter', 'skipped_inphared', 'missing_pmfdb'])
def test_full_benchmark_rejects_incomplete_searches(tmp_path, failure):
    module = runner()
    adapters = [{'adapter': name, 'status': 'REAL'} for name in module.REQUIRED_ADAPTERS]
    comparative = {'inphared': {'status': 'COMPLETE'},
                   'pmfdb': {'reference_type': 'PMFDB', 'command': ['mmseqs', 'search']}}
    if failure == 'absent_adapter': adapters.pop()
    if failure == 'unavailable_adapter': adapters[0]['status'] = 'UNAVAILABLE'
    if failure == 'skipped_inphared': comparative['inphared']['status'] = 'SKIPPED'
    if failure == 'missing_pmfdb': comparative['pmfdb'] = {}
    (tmp_path / 'run_manifest.json').write_text(json.dumps({'evidence_adapters': adapters,
                                                          'comparative_analysis': comparative}))
    with pytest.raises(ValueError): module.verify_full_run(tmp_path)
