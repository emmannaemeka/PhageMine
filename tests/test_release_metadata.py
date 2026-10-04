import re
import json
from pathlib import Path

import phagemine


ROOT = Path(__file__).resolve().parents[1]


def test_release_versions_are_synchronized():
    pyproject = (ROOT / "pyproject.toml").read_text()
    citation = (ROOT / "CITATION.cff").read_text()
    project_version = re.search(r'^version = "([^"]+)"$', pyproject, re.MULTILINE)
    citation_version = re.search(r"^version: ([^\s]+)$", citation, re.MULTILINE)
    assert project_version is not None
    assert citation_version is not None
    assert project_version.group(1) == phagemine.__version__ == citation_version.group(1)
    assert json.loads((ROOT / '.zenodo.json').read_text())['version'] == phagemine.__version__
    recipe = (ROOT / 'packaging/bioconda/recipe/meta.yaml').read_text()
    assert f'set version = "{phagemine.__version__}"' in recipe
    if re.fullmatch(r'\d+\.\d+\.\d+', phagemine.__version__):
        assert (ROOT / 'docs/releases' / f'v{phagemine.__version__}.md').is_file()
