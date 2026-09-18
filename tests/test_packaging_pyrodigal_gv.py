import subprocess
import sys
import tomllib
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[1]


def test_automatic_provider_dependencies_are_core_project_requirements():
    project = tomllib.loads((ROOT / "pyproject.toml").read_text())
    requirements = set(project["project"]["dependencies"])
    assert "pyrodigal==3.7.1" in requirements
    assert "pyrodigal-gv==0.3.2" in requirements
    assert "pyrodigal-rv==0.1.0" in requirements
    assert "pyrodigal==3.7.1" not in set(project["project"]["optional-dependencies"]["test"])


@pytest.mark.integration
def test_clean_runtime_imports_for_automatic_provider_dependencies():
    for module in ("pyrodigal", "pyrodigal_gv"):
        result = subprocess.run(
            [sys.executable, "-c", f"import {module}; print({module}.__version__)"],
            cwd=ROOT, capture_output=True, text=True, check=False,
        )
        assert result.returncode == 0, result.stderr
        assert result.stdout.strip()


@pytest.mark.integration
def test_fresh_annotation_has_operational_pyrodigal_gv_observer(tmp_path):
    output = tmp_path / "run"
    command = [
        sys.executable, "-c",
        "from phagemine.cli import main; raise SystemExit(main(__import__('sys').argv[1:]))",
        "annotate", str(ROOT / "examples" / "demo_phage.fasta"),
        "--output", str(output),
    ]
    result = subprocess.run(command, cwd=ROOT, capture_output=True, text=True, check=False)
    assert result.returncode == 0, result.stderr
    combined = result.stdout + result.stderr
    assert "Prodigal-gv is required for the prodigal_gv provider" not in combined
    assert "Pyrodigal-gv observational gene prediction" in combined
    manifest = output / "gene_calls" / "gene_call_manifest.json"
    assert manifest.is_file()
    assert '"provider_id": "prodigal_gv"' in manifest.read_text()
