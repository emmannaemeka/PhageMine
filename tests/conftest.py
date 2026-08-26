"""Test isolation for persistent PhageMine state.

Every test receives a private resource registry.  This is especially important
on macOS, where XDG_CONFIG_HOME is not used by the production registry path.
"""
from __future__ import annotations

import pytest


@pytest.fixture(autouse=True)
def isolate_phagemine_registry(tmp_path, monkeypatch):
    registry = tmp_path / "phagemine-test-state" / "resources.json"
    monkeypatch.setenv("PHAGEMINE_REGISTRY_PATH", str(registry))
    return registry
