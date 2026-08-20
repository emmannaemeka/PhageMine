"""Tests for release-candidate desktop bundle validation."""
from __future__ import annotations

import importlib.util
from pathlib import Path


SCRIPT = Path(__file__).parents[1] / "scripts" / "audit_macos_bundle.py"
SPEC = importlib.util.spec_from_file_location("audit_macos_bundle", SCRIPT)
assert SPEC and SPEC.loader
AUDIT = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(AUDIT)


def test_macos_minimum_parser_ignores_dylib_compatibility_versions():
    output = """
Load command 1
      cmd LC_ID_DYLIB
  cmdsize 48
  current version 99.0.0
compatibility version 99.0.0
Load command 2
      cmd LC_BUILD_VERSION
  cmdsize 32
 platform 1
    minos 12.0
      sdk 15.0
Load command 3
      cmd LC_VERSION_MIN_MACOSX
  cmdsize 16
  version 10.13
      sdk 14.0
"""
    assert AUDIT.parse_minimum_versions(output) == [(12, 0), (10, 13)]
