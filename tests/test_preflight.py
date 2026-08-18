import json
import tempfile
import unittest
from pathlib import Path

from phagemine.evidence import EvidenceAdapterResult, EvidenceState
from phagemine.preflight import preflight_resources
from phagemine.resources import EvidenceResourceManager, ResourceType


class PreflightTests(unittest.TestCase):
    def test_evidence_state_normalizes_legacy_real_results(self):
        self.assertEqual(EvidenceAdapterResult("x", "REAL").state, EvidenceState.SUCCESS_NO_HIT)

    def test_explicit_missing_resource_fails_closed(self):
        with self.assertRaises(RuntimeError):
            preflight_resources({"PFAM": {"path": "/definitely/missing/Pfam-A.hmm", "required_tools": ["hmmscan"]}})

    def test_resource_validation_preserves_registered_status(self):
        with tempfile.TemporaryDirectory() as temp:
            registry = Path(temp) / "resources.json"
            db = Path(temp) / "fake.db"
            db.write_text("x")
            manager = EvidenceResourceManager(registry)
            manager.register("CUSTOM", ResourceType.CUSTOM, db)
            checked = manager.validate("CUSTOM")
            self.assertEqual(checked["status"], "READY")


if __name__ == "__main__":
    unittest.main()
