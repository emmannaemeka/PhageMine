import json
import tempfile
import unittest
from unittest.mock import patch
from pathlib import Path

from phagemine.evidence import EvidenceAdapterResult, EvidenceState
from phagemine.preflight import doctor, doctor_text, executable_status, preflight_resources
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

    def test_doctor_text_gives_actionable_database_install_commands(self):
        payload = {
            "phagemine_version": "1.0.0",
            "executables": [],
            "resources": [],
            "capabilities": {},
            "recommendations": [{
                "action": "INSTALL_EVIDENCE_DATABASES",
                "commands": ["phagemine databases install pfam"],
                "all_command": "phagemine databases install --all",
                "verify_command": "phagemine doctor",
            }],
        }
        text = doctor_text(payload)
        self.assertIn("phagemine databases install pfam", text)
        self.assertIn("phagemine databases install --all", text)
        self.assertIn("phagemine doctor", text)

    @patch("phagemine.preflight.os.access", return_value=True)
    @patch("phagemine.preflight.shutil.which", return_value="/usr/local/bin/diamond")
    @patch("phagemine.preflight.subprocess.run")
    def test_executable_linker_failure_is_broken(self, run, _which, _access):
        run.return_value.returncode = 1
        run.return_value.stdout = ""
        run.return_value.stderr = "dyld: Symbol not found: broken_symbol"
        result = executable_status("diamond")
        self.assertEqual(result["status"], "BROKEN")
        self.assertIn("Symbol not found", result["diagnostic"])

    @patch("phagemine.preflight.os.access", return_value=True)
    @patch("phagemine.preflight.shutil.which", return_value="/usr/local/bin/diamond")
    @patch("phagemine.preflight.subprocess.run")
    def test_executable_requires_successful_probe(self, run, _which, _access):
        run.return_value.returncode = 0
        run.return_value.stdout = "diamond version 2.2.5"
        run.return_value.stderr = ""
        result = executable_status("diamond")
        self.assertEqual(result["status"], "READY")
        self.assertEqual(result["version"], "diamond version 2.2.5")

    def test_unknown_future_resource_type_does_not_crash_registry(self):
        with tempfile.TemporaryDirectory() as temp:
            registry = Path(temp) / "resources.json"
            db = Path(temp) / "future.db"
            db.write_text("x")
            registry.write_text(json.dumps({
                "format_version": "2",
                "resources": {
                    "Future": {"name": "Future", "resource_type": "FUTURE_DB", "path": str(db)}
                },
            }))
            rows = EvidenceResourceManager(registry).list()
            self.assertEqual(rows[0]["resource_type"], "FUTURE_DB")
            self.assertEqual(rows[0]["status"], "READY")

    def test_doctor_text_reports_broken_tool_and_does_not_reinstall_registered_database(self):
        payload = {
            "phagemine_version": "1.0.1",
            "executables": [{"name": "mash", "status": "MISSING", "version": None,
                             "diagnostic": None}],
            "resources": [{"resource_type": "INPHARED_GENOMES", "status": "INVALID"}],
            "capabilities": {"WHOLE_GENOME_REFERENCE_COMPARISON": "UNAVAILABLE"},
            "recommendations": [{
                "action": "INSTALL_OR_REPAIR_EXECUTABLES",
                "command": "conda install --channel conda-forge --channel bioconda --strict-channel-priority mash",
                "verify_command": "phagemine doctor",
            }],
        }
        output = doctor_text(payload)
        self.assertIn("INPHARED genomes  BLOCKED/INVALID", output)
        self.assertIn("conda install", output)
        self.assertNotIn("phagemine databases install inphared", output)

    @patch("phagemine.preflight.executable_status")
    @patch("phagemine.preflight.EvidenceResourceManager.validate_all")
    def test_doctor_recommends_mash_not_inphared_redownload(self, validate_all, status):
        validate_all.return_value = [{
            "name": "INPHARED-Genomes",
            "resource_type": "INPHARED_GENOMES",
            "status": "INVALID",
            "validation_errors": ["required executable unavailable: mash"],
        }]

        def tool(name):
            return {"name": name, "status": "MISSING" if name == "mash" else "READY",
                    "path": None, "version": "test" if name != "mash" else None,
                    "diagnostic": None}

        status.side_effect = tool
        actions = {item["action"]: item for item in doctor()["recommendations"]}
        self.assertIn("INSTALL_OR_REPAIR_EXECUTABLES", actions)
        self.assertIn("mash", actions["INSTALL_OR_REPAIR_EXECUTABLES"]["command"])
        self.assertNotIn("REPAIR_EVIDENCE_DATABASES", actions)


if __name__ == "__main__":
    unittest.main()
