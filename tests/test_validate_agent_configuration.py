#!/usr/bin/env python3
"""Behavioral tests for the read-only agent configuration validator."""

from __future__ import annotations

import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest


REPO_ROOT = Path(__file__).resolve().parents[1]
VALIDATOR = REPO_ROOT / "scripts" / "validate_agent_configuration.py"


class AgentConfigurationValidatorTests(unittest.TestCase):
    maxDiff = None

    def run_validator(self, *arguments: str) -> subprocess.CompletedProcess[str]:
        return subprocess.run(
            [sys.executable, str(VALIDATOR), *arguments, "--format", "json"],
            cwd=REPO_ROOT,
            text=True,
            capture_output=True,
            check=False,
        )

    def test_minimal_example_passes_with_pending_decision_warnings(self) -> None:
        result = self.run_validator(
            "--samples",
            "agent/examples/minimal_samples.tsv",
            "--parameters",
            "agent/examples/minimal_parameters.json",
            "--status",
            "agent/examples/minimal_status.json",
        )
        self.assertEqual(result.returncode, 0, result.stderr or result.stdout)
        payload = json.loads(result.stdout)
        self.assertEqual(payload["status"], "PASS_WITH_WARNINGS")
        self.assertEqual(payload["counts"]["errors"], 0)
        self.assertEqual(payload["counts"]["warnings"], 8)
        self.assertEqual(payload["counts"]["passes"], 4)

    def test_strict_mode_rejects_pending_decision_warnings(self) -> None:
        result = self.run_validator(
            "--samples",
            "agent/examples/minimal_samples.tsv",
            "--parameters",
            "agent/examples/minimal_parameters.json",
            "--status",
            "agent/examples/minimal_status.json",
            "--strict",
        )
        self.assertEqual(result.returncode, 1)
        payload = json.loads(result.stdout)
        self.assertEqual(payload["status"], "PASS_WITH_WARNINGS")
        self.assertEqual(payload["counts"]["errors"], 0)

    def test_unfilled_templates_are_rejected(self) -> None:
        result = self.run_validator(
            "--samples",
            "config/templates/samples.template.tsv",
            "--parameters",
            "config/templates/analysis_parameters.agent.template.json",
            "--status",
            "agent/analysis_status.template.json",
        )
        self.assertEqual(result.returncode, 1)
        payload = json.loads(result.stdout)
        self.assertEqual(payload["status"], "FAIL")
        self.assertGreater(payload["counts"]["errors"], 0)
        self.assertTrue(
            any("placeholder" in message for message in payload["errors"])
        )

    def test_file_check_rejects_missing_example_fastq_files(self) -> None:
        result = self.run_validator(
            "--samples",
            "agent/examples/minimal_samples.tsv",
            "--parameters",
            "agent/examples/minimal_parameters.json",
            "--check-files",
        )
        self.assertEqual(result.returncode, 1)
        payload = json.loads(result.stdout)
        self.assertEqual(payload["status"], "FAIL")
        self.assertTrue(
            any("missing or empty" in message for message in payload["errors"])
        )

    def test_complete_stage_requires_pass_validation_evidence(self) -> None:
        source = REPO_ROOT / "agent" / "examples" / "minimal_status.json"
        status = json.loads(source.read_text(encoding="utf-8"))
        status["human_decisions"]["study_scope"] = {
            "status": "approved",
            "approved_by": "example_reviewer",
            "approved_at": "2026-01-01T00:00:00Z",
            "decision": "Configuration-only test decision.",
        }
        status["stages"][0]["status"] = "complete"

        with tempfile.TemporaryDirectory() as temporary_directory:
            status_path = Path(temporary_directory) / "status.json"
            status_path.write_text(json.dumps(status), encoding="utf-8")
            result = self.run_validator(
                "--samples",
                "agent/examples/minimal_samples.tsv",
                "--parameters",
                "agent/examples/minimal_parameters.json",
                "--status",
                str(status_path),
            )

        self.assertEqual(result.returncode, 1)
        payload = json.loads(result.stdout)
        self.assertTrue(
            any(
                "no recorded PASS validation" in message
                for message in payload["errors"]
            )
        )


if __name__ == "__main__":
    unittest.main()
