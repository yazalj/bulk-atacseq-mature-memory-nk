from __future__ import annotations

import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "run_synthetic_demo.py"
EXPECTED = ROOT / "demo" / "expected"
INPUT = ROOT / "demo" / "input"


class SyntheticDemoTests(unittest.TestCase):
    def test_read_only_check_matches_expected_outputs(self) -> None:
        result = subprocess.run(
            [sys.executable, str(SCRIPT), "--check"],
            cwd=ROOT,
            check=True,
            capture_output=True,
            text=True,
        )
        self.assertIn("Synthetic demo: PASS (5 of 6 peaks retained)", result.stdout)

    def test_output_attempt_is_byte_identical_to_expected(self) -> None:
        with tempfile.TemporaryDirectory(prefix="bulk_atacseq_demo_test_") as directory:
            output = Path(directory) / "attempt1"
            subprocess.run(
                [sys.executable, str(SCRIPT), "--output-dir", str(output)],
                cwd=ROOT,
                check=True,
                capture_output=True,
                text=True,
            )
            for expected in EXPECTED.iterdir():
                self.assertEqual(expected.read_bytes(), (output / expected.name).read_bytes())

    def test_reordered_count_columns_are_rejected(self) -> None:
        with tempfile.TemporaryDirectory(prefix="bulk_atacseq_demo_test_") as directory:
            temporary = Path(directory)
            bad_counts = temporary / "reordered_counts.tsv"
            original = (INPUT / "synthetic_counts.tsv").read_text(encoding="utf-8")
            bad_counts.write_text(
                original.replace(
                    "peak_id\tD1_A\tD1_B\tD2_A\tD2_B\tD3_A\tD3_B",
                    "peak_id\tD1_B\tD1_A\tD2_A\tD2_B\tD3_A\tD3_B",
                    1,
                ),
                encoding="utf-8",
            )
            result = subprocess.run(
                [
                    sys.executable,
                    str(SCRIPT),
                    "--counts",
                    str(bad_counts),
                    "--output-dir",
                    str(temporary / "attempt1"),
                ],
                cwd=ROOT,
                capture_output=True,
                text=True,
            )
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("exactly match ordered sample IDs", result.stderr)


if __name__ == "__main__":
    unittest.main()
