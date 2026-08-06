#!/usr/bin/env python3
"""Tests for the versioned NeoQC technical QC decision engine."""

from __future__ import annotations

import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "scripts"
sys.path.insert(0, str(SCRIPTS))

from qc_rules import evaluate_directory, write_evaluation
try:
    from test_plot_results import write_fixture_set
except ModuleNotFoundError:
    from tests.test_plot_results import write_fixture_set


class QcRulesTest(unittest.TestCase):
    def test_complete_fixture_produces_explainable_statuses(self) -> None:
        with tempfile.TemporaryDirectory(prefix="neoqc-rules-") as temporary:
            input_dir = Path(temporary)
            write_fixture_set(input_dir, ("R1",))

            result = evaluate_directory(input_dir)
            self.assertEqual(result["schema_version"], 1)
            self.assertEqual(result["ruleset"]["id"], "fastqc-compatible-v1")
            self.assertEqual(len(result["ruleset"]["sha256"]), 64)
            self.assertEqual(len(result["evaluations"]), 8)
            by_metric = {item["metric_id"]: item for item in result["evaluations"]}
            self.assertEqual(by_metric["per_base_quality"]["qc_status"], "pass")
            self.assertEqual(by_metric["sequence_length_distribution"]["qc_status"], "warning")
            duplication = by_metric["sequence_duplication_levels"]
            self.assertEqual(duplication["qc_status"], "fail")
            self.assertAlmostEqual(
                duplication["observations"]["deduplicated_remaining_percent"],
                100 * 13.4 / 51.5,
                places=6,
            )
            self.assertIn("threshold", duplication["reasons"][0])
            self.assertEqual(result["summary"]["overall_status"], "fail")

    def test_old_mean_only_quality_tsv_is_not_evaluated(self) -> None:
        with tempfile.TemporaryDirectory(prefix="neoqc-old-quality-") as temporary:
            input_dir = Path(temporary)
            (input_dir / "per_cycle_R1.tsv").write_text(
                "cycle\tmean_quality\n1\t35\n2\t34\n", encoding="utf-8"
            )

            result = evaluate_directory(input_dir)
            quality = next(
                item for item in result["evaluations"] if item["metric_id"] == "per_base_quality"
            )
            self.assertEqual(quality["qc_status"], "not_evaluated")
            self.assertEqual(quality["reasons"][0]["code"], "evaluation.data_invalid")
            self.assertNotEqual(result["summary"]["overall_status"], "pass")

    def test_per_sequence_quality_threshold_boundaries(self) -> None:
        expected = ((27, "pass"), (26, "warning"), (20, "warning"), (19, "fail"))
        for mode, status in expected:
            with self.subTest(mode=mode), tempfile.TemporaryDirectory(
                prefix="neoqc-boundary-"
            ) as temporary:
                input_dir = Path(temporary)
                (input_dir / "per_sequence_quality_R1.tsv").write_text(
                    f"mean_quality\tread_count\n{mode}\t100\n", encoding="utf-8"
                )
                result = evaluate_directory(input_dir)
                quality = next(
                    item
                    for item in result["evaluations"]
                    if item["metric_id"] == "per_sequence_quality"
                )
                self.assertEqual(quality["qc_status"], status)

    def test_writer_is_atomic_and_cli_is_standalone(self) -> None:
        with tempfile.TemporaryDirectory(prefix="neoqc-evaluation-cli-") as temporary:
            input_dir = Path(temporary)
            write_fixture_set(input_dir, ("R1",))
            output = input_dir / "custom evaluation.json"
            output.write_text("old", encoding="utf-8")

            written = write_evaluation(input_dir, output)
            self.assertEqual(written, output.resolve())
            self.assertEqual(json.loads(output.read_text(encoding="utf-8"))["schema_version"], 1)
            result = subprocess.run(
                [
                    sys.executable,
                    str(SCRIPTS / "evaluate_qc.py"),
                    str(input_dir),
                    "--output",
                    str(output),
                ],
                cwd=ROOT,
                text=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                check=False,
            )
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertIn("QC evaluation:", result.stdout)


if __name__ == "__main__":
    unittest.main()
