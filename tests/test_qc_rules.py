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
    @staticmethod
    def write_native_duplication_summary(
        directory: Path,
        *,
        remaining: float,
    ) -> None:
        (directory / "sequence_duplication_summary_R1.tsv").write_text(
            "source_kind\talgorithm\tsource_fastq\tprefix_length\ttotal_reads\t"
            "unique_sequences\tdeduplicated_remaining_percent\n"
            "native_fastq\tneoqc-exact-prefix-v1\tsample_R1.fastq.gz\t50\t"
            f"10\t6\t{remaining}\n",
            encoding="utf-8",
        )

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

    def test_native_duplication_summary_supports_profiles_without_level_one(self) -> None:
        with tempfile.TemporaryDirectory(prefix="neoqc-native-duplication-") as temporary:
            input_dir = Path(temporary)
            (input_dir / "sequence_duplication_levels_R1.tsv").write_text(
                "duplication_level\ttotal_sequences_percent\t"
                "deduplicated_sequences_percent\n"
                "2\t40\t70\n4\t60\t30\n",
                encoding="utf-8",
            )
            self.write_native_duplication_summary(input_dir, remaining=60)

            result = evaluate_directory(input_dir)
            duplication = next(
                item
                for item in result["evaluations"]
                if item["metric_id"] == "sequence_duplication_levels"
            )
            self.assertEqual(duplication["qc_status"], "warning")
            self.assertEqual(
                duplication["observations"]["deduplicated_remaining_percent"], 60
            )

    def test_incomplete_native_duplication_transaction_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory(prefix="neoqc-incomplete-duplication-") as temporary:
            input_dir = Path(temporary)
            (input_dir / "sequence_duplication_levels_R1.tsv").write_text(
                "duplication_level\ttotal_sequences_percent\t"
                "deduplicated_sequences_percent\n1\t100\t100\n",
                encoding="utf-8",
            )
            (input_dir / "sequence_duplication_R1.incomplete").write_text(
                "incomplete\n", encoding="utf-8"
            )

            result = evaluate_directory(input_dir)
            duplication = next(
                item
                for item in result["evaluations"]
                if item["metric_id"] == "sequence_duplication_levels"
            )
            self.assertEqual(duplication["qc_status"], "not_evaluated")
            self.assertEqual(duplication["reasons"][0]["code"], "evaluation.data_invalid")
            self.assertIn("incomplete", duplication["reasons"][0]["message"])

    def test_bounded_prototype_summary_remains_readable(self) -> None:
        with tempfile.TemporaryDirectory(prefix="neoqc-bounded-compatibility-") as temporary:
            input_dir = Path(temporary)
            (input_dir / "sequence_duplication_levels_R1.tsv").write_text(
                "duplication_level\ttotal_sequences_percent\t"
                "deduplicated_sequences_percent\n1\t50\t100\n2\t50\t0\n",
                encoding="utf-8",
            )
            (input_dir / "sequence_duplication_summary_R1.tsv").write_text(
                "source_kind\talgorithm\tsource_fastq\tprefix_length\t"
                "max_tracked_unique\ttotal_reads\ttracked_unique_sequences\t"
                "count_at_unique_limit\tsampling_limited\t"
                "deduplicated_remaining_percent\n"
                "native_fastq\tfastqc-compatible-bounded-v1\tsample_R1.fastq.gz\t"
                "50\t100000\t200000\t100000\t100250\ttrue\t50\n",
                encoding="utf-8",
            )

            result = evaluate_directory(input_dir)
            duplication = next(
                item
                for item in result["evaluations"]
                if item["metric_id"] == "sequence_duplication_levels"
            )
            self.assertEqual(duplication["qc_status"], "warning")
            self.assertEqual(
                duplication["observations"]["deduplicated_remaining_percent"], 50
            )

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
