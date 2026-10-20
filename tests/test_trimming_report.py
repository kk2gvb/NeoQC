#!/usr/bin/env python3
"""Tests for the standalone before/after trimming HTML report."""

from __future__ import annotations

import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch


ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "scripts"
sys.path.insert(0, str(SCRIPTS))

from trimming_report import REPORT_FILENAME, _render_parameters, generate_trimming_report
from trimming_report import _plain_value
import plot_results


SVG = (
    '<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 20 10">'
    '<defs><path id="curve" d="M0 5h20"/></defs><use href="#curve"/></svg>'
)


def make_result(root: Path) -> None:
    for stage, reads, status in (("before", "1200", "warning"), ("after", "1180", "pass")):
        result_dir = root / "qc" / stage
        plot_dir = root / "plots" / stage
        result_dir.mkdir(parents=True)
        plot_dir.mkdir(parents=True)
        (result_dir / "sample_R1_summary.txt").write_text(
            "=== sample_R1 Summary ===\n"
            f"Processed reads : {reads}\n"
            "Total bases : 180000\n"
            "GC content : 49.50%\n",
            encoding="utf-8",
        )
        (plot_dir / "quality_R1.svg").write_text(SVG, encoding="utf-8")
        (plot_dir / "plots_manifest.json").write_text(
            json.dumps(
                {
                    "schema_version": 1,
                    "plots": [
                        {
                            "id": "per_base_quality",
                            "read": "R1",
                            "title": "Per base sequence quality",
                            "status": "generated",
                            "svg": "quality_R1.svg",
                            "alt_text": "Quality by position",
                        }
                    ],
                }
            ),
            encoding="utf-8",
        )
        (result_dir / "qc_evaluation.json").write_text(
            json.dumps(
                {
                    "schema_version": 1,
                    "ruleset": {"id": "test-rules", "version": "1"},
                    "evaluations": [
                        {
                            "metric_id": "per_base_quality",
                            "read": "R1",
                            "qc_status": status,
                            "observations": {"minimum_median": 30.0},
                            "checks": [
                                {
                                    "observation": "minimum_median",
                                    "label": "Minimum median",
                                    "unit": "Phred",
                                }
                            ],
                            "reasons": [],
                        }
                    ],
                }
            ),
            encoding="utf-8",
        )
    (root / "trimming_report.json").write_text(
        json.dumps(
            {
                "reads": {"total": 1200, "passed": 1180, "discarded": 20},
                "bases": {"before": 180000, "after": 170000, "trimmed": 10000},
                "trimming": {
                    "quality_trimmed_reads": 30,
                    "adapter_trimmed_reads": 100,
                    "polyG_trimmed_reads": 2,
                    "polyX_trimmed_reads": 1,
                    "corrected_reads": 5,
                    "corrected_bases": 8,
                    "merged_pairs": 4,
                },
                "filtering": {"too_short_reads": 20},
                "output": {"reads": 1180, "discarded_pairs": 0},
                "adapter_positions": {"42": 71},
            }
        ),
        encoding="utf-8",
    )
    (root / "run_manifest.json").write_text(
        json.dumps(
            {
                "schema_version": 1,
                "trimming": {
                    "enabled": True,
                    "trim_front": 0,
                    "trim_tail": 0,
                    "cut_front": True,
                    "cut_tail": True,
                    "cut_right": True,
                    "quality_threshold": 20,
                    "window_size": 4,
                    "adapter_trimming": True,
                    "adapter_count": 4,
                    "min_adapter_match": 6,
                    "min_length": 36,
                    "trim_poly_g": True,
                    "trim_poly_x": True,
                    "poly_g_min_length": 10,
                    "poly_x_min_length": 10,
                    "overlap_correction": True,
                    "merge_reads": True,
                    "umi_enabled": False,
                    "umi_length": 0,
                },
            }
        ),
        encoding="utf-8",
    )
    trimmed_dir = root / "trimmed"
    trimmed_dir.mkdir()
    (trimmed_dir / "R1.trimmed.fastq.gz").touch()


class TrimmingReportTest(unittest.TestCase):
    def test_statistic_values_do_not_use_scientific_notation(self) -> None:
        self.assertEqual(_plain_value("1.234567e+06"), "1,234,567")
        self.assertEqual(_plain_value("-2.5e-04%"), "-0.0003%")

    def test_parameters_hide_inactive_optional_operations(self) -> None:
        parameters = _render_parameters(
            {
                "trimming": {
                    "overlap_correction": False,
                    "merge_reads": False,
                    "umi_enabled": False,
                    "umi_length": 0,
                }
            },
            ("R1", "R2"),
        )
        self.assertNotIn("PE overlap correction", parameters)
        self.assertNotIn("Merge overlapping pairs", parameters)
        self.assertNotIn("5′ UMI extraction", parameters)
        self.assertNotIn("UMI length (bp)", parameters)

    def test_report_compares_both_stages_and_is_self_contained(self) -> None:
        with tempfile.TemporaryDirectory(prefix="neoqc-trimming-report-") as temporary:
            root = Path(temporary)
            make_result(root)
            report = generate_trimming_report(root)
            self.assertEqual(report.name, REPORT_FILENAME)
            document = report.read_text(encoding="utf-8")
            self.assertIn("NeoQC trimming report — sample", document)
            self.assertIn("Before trimming", document)
            self.assertIn("After trimming", document)
            self.assertIn("Processed reads", document)
            self.assertIn("1,200", document)
            self.assertIn("1,180", document)
            self.assertIn('id="trim-accounting"', document)
            self.assertIn('id="trim-parameters"', document)
            self.assertIn('id="trim-compare"', document)
            self.assertIn("10,000", document)
            self.assertIn("trimmed/R1.trimmed.fastq.gz", document)
            self.assertIn('class="delta neutral">-20</span>', document)
            self.assertIn('id="before-per_base_quality-R1-en-', document)
            self.assertIn('id="after-per_base_quality-R1-en-', document)
            self.assertIn('<span class="l-ru">До тримминга</span>', document)
            self.assertNotIn('src="http', document)
            self.assertNotIn('href="http', document)

    def test_cli_writes_a_requested_output_path(self) -> None:
        with tempfile.TemporaryDirectory(prefix="neoqc-trimming-report-cli-") as temporary:
            root = Path(temporary)
            make_result(root)
            output = root / "custom.html"
            result = subprocess.run(
                [sys.executable, str(SCRIPTS / "generate_trimming_report.py"), str(root), "--output", str(output)],
                text=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                check=False,
            )
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertTrue(output.is_file())

    def test_root_plot_command_creates_the_combined_report(self) -> None:
        with tempfile.TemporaryDirectory(prefix="neoqc-trimming-report-plot-") as temporary:
            root = Path(temporary)
            make_result(root)
            with patch.object(plot_results, "_plot_one", return_value=0) as plot_one:
                result = plot_results.main([str(root), str(root / "plots")])
            self.assertEqual(result, 0)
            self.assertEqual(plot_one.call_count, 2)
            self.assertTrue((root / REPORT_FILENAME).is_file())


if __name__ == "__main__":
    unittest.main()
