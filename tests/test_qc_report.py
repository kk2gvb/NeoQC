#!/usr/bin/env python3
"""Tests for the standalone, self-contained NeoQC QC report."""

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

from qc_report import QcReportError, generate_qc_report, load_report_model


def write_summary(path: Path, sample: str, read: str, injected: str = "") -> None:
    path.write_text(
        f"\n=== {sample}_{read} Summary ===\n"
        "Processed reads : 1200\n"
        "Total bases     : 180000\n"
        "Min length      : 150\n"
        "Max length      : 150\n"
        "Avg length      : 150.00\n\n"
        f"GC content      : 49.50%{injected}\n"
        "%N              : 0.01%\n",
        encoding="utf-8",
    )


def write_manifest(plot_dir: Path, *, unsafe: bool = False) -> None:
    plot_dir.mkdir(parents=True, exist_ok=True)
    (plot_dir / "quality_R1.svg").write_text(
        '<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 20 10"><path d="M0 5h20"/></svg>',
        encoding="utf-8",
    )
    manifest = {
        "schema_version": 1,
        "theme": "neo-report",
        "plots": [
            {
                "id": "per_base_quality",
                "read": "R1",
                "title": "Per base sequence quality",
                "source": "per_cycle_R1.tsv",
                "status": "generated",
                "svg": "../../outside.svg" if unsafe else "quality_R1.svg",
                "alt_text": "Mean quality for R1",
            },
            {
                "id": "adapter_content",
                "read": "R1",
                "title": "Adapter content",
                "source": "adapter_content_R1.tsv",
                "status": "skipped",
                "reason": "adapter_analysis_disabled",
            },
            {
                "id": "per_base_quality",
                "read": "R2",
                "title": "Per base sequence quality",
                "source": "per_cycle_R2.tsv",
                "status": "skipped",
                "reason": "source_not_found",
            },
        ],
    }
    (plot_dir / "plots_manifest.json").write_text(
        json.dumps(manifest), encoding="utf-8"
    )


class QcReportTest(unittest.TestCase):
    def test_self_contained_single_read_report(self) -> None:
        with tempfile.TemporaryDirectory(prefix="neoqc-report-") as temporary:
            result_dir = Path(temporary) / "sample result"
            plot_dir = result_dir / "plots"
            result_dir.mkdir(parents=True)
            write_summary(result_dir / "sample_A_R1_summary.txt", "sample_A", "R1")
            write_manifest(plot_dir)

            output = generate_qc_report(result_dir)
            document = output.read_text(encoding="utf-8")
            self.assertIn("<!doctype html>", document)
            self.assertIn("sample_A", document)
            self.assertIn("Processed reads", document)
            self.assertIn("data:image/svg+xml;base64,", document)
            self.assertNotIn("quality_R1.svg\"", document)
            self.assertIn("Adapter analysis was disabled", document)
            self.assertNotIn("source_not_found", document)
            self.assertNotIn("MODULE 00", document)
            self.assertNotRegex(document, r"MODULE [0-9]{2}")
            self.assertRegex(
                document,
                r"Generated</small><strong>\d{2}\.\d{2}\.\d{4}, \d{2}:\d{2} МСК",
            )
            self.assertIn('class="nav-index">00</span>', document)
            self.assertNotIn("IntersectionObserver", document)
            self.assertIn("window.scrollTo", document)
            self.assertIn("requestAnimationFrame", document)
            model = load_report_model(result_dir)
            self.assertEqual(model.reads, ("R1",))
            self.assertEqual(model.generated_plots, 1)

    def test_summary_values_are_html_escaped(self) -> None:
        with tempfile.TemporaryDirectory(prefix="neoqc-report-escape-") as temporary:
            result_dir = Path(temporary)
            write_summary(
                result_dir / "safe_R1_summary.txt",
                "safe",
                "R1",
                injected="<script>alert(1)</script>",
            )
            write_manifest(result_dir / "plots")
            document = generate_qc_report(result_dir).read_text(encoding="utf-8")
            self.assertNotIn("<script>alert(1)</script>", document)
            self.assertIn("&lt;script&gt;alert(1)&lt;/script&gt;", document)

    def test_asset_path_cannot_escape_plot_directory(self) -> None:
        with tempfile.TemporaryDirectory(prefix="neoqc-report-safe-path-") as temporary:
            result_dir = Path(temporary)
            write_summary(result_dir / "safe_R1_summary.txt", "safe", "R1")
            write_manifest(result_dir / "plots", unsafe=True)
            document = generate_qc_report(result_dir).read_text(encoding="utf-8")
            self.assertIn("asset path escapes the plot directory", document)
            self.assertNotIn("data:image/svg+xml;base64,", document)

    def test_standalone_cli(self) -> None:
        with tempfile.TemporaryDirectory(prefix="neoqc-report-cli-") as temporary:
            result_dir = Path(temporary)
            write_summary(result_dir / "cli_R1_summary.txt", "cli", "R1")
            write_manifest(result_dir / "plots")
            output = result_dir / "custom report.html"
            result = subprocess.run(
                [
                    sys.executable,
                    str(SCRIPTS / "generate_qc_report.py"),
                    str(result_dir),
                    "--output",
                    str(output),
                ],
                text=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                check=False,
            )
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertTrue(output.is_file())


if __name__ == "__main__":
    unittest.main()
