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


def add_duplication_plot(result_dir: Path) -> None:
    plot_dir = result_dir / "plots"
    (plot_dir / "duplication_R1.svg").write_text(
        '<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 20 10"></svg>',
        encoding="utf-8",
    )
    manifest_path = plot_dir / "plots_manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest["plots"].append(
        {
            "id": "sequence_duplication_levels",
            "read": "R1",
            "title": "Sequence duplication levels",
            "source": "sequence_duplication_levels_R1.tsv",
            "status": "generated",
            "svg": "duplication_R1.svg",
            "alt_text": "Sequence duplication levels for R1",
        }
    )
    manifest_path.write_text(json.dumps(manifest), encoding="utf-8")


def write_evaluation(result_dir: Path) -> None:
    evaluation = {
        "schema_version": 1,
        "ruleset": {"id": "test-rules", "version": "1.2.3"},
        "evaluations": [
            {
                "metric_id": "per_base_quality",
                "read": "R1",
                "qc_status": "pass",
                "observations": {"minimum_median": 31.0},
                "checks": [
                    {
                        "observation": "minimum_median",
                        "label": "Minimum median",
                        "unit": "Phred",
                    }
                ],
                "reasons": [
                    {
                        "code": "quality.within_thresholds",
                        "message": "All observations are within thresholds.",
                    }
                ],
            },
            {
                "metric_id": "adapter_content",
                "read": "R1",
                "qc_status": "warning",
                "observations": {"maximum_adapter_percent": 7.0},
                "checks": [
                    {
                        "observation": "maximum_adapter_percent",
                        "label": "Maximum adapter content",
                        "unit": "%",
                    }
                ],
                "reasons": [
                    {
                        "code": "adapter.warning",
                        "message": "Adapter content exceeds 5%.",
                    }
                ],
            },
        ],
    }
    (result_dir / "qc_evaluation.json").write_text(
        json.dumps(evaluation), encoding="utf-8"
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

    def test_qc_distribution_matrix_and_card_decisions(self) -> None:
        with tempfile.TemporaryDirectory(prefix="neoqc-report-status-") as temporary:
            result_dir = Path(temporary)
            write_summary(result_dir / "status_R1_summary.txt", "status", "R1")
            write_manifest(result_dir / "plots")
            write_evaluation(result_dir)

            document = generate_qc_report(result_dir).read_text(encoding="utf-8")
            self.assertIn("PASS / WARNING / FAIL distribution", document)
            self.assertIn("test-rules · v1.2.3", document)
            self.assertIn("PASS 1, WARNING 1, FAIL 0, NOT EVALUATED 0", document)
            self.assertIn("✓</span> PASS", document)
            self.assertIn("▲</span> WARNING", document)
            self.assertIn("Minimum median", document)
            self.assertIn("31 Phred", document)
            self.assertIn("Adapter content exceeds 5%.", document)

    def test_asset_path_cannot_escape_plot_directory(self) -> None:
        with tempfile.TemporaryDirectory(prefix="neoqc-report-safe-path-") as temporary:
            result_dir = Path(temporary)
            write_summary(result_dir / "safe_R1_summary.txt", "safe", "R1")
            write_manifest(result_dir / "plots", unsafe=True)
            document = generate_qc_report(result_dir).read_text(encoding="utf-8")
            self.assertIn("asset path escapes the plot directory", document)
            self.assertNotIn("data:image/svg+xml;base64,", document)

    def test_overrepresented_sequences_table_is_embedded_and_escaped(self) -> None:
        with tempfile.TemporaryDirectory(prefix="neoqc-report-overrepresented-") as temporary:
            result_dir = Path(temporary)
            write_summary(result_dir / "sample_R1_summary.txt", "sample", "R1")
            write_manifest(result_dir / "plots")
            add_duplication_plot(result_dir)
            sequence = "T" * 50
            (result_dir / "overrepresented_sequences_R1.tsv").write_text(
                "sequence\tcount\tpercentage\tpossible_source\n"
                f"{sequence}\t269055\t0.4065003125\tNo Hit <script>alert(1)</script>\n",
                encoding="utf-8",
            )

            document = generate_qc_report(result_dir).read_text(encoding="utf-8")
            self.assertIn("Overrepresented sequences", document)
            self.assertIn(sequence, document)
            self.assertIn("269,055", document)
            self.assertIn("0.4065%", document)
            self.assertIn("No Hit &lt;script&gt;alert(1)&lt;/script&gt;", document)
            self.assertNotIn("No Hit <script>alert(1)</script>", document)

    def test_empty_overrepresented_sequences_table_has_explicit_state(self) -> None:
        with tempfile.TemporaryDirectory(prefix="neoqc-report-overrepresented-empty-") as temporary:
            result_dir = Path(temporary)
            write_summary(result_dir / "sample_R1_summary.txt", "sample", "R1")
            write_manifest(result_dir / "plots")
            add_duplication_plot(result_dir)
            (result_dir / "overrepresented_sequences_R1.tsv").write_text(
                "sequence\tcount\tpercentage\tpossible_source\n", encoding="utf-8"
            )

            document = generate_qc_report(result_dir).read_text(encoding="utf-8")
            self.assertIn("No sequences exceeded the reporting threshold.", document)
            self.assertIn('class="table-count">0</span>', document)

    def test_invalid_overrepresented_percentage_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory(prefix="neoqc-report-overrepresented-invalid-") as temporary:
            result_dir = Path(temporary)
            write_summary(result_dir / "sample_R1_summary.txt", "sample", "R1")
            write_manifest(result_dir / "plots")
            add_duplication_plot(result_dir)
            (result_dir / "overrepresented_sequences_R1.tsv").write_text(
                "sequence\tcount\tpercentage\tpossible_source\n"
                "ACGT\t12\t101\tNo Hit\n",
                encoding="utf-8",
            )

            with self.assertRaisesRegex(QcReportError, "percentage must be between 0 and 100"):
                load_report_model(result_dir)

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
