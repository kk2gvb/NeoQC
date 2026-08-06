from __future__ import annotations

import json
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "python"))

from neo_mrna_vax_report import (
    ReportData,
    ReportValidationError,
    attach_qc_evaluation,
    render_html_report,
    write_html_report,
)
from neo_mrna_vax_report.cli import main


def valid_data() -> dict:
    return {
        "title": "Итоговый отчёт программы neo-mRNA-vax",
        "report_version": "0.1",
        "case_id": "CASE-001",
        "generated_at": "2026-08-03T12:00:00Z",
        "organization": "Клиника",
        "metadata": [{"label": "Референс", "value": "GRCh38"}],
        "sections": [
            {
                "id": "sequencing-quality-control",
                "title": "Контроль качества исходных данных",
                "status": "passed",
                "metrics": [
                    {
                        "label": "Доля Q30",
                        "value": 94.8,
                        "unit": "%",
                        "reference_range": ">= 85",
                        "status": "passed",
                    }
                ],
                "tables": [
                    {
                        "title": "Образцы",
                        "columns": ["Образец", "Роль"],
                        "rows": [["T-01", "Опухоль"]],
                    }
                ],
            }
        ],
        "conclusions": [{"text": "Материал пригоден.", "status": "passed"}],
        "disclaimer": "Для профессионального использования.",
    }


def qc_evaluation_data() -> dict:
    return {
        "schema_version": 1,
        "ruleset": {"id": "fastqc-compatible-v1", "version": "1.0.0"},
        "evaluations": [
            {
                "metric_id": "per_base_n_content",
                "title": "Per base N content",
                "read": "R1",
                "qc_status": "pass",
                "observations": {"maximum_n_percent": 0.2},
                "checks": [
                    {
                        "observation": "maximum_n_percent",
                        "label": "Maximum N content",
                        "unit": "%",
                    }
                ],
                "reasons": [
                    {"code": "n.pass", "message": "Within configured thresholds."}
                ],
            },
            {
                "metric_id": "adapter_content",
                "title": "Adapter content",
                "read": "R2",
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
                    {"code": "adapter.warning", "message": "Adapter content exceeds 5%."}
                ],
            },
        ],
    }


class HtmlReportTest(unittest.TestCase):
    def test_render_and_escape_untrusted_values(self) -> None:
        data = valid_data()
        data["metadata"].append(
            {"label": "Проверка", "value": "<script>alert('x')</script> & test"}
        )
        report = ReportData.from_dict(data)

        html = render_html_report(report, include_print_button=False)

        self.assertIn("<!doctype html>", html)
        self.assertIn("CASE-001", html)
        self.assertNotIn("<script>alert", html)
        self.assertIn("&lt;script&gt;alert(&#x27;x&#x27;)&lt;/script&gt; &amp; test", html)
        self.assertNotIn('onclick="window.print()"', html)

    def test_render_uses_approved_brand_palette(self) -> None:
        html = render_html_report(ReportData.from_dict(valid_data())).lower()

        for approved_color in ("#0a132d", "#2947a0", "#192f70", "#031e4f", "#539d96"):
            self.assertIn(approved_color, html)
        self.assertNotIn("#104a76", html)
        self.assertNotIn("#155e75", html)

    def test_invalid_table_is_rejected_with_field_path(self) -> None:
        data = valid_data()
        data["sections"][0]["tables"][0]["rows"] = [["only one cell"]]

        with self.assertRaisesRegex(ReportValidationError, r"rows\[0\].*expected 2"):
            ReportData.from_dict(data)

    def test_unknown_status_is_rejected(self) -> None:
        data = valid_data()
        data["sections"][0]["status"] = "maybe"

        with self.assertRaisesRegex(ReportValidationError, "must be one of"):
            ReportData.from_dict(data)

    def test_atomic_writer_replaces_existing_report(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            output = Path(temporary_directory) / "case" / "report.html"
            output.parent.mkdir()
            output.write_text("stale", encoding="utf-8")

            returned_path = write_html_report(ReportData.from_dict(valid_data()), output)

            self.assertEqual(returned_path, output.resolve())
            self.assertIn("CASE-001", output.read_text(encoding="utf-8"))
            self.assertEqual(list(output.parent.glob("*.tmp")), [])

    def test_cli_generates_report_from_json(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            source = root / "report.json"
            output = root / "report.html"
            source.write_text(json.dumps(valid_data(), ensure_ascii=False), encoding="utf-8")

            exit_code = main(
                ["--input", str(source), "--output", str(output), "--no-print-button"]
            )

            self.assertEqual(exit_code, 0)
            self.assertTrue(output.is_file())
            self.assertNotIn('onclick="window.print()"', output.read_text(encoding="utf-8"))

    def test_qc_evaluation_replaces_full_report_quality_section(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            evaluation = root / "qc_evaluation.json"
            evaluation.write_text(
                json.dumps(qc_evaluation_data()), encoding="utf-8"
            )
            report = attach_qc_evaluation(ReportData.from_dict(valid_data()), evaluation)

            section = next(
                item for item in report.sections if item.section_id == "sequencing-quality-control"
            )
            self.assertEqual(section.status.value, "warning")
            self.assertEqual([metric.value for metric in section.metrics], ["1", "1", "0", "0"])
            html = render_html_report(report)
            self.assertIn("Матрица проверок NeoQC", html)
            self.assertIn("fastqc-compatible-v1 v1.0.0", html)
            self.assertIn("Adapter content exceeds 5%.", html)
            self.assertIn("Maximum adapter content: 7 %", html)

    def test_cli_accepts_qc_evaluation(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            source = root / "report.json"
            evaluation = root / "qc_evaluation.json"
            output = root / "report.html"
            source.write_text(json.dumps(valid_data()), encoding="utf-8")
            evaluation.write_text(json.dumps(qc_evaluation_data()), encoding="utf-8")

            exit_code = main(
                [
                    "--input",
                    str(source),
                    "--qc-evaluation",
                    str(evaluation),
                    "--output",
                    str(output),
                ]
            )
            self.assertEqual(exit_code, 0)
            self.assertIn("Матрица проверок NeoQC", output.read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
