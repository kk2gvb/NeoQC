from __future__ import annotations

import json
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "python"))

from neoqc_report import ReportData, ReportValidationError, render_html_report, write_html_report
from neoqc_report.cli import main


def valid_data() -> dict:
    return {
        "title": "Отчёт по персонализированной мРНК-вакцине",
        "report_version": "0.1",
        "case_id": "CASE-001",
        "generated_at": "2026-08-03T12:00:00Z",
        "organization": "Клиника",
        "metadata": [{"label": "Референс", "value": "GRCh38"}],
        "sections": [
            {
                "id": "sequencing-qc",
                "title": "Контроль качества",
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


if __name__ == "__main__":
    unittest.main()
