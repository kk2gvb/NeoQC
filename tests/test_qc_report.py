#!/usr/bin/env python3
"""Tests for the standalone, self-contained NeoQC QC report."""

from __future__ import annotations

import json
import re
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "scripts"
sys.path.insert(0, str(SCRIPTS))

from neoqc_i18n import chart_ru, ru
from neoqc_theme import CHART_COLORS, LIGHT_TO_VAR
from qc_report import QcReportError, generate_qc_report, load_report_model, render_qc_report


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


MALICIOUS_SVG = """<?xml version="1.0" encoding="utf-8" standalone="no"?>
<!DOCTYPE svg PUBLIC "-//W3C//DTD SVG 1.1//EN" "http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd">
<svg xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink"
     width="576pt" height="324pt" viewBox="0 0 576 324" onload="alert(1)">
 <style type="text/css">*{stroke-linecap: butt}</style>
 <script>alert(2)</script>
 <foreignObject><div xmlns="http://www.w3.org/1999/xhtml">x</div></foreignObject>
 <defs><clipPath id="p1"><rect width="10" height="10"/></clipPath>
   <path id="m1" d="M0 0h1" style="stroke: #1f5fd1"/></defs>
 <g clip-path="url(#p1)">
  <path d="M0 5h20" style="fill: #ffffff; stroke: #d8343d; clip-path: url(#p1)" onclick="alert(3)"/>
  <use xlink:href="#m1" x="1" y="1" fill="#12a150"/>
  <use xlink:href="https://example.org/evil.svg#x"/>
  <path d="M0 0" style="fill: url(https://example.org/x.svg#y)"/>
  <text style="font-size: 10px; font-family: 'IBM Plex Sans'; fill: #1b2430">Mean</text>
 </g>
</svg>"""


def set_plot_svg(result_dir: Path, filename: str, content: str, *, localized: dict | None = None) -> None:
    plot_dir = result_dir / "plots"
    (plot_dir / filename).write_text(content, encoding="utf-8")
    manifest_path = plot_dir / "plots_manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    entry = manifest["plots"][0]
    entry["svg"] = filename
    if localized is not None:
        entry["localized"] = localized
    manifest_path.write_text(json.dumps(manifest), encoding="utf-8")


def make_report(temporary: str) -> Path:
    result_dir = Path(temporary)
    write_summary(result_dir / "sample_R1_summary.txt", "sample", "R1")
    write_manifest(result_dir / "plots")
    return result_dir


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
            self.assertIn('class="chart-svg"', document)
            self.assertNotIn("data:image/svg+xml;base64,", document)
            self.assertNotIn("quality_R1.svg\"", document)
            self.assertIn("Adapter analysis was disabled", document)
            self.assertNotIn("source_not_found", document)
            self.assertNotIn("MODULE 00", document)
            self.assertNotRegex(document, r"MODULE [0-9]{2}")
            self.assertRegex(
                document,
                r'Generated</span><span class="l-ru">Создан</span></small>'
                r"<strong>\d{2}\.\d{2}\.\d{4}, \d{2}:\d{2} GMT\+3</strong>",
            )
            self.assertIn('<span class="i">00</span>', document)
            self.assertNotIn("IntersectionObserver", document)
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
            self.assertNotIn('class="chart-svg"', document)

    def test_overrepresented_sequences_table_is_embedded_and_escaped(self) -> None:
        with tempfile.TemporaryDirectory(prefix="neoqc-report-overrepresented-") as temporary:
            result_dir = Path(temporary)
            write_summary(result_dir / "sample_R1_summary.txt", "sample", "R1")
            write_manifest(result_dir / "plots")
            add_duplication_plot(result_dir)
            sequence = "T" * 50
            (result_dir / "overrepresented_sequences_R1.tsv").write_text(
                "sequence\tcount\tpercentage\n"
                f"{sequence}\t269055\t0.4065003125\n",
                encoding="utf-8",
            )

            document = generate_qc_report(result_dir).read_text(encoding="utf-8")
            self.assertIn("Overrepresented sequences", document)
            self.assertIn(sequence, document)
            self.assertIn("269,055", document)
            self.assertIn("0.4065%", document)

    def test_empty_overrepresented_sequences_table_has_explicit_state(self) -> None:
        with tempfile.TemporaryDirectory(prefix="neoqc-report-overrepresented-empty-") as temporary:
            result_dir = Path(temporary)
            write_summary(result_dir / "sample_R1_summary.txt", "sample", "R1")
            write_manifest(result_dir / "plots")
            add_duplication_plot(result_dir)
            (result_dir / "overrepresented_sequences_R1.tsv").write_text(
                "sequence\tcount\tpercentage\n", encoding="utf-8"
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
                "sequence\tcount\tpercentage\n"
                "ACGT\t12\t101\n",
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

    def test_inline_svg_is_sanitized_themed_and_isolated(self) -> None:
        with tempfile.TemporaryDirectory(prefix="neoqc-report-svg-") as temporary:
            result_dir = make_report(temporary)
            set_plot_svg(result_dir, "quality_R1.svg", MALICIOUS_SVG)
            document = generate_qc_report(result_dir).read_text(encoding="utf-8")
            body = document[document.index("<body"):]
            chart = re.search(r'<svg[^>]*class="chart-svg".*?</svg>', body, re.S).group(0)
            for forbidden in ("<script", "onload", "onclick", "foreignObject", "<style", "example.org", "alert("):
                self.assertNotIn(forbidden, chart)
            # Theme colours become CSS variables, including presentation attributes.
            self.assertIn("var(--c-bg)", chart)
            self.assertIn("var(--c-danger)", chart)
            self.assertIn("fill: var(--c-a)", chart)
            self.assertNotRegex(chart, r"#[0-9a-fA-F]{6}\b")
            # Ids and references are prefixed so several charts can share one page.
            self.assertIn('id="per_base_quality-R1-en-p1"', chart)
            self.assertIn("url(#per_base_quality-R1-en-p1)", chart)
            self.assertIn('href="#per_base_quality-R1-en-m1"', chart)
            self.assertNotIn(' width="576pt"', chart)
            self.assertIn("'IBM Plex Sans', 'Segoe UI', Arial, sans-serif", chart)

    def test_svg_with_entity_declarations_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory(prefix="neoqc-report-entity-") as temporary:
            result_dir = make_report(temporary)
            set_plot_svg(
                result_dir,
                "quality_R1.svg",
                '<?xml version="1.0"?><!DOCTYPE svg [<!ENTITY a "aaaa">]>'
                '<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 1 1"><text>&a;</text></svg>',
            )
            document = generate_qc_report(result_dir).read_text(encoding="utf-8")
            self.assertIn("SVG with entity declarations is not accepted", document)
            self.assertNotIn('class="chart-svg"', document)

    def test_localized_chart_variants_and_fallback(self) -> None:
        svg = '<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 20 10"><text>{}</text></svg>'
        with tempfile.TemporaryDirectory(prefix="neoqc-report-l10n-") as temporary:
            result_dir = make_report(temporary)
            (result_dir / "plots" / "quality_R1.ru.svg").write_text(svg.format("Среднее"), encoding="utf-8")
            set_plot_svg(result_dir, "quality_R1.svg", svg.format("Mean"),
                         localized={"ru": {"svg": "quality_R1.ru.svg"}})
            document = generate_qc_report(result_dir).read_text(encoding="utf-8")
            self.assertIn('<span class="chart l-en">', document)
            self.assertIn('<span class="chart l-ru">', document)
            self.assertIn("Среднее", document)

            # A missing translation keeps the English chart visible in both languages.
            (result_dir / "plots" / "quality_R1.ru.svg").unlink()
            document = generate_qc_report(result_dir).read_text(encoding="utf-8")
            self.assertIn('<span class="chart"><svg', document)
            self.assertNotIn('<span class="chart l-ru">', document)

            set_plot_svg(result_dir, "quality_R1.svg", svg.format("Mean"), localized={"de": {"svg": "x.svg"}})
            with self.assertRaisesRegex(QcReportError, "unsupported locale"):
                load_report_model(result_dir)

    def test_theme_language_controls_and_print_are_self_contained(self) -> None:
        with tempfile.TemporaryDirectory(prefix="neoqc-report-ui-") as temporary:
            document = generate_qc_report(make_report(temporary)).read_text(encoding="utf-8")
            head = document[: document.index("</head>")]
            # Saved theme/language are applied in <head>, before the first paint.
            self.assertIn("neoqc.report.theme", head)
            self.assertIn("prefers-color-scheme: dark", head)
            for control in ('data-theme-btn="light"', 'data-theme-btn="dark"',
                            'data-lang-btn="en"', 'data-lang-btn="ru"',
                            'data-order-btn="triage"', 'data-order-btn="sequential"', 'id="fold-all"'):
                self.assertIn(control, document)
            self.assertIn('<span class="l-ru">Основная статистика</span>', document)
            self.assertIn(':root[data-theme="dark"]', document)
            # Printing always uses the light palette.
            self.assertRegex(document, r'@media print\{:root,:root\[data-theme="dark"\]\{color-scheme:light')
            for name, (light, dark) in CHART_COLORS.items():
                self.assertIn(f"{name}:{light}", document)
                self.assertIn(f"{name}:{dark}", document)
            # No external resources: the file must work offline.
            self.assertNotRegex(document, r'(src|href)="https?://')
            self.assertNotIn("@import", document)

    def test_embedded_fonts_cover_latin_and_cyrillic(self) -> None:
        with tempfile.TemporaryDirectory(prefix="neoqc-report-fonts-") as temporary:
            document = generate_qc_report(make_report(temporary)).read_text(encoding="utf-8")
            faces = re.findall(r"@font-face\{[^}]*\}", document)
            self.assertEqual(len(faces), 14)  # 7 weights x (latin + cyrillic)
            self.assertTrue(any("font-family:'IBM Plex Sans'" in face and "U+0400-045F" in face for face in faces))
            self.assertTrue(any("font-family:'IBM Plex Mono'" in face and "U+0000-00FF" in face for face in faces))

    def test_triage_order_keys_put_problems_first(self) -> None:
        with tempfile.TemporaryDirectory(prefix="neoqc-report-order-") as temporary:
            result_dir = make_report(temporary)
            write_evaluation(result_dir)
            document = generate_qc_report(result_dir).read_text(encoding="utf-8")

            def section_keys(section_id: str) -> tuple[int, int, str]:
                match = re.search(
                    rf'<section class="panel[^"]*" id="{section_id}"[^>]*?(?:data-sev="([a-z_]+)" )?'
                    r'data-seq="(\d+)" data-tri="(\d+)"',
                    document,
                )
                self.assertIsNotNone(match, section_id)
                return int(match.group(2)), int(match.group(3)), match.group(1) or ""

            overview = section_keys("qc-summary")
            basic = section_keys("basic-statistics")
            quality = section_keys("module-per_base_quality")   # PASS
            adapter = section_keys("module-adapter_content")    # WARNING
            self.assertEqual(overview[:2], (0, 0))
            self.assertEqual(basic[:2], (1, 1))                  # INFO always follows the overview
            self.assertLess(quality[0], adapter[0])              # pipeline order
            self.assertLess(adapter[1], quality[1])              # triage: problems first
            self.assertEqual(adapter[2], "warning")
            self.assertEqual(quality[2], "ok")                   # folded in triage mode
            self.assertIn("Needs attention", document)

    def test_report_javascript_is_valid(self) -> None:
        node = shutil.which("node")
        if node is None:
            self.skipTest("node is not installed")
        with tempfile.TemporaryDirectory(prefix="neoqc-report-js-") as temporary:
            document = generate_qc_report(make_report(temporary)).read_text(encoding="utf-8")
            scripts = re.findall(r"<script>(.*?)</script>", document, re.S)
            self.assertEqual(len(scripts), 2)
            for index, script in enumerate(scripts):
                path = Path(temporary) / f"script{index}.js"
                path.write_text(script, encoding="utf-8")
                result = subprocess.run([node, "--check", str(path)], capture_output=True, text=True, check=False)
                self.assertEqual(result.returncode, 0, result.stderr)


class TranslationTest(unittest.TestCase):
    def test_reason_messages_are_translated_with_values(self) -> None:
        self.assertEqual(
            ru("Maximum adapter content: 11.95 %; FAIL threshold > 10 %."),
            "Макс. доля адаптеров: 11.95 %; порог FAIL > 10 %.",
        )
        self.assertEqual(ru("All evaluated observations are within configured thresholds."),
                         "Все наблюдения в пределах заданных порогов.")
        self.assertEqual(ru("Custom ruleset message"), "Custom ruleset message")

    def test_chart_labels_and_patterns(self) -> None:
        self.assertEqual(chart_ru("Position in read (bp)"), "Позиция в риде (п.н.)")
        self.assertEqual(chart_ru("Mean GC: 47.9"), "Среднее GC: 47.9")
        self.assertEqual(chart_ru("Peak 0.26%"), "Пик 0.26%")
        self.assertEqual(chart_ru("1.0M reads"), "1.0M ридов")
        self.assertEqual(chart_ru("TruSeq_R1"), "TruSeq_R1")

    def test_theme_tokens_map_back_to_css_variables(self) -> None:
        self.assertEqual(len(LIGHT_TO_VAR), len(CHART_COLORS))
        for name, (light, dark) in CHART_COLORS.items():
            self.assertEqual(LIGHT_TO_VAR[light.lower()], name)
            self.assertRegex(dark, r"^#[0-9a-f]{6}$")

    def test_render_does_not_mutate_model(self) -> None:
        with tempfile.TemporaryDirectory(prefix="neoqc-report-pure-") as temporary:
            result_dir = make_report(temporary)
            model = load_report_model(result_dir)
            first = render_qc_report(model, result_dir / "plots")
            second = render_qc_report(model, result_dir / "plots")
            self.assertEqual(first, second)


if __name__ == "__main__":
    unittest.main()
