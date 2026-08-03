"""Safe, dependency-free rendering of standalone clinical HTML reports."""

from __future__ import annotations

import os
import tempfile
from html import escape
from pathlib import Path

from .models import ReportData, Status


_STATUS_LABELS = {
    Status.NOT_EVALUATED: "Не оценено",
    Status.PASSED: "Пройдено",
    Status.WARNING: "Требует внимания",
    Status.FAILED: "Не пройдено",
}

_CSS = """
:root{--ink:#18212b;--muted:#637080;--line:#dbe2e8;--panel:#f7f9fb;--brand:#155e75;--ok:#166534;--ok-bg:#dcfce7;--warn:#92400e;--warn-bg:#fef3c7;--fail:#991b1b;--fail-bg:#fee2e2;--neutral:#475569;--neutral-bg:#e2e8f0}
*{box-sizing:border-box}body{margin:0;background:#eef2f5;color:var(--ink);font:15px/1.55 -apple-system,BlinkMacSystemFont,"Segoe UI",Arial,sans-serif}.page{max-width:1080px;margin:32px auto;background:#fff;box-shadow:0 8px 28px #2634421a}.header{padding:40px 48px 32px;border-top:8px solid var(--brand);border-bottom:1px solid var(--line)}h1{font-size:30px;line-height:1.2;margin:0 0 10px}.subtitle{color:var(--muted);margin:0}.meta{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:0 32px;margin-top:28px}.meta-item{display:flex;justify-content:space-between;gap:20px;padding:9px 0;border-bottom:1px solid var(--line)}.meta-label{color:var(--muted)}main{padding:12px 48px 40px}.section{padding:28px 0;border-bottom:1px solid var(--line)}.section:last-child{border-bottom:0}.section-head{display:flex;align-items:center;justify-content:space-between;gap:20px;margin-bottom:12px}h2{font-size:21px;margin:0}h3{font-size:16px;margin:22px 0 8px}.summary{color:#334155;max-width:84ch;white-space:pre-line}.status{display:inline-block;border-radius:999px;padding:4px 10px;font-size:12px;font-weight:700;white-space:nowrap}.passed{color:var(--ok);background:var(--ok-bg)}.warning{color:var(--warn);background:var(--warn-bg)}.failed{color:var(--fail);background:var(--fail-bg)}.not_evaluated{color:var(--neutral);background:var(--neutral-bg)}.metrics{display:grid;grid-template-columns:repeat(3,minmax(0,1fr));gap:12px;margin-top:18px}.metric{padding:16px;border:1px solid var(--line);border-radius:8px;background:var(--panel)}.metric-label{color:var(--muted);font-size:13px}.metric-value{font-size:22px;font-weight:700;margin:5px 0}.unit{font-size:13px;font-weight:400;color:var(--muted)}.reference{font-size:12px;color:var(--muted);min-height:19px}.metric .status{margin-top:9px}table{width:100%;border-collapse:collapse;font-size:14px}th,td{text-align:left;vertical-align:top;padding:10px 12px;border:1px solid var(--line)}th{background:var(--panel)}.conclusions{margin:14px 0 0;padding:0;list-style:none}.conclusions li{display:flex;gap:12px;align-items:flex-start;padding:12px 0;border-bottom:1px solid var(--line)}.conclusions li:last-child{border:0}.disclaimer{padding:22px 48px;background:var(--panel);border-top:1px solid var(--line);color:var(--muted);font-size:12px;white-space:pre-line}.print-button{position:fixed;right:24px;bottom:24px;border:0;border-radius:8px;padding:11px 16px;background:var(--brand);color:#fff;font-weight:700;cursor:pointer}@media(max-width:700px){.page{margin:0}.header,main{padding-left:22px;padding-right:22px}.meta,.metrics{grid-template-columns:1fr}.disclaimer{padding:20px 22px}.table-wrap{overflow-x:auto}}@media print{body{background:#fff}.page{margin:0;max-width:none;box-shadow:none}.print-button{display:none}.section,.metric,table{break-inside:avoid}}
""".strip()


def _status(status: Status) -> str:
    return (
        f'<span class="status {status.value}">'
        f"{_STATUS_LABELS[status]}</span>"
    )


def render_html_report(
    report: ReportData, *, language: str = "ru", include_print_button: bool = True
) -> str:
    """Return a standalone UTF-8 HTML document for validated report data."""
    report.validate()
    parts = [
        "<!doctype html>",
        f'<html lang="{escape(language, quote=True)}">',
        "<head>",
        '<meta charset="utf-8">',
        '<meta name="viewport" content="width=device-width, initial-scale=1">',
        f"<title>{escape(report.title)} — {escape(report.case_id)}</title>",
        f"<style>{_CSS}</style>",
        "</head>",
        "<body>",
        '<article class="page">',
        '<header class="header">',
        f"<h1>{escape(report.title)}</h1>",
        '<p class="subtitle">Версия отчёта '
        f"{escape(report.report_version)}"
        + (f" · {escape(report.organization)}" if report.organization else "")
        + "</p>",
        '<div class="meta">',
        '<div class="meta-item"><span class="meta-label">Идентификатор случая</span>'
        f"<strong>{escape(report.case_id)}</strong></div>",
        '<div class="meta-item"><span class="meta-label">Дата формирования</span>'
        f"<strong>{escape(report.generated_at)}</strong></div>",
    ]
    for item in report.metadata:
        parts.append(
            '<div class="meta-item"><span class="meta-label">'
            f"{escape(item.label)}</span><strong>{escape(item.value)}</strong></div>"
        )
    parts.extend(["</div>", "</header>", "<main>"])

    for section in report.sections:
        parts.extend(
            [
                f'<section class="section" id="{escape(section.section_id, quote=True)}">',
                f'<div class="section-head"><h2>{escape(section.title)}</h2>'
                f"{_status(section.status)}</div>",
            ]
        )
        if section.summary:
            parts.append(f'<p class="summary">{escape(section.summary)}</p>')
        if section.metrics:
            parts.append('<div class="metrics">')
            for metric in section.metrics:
                unit = (
                    f' <span class="unit">{escape(metric.unit)}</span>'
                    if metric.unit
                    else ""
                )
                reference = (
                    f"Референс: {escape(metric.reference_range)}"
                    if metric.reference_range
                    else ""
                )
                parts.append(
                    '<div class="metric"><div class="metric-label">'
                    f'{escape(metric.label)}</div><div class="metric-value">'
                    f'{escape(metric.value)}{unit}</div><div class="reference">'
                    f"{reference}</div>{_status(metric.status)}</div>"
                )
            parts.append("</div>")
        for table in section.tables:
            if table.title:
                parts.append(f"<h3>{escape(table.title)}</h3>")
            headers = "".join(
                f'<th scope="col">{escape(column)}</th>' for column in table.columns
            )
            rows = "".join(
                "<tr>" + "".join(f"<td>{escape(cell)}</td>" for cell in row) + "</tr>"
                for row in table.rows
            )
            parts.append(
                '<div class="table-wrap"><table><thead><tr>'
                f"{headers}</tr></thead><tbody>{rows}</tbody></table></div>"
            )
        parts.append("</section>")

    if report.conclusions:
        parts.append(
            '<section class="section"><div class="section-head">'
            "<h2>Итоговые выводы</h2></div><ul class=\"conclusions\">"
        )
        parts.extend(
            f"<li>{_status(item.status)}<span>{escape(item.text)}</span></li>"
            for item in report.conclusions
        )
        parts.append("</ul></section>")
    parts.extend(["</main>"])
    if report.disclaimer:
        parts.append(f'<footer class="disclaimer">{escape(report.disclaimer)}</footer>')
    parts.append("</article>")
    if include_print_button:
        parts.append(
            '<button class="print-button" type="button" '
            'onclick="window.print()">Печать / PDF</button>'
        )
    parts.extend(["</body>", "</html>", ""])
    return "\n".join(parts)


def write_html_report(
    report: ReportData,
    output_path: str | os.PathLike[str],
    *,
    language: str = "ru",
    include_print_button: bool = True,
) -> Path:
    """Atomically write a standalone report and return its resolved path."""
    destination = Path(output_path).resolve()
    destination.parent.mkdir(parents=True, exist_ok=True)
    html = render_html_report(
        report, language=language, include_print_button=include_print_button
    )

    temporary_path: str | None = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w",
            encoding="utf-8",
            dir=destination.parent,
            prefix=f".{destination.name}.",
            suffix=".tmp",
            delete=False,
        ) as temporary:
            temporary_path = temporary.name
            temporary.write(html)
            temporary.flush()
            os.fsync(temporary.fileno())
        os.replace(temporary_path, destination)
    except Exception:
        if temporary_path is not None:
            Path(temporary_path).unlink(missing_ok=True)
        raise
    return destination
