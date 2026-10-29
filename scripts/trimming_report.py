#!/usr/bin/env python3
"""Render a self-contained NeoQC before/after trimming HTML report."""

from __future__ import annotations

from dataclasses import replace
from html import escape
import json
from pathlib import Path
import re
from typing import Mapping

import qc_report as qc


REPORT_FILENAME = "neoqc_trimming_report.html"


def _by_read(model: qc.QcReportModel) -> dict[str, qc.BasicStatistics]:
    return {summary.read: summary for summary in model.basic_statistics}


def _cards(module: qc.PlotModule | None) -> dict[str, qc.PlotCard]:
    return {card.read: card for card in module.cards} if module else {}


def _modules(model: qc.QcReportModel) -> dict[str, qc.PlotModule]:
    return {module.metric_id: module for module in model.modules}


def _status(card: qc.PlotCard | None) -> str:
    return card.qc_status if card else "not_evaluated"


def _stage_chart(card: qc.PlotCard | None, plot_dir: Path, stage: str) -> str:
    if card is None:
        return f'<div class="empty">{qc._t("No result")}</div>'
    # The SVG sanitizer uses metric/read as an ID prefix.  Add the stage so
    # before and after charts cannot collide when embedded in one document.
    return qc._chart_markup(replace(card, metric_id=f"{stage}-{card.metric_id}"), plot_dir)


def _value(card: qc.PlotCard | None, label: str) -> str:
    if card is None:
        return "—"
    for observed_label, number, unit in card.qc_observations:
        if observed_label == label:
            return f"{number:.4g}" + (f" {unit}" if unit else "")
    return "—"


_NUMBER = re.compile(r"^\s*([+-]?(?:[0-9][0-9,]*(?:\.[0-9]+)?|\.[0-9]+)(?:[eE][+-]?[0-9]+)?)\s*(%?)\s*$")


def _plain_number(value: float, *, signed: bool = False) -> str:
    """Render numeric values for people, never using scientific notation."""
    sign = "+" if signed else ""
    if abs(value) >= 1000 or value.is_integer():
        return f"{value:{sign},.0f}"
    return f"{value:{sign},.4f}".rstrip("0").rstrip(".")


def _plain_value(value: str) -> str:
    match = _NUMBER.match(value)
    if not match:
        return value
    return _plain_number(float(match.group(1).replace(",", ""))) + match.group(2)


def _delta(before: str, after: str) -> str:
    """Format a human-readable difference when two displayed values are numeric."""
    before_match, after_match = _NUMBER.match(before), _NUMBER.match(after)
    if not before_match or not after_match:
        return "=" if before == after else "—"
    first = float(before_match.group(1).replace(",", ""))
    second = float(after_match.group(1).replace(",", ""))
    difference = second - first
    if abs(difference) < 1e-10:
        return "="
    suffix = " pp" if before_match.group(2) == after_match.group(2) == "%" else ""
    return _plain_number(difference, signed=True) + suffix


def _change(before: str, after: str) -> str:
    delta = _delta(before, after)
    css = "same" if delta == "=" else "neutral"
    return f'<span class="delta {css}">{escape(delta)}</span>'


def _load_json_object(path: Path) -> Mapping[str, object] | None:
    if not path.is_file():
        return None
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise qc.QcReportError(f"cannot read {path.name}: {error}") from error
    if not isinstance(value, Mapping):
        raise qc.QcReportError(f"{path.name} must contain a JSON object")
    return value


def _number(value: object) -> int:
    return int(value) if isinstance(value, int) and not isinstance(value, bool) else 0


def _nested(data: Mapping[str, object], key: str) -> Mapping[str, object]:
    value = data.get(key)
    return value if isinstance(value, Mapping) else {}


def _sample_id_from_run_manifest(manifest: Mapping[str, object] | None) -> str | None:
    """Return the authoritative sample ID recorded for the trimming run."""
    sample_id = _nested(manifest or {}, "run").get("sample_id")
    return sample_id.strip() if isinstance(sample_id, str) and sample_id.strip() else None


def _with_root_sample_id_for_missing_summaries(
    model: qc.QcReportModel, root_sample_id: str
) -> qc.QcReportModel:
    """Avoid treating ``before`` and ``after`` directory names as sample IDs.

    ``qc.load_report_model`` uses its input directory name if no summary file is
    available.  That is a suitable fallback for a standalone QC report, but it
    makes the two stages of one trimming run look like different samples.
    """
    return replace(model, sample_id=root_sample_id) if not model.basic_statistics else model


def _render_statistics(before: qc.QcReportModel, after: qc.QcReportModel) -> str:
    before_stats, after_stats = _by_read(before), _by_read(after)
    reads = tuple(read for read in ("R1", "R2") if read in before.reads or read in after.reads)
    labels: list[str] = []
    for stats in (*before_stats.values(), *after_stats.values()):
        for label, _ in stats.items:
            if label not in labels:
                labels.append(label)
    if not labels:
        body = f'<div class="empty">{qc._t("Summary files were not found.")}</div>'
    else:
        head = "".join(
            f'<th colspan="3" scope="colgroup" class="n">{escape(read)}</th>' for read in reads
        )
        stages = "".join(
            f'<th scope="col" class="n">{qc._t("Before")}</th><th scope="col" class="n">{qc._t("After")}</th><th scope="col" class="n">{qc._t("Change")}</th>'
            for _ in reads
        )
        values_before = {read: dict(stats.items) for read, stats in before_stats.items()}
        values_after = {read: dict(stats.items) for read, stats in after_stats.items()}
        rows = "".join(
            f'<tr><th scope="row" class="row-h">{qc._t(label)}</th>'
            + "".join(
                f'<td class="n">{escape(_plain_value(before_value := values_before.get(read, {}).get(label, "—")))}</td>'
                f'<td class="n">{escape(_plain_value(after_value := values_after.get(read, {}).get(label, "—")))}</td>'
                f'<td class="n">{_change(before_value, after_value)}</td>'
                for read in reads
            )
            + "</tr>"
            for label in labels
        )
        body = (
            '<div class="table-wrap"><table class="grid"><thead><tr><th rowspan="2" scope="col">'
            f'{qc._t("Metric")}</th>{head}</tr><tr>{stages}</tr></thead><tbody>{rows}</tbody></table></div>'
        )
    return (
        '<section class="panel foldable" id="trim-basic" data-section data-sev="info" '
        f'{qc._keys(3, 3)}><div class="phead" data-fold>{qc.CHEVRON}'
        f'<h2>{qc._t("Basic statistics: before / after")}</h2><span class="sp"></span>{qc._pill("info", "INFO")}'
        f'</div>{body}</section>'
    )


def _yes_no(value: object) -> str:
    return qc._t("Yes") if value is True else qc._t("No")


def _render_parameters(
    manifest: Mapping[str, object] | None,
    reads: tuple[str, ...],
    output_fastqs: Mapping[str, Path] | None = None,
) -> str:
    trim = _nested(manifest or {}, "trimming")
    if not trim:
        body = f'<div class="empty">{qc._t("Trimming parameters were not recorded for this run.")}</div>'
    else:
        entries = [
            ("Library layout", qc._t("Paired") if len(reads) == 2 else qc._t("Single")),
            ("Fixed trim at 5′ end (bp)", str(_number(trim.get("trim_front")))),
            ("Fixed trim at 3′ end (bp)", str(_number(trim.get("trim_tail")))),
            ("Quality trim from 5′ end", _yes_no(trim.get("cut_front"))),
            ("Quality trim from 3′ end", _yes_no(trim.get("cut_tail"))),
            ("Sliding-window quality trim", _yes_no(trim.get("cut_right"))),
            ("Quality threshold (Phred)", str(_number(trim.get("quality_threshold")))),
            ("Window size (bp)", str(_number(trim.get("window_size")))),
            ("Adapter trimming", _yes_no(trim.get("adapter_trimming"))),
            ("Adapter sequences configured", str(_number(trim.get("adapter_count")))),
            ("Minimum adapter match (bp)", str(_number(trim.get("min_adapter_match")))),
            ("Minimum read length (bp)", str(_number(trim.get("min_length")))),
            ("Terminal polyG trimming", _yes_no(trim.get("trim_poly_g"))),
            ("Terminal polyX trimming", _yes_no(trim.get("trim_poly_x"))),
            ("Minimum polyG length (bp)", str(_number(trim.get("poly_g_min_length")))),
            ("Minimum polyX length (bp)", str(_number(trim.get("poly_x_min_length")))),
        ]
        # These operations are optional and should be listed only when they
        # actually ran; displaying a row with "No" makes the report look as
        # though they are part of the selected trimming workflow.
        if trim.get("overlap_correction") is True:
            entries.append(("PE overlap correction", _yes_no(True)))
        if trim.get("merge_reads") is True:
            entries.append(("Merge overlapping pairs", _yes_no(True)))
        if trim.get("umi_enabled") is True:
            entries.extend((
                ("5′ UMI extraction", _yes_no(True)),
                ("UMI length (bp)", str(_number(trim.get("umi_length")))),
            ))
        for read in ("R1", "R2", "merged"):
            path = (output_fastqs or {}).get(read)
            if path is not None:
                label = "Merged FASTQ.gz" if read == "merged" else f"Trimmed {read} FASTQ.gz"
                entries.append((label, f'<code class="output-path">{escape(str(path))}</code>'))
        rows = "".join(
            f'<tr><th scope="row" class="row-h">{qc._t(label)}</th><td class="v">{value}</td></tr>'
            for label, value in entries
        )
        body = f'<div class="table-wrap"><table class="grid params"><tbody>{rows}</tbody></table></div>'
    return (
        '<section class="panel foldable" id="trim-parameters" data-section data-sev="info" '
        f'{qc._keys(1, 1)}><div class="phead" data-fold>{qc.CHEVRON}<h2>{qc._t("Trimming parameters")}</h2>'
        f'<span class="sp"></span>{qc._pill("info", "INFO")}</div>{body}</section>'
    )


def _account_value(value: int, total: int) -> str:
    suffix = f' <small>{value / total * 100:.2f}%</small>' if total else ""
    return f"{value:,}{suffix}"


def _render_accounting(stats: Mapping[str, object] | None) -> str:
    if not stats:
        body = f'<div class="empty">{qc._t("Trimming statistics were not found.")}</div>'
    else:
        reads, bases = _nested(stats, "reads"), _nested(stats, "bases")
        trimming, filtering, output = (_nested(stats, "trimming"), _nested(stats, "filtering"), _nested(stats, "output"))
        input_reads, kept_reads, removed_reads = (_number(reads.get("total")), _number(reads.get("passed")), _number(reads.get("discarded")))
        input_bases, kept_bases, removed_bases = (_number(bases.get("before")), _number(bases.get("after")), _number(bases.get("trimmed")))
        entries = (
            ("Reads before", _account_value(input_reads, input_reads)),
            ("Reads after", _account_value(kept_reads, input_reads)),
            ("Reads discarded", _account_value(removed_reads, input_reads)),
            ("Too short after trimming", _account_value(_number(filtering.get("too_short_reads")), input_reads)),
            ("Adapter-trimmed reads", _account_value(_number(trimming.get("adapter_trimmed_reads")), input_reads)),
            ("Quality-trimmed reads", _account_value(_number(trimming.get("quality_trimmed_reads")), input_reads)),
            ("PolyG-trimmed reads", _account_value(_number(trimming.get("polyG_trimmed_reads")), input_reads)),
            ("PolyX-trimmed reads", _account_value(_number(trimming.get("polyX_trimmed_reads")), input_reads)),
            ("Bases before", _account_value(input_bases, input_bases)),
            ("Bases after", _account_value(kept_bases, input_bases)),
            ("Bases removed", _account_value(removed_bases, input_bases)),
            ("Corrected reads", f'{_number(trimming.get("corrected_reads")):,}'),
            ("Corrected bases", f'{_number(trimming.get("corrected_bases")):,}'),
            ("Merged pairs", f'{_number(trimming.get("merged_pairs")):,}'),
            ("Output reads", f'{_number(output.get("reads")):,}'),
            ("Discarded pairs", f'{_number(output.get("discarded_pairs")):,}'),
        )
        rows = "".join(
            f'<tr><th scope="row" class="row-h">{qc._t(label)}</th><td class="n">{value}</td></tr>'
            for label, value in entries
        )
        body = f'<div class="table-wrap"><table class="grid acc"><thead><tr><th>{qc._t("Metric")}</th><th class="n">{qc._t("All reads")}</th></tr></thead><tbody>{rows}</tbody></table></div>'
    return (
        '<section class="panel foldable" id="trim-accounting" data-section data-sev="info" '
        f'{qc._keys(2, 2)}><div class="phead" data-fold>{qc.CHEVRON}<h2>{qc._t("Reads and bases")}</h2>'
        f'<span class="sp"></span>{qc._pill("info", "INFO")}</div>{body}</section>'
    )


def _status_change(before: str, after: str) -> str:
    severity = {"not_evaluated": 0, "pass": 1, "warning": 2, "fail": 3}
    if before == after:
        return '<span class="chg same" title="unchanged">=</span>'
    return ('<span class="chg good" title="improved">▲</span>' if severity[after] < severity[before]
            else '<span class="chg bad" title="worsened">▼</span>')


def _render_overview(
    before: qc.QcReportModel,
    after: qc.QcReportModel,
    modules: list[str],
    before_modules: Mapping[str, qc.PlotModule],
    after_modules: Mapping[str, qc.PlotModule],
    metric_ids: list[str],
    reads: tuple[str, ...],
    stats: Mapping[str, object] | None,
) -> str:
    statistics_before, statistics_after = _by_read(before), _by_read(after)
    read_counts = _nested(stats or {}, "reads")
    base_counts = _nested(stats or {}, "bases")
    total_reads, kept_reads, dropped_reads = (_number(read_counts.get("total")), _number(read_counts.get("passed")), _number(read_counts.get("discarded")))
    total_bases, removed_bases = (_number(base_counts.get("before")), _number(base_counts.get("trimmed")))
    kpis = [
        ("Input reads", f"{total_reads:,}", ""),
        ("Kept reads", f"{kept_reads:,}", f"{kept_reads / total_reads * 100:.2f}%" if total_reads else ""),
        ("Discarded reads", f"{dropped_reads:,}", f"{dropped_reads / total_reads * 100:.2f}%" if total_reads else ""),
        ("Bases removed", f"{removed_bases:,}", f"{removed_bases / total_bases * 100:.2f}%" if total_bases else ""),
    ]
    for read in reads:
        before_length = dict(statistics_before.get(read, qc.BasicStatistics(read, "", ())).items).get("Avg length", "—")
        after_length = dict(statistics_after.get(read, qc.BasicStatistics(read, "", ())).items).get("Avg length", "—")
        kpis.append((f"Mean length {read}", f"{before_length} → {after_length}", "bp"))
    kpi_html = "".join(
        f'<div class="kpi"><small>{qc._t(label)}</small><strong>{escape(value)}</strong><em>{escape(note)}</em></div>'
        for label, value, note in kpis
    )
    before_counts, after_counts = dict(before.qc_counts), dict(after.qc_counts)
    count_rows = "".join(
        f'<tr><th scope="row">{qc._pill(status)}</th><td class="n">{before_counts.get(status, 0)}</td><td class="n">{after_counts.get(status, 0)}</td></tr>'
        for status in qc.STATUS_ORDER
    )
    status_rows = []
    for metric_id in metric_ids:
        title = (before_modules.get(metric_id) or after_modules[metric_id]).title
        before_cards, after_cards = _cards(before_modules.get(metric_id)), _cards(after_modules.get(metric_id))
        cells = "".join(
            f'<td class="c">{qc._cell(_status(before_cards.get(read)), qc.STATUS_LABELS[_status(before_cards.get(read))])}</td>'
            f'<td class="c">{qc._cell(_status(after_cards.get(read)), qc.STATUS_LABELS[_status(after_cards.get(read))])}</td>'
            f'<td class="c">{_status_change(_status(before_cards.get(read)), _status(after_cards.get(read)))}</td>'
            for read in reads
        )
        status_rows.append(
            f'<tr><th scope="row" class="row-h"><a href="#cmp-{escape(metric_id)}">{qc._t(title)}</a></th>{cells}</tr>'
        )
    read_head = "".join(f'<th colspan="3" class="c">{escape(read)}</th>' for read in reads)
    stage_head = "".join(
        f'<th class="c">{qc._t("Before")}</th><th class="c">{qc._t("After")}</th><th class="c">{qc._t("Change")}</th>'
        for _ in reads
    )
    return (
        f'<section class="panel" id="trim-overview" data-section {qc._keys(0, 0)}><div class="hero"><div class="id">'
        f'<div class="eyebrow">{qc._t("Trimming report")}</div><h1>{escape(before.sample_id)}</h1><div class="kpis">{kpi_html}</div>'
        f'<div class="notice">{qc._t("Trimming removes adapters, low-quality ends and too short reads. QC after trimming is calculated on the trimmed FASTQ with the same ruleset; a status change shows the effect of trimming.")}</div></div>'
        f'<div class="qcbox"><div class="top"><div><h2>{qc._t("QC status before and after trimming")}</h2><div class="rs">{qc._t("Ruleset")}: {escape(after.ruleset_label)}</div></div>'
        f'<div class="overall">{qc._pill(before.overall_qc_status)} <span class="arrow">→</span> {qc._pill(after.overall_qc_status)}</div></div>'
        f'<div class="table-wrap"><table class="grid counts-cmp"><thead><tr><th>{qc._t("Status")}</th><th class="n">{qc._t("Before trimming")}</th><th class="n">{qc._t("After trimming")}</th></tr></thead><tbody>{count_rows}</tbody></table></div></div></div>'
        f'<div class="table-wrap"><table class="grid hm status-cmp"><thead><tr><th rowspan="2">{qc._t("QC module")}</th>{read_head}</tr><tr>{stage_head}</tr></thead><tbody>{"".join(status_rows)}</tbody></table></div></section>'
    )


def _render_key_comparisons(
    before_modules: Mapping[str, qc.PlotModule],
    after_modules: Mapping[str, qc.PlotModule],
    before_plots: Path,
    after_plots: Path,
    reads: tuple[str, ...],
) -> str:
    selected = ("per_base_quality", "adapter_content", "sequence_length_distribution", "per_sequence_quality")
    groups = []
    for metric_id in selected:
        before_module, after_module = before_modules.get(metric_id), after_modules.get(metric_id)
        if not before_module and not after_module:
            continue
        title = (before_module or after_module).title
        cards_before, cards_after = _cards(before_module), _cards(after_module)
        lanes = "".join(
            _render_lane(read, cards_before.get(read), cards_after.get(read), before_plots, after_plots)
            for read in reads
        )
        groups.append(f'<div class="cmp-group"><h3>{qc._t(title)}</h3><div class="tracks compare-tracks">{lanes}</div></div>')
    body = "".join(groups) if groups else f'<div class="empty">{qc._t("No comparison charts are available.")}</div>'
    return (
        '<section class="panel foldable" id="trim-compare" data-section data-sev="info" '
        f'{qc._keys(4, 4)}><div class="phead" data-fold>{qc.CHEVRON}<h2>{qc._t("Before / after comparison")}</h2>'
        f'<span class="sp"></span>{qc._pill("info", "INFO")}</div><div class="cmp-legend">{qc._t("Before trimming")} <span class="arrow">→</span> {qc._t("After trimming")}</div>{body}</section>'
    )


def _render_lane(
    read: str,
    before_card: qc.PlotCard | None,
    after_card: qc.PlotCard | None,
    before_plots: Path,
    after_plots: Path,
) -> str:
    labels: list[str] = []
    for card in (before_card, after_card):
        if card:
            for label, _, _ in card.qc_observations:
                if label not in labels:
                    labels.append(label)
    observations = "".join(
        f'<tr><td>{qc._t(label)}</td><td class="n">{escape(_value(before_card, label))}</td>'
        f'<td class="n">{escape(_value(after_card, label))}</td></tr>'
        for label in labels
    )
    metrics = (
        '<div class="metrics"><table><thead><tr><th></th>'
        f'<th class="n">{qc._t("Before")}</th><th class="n">{qc._t("After")}</th></tr></thead><tbody>{observations}</tbody></table></div>'
        if observations else ""
    )
    after_reasons = "".join(
        f'<div class="why">{qc._t(reason)}</div>' for reason in (after_card.qc_reasons if after_card else ())
    )
    return (
        '<div class="lane compare-lane"><div class="lane-h"><span>' + escape(read) + '</span>'
        + qc._pill(_status(before_card)) + '<span class="arrow">→</span>' + qc._pill(_status(after_card))
        + '</div><div class="stage-charts"><div><h3>' + qc._t("Before trimming") + '</h3>'
        + _stage_chart(before_card, before_plots, "before") + '</div><div><h3>' + qc._t("After trimming") + '</h3>'
        + _stage_chart(after_card, after_plots, "after") + f'</div></div>{metrics}{after_reasons}</div>'
    )


def _render_module(
    metric_id: str,
    title: str,
    before_module: qc.PlotModule | None,
    after_module: qc.PlotModule | None,
    before_plots: Path,
    after_plots: Path,
    index: int,
    reads: tuple[str, ...],
) -> str:
    before_cards, after_cards = _cards(before_module), _cards(after_module)
    statuses = [_status(after_cards.get(read)) for read in reads]
    severity = max(statuses, key={"not_evaluated": 0, "pass": 1, "warning": 2, "fail": 3}.get)
    panel_severity = severity if qc._needs_attention(severity) else "ok"
    lanes = "".join(
        _render_lane(read, before_cards.get(read), after_cards.get(read), before_plots, after_plots)
        for read in reads
    )
    return (
        f'<section class="panel foldable" id="cmp-{escape(metric_id)}" data-section data-sev="{panel_severity}" '
        f'{qc._keys(index, qc._triage_rank(severity, index))}><div class="phead" data-fold>{qc.CHEVRON}'
        f'<h2>{qc._t(title)}</h2><span class="sp"></span>{qc._pill(severity)}</div><div class="tracks compare-tracks">{lanes}</div></section>'
    )


_TRIMMING_CSS = """
.kpis{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));border-top:1px solid var(--line);border-left:1px solid var(--line);margin-top:14px}
.kpi{padding:9px;border-right:1px solid var(--line);border-bottom:1px solid var(--line)}.kpi small,.kpi em{display:block;font:11px var(--mono);color:var(--muted);font-style:normal}.kpi strong{display:block;font:600 16px var(--sans);margin:3px 0}
.qcbox{padding:14px}.qcbox .top{display:flex;justify-content:space-between;gap:12px;align-items:flex-start}.qcbox h2{font:600 15px var(--sans);margin:0 0 3px}.qcbox .rs{font:11px var(--mono);color:var(--muted)}.overall{white-space:nowrap}.status-cmp .c,.grid .c{text-align:center}.status-cmp .cell{margin:auto}.chg{font:600 14px var(--mono)}.chg.good{color:var(--pass)}.chg.bad{color:var(--fail)}.chg.same{color:var(--muted)}
.cmp-legend{padding:10px 12px;border-bottom:1px solid var(--line);font:12px var(--mono);color:var(--muted)}.cmp-group{padding:12px;border-bottom:1px solid var(--line)}.cmp-group:last-child{border-bottom:0}.cmp-group>h3{margin:0 0 8px;font:600 14px var(--sans)}
.stage-charts{display:grid;grid-template-columns:minmax(0,1fr) minmax(0,1fr);gap:12px;padding:10px 0}
.stage-charts>div+div{border-left:1px solid var(--line);padding-left:12px}
.stage-charts h3{font:600 11px var(--mono);text-transform:uppercase;letter-spacing:.035em;color:var(--muted);margin:0 0 7px}
.compare-lane{grid-column:span 2}.compare-tracks{grid-template-columns:minmax(0,1fr)}
.compare-lane .metrics table{width:100%}.compare-lane .metrics th,.compare-lane .metrics td{padding:5px 0}
.params .output-path{font:11px var(--mono);overflow-wrap:anywhere}
.stage-report-links{padding:12px 10px 8px;border-bottom:1px solid var(--line)}.stage-report-links small{display:block;margin:0 6px 6px;font:600 10px var(--mono);letter-spacing:.05em;text-transform:uppercase;color:var(--muted)}.stage-report-links a{display:flex;align-items:center;justify-content:space-between;gap:8px;padding:7px 6px;border-radius:3px;color:var(--ink);font:12px var(--sans);text-decoration:none}.stage-report-links a:hover{background:var(--row-hover);color:var(--brand)}.stage-report-links a::after{content:"↗";font:13px var(--mono);color:var(--muted)}
@media(max-width:700px){.kpis{grid-template-columns:minmax(0,1fr)}.qcbox .top{display:block}.overall{margin-top:8px}.stage-charts{grid-template-columns:minmax(0,1fr)}.stage-charts>div+div{border-left:0;border-top:1px solid var(--line);padding-left:0;padding-top:12px}.compare-lane{grid-column:auto}}
"""


def render_trimming_report(
    before: qc.QcReportModel,
    after: qc.QcReportModel,
    before_plots: Path,
    after_plots: Path,
    trimming_stats: Mapping[str, object] | None = None,
    run_manifest: Mapping[str, object] | None = None,
    output_fastqs: Mapping[str, Path] | None = None,
    stage_report_links: Mapping[str, Path] | None = None,
) -> str:
    if before.sample_id != after.sample_id:
        raise qc.QcReportError("before and after reports have different sample identifiers")
    reads = tuple(read for read in ("R1", "R2") if read in before.reads or read in after.reads)
    before_modules, after_modules = _modules(before), _modules(after)
    metric_ids = list(before_modules)
    metric_ids.extend(metric_id for metric_id in after_modules if metric_id not in before_modules)

    def after_status(metric_id: str) -> str:
        statuses = [_status(_cards(after_modules.get(metric_id)).get(read)) for read in reads]
        return max(statuses, key={"not_evaluated": 0, "pass": 1, "warning": 2, "fail": 3}.get)

    modules = [
        _render_module(
            metric_id,
            (before_modules.get(metric_id) or after_modules[metric_id]).title,
            before_modules.get(metric_id), after_modules.get(metric_id), before_plots, after_plots, index + 5, reads,
        )
        for index, metric_id in enumerate(metric_ids)
    ]
    tracks = [
        f'<a class="track" data-nav href="#trim-overview" {qc._keys(0, 0)}>'
        f'<span class="l">{qc._t("Trimming overview")}</span></a>',
        f'<a class="track" data-nav href="#trim-parameters" {qc._keys(1, 1)}><span class="l">{qc._t("Trimming parameters")}</span></a>',
        f'<a class="track" data-nav href="#trim-accounting" {qc._keys(2, 2)}><span class="l">{qc._t("Reads and bases")}</span></a>',
        f'<a class="track" data-nav href="#trim-basic" {qc._keys(3, 3)}><span class="l">{qc._t("Basic statistics: before / after")}</span></a>',
        f'<a class="track" data-nav href="#trim-compare" {qc._keys(4, 4)}><span class="l">{qc._t("Before / after comparison")}</span></a>',
    ]
    tracks.extend(
        f'<a class="track" data-nav href="#cmp-{escape(metric_id)}" '
        f'{qc._keys(index + 5, qc._triage_rank(after_status(metric_id), index + 5))}>'
        f'<span class="l">{qc._t((before_modules.get(metric_id) or after_modules[metric_id]).title)}</span></a>'
        for index, metric_id in enumerate(metric_ids)
    )
    reads_head = "".join(f"<b>{escape(read)}</b>" for read in reads)
    legend = "".join(
        f'<span>{qc._cell(status, qc.STATUS_LABELS[status])}{qc.STATUS_LABELS[status]}</span>'
        for status in qc.STATUS_ORDER
    )
    prefs = (
        f'<div class="prefs"><div class="pref"><small>{qc._t("Theme")}</small><div class="seg" role="group">'
        f'<button type="button" data-theme-btn="light" aria-pressed="true" title="Light">{qc.SUN}</button>'
        f'<button type="button" data-theme-btn="dark" aria-pressed="false" title="Dark">{qc.MOON}</button></div></div>'
        f'<div class="pref"><small>{qc._t("Language")}</small><div class="seg" role="group">'
        '<button type="button" data-lang-btn="en" aria-pressed="true">EN</button>'
        '<button type="button" data-lang-btn="ru" aria-pressed="false">RU</button></div></div></div>'
    )
    report_links = "".join(
        f'<a href="{escape(path.as_posix())}">{qc._t(label)}</a>'
        for stage, label in (
            ("before", "QC before trimming"),
            ("after", "QC after trimming"),
        )
        if (path := (stage_report_links or {}).get(stage)) is not None
    )
    stage_reports = (
        '<div class="stage-report-links" aria-label="QC reports">'
        f'<small>{qc._t("QC reports")}</small>{report_links}</div>'
        if report_links else ""
    )
    rail = (
        f'<aside class="rail">{prefs}{stage_reports}<button class="rail-h" type="button" aria-expanded="false" aria-controls="tracklist">'
        f'<span>{qc._t("Report sections")}</span><em>{reads_head}</em></button><nav class="tracklist" id="tracklist">'
        f'{"".join(tracks)}</nav><div class="legend">{legend}</div></aside>'
    )
    appbar = (
        '<header class="appbar"><div class="logo"><i aria-hidden="true"><b></b><b></b><b></b><b></b></i>NeoQC Trimming Report</div>'
        '<div class="crumbs">'
        f'<span class="cr"><small>{qc._t("Sample")}</small><strong>{escape(before.sample_id)}</strong></span>'
        f'<span class="cr"><small>{qc._t("Reads")}</small><strong>{escape(" / ".join(reads))}</strong></span>'
        f'<span class="cr"><small>{qc._t("Generated")}</small><strong>{escape(after.generated_at)}</strong></span>'
        f'<span class="cr"><small>{qc._t("Ruleset")}</small><strong>{escape(after.ruleset_label)}</strong></span></div>'
        f'<div class="seg" role="group"><small>{qc._t("Order")}</small>'
        f'<button type="button" data-order-btn="triage" aria-pressed="true">{qc._t("Triage")}</button>'
        f'<button type="button" data-order-btn="sequential" aria-pressed="false">{qc._t("Sequential")}</button></div>'
        f'<button class="btn ghost" id="fold-all" type="button" data-state="collapse"><span class="fa-c">{qc._t("Collapse all")}</span>'
        f'<span class="fa-e">{qc._t("Expand all")}</span></button><button class="btn print" type="button" onclick="window.print()">{qc._t("Print / Save PDF")}</button></header>'
    )
    overview = _render_overview(before, after, modules, before_modules, after_modules, metric_ids, reads, trimming_stats)
    main = (
        overview
        + _render_parameters(run_manifest, reads, output_fastqs)
        + _render_accounting(trimming_stats)
        + _render_statistics(before, after)
        + _render_key_comparisons(before_modules, after_modules, before_plots, after_plots, reads)
        + "".join(modules)
        + f'<footer {qc._keys(99999, 99999)}>{qc._t("Generated by NeoQC")}</footer>'
    )
    dialog = (
        '<dialog id="chart-dialog" aria-label="Chart"><div class="dialog-head">'
        f'<button class="dialog-close" id="dialog-close" type="button">{qc._t("Close")}</button></div><div id="dialog-body"></div></dialog>'
    )
    return (
        '<!doctype html>\n<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">'
        f'<title>NeoQC trimming report — {escape(before.sample_id)}</title><style>{qc._font_face_css()}\n{qc._CSS}\n{_TRIMMING_CSS}</style>'
        f'<script>{qc._HEAD_JS}</script></head><body>{appbar}<div class="shell">{rail}<main class="main">{main}</main></div>{dialog}'
        f'<script>{qc._JS}</script></body></html>'
    )


def generate_trimming_report(result_dir: Path, output_path: Path | None = None) -> Path:
    result_dir = result_dir.resolve()
    before_dir, after_dir = result_dir / "qc" / "before", result_dir / "qc" / "after"
    before_plots, after_plots = result_dir / "plots" / "before", result_dir / "plots" / "after"
    before = qc.load_report_model(before_dir, before_plots)
    after = qc.load_report_model(after_dir, after_plots)
    trimming_stats = _load_json_object(result_dir / "trimming_report.json")
    run_manifest = _load_json_object(result_dir / "run_manifest.json")
    root_sample_id = _sample_id_from_run_manifest(run_manifest) or result_dir.name
    before = _with_root_sample_id_for_missing_summaries(before, root_sample_id)
    after = _with_root_sample_id_for_missing_summaries(after, root_sample_id)
    trimmed_dir = result_dir / "trimmed"
    output_fastqs = {
        read: path.relative_to(result_dir)
        for read, path in (
            ("R1", trimmed_dir / "R1.trimmed.fastq.gz"),
            ("R2", trimmed_dir / "R2.trimmed.fastq.gz"),
            ("merged", trimmed_dir / "merged.fastq.gz"),
        )
        if path.is_file()
    }
    stage_report_links = {
        stage: report_path.relative_to(result_dir)
        for stage, report_path in (
            ("before", before_dir / qc.REPORT_FILENAME),
            ("after", after_dir / qc.REPORT_FILENAME),
        )
        if report_path.is_file()
    }
    output_path = (output_path or result_dir / REPORT_FILENAME).resolve()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(
        render_trimming_report(
            before, after, before_plots, after_plots, trimming_stats, run_manifest, output_fastqs,
            stage_report_links,
        ),
        encoding="utf-8",
    )
    return output_path
