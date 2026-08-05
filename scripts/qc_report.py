"""Self-contained, FastQC-style report model and renderer for NeoQC charts."""

from __future__ import annotations

import base64
import json
from collections import OrderedDict
from dataclasses import dataclass
from datetime import datetime, timedelta, timezone
from html import escape
from pathlib import Path
from typing import Mapping, Sequence


REPORT_FILENAME = "neoqc_qc_report.html"
MOSCOW_TIME = timezone(timedelta(hours=3), name="МСК")


class QcReportError(ValueError):
    """Raised when plot artifacts cannot form a consistent QC report."""


@dataclass(frozen=True)
class BasicStatistics:
    read: str
    sample_label: str
    items: tuple[tuple[str, str], ...]


@dataclass(frozen=True)
class PlotCard:
    metric_id: str
    read: str
    title: str
    status: str
    reason: str = ""
    alt_text: str = ""
    svg: str = ""
    png: str = ""


@dataclass(frozen=True)
class PlotModule:
    metric_id: str
    title: str
    cards: tuple[PlotCard, ...]

    @property
    def availability(self) -> str:
        statuses = {card.status for card in self.cards}
        if "error" in statuses:
            return "error"
        if "generated" in statuses:
            return "ready"
        return "not_run"


@dataclass(frozen=True)
class QcReportModel:
    sample_id: str
    generated_at: str
    reads: tuple[str, ...]
    basic_statistics: tuple[BasicStatistics, ...]
    modules: tuple[PlotModule, ...]
    generated_plots: int
    errors: int


def _text(value: object, path: str, *, required: bool = False) -> str:
    if value is None:
        result = ""
    elif isinstance(value, (str, int, float, bool)):
        result = str(value)
    else:
        raise QcReportError(f"{path} must be a scalar value")
    if required and not result.strip():
        raise QcReportError(f"{path} is required")
    return result


def _mapping(value: object, path: str) -> Mapping[str, object]:
    if not isinstance(value, Mapping):
        raise QcReportError(f"{path} must be an object")
    return value


def _sequence(value: object, path: str) -> Sequence[object]:
    if isinstance(value, (str, bytes)) or not isinstance(value, Sequence):
        raise QcReportError(f"{path} must be an array")
    return value


def _parse_summary(path: Path, read: str) -> BasicStatistics:
    try:
        lines = path.read_text(encoding="utf-8").splitlines()
    except OSError as error:
        raise QcReportError(f"cannot read {path.name}: {error}") from error

    sample_label = path.name.removesuffix(f"_{read}_summary.txt")
    items: list[tuple[str, str]] = []
    in_base_composition = False
    for line in lines:
        stripped = line.strip()
        if not stripped or stripped.startswith("==="):
            continue
        if stripped == "Base composition":
            in_base_composition = True
            continue
        if ":" not in stripped:
            continue
        label, value = (part.strip() for part in stripped.split(":", 1))
        if in_base_composition and label in {"A", "C", "G", "T", "N"}:
            label = f"Base {label}"
        items.append((label, value))
    return BasicStatistics(read=read, sample_label=sample_label, items=tuple(items))


def _find_summary(result_dir: Path, read: str) -> BasicStatistics | None:
    candidates = sorted(result_dir.glob(f"*_{read}_summary.txt"))
    if not candidates:
        return None
    if len(candidates) > 1:
        raise QcReportError(f"multiple {read} summary files found in {result_dir}")
    return _parse_summary(candidates[0], read)


def load_report_model(result_dir: Path, plot_dir: Path | None = None) -> QcReportModel:
    result_dir = result_dir.resolve()
    plot_dir = (plot_dir or result_dir / "plots").resolve()
    manifest_path = plot_dir / "plots_manifest.json"
    try:
        raw_manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    except OSError as error:
        raise QcReportError(f"cannot read {manifest_path}: {error}") from error
    except json.JSONDecodeError as error:
        raise QcReportError(f"invalid JSON in {manifest_path}: {error}") from error

    manifest = _mapping(raw_manifest, "manifest")
    if manifest.get("schema_version") != 1:
        raise QcReportError("unsupported plots manifest schema")

    raw_plots = _sequence(manifest.get("plots"), "manifest.plots")
    cards: list[PlotCard] = []
    for index, raw_entry in enumerate(raw_plots):
        entry = _mapping(raw_entry, f"manifest.plots[{index}]")
        status = _text(entry.get("status"), f"manifest.plots[{index}].status", required=True)
        if status not in {"generated", "skipped", "error"}:
            raise QcReportError(f"manifest.plots[{index}].status is invalid")
        read = _text(entry.get("read"), f"manifest.plots[{index}].read", required=True)
        if read not in {"R1", "R2"}:
            raise QcReportError(f"manifest.plots[{index}].read must be R1 or R2")
        cards.append(
            PlotCard(
                metric_id=_text(entry.get("id"), f"manifest.plots[{index}].id", required=True),
                read=read,
                title=_text(entry.get("title"), f"manifest.plots[{index}].title", required=True),
                status=status,
                reason=_text(entry.get("reason"), f"manifest.plots[{index}].reason"),
                alt_text=_text(entry.get("alt_text"), f"manifest.plots[{index}].alt_text"),
                svg=_text(entry.get("svg"), f"manifest.plots[{index}].svg"),
                png=_text(entry.get("png"), f"manifest.plots[{index}].png"),
            )
        )

    summaries = tuple(
        summary
        for read in ("R1", "R2")
        if (summary := _find_summary(result_dir, read)) is not None
    )
    active_reads = {
        card.read
        for card in cards
        if card.status != "skipped" or card.reason != "source_not_found"
    }
    active_reads.update(summary.read for summary in summaries)
    reads = tuple(read for read in ("R1", "R2") if read in active_reads)
    if not reads:
        raise QcReportError("report contains neither R1 nor R2 results")

    grouped: OrderedDict[str, list[PlotCard]] = OrderedDict()
    for card in cards:
        if card.read not in reads:
            continue
        grouped.setdefault(card.metric_id, []).append(card)
    modules = tuple(
        PlotModule(metric_id=metric_id, title=module_cards[0].title, cards=tuple(module_cards))
        for metric_id, module_cards in grouped.items()
    )

    sample_ids = {summary.sample_label for summary in summaries if summary.sample_label}
    sample_id = next(iter(sample_ids)) if len(sample_ids) == 1 else result_dir.name
    generated_plots = sum(card.status == "generated" for card in cards if card.read in reads)
    errors = sum(card.status == "error" for card in cards if card.read in reads)
    return QcReportModel(
        sample_id=sample_id,
        generated_at=datetime.now(MOSCOW_TIME).strftime(
            "%d.%m.%Y, %H:%M МСК"
        ),
        reads=reads,
        basic_statistics=summaries,
        modules=modules,
        generated_plots=generated_plots,
        errors=errors,
    )


def _safe_asset(plot_dir: Path, filename: str, expected_suffix: str) -> Path:
    if not filename:
        raise QcReportError("asset filename is missing")
    root = plot_dir.resolve()
    candidate = (root / filename).resolve()
    try:
        candidate.relative_to(root)
    except ValueError as error:
        raise QcReportError("asset path escapes the plot directory") from error
    if candidate.suffix.lower() != expected_suffix:
        raise QcReportError(f"unexpected asset type: {candidate.suffix}")
    if not candidate.is_file():
        raise QcReportError(f"asset not found: {filename}")
    return candidate


def _asset_data_uri(plot_dir: Path, card: PlotCard) -> str:
    candidates = ((card.svg, ".svg", "image/svg+xml"), (card.png, ".png", "image/png"))
    errors: list[str] = []
    for filename, suffix, mime_type in candidates:
        if not filename:
            continue
        try:
            path = _safe_asset(plot_dir, filename, suffix)
            encoded = base64.b64encode(path.read_bytes()).decode("ascii")
            return f"data:{mime_type};base64,{encoded}"
        except (OSError, QcReportError) as error:
            errors.append(str(error))
    raise QcReportError("; ".join(errors) or "no chart asset is available")


def _badge(status: str, label: str | None = None) -> str:
    labels = {"ready": "READY", "not_run": "NOT RUN", "error": "ERROR", "info": "INFO"}
    return f'<span class="badge {escape(status)}">{escape(label or labels[status])}</span>'


def _render_basic_statistics(model: QcReportModel) -> str:
    cards: list[str] = []
    for stats in model.basic_statistics:
        rows = "".join(
            f"<tr><th>{escape(label)}</th><td>{escape(value)}</td></tr>"
            for label, value in stats.items
        )
        cards.append(
            f'<article class="stats-card"><h3>{escape(stats.read)}</h3>'
            f'<table><tbody>{rows}</tbody></table></article>'
        )
    if not cards:
        cards.append('<p class="empty">Summary files were not found.</p>')
    return (
        '<section class="module" id="basic-statistics">'
        '<div class="module-head"><h2>Basic Statistics</h2>' + _badge("info") + '</div>'
        '<div class="stats-grid">' + "".join(cards) + '</div></section>'
    )


def _render_plot_card(card: PlotCard, plot_dir: Path) -> str:
    if card.status == "generated":
        try:
            source = _asset_data_uri(plot_dir, card)
            image = (
                f'<button class="chart-button" type="button" aria-label="Open {escape(card.title)} {escape(card.read)}">'
                f'<img class="chart-image" src="{source}" alt="{escape(card.alt_text or card.title)}" loading="lazy">'
                '<span class="expand-hint">Open full size</span></button>'
            )
            state = _badge("ready")
        except QcReportError as error:
            image = f'<div class="empty error-box">Chart asset unavailable: {escape(str(error))}</div>'
            state = _badge("error")
    else:
        reasons = {
            "adapter_analysis_disabled": "Adapter analysis was disabled for this run.",
            "source_not_found": "The source TSV was not produced.",
        }
        message = reasons.get(card.reason, card.reason or "Chart was not generated.")
        image = f'<div class="empty">{escape(message)}</div>'
        state = _badge("error" if card.status == "error" else "not_run")
    return (
        '<article class="plot-card">'
        f'<div class="plot-card-head"><h3>{escape(card.read)}</h3>{state}</div>'
        f'{image}</article>'
    )


def _render_module(module: PlotModule, plot_dir: Path) -> str:
    cards = "".join(_render_plot_card(card, plot_dir) for card in module.cards)
    return (
        f'<section class="module" id="module-{escape(module.metric_id)}">'
        '<div class="module-head">'
        f'<h2>{escape(module.title)}</h2>{_badge(module.availability)}</div>'
        f'<div class="plot-grid">{cards}</div></section>'
    )


_CSS = """
:root{--ink:#0a132d;--muted:#657285;--line:#dce3ea;--panel:#f4f7f9;--brand:#2947a0;--brand-dark:#192f70;--accent:#539d96;--white:#fff;--danger:#c84a5a;--warn:#e69f00}
*{box-sizing:border-box}html{scroll-behavior:auto;scroll-padding-top:20px}body{margin:0;background:#edf1f5;color:var(--ink);font:15px/1.55 -apple-system,BlinkMacSystemFont,"Segoe UI",Arial,sans-serif}.layout{display:grid;grid-template-columns:280px minmax(0,1fr);max-width:1540px;margin:0 auto;min-height:100vh}.sidebar{position:sticky;top:0;height:100vh;padding:28px 22px;color:#fff;overflow:auto;background-color:var(--ink);background-image:radial-gradient(circle at 18% 8%,#365fc03d 0,transparent 27%),linear-gradient(#ffffff08 1px,transparent 1px),linear-gradient(90deg,#ffffff08 1px,transparent 1px);background-size:auto,28px 28px,28px 28px;background-position:0 0,-1px -1px,-1px -1px;border-top:3px solid #7190e2}.brand{display:flex;align-items:center;gap:11px;font-size:19px;font-weight:800;letter-spacing:.01em}.brand-mark{display:grid;place-items:center;width:34px;height:34px;border-radius:9px;background:linear-gradient(145deg,#4268c9,#2947a0);color:#fff;box-shadow:0 8px 22px #0004}.sidebar-meta{margin:24px 0;padding:16px;border:1px solid #ffffff24;border-radius:11px;background:#101b38cc;box-shadow:inset 0 1px #ffffff0d,0 12px 28px #0002;backdrop-filter:blur(4px)}.sidebar-meta small{display:block;color:#aebbd1}.sidebar-meta strong{display:block;margin:3px 0 12px;overflow-wrap:anywhere}.sidebar-meta strong:last-child{margin-bottom:0}.nav-title{margin:22px 8px 8px;color:#91a9de;font-size:11px;font-weight:800;letter-spacing:.13em;text-transform:uppercase}.nav{list-style:none;margin:0;padding:0}.nav a{display:grid;grid-template-columns:24px minmax(0,1fr) 8px;align-items:center;gap:8px;margin:2px 0;padding:9px 10px;border:1px solid transparent;border-radius:8px;color:#dce3ef;text-decoration:none;font-size:13px;transition:background .12s,border-color .12s}.nav a:hover,.nav a:focus{background:#ffffff10;border-color:#ffffff13;color:#fff}.nav a.active{background:linear-gradient(90deg,#3153a650,#ffffff0b);border-color:#7894db42;color:#fff;box-shadow:inset 3px 0 #7190e2}.nav-index{color:#8395b4;font:700 10px/1.2 ui-monospace,SFMono-Regular,Consolas,monospace;letter-spacing:.04em}.nav a.active .nav-index{color:#aac0fa}.nav-label{min-width:0}.nav-dot{width:7px;height:7px;border-radius:50%;background:var(--accent);box-shadow:0 0 0 3px #539d961c;flex:0 0 auto}.nav-dot.info{background:#7190e2}.nav-dot.not_run{background:#77849a}.nav-dot.error{background:var(--danger)}.main{min-width:0;padding:34px}.report-head,.module{background:#fff;border:1px solid var(--line);border-radius:14px;box-shadow:0 8px 28px #0a132d0d}.report-head{padding:34px 38px;border-top:6px solid var(--brand)}.eyebrow{margin:0 0 7px;color:var(--brand);font-size:11px;font-weight:800;letter-spacing:.12em;text-transform:uppercase}h1{margin:0;font-size:30px;line-height:1.2}h2{margin:0;font-size:21px}h3{margin:0;font-size:15px}.subtitle{max-width:80ch;margin:12px 0 0;color:var(--muted)}.overview{display:grid;grid-template-columns:repeat(3,minmax(0,1fr));gap:12px;margin-top:24px}.overview-card{padding:15px 17px;border:1px solid var(--line);border-radius:9px;background:var(--panel)}.overview-card span{display:block;color:var(--muted);font-size:12px}.overview-card strong{display:block;margin-top:3px;font-size:21px}.notice{margin-top:18px;padding:12px 14px;border-left:3px solid var(--accent);background:#e8f3f1;color:var(--brand-dark);font-size:13px}.module{margin-top:18px;padding:28px 30px;scroll-margin-top:20px}.module-head,.plot-card-head{display:flex;align-items:center;justify-content:space-between;gap:18px}.badge{display:inline-flex;align-items:center;border-radius:999px;padding:4px 9px;font-size:10px;font-weight:800;letter-spacing:.06em;white-space:nowrap}.badge.ready{background:#e8f3f1;color:#276c66}.badge.info{background:#eef1f8;color:var(--brand-dark)}.badge.not_run{background:#f1f3f5;color:var(--muted)}.badge.error{background:#fae9ec;color:#9f3041}.stats-grid,.plot-grid{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:16px;margin-top:20px}.stats-card,.plot-card{min-width:0;border:1px solid var(--line);border-radius:10px;overflow:hidden}.stats-card h3,.plot-card-head{padding:12px 15px;background:var(--panel);border-bottom:1px solid var(--line)}table{width:100%;border-collapse:collapse;font-size:13px}th,td{padding:8px 13px;border-bottom:1px solid var(--line);text-align:left}th{width:58%;color:var(--muted);font-weight:500}tr:last-child th,tr:last-child td{border-bottom:0}.chart-button{position:relative;display:block;width:100%;padding:0;border:0;background:#fff;cursor:zoom-in}.chart-image{display:block;width:100%;height:auto}.expand-hint{position:absolute;right:10px;bottom:10px;padding:5px 8px;border-radius:6px;background:#0a132ddd;color:#fff;font-size:10px;opacity:0;transition:opacity .15s}.chart-button:hover .expand-hint,.chart-button:focus .expand-hint{opacity:1}.empty{display:grid;place-items:center;min-height:210px;padding:28px;text-align:center;color:var(--muted);background:var(--panel)}.error-box{color:#9f3041}.footer{padding:24px 6px;color:var(--muted);font-size:12px;text-align:center}.print-button{position:fixed;right:22px;bottom:22px;border:0;border-radius:9px;padding:11px 16px;background:var(--brand);color:#fff;font-weight:700;box-shadow:0 8px 22px #0a132d35;cursor:pointer}dialog{width:min(96vw,1400px);max-height:94vh;padding:18px;border:0;border-radius:14px;box-shadow:0 20px 70px #0008}dialog::backdrop{background:#071127c9}.dialog-head{display:flex;justify-content:flex-end;margin-bottom:8px}.dialog-close{border:0;border-radius:7px;padding:8px 12px;background:var(--ink);color:#fff;cursor:pointer}#dialog-image{display:block;max-width:100%;max-height:82vh;margin:auto}
@media(max-width:980px){.layout{grid-template-columns:1fr}.sidebar{position:relative;height:auto}.nav{display:grid;grid-template-columns:repeat(2,minmax(0,1fr))}.main{padding:20px}.stats-grid,.plot-grid{grid-template-columns:1fr}}
@media(max-width:620px){.main{padding:0}.report-head,.module{border-radius:0;border-left:0;border-right:0}.report-head,.module{padding:22px}.overview,.nav{grid-template-columns:1fr}.sidebar{padding:22px}h1{font-size:25px}}
@media print{body{background:#fff}.layout{display:block;max-width:none}.sidebar,.print-button,dialog{display:none!important}.main{padding:0}.report-head,.module{box-shadow:none;border-radius:0;break-inside:avoid}.module{margin-top:12px}.chart-button{cursor:default}.expand-hint{display:none}.plot-grid{grid-template-columns:repeat(2,minmax(0,1fr))}}
""".strip()


_JS = """
const dialog=document.getElementById('chart-dialog');
const dialogImage=document.getElementById('dialog-image');
document.querySelectorAll('.chart-button').forEach((button)=>{
  button.addEventListener('click',()=>{
    const image=button.querySelector('img');
    dialogImage.src=image.src; dialogImage.alt=image.alt; dialog.showModal();
  });
});
document.getElementById('dialog-close').addEventListener('click',()=>dialog.close());
dialog.addEventListener('click',(event)=>{if(event.target===dialog)dialog.close();});
const reportSections=[...document.querySelectorAll('main section[id]')];
const navLinks=[...document.querySelectorAll('.nav a')];
const activateNav=(activeLink)=>navLinks.forEach((link)=>{
  const active=link===activeLink;
  link.classList.toggle('active',active);
  if(active)link.setAttribute('aria-current','true'); else link.removeAttribute('aria-current');
});
const updateActiveNav=()=>{
  if(!reportSections.length)return;
  const marker=24;
  let activeSection=reportSections[0];
  for(const section of reportSections){
    if(section.getBoundingClientRect().top<=marker)activeSection=section; else break;
  }
  if(window.scrollY+window.innerHeight>=document.documentElement.scrollHeight-2){
    activeSection=reportSections[reportSections.length-1];
  }
  activateNav(navLinks.find((link)=>link.hash===`#${activeSection.id}`));
};
navLinks.forEach((link)=>link.addEventListener('click',(event)=>{
  const target=document.querySelector(link.hash);
  if(!target)return;
  event.preventDefault();
  activateNav(link);
  const top=window.scrollY+target.getBoundingClientRect().top-20;
  window.scrollTo({top:Math.max(0,top),behavior:'auto'});
}));
let scrollFrame=0;
window.addEventListener('scroll',()=>{
  if(scrollFrame)return;
  scrollFrame=requestAnimationFrame(()=>{scrollFrame=0;updateActiveNav();});
},{passive:true});
updateActiveNav();
""".strip()


def render_qc_report(model: QcReportModel, plot_dir: Path) -> str:
    modules = "".join(_render_module(module, plot_dir) for module in model.modules)
    nav_items = [
        '<li><a href="#basic-statistics"><span class="nav-index">00</span>'
        '<span class="nav-label">Basic Statistics</span><i class="nav-dot info"></i></a></li>'
    ]
    nav_items.extend(
        f'<li><a href="#module-{escape(module.metric_id)}"><span class="nav-index">{index:02d}</span>'
        f'<span class="nav-label">{escape(module.title)}</span>'
        f'<i class="nav-dot {escape(module.availability)}"></i></a></li>'
        for index, module in enumerate(model.modules, start=1)
    )
    read_label = " / ".join(model.reads)
    error_label = str(model.errors) if model.errors else "None"
    return f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>NeoQC report — {escape(model.sample_id)}</title><style>{_CSS}</style></head>
<body><div class="layout"><aside class="sidebar"><div class="brand"><span class="brand-mark">NQ</span>NeoQC Report</div>
<div class="sidebar-meta"><small>Sample</small><strong>{escape(model.sample_id)}</strong><small>Generated</small><strong>{escape(model.generated_at)}</strong><small>Reads</small><strong>{escape(read_label)}</strong></div>
<p class="nav-title">Report sections</p><ul class="nav">{''.join(nav_items)}</ul></aside>
<main class="main"><header class="report-head"><p class="eyebrow">Sequencing quality control</p><h1>{escape(model.sample_id)}</h1>
<p class="subtitle">Compact, self-contained NeoQC report with all available read-quality charts.</p>
<div class="overview"><div class="overview-card"><span>Read sets</span><strong>{escape(read_label)}</strong></div><div class="overview-card"><span>Charts available</span><strong>{model.generated_plots}</strong></div><div class="overview-card"><span>Rendering errors</span><strong>{escape(error_label)}</strong></div></div>
<div class="notice">READY / NOT RUN / ERROR describe report artifact availability. Biological QC thresholds will be reported separately when the FastQC-compatible decision engine is implemented.</div></header>
{_render_basic_statistics(model)}{modules}<footer class="footer">Generated by NeoQC • Self-contained QC report</footer></main></div>
<button class="print-button" type="button" onclick="window.print()">Print / Save PDF</button>
<dialog id="chart-dialog"><div class="dialog-head"><button class="dialog-close" id="dialog-close" type="button">Close</button></div><img id="dialog-image" alt=""></dialog>
<script>{_JS}</script></body></html>"""


def generate_qc_report(
    result_dir: Path,
    plot_dir: Path | None = None,
    output_path: Path | None = None,
) -> Path:
    result_dir = result_dir.resolve()
    plot_dir = (plot_dir or result_dir / "plots").resolve()
    output_path = (output_path or result_dir / REPORT_FILENAME).resolve()
    model = load_report_model(result_dir, plot_dir)
    document = render_qc_report(model, plot_dir)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(document, encoding="utf-8")
    return output_path
