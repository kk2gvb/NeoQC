"""Self-contained NeoQC HTML report: data model and renderer.

The report is a single HTML file with embedded fonts, inline chart SVG, CSS and
JavaScript. It supports light/dark themes, English/Russian text, a triage or
sequential section order and foldable panels. See docs/html-report.md.
"""

from __future__ import annotations

import base64
import csv
import json
import math
import re
import xml.etree.ElementTree as ET
from collections import OrderedDict
from dataclasses import dataclass
from datetime import datetime, timedelta, timezone
from html import escape
from pathlib import Path
from typing import Mapping, Sequence

from neoqc_i18n import ru
from neoqc_theme import CHART_COLORS, LIGHT_TO_VAR


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
class OverrepresentedSequence:
    sequence: str
    count: int
    percentage: float


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
    qc_status: str = "not_evaluated"
    qc_reasons: tuple[str, ...] = ()
    qc_observations: tuple[tuple[str, float, str], ...] = ()
    overrepresented_sequences: tuple[OverrepresentedSequence, ...] | None = None
    # (locale, svg filename) pairs for translated chart variants, e.g. (("ru", "x_R1.ru.svg"),)
    localized_svg: tuple[tuple[str, str], ...] = ()


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

    @property
    def qc_status(self) -> str:
        severity = {"not_evaluated": -1, "pass": 0, "warning": 1, "fail": 2}
        statuses = [card.qc_status for card in self.cards]
        evaluated = [status for status in statuses if status != "not_evaluated"]
        return max(evaluated, key=severity.__getitem__) if evaluated else "not_evaluated"


@dataclass(frozen=True)
class QcReportModel:
    sample_id: str
    generated_at: str
    reads: tuple[str, ...]
    basic_statistics: tuple[BasicStatistics, ...]
    modules: tuple[PlotModule, ...]
    generated_plots: int
    errors: int
    qc_counts: tuple[tuple[str, int], ...]
    overall_qc_status: str
    ruleset_label: str


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


def _load_overrepresented_sequences(
    result_dir: Path, read: str
) -> tuple[OverrepresentedSequence, ...] | None:
    path = result_dir / f"overrepresented_sequences_{read}.tsv"
    if not path.is_file():
        return None
    expected_header = ["sequence", "count", "percentage"]
    try:
        with path.open(encoding="utf-8", newline="") as stream:
            reader = csv.DictReader(stream, delimiter="\t")
            if reader.fieldnames != expected_header:
                raise QcReportError(
                    f"{path.name} must use columns: " + ", ".join(expected_header)
                )
            rows: list[OverrepresentedSequence] = []
            seen_sequences: set[str] = set()
            for line_number, row in enumerate(reader, start=2):
                if None in row or any(row[column] is None for column in expected_header):
                    raise QcReportError(f"{path.name}:{line_number} has an invalid column count")
                sequence = row["sequence"].strip()
                if not sequence:
                    raise QcReportError(f"{path.name}:{line_number} has an empty sequence")
                if sequence in seen_sequences:
                    raise QcReportError(
                        f"{path.name}:{line_number} repeats sequence {sequence!r}"
                    )
                try:
                    count = int(row["count"])
                    percentage = float(row["percentage"])
                except ValueError as error:
                    raise QcReportError(
                        f"{path.name}:{line_number} contains a non-numeric value"
                    ) from error
                if count <= 0:
                    raise QcReportError(f"{path.name}:{line_number} count must be positive")
                if not math.isfinite(percentage) or not 0.0 <= percentage <= 100.0:
                    raise QcReportError(
                        f"{path.name}:{line_number} percentage must be between 0 and 100"
                    )
                seen_sequences.add(sequence)
                rows.append(
                    OverrepresentedSequence(
                        sequence=sequence,
                        count=count,
                        percentage=percentage,
                    )
                )
    except OSError as error:
        raise QcReportError(f"cannot read {path.name}: {error}") from error
    return tuple(rows)


def _load_qc_evaluations(
    result_dir: Path,
) -> tuple[
    dict[tuple[str, str], tuple[str, tuple[str, ...], tuple[tuple[str, float, str], ...]]],
    str,
]:
    path = result_dir / "qc_evaluation.json"
    if not path.is_file():
        return {}, "No QC ruleset evaluation"
    try:
        raw = json.loads(path.read_text(encoding="utf-8"))
    except OSError as error:
        raise QcReportError(f"cannot read {path}: {error}") from error
    except json.JSONDecodeError as error:
        raise QcReportError(f"invalid JSON in {path}: {error}") from error
    document = _mapping(raw, "qc_evaluation")
    if document.get("schema_version") != 1:
        raise QcReportError("unsupported QC evaluation schema")
    ruleset = _mapping(document.get("ruleset"), "qc_evaluation.ruleset")
    ruleset_id = _text(ruleset.get("id"), "qc_evaluation.ruleset.id", required=True)
    ruleset_version = _text(
        ruleset.get("version"), "qc_evaluation.ruleset.version", required=True
    )
    decisions: dict[
        tuple[str, str], tuple[str, tuple[str, ...], tuple[tuple[str, float, str], ...]]
    ] = {}
    for index, raw_evaluation in enumerate(
        _sequence(document.get("evaluations"), "qc_evaluation.evaluations")
    ):
        evaluation = _mapping(raw_evaluation, f"qc_evaluation.evaluations[{index}]")
        metric_id = _text(
            evaluation.get("metric_id"),
            f"qc_evaluation.evaluations[{index}].metric_id",
            required=True,
        )
        read = _text(
            evaluation.get("read"),
            f"qc_evaluation.evaluations[{index}].read",
            required=True,
        )
        status = _text(
            evaluation.get("qc_status"),
            f"qc_evaluation.evaluations[{index}].qc_status",
            required=True,
        )
        if read not in {"R1", "R2"} or status not in {
            "pass",
            "warning",
            "fail",
            "not_evaluated",
        }:
            raise QcReportError(f"qc_evaluation.evaluations[{index}] is invalid")
        key = (metric_id, read)
        if key in decisions:
            raise QcReportError(f"duplicate QC evaluation for {metric_id}/{read}")
        raw_reasons = _sequence(
            evaluation.get("reasons"), f"qc_evaluation.evaluations[{index}].reasons"
        )
        reasons = tuple(
            _text(
                _mapping(reason, f"qc_evaluation.evaluations[{index}].reasons[{reason_index}]").get("message"),
                f"qc_evaluation.evaluations[{index}].reasons[{reason_index}].message",
                required=True,
            )
            for reason_index, reason in enumerate(raw_reasons)
        )
        raw_observations = _mapping(
            evaluation.get("observations"),
            f"qc_evaluation.evaluations[{index}].observations",
        )
        check_labels: dict[str, tuple[str, str]] = {}
        for check_index, raw_check in enumerate(
            _sequence(evaluation.get("checks"), f"qc_evaluation.evaluations[{index}].checks")
        ):
            check = _mapping(
                raw_check, f"qc_evaluation.evaluations[{index}].checks[{check_index}]"
            )
            observation = _text(
                check.get("observation"),
                f"qc_evaluation.evaluations[{index}].checks[{check_index}].observation",
                required=True,
            )
            check_labels[observation] = (
                _text(check.get("label"), f"check[{check_index}].label", required=True),
                _text(check.get("unit"), f"check[{check_index}].unit"),
            )
        observations: list[tuple[str, float, str]] = []
        for name, value in raw_observations.items():
            if isinstance(value, bool) or not isinstance(value, (int, float)):
                raise QcReportError(f"observation {name} must be numeric")
            if name not in check_labels:
                continue
            label, unit = check_labels[name]
            observations.append((label, float(value), unit))
        decisions[key] = (status, reasons, tuple(observations))
    return decisions, f"{ruleset_id} · v{ruleset_version}"


SUPPORTED_LOCALES = ("en", "ru")


def _localized_svg(value: object, path: str) -> tuple[tuple[str, str], ...]:
    if value is None:
        return ()
    localized = _mapping(value, path)
    variants: list[tuple[str, str]] = []
    for locale, raw_variant in localized.items():
        if locale not in SUPPORTED_LOCALES or locale == "en":
            raise QcReportError(f"{path} has an unsupported locale: {locale!r}")
        variant = _mapping(raw_variant, f"{path}.{locale}")
        variants.append((locale, _text(variant.get("svg"), f"{path}.{locale}.svg", required=True)))
    return tuple(sorted(variants))


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

    decisions, ruleset_label = _load_qc_evaluations(result_dir)
    raw_plots = _sequence(manifest.get("plots"), "manifest.plots")
    cards: list[PlotCard] = []
    overrepresented_by_read: dict[
        str, tuple[OverrepresentedSequence, ...] | None
    ] = {}
    for index, raw_entry in enumerate(raw_plots):
        entry = _mapping(raw_entry, f"manifest.plots[{index}]")
        status = _text(entry.get("status"), f"manifest.plots[{index}].status", required=True)
        if status not in {"generated", "skipped", "error"}:
            raise QcReportError(f"manifest.plots[{index}].status is invalid")
        read = _text(entry.get("read"), f"manifest.plots[{index}].read", required=True)
        if read not in {"R1", "R2"}:
            raise QcReportError(f"manifest.plots[{index}].read must be R1 or R2")
        metric_id = _text(entry.get("id"), f"manifest.plots[{index}].id", required=True)
        decision = decisions.get((metric_id, read))
        if metric_id == "sequence_duplication_levels" and read not in overrepresented_by_read:
            overrepresented_by_read[read] = _load_overrepresented_sequences(result_dir, read)
        cards.append(
            PlotCard(
                metric_id=metric_id,
                read=read,
                title=_text(entry.get("title"), f"manifest.plots[{index}].title", required=True),
                status=status,
                reason=_text(entry.get("reason"), f"manifest.plots[{index}].reason"),
                alt_text=_text(entry.get("alt_text"), f"manifest.plots[{index}].alt_text"),
                svg=_text(entry.get("svg"), f"manifest.plots[{index}].svg"),
                png=_text(entry.get("png"), f"manifest.plots[{index}].png"),
                qc_status=decision[0] if decision else "not_evaluated",
                qc_reasons=decision[1] if decision else ("QC evaluation is not available.",),
                qc_observations=decision[2] if decision else (),
                overrepresented_sequences=(
                    overrepresented_by_read[read]
                    if metric_id == "sequence_duplication_levels"
                    else None
                ),
                localized_svg=_localized_svg(entry.get("localized"), f"manifest.plots[{index}].localized"),
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
    counts = {status: 0 for status in ("pass", "warning", "fail", "not_evaluated")}
    for card in cards:
        if card.read in reads:
            counts[card.qc_status] += 1
    severity = {"not_evaluated": -1, "pass": 0, "warning": 1, "fail": 2}
    evaluated = [status for status in ("pass", "warning", "fail") if counts[status]]
    overall = max(evaluated, key=severity.__getitem__) if evaluated else "not_evaluated"
    return QcReportModel(
        sample_id=sample_id,
        generated_at=datetime.now(MOSCOW_TIME).strftime(
            "%d.%m.%Y, %H:%M GMT+3"
        ),
        reads=reads,
        basic_statistics=summaries,
        modules=modules,
        generated_plots=generated_plots,
        errors=errors,
        qc_counts=tuple(counts.items()),
        overall_qc_status=overall,
        ruleset_label=ruleset_label,
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


STATUS_LABELS = {
    "pass": "PASS",
    "warning": "WARNING",
    "fail": "FAIL",
    "not_evaluated": "NOT EVALUATED",
}
STATUS_SYMBOLS = {"pass": "✓", "warning": "▲", "fail": "✕", "not_evaluated": "—"}
STATUS_ORDER = ("pass", "warning", "fail", "not_evaluated")
SKIP_REASONS = {
    "adapter_analysis_disabled": "Adapter analysis was disabled for this run.",
    "source_not_found": "The source TSV was not produced.",
}
NOTICE = (
    "PASS / WARNING / FAIL are technical QC flags from the displayed versioned ruleset, "
    "not clinical conclusions. Plot availability errors are reported separately and never "
    "converted into biological FAIL results."
)
FONT_DIR = Path(__file__).resolve().parents[1] / "assets" / "fonts" / "ibm-plex"
FONT_FACES = (
    ("IBM Plex Sans", "IBMPlexSans", (400, 500, 600, 700)),
    ("IBM Plex Mono", "IBMPlexMono", (400, 500, 600)),
)


# ---------------------------------------------------------------------------
# Bilingual text
# ---------------------------------------------------------------------------
def _t(text: str) -> str:
    """Escaped bilingual text; CSS shows the span that matches <html lang>."""
    translated = ru(text)
    if translated == text:
        return escape(text)
    return f'<span class="l-en">{escape(text)}</span><span class="l-ru">{escape(translated)}</span>'


# ---------------------------------------------------------------------------
# Inline chart SVG
# ---------------------------------------------------------------------------
_SVG_NS = "http://www.w3.org/2000/svg"
_XLINK_NS = "http://www.w3.org/1999/xlink"
ET.register_namespace("", _SVG_NS)
ET.register_namespace("xlink", _XLINK_NS)
# Everything Matplotlib emits for NeoQC charts; anything else (script,
# foreignObject, image, a, style, metadata, ...) is dropped before inlining.
_SVG_TAGS = frozenset({
    "svg", "g", "defs", "path", "rect", "circle", "ellipse", "line", "polyline", "polygon",
    "text", "tspan", "clipPath", "use", "title", "desc", "linearGradient", "radialGradient",
    "stop", "pattern", "mask", "symbol", "marker",
})
_COLOR_ATTRIBUTES = ("fill", "stroke", "stop-color", "color")
_HEX_COLOR = re.compile(r"#[0-9a-fA-F]{6}\b")
_URL_REF = re.compile(r"url\(\s*#([^)\s]+)\s*\)")
_UNSAFE_STYLE = re.compile(r"url\(\s*(?!#)|expression\s*\(|@import|javascript:", re.IGNORECASE)
_FONT_FALLBACK = re.compile(r"(font(?:-family)?:[^;\"]*?)'(IBM Plex Sans|IBM Plex Mono)'(?!\s*,)")


def _local_name(name: str) -> str:
    return name.rsplit("}", 1)[-1]


def _themed(value: str) -> str:
    return _HEX_COLOR.sub(
        lambda match: f"var({LIGHT_TO_VAR[match[0].lower()]})" if match[0].lower() in LIGHT_TO_VAR else match[0],
        value,
    )


def _with_font_fallback(style: str) -> str:
    def fallback(match: re.Match[str]) -> str:
        generic = "Consolas, monospace" if match[2] == "IBM Plex Mono" else "'Segoe UI', Arial, sans-serif"
        return f"{match[1]}'{match[2]}', {generic}"

    return _FONT_FALLBACK.sub(fallback, style)


def _sanitize_svg_element(element: ET.Element, prefix: str) -> None:
    for child in list(element):
        if not isinstance(child.tag, str) or _local_name(child.tag) not in _SVG_TAGS:
            element.remove(child)
        else:
            _sanitize_svg_element(child, prefix)
    style_parts: list[str] = []
    for name, value in list(element.attrib.items()):
        local = _local_name(name)
        lowered = value.strip().lower()
        if local.lower().startswith("on") or "javascript:" in lowered:
            del element.attrib[name]
        elif local == "href":
            if value.startswith("#"):
                element.attrib[name] = f"#{prefix}{value[1:]}"
            else:  # external references are never followed
                del element.attrib[name]
        elif local == "id":
            element.attrib[name] = f"{prefix}{value}"
        elif local == "style":
            if _UNSAFE_STYLE.search(value):
                del element.attrib[name]
            else:
                element.attrib[name] = _with_font_fallback(_themed(_URL_REF.sub(rf"url(#{prefix}\1)", value)))
        elif local in _COLOR_ATTRIBUTES and _HEX_COLOR.search(value):
            # CSS variables are not valid in presentation attributes; move the colour into style.
            del element.attrib[name]
            style_parts.append(f"{local}: {_themed(value)}")
        elif "url(" in lowered:
            if _UNSAFE_STYLE.search(value):
                del element.attrib[name]
            else:
                element.attrib[name] = _URL_REF.sub(rf"url(#{prefix}\1)", value)
    if style_parts:
        existing = element.attrib.get("style", "").strip().rstrip(";")
        element.attrib["style"] = "; ".join(filter(None, [existing, *style_parts]))


def _inline_svg(path: Path, prefix: str, label: str) -> str:
    """Sanitize a chart SVG and return inline markup whose colours follow the report theme."""
    text = path.read_text(encoding="utf-8")
    # Matplotlib writes a plain SVG 1.1 DOCTYPE. Entity declarations or an
    # internal DTD subset could expand content, so they are rejected.
    if re.search(r"<!ENTITY|<!DOCTYPE[^>]*\[", text, re.IGNORECASE):
        raise QcReportError(f"{path.name}: SVG with entity declarations is not accepted")
    text = re.sub(r"<!DOCTYPE[^>]*>", "", text, count=1, flags=re.IGNORECASE)
    try:
        root = ET.fromstring(text)
    except ET.ParseError as error:
        raise QcReportError(f"{path.name}: invalid SVG ({error})") from error
    if _local_name(root.tag) != "svg":
        raise QcReportError(f"{path.name}: root element is not <svg>")
    _sanitize_svg_element(root, prefix)
    if "viewBox" in root.attrib:  # let CSS size the chart
        root.attrib.pop("width", None)
        root.attrib.pop("height", None)
    root.attrib.update({"class": "chart-svg", "role": "img", "aria-label": label, "focusable": "false"})
    return ET.tostring(root, encoding="unicode", short_empty_elements=True)


def _chart_markup(card: PlotCard, plot_dir: Path) -> str:
    label = card.alt_text or f"{card.title} {card.read}"
    if card.status != "generated":
        message = SKIP_REASONS.get(card.reason, card.reason or "Chart was not generated.")
        return f'<div class="empty">{_t(message)}</div>'
    try:
        if card.svg:
            prefix = f"{card.metric_id}-{card.read}-"
            primary = _inline_svg(_safe_asset(plot_dir, card.svg, ".svg"), f"{prefix}en-", label)
            translated = []
            for locale, filename in card.localized_svg:
                try:
                    svg = _inline_svg(_safe_asset(plot_dir, filename, ".svg"), f"{prefix}{locale}-", label)
                except (OSError, QcReportError):
                    continue  # a missing translation falls back to the English chart
                translated.append(f'<span class="chart l-{locale}">{svg}</span>')
            if translated:
                inner = f'<span class="chart l-en">{primary}</span>' + "".join(translated)
            else:
                inner = f'<span class="chart">{primary}</span>'
        else:  # PNG-only runs: a static image that does not follow the dark theme
            inner = f'<img class="chart-img" src="{_asset_data_uri(plot_dir, card)}" alt="{escape(label)}">'
    except (OSError, QcReportError) as error:
        return (f'<div class="empty error-box">{_t("Chart asset unavailable")}: '
                f'{escape(str(error))}</div>')
    return (f'<button class="chart-button" type="button" aria-label="{escape(label)} — open full size">'
            f"{inner}</button>")


def _font_face_css() -> str:
    """Embed the bundled IBM Plex web fonts (latin + cyrillic); missing files fall back to system fonts."""
    try:
        ranges = json.loads((FONT_DIR / "unicode_ranges.json").read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return ""
    rules = []
    for family, slug, weights in FONT_FACES:
        for weight in weights:
            for subset, unicode_range in ranges.items():
                path = FONT_DIR / f"{slug}-{weight}.{subset}.woff2"
                if not path.is_file():
                    continue
                data = base64.b64encode(path.read_bytes()).decode("ascii")
                rules.append(
                    f"@font-face{{font-family:'{family}';font-style:normal;font-weight:{weight};"
                    f"font-display:swap;src:url(data:font/woff2;base64,{data}) format('woff2');"
                    f"unicode-range:{unicode_range}}}"
                )
    return "\n".join(rules)


# ---------------------------------------------------------------------------
# Markup helpers
# ---------------------------------------------------------------------------
def _pill(status: str, label: str | None = None) -> str:
    symbol = STATUS_SYMBOLS.get(status)
    text = escape(label or STATUS_LABELS.get(status, status.upper()))
    mark = f'<span aria-hidden="true">{symbol}</span> ' if symbol else ""
    return f'<span class="pill {escape(status)}">{mark}{text}</span>'


def _cell(status: str, title: str) -> str:
    return f'<i class="cell {escape(status)}" title="{escape(title)}"></i>'


def _keys(seq: int, tri: int) -> str:
    return f'data-seq="{seq}" data-tri="{tri}"'


def _triage_rank(status: str, index: int) -> int:
    """FAIL, WARNING and NOT EVALUATED need attention; PASS modules follow folded."""
    return {"fail": 100, "warning": 200, "not_evaluated": 300}.get(status, 6000) + index


def _needs_attention(status: str) -> bool:
    return status in ("fail", "warning", "not_evaluated")


CHEVRON = (
    '<button class="fold" type="button" aria-label="Collapse or expand section" aria-expanded="true">'
    '<svg width="12" height="12" viewBox="0 0 12 12" aria-hidden="true"><path d="M2.5 4.5 6 8l3.5-3.5" '
    'fill="none" stroke="currentColor" stroke-width="1.6" stroke-linecap="round" stroke-linejoin="round"/>'
    "</svg></button>"
)
SUN = (
    '<svg width="13" height="13" viewBox="0 0 16 16" aria-hidden="true"><circle cx="8" cy="8" r="3.2" '
    'fill="none" stroke="currentColor" stroke-width="1.5"/><path d="M8 1v2M8 13v2M1 8h2M13 8h2M3 3l1.4 '
    '1.4M11.6 11.6 13 13M3 13l1.4-1.4M11.6 4.4 13 3" stroke="currentColor" stroke-width="1.5" '
    'stroke-linecap="round"/></svg>'
)
MOON = (
    '<svg width="13" height="13" viewBox="0 0 16 16" aria-hidden="true"><path d="M13.5 10.2A6 6 0 0 1 '
    '5.8 2.5a6 6 0 1 0 7.7 7.7Z" fill="none" stroke="currentColor" stroke-width="1.5" '
    'stroke-linejoin="round"/></svg>'
)


def _render_overview(model: QcReportModel, keys: str) -> str:
    counts = dict(model.qc_counts)
    total = sum(counts.values()) or 1
    counters = "".join(
        f'<div class="count {status}"><small>{STATUS_LABELS[status]}</small><strong>{counts[status]}</strong></div>'
        for status in STATUS_ORDER
    )
    segments = "".join(
        f'<i class="{status}" style="width:{counts[status] / total * 100:.4f}%" '
        f'title="{STATUS_LABELS[status]}: {counts[status]}"></i>'
        for status in STATUS_ORDER
        if counts[status]
    )
    bar_label = ", ".join(f"{STATUS_LABELS[status]} {counts[status]}" for status in STATUS_ORDER)
    rows = []
    for index, module in enumerate(model.modules, start=2):
        by_read = {card.read: card for card in module.cards}
        cells = "".join(
            f"<td>{_cell(by_read[read].qc_status, STATUS_LABELS[by_read[read].qc_status])}</td>"
            if read in by_read else "<td>—</td>"
            for read in model.reads
        )
        rows.append(
            f'<tr {_keys(index, _triage_rank(module.qc_status, index))}><th scope="row" class="row-h">'
            f"{_t(module.title)}</th>{cells}</tr>"
        )
    headers = "".join(f'<th scope="col" class="c">{escape(read)}</th>' for read in model.reads)
    error_label = str(model.errors) if model.errors else "None"
    facts = (
        ("Read sets", " / ".join(model.reads)),
        ("Charts available", str(model.generated_plots)),
        ("Rendering errors", error_label),
    )
    facts_html = "".join(f"<div><small>{_t(k)}</small><strong>{_t(v)}</strong></div>" for k, v in facts)
    return (
        f'<section class="panel" id="qc-summary" data-section {keys}><div class="hero"><div class="id">'
        f'<div class="eyebrow">{_t("Sequencing quality control")}</div><h1>{escape(model.sample_id)}</h1>'
        f'<div class="kv">{facts_html}</div>'
        f'<p class="print-meta">{_t("Generated")}: {escape(model.generated_at)}</p>'
        f'<div class="notice">{_t(NOTICE)}</div></div>'
        f'<div class="qcbox"><div class="top"><div><h2>{_t("Technical QC overview · PASS / WARNING / FAIL distribution")}</h2>'
        f'<div class="rs">{_t("Ruleset")}: {escape(model.ruleset_label)}</div></div>{_pill(model.overall_qc_status)}</div>'
        f'<div class="counts">{counters}</div><div class="sbar" role="img" aria-label="{bar_label}">{segments}</div>'
        f'</div></div><div class="table-wrap"><table class="grid hm"><thead><tr><th scope="col">{_t("QC module")}</th>'
        f'{headers}</tr></thead><tbody>{"".join(rows)}</tbody></table></div></section>'
    )


def _render_basic_statistics(model: QcReportModel, keys: str) -> str:
    if model.basic_statistics:
        labels: list[str] = []
        for stats in model.basic_statistics:
            labels.extend(label for label, _ in stats.items if label not in labels)
        values = {stats.read: dict(stats.items) for stats in model.basic_statistics}
        head = "".join(f'<th scope="col" class="n">{escape(stats.read)}</th>' for stats in model.basic_statistics)
        rows = "".join(
            f'<tr><th scope="row" class="row-h">{_t(label)}</th>' + "".join(
                f'<td class="n">{escape(values[stats.read].get(label, "—"))}</td>'
                for stats in model.basic_statistics
            ) + "</tr>"
            for label in labels
        )
        body = (f'<div class="table-wrap"><table class="grid"><thead><tr><th scope="col">{_t("Metric")}</th>{head}'
                f"</tr></thead><tbody>{rows}</tbody></table></div>")
        first = labels[0] if labels else ""
        peek = "".join(
            f"<span>{escape(stats.read)} <b>{escape(values[stats.read].get(first, '—'))}</b></span>"
            for stats in model.basic_statistics
        )
    else:
        body = f'<div class="empty">{_t("Summary files were not found.")}</div>'
        peek = ""
    return (
        f'<section class="panel foldable" id="basic-statistics" data-section data-sev="info" {keys}>'
        f'<div class="phead" data-fold>{CHEVRON}<h2>{_t("Basic Statistics")}</h2><span class="peek">{peek}</span>'
        f'<span class="sp"></span>{_pill("info", "INFO")}</div>{body}</section>'
    )


def _render_overrepresented_sequences(card: PlotCard) -> str:
    rows = card.overrepresented_sequences
    if rows is None:
        return ""
    if rows:
        body = "".join(
            f'<tr><td><code>{escape(row.sequence)}</code></td><td class="n">{row.count:,}</td>'
            f'<td class="n">{row.percentage:.4f}%</td></tr>'
            for row in rows
        )
        table = (f'<div class="table-wrap"><table class="grid"><thead><tr><th scope="col">{_t("Sequence")}</th>'
                 f'<th scope="col" class="n">{_t("Count")}</th><th scope="col" class="n">{_t("Percentage")}</th>'
                 f"</tr></thead><tbody>{body}</tbody></table></div>")
    else:
        table = f'<p class="none">{_t("No sequences exceeded the reporting threshold.")}</p>'
    return (
        f'<div class="orep" aria-label="Overrepresented sequences"><div class="orep-h">'
        f'<span>{_t("Sequence screen · Overrepresented sequences")}</span><span class="table-count">{len(rows)}</span></div>{table}</div>'
    )


def _render_lane(card: PlotCard, plot_dir: Path) -> str:
    observations = "".join(
        f"<tr><td>{_t(label)}</td><td>{escape(f'{value:.4g}' + (' ' + unit if unit else ''))}</td></tr>"
        for label, value, unit in card.qc_observations
    )
    reasons = "".join(f'<div class="why">{_t(reason)}</div>' for reason in card.qc_reasons)
    if card.status == "generated":
        state = _pill(card.qc_status)
    else:
        state = _pill("not_evaluated", "ERROR" if card.status == "error" else "NOT RUN")
    return (
        f'<div class="lane"><div class="lane-h"><span>{escape(card.read)}</span>{state}</div>'
        f"{_chart_markup(card, plot_dir)}"
        f'<div class="metrics">{"<table>" + observations + "</table>" if observations else ""}{reasons}</div>'
        f"{_render_overrepresented_sequences(card)}</div>"
    )


def _render_module(module: PlotModule, index: int, plot_dir: Path) -> str:
    peek = "".join(
        f"<span>{escape(card.read)} <b>{escape(f'{card.qc_observations[0][1]:.4g}' + (' ' + card.qc_observations[0][2] if card.qc_observations[0][2] else ''))}</b></span>"
        for card in module.cards
        if card.qc_observations
    )
    severity = module.qc_status if _needs_attention(module.qc_status) else "ok"
    lanes = "".join(_render_lane(card, plot_dir) for card in module.cards)
    return (
        f'<section class="panel foldable" id="module-{escape(module.metric_id)}" data-section data-sev="{severity}" '
        f'{_keys(index, _triage_rank(module.qc_status, index))}>'
        f'<div class="phead" data-fold>{CHEVRON}<h2>{_t(module.title)}</h2><span class="peek">{peek}</span>'
        f'<span class="sp"></span>{_pill(module.qc_status)}</div><div class="tracks">{lanes}</div></section>'
    )


# ---------------------------------------------------------------------------
# Styles and behaviour
# ---------------------------------------------------------------------------
_LIGHT_VARS = (
    "color-scheme:light;--bg:#eceff3;--panel:#fff;--sub:#f6f7f9;--line:#d5dbe3;--line2:#e6eaef;--ink:#1b2430;"
    "--muted:#5d6878;--faint:#626d7d;--brand:#1f5fd1;--pass:#0e9f6e;--warning:#e09b00;--fail:#d8343d;--ne:#9aa3b0;"
    "--info:#1f5fd1;--pill-pass:#0a7f57;--pill-warning:#e09b00;--pill-fail:#d8343d;--pill-ne:#6b7686;"
    "--pill-info:#1f5fd1;--on-pill:#fff;--on-pill-warning:#1b2430;--fail-text:#b52a33;--hover:#e9edf2;"
    "--active:#e3ebfb;--row-hover:#f7f9fc;--head-hover:#eef1f5;--fold-hover:#dfe4ea;--notice-bg:#f3f7ff;"
    "--notice-line:#d6e2fb;--notice-ink:#34425a;--backdrop:rgba(20,28,40,.5);"
    + ";".join(f"{name}:{light}" for name, (light, _) in CHART_COLORS.items())
)
_DARK_VARS = (
    "color-scheme:dark;--bg:#0e1319;--panel:#151b23;--sub:#1a212b;--line:#2b3441;--line2:#222a35;--ink:#e3e8ef;"
    "--muted:#9aa6b6;--faint:#8591a3;--brand:#6c9bff;--pass:#1fae7a;--warning:#e0a52a;--fail:#e5535b;--ne:#5f6a7a;"
    "--info:#5b8cff;--pill-pass:#1fae7a;--pill-warning:#e0a52a;--pill-fail:#e5535b;--pill-ne:#8591a3;"
    "--pill-info:#5b8cff;--on-pill:#0e1319;--on-pill-warning:#0e1319;--fail-text:#ff7b82;--hover:#1f2833;"
    "--active:#1c2a47;--row-hover:#1a2230;--head-hover:#1d2530;--fold-hover:#2a3340;--notice-bg:#16213a;"
    "--notice-line:#26375e;--notice-ink:#b8c6e0;--backdrop:rgba(0,0,0,.62);"
    + ";".join(f"{name}:{dark}" for name, (_, dark) in CHART_COLORS.items())
)

_CSS = (
    ":root{" + _LIGHT_VARS + ";--sans:'IBM Plex Sans','Segoe UI',Arial,sans-serif;"
    "--mono:'IBM Plex Mono',Consolas,monospace}\n"
    ':root[data-theme="dark"]{' + _DARK_VARS + "}\n"
    "@media print{:root,:root[data-theme=\"dark\"]{" + _LIGHT_VARS + "}}\n"
    + """
html[lang="ru"] .l-en,html:not([lang="ru"]) .l-ru{display:none!important}
*{box-sizing:border-box}html{scroll-padding-top:62px}
body{margin:0;background:var(--bg);color:var(--ink);font:13px/1.5 var(--sans)}
:focus-visible{outline:2px solid var(--brand);outline-offset:2px}
.appbar{position:sticky;top:0;z-index:10;display:flex;align-items:stretch;height:46px;background:var(--panel);border-bottom:1px solid var(--line)}
.logo{display:flex;align-items:center;gap:8px;padding:0 16px;border-right:1px solid var(--line);font:600 14px var(--sans);width:268px;flex:none}
.logo i{display:grid;grid-template-columns:repeat(2,7px);gap:2px}.logo i b{width:7px;height:7px;display:block}
.logo i b:nth-child(1){background:var(--c-a)}.logo i b:nth-child(2){background:var(--c-c)}.logo i b:nth-child(3){background:var(--c-warning)}.logo i b:nth-child(4){background:var(--c-danger)}
.crumbs{display:flex;align-items:center;flex:1;min-width:0;overflow:hidden}
.crumbs .cr{display:flex;flex-direction:column;justify-content:center;padding:0 16px;height:100%;border-right:1px solid var(--line2);white-space:nowrap}
.crumbs small{font:500 10px/1.2 var(--mono);color:var(--faint);text-transform:uppercase;letter-spacing:.04em}
.crumbs strong{font:500 12.5px/1.3 var(--mono)}
.btn{border:1px solid var(--line);background:var(--panel);color:var(--ink);padding:0 12px;font:500 12px var(--sans);cursor:pointer;border-radius:3px;white-space:nowrap}
.btn:hover{border-color:var(--brand);color:var(--brand)}
.appbar .btn{margin:8px 12px 8px 0}.appbar .btn.ghost{margin-right:8px}
.seg{display:flex;align-items:stretch;border:1px solid var(--line);border-radius:3px;overflow:hidden;flex:none}
.appbar .seg{margin:8px 8px 8px auto}
.seg small{display:flex;align-items:center;padding:0 9px;font:500 10px var(--mono);color:var(--faint);text-transform:uppercase;letter-spacing:.06em;border-right:1px solid var(--line);background:var(--sub)}
.seg button{border:0;background:var(--panel);padding:0 12px;font:500 12px var(--sans);color:var(--muted);cursor:pointer;display:flex;align-items:center;gap:6px}
.seg button+button{border-left:1px solid var(--line)}
.seg button[aria-pressed="true"]{background:var(--ink);color:var(--panel)}
.shell{display:grid;grid-template-columns:268px minmax(0,1fr)}
.rail{position:sticky;top:46px;height:calc(100vh - 46px);overflow:auto;background:var(--sub);border-right:1px solid var(--line);padding:12px 0}
.prefs{display:grid;grid-template-columns:1fr 1fr;gap:10px;margin:0 16px 12px;padding-bottom:12px;border-bottom:1px solid var(--line)}
.pref small{display:block;margin-bottom:5px;font:600 10px/1 var(--mono);color:var(--faint);text-transform:uppercase;letter-spacing:.08em}
.pref .seg{height:28px}.pref .seg button{flex:1;justify-content:center;padding:0 6px;font:600 11.5px var(--mono)}
.pref svg{display:block}
.rail-h{width:calc(100% - 32px);border:0;background:none;padding:0;cursor:default;margin:8px 16px 6px;font:600 10px/1 var(--mono);color:var(--faint);text-transform:uppercase;letter-spacing:.08em;display:flex;justify-content:space-between}
.rail-h em{font-style:normal;display:flex;gap:4px}.rail-h em b{width:18px;text-align:center;font-weight:600}
.track{display:flex;align-items:center;gap:8px;padding:6px 16px;text-decoration:none;color:var(--ink);border-left:3px solid transparent}
.track:hover{background:var(--hover)}.track.active{background:var(--active);border-left-color:var(--brand)}
.track .i{font:500 11px var(--mono);color:var(--faint);width:18px;flex:none}.track .l{flex:1;min-width:0;overflow:hidden;text-overflow:ellipsis;white-space:nowrap}
.track .cells{display:flex;gap:4px}
.cell{width:18px;height:14px;display:grid;place-items:center;border-radius:2px;font:700 9px/1 var(--sans);font-style:normal;color:var(--on-pill)}
.cell.pass{background:var(--pass)}.cell.warning{background:var(--warning);color:var(--on-pill-warning)}.cell.fail{background:var(--fail)}.cell.not_evaluated{background:var(--ne)}
.cell.info{background:var(--panel);border:1px solid var(--line)}
.cell.pass::after{content:"✓"}.cell.warning::after{content:"▲";font-size:7px}.cell.fail::after{content:"✕"}.cell.not_evaluated::after{content:"—"}
.legend{margin:14px 16px 0;padding-top:12px;border-top:1px solid var(--line);display:grid;grid-template-columns:1fr 1fr;gap:6px;font:11px var(--mono);color:var(--muted)}
.legend span{display:flex;align-items:center;gap:6px;white-space:nowrap;font-size:10.5px}
.main{padding:16px;display:grid;grid-template-columns:minmax(0,1fr);gap:12px;min-width:0}
.panel{background:var(--panel);border:1px solid var(--line);border-radius:3px;min-width:0}
.phead{display:flex;align-items:center;gap:10px;padding:8px 12px;border-bottom:1px solid var(--line);background:var(--sub)}
.phead h2{margin:0;font:600 13.5px var(--sans)}.phead .sp{flex:1}
.pill{display:inline-flex;align-items:center;gap:4px;padding:1px 7px;border-radius:2px;font:600 10.5px/1.6 var(--mono);letter-spacing:.03em;color:var(--on-pill);white-space:nowrap}
.pill.pass{background:var(--pill-pass)}.pill.warning{background:var(--pill-warning);color:var(--on-pill-warning)}.pill.fail{background:var(--pill-fail)}
.pill.not_evaluated{background:var(--pill-ne)}.pill.info{background:var(--pill-info)}
.hero{display:grid;grid-template-columns:minmax(0,1.1fr) minmax(0,1fr)}
.hero .id{padding:16px 16px 14px;border-right:1px solid var(--line);min-width:0}
.eyebrow{font:500 11px var(--mono);color:var(--brand);text-transform:uppercase;letter-spacing:.06em}
.hero h1{margin:4px 0 2px;font:600 30px/1.1 var(--mono);letter-spacing:-.02em;overflow-wrap:anywhere}
.kv{display:grid;grid-template-columns:repeat(3,1fr);border-top:1px solid var(--line2);margin-top:14px}
.kv div{padding:8px 0 0}.kv small{display:block;font:500 10px var(--mono);color:var(--faint);text-transform:uppercase}.kv strong{font:500 15px var(--mono)}
.print-meta{display:none;margin:10px 0 0;font:12px var(--mono);color:var(--muted)}
.notice{margin:12px 0 0;padding:8px 10px;background:var(--notice-bg);border:1px solid var(--notice-line);border-left:3px solid var(--brand);font-size:12px;color:var(--notice-ink)}
.qcbox{padding:14px 16px;min-width:0}
.qcbox .top{display:flex;justify-content:space-between;align-items:center;gap:10px}
.qcbox h2{margin:0;font:600 13px var(--sans)}.qcbox .rs{font:11px var(--mono);color:var(--faint);overflow-wrap:anywhere}
.counts{display:grid;grid-template-columns:repeat(4,minmax(0,1fr));gap:8px;margin:12px 0 10px}
.count{border:1px solid var(--line);border-top:3px solid;padding:6px 8px;border-radius:2px;min-width:0}
.count small{display:block;font:500 10px var(--mono);color:var(--muted);white-space:nowrap;overflow:hidden;text-overflow:ellipsis}.count strong{font:500 22px/1.2 var(--mono)}
.count.pass{border-top-color:var(--pass)}.count.warning{border-top-color:var(--warning)}.count.fail{border-top-color:var(--fail)}.count.not_evaluated{border-top-color:var(--ne)}
.sbar{display:flex;height:8px;border-radius:2px;overflow:hidden;background:var(--line2)}.sbar i{display:block}
.sbar .pass{background:var(--pass)}.sbar .warning{background:var(--warning)}.sbar .fail{background:var(--fail)}.sbar .not_evaluated{background:var(--ne)}
.table-wrap{overflow-x:auto}
table{border-collapse:collapse;width:100%}
.grid th,.grid td{border-bottom:1px solid var(--line2);padding:5px 12px;text-align:left;font-weight:400;white-space:nowrap}
.grid thead th{font:600 10.5px var(--mono);color:var(--muted);text-transform:uppercase;letter-spacing:.04em;background:var(--sub);border-bottom:1px solid var(--line)}
.grid td.n{text-align:right;font:12.5px var(--mono)}.grid thead th.n{text-align:right}.grid th.c{text-align:center}
.grid tbody tr:hover{background:var(--row-hover)}
.grid .row-h{color:var(--muted);white-space:normal}
.hm td{text-align:center}.hm td .cell{margin:0 auto;width:44px;height:18px;font-size:11px}
.hm td .cell.warning::after{font-size:9px}
.tracks{display:grid;grid-template-columns:repeat(2,minmax(0,1fr))}
.lane{min-width:0}.lane+.lane{border-left:1px solid var(--line)}
.lane-h{display:flex;align-items:center;justify-content:space-between;padding:6px 12px;border-bottom:1px solid var(--line2);font:600 12px var(--mono)}
.chart-button{display:block;width:100%;padding:6px 6px 0;border:0;background:var(--panel);cursor:zoom-in;color:inherit;font:inherit}
.chart{display:block}
.chart-svg,.chart-img{display:block;width:100%;height:auto}
.chart-svg path,.chart-svg use{stroke-linejoin:round;stroke-linecap:butt}
.empty{margin:12px;padding:30px 12px;border:1px dashed var(--line);text-align:center;color:var(--muted);font:12px var(--mono);overflow-wrap:anywhere}
.error-box{border-color:var(--fail);color:var(--fail-text)}
.metrics{border-top:1px solid var(--line2)}
.metrics table td{padding:4px 12px;border-bottom:1px solid var(--line2)}
.metrics table td:last-child{text-align:right;font:500 12.5px var(--mono);white-space:nowrap}
.why{padding:6px 12px 10px;font:12px/1.45 var(--mono);color:var(--muted)}
.why::before{content:"› ";color:var(--brand)}
.orep{border-top:1px solid var(--line)}
.orep-h{display:flex;justify-content:space-between;align-items:center;padding:6px 12px;background:var(--sub);border-bottom:1px solid var(--line2);font:600 11px var(--mono);text-transform:uppercase;color:var(--muted)}
.table-count{display:inline-grid;place-items:center;min-width:22px;height:18px;padding:0 6px;border:1px solid var(--line);border-radius:2px;background:var(--panel);color:var(--ink);font:600 11px var(--mono)}
.orep code{font:11.5px/1.4 var(--mono);word-break:break-all}
.orep td:first-child{white-space:normal}
.orep .none{padding:8px 12px;margin:0;color:var(--muted);font:12px var(--mono)}
.divider{display:none;align-items:center;gap:10px;margin:6px 2px -2px;font:600 11px var(--mono);text-transform:uppercase;letter-spacing:.06em;color:var(--muted)}
.divider::after{content:"";flex:1;height:1px;background:var(--line)}
.divider b{display:inline-grid;place-items:center;min-width:20px;height:18px;padding:0 5px;border-radius:2px;background:var(--panel);border:1px solid var(--line);color:var(--ink);font-weight:600}
.divider.attention{color:var(--fail-text)}
.tracklist .divider{margin:10px 16px 4px;font-size:10px}
body[data-order="triage"] .divider{display:flex}
body[data-order="triage"] .panel[data-sev="fail"]{border-left:3px solid var(--fail)}
body[data-order="triage"] .panel[data-sev="warning"]{border-left:3px solid var(--warning)}
body[data-order="triage"] .panel[data-sev="not_evaluated"]{border-left:3px solid var(--ne)}
.phead[data-fold]{cursor:pointer;user-select:none}
.phead[data-fold]:hover{background:var(--head-hover)}
.fold{flex:none;width:22px;height:22px;margin-left:-4px;border:0;border-radius:2px;background:none;display:grid;place-items:center;color:var(--muted);cursor:pointer;transition:transform .15s}
.fold:hover{background:var(--fold-hover)}
.fold svg{display:block}
.panel.collapsed .fold{transform:rotate(-90deg)}
.panel.collapsed>:not(.phead){display:none}
.panel.collapsed>.phead{border-bottom:0}
.peek{display:none;gap:14px;font:12px var(--mono);color:var(--muted);white-space:nowrap;overflow:hidden;min-width:0}
.peek b{font-weight:500;color:var(--ink)}
.panel.collapsed .peek{display:flex}
#fold-all[data-state="collapse"] .fa-e,#fold-all[data-state="expand"] .fa-c{display:none}
footer{padding:8px 4px 24px;font:11px var(--mono);color:var(--faint)}
dialog{border:1px solid var(--line);border-radius:3px;padding:0;width:min(1200px,94vw);background:var(--panel);color:var(--ink)}
dialog::backdrop{background:var(--backdrop)}
.dialog-head{display:flex;justify-content:flex-end;padding:6px;border-bottom:1px solid var(--line);background:var(--sub)}
.dialog-close{border:1px solid var(--line);background:var(--panel);color:var(--ink);padding:3px 10px;font:12px var(--sans);cursor:pointer}
#dialog-body{padding:8px}
@media(prefers-reduced-motion:reduce){.fold{transition:none}}
@media(max-width:1000px){.shell{grid-template-columns:minmax(0,1fr)}.rail{position:relative;top:0;height:auto}.logo{width:auto}
.crumbs .cr:nth-child(n+3){display:none}.hero{grid-template-columns:minmax(0,1fr)}.hero .id{border-right:0;border-bottom:1px solid var(--line)}
.appbar .seg small{display:none}.appbar .btn.print{display:none}
.rail-h{cursor:pointer;padding:8px 0;margin-top:0}.rail-h>span:first-child::after{content:" ▾"}.rail.open .rail-h>span:first-child::after{content:" ▴"}
.rail:not(.open) .tracklist,.rail:not(.open) .legend{display:none}
#dialog-body{overflow-x:auto}#dialog-body .chart-svg,#dialog-body .chart-img{min-width:720px}}
@media(max-width:700px){.hm td .cell{width:28px}.grid th,.grid td{padding:5px 8px}.tracks{grid-template-columns:minmax(0,1fr)}
.lane+.lane{border-left:0;border-top:1px solid var(--line)}.counts{grid-template-columns:repeat(2,minmax(0,1fr))}.main{padding:8px}
.peek{display:none!important}.phead h2{min-width:0;flex:1 1 auto}.phead{flex-wrap:wrap}.appbar{height:auto;flex-wrap:wrap}
.logo{border-right:0;height:44px}.crumbs{display:none}.appbar .seg{margin:0 8px 8px 12px}.kv{grid-template-columns:repeat(2,minmax(0,1fr))}}
@media print{.appbar,.rail,.fold,dialog{display:none!important}.print-meta{display:block}.shell{display:block}body{background:#fff}.main{padding:0}
.panel{break-inside:avoid}.panel.collapsed>:not(.phead){display:revert!important}.panel.collapsed .tracks{display:grid!important}
.panel.collapsed .peek{display:none!important}.chart-button{cursor:default}}
"""
)

# Runs in <head> so the saved theme and language apply before the first paint.
_HEAD_JS = """
(function(){var d=document.documentElement,t=null,l=null;
try{t=localStorage.getItem('neoqc.report.theme');l=localStorage.getItem('neoqc.report.lang');}catch(e){}
if(t!=='light'&&t!=='dark')t=(window.matchMedia&&matchMedia('(prefers-color-scheme: dark)').matches)?'dark':'light';
if(l!=='en'&&l!=='ru')l=((navigator.language||'').toLowerCase().indexOf('ru')===0)?'ru':'en';
d.setAttribute('data-theme',t);d.setAttribute('lang',l);})();
""".strip()

_JS = """
const root=document.documentElement;
const store=(key,value)=>{try{localStorage.setItem(key,value);}catch(e){}};
const press=(attr,value)=>document.querySelectorAll('['+attr+']').forEach((b)=>b.setAttribute('aria-pressed',String(b.getAttribute(attr)===value)));

// Scroll-spy: highlight the rail entry of the section at the top of the viewport.
const navLinks=[...document.querySelectorAll('[data-nav]')];
let clickedHash=null;
const spy=()=>{
  const sections=[...document.querySelectorAll('[data-section]')].filter((s)=>s.offsetParent!==null);
  if(!sections.length)return;
  let active=sections[0];
  for(const section of sections){if(section.getBoundingClientRect().top<=120)active=section;else break;}
  if(window.scrollY+window.innerHeight>=document.documentElement.scrollHeight-2){
    // At the bottom several short sections are visible at once: keep the one the user picked.
    const picked=clickedHash&&sections.find((s)=>'#'+s.id===clickedHash);
    active=picked&&picked.getBoundingClientRect().top<window.innerHeight?picked:sections[sections.length-1];
  }
  navLinks.forEach((link)=>{const on=link.getAttribute('href')==='#'+active.id;link.classList.toggle('active',on);
    if(on)link.setAttribute('aria-current','true');else link.removeAttribute('aria-current');});
};
let spyFrame=0;
window.addEventListener('scroll',()=>{if(spyFrame)return;spyFrame=requestAnimationFrame(()=>{spyFrame=0;spy();});},{passive:true});

// Theme and language.
const setTheme=(theme,save)=>{root.setAttribute('data-theme',theme);press('data-theme-btn',theme);if(save)store('neoqc.report.theme',theme);};
const setLang=(lang,save)=>{root.setAttribute('lang',lang);press('data-lang-btn',lang);if(save)store('neoqc.report.lang',lang);spy();};
document.querySelectorAll('[data-theme-btn]').forEach((b)=>b.addEventListener('click',()=>setTheme(b.dataset.themeBtn,true)));
document.querySelectorAll('[data-lang-btn]').forEach((b)=>b.addEventListener('click',()=>setLang(b.dataset.langBtn,true)));
setTheme(root.getAttribute('data-theme')==='dark'?'dark':'light',false);
setLang(root.getAttribute('lang')==='ru'?'ru':'en',false);

// Section order (triage / sequential) and folding.
const ORDER_KEY='neoqc.report.order';
const main=document.querySelector('.main'),tracklist=document.querySelector('.tracklist'),heatmap=document.querySelector('.hm tbody');
const foldables=[...document.querySelectorAll('.panel.foldable')];
const foldAll=document.getElementById('fold-all');
const sortBy=(box,key)=>{if(!box)return;[...box.children].sort((a,b)=>(+a.dataset[key])-(+b.dataset[key])).forEach((e)=>box.appendChild(e));};
const setFold=(panel,folded)=>{panel.classList.toggle('collapsed',folded);
  const button=panel.querySelector('.fold');if(button)button.setAttribute('aria-expanded',String(!folded));};
const syncFoldAll=()=>{if(foldAll)foldAll.dataset.state=foldables.some((p)=>!p.classList.contains('collapsed'))?'collapse':'expand';};
const applyOrder=(mode,resetFold)=>{
  document.body.dataset.order=mode;
  [main,tracklist,heatmap].forEach((box)=>sortBy(box,mode==='triage'?'tri':'seq'));
  press('data-order-btn',mode);
  if(resetFold)foldables.forEach((p)=>setFold(p,mode==='triage'&&p.dataset.sev==='ok'));
  syncFoldAll();spy();store(ORDER_KEY,mode);
};
document.querySelectorAll('[data-order-btn]').forEach((b)=>b.addEventListener('click',()=>applyOrder(b.dataset.orderBtn,true)));
document.querySelectorAll('.phead[data-fold]').forEach((head)=>head.addEventListener('click',()=>{
  const panel=head.parentElement;setFold(panel,!panel.classList.contains('collapsed'));syncFoldAll();spy();}));
if(foldAll)foldAll.addEventListener('click',()=>{const fold=foldables.some((p)=>!p.classList.contains('collapsed'));
  foldables.forEach((p)=>setFold(p,fold));syncFoldAll();spy();});
const reveal=(hash)=>{let target=null;try{target=hash&&document.querySelector(hash);}catch(e){}
  if(target&&target.classList.contains('collapsed')){setFold(target,false);syncFoldAll();}};
navLinks.forEach((link)=>link.addEventListener('click',()=>{clickedHash=link.getAttribute('href');reveal(clickedHash);
  requestAnimationFrame(spy);}));

// Narrow screens: the section list is folded behind the "Report sections" header.
const rail=document.querySelector('.rail'),railToggle=document.querySelector('.rail-h');
if(railToggle)railToggle.addEventListener('click',()=>{const open=!rail.classList.contains('open');
  rail.classList.toggle('open',open);railToggle.setAttribute('aria-expanded',String(open));});

// Full-size chart dialog.
const dialog=document.getElementById('chart-dialog'),dialogBody=document.getElementById('dialog-body');
document.querySelectorAll('.chart-button').forEach((b)=>b.addEventListener('click',()=>{dialogBody.innerHTML=b.innerHTML;dialog.showModal();}));
document.getElementById('dialog-close').addEventListener('click',()=>dialog.close());
dialog.addEventListener('click',(event)=>{if(event.target===dialog)dialog.close();});
dialog.addEventListener('close',()=>{dialogBody.innerHTML='';});

// Printing always uses the light theme and expands every section.
let printFolds=null,printTheme=null;
window.addEventListener('beforeprint',()=>{printFolds=foldables.map((p)=>p.classList.contains('collapsed'));
  foldables.forEach((p)=>setFold(p,false));printTheme=root.getAttribute('data-theme');root.setAttribute('data-theme','light');});
window.addEventListener('afterprint',()=>{if(printFolds)foldables.forEach((p,i)=>setFold(p,printFolds[i]));printFolds=null;syncFoldAll();
  if(printTheme)root.setAttribute('data-theme',printTheme);printTheme=null;});

let savedOrder=null;try{savedOrder=localStorage.getItem(ORDER_KEY);}catch(e){}
applyOrder(savedOrder==='sequential'?'sequential':'triage',true);
reveal(location.hash);
""".strip()


def render_qc_report(model: QcReportModel, plot_dir: Path) -> str:
    attention = [module for module in model.modules if _needs_attention(module.qc_status)]
    passed = [module for module in model.modules if not _needs_attention(module.qc_status)]
    k_overview, k_basic, k_footer = _keys(0, 0), _keys(1, 1), _keys(99999, 99999)
    k_attention, k_passed = _keys(99990, 50), _keys(99991, 5000)

    def divider(kind: str, label: str, count: int, keys: str) -> str:
        if not count:
            return ""
        return f'<div class="divider {kind}" {keys} aria-hidden="true">{_t(label)} <b>{count}</b></div>'

    def cells(module: PlotModule) -> str:
        by_read = {card.read: card for card in module.cards}
        return "".join(
            _cell(by_read[read].qc_status, f"{read}: {STATUS_LABELS[by_read[read].qc_status]}")
            if read in by_read else _cell("not_evaluated", f"{read}: —")
            for read in model.reads
        )

    tracks = [
        f'<a class="track" data-nav href="#qc-summary" {k_overview}><span class="i">00</span>'
        f'<span class="l">{_t("QC overview")}</span><span class="cells">'
        f'{_cell(model.overall_qc_status, STATUS_LABELS[model.overall_qc_status])}</span></a>',
        f'<a class="track" data-nav href="#basic-statistics" {k_basic}><span class="i">01</span>'
        f'<span class="l">{_t("Basic Statistics")}</span><span class="cells">{_cell("info", "INFO")}</span></a>',
        divider("attention", "Needs attention", len(attention), k_attention),
        divider("", "Passed", len(passed), k_passed),
    ]
    tracks += [
        f'<a class="track" data-nav href="#module-{escape(module.metric_id)}" '
        f'{_keys(index, _triage_rank(module.qc_status, index))}><span class="i">{index:02d}</span>'
        f'<span class="l">{_t(module.title)}</span><span class="cells">{cells(module)}</span></a>'
        for index, module in enumerate(model.modules, start=2)
    ]
    reads_head = "".join(f"<b>{escape(read)}</b>" for read in model.reads)
    legend = "".join(f"<span>{_cell(status, STATUS_LABELS[status])}{STATUS_LABELS[status]}</span>" for status in STATUS_ORDER)
    prefs = (
        f'<div class="prefs"><div class="pref"><small>{_t("Theme")}</small><div class="seg" role="group" aria-label="Theme">'
        f'<button type="button" data-theme-btn="light" aria-pressed="true" aria-label="Light theme" title="Light">{SUN}</button>'
        f'<button type="button" data-theme-btn="dark" aria-pressed="false" aria-label="Dark theme" title="Dark">{MOON}</button>'
        f'</div></div><div class="pref"><small>{_t("Language")}</small><div class="seg" role="group" aria-label="Language">'
        '<button type="button" data-lang-btn="en" aria-pressed="true" lang="en">EN</button>'
        '<button type="button" data-lang-btn="ru" aria-pressed="false" lang="ru">RU</button></div></div></div>'
    )
    rail = (
        f'<aside class="rail">{prefs}<button class="rail-h" type="button" aria-expanded="false" aria-controls="tracklist">'
        f'<span>{_t("Report sections")}</span><em>{reads_head}</em></button>'
        f'<nav class="tracklist" id="tracklist" aria-label="Report sections">{"".join(tracks)}</nav><div class="legend">{legend}</div></aside>'
    )
    appbar = (
        '<header class="appbar"><div class="logo"><i aria-hidden="true"><b></b><b></b><b></b><b></b></i>NeoQC Report</div>'
        '<div class="crumbs">'
        f'<span class="cr"><small>{_t("Sample")}</small><strong>{escape(model.sample_id)}</strong></span>'
        f'<span class="cr"><small>{_t("Reads")}</small><strong>{escape(" / ".join(model.reads))}</strong></span>'
        f'<span class="cr"><small>{_t("Generated")}</small><strong>{escape(model.generated_at)}</strong></span>'
        f'<span class="cr"><small>{_t("Ruleset")}</small><strong>{escape(model.ruleset_label)}</strong></span></div>'
        f'<div class="seg" role="group" aria-label="Section order"><small>{_t("Order")}</small>'
        f'<button type="button" data-order-btn="triage" aria-pressed="true">{_t("Triage")}</button>'
        f'<button type="button" data-order-btn="sequential" aria-pressed="false">{_t("Sequential")}</button></div>'
        f'<button class="btn ghost" id="fold-all" type="button" data-state="collapse"><span class="fa-c">{_t("Collapse all")}</span>'
        f'<span class="fa-e">{_t("Expand all")}</span></button>'
        f'<button class="btn print" type="button" onclick="window.print()">{_t("Print / Save PDF")}</button></header>'
    )
    main = (
        _render_overview(model, k_overview)
        + _render_basic_statistics(model, k_basic)
        + divider("attention", "Needs attention", len(attention), k_attention)
        + divider("", "Passed · folded", len(passed), k_passed)
        + "".join(_render_module(module, index, plot_dir) for index, module in enumerate(model.modules, start=2))
        + f'<footer {k_footer}>{_t("Generated by NeoQC")}</footer>'
    )
    dialog = (
        '<dialog id="chart-dialog" aria-label="Chart"><div class="dialog-head">'
        f'<button class="dialog-close" id="dialog-close" type="button">{_t("Close")}</button></div>'
        '<div id="dialog-body"></div></dialog>'
    )
    return (
        '<!doctype html>\n<html lang="en"><head><meta charset="utf-8">'
        '<meta name="viewport" content="width=device-width,initial-scale=1">'
        f"<title>NeoQC report — {escape(model.sample_id)}</title>"
        f"<style>{_font_face_css()}\n{_CSS}</style><script>{_HEAD_JS}</script></head>"
        f'<body>{appbar}<div class="shell">{rail}<main class="main">{main}</main></div>{dialog}'
        f"<script>{_JS}</script></body></html>"
    )


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
