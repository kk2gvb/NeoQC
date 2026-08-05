"""Validated TSV-to-chart rendering for NeoQC."""

from __future__ import annotations

import csv
import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Mapping, Sequence

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter

from plot_style import (
    ACCENT,
    BASE_COLORS,
    BRAND,
    BRAND_DARK,
    DANGER,
    FIGURE_SIZE,
    GRID,
    INK,
    LINE_STYLES,
    MARKERS,
    MUTED,
    PANEL,
    PNG_DPI,
    SERIES_COLORS,
    WARNING,
    WHITE,
    apply_theme,
    compact_number,
    finish_figure,
    setup_axes,
    use_compact_y_axis,
)


class PlotDataError(ValueError):
    """Raised when a NeoQC TSV cannot be rendered safely."""


Rows = list[dict[str, float]]
PlotFunction = Callable[[Rows, str], tuple[plt.Figure, str]]


@dataclass(frozen=True)
class MetricSpec:
    metric_id: str
    source_prefix: str
    output_prefix: str
    title: str
    required_columns: tuple[str, ...]
    plot: PlotFunction
    adapters_only: bool = False
    variable_series: bool = False

    def source_name(self, read: str) -> str:
        return f"{self.source_prefix}_{read}.tsv"

    def output_name(self, read: str, suffix: str) -> str:
        return f"{self.output_prefix}_{read}.{suffix}"


def _read_numeric_tsv(path: Path, spec: MetricSpec) -> Rows:
    try:
        handle = path.open("r", encoding="utf-8", newline="")
    except OSError as error:
        raise PlotDataError(f"cannot read {path.name}: {error}") from error

    with handle:
        reader = csv.DictReader(handle, delimiter="\t")
        columns = tuple(reader.fieldnames or ())
        if not columns or any(not column.strip() for column in columns):
            raise PlotDataError(f"{path.name}: contains an empty column name")
        if len(set(columns)) != len(columns):
            raise PlotDataError(f"{path.name}: contains duplicate column names")
        missing = [column for column in spec.required_columns if column not in columns]
        if missing:
            raise PlotDataError(
                f"{path.name}: missing required column(s): {', '.join(missing)}"
            )
        if spec.variable_series and len(columns) <= len(spec.required_columns):
            raise PlotDataError(f"{path.name}: contains no data series")
        if not spec.variable_series:
            unexpected = [column for column in columns if column not in spec.required_columns]
            if unexpected:
                raise PlotDataError(
                    f"{path.name}: unexpected column(s): {', '.join(unexpected)}"
                )

        rows: Rows = []
        for line_number, raw in enumerate(reader, start=2):
            if None in raw:
                raise PlotDataError(f"{path.name}:{line_number}: too many fields")
            converted: dict[str, float] = {}
            for column in columns:
                value = (raw.get(column) or "").strip()
                try:
                    number = float(value)
                except ValueError as error:
                    raise PlotDataError(
                        f"{path.name}:{line_number}: {column} is not numeric"
                    ) from error
                if not math.isfinite(number):
                    raise PlotDataError(
                        f"{path.name}:{line_number}: {column} is not finite"
                    )
                if number < 0:
                    raise PlotDataError(
                        f"{path.name}:{line_number}: {column} must not be negative"
                    )
                converted[column] = number
            rows.append(converted)

    if not rows:
        raise PlotDataError(f"{path.name}: contains no data rows")
    return rows


def _values(rows: Rows, column: str) -> list[float]:
    return [row[column] for row in rows]


def _ensure_range(values: Sequence[float], label: str, lower: float, upper: float) -> None:
    if any(value < lower or value > upper for value in values):
        raise PlotDataError(f"{label} must be between {lower:g} and {upper:g}")


def _weighted_summary(x: Sequence[float], weights: Sequence[float]) -> tuple[float, float, float]:
    total = sum(weights)
    if total <= 0:
        return 0.0, 0.0, x[0] if x else 0.0
    mean = sum(value * weight for value, weight in zip(x, weights)) / total
    half = total / 2.0
    cumulative = 0.0
    median = x[-1]
    for value, weight in sorted(zip(x, weights)):
        cumulative += weight
        if cumulative >= half:
            median = value
            break
    mode = x[max(range(len(weights)), key=weights.__getitem__)]
    return mean, median, mode


def _mark_statistic(ax, value: float, label: str, color: str, linestyle: str) -> None:
    ax.axvline(value, color=color, linestyle=linestyle, linewidth=1.3, label=f"{label}: {value:.1f}")


def _line_marker_stride(values: Sequence[float]) -> int:
    return max(1, len(values) // 16)


def _fastqc_theoretical_gc(
    x: Sequence[float], counts: Sequence[float]
) -> tuple[list[int], list[float], float, float]:
    """Reproduce FastQC's mode-centred normal model on the 0..100 GC scale."""

    observed = [0.0] * 101
    for percentage, count in zip(x, counts):
        index = int(round(percentage))
        if not math.isclose(percentage, index, abs_tol=1e-6):
            raise PlotDataError("GC percentage values must be whole numbers")
        observed[index] += count

    total = sum(observed)
    if total <= 1:
        raise PlotDataError("GC distribution needs at least two reads")

    first_mode = max(range(101), key=observed.__getitem__)
    threshold = observed[first_mode] * 0.9
    modal_bins = [first_mode]
    for index in range(first_mode + 1, 101):
        if observed[index] <= threshold:
            break
        modal_bins.append(index)
    for index in range(first_mode - 1, -1, -1):
        if observed[index] <= threshold:
            break
        modal_bins.append(index)

    # FastQC keeps the first mode when the high plateau touches either edge.
    touches_edge = min(modal_bins) == 0 or max(modal_bins) == 100
    centre = float(first_mode) if touches_edge else sum(modal_bins) / len(modal_bins)
    variance = sum(((index - centre) ** 2) * count for index, count in enumerate(observed))
    stdev = math.sqrt(variance / (total - 1))
    if stdev <= 0:
        theoretical = [0.0] * 101
        theoretical[int(round(centre))] = total
    else:
        scale = total / (math.sqrt(2.0 * math.pi) * stdev)
        theoretical = [
            scale * math.exp(-((index - centre) ** 2) / (2.0 * stdev * stdev))
            for index in range(101)
        ]
    deviation = sum(abs(model - actual) for model, actual in zip(theoretical, observed)) / total * 100.0
    return list(range(101)), theoretical, centre, deviation


def plot_per_base_quality(rows: Rows, read: str) -> tuple[plt.Figure, str]:
    x = _values(rows, "cycle")
    y = _values(rows, "mean_quality")
    fig, ax = plt.subplots(figsize=FIGURE_SIZE)
    setup_axes(ax, "Per base sequence quality", "Position in read (bp)", "Mean Phred quality", read)
    upper = max(42.0, max(y) + 3.0)
    ax.axhspan(0, 20, color=DANGER, alpha=0.09, linewidth=0)
    ax.axhspan(20, 30, color=WARNING, alpha=0.10, linewidth=0)
    ax.axhspan(30, upper, color=ACCENT, alpha=0.09, linewidth=0)
    ax.axhline(20, color=DANGER, linewidth=0.9, linestyle="--", alpha=0.8)
    ax.axhline(30, color=ACCENT, linewidth=0.9, linestyle="--", alpha=0.8)
    ax.plot(x, y, color=BRAND, marker="o", markevery=_line_marker_stride(x), markersize=3.2, zorder=3)
    ax.set_xlim(min(x), max(x) if len(x) > 1 else min(x) + 1)
    ax.set_ylim(0, upper)
    label_x = min(x) + (max(x) - min(x)) * 0.01
    ax.text(label_x, 30.4, "Q30", color=ACCENT, fontsize=7, va="bottom")
    ax.text(label_x, 20.4, "Q20", color=WARNING, fontsize=7, va="bottom")
    finish_figure(fig)
    return fig, f"Mean Phred quality across {len(x)} {read} read positions."


def plot_per_sequence_quality(rows: Rows, read: str) -> tuple[plt.Figure, str]:
    x = _values(rows, "mean_quality")
    counts = _values(rows, "read_count")
    mean, median, mode = _weighted_summary(x, counts)
    total = sum(counts)
    if total <= 0:
        raise PlotDataError("Per sequence quality scores contains no observations")
    shares = [count / total * 100.0 for count in counts]
    fig, ax = plt.subplots(figsize=FIGURE_SIZE)
    setup_axes(
        ax,
        "Per sequence quality scores",
        "Mean Phred quality per read",
        "Reads (%)",
        read,
    )
    ax.axvspan(min(x), 20, color=DANGER, alpha=0.055, linewidth=0)
    ax.axvspan(20, 30, color=WARNING, alpha=0.065, linewidth=0)
    ax.axvspan(30, max(x), color=ACCENT, alpha=0.055, linewidth=0)
    ax.fill_between(x, shares, color=ACCENT, alpha=0.18, linewidth=0)
    ax.plot(x, shares, color=BRAND_DARK, linewidth=2.2, marker="o", markersize=3.2, zorder=3)
    ax.set_ylim(0, min(100.0, max(shares) * 1.18 + 1.0))
    ax.yaxis.set_major_formatter(FuncFormatter(lambda value, _: f"{value:g}%"))

    x_padding = max(0.8, (max(x) - min(x)) * 0.018)
    ax.set_xlim(min(x) - x_padding, max(x) + x_padding)
    ax.text(
        0.02,
        0.96,
        f"Mean Q{mean:.1f}  •  Median Q{median:.0f}  •  Mode Q{mode:.0f}",
        transform=ax.transAxes,
        ha="left",
        va="top",
        color=INK,
        fontsize=8,
        fontweight=600,
        bbox={"boxstyle": "round,pad=0.45", "facecolor": WHITE, "edgecolor": GRID, "alpha": 0.94},
    )
    finish_figure(fig)
    return fig, (
        f"Per sequence quality scores for {read}, shown as normalized percentages; "
        f"mean Q{mean:.1f}, median Q{median:.1f}, mode Q{mode:.0f}."
    )


def plot_adapter_content(rows: Rows, read: str) -> tuple[plt.Figure, str]:
    x = _values(rows, "pos")
    series = [column for column in rows[0] if column != "pos"]
    fig, ax = plt.subplots(figsize=FIGURE_SIZE)
    setup_axes(ax, "Adapter content", "Position in read (bp)", "Reads with adapter (%)", read)
    plotted = 0
    peak = 0.0
    for index, name in enumerate(series):
        y = _values(rows, name)
        _ensure_range(y, f"adapter series {name}", 0, 100)
        peak = max(peak, max(y))
        if max(y) <= 0:
            continue
        ax.plot(
            x,
            y,
            color=SERIES_COLORS[index % len(SERIES_COLORS)],
            linestyle=LINE_STYLES[index % len(LINE_STYLES)],
            marker=MARKERS[index % len(MARKERS)],
            markevery=_line_marker_stride(x),
            markersize=3,
            label=name,
        )
        plotted += 1
    if plotted:
        ax.legend(loc="upper left", ncols=min(2, plotted))
    else:
        ax.text(
            0.5,
            0.52,
            "No adapters detected",
            transform=ax.transAxes,
            ha="center",
            va="center",
            color=MUTED,
            fontsize=12,
            fontweight=600,
            bbox={"boxstyle": "round,pad=0.7", "facecolor": PANEL, "edgecolor": GRID},
        )
    ax.set_xlim(min(x), max(x) if len(x) > 1 else min(x) + 1)
    ax.set_ylim(0, max(5.0, min(100.0, peak * 1.2 + 1.0)))
    finish_figure(fig)
    return fig, f"Adapter content by {read} read position; peak {peak:.2f}%."


def plot_base_content(rows: Rows, read: str) -> tuple[plt.Figure, str]:
    x = _values(rows, "position")
    fig, ax = plt.subplots(figsize=FIGURE_SIZE)
    setup_axes(ax, "Per base sequence content", "Position in read (bp)", "Base content (%)", read)
    for index, base in enumerate(("A", "C", "G", "T", "N")):
        y = _values(rows, base)
        _ensure_range(y, f"base content {base}", 0, 100)
        ax.plot(
            x,
            y,
            color=BASE_COLORS[base],
            linestyle=LINE_STYLES[index % len(LINE_STYLES)],
            marker=MARKERS[index % len(MARKERS)],
            markevery=_line_marker_stride(x),
            markersize=2.8,
            label=base,
        )
    ax.axhline(25, color=GRID, linewidth=1.0, linestyle="--", zorder=0)
    ax.set_xlim(min(x), max(x) if len(x) > 1 else min(x) + 1)
    ax.set_ylim(0, 100)
    ax.legend(loc="upper left", ncols=5)
    finish_figure(fig)
    return fig, f"Relative A, C, G, T and N content across {len(x)} {read} positions."


def plot_gc_content(rows: Rows, read: str) -> tuple[plt.Figure, str]:
    x = _values(rows, "gc_percent")
    _ensure_range(x, "GC percentage", 0, 100)
    counts = _values(rows, "reads")
    mean, median, mode = _weighted_summary(x, counts)
    theoretical_x, theoretical, theoretical_centre, deviation = _fastqc_theoretical_gc(x, counts)
    fig, ax = plt.subplots(figsize=FIGURE_SIZE)
    setup_axes(ax, "Per sequence GC content", "GC content (%)", "Reads", read)
    ax.fill_between(x, counts, color=ACCENT, alpha=0.20, linewidth=0)
    ax.plot(x, counts, color=BRAND_DARK, linewidth=2.2, label="Observed")
    ax.plot(
        theoretical_x,
        theoretical,
        color=WARNING,
        linewidth=2.0,
        linestyle="--",
        label="Theoretical distribution",
    )
    _mark_statistic(ax, mean, "Mean GC", ACCENT, "--")
    ax.set_xlim(0, 100)
    ax.set_title("Per sequence GC content", loc="left", color=INK, pad=46)
    if ax.texts:
        ax.texts[0].set_y(1.13)
    ax.legend(
        loc="lower left",
        bbox_to_anchor=(0.0, 1.015),
        borderaxespad=0,
        ncols=3,
    )
    use_compact_y_axis(ax)
    finish_figure(fig)
    return fig, (
        f"Observed and theoretical GC distributions for {read}; mean {mean:.1f}%, "
        f"median {median:.1f}%, mode {mode:.0f}%, theoretical centre "
        f"{theoretical_centre:.1f}%, deviation {deviation:.1f}%."
    )


def plot_n_content(rows: Rows, read: str) -> tuple[plt.Figure, str]:
    x = _values(rows, "position")
    y = _values(rows, "N_percent")
    _ensure_range(y, "N percentage", 0, 100)
    peak = max(y)
    fig, ax = plt.subplots(figsize=FIGURE_SIZE)
    setup_axes(ax, "Per base N content", "Position in read (bp)", "N content (%)", read)
    ax.fill_between(x, y, color=ACCENT, alpha=0.18, linewidth=0)
    ax.plot(x, y, color=BRAND, marker="o", markevery=_line_marker_stride(x), markersize=3)
    peak_index = y.index(peak)
    ax.scatter([x[peak_index]], [peak], color=DANGER if peak > 0 else ACCENT, s=24, zorder=4)
    ax.annotate(
        f"Peak {peak:.2f}%",
        (x[peak_index], peak),
        xytext=(8, 10),
        textcoords="offset points",
        color=INK,
        fontsize=8,
    )
    ax.set_xlim(min(x), max(x) if len(x) > 1 else min(x) + 1)
    ax.set_ylim(0, max(1.0, min(100.0, peak * 1.25 + 0.2)))
    finish_figure(fig)
    return fig, f"N content by {read} position; peak {peak:.2f}%."


def plot_length_distribution(rows: Rows, read: str) -> tuple[plt.Figure, str]:
    x = _values(rows, "length")
    counts = _values(rows, "reads")
    mean, median, mode = _weighted_summary(x, counts)
    fig, ax = plt.subplots(figsize=FIGURE_SIZE)
    setup_axes(ax, "Sequence length distribution", "Read length (bp)", "Reads", read)
    by_length: dict[float, float] = {}
    for length, count in zip(x, counts):
        by_length[length] = by_length.get(length, 0.0) + count
    lengths = sorted(by_length)
    values = [by_length[length] for length in lengths]

    if len(lengths) == 1:
        length = lengths[0]
        total = values[0]
        plot_x = [length - 1.0, length, length + 1.0]
        plot_y = [0.0, total, 0.0]
        ax.axvspan(length - 0.5, length + 0.5, color=GRID, alpha=0.52, linewidth=0, zorder=0)
        ax.fill_between(plot_x, plot_y, color=ACCENT, alpha=0.10, linewidth=0, zorder=1)
        ax.plot(plot_x, plot_y, color=BRAND_DARK, linewidth=2.2, marker="o", markersize=4, zorder=3)
        ax.set_xlim(length - 1.5, length + 1.5)
        ax.set_ylim(0, total * 1.18 if total > 0 else 1)
        ax.text(
            0.02,
            0.94,
            f"Fixed length  •  {length:g} bp",
            transform=ax.transAxes,
            ha="left",
            va="top",
            color=INK,
            fontsize=8.5,
            fontweight=600,
            bbox={"boxstyle": "round,pad=0.45", "facecolor": PANEL, "edgecolor": GRID},
        )
        ax.annotate(
            f"{compact_number(total)} reads",
            (length, total),
            xytext=(0, 8),
            textcoords="offset points",
            ha="center",
            va="bottom",
            color=INK,
            fontsize=8,
            fontweight=600,
        )
        use_compact_y_axis(ax)
        finish_figure(fig)
        return fig, (
            f"FastQC-style sequence length line for {read}; all {total:g} reads "
            f"have a fixed length of {length:g} bp."
        )

    gaps = [right - left for left, right in zip(lengths, lengths[1:]) if right > left]
    endpoint_gap = min(gaps) if gaps else 1.0
    plot_x = [lengths[0] - endpoint_gap, *lengths, lengths[-1] + endpoint_gap]
    plot_y = [0.0, *values, 0.0]
    if len(lengths) <= 50:
        band_width = endpoint_gap / 2.0
        for index, length in enumerate(lengths):
            if index % 2 == 0:
                ax.axvspan(length - band_width, length + band_width, color=GRID, alpha=0.34, linewidth=0)
    ax.fill_between(plot_x, plot_y, color=ACCENT, alpha=0.10, linewidth=0, zorder=1)
    ax.plot(plot_x, plot_y, color=BRAND_DARK, linewidth=2.1, marker="o", markersize=3.5, zorder=3)
    ax.set_xlim(plot_x[0], plot_x[-1])
    ax.set_ylim(0, max(values) * 1.18 if max(values) > 0 else 1)
    ax.text(
        0.02,
        0.96,
        f"Mean {mean:.1f} bp  •  Median {median:.0f} bp  •  Mode {mode:.0f} bp",
        transform=ax.transAxes,
        ha="left",
        va="top",
        color=INK,
        fontsize=8,
        fontweight=600,
        bbox={"boxstyle": "round,pad=0.45", "facecolor": WHITE, "edgecolor": GRID, "alpha": 0.94},
    )
    use_compact_y_axis(ax)
    finish_figure(fig)
    return fig, (
        f"FastQC-style sequence length line for {read}, from {lengths[0]:.0f} to "
        f"{lengths[-1]:.0f} bp; mode {mode:.0f} bp."
    )


METRICS: tuple[MetricSpec, ...] = (
    MetricSpec("per_base_quality", "per_cycle", "per_base_quality", "Per base sequence quality", ("cycle", "mean_quality"), plot_per_base_quality),
    MetricSpec("adapter_content", "adapter_content", "adapter_content", "Adapter content", ("pos",), plot_adapter_content, adapters_only=True, variable_series=True),
    MetricSpec("per_base_sequence_content", "per_base_sequence_content", "per_base_sequence_content", "Per base sequence content", ("position", "A", "C", "G", "T", "N"), plot_base_content),
    MetricSpec("per_sequence_gc_content", "per_sequence_gc_content", "per_sequence_gc_content", "Per sequence GC content", ("gc_percent", "reads"), plot_gc_content),
    MetricSpec("per_base_n_content", "per_base_n_content", "per_base_n_content", "Per base N content", ("position", "N_percent"), plot_n_content),
    MetricSpec("sequence_length_distribution", "sequence_length_distribution", "sequence_length_distribution", "Sequence length distribution", ("length", "reads"), plot_length_distribution),
    MetricSpec("per_sequence_quality", "per_sequence_quality", "per_sequence_quality", "Per sequence quality scores", ("mean_quality", "read_count"), plot_per_sequence_quality),
)

RETIRED_PLOT_PREFIXES = ("quality_distribution",)


def _save_figure(fig: plt.Figure, output_dir: Path, spec: MetricSpec, read: str, formats: Sequence[str]) -> Mapping[str, str]:
    files: dict[str, str] = {}
    metadata = {"Creator": "NeoQC", "Title": f"{spec.title} — {read}"}
    for output_format in formats:
        filename = spec.output_name(read, output_format)
        path = output_dir / filename
        if output_format == "png":
            fig.savefig(path, format="png", dpi=PNG_DPI, metadata=metadata)
        else:
            fig.savefig(path, format="svg", metadata=metadata)
        files[output_format] = filename
    return files


def generate_plots(
    input_dir: Path,
    output_dir: Path,
    *,
    include_adapters: bool = True,
    formats: Sequence[str] = ("svg", "png"),
) -> dict[str, object]:
    """Render all recognized NeoQC TSV files and return the report contract."""

    apply_theme()
    input_dir = input_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    for prefix in RETIRED_PLOT_PREFIXES:
        for read in ("R1", "R2"):
            for suffix in ("svg", "png"):
                retired_output = output_dir / f"{prefix}_{read}.{suffix}"
                if retired_output.is_file():
                    retired_output.unlink()
    formats = tuple(dict.fromkeys(formats))
    invalid_formats = sorted(set(formats) - {"svg", "png"})
    if invalid_formats:
        raise ValueError(f"unsupported output format(s): {', '.join(invalid_formats)}")
    if not formats:
        raise ValueError("at least one output format is required")

    recognized_sources = {
        spec.source_name(read)
        for spec in METRICS
        for read in ("R1", "R2")
    }
    if not any((input_dir / name).is_file() for name in recognized_sources):
        raise PlotDataError(f"no recognized NeoQC TSV files found in {input_dir}")

    entries: list[dict[str, object]] = []
    for read in ("R1", "R2"):
        for spec in METRICS:
            source_name = spec.source_name(read)
            source = input_dir / source_name
            entry: dict[str, object] = {
                "id": spec.metric_id,
                "read": read,
                "title": spec.title,
                "source": source_name,
                "status": "skipped",
            }
            # A rerun with fewer reads or disabled adapters must not leave stale
            # report assets that contradict the new manifest.
            for suffix in ("svg", "png"):
                stale_output = output_dir / spec.output_name(read, suffix)
                if stale_output.is_file():
                    stale_output.unlink()
            if spec.adapters_only and not include_adapters:
                entry["reason"] = "adapter_analysis_disabled"
                entries.append(entry)
                continue
            if not source.is_file():
                entry["reason"] = "source_not_found"
                entries.append(entry)
                continue

            try:
                rows = _read_numeric_tsv(source, spec)
                fig, alt_text = spec.plot(rows, read)
                try:
                    files = _save_figure(fig, output_dir, spec, read, formats)
                finally:
                    plt.close(fig)
                entry.update(files)
                entry["alt_text"] = alt_text
                entry["status"] = "generated"
            except Exception as error:  # retain failures in the machine-readable contract
                entry["status"] = "error"
                entry["reason"] = str(error)
            entries.append(entry)

    generated = sum(entry["status"] == "generated" for entry in entries)
    errors = sum(entry["status"] == "error" for entry in entries)
    manifest: dict[str, object] = {
        "schema_version": 1,
        "theme": "neo-report",
        "formats": list(formats),
        "figure": {
            "width_inches": FIGURE_SIZE[0],
            "height_inches": FIGURE_SIZE[1],
            "aspect_ratio": "16:9",
            "png_width_px": round(FIGURE_SIZE[0] * PNG_DPI),
            "png_height_px": round(FIGURE_SIZE[1] * PNG_DPI),
            "png_dpi": PNG_DPI,
        },
        "summary": {"generated": generated, "errors": errors, "skipped": len(entries) - generated - errors},
        "plots": entries,
    }
    manifest_path = output_dir / "plots_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    return manifest
