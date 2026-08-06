"""Shared visual language for every NeoQC quality-control chart."""

from __future__ import annotations

from matplotlib import rcParams
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from matplotlib.ticker import FuncFormatter, MaxNLocator


INK = "#0A132D"
MUTED = "#657285"
BRAND = "#2947A0"
BRAND_DARK = "#192F70"
NAVY = "#031E4F"
ACCENT = "#539D96"
PANEL = "#F4F7F9"
GRID = "#DCE3EA"
WHITE = "#FFFFFF"
WARNING = "#E69F00"
DANGER = "#C84A5A"

SERIES_COLORS = (BRAND, ACCENT, WARNING, DANGER, "#7A5AA6", "#2A839A")
BASE_COLORS = {
    "A": BRAND,
    "C": ACCENT,
    "G": WARNING,
    "T": "#D55E00",
    "N": MUTED,
}
LINE_STYLES = ("-", "--", "-.", ":")
MARKERS = ("o", "s", "^", "D", "P", "X")

FIGURE_SIZE = (8.0, 4.5)
PNG_DPI = 300


def apply_theme() -> None:
    """Install deterministic, report-aligned Matplotlib defaults."""

    rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 9,
            "axes.titlesize": 14,
            "axes.titleweight": 700,
            "axes.labelsize": 10,
            "axes.labelcolor": INK,
            "axes.edgecolor": GRID,
            "axes.linewidth": 0.8,
            "axes.facecolor": WHITE,
            "figure.facecolor": WHITE,
            "xtick.color": MUTED,
            "ytick.color": MUTED,
            "xtick.labelsize": 8,
            "ytick.labelsize": 8,
            "grid.color": GRID,
            "grid.linewidth": 0.7,
            "grid.alpha": 0.75,
            "legend.frameon": False,
            "legend.fontsize": 8,
            "lines.linewidth": 2.0,
            "savefig.facecolor": WHITE,
            "savefig.edgecolor": WHITE,
            "svg.fonttype": "none",
        }
    )


def compact_number(value: float, _position: float | None = None) -> str:
    absolute = abs(value)
    if absolute >= 1_000_000_000:
        return f"{value / 1_000_000_000:.1f}B"
    if absolute >= 1_000_000:
        return f"{value / 1_000_000:.1f}M"
    if absolute >= 1_000:
        return f"{value / 1_000:.1f}k"
    return f"{value:g}"


def setup_axes(ax: Axes, title: str, xlabel: str, ylabel: str, read: str) -> None:
    ax.set_title(title, loc="left", color=INK, pad=14)
    ax.set_xlabel(xlabel, labelpad=8)
    ax.set_ylabel(ylabel, labelpad=8)
    ax.grid(axis="y", zorder=0)
    ax.grid(axis="x", visible=False)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(length=0, pad=5)
    ax.xaxis.set_major_locator(MaxNLocator(nbins=10, integer=True))
    ax.text(
        1.0,
        1.025,
        f"NeoQC  •  {read}",
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        color=MUTED,
        fontsize=7.5,
        fontweight=700,
    )


def use_compact_y_axis(ax: Axes) -> None:
    ax.yaxis.set_major_formatter(FuncFormatter(compact_number))


def finish_figure(fig: Figure) -> None:
    fig.tight_layout(pad=1.6)
