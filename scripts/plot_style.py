"""Shared visual language for every NeoQC quality-control chart."""

from __future__ import annotations

import re
from pathlib import Path

from matplotlib import font_manager, rcParams
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from matplotlib.ticker import FuncFormatter, MaxNLocator

from neoqc_theme import LIGHT


INK = LIGHT["--c-ink"]
MUTED = LIGHT["--c-muted"]
BRAND = LIGHT["--c-brand"]
BRAND_DARK = LIGHT["--c-brand-dark"]
ACCENT = LIGHT["--c-accent"]
PANEL = LIGHT["--c-panel"]
GRID = LIGHT["--c-grid"]
SPINE = LIGHT["--c-spine"]
WHITE = LIGHT["--c-bg"]
WARNING = LIGHT["--c-warning"]
DANGER = LIGHT["--c-danger"]
NEUTRAL = LIGHT["--c-n"]

SERIES_COLORS = (BRAND, ACCENT, WARNING, DANGER, LIGHT["--c-violet"], LIGHT["--c-teal"])
# Classic genome-browser nucleotide colours: A green, C blue, G amber, T red.
BASE_COLORS = {
    "A": LIGHT["--c-a"],
    "C": LIGHT["--c-c"],
    "G": WARNING,
    "T": DANGER,
    "N": NEUTRAL,
}
LINE_STYLES = ("-", "--", "-.", ":")
MARKERS = ("o", "s", "^", "D", "P", "X")

FIGURE_SIZE = (8.0, 4.5)
PNG_DPI = 300

FONT_DIR = Path(__file__).resolve().parents[1] / "assets" / "fonts" / "ibm-plex"
FONT_FILES = ("IBMPlexSans-Regular.ttf", "IBMPlexSans-Bold.ttf", "IBMPlexMono-Regular.ttf")
SANS_FALLBACK = "'Segoe UI', Arial, sans-serif"
MONO_FALLBACK = "Consolas, monospace"


def _register_fonts() -> tuple[str, str]:
    """Register bundled IBM Plex fonts; fall back to DejaVu when they are absent."""

    registered = set()
    for name in FONT_FILES:
        path = FONT_DIR / name
        if path.is_file():
            font_manager.fontManager.addfont(str(path))
            registered.add(name)
    sans = "IBM Plex Sans" if {"IBMPlexSans-Regular.ttf", "IBMPlexSans-Bold.ttf"} <= registered else "DejaVu Sans"
    mono = "IBM Plex Mono" if "IBMPlexMono-Regular.ttf" in registered else "DejaVu Sans Mono"
    return sans, mono


SANS_FONT, MONO_FONT = _register_fonts()


def apply_theme() -> None:
    """Install deterministic, report-aligned Matplotlib defaults."""

    rcParams.update(
        {
            "font.family": SANS_FONT,
            "font.size": 11,
            "text.color": INK,
            "axes.titlesize": 15,
            "axes.titleweight": 700,
            "axes.labelsize": 12,
            "axes.labelcolor": INK,
            "axes.edgecolor": SPINE,
            "axes.linewidth": 0.8,
            "axes.facecolor": WHITE,
            "figure.facecolor": WHITE,
            "xtick.color": MUTED,
            "ytick.color": MUTED,
            "xtick.labelsize": 10,
            "ytick.labelsize": 10,
            "grid.color": GRID,
            "grid.linewidth": 0.7,
            "grid.alpha": 1.0,
            "legend.frameon": False,
            "legend.fontsize": 10,
            "legend.labelcolor": INK,
            "lines.linewidth": 1.8,
            "savefig.facecolor": WHITE,
            "savefig.edgecolor": WHITE,
            # Keep text as text: the report supplies the embedded web fonts and
            # translated SVG remains searchable and selectable.
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
    """Common axes styling.

    The chart title and read are shown by the report panel around the chart, so
    they are not repeated inside the figure; file metadata still carries both.
    """

    del title, read
    ax.set_xlabel(xlabel, labelpad=8)
    ax.set_ylabel(ylabel, labelpad=8)
    ax.grid(axis="both", zorder=0)
    ax.set_axisbelow(True)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(length=0, pad=5)
    try:  # Matplotlib >= 3.8; finish_figure() covers older releases
        ax.tick_params(labelfontfamily=MONO_FONT)
    except (TypeError, ValueError, AttributeError):
        pass
    ax.xaxis.set_major_locator(MaxNLocator(nbins=10, integer=True))


def use_compact_y_axis(ax: Axes) -> None:
    ax.yaxis.set_major_formatter(FuncFormatter(compact_number))


def finish_figure(fig: Figure) -> None:
    for ax in fig.axes:
        for label in ax.get_xticklabels() + ax.get_yticklabels():
            label.set_fontfamily(MONO_FONT)
    fig.tight_layout(pad=1.6)


_FONT_FAMILY = re.compile(r"font-family: '?(IBM Plex Sans|IBM Plex Mono|DejaVu Sans Mono|DejaVu Sans)'?")


def add_svg_font_fallbacks(svg: str) -> str:
    """Give SVG text a generic fallback when the chart is opened outside the report."""

    def fallback(match: re.Match[str]) -> str:
        family = match.group(1)
        generic = MONO_FALLBACK if "Mono" in family else SANS_FALLBACK
        return f"font-family: '{family}', {generic}"

    return _FONT_FAMILY.sub(fallback, svg)
