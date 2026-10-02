"""Colour tokens shared by NeoQC charts and the HTML report.

Charts are drawn with the light values. The HTML report inlines chart SVG and
replaces every light value with its CSS variable, so the same chart follows the
report's light or dark theme. Each light value must therefore be unique.
This module must stay free of Matplotlib so the report can be regenerated
without the plotting stack.
"""

from __future__ import annotations

# CSS variable -> (light, dark)
CHART_COLORS: dict[str, tuple[str, str]] = {
    "--c-bg": ("#ffffff", "#151b23"),
    "--c-ink": ("#1b2430", "#e3e8ef"),
    "--c-muted": ("#5d6878", "#9aa6b6"),
    "--c-grid": ("#e4e8ee", "#262f3b"),
    "--c-spine": ("#c9d0d9", "#3a4555"),
    "--c-panel": ("#f2f4f7", "#1d2530"),
    "--c-brand": ("#1f5fd1", "#6c9bff"),
    "--c-brand-dark": ("#163f8c", "#8fb2ff"),
    "--c-accent": ("#0e9f6e", "#2cc28f"),
    "--c-warning": ("#e09b00", "#f0b23a"),
    "--c-danger": ("#d8343d", "#ff646b"),
    "--c-a": ("#12a150", "#3ccf7c"),
    "--c-c": ("#1f6fd1", "#58a0ff"),
    "--c-n": ("#8a94a3", "#7d8898"),
    "--c-violet": ("#7d56c2", "#a889f0"),
    "--c-teal": ("#138a9e", "#37b9cc"),
}

LIGHT = {name: light for name, (light, _) in CHART_COLORS.items()}

# Matplotlib writes colours as lower-case #rrggbb.
LIGHT_TO_VAR = {light.lower(): name for name, (light, _) in CHART_COLORS.items()}

if len(LIGHT_TO_VAR) != len(CHART_COLORS):  # pragma: no cover - guarded at import time
    raise RuntimeError("chart colour tokens must have unique light values")
