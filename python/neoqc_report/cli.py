"""Command-line boundary for generating a report from JSON."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Sequence

from .html_report import write_html_report
from .models import ReportData, ReportValidationError


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Generate a standalone NeoQC HTML report from JSON data."
    )
    parser.add_argument("--input", required=True, type=Path, help="Input UTF-8 JSON file")
    parser.add_argument("--output", required=True, type=Path, help="Output HTML file")
    parser.add_argument(
        "--no-print-button", action="store_true", help="Do not include the print/PDF button"
    )
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        with args.input.open(encoding="utf-8") as source:
            raw_data = json.load(source)
        report = ReportData.from_dict(raw_data)
        output = write_html_report(
            report,
            args.output,
            include_print_button=not args.no_print_button,
        )
    except (OSError, json.JSONDecodeError, ReportValidationError) as error:
        print(f"HTML report error: {error}", file=sys.stderr)
        return 1
    print(f"HTML report: {output}")
    return 0
