#!/usr/bin/env python3
"""Generate a complete, report-ready NeoQC plot set."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input_dir", type=Path, help="directory containing NeoQC TSV files")
    parser.add_argument("output_dir", type=Path, help="directory for charts and plots_manifest.json")
    parser.add_argument(
        "--formats",
        nargs="+",
        choices=("svg", "png"),
        default=("svg", "png"),
        help="output formats (default: svg png)",
    )
    parser.add_argument(
        "--skip-adapters",
        action="store_true",
        help="record adapter charts as skipped even if adapter TSV files exist",
    )
    parser.add_argument(
        "--strict",
        action="store_true",
        help="return a failure status when any present TSV cannot be rendered",
    )
    parser.add_argument(
        "--no-html-report",
        action="store_true",
        help="generate chart assets only, without neoqc_qc_report.html",
    )
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    if not args.input_dir.is_dir():
        print(f"plot error: input directory does not exist: {args.input_dir}", file=sys.stderr)
        return 2

    try:
        from plot_metrics import generate_plots
    except ImportError as error:
        print(f"plot error: plotting dependency unavailable: {error}", file=sys.stderr)
        return 2

    try:
        manifest = generate_plots(
            args.input_dir,
            args.output_dir,
            include_adapters=not args.skip_adapters,
            formats=args.formats,
        )
    except (OSError, ValueError) as error:
        print(f"plot error: {error}", file=sys.stderr)
        return 2

    summary = manifest["summary"]
    print(
        "NeoQC plots: "
        f"generated={summary['generated']}, "
        f"skipped={summary['skipped']}, "
        f"errors={summary['errors']}"
    )
    print(f"Manifest: {args.output_dir / 'plots_manifest.json'}")
    if not args.no_html_report:
        try:
            from qc_report import generate_qc_report

            report_path = generate_qc_report(args.input_dir, args.output_dir)
        except (ImportError, OSError, ValueError) as error:
            print(f"report error: {error}", file=sys.stderr)
            return 2
        print(f"QC report: {report_path}")
    if args.strict and summary["errors"]:
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
