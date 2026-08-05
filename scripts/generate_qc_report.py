#!/usr/bin/env python3
"""Build a compact self-contained NeoQC HTML report from plot artifacts."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from qc_report import QcReportError, generate_qc_report


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("result_dir", type=Path, help="NeoQC sample result directory")
    parser.add_argument("--plot-dir", type=Path, help="plot directory (default: <result>/plots)")
    parser.add_argument("--output", type=Path, help="HTML output path")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        output = generate_qc_report(args.result_dir, args.plot_dir, args.output)
    except (OSError, QcReportError) as error:
        print(f"report error: {error}", file=sys.stderr)
        return 2
    print(f"QC report: {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
