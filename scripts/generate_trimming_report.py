#!/usr/bin/env python3
"""Build a self-contained NeoQC before/after trimming HTML report."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from qc_report import QcReportError
from trimming_report import generate_trimming_report


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("result_dir", type=Path, help="NeoQC trimming result directory")
    parser.add_argument("--output", type=Path, help="HTML output path")
    args = parser.parse_args(argv)
    try:
        output = generate_trimming_report(args.result_dir, args.output)
    except (OSError, QcReportError) as error:
        print(f"trimming report error: {error}", file=sys.stderr)
        return 2
    print(f"Trimming report: {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
