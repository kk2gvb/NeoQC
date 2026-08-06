#!/usr/bin/env python3
"""Evaluate NeoQC TSV metrics with a versioned PASS/WARNING/FAIL ruleset."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from qc_rules import DEFAULT_RULESET, QcRuleError, write_evaluation


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input_dir", type=Path, help="directory containing NeoQC TSV files")
    parser.add_argument("--output", type=Path, help="qc_evaluation.json output path")
    parser.add_argument(
        "--ruleset",
        type=Path,
        default=DEFAULT_RULESET,
        help="versioned JSON ruleset",
    )
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    if not args.input_dir.is_dir():
        print(f"QC evaluation error: input directory does not exist: {args.input_dir}", file=sys.stderr)
        return 2
    try:
        output = write_evaluation(args.input_dir, args.output, args.ruleset)
    except (OSError, ValueError, QcRuleError) as error:
        print(f"QC evaluation error: {error}", file=sys.stderr)
        return 2
    print(f"QC evaluation: {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
