#!/usr/bin/env python3
"""Generate report-ready NeoQC plots for plain QC or trimming results."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path


TRIMMING_STAGES = ("before", "after")


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "input_dir",
        type=Path,
        help=(
            "directory containing NeoQC TSV files, or a trimming result root "
            "containing qc/before and qc/after"
        ),
    )
    parser.add_argument(
        "output_dir",
        type=Path,
        help=(
            "directory for charts and plots_manifest.json; trimming root mode "
            "uses before/ and after/ subdirectories"
        ),
    )
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
    parser.add_argument(
        "--ruleset",
        type=Path,
        help="versioned QC ruleset JSON (default: fastqc-compatible-v1)",
    )
    parser.add_argument(
        "--run-manifest",
        type=Path,
        help="override the automatically selected run_manifest.json",
    )
    return parser.parse_args(argv)


def _run_manifest_path(input_dir: Path, explicit_path: Path | None) -> Path:
    """Resolve the publication manifest for plain and nested trimming layouts."""
    if explicit_path is not None:
        return explicit_path
    if input_dir.name in TRIMMING_STAGES and input_dir.parent.name == "qc":
        return input_dir.parent.parent / "run_manifest.json"
    return input_dir / "run_manifest.json"


def _plot_jobs(
    input_dir: Path,
    output_dir: Path,
    explicit_manifest: Path | None,
) -> list[tuple[Path, Path, Path]]:
    """Return one plain-QC job or two jobs for a trimming result root."""
    trimming_qc_dir = input_dir / "qc"
    trimming_inputs = [trimming_qc_dir / stage for stage in TRIMMING_STAGES]
    trimming_stages_present = [path.is_dir() for path in trimming_inputs]
    if any(trimming_stages_present) and not all(trimming_stages_present):
        missing = [
            stage
            for stage, present in zip(TRIMMING_STAGES, trimming_stages_present)
            if not present
        ]
        raise ValueError(
            "incomplete trimming QC layout; missing: "
            + ", ".join(f"qc/{stage}" for stage in missing)
        )
    if all(trimming_stages_present):
        manifest_path = explicit_manifest or input_dir / "run_manifest.json"
        return [
            (stage_input, output_dir / stage, manifest_path)
            for stage, stage_input in zip(TRIMMING_STAGES, trimming_inputs)
        ]
    return [
        (
            input_dir,
            output_dir,
            _run_manifest_path(input_dir, explicit_manifest),
        )
    ]


def _plot_one(
    args: argparse.Namespace,
    input_dir: Path,
    output_dir: Path,
    run_manifest_path: Path,
) -> int:
    try:
        from plot_metrics import generate_plots
    except ImportError as error:
        print(f"plot error: plotting dependency unavailable: {error}", file=sys.stderr)
        return 2

    print(f"NeoQC plot input: {input_dir}")

    try:
        manifest = generate_plots(
            input_dir,
            output_dir,
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
    print(f"Manifest: {output_dir / 'plots_manifest.json'}")
    try:
        from qc_rules import DEFAULT_RULESET, write_evaluation

        evaluation_path = write_evaluation(
            input_dir,
            ruleset_path=args.ruleset or DEFAULT_RULESET,
            run_manifest_path=run_manifest_path,
        )
    except (ImportError, OSError, ValueError) as error:
        print(f"QC evaluation error: {error}", file=sys.stderr)
        return 2
    print(f"QC evaluation: {evaluation_path}")
    if not args.no_html_report:
        try:
            from qc_report import generate_qc_report

            report_path = generate_qc_report(input_dir, output_dir)
        except (ImportError, OSError, ValueError) as error:
            print(f"report error: {error}", file=sys.stderr)
            return 2
        print(f"QC report: {report_path}")
    if args.strict and summary["errors"]:
        return 1
    return 0


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    if not args.input_dir.is_dir():
        print(f"plot error: input directory does not exist: {args.input_dir}", file=sys.stderr)
        return 2

    try:
        jobs = _plot_jobs(args.input_dir, args.output_dir, args.run_manifest)
    except ValueError as error:
        print(f"plot error: {error}", file=sys.stderr)
        return 2
    return_code = 0
    for input_dir, output_dir, run_manifest_path in jobs:
        return_code = max(
            return_code,
            _plot_one(args, input_dir, output_dir, run_manifest_path),
        )
    return return_code


if __name__ == "__main__":
    raise SystemExit(main())
