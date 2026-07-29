#!/usr/bin/env python3
"""Benchmark NeoQC on one single-end or paired-end FASTQ sample.

The script runs NeoQC several times, samples its resident memory with psutil,
and stores a reproducible, human-readable report in benchmarks/<sample>.txt.
"""

from __future__ import annotations

import argparse
import os
import platform
import re
import statistics
import subprocess
import sys
import tempfile
import threading
import time
from dataclasses import dataclass
from pathlib import Path

try:
    import psutil
except ImportError as exc:  # pragma: no cover - depends on local environment
    raise SystemExit("psutil is required: install it with 'python3 -m pip install psutil'.") from exc


SAMPLE_INTERVAL_SECONDS = 0.02


@dataclass
class Summary:
    reads: int
    min_length: int
    max_length: int
    avg_length: float


@dataclass
class RunResult:
    seconds: float
    peak_rss_bytes: int
    r1: Summary
    r2: Summary | None


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--r1", required=True, type=Path, help="R1 FASTQ file")
    parser.add_argument("--r2", type=Path, help="R2 FASTQ file (paired-end)")
    parser.add_argument("--sample-id", required=True, help="Sample identifier")
    parser.add_argument(
        "--build-dir", type=Path, default=Path("build"),
        help="Directory containing the neoqc executable (default: build)",
    )
    parser.add_argument(
        "--runs", type=int, default=3,
        help="Number of benchmark runs (default: 3)",
    )
    parser.add_argument("--plot", action="store_true", help="Pass --plot to NeoQC")
    args = parser.parse_args()
    if args.runs < 1:
        parser.error("--runs must be at least 1")
    for option, path in (("--r1", args.r1), ("--r2", args.r2)):
        if path is not None and not path.is_file():
            parser.error(f"{option}: file not found: {path}")
    return args


def find_neoqc(build_dir: Path) -> Path:
    executable = build_dir / ("neoqc.exe" if os.name == "nt" else "neoqc")
    if not executable.is_file():
        raise FileNotFoundError(
            f"NeoQC executable not found: {executable}. Build the project or pass --build-dir."
        )
    return executable.resolve()


def parse_summary(path: Path) -> Summary:
    """Read the fields NeoQC emits in *_summary.txt."""
    try:
        text = path.read_text(encoding="utf-8")
    except FileNotFoundError as exc:
        raise RuntimeError(f"NeoQC did not create expected summary: {path}") from exc

    fields = {
        "reads": r"^Processed reads\s*:\s*(\d+)\s*$",
        "min_length": r"^Min length\s*:\s*(\d+)\s*$",
        "max_length": r"^Max length\s*:\s*(\d+)\s*$",
        "avg_length": r"^Avg length\s*:\s*([0-9]+(?:\.[0-9]+)?)\s*$",
    }
    values: dict[str, str] = {}
    for line in text.splitlines():
        for name, pattern in fields.items():
            match = re.match(pattern, line)
            if match:
                values[name] = match.group(1)

    missing = set(fields) - set(values)
    if missing:
        raise RuntimeError(f"Cannot parse {path}; missing: {', '.join(sorted(missing))}")
    return Summary(
        reads=int(values["reads"]), min_length=int(values["min_length"]),
        max_length=int(values["max_length"]), avg_length=float(values["avg_length"]),
    )


def process_tree_rss(process: psutil.Process) -> int:
    """Return RSS of NeoQC and any currently running child processes."""
    total = 0
    try:
        candidates = [process, *process.children(recursive=True)]
    except (psutil.NoSuchProcess, psutil.AccessDenied):
        return total
    for candidate in candidates:
        try:
            total += candidate.memory_info().rss
        except (psutil.NoSuchProcess, psutil.AccessDenied):
            pass
    return total


def run_neoqc(command: list[str], output_dir: Path, sample_id: str, paired: bool) -> RunResult:
    peak_rss = 0
    stop_sampling = threading.Event()
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    root_process = psutil.Process(process.pid)

    def sample_memory() -> None:
        nonlocal peak_rss
        while not stop_sampling.is_set():
            peak_rss = max(peak_rss, process_tree_rss(root_process))
            stop_sampling.wait(SAMPLE_INTERVAL_SECONDS)
        peak_rss = max(peak_rss, process_tree_rss(root_process))

    sampler = threading.Thread(target=sample_memory, daemon=True)
    started = time.perf_counter()
    sampler.start()
    stdout, stderr = process.communicate()
    elapsed = time.perf_counter() - started
    stop_sampling.set()
    sampler.join()

    if process.returncode:
        details = stderr.strip() or stdout.strip() or "no output"
        raise RuntimeError(f"NeoQC exited with code {process.returncode}: {details}")

    r1 = parse_summary(output_dir / f"{sample_id}_R1_summary.txt")
    r2 = parse_summary(output_dir / f"{sample_id}_R2_summary.txt") if paired else None
    return RunResult(elapsed, peak_rss, r1, r2)


def file_size(path: Path) -> str:
    size = path.stat().st_size
    units = ("B", "KiB", "MiB", "GiB", "TiB")
    for unit in units:
        if size < 1024 or unit == units[-1]:
            return f"{size:.2f} {unit}" if unit != "B" else f"{size} B"
        size /= 1024
    raise AssertionError("unreachable")


def bytes_to_mib(value: float) -> str:
    return f"{value / 1024 / 1024:.2f} MiB"


def disk_models() -> str:
    if sys.platform.startswith("linux"):
        models: list[str] = []
        for device in sorted(Path("/sys/block").glob("*")):
            model_file = device / "device" / "model"
            try:
                model = model_file.read_text().strip()
            except OSError:
                continue
            if model and model not in models:
                models.append(model)
        return ", ".join(models) if models else "unknown"
    return "unavailable on this OS"


def system_info() -> list[tuple[str, str]]:
    memory = psutil.virtual_memory()
    return [
        ("CPU", platform.processor() or platform.machine() or "unknown"),
        ("CPU logical cores", str(psutil.cpu_count(logical=True) or "unknown")),
        ("RAM", bytes_to_mib(memory.total)),
        ("OS", f"{platform.system()} {platform.release()} ({platform.machine()})"),
        ("Storage model", disk_models()),
    ]


def render_report(args: argparse.Namespace, results: list[RunResult]) -> str:
    first = results[0]
    times = [item.seconds for item in results]
    memory = [item.peak_rss_bytes for item in results]
    pairs = first.r1.reads
    if first.r2 is not None and first.r2.reads != pairs:
        raise RuntimeError("R1 and R2 have different read counts in the benchmark summary")
    throughput = [pairs / seconds if seconds else 0.0 for seconds in times]
    throughput_unit = "pairs/s" if first.r2 else "reads/s"

    lines = ["NeoQC benchmark", "=" * 68, "", "Input"]
    lines += [f"  Sample ID: {args.sample_id}", f"  R1: {args.r1} ({file_size(args.r1)})"]
    if args.r2:
        lines.append(f"  R2: {args.r2} ({file_size(args.r2)})")
    lines += [f"  Layout: {'paired-end' if args.r2 else 'single-end'}", "", "Read statistics"]
    lines += [
        f"  R1 reads: {first.r1.reads:,}",
        f"  R1 length: {first.r1.min_length}-{first.r1.max_length} bp (mean {first.r1.avg_length:.2f} bp)",
    ]
    if first.r2:
        lines += [
            f"  R2 reads: {first.r2.reads:,}",
            f"  R2 length: {first.r2.min_length}-{first.r2.max_length} bp (mean {first.r2.avg_length:.2f} bp)",
            f"  Read pairs: {pairs:,}",
        ]
    else:
        lines.append(f"  Reads: {pairs:,}")

    lines += ["", "Performance"]
    for index, item in enumerate(results, start=1):
        lines.append(
            f"  Run {index}: {item.seconds:.3f} s | peak RSS {bytes_to_mib(item.peak_rss_bytes)} | {throughput[index - 1]:,.0f} {throughput_unit}"
        )
    lines += [
        f"  Time, median: {statistics.median(times):.3f} s",
        f"  Time, mean:   {statistics.mean(times):.3f} s",
        f"  Peak RSS, max:  {bytes_to_mib(max(memory))}",
        f"  Throughput, median: {statistics.median(throughput):,.0f} {throughput_unit}",
        f"  Throughput, mean:   {statistics.mean(throughput):,.0f} {throughput_unit}",
        "", "System",
    ]
    lines += [f"  {name}: {value}" for name, value in system_info()]
    return "\n".join(lines) + "\n"


def main() -> int:
    args = parse_args()
    executable = find_neoqc(args.build_dir)
    results: list[RunResult] = []
    with tempfile.TemporaryDirectory(prefix="neoqc-benchmark-") as temporary:
        temporary_dir = Path(temporary)
        for run in range(1, args.runs + 1):
            output_dir = temporary_dir / f"run-{run}"
            command = [str(executable), "--r1", str(args.r1.resolve()), "--sample-id", args.sample_id, "--out", str(output_dir)]
            if args.r2:
                command.extend(["--r2", str(args.r2.resolve())])
            if args.plot:
                command.append("--plot")
            print(f"Running NeoQC: {run}/{args.runs}", file=sys.stderr)
            results.append(run_neoqc(command, output_dir, args.sample_id, args.r2 is not None))

    report = render_report(args, results)
    report_dir = Path("benchmarks")
    report_dir.mkdir(parents=True, exist_ok=True)
    report_path = report_dir / f"{args.sample_id}.txt"
    report_path.write_text(report, encoding="utf-8")
    print(report, end="")
    print(f"Report saved to {report_path}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (FileNotFoundError, RuntimeError) as exc:
        print(f"benchmark.py: error: {exc}", file=sys.stderr)
        raise SystemExit(1) from exc
