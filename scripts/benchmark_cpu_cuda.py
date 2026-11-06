#!/usr/bin/env python3
"""Reproducible, sequential CPU/CUDA benchmark runner for NeoQC.

The runner deliberately treats CUDA fallback as a non-successful CUDA benchmark.
It never removes input FASTQ files and keeps each NeoQC invocation in its own
directory so that an interrupted benchmark can be resumed safely.
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import gzip
import json
import os
import platform
import re
import shutil
import signal
import subprocess
import sys
import time
from dataclasses import asdict, dataclass
from pathlib import Path
from statistics import median
from typing import Any, Iterable


DEFAULT_CONFIG = Path(__file__).with_name("benchmark_config.json")
DEFAULT_EXCLUDED = {"run_manifest.json"}
CANCELLED = False

RESULT_FIELDS = [
    "run_key", "sample", "material", "assay", "backend", "mode", "repeat",
    "status", "comparison_status", "fallback", "exit_code", "wall_time_s",
    "reads", "bases", "read_fastq_s", "parallel_processing_s", "write_fastq_s",
    "merge_analyzers_s", "statistics_s", "reports_s", "cuda_packing_s", "h2d_s",
    "kernel_s", "d2h_s", "gpu_qc_total_s", "cpu_remaining_qc_s", "peak_rss_gib",
    "peak_gpu_memory_mib", "output_dir", "stdout", "stderr", "comparison_file",
]

TIMERS = {
    "Read FASTQ": "read_fastq_s",
    "Parallel processing": "parallel_processing_s",
    "Write FASTQ": "write_fastq_s",
    "Merge analyzers": "merge_analyzers_s",
    "Statistics": "statistics_s",
    "Reports": "reports_s",
    "Packing": "cuda_packing_s",
    "Host to Device": "h2d_s",
    "Kernel": "kernel_s",
    "Device to Host": "d2h_s",
    "GPU QC total": "gpu_qc_total_s",
    "CPU remaining QC": "cpu_remaining_qc_s",
}


@dataclass(frozen=True)
class Dataset:
    sample: str
    material: str
    assay: str
    r1: Path
    r2: Path

    @property
    def dataset_id(self) -> str:
        return f"{self.sample}/{self.assay}"

    @property
    def safe_id(self) -> str:
        return f"{self.sample}_{self.assay}"


def die(message: str) -> None:
    raise RuntimeError(message)


def path_is_within(path: Path, parent: Path) -> bool:
    try:
        path.resolve().relative_to(parent.resolve())
        return True
    except ValueError:
        return False


def ensure_safe_output(output: Path, data: Path) -> None:
    resolved = output.resolve()
    data_resolved = data.resolve()
    if resolved == data_resolved or path_is_within(resolved, data_resolved):
        die(f"Refusing output inside input data directory: {resolved}")
    if resolved == Path(resolved.anchor):
        die("Refusing filesystem root as benchmark output")


def load_config(path: Path) -> dict[str, Any]:
    try:
        with path.open(encoding="utf-8") as handle:
            config = json.load(handle)
    except (OSError, json.JSONDecodeError) as exc:
        die(f"Cannot read benchmark config {path}: {exc}")
    if not isinstance(config, dict):
        die("Benchmark config must be a JSON object")
    return config


def pair_for_r1(r1: Path, suffix: str) -> Path:
    if not r1.name.endswith(suffix):
        die(f"Internal discovery error: {r1.name} does not end in {suffix}")
    return r1.with_name(r1.name[: -len(suffix)] + suffix.replace("R1", "R2", 1))


def discover_datasets(data: Path, expected_samples: Iterable[str]) -> list[Dataset]:
    """Discover exactly GED, SOM and RNA pairs for every configured sample."""
    if not data.is_dir():
        die(f"Input data directory does not exist: {data}")

    datasets: list[Dataset] = []
    used_paths: set[Path] = set()
    expected = list(expected_samples)
    for sample in expected:
        sample_root = data / sample
        if not sample_root.is_dir():
            die(f"Missing sample directory: {sample_root}")

        candidates = [
            ("DNA", "GED", sample_root / "dna", "*_GED_R1.fastq.gz", "_R1.fastq.gz"),
            ("DNA", "SOM", sample_root / "dna", "*_SOM_R1.fastq.gz", "_R1.fastq.gz"),
            ("RNA", "RNA", sample_root / "rna", "*_R1_001.fastq.gz", "_R1_001.fastq.gz"),
        ]
        for material, assay, directory, pattern, suffix in candidates:
            if not directory.is_dir():
                die(f"Missing {material} directory for {sample}: {directory}")
            r1_matches = sorted(directory.glob(pattern))
            if len(r1_matches) != 1:
                die(
                    f"Expected exactly one {sample}/{assay} R1 matching "
                    f"{directory / pattern}, found {len(r1_matches)}"
                )
            r1 = r1_matches[0].resolve()
            r2 = pair_for_r1(r1, suffix).resolve()
            if not r2.is_file():
                die(f"Missing mate for {sample}/{assay}: {r2}")
            if r1 == r2 or r1 in used_paths or r2 in used_paths:
                die(f"Duplicate or ambiguous FASTQ mapping for {sample}/{assay}")
            used_paths.update({r1, r2})
            datasets.append(Dataset(sample, material, assay, r1, r2))

    identifiers = [item.dataset_id for item in datasets]
    if len(identifiers) != len(set(identifiers)):
        die("Duplicate dataset identifiers discovered")
    return datasets


def stream_files_equal(left: Path, right: Path, chunk_size: int) -> tuple[bool, dict[str, Any]]:
    """Byte-compare files; gzip inputs are compared after streaming decompression."""
    if not left.exists() or not right.exists():
        return False, {"reason": "missing_file", "left_exists": left.exists(), "right_exists": right.exists()}
    use_gzip = left.suffix == ".gz" or right.suffix == ".gz"
    left_open = gzip.open if use_gzip else open
    right_open = gzip.open if use_gzip else open
    offset = 0
    try:
        with left_open(left, "rb") as lhs, right_open(right, "rb") as rhs:
            while True:
                a = lhs.read(chunk_size)
                b = rhs.read(chunk_size)
                if a == b:
                    if not a:
                        return True, {"method": "gzip-stream" if use_gzip else "byte-stream", "bytes": offset}
                    offset += len(a)
                    continue
                common = min(len(a), len(b))
                first = next((index for index in range(common) if a[index] != b[index]), common)
                return False, {
                    "method": "gzip-stream" if use_gzip else "byte-stream",
                    "reason": "content_differs",
                    "offset": offset + first,
                    "left_byte": a[first] if first < len(a) else None,
                    "right_byte": b[first] if first < len(b) else None,
                }
    except (OSError, EOFError, gzip.BadGzipFile) as exc:
        return False, {"reason": "read_error", "error": str(exc), "offset": offset}


def compare_outputs(cpu_dir: Path, cuda_dir: Path, config: dict[str, Any]) -> dict[str, Any]:
    excluded = set(config.get("excluded_artifacts", DEFAULT_EXCLUDED)) | DEFAULT_EXCLUDED
    chunk_size = int(config.get("comparison_chunk_bytes", 1024 * 1024))
    if chunk_size <= 0:
        die("comparison_chunk_bytes must be positive")

    def files(root: Path) -> set[Path]:
        return {
            path.relative_to(root)
            for path in root.rglob("*")
            if path.is_file() and path.name not in excluded
        }

    if not cpu_dir.is_dir() or not cuda_dir.is_dir():
        return {"status": "NOT_RUN", "reason": "missing_output_directory"}
    cpu_files, cuda_files = files(cpu_dir), files(cuda_dir)
    if cpu_files != cuda_files:
        return {
            "status": "MISMATCH", "reason": "artifact_set_differs",
            "cpu_only": sorted(str(p) for p in cpu_files - cuda_files),
            "cuda_only": sorted(str(p) for p in cuda_files - cpu_files),
        }
    for relative in sorted(cpu_files):
        equal, detail = stream_files_equal(cpu_dir / relative, cuda_dir / relative, chunk_size)
        if not equal:
            return {"status": "MISMATCH", "artifact": str(relative), **detail}
    return {"status": "MATCH", "artifacts": len(cpu_files)}


def parse_metrics(stdout_path: Path, mode: str) -> dict[str, Any]:
    metrics: dict[str, Any] = {}
    if not stdout_path.exists():
        return metrics
    text = stdout_path.read_text(encoding="utf-8", errors="replace")
    for label, key in TIMERS.items():
        match = re.search(rf"(?m)^\s*{re.escape(label)}:\s*([0-9]+(?:\.[0-9]+)?)\s*s\s*$", text)
        if match:
            metrics[key] = float(match.group(1))
    rss = re.search(r"(?m)^\s*Peak RSS \(VmHWM\):\s*([0-9]+(?:\.[0-9]+)?)\s*GiB\s*$", text)
    if rss:
        metrics["peak_rss_gib"] = float(rss.group(1))
    return metrics


def parse_read_counts(output_dir: Path, mode: str) -> tuple[int | str, int | str]:
    report_root = output_dir / "qc" / "before" if mode == "trim" else output_dir
    if not report_root.is_dir():
        return "", ""
    reads = bases = 0
    summaries = sorted(report_root.glob("*_summary.txt"))
    if not summaries:
        return "", ""
    for summary in summaries:
        text = summary.read_text(encoding="utf-8", errors="replace")
        read = re.search(r"(?im)^\s*(?:Processed |Total )?reads\s*:\s*([0-9]+)\s*$", text)
        base = re.search(r"(?im)^\s*Total bases\s*:\s*([0-9]+)\s*$", text)
        if not read or not base:
            return "", ""
        reads += int(read.group(1))
        bases += int(base.group(1))
    return reads, bases


def contains_fallback(*paths: Path) -> bool:
    expression = re.compile(r"CUDA QC batch failed|Falling back to CPU", re.IGNORECASE)
    return any(path.exists() and expression.search(path.read_text(encoding="utf-8", errors="replace")) for path in paths)


def run_command(command: list[str], stdout_path: Path, stderr_path: Path, timeout_seconds: int, env: dict[str, str]) -> tuple[int, float, bool]:
    global CANCELLED
    start = time.perf_counter()
    with stdout_path.open("w", encoding="utf-8") as stdout, stderr_path.open("w", encoding="utf-8") as stderr:
        process = subprocess.Popen(command, stdout=stdout, stderr=stderr, env=env, start_new_session=True)
        timed_out = False
        while True:
            try:
                exit_code = process.wait(timeout=0.1)
                break
            except subprocess.TimeoutExpired:
                if CANCELLED or time.perf_counter() - start >= timeout_seconds:
                    timed_out = not CANCELLED
                    os.killpg(process.pid, signal.SIGTERM)
                    try:
                        exit_code = process.wait(timeout=10)
                    except subprocess.TimeoutExpired:
                        os.killpg(process.pid, signal.SIGKILL)
                        exit_code = process.wait()
                    break
    return exit_code, time.perf_counter() - start, timed_out


def atomic_json(path: Path, value: Any) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    temporary.replace(path)


def write_results(output: Path, records: dict[str, dict[str, Any]]) -> None:
    temporary = output / "results.csv.tmp"
    with temporary.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=RESULT_FIELDS, extrasaction="ignore")
        writer.writeheader()
        for key in sorted(records):
            writer.writerow({field: records[key].get(field, "") for field in RESULT_FIELDS})
    temporary.replace(output / "results.csv")


def load_existing(output: Path) -> dict[str, dict[str, Any]]:
    results = output / "results.csv"
    if not results.is_file():
        return {}
    with results.open(newline="", encoding="utf-8") as handle:
        return {row["run_key"]: row for row in csv.DictReader(handle) if row.get("run_key")}


def command_for(binary: Path, dataset: Dataset, backend: str, mode: str, output: Path) -> list[str]:
    command = [str(binary), "--r1", str(dataset.r1), "--r2", str(dataset.r2), "--sample-id", dataset.safe_id,
               "--out", str(output), "--backend", backend, "--timing"]
    if mode == "trim":
        command.append("--full-trim")
    return command


def make_run_key(dataset: Dataset, mode: str, repeat: int, backend: str) -> str:
    return f"{dataset.safe_id}:{mode}:repeat-{repeat:03d}:{backend}"


def run_one(binary: Path, dataset: Dataset, mode: str, repeat: int, backend: str, root: Path,
            config: dict[str, Any], threads: int | None, dry_run: bool) -> dict[str, Any]:
    run_dir = root / "runs" / dataset.safe_id / mode / f"repeat-{repeat:03d}" / backend
    neoqc_out = run_dir / "neoqc-output"
    run_dir.mkdir(parents=True, exist_ok=True)
    command = command_for(binary, dataset, backend, mode, neoqc_out)
    record: dict[str, Any] = {field: "" for field in RESULT_FIELDS}
    record.update({
        "run_key": make_run_key(dataset, mode, repeat, backend), "sample": dataset.sample,
        "material": dataset.material, "assay": dataset.assay, "backend": backend, "mode": mode,
        "repeat": repeat, "status": "PLANNED" if dry_run else "ERROR", "comparison_status": "PENDING",
        "fallback": False, "exit_code": "", "wall_time_s": "", "reads": "", "bases": "",
        "peak_gpu_memory_mib": "", "output_dir": str(neoqc_out), "stdout": str(run_dir / "stdout.log"),
        "stderr": str(run_dir / "stderr.log"), "comparison_file": "",
    })
    if dry_run:
        atomic_json(run_dir / "command.json", {"command": command, "environment_threads": threads})
        atomic_json(run_dir / "result.json", record)
        return record

    if mode == "trim":
        reserve = float(config.get("minimum_free_disk_gib", 20)) * 1024 ** 3
        free = shutil.disk_usage(root).free
        if free < reserve:
            record.update(status="ERROR", comparison_status="NOT_RUN")
            atomic_json(run_dir / "result.json", record)
            die(f"Only {free / 1024 ** 3:.1f} GiB free; full-trim requires reserve of {reserve / 1024 ** 3:.1f} GiB")

    env = os.environ.copy()
    if threads is not None:
        env["OMP_NUM_THREADS"] = str(threads)
        env["OMP_DYNAMIC"] = "FALSE"
    exit_code, wall_time, timed_out = run_command(
        command, run_dir / "stdout.log", run_dir / "stderr.log", int(config.get("timeout_seconds", 7200)), env)
    record.update(exit_code=exit_code, wall_time_s=round(wall_time, 6))
    record.update(parse_metrics(run_dir / "stdout.log", mode))
    record["reads"], record["bases"] = parse_read_counts(neoqc_out, mode)
    record["fallback"] = backend == "cuda" and contains_fallback(run_dir / "stdout.log", run_dir / "stderr.log")
    if CANCELLED:
        record["status"] = "INTERRUPTED"
        record["comparison_status"] = "NOT_RUN"
    elif timed_out:
        record["status"] = "TIMEOUT"
        record["comparison_status"] = "NOT_RUN"
    elif exit_code == 0:
        record["status"] = "FALLBACK" if record["fallback"] else "OK"
        record["comparison_status"] = "PENDING"
    else:
        record["status"] = "ERROR"
        record["comparison_status"] = "NOT_RUN"
    atomic_json(run_dir / "command.json", {"command": command, "environment_threads": threads})
    atomic_json(run_dir / "result.json", record)
    return record


def compare_pair(records: dict[str, dict[str, Any]], dataset: Dataset, mode: str, repeat: int, root: Path,
                 config: dict[str, Any]) -> None:
    cpu_key, cuda_key = make_run_key(dataset, mode, repeat, "cpu"), make_run_key(dataset, mode, repeat, "cuda")
    cpu, cuda = records.get(cpu_key), records.get(cuda_key)
    if not cpu or not cuda or cpu["status"] not in {"OK", "FALLBACK"} or cuda["status"] not in {"OK", "FALLBACK"}:
        return
    result = compare_outputs(Path(cpu["output_dir"]), Path(cuda["output_dir"]), config)
    comparison_path = root / "runs" / dataset.safe_id / mode / f"repeat-{repeat:03d}" / "comparison.json"
    atomic_json(comparison_path, result)
    for record in (cpu, cuda):
        record["comparison_status"] = result["status"]
        record["comparison_file"] = str(comparison_path)
        if result["status"] == "MISMATCH":
            record["status"] = "MISMATCH"
        run_result = Path(record["stdout"]).parent / "result.json"
        atomic_json(run_result, record)


def format_number(value: Any, digits: int = 3) -> str:
    try:
        return f"{float(value):.{digits}f}"
    except (TypeError, ValueError):
        return "—"


def write_summary(output: Path, records: dict[str, dict[str, Any]]) -> None:
    rows = list(records.values())
    lines = ["# NeoQC CPU/CUDA benchmark", "", "Results are saved after every run. `FALLBACK` means NeoQC reported CPU fallback for at least one CUDA batch; it is not a successful GPU measurement.", "", "## Runs", "", "| Dataset | Type | Mode | Repeat | CPU wall, s | CUDA wall, s | CPU/CUDA | CUDA status | Comparison |", "|---|---|---|---:|---:|---:|---:|---|---|"]
    pairs: dict[tuple[str, str, str], dict[str, dict[str, Any]]] = {}
    for row in rows:
        pairs.setdefault((row["sample"] + "/" + row["assay"], row["mode"], str(row["repeat"])), {})[row["backend"]] = row
    for (dataset, mode, repeat), pair in sorted(pairs.items()):
        cpu, cuda = pair.get("cpu", {}), pair.get("cuda", {})
        try:
            speedup = float(cpu["wall_time_s"]) / float(cuda["wall_time_s"])
            speedup_text = f"{speedup:.3f}x"
        except (KeyError, TypeError, ValueError, ZeroDivisionError):
            speedup_text = "—"
        material = (cpu or cuda).get("material", "")
        lines.append(f"| {dataset} | {material} | {mode} | {repeat} | {format_number(cpu.get('wall_time_s'))} | {format_number(cuda.get('wall_time_s'))} | {speedup_text} | {cuda.get('status', '—')} | {cuda.get('comparison_status', cpu.get('comparison_status', '—'))} |")

    lines += ["", "## Median wall time by dataset", "", "| Dataset | Mode | CPU median, s | CUDA median, s | CPU/CUDA |", "|---|---|---:|---:|---:|"]
    grouped: dict[tuple[str, str], dict[str, list[float]]] = {}
    for row in rows:
        try:
            value = float(row["wall_time_s"])
        except (TypeError, ValueError):
            continue
        grouped.setdefault((row["sample"] + "/" + row["assay"], row["mode"]), {}).setdefault(row["backend"], []).append(value)
    for key, values in sorted(grouped.items()):
        cpu = median(values.get("cpu", [])) if values.get("cpu") else None
        cuda = median(values.get("cuda", [])) if values.get("cuda") else None
        ratio = f"{cpu / cuda:.3f}x" if cpu is not None and cuda else "—"
        lines.append(f"| {key[0]} | {key[1]} | {format_number(cpu)} | {format_number(cuda)} | {ratio} |")

    counts: dict[str, int] = {}
    for row in rows:
        counts[row.get("status", "unknown")] = counts.get(row.get("status", "unknown"), 0) + 1
    lines += ["", "## Status", "", ", ".join(f"{status}: {count}" for status, count in sorted(counts.items())) or "No runs yet.", "", "## CUDA phase timings", "", "| Phase | Samples | Median, s |", "|---|---:|---:|"]
    for field, title in (("cuda_packing_s", "Packing"), ("h2d_s", "H2D"), ("kernel_s", "Kernel"), ("d2h_s", "D2H"), ("gpu_qc_total_s", "GPU QC total"), ("cpu_remaining_qc_s", "CPU remaining QC"), ("peak_rss_gib", "Peak RSS, GiB")):
        values = []
        for row in rows:
            if row.get("backend") != "cuda":
                continue
            try:
                values.append(float(row[field]))
            except (TypeError, ValueError):
                pass
        lines.append(f"| {title} | {len(values)} | {format_number(median(values) if values else None)} |")
    lines += ["", "Peak GPU memory is intentionally blank: sampling `nvidia-smi` before/after a short process is not a reliable peak measurement.", ""]
    (output / "summary.md").write_text("\n".join(lines), encoding="utf-8")


def capture_system_info(output: Path, threads: int | None) -> None:
    lines = [
        f"timestamp_utc: {dt.datetime.now(dt.timezone.utc).isoformat()}", f"platform: {platform.platform()}",
        f"python: {sys.version}", f"hostname: {platform.node()}", f"cpu_count: {os.cpu_count()}",
        f"OMP_NUM_THREADS: {threads if threads is not None else os.environ.get('OMP_NUM_THREADS', 'unset')}", "",
    ]
    try:
        result = subprocess.run(["nvidia-smi", "--query-gpu=name,memory.total,driver_version", "--format=csv,noheader"], text=True, capture_output=True, timeout=2, check=False)
        lines += ["nvidia-smi:", result.stdout.strip() or result.stderr.strip() or f"exit code {result.returncode}"]
    except (OSError, subprocess.TimeoutExpired) as exc:
        lines += ["nvidia-smi:", f"unavailable: {exc}"]
    (output / "system_info.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")


def planned_jobs(datasets: list[Dataset], mode: str, repeats: int, trim_set: set[str]) -> list[tuple[Dataset, str, int, str]]:
    jobs = []
    if mode in {"qc", "all"}:
        jobs.extend((dataset, "qc", repeat, backend) for dataset in datasets for repeat in range(1, repeats + 1) for backend in ("cpu", "cuda"))
    if mode in {"trim", "all"}:
        chosen = [dataset for dataset in datasets if dataset.dataset_id in trim_set]
        if len(chosen) != len(trim_set):
            missing = sorted(trim_set - {dataset.dataset_id for dataset in chosen})
            die(f"Configured full-trim dataset was not discovered: {', '.join(missing)}")
        jobs.extend((dataset, "trim", repeat, backend) for dataset in chosen for repeat in range(1, repeats + 1) for backend in ("cpu", "cuda"))
    return jobs


def on_signal(signum: int, _frame: Any) -> None:
    global CANCELLED
    CANCELLED = True
    print(f"Received signal {signum}; current NeoQC process will be stopped and results retained.", file=sys.stderr)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", type=Path, required=True, help="NeoQC executable")
    parser.add_argument("--data", type=Path, required=True, help="Root of BAI/.../VOR data directories")
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--mode", choices=("qc", "trim", "all"), default="all")
    parser.add_argument("--repeats", type=int, default=1)
    parser.add_argument("--resume", action="store_true", help="Skip runs already recorded in results.csv")
    parser.add_argument("--dry-run", action="store_true", help="Discover data and save the planned commands without launching NeoQC")
    parser.add_argument("--output", type=Path, help="Result directory (default: benchmark-results/<timestamp>)")
    parser.add_argument("--threads", type=int, help="Set identical OMP_NUM_THREADS for CPU and CUDA runs")
    parser.add_argument("--clean-output", action="store_true", help="Explicitly remove an existing --output before starting; input data is always protected")
    args = parser.parse_args(argv)
    if args.repeats < 1:
        parser.error("--repeats must be positive")
    if args.threads is not None and args.threads < 1:
        parser.error("--threads must be positive")
    if args.resume and args.clean_output:
        parser.error("--resume cannot be combined with --clean-output")
    if args.clean_output and args.output is None:
        parser.error("--clean-output requires explicit --output")
    return args


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    data, binary = args.data.resolve(), args.binary.resolve()
    config = load_config(args.config.resolve())
    expected = config.get("expected_samples", [])
    if not isinstance(expected, list) or not expected:
        die("Config expected_samples must be a non-empty list")
    datasets = discover_datasets(data, expected)
    timestamp = dt.datetime.now().strftime("%Y%m%d-%H%M%S")
    output = (args.output or Path("benchmark-results") / timestamp).resolve()
    ensure_safe_output(output, data)
    if args.clean_output and output.exists():
        shutil.rmtree(output)
    if output.exists() and not args.resume and any(output.iterdir()):
        die(f"Output already exists and is not empty: {output}; use --resume or an empty/new --output")
    output.mkdir(parents=True, exist_ok=True)
    if not args.dry_run and not binary.is_file():
        die(f"NeoQC binary does not exist: {binary}")
    serializable_args = {key: str(value) if isinstance(value, Path) else value for key, value in vars(args).items()}
    serializable_datasets = [
        {**asdict(item), "r1": str(item.r1), "r2": str(item.r2)}
        for item in datasets
    ]
    atomic_json(output / "config.json", {
        "benchmark_config": config,
        "arguments": serializable_args,
        "datasets": serializable_datasets,
    })
    capture_system_info(output, args.threads)
    signal.signal(signal.SIGINT, on_signal)
    signal.signal(signal.SIGTERM, on_signal)

    trim_set = set(config.get("full_trim_datasets", []))
    jobs = planned_jobs(datasets, args.mode, args.repeats, trim_set)
    records = load_existing(output) if args.resume else {}
    for dataset, mode, repeat, backend in jobs:
        if CANCELLED:
            break
        key = make_run_key(dataset, mode, repeat, backend)
        if args.resume and key in records:
            continue
        record = run_one(binary, dataset, mode, repeat, backend, output, config, args.threads, args.dry_run)
        records[key] = record
        write_results(output, records)
        # CPU and CUDA are adjacent by construction, so compare immediately after CUDA.
        if backend == "cuda":
            compare_pair(records, dataset, mode, repeat, output, config)
            write_results(output, records)
        write_summary(output, records)
    write_results(output, records)
    write_summary(output, records)
    return 130 if CANCELLED else 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except RuntimeError as exc:
        print(f"benchmark error: {exc}", file=sys.stderr)
        raise SystemExit(2)
