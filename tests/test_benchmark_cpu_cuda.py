#!/usr/bin/env python3
"""Unit tests for the standard-library CPU/CUDA benchmark runner."""

import gzip
import importlib.util
import json
import sys
import tempfile
import unittest
from unittest import mock
from pathlib import Path


SCRIPT = Path(__file__).parents[1] / "scripts" / "benchmark_cpu_cuda.py"
SPEC = importlib.util.spec_from_file_location("benchmark_cpu_cuda", SCRIPT)
assert SPEC and SPEC.loader
benchmark = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = benchmark
SPEC.loader.exec_module(benchmark)


def write_gzip(path: Path, content: bytes = b"@r\nACGT\n+\nIIII\n") -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wb") as handle:
        handle.write(content)


class BenchmarkRunnerTests(unittest.TestCase):
    def make_data(self, root: Path) -> Path:
        data = root / "data"
        for sample in ("ONE",):
            write_gzip(data / sample / "dna" / "I_ONE_GED_R1.fastq.gz")
            write_gzip(data / sample / "dna" / "I_ONE_GED_R2.fastq.gz")
            write_gzip(data / sample / "dna" / "I_ONE_SOM_R1.fastq.gz")
            write_gzip(data / sample / "dna" / "I_ONE_SOM_R2.fastq.gz")
            write_gzip(data / sample / "rna" / "ONE_S1_L001_R1_001.fastq.gz")
            write_gzip(data / sample / "rna" / "ONE_S1_L001_R2_001.fastq.gz")
        return data

    def test_discovery_and_missing_mate(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            data = self.make_data(Path(temporary))
            datasets = benchmark.discover_datasets(data, ["ONE"])
            self.assertEqual([item.dataset_id for item in datasets], ["ONE/GED", "ONE/SOM", "ONE/RNA"])
            (data / "ONE" / "dna" / "I_ONE_GED_R2.fastq.gz").unlink()
            with self.assertRaisesRegex(RuntimeError, "Missing mate"):
                benchmark.discover_datasets(data, ["ONE"])

    def test_streamed_gzip_comparison_reports_first_difference(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            left, right = root / "left.fastq.gz", root / "right.fastq.gz"
            write_gzip(left, b"abcdef")
            write_gzip(right, b"abcxef")
            equal, detail = benchmark.stream_files_equal(left, right, 2)
            self.assertFalse(equal)
            self.assertEqual(detail["offset"], 3)
            write_gzip(right, b"abcdef")
            self.assertTrue(benchmark.stream_files_equal(left, right, 2)[0])

    def test_output_protection(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            data = root / "data"
            data.mkdir()
            with self.assertRaisesRegex(RuntimeError, "inside input"):
                benchmark.ensure_safe_output(data / "results", data)
            benchmark.ensure_safe_output(root / "benchmark-results", data)

    def test_dry_run_and_resume_record(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            data, output = self.make_data(root), root / "results"
            config = root / "config.json"
            config.write_text(json.dumps({
                "expected_samples": ["ONE"], "full_trim_datasets": ["ONE/GED"],
                "timeout_seconds": 1, "minimum_free_disk_gib": 0,
            }), encoding="utf-8")
            with mock.patch.object(benchmark, "capture_system_info"):
                result = benchmark.main([
                    "--binary", str(root / "not-needed-for-dry-run"), "--data", str(data),
                    "--config", str(config), "--mode", "qc", "--dry-run", "--output", str(output),
                ])
            self.assertEqual(result, 0)
            with (output / "results.csv").open(encoding="utf-8") as handle:
                records = list(__import__("csv").DictReader(handle))
            self.assertEqual(len(records), 6)
            self.assertTrue(all(item["status"] == "PLANNED" for item in records))
            self.assertTrue((output / "runs" / "ONE_GED" / "qc" / "repeat-001" / "cpu" / "command.json").is_file())
            with mock.patch.object(benchmark, "capture_system_info"):
                self.assertEqual(benchmark.main([
                    "--binary", str(root / "not-needed-for-dry-run"), "--data", str(data),
                    "--config", str(config), "--mode", "qc", "--dry-run", "--output", str(output), "--resume",
                ]), 0)

    def write_fake_neoqc(self, root: Path) -> Path:
        binary = root / "fake_neoqc.py"
        binary.write_text("""#!/usr/bin/env python3
import os
import pathlib
import sys

args = sys.argv[1:]
def value(name): return args[args.index(name) + 1]
out = pathlib.Path(value('--out'))
sample = value('--sample-id')
backend = value('--backend')
out.mkdir(parents=True, exist_ok=True)
suffix = 'different\\n' if backend == 'cuda' and os.environ.get('BENCH_FAKE_MISMATCH') else 'same\\n'
for read in ('R1', 'R2'):
    (out / f'{sample}_{read}_summary.txt').write_text('Processed reads : 2\\nTotal bases     : 8\\n')
    (out / f'per_cycle_{read}.tsv').write_text(suffix)
(out / 'run_manifest.json').write_text('nondeterministic\\n')
print('  Read FASTQ:             0.01 s')
print('  Parallel processing:    0.02 s')
print('  Merge analyzers:        0.01 s')
print('  Statistics:             0.01 s')
print('  Reports:                0.01 s')
print('  Peak RSS (VmHWM):       0.10 GiB')
""", encoding="utf-8")
        binary.chmod(0o755)
        return binary

    def test_runner_records_matching_cpu_cuda_outputs(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            data, binary = self.make_data(root), self.write_fake_neoqc(root)
            config = root / "config.json"
            config.write_text(json.dumps({
                "expected_samples": ["ONE"], "full_trim_datasets": ["ONE/GED"],
                "timeout_seconds": 10, "minimum_free_disk_gib": 0,
            }), encoding="utf-8")
            matching = root / "matching"
            dataset = benchmark.discover_datasets(data, ["ONE"])[0]
            settings = benchmark.load_config(config)
            cpu = benchmark.run_one(binary, dataset, "qc", 1, "cpu", matching, settings, 2, False)
            cuda = benchmark.run_one(binary, dataset, "qc", 1, "cuda", matching, settings, 2, False)
            records = {cpu["run_key"]: cpu, cuda["run_key"]: cuda}
            benchmark.compare_pair(records, dataset, "qc", 1, matching, settings)
            self.assertEqual({item["status"] for item in records.values()}, {"OK"})
            self.assertEqual({item["comparison_status"] for item in records.values()}, {"MATCH"})
            self.assertTrue(Path(cuda["stdout"]).is_file())
            self.assertTrue(Path(cuda["stderr"]).is_file())

            cpu_output = Path(records["ONE_GED:qc:repeat-001:cpu"]["output_dir"])
            cuda_output = Path(records["ONE_GED:qc:repeat-001:cuda"]["output_dir"])
            (cuda_output / "per_cycle_R1.tsv").write_text("changed\n", encoding="utf-8")
            mismatch = benchmark.compare_outputs(cpu_output, cuda_output, {})
            self.assertEqual(mismatch["status"], "MISMATCH")
            self.assertEqual(mismatch["artifact"], "per_cycle_R1.tsv")


if __name__ == "__main__":
    unittest.main()
