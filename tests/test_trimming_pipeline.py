#!/usr/bin/env python3

import gzip
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import sys
import tempfile


NEOQC = Path(sys.argv[1]).resolve()
ROOT = Path(__file__).resolve().parents[1]


def write_fastq(path: Path, records: list[tuple[str, str, str]]) -> None:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "wt", encoding="ascii", newline="\n") as output:
        for header, sequence, quality in records:
            output.write(f"{header}\n{sequence}\n+\n{quality}\n")


def read_fastq(path: Path) -> list[tuple[str, str, str, str]]:
    with gzip.open(path, "rt", encoding="ascii") as source:
        lines = source.read().splitlines()
    if len(lines) % 4 != 0:
        raise AssertionError(f"malformed FASTQ output: {path}")
    records = []
    for index in range(0, len(lines), 4):
        header, sequence, separator, quality = lines[index:index + 4]
        assert header.startswith("@")
        assert separator.startswith("+")
        assert len(sequence) == len(quality)
        records.append((header, sequence, separator, quality))
    return records


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run(arguments: list[str], expect_success: bool = True) -> subprocess.CompletedProcess[str]:
    environment = os.environ.copy()
    environment["OMP_NUM_THREADS"] = "2"
    result = subprocess.run(
        [str(NEOQC), *arguments],
        cwd=ROOT,
        env=environment,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )
    if expect_success and result.returncode != 0:
        raise AssertionError(
            f"NeoQC failed ({result.returncode}):\n{result.stdout}\n{result.stderr}")
    if not expect_success and result.returncode == 0:
        raise AssertionError(f"NeoQC unexpectedly succeeded: {' '.join(arguments)}")
    return result


def summary_read_count(path: Path) -> int:
    match = re.search(r"Processed reads\s*:\s*(\d+)", path.read_text())
    if match is None:
        raise AssertionError(f"read count missing from {path}")
    return int(match.group(1))


def test_cli(root: Path) -> None:
    input_path = root / "cli.fastq.gz"
    write_fastq(input_path, [("@CLI", "ACGTAC", "IIIIII")])
    base = ["--r1", str(input_path), "--sample-id", "cli"]
    adapter_fasta = root / "adapters.fa"
    adapter_fasta.write_text(">one\nAGAT\nCGGA\n>two\nTTTTTT\n", encoding="ascii")
    positive = [
        [],
        ["--trim-front", "1"],
        ["--trim-tail", "1"],
        ["--cut-front", "--quality-threshold", "20"],
        ["--cut-tail", "--quality-threshold", "20"],
        ["--cut-right", "--quality-threshold", "20", "--window-size", "2"],
        ["--adapter-sequence", "TTTTTT"],
        ["--adapter-fasta", str(adapter_fasta)],
        ["--umi", "--umi-length", "1"],
        ["--min-length", "1"],
        ["--trim-polyg", "--polyg-min-length", "4"],
        ["--trim-polyx", "--polyx-min-length", "4"],
    ]
    for index, options in enumerate(positive):
        output = root / f"cli-positive-{index}"
        run([*base, "--out", str(output), "--trim", *options])
        assert (output / "trimmed" / "R1.trimmed.fastq.gz").is_file()

    invalid = [
        (["--trim-front", "1"], "require --trim"),
        (["--trim", "--trim-front", "invalid"], "non-negative integer"),
        (["--trim", "--cut-right", "--window-size", "0"], "greater than zero"),
        (["--trim", "--trim-front"], "requires a value"),
        (["--trim", "--adapter-sequence", "ACGTZQ"], "only A, C, G, T, or N"),
        (["--trim", "--unknown-trim-option"], "Unknown argument"),
        (["--trim", "--adapter-fasta", str(root / "missing.fa")], "Cannot open adapter FASTA"),
        (["--trim", "--overlap-correction"], "require paired-end input"),
        (["--trim", "--merge"], "require paired-end input"),
        (["--trim", "--umi"], "greater than zero"),
        (["--trim", "--umi-length", "2"], "requires --umi"),
        (["--trim", "--trim-polyg", "--polyg-min-length", "0"], "greater than zero"),
        (["--trim", "--trim-polyx", "--polyx-min-length", "0"], "greater than zero"),
    ]
    for index, (options, message) in enumerate(invalid):
        output = root / f"cli-invalid-{index}"
        result = run([*base, "--out", str(output), *options], expect_success=False)
        assert message in result.stderr, (message, result.stderr)

    unsafe = run(
        [*base, "--out", str(root), "--trim"],
        expect_success=False,
    )
    assert "must not contain an input FASTQ" in unsafe.stderr

    no_trim_output = root / "no-trim-regression"
    run([*base, "--out", str(no_trim_output)])
    assert not (no_trim_output / "trimmed").exists()
    assert not (no_trim_output / "trimming_report.json").exists()


def common_options() -> list[str]:
    return [
        "--trim",
        "--cut-tail",
        "--quality-threshold", "20",
        "--adapter-sequence", "TTTTTT",
        "--trim-polyg", "--polyg-min-length", "4",
        "--trim-polyx", "--polyx-min-length", "4",
        "--min-length", "4",
    ]


def test_single_end(root: Path) -> None:
    input_path = root / "single.fastq.gz"
    records = [
        ("@clean metadata", "ACGTACGT", "IIIIIIII"),
        ("@adapter", "ACGTACTTTTTT", "IIIIIIIIIIII"),
        ("@quality", "ACGTACGT", "IIIIII!!"),
        ("@polyg", "ACGTACGGGG", "IIIIIIIIII"),
        ("@polyx", "ACGTACAAAA", "IIIIIIIIII"),
        ("@short", "ACGGGG", "IIIIII"),
    ]
    write_fastq(input_path, records)
    original_digest = digest(input_path)
    output = root / "single-output"
    run([
        "--r1", str(input_path),
        "--sample-id", "single_e2e",
        "--out", str(output),
        *common_options(),
    ])

    assert digest(input_path) == original_digest
    trimmed = read_fastq(output / "trimmed" / "R1.trimmed.fastq.gz")
    assert [record[0].split()[0] for record in trimmed] == [
        "@clean", "@adapter", "@quality", "@polyg", "@polyx"
    ]
    assert [record[1] for record in trimmed] == [
        "ACGTACGT", "ACGTAC", "ACGTAC", "ACGTAC", "ACGTAC"
    ]

    report = json.loads((output / "trimming_report.json").read_text())
    assert report["reads"] == {"total": 6, "passed": 5, "discarded": 1}
    assert report["bases"] == {"before": 54, "after": 34, "trimmed": 20}
    assert report["trimming"]["quality_trimmed_reads"] == 1
    assert report["trimming"]["adapter_trimmed_reads"] == 1
    assert report["trimming"]["polyG_trimmed_reads"] == 2
    assert report["trimming"]["polyX_trimmed_reads"] == 1
    assert report["filtering"]["too_short_reads"] == 1
    assert report["output"]["reads"] == 5
    assert summary_read_count(
        output / "qc" / "before" / "single_e2e_R1_summary.txt") == 6
    assert summary_read_count(
        output / "qc" / "after" / "single_e2e_R1_summary.txt") == 5


def normalized(header: str) -> str:
    identifier = header.split()[0].removeprefix("@")
    return identifier[:-2] if identifier.endswith(("/1", "/2")) else identifier


def reverse_complement(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


def test_advanced_paired_end(root: Path) -> None:
    overlap = "ACGTTGCAACGATCGTACCTGATCGTTAACGT"
    altered = "T" + overlap[1:]
    r1_sequence = "GGGG" + overlap
    r2_reverse = altered + "CCCC"
    r2_sequence = reverse_complement(r2_reverse)
    r1_quality = "I" * 4 + "!" + "I" * (len(overlap) - 1)
    r2_quality = "I" * len(r2_sequence)
    r1_path = root / "advanced_R1.fastq.gz"
    r2_path = root / "advanced_R2.fastq.gz"
    write_fastq(r1_path, [("@ADVANCED/1", r1_sequence, r1_quality)])
    write_fastq(r2_path, [("@ADVANCED/2", r2_sequence, r2_quality)])
    original = digest(r1_path), digest(r2_path)
    output = root / "advanced-output"
    run([
        "--r1", str(r1_path), "--r2", str(r2_path),
        "--sample-id", "advanced", "--out", str(output),
        "--trim", "--overlap-correction", "--merge",
    ])

    assert (digest(r1_path), digest(r2_path)) == original
    trimmed_r1 = read_fastq(output / "trimmed" / "R1.trimmed.fastq.gz")
    trimmed_r2 = read_fastq(output / "trimmed" / "R2.trimmed.fastq.gz")
    merged = read_fastq(output / "trimmed" / "merged.fastq.gz")
    assert len(trimmed_r1) == len(trimmed_r2) == len(merged) == 1
    assert trimmed_r1[0][0] == "@ADVANCED/1"
    assert trimmed_r2[0][0] == "@ADVANCED/2"
    assert merged[0][0] == "@ADVANCED merged"
    assert len(merged[0][1]) == len(merged[0][3])
    report = json.loads((output / "trimming_report.json").read_text())
    assert report["trimming"]["corrected_reads"] == 1
    assert report["trimming"]["corrected_bases"] == 1
    assert report["trimming"]["merged_pairs"] == 1

    mismatch_r1 = root / "count_mismatch_R1.fastq.gz"
    mismatch_r2 = root / "count_mismatch_R2.fastq.gz"
    write_fastq(mismatch_r1, [
        ("@COUNT_1/1", "ACGT", "IIII"),
        ("@COUNT_2/1", "ACGT", "IIII"),
    ])
    write_fastq(mismatch_r2, [("@COUNT_1/2", "ACGT", "IIII")])
    failed_output = root / "count-mismatch-output"
    failure = run([
        "--r1", str(mismatch_r1), "--r2", str(mismatch_r2),
        "--sample-id", "count_mismatch", "--out", str(failed_output),
        "--trim",
    ], expect_success=False)
    assert "R2 contains fewer reads than R1" in failure.stderr
    assert not failed_output.exists()


def test_paired_end(root: Path) -> None:
    r1_path = root / "paired_R1.fastq.gz"
    r2_path = root / "paired_R2.fastq.gz"
    r1_records = [
        ("@PAIR_001/1", "ACGTAC", "IIIIII"),
        ("@PAIR_002/1", "ACGTACTTTTTT", "IIIIIIIIIIII"),
        ("@PAIR_003/1", "ACGTACGT", "IIIIII!!"),
        ("@PAIR_004/1", "ACGTACGGGG", "IIIIIIIIII"),
        ("@PAIR_005/1", "ACGGGG", "IIIIII"),
    ]
    r2_records = [
        ("@PAIR_001/2", "TGCATG", "IIIIII"),
        ("@PAIR_002/2", "TGCAACTTTTTT", "IIIIIIIIIIII"),
        ("@PAIR_003/2", "TGCATGCA", "IIIIII!!"),
        ("@PAIR_004/2", "TGCAACAAAA", "IIIIIIIIII"),
        ("@PAIR_005/2", "TGCATG", "IIIIII"),
    ]
    write_fastq(r1_path, r1_records)
    write_fastq(r2_path, r2_records)
    original_digests = (digest(r1_path), digest(r2_path))
    output = root / "paired-output"
    run([
        "--r1", str(r1_path), "--r2", str(r2_path),
        "--sample-id", "paired_e2e",
        "--out", str(output),
        *common_options(),
    ])

    assert (digest(r1_path), digest(r2_path)) == original_digests
    trimmed_r1 = read_fastq(output / "trimmed" / "R1.trimmed.fastq.gz")
    trimmed_r2 = read_fastq(output / "trimmed" / "R2.trimmed.fastq.gz")
    assert len(trimmed_r1) == len(trimmed_r2) == 4
    assert [normalized(record[0]) for record in trimmed_r1] == [
        "PAIR_001", "PAIR_002", "PAIR_003", "PAIR_004"
    ]
    assert [normalized(record[0]) for record in trimmed_r1] == [
        normalized(record[0]) for record in trimmed_r2
    ]

    report = json.loads((output / "trimming_report.json").read_text())
    assert report["reads"] == {"total": 10, "passed": 9, "discarded": 1}
    assert report["output"] == {
        "reads": 8,
        "mate_discarded_reads": 1,
        "total_pairs": 5,
        "passed_pairs": 4,
        "discarded_pairs": 1,
    }
    for read in ("R1", "R2"):
        assert summary_read_count(
            output / "qc" / "before" / f"paired_e2e_{read}_summary.txt") == 5
        assert summary_read_count(
            output / "qc" / "after" / f"paired_e2e_{read}_summary.txt") == 4


def main() -> int:
    with tempfile.TemporaryDirectory(prefix="neoqc-trimming-e2e-") as directory:
        root = Path(directory)
        test_cli(root)
        test_single_end(root)
        test_paired_end(root)
        test_advanced_paired_end(root)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
