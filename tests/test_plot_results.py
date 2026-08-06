#!/usr/bin/env python3
"""Contract and rendering tests for the report-ready NeoQC charts."""

from __future__ import annotations

import json
import struct
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "plot_results.py"


FIXTURES = {
    "per_cycle": (
        "cycle\tmean_quality\tlower_quartile\tmedian\n"
        "1\t35.2\t33\t35\n2\t34.8\t32\t35\n3\t31.4\t28\t31\n4\t28.6\t24\t28\n"
    ),
    "adapter_content": "pos\tTruSeq\tNextera\n1\t0\t0\n2\t0.5\t0\n3\t2.5\t0.25\n4\t7.0\t1.5\n",
    "per_base_sequence_content": (
        "position\tA\tC\tG\tT\tN\n"
        "1\t25\t25\t25\t25\t0\n2\t27\t24\t24\t24\t1\n3\t29\t23\t23\t23\t2\n"
    ),
    "per_sequence_gc_content": "gc_percent\treads\n30\t5\n40\t30\n50\t80\n60\t25\n70\t4\n",
    "per_base_n_content": "position\tN_percent\n1\t0\n2\t0.2\n3\t1.5\n4\t0.4\n",
    "sequence_length_distribution": "length\treads\n75\t5\n100\t15\n150\t80\n",
    "sequence_duplication_levels": (
        "duplication_level\ttotal_sequences_percent\tdeduplicated_sequences_percent\n"
        "1\t13.4\t51.5\n2\t10.2\t20.4\n3\t8.1\t10.7\n"
        "4\t6.0\t6.1\n5\t3.8\t3.0\n>10\t14.0\t2.9\n>50\t5.2\t0.3\n"
    ),
    "per_sequence_quality": "mean_quality\tread_count\n20\t4\n30\t45\n35\t70\n40\t10\n",
}


def write_fixture_set(directory: Path, reads: tuple[str, ...] = ("R1", "R2")) -> None:
    directory.mkdir(parents=True, exist_ok=True)
    for read in reads:
        for prefix, content in FIXTURES.items():
            (directory / f"{prefix}_{read}.tsv").write_text(content, encoding="utf-8")


def run_plotter(input_dir: Path, output_dir: Path, *arguments: str) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [sys.executable, str(SCRIPT), str(input_dir), str(output_dir), *arguments],
        cwd=ROOT,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )


def png_metadata(path: Path) -> tuple[int, int, int | None]:
    data = path.read_bytes()
    if not data.startswith(b"\x89PNG\r\n\x1a\n"):
        raise AssertionError(f"{path} is not a PNG")
    offset = 8
    width = height = 0
    pixels_per_metre = None
    while offset < len(data):
        length = struct.unpack(">I", data[offset : offset + 4])[0]
        chunk_type = data[offset + 4 : offset + 8]
        chunk = data[offset + 8 : offset + 8 + length]
        if chunk_type == b"IHDR":
            width, height = struct.unpack(">II", chunk[:8])
        elif chunk_type == b"pHYs" and chunk[8] == 1:
            pixels_per_metre = struct.unpack(">I", chunk[:4])[0]
        offset += 12 + length
        if chunk_type == b"IEND":
            break
    return width, height, pixels_per_metre


class PlotResultsTest(unittest.TestCase):
    def test_complete_r1_r2_svg_contract(self) -> None:
        with tempfile.TemporaryDirectory(prefix="neo qc plots ") as temporary:
            root = Path(temporary)
            input_dir = root / "input data"
            output_dir = root / "report plots"
            write_fixture_set(input_dir)
            output_dir.mkdir(parents=True)
            (output_dir / "quality_distribution_R1.svg").write_text("stale", encoding="utf-8")
            (output_dir / "quality_distribution_R1.png").write_text("stale", encoding="utf-8")

            result = run_plotter(input_dir, output_dir, "--formats", "svg", "--strict")
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertNotIn("findfont:", result.stderr)

            manifest = json.loads((output_dir / "plots_manifest.json").read_text(encoding="utf-8"))
            self.assertEqual(manifest["schema_version"], 1)
            self.assertEqual(manifest["theme"], "neo-report")
            self.assertEqual(manifest["figure"]["aspect_ratio"], "16:9")
            self.assertEqual(manifest["figure"]["png_width_px"], 2400)
            self.assertEqual(manifest["figure"]["png_height_px"], 1350)
            self.assertEqual(manifest["figure"]["png_dpi"], 300)
            self.assertEqual(manifest["summary"], {"generated": 16, "errors": 0, "skipped": 0})
            self.assertEqual(len(manifest["plots"]), 16)
            for entry in manifest["plots"]:
                self.assertEqual(entry["status"], "generated")
                self.assertFalse(Path(entry["svg"]).is_absolute())
                svg = output_dir / entry["svg"]
                self.assertTrue(svg.is_file(), svg)
                self.assertIn("viewBox", svg.read_text(encoding="utf-8")[:1000])
                self.assertTrue(entry["alt_text"])
            report = input_dir / "neoqc_qc_report.html"
            self.assertTrue(report.is_file())
            report_text = report.read_text(encoding="utf-8")
            self.assertIn("data:image/svg+xml;base64,", report_text)
            self.assertIn("Per base sequence quality", report_text)
            self.assertIn("R1", report_text)
            self.assertIn("R2", report_text)
            self.assertNotIn("Base quality distribution", report_text)
            self.assertFalse((output_dir / "quality_distribution_R1.svg").exists())
            self.assertFalse((output_dir / "quality_distribution_R1.png").exists())
            gc_svg = (output_dir / "per_sequence_gc_content_R1.svg").read_text(encoding="utf-8")
            self.assertIn("Observed", gc_svg)
            self.assertIn("Theoretical distribution", gc_svg)
            duplication_svg = (output_dir / "sequence_duplication_levels_R1.svg").read_text(
                encoding="utf-8"
            )
            self.assertIn("Total sequences", duplication_svg)
            self.assertIn("Deduplicated sequences", duplication_svg)
            self.assertIn("&gt;50", duplication_svg)
            self.assertIn("Sequence duplication levels", report_text)
            self.assertIn("PASS / WARNING / FAIL distribution", report_text)
            self.assertIn("fastqc-compatible-v1", report_text)
            self.assertTrue((input_dir / "qc_evaluation.json").is_file())

    def test_duplication_percentages_outside_range_are_rejected(self) -> None:
        with tempfile.TemporaryDirectory(prefix="neoqc-duplication-invalid-") as temporary:
            root = Path(temporary)
            input_dir = root / "input"
            output_dir = root / "output"
            input_dir.mkdir()
            (input_dir / "sequence_duplication_levels_R1.tsv").write_text(
                "duplication_level\ttotal_sequences_percent\tdeduplicated_sequences_percent\n"
                "1\t101\t50\n",
                encoding="utf-8",
            )

            result = run_plotter(input_dir, output_dir, "--formats", "svg", "--strict")
            self.assertEqual(result.returncode, 1)
            manifest = json.loads((output_dir / "plots_manifest.json").read_text(encoding="utf-8"))
            failed = next(
                entry
                for entry in manifest["plots"]
                if entry["id"] == "sequence_duplication_levels" and entry["read"] == "R1"
            )
            self.assertEqual(failed["status"], "error")
            self.assertIn("between 0 and 100", failed["reason"])

    def test_fixed_length_chart_uses_readable_single_value_layout(self) -> None:
        with tempfile.TemporaryDirectory(prefix="neoqc-fixed-length-") as temporary:
            root = Path(temporary)
            input_dir = root / "input"
            output_dir = root / "output"
            input_dir.mkdir()
            (input_dir / "sequence_length_distribution_R1.tsv").write_text(
                "length\treads\n100\t66188141\n", encoding="utf-8"
            )

            result = run_plotter(input_dir, output_dir, "--formats", "svg", "--strict")
            self.assertEqual(result.returncode, 0, result.stderr)
            svg = (output_dir / "sequence_length_distribution_R1.svg").read_text(encoding="utf-8")
            self.assertIn("Fixed length", svg)
            self.assertIn("66.2M reads", svg)
            self.assertNotIn("Mean:", svg)
            self.assertNotIn("Mode:", svg)

            manifest = json.loads((output_dir / "plots_manifest.json").read_text(encoding="utf-8"))
            entry = next(item for item in manifest["plots"] if item["status"] == "generated")
            self.assertIn("FastQC-style sequence length line", entry["alt_text"])
            self.assertIn("fixed length of 100 bp", entry["alt_text"])

    def test_png_is_2400_by_1350_at_300_dpi(self) -> None:
        with tempfile.TemporaryDirectory(prefix="neoqc-png-") as temporary:
            root = Path(temporary)
            input_dir = root / "input"
            output_dir = root / "output"
            input_dir.mkdir()
            (input_dir / "per_cycle_R1.tsv").write_text(FIXTURES["per_cycle"], encoding="utf-8")

            result = run_plotter(input_dir, output_dir, "--formats", "png", "--strict")
            self.assertEqual(result.returncode, 0, result.stderr)
            width, height, pixels_per_metre = png_metadata(output_dir / "per_base_quality_R1.png")
            self.assertEqual((width, height), (2400, 1350))
            self.assertIsNotNone(pixels_per_metre)
            self.assertAlmostEqual(pixels_per_metre or 0, 11_811, delta=2)

    def test_skip_adapters_removes_stale_assets_and_records_reason(self) -> None:
        with tempfile.TemporaryDirectory(prefix="neoqc-skip-") as temporary:
            root = Path(temporary)
            input_dir = root / "input"
            output_dir = root / "output"
            input_dir.mkdir()
            (input_dir / "per_cycle_R1.tsv").write_text(FIXTURES["per_cycle"], encoding="utf-8")
            (input_dir / "adapter_content_R1.tsv").write_text(FIXTURES["adapter_content"], encoding="utf-8")
            output_dir.mkdir()
            (output_dir / "adapter_content_R1.svg").write_text("stale", encoding="utf-8")

            result = run_plotter(input_dir, output_dir, "--formats", "svg", "--skip-adapters", "--strict")
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertFalse((output_dir / "adapter_content_R1.svg").exists())
            manifest = json.loads((output_dir / "plots_manifest.json").read_text(encoding="utf-8"))
            adapter = next(
                entry for entry in manifest["plots"]
                if entry["id"] == "adapter_content" and entry["read"] == "R1"
            )
            self.assertEqual(adapter["status"], "skipped")
            self.assertEqual(adapter["reason"], "adapter_analysis_disabled")
            self.assertTrue((output_dir / "per_base_quality_R1.svg").is_file())

    def test_malformed_tsv_is_reported_in_manifest(self) -> None:
        with tempfile.TemporaryDirectory(prefix="neoqc-invalid-") as temporary:
            root = Path(temporary)
            input_dir = root / "input"
            output_dir = root / "output"
            input_dir.mkdir()
            (input_dir / "per_cycle_R1.tsv").write_text("cycle\twrong\n1\t30\n", encoding="utf-8")

            result = run_plotter(input_dir, output_dir, "--formats", "svg", "--strict")
            self.assertEqual(result.returncode, 1)
            manifest = json.loads((output_dir / "plots_manifest.json").read_text(encoding="utf-8"))
            self.assertEqual(manifest["summary"]["errors"], 1)
            failed = next(entry for entry in manifest["plots"] if entry["status"] == "error")
            self.assertIn("missing required column", failed["reason"])


if __name__ == "__main__":
    unittest.main()
