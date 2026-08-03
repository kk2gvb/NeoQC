#!/usr/bin/env python3
"""Repository-local entry point for the NeoQC HTML report module."""

from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "python"))

from neoqc_report.cli import main  # noqa: E402


raise SystemExit(main())
