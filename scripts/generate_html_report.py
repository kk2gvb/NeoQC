#!/usr/bin/env python3
"""Repository-local entry point for complete neo-mRNA-vax reports."""

from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "python"))

from neo_mrna_vax_report.cli import main  # noqa: E402


raise SystemExit(main())
