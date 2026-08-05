# NeoQC plot contract

NeoQC renders quality-control metrics through a single report-oriented plotting
pipeline. C++ writes stable TSV files; `scripts/plot_results.py` validates and
renders every available metric in one process.

## Output

With `--plot`, files are written to `<result>/plots/` in two formats:

- SVG is the primary responsive asset for the self-contained HTML report;
- PNG is a 2400 x 1350 pixel, 300 dpi fallback for export and print workflows.

`plots_manifest.json` is the integration boundary for downstream consumers.
It contains one entry for every supported metric and read direction, including
generated, skipped and failed plots. File paths in the manifest are relative to
the plot directory.

NeoQC also writes `<result>/neoqc_qc_report.html`. This is a compact,
self-contained sequencing QC report, separate from the complete clinical HTML
report. It groups R1 and R2 by metric, embeds chart assets as data URIs, includes
basic statistics, provides print/PDF styling and can open charts at full size.

Supported plot identifiers are:

- `per_base_quality`;
- `adapter_content`;
- `per_base_sequence_content`;
- `per_sequence_gc_content`;
- `per_base_n_content`;
- `sequence_length_distribution`;
- `per_sequence_quality`.

Missing R2 inputs are recorded as `source_not_found`. When adapter analysis is
disabled, adapter entries are recorded as `adapter_analysis_disabled`; all
other available plots are still generated.

## Standalone use

```bash
python3 scripts/plot_results.py results/sample01 results/sample01/plots
```

Generate only one format:

```bash
python3 scripts/plot_results.py results/sample01 results/sample01/plots \
  --formats svg --strict
```

The plotting process reports failures without invalidating TSV and summary
outputs. C++ treats a plotting failure as a warning.

Regenerate only the HTML report from existing artifacts:

```bash
python3 scripts/generate_qc_report.py results/sample01
```

Choose a custom output path:

```bash
python3 scripts/generate_qc_report.py results/sample01 \
  --output results/sample01/sample01_fastqc.html
```

The report badges `READY`, `NOT RUN` and `ERROR` describe artifact availability.
They intentionally do not claim FastQC-compatible biological
`PASS`/`WARNING`/`FAIL` decisions; that requires the separate metric threshold
engine.

## Visual language

All figures share `scripts/plot_style.py`, which matches the HTML report brand
palette, typography, spacing, grid, number formatting and accessible series
styles. Individual metric renderers must not define their own global theme.

The HTML renderer should consume SVG first, use PNG as a fallback, and embed
the chosen asset into the self-contained report. It should use `alt_text` and
`title` from the manifest instead of reconstructing chart metadata from file
names.
