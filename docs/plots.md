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

NeoQC also writes `<result>/neoqc_qc_report.html`. This is a self-contained
sequencing QC report. It groups R1 and R2 by metric, embeds chart assets as data
URIs, includes basic statistics, provides print/PDF styling and can open charts
at full size.

NeoQC also writes `<result>/qc_evaluation.json`. It contains explainable,
versioned technical QC decisions for every metric/read pair. Plot artifact
status remains independent from QC status.

Supported plot identifiers are:

- `per_base_quality`;
- `adapter_content`;
- `per_base_sequence_content`;
- `per_sequence_gc_content`;
- `per_base_n_content`;
- `sequence_length_distribution`;
- `sequence_duplication_levels`;
- `per_sequence_quality`.

`sequence_duplication_levels_R1.tsv` and
`sequence_duplication_levels_R2.tsv` use the columns
`duplication_level`, `total_sequences_percent` and
`deduplicated_sequences_percent`. Duplication levels are labels (for example,
`1`, `2`, `>10`), while both percentage columns must contain values from 0 to
100.

Native FASTQ analysis also writes
`sequence_duplication_summary_R1.tsv`/`R2.tsv`. This one-row provenance record
contains the algorithm identifier, source FASTQ filename, prefix length,
total/unique counts and the exact `deduplicated_remaining_percent` used by the QC
decision engine. The level-1 ratio is retained only as a compatibility fallback
for previously imported FastQC two-series TSV files.

Duplication outputs are a small transaction. NeoQC removes stale artifacts at
the beginning of a run, publishes each file atomically, publishes the summary
last and removes `sequence_duplication_R1.incomplete`/`R2.incomplete` only after
success. The evaluator rejects a set while this marker exists, so a crashed run
cannot silently reuse a partial or stale duplication profile. The full method
and file semantics are documented in
[`sequence-duplication.md`](sequence-duplication.md).

`per_cycle_R1.tsv` and `per_cycle_R2.tsv` use `cycle`, `mean_quality`,
`lower_quartile` and `median`. The latter two columns are required for the
FastQC-compatible per-base quality decision. Older mean-only files still plot,
but are reported as `NOT EVALUATED` rather than receiving an inferred PASS.

`per_sequence_quality_R1.tsv` and `per_sequence_quality_R2.tsv` use
`mean_quality`, `read_count` and `read_count_truncate`. The rounded distribution
is the native NeoQC view; the truncated distribution matches FastQC binning and
is used by the FastQC-compatible decision profile. Legacy two-column files
without `read_count_truncate` remain readable and fall back to `read_count`.

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

Evaluate existing TSV files without rendering charts:

```bash
python3 scripts/evaluate_qc.py results/sample01
```

Select an explicit ruleset when plotting:

```bash
python3 scripts/plot_results.py results/sample01 results/sample01/plots \
  --ruleset config/qc_rules/fastqc-compatible-v1.json
```

The report displays technical `PASS`, `WARNING`, `FAIL` and `NOT EVALUATED`
decisions from the named ruleset. Artifact `generated`, `skipped` and `error`
states remain separate; a rendering failure is never converted into QC FAIL.
The data contract, ruleset strategy and aggregation are described in
[`qc-status-engine.md`](qc-status-engine.md).

## Visual language

All figures share `scripts/plot_style.py`, which matches the HTML report brand
palette, typography, spacing, grid, number formatting and accessible series
styles. Individual metric renderers must not define their own global theme.

The HTML renderer should consume SVG first, use PNG as a fallback, and embed
the chosen asset into the self-contained report. It should use `alt_text` and
`title` from the manifest instead of reconstructing chart metadata from file
names.
