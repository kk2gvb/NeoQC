# Sequence duplication calculation

NeoQC calculates sequence duplication during the same FASTQ pass as the other
quality metrics. Every distinct 50-nucleotide prefix is retained, so the result
is exact for this key definition and has no unique-sequence sampling limit.

## Method

1. The first 50 bases of each read form the duplication key.
2. A/C/G/T/N symbols are encoded with three bits each in a fixed 150-bit key.
3. Every unique key is counted through the end of the file.
4. Exact counts are grouped into 16
   stable bins: `1` through `9`, `>10`, `>50`, `>100`, `>500`, `>1k`, `>5k`
   and `>10k+`.

NeoQC retains the historical two-series TSV needed by its report: percentage of
total sequences and percentage after deduplication. Unlike FastQC's bounded
sampling strategy, NeoQC does not stop admitting new keys after a threshold.
Consequently memory grows with the number of distinct prefixes. The compact key
avoids per-key sequence strings, but very large high-complexity paired libraries
still require capacity planning.

## Output contract

`sequence_duplication_levels_R*.tsv` always has exactly these columns:

```text
duplication_level	total_sequences_percent	deduplicated_sequences_percent
```

Native runs additionally produce:

- `sequence_duplication_summary_R*.tsv` — source filename, algorithm ID,
  total/unique counts and the headline deduplicated percentage;
- `overrepresented_sequences_R*.tsv` — tracked sequences above 0.1% of all
  reads, sorted deterministically;
- `sequence_duplication_R*.incomplete` — a temporary transaction marker which
remains only when calculation or publication did not finish.

When the compact HTML report is generated, each duplication chart is followed
by its read-specific `Overrepresented sequences` table. The report shows the
sequence, exact count, percentage and possible source; a valid header-only TSV
is rendered explicitly as “No sequences exceeded the reporting threshold.”

The summary is published last. The QC evaluator rejects native artifacts while
the incomplete marker exists and validates the summary against the percentage
table. A PNG/SVG rendering error remains an artifact error and does not alter
the biological QC status.

## Compatibility

Older result directories may contain a two-series TSV imported from a FastQC
archive without a NeoQC summary, or a summary created by the earlier bounded
NeoQC prototype. Both remain readable for report compatibility, but new native
runs use `neoqc-exact-prefix-v1` and contain no sampling-limit fields.

Reference implementation: [FastQC `DuplicationLevel.java`](https://github.com/s-andrews/FastQC/blob/master/uk/ac/babraham/FastQC/Modules/DuplicationLevel.java)
and [`OverRepresentedSeqs.java`](https://github.com/s-andrews/FastQC/blob/master/uk/ac/babraham/FastQC/Modules/OverRepresentedSeqs.java).
