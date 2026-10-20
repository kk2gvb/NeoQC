# Deterministic trimming fixtures

All qualities use Phred+33. Expected behavior with the named operation enabled:

| File | Input | Operation | Expected result |
|---|---|---|---|
| `clean.fastq` | `ACGTACGT` / `IIIIIIII` | none | PASS, unchanged |
| `adapter.fastq` | `ACGTAGATCGGA` | adapter `AGATCGGA` | PASS, `ACGT` |
| `low_quality.fastq` | `ACGTACGT` / `IIIIII!!` | quality tail Q20 | PASS, `ACGTAC` |
| `polyG.fastq` | `ACGTGGGG` | polyG threshold 4 | PASS, `ACGT` |
| `polyX.fastq` | `ACGTAAAA` | polyX threshold 4 | PASS, `ACGT` |
| `short.fastq` | `ACG` | minimum length 4 | DISCARD |
| `PE_clean_*` | insert-spanning pair | no trimming | PASS, unchanged and synchronized |
| `PE_adapter_*` | `ACGTAC` plus 3' tails | PE read-through | PASS, both reads become the six-base insert orientations |

