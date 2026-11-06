# CPU/CUDA benchmark NeoQC

`scripts/benchmark_cpu_cuda.py` runs the same paired FASTQ input through NeoQC
sequentially with `--backend cpu` and `--backend cuda`.  It does not change QC,
CUDA, or trimming logic.  CUDA currently owns only the implemented per-base
metrics; the rest of QC and all trimming remain CPU/OpenMP work.

## Build and run

Build both variants in separate directories. The CPU executable is useful for
the CPU reference; the CUDA executable is required for the CUDA side. `build.sh`
supports `--cuda` and `--no-cuda` but intentionally uses its fixed `build/`
directory, so use CMake directly when both artifacts must coexist.

```bash
cmake -S . -B build-cpu -DNEOQC_ENABLE_CUDA=OFF -DCMAKE_BUILD_TYPE=Release
cmake --build build-cpu --parallel
cmake -S . -B build-cuda -DNEOQC_ENABLE_CUDA=ON -DCMAKE_CUDA_ARCHITECTURES=86 -DCMAKE_BUILD_TYPE=Release
cmake --build build-cuda --parallel

python3 scripts/benchmark_cpu_cuda.py \
  --binary ./build-cuda/neoqc \
  --data ./data \
  --config scripts/benchmark_config.json \
  --mode all --repeats 1 --threads 20
```

The CUDA-built binary supports both requested backends, so one executable is
used in the example.  Before an overnight run, check discovery and generated
commands only:

```bash
python3 scripts/benchmark_cpu_cuda.py \
  --binary ./build-cuda/neoqc --data ./data --mode all --dry-run
```

Use a known result directory to continue after interruption:

```bash
python3 scripts/benchmark_cpu_cuda.py \
  --binary ./build-cuda/neoqc --data ./data \
  --output benchmark-results/20260101-120000 --resume --threads 20
```

## Dataset contract

The default configuration expects 27 paired datasets: GED, SOM, and RNA for
each of `BAI`, `CHER`, `CHU`, `IVA`, `KLU`, `KUD`, `LAP`, `SUN`, and `VOR`.
For each sample it requires exactly one R1 file for each pattern and derives a
single R2 mate:

```text
data/<sample>/dna/*_GED_R1.fastq.gz  -> *_GED_R2.fastq.gz
data/<sample>/dna/*_SOM_R1.fastq.gz  -> *_SOM_R2.fastq.gz
data/<sample>/rna/*_R1_001.fastq.gz  -> *_R2_001.fastq.gz
```

Ambiguous patterns, missing sample directories, missing mates, and reused paths
are errors.  This deliberately prevents mixing GED, SOM, and RNA inputs.

`--mode qc` covers all 27 inputs. `--mode trim` runs only CHU/GED, BAI/RNA, and
LAP/SOM as declared in `benchmark_config.json`; `--mode all` runs both.

## Outputs and safety

The default output is `benchmark-results/<timestamp>/` and contains:

- `results.csv`: one durable record per process invocation;
- `summary.md`: CPU/CUDA table, medians, status counts, and CUDA-phase medians;
- `config.json` and `system_info.txt`;
- `runs/<dataset>/<mode>/repeat-XXX/{cpu,cuda}/`: command, separate stdout and
  stderr, NeoQC output, and the per-run JSON record. `comparison.json` is next
  to the CPU/CUDA pair.

The runner writes CSV and summary after every invocation.  It runs CPU then CUDA
for a dataset, never concurrently. `--threads N` supplies the same
`OMP_NUM_THREADS=N` and `OMP_DYNAMIC=FALSE` to both processes. Each process has
the configured timeout (default two hours). Full trim checks the configured
free-space reserve before launch (default 20 GiB).

Input `data/` is protected: an output directory inside it is rejected. Results
are retained by default. The only removal option is `--clean-output`, which
requires an explicit `--output`, cannot be combined with `--resume`, and still
rejects `data/` and its descendants.

## Correctness checks and interpretation

All output files except `run_manifest.json` are compared exactly. For `.gz`
files, comparison is streamed after decompression, so trimmed FASTQ content is
compared without loading the file into memory. The first missing or differing
artifact and byte offset are saved in `comparison.json`; runs are marked
`MISMATCH` and the remaining plan continues.

The runner separately detects current NeoQC messages `CUDA QC batch failed` and
`Falling back to CPU`. A CUDA run that emits either is `FALLBACK`, never a
successful CUDA result. This is reliable only for fallback paths that currently
log those messages; lack of a message is not proof that every batch ran on GPU.
A future machine-readable per-batch backend counter would make this audit
stronger.

Wall time is measured externally with `time.perf_counter()`. NeoQC internal
timers are parsed separately and are not summed. Peak RSS is taken from NeoQC's
`--timing` output. Peak GPU memory is left empty because sampling `nvidia-smi`
before/after a short process is not a reliable peak measurement.

The runner intentionally does not claim a speedup: inspect `OK` CUDA rows with
`MATCH` comparisons and no fallback before interpreting `CPU/CUDA` ratios.
