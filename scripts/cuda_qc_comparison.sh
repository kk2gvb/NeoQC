#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 4 || $# -gt 5 ]]; then
    echo "Usage: $0 <neoqc-binary> <r1.fastq[.gz]> <sample-id> <output-root> [r2.fastq[.gz]]" >&2
    exit 2
fi

binary=$1
r1=$2
sample_id=$3
output_root=$4
r2=${5:-}

cpu_out="$output_root/cpu"
cuda_out="$output_root/cuda"
common=(--r1 "$r1" --sample-id "$sample_id" --timing)
if [[ -n "$r2" ]]; then
    common+=(--r2 "$r2")
fi

mkdir -p "$output_root"

echo "=== CPU reference ==="
"$binary" "${common[@]}" --backend cpu --out "$cpu_out" | tee "$output_root/cpu.log"

echo "=== CUDA per-base QC ==="
"$binary" "${common[@]}" --backend cuda --out "$cuda_out" | tee "$output_root/cuda.log"

echo "=== Exact TSV comparison for CUDA-owned metrics ==="
for read in R1 R2; do
    for report in "per_cycle_${read}.tsv" \
                  "per_base_sequence_content_${read}.tsv" \
                  "per_base_n_content_${read}.tsv"; do
        cpu_report="$cpu_out/$report"
        cuda_report="$cuda_out/$report"
        if [[ -e "$cpu_report" || -e "$cuda_report" ]]; then
            diff -u "$cpu_report" "$cuda_report"
        fi
    done
done

echo "CPU/CUDA per-base TSV outputs match exactly."
