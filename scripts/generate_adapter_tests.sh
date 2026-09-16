#!/usr/bin/env bash

set -euo pipefail

OUT_DIR="tests/data"
mkdir -p "$OUT_DIR"

# Built-in NeoQC adapter set.
TRUSEQ_R1="AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
TRUSEQ_R2="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"
SMALLRNA3="TGGAATTCTCGGGTGCCAAGG"
SMALLRNA5="GATCGTCGGACTGTAGAACTCTGAAC"
NEXTERA="CTGTCTCTTATACACATCT"
POLYA="AAAAAAAAAAAA"
POLYG="GGGGGGGGGGGG"

# Exactly 30 nt, so adapter detection starts at position 31.
PREFIX="ACGTACGTACGTACGTACGTACGTACGTAC"

echo "Generating adapter test datasets..."

########################################
# Helpers
########################################

write_record() {
    local file="$1"
    local name="$2"
    local sequence="$3"

    local quality
    quality=$(printf 'I%.0s' $(seq 1 "${#sequence}"))

    {
        printf '@%s\n' "$name"
        printf '%s\n' "$sequence"
        printf '+\n'
        printf '%s\n' "$quality"
    } >> "$file"
}

reset_file() {
    : > "$1"
}

########################################
# adapter_profile
#
# Full TruSeq adapter on every read.
########################################

reset_file "${OUT_DIR}/adapter_profile_R1.fq"
reset_file "${OUT_DIR}/adapter_profile_R2.fq"

write_record \
    "${OUT_DIR}/adapter_profile_R1.fq" \
    "READ_001/1" \
    "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA${TRUSEQ_R1}"

write_record \
    "${OUT_DIR}/adapter_profile_R1.fq" \
    "READ_002/1" \
    "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCC${TRUSEQ_R1}"

write_record \
    "${OUT_DIR}/adapter_profile_R2.fq" \
    "READ_001/2" \
    "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT${TRUSEQ_R2}"

write_record \
    "${OUT_DIR}/adapter_profile_R2.fq" \
    "READ_002/2" \
    "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGG${TRUSEQ_R2}"

########################################
# partial_adapter
#
# One positive and one negative read.
########################################

reset_file "${OUT_DIR}/partial_adapter_R1.fq"
reset_file "${OUT_DIR}/partial_adapter_R2.fq"

write_record \
    "${OUT_DIR}/partial_adapter_R1.fq" \
    "READ_001/1" \
    "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA${TRUSEQ_R1}"

write_record \
    "${OUT_DIR}/partial_adapter_R1.fq" \
    "READ_002/1" \
    "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"

write_record \
    "${OUT_DIR}/partial_adapter_R2.fq" \
    "READ_001/2" \
    "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT${TRUSEQ_R2}"

write_record \
    "${OUT_DIR}/partial_adapter_R2.fq" \
    "READ_002/2" \
    "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"

########################################
# adapter_all
#
# One read per built-in adapter.
#
# R1/R2 are paired by identical read IDs.
# Every adapter starts at position 31.
########################################

reset_file "${OUT_DIR}/adapter_all_R1.fq"
reset_file "${OUT_DIR}/adapter_all_R2.fq"

write_record \
    "${OUT_DIR}/adapter_all_R1.fq" \
    "READ_001/1" \
    "${PREFIX}${TRUSEQ_R1}"

write_record \
    "${OUT_DIR}/adapter_all_R2.fq" \
    "READ_001/2" \
    "${PREFIX}${TRUSEQ_R2}"

write_record \
    "${OUT_DIR}/adapter_all_R1.fq" \
    "READ_002/1" \
    "${PREFIX}${SMALLRNA3}"

write_record \
    "${OUT_DIR}/adapter_all_R2.fq" \
    "READ_002/2" \
    "${PREFIX}${SMALLRNA3}"

write_record \
    "${OUT_DIR}/adapter_all_R1.fq" \
    "READ_003/1" \
    "${PREFIX}${SMALLRNA5}"

write_record \
    "${OUT_DIR}/adapter_all_R2.fq" \
    "READ_003/2" \
    "${PREFIX}${SMALLRNA5}"

write_record \
    "${OUT_DIR}/adapter_all_R1.fq" \
    "READ_004/1" \
    "${PREFIX}${NEXTERA}"

write_record \
    "${OUT_DIR}/adapter_all_R2.fq" \
    "READ_004/2" \
    "${PREFIX}${NEXTERA}"

write_record \
    "${OUT_DIR}/adapter_all_R1.fq" \
    "READ_005/1" \
    "${PREFIX}${POLYA}"

write_record \
    "${OUT_DIR}/adapter_all_R2.fq" \
    "READ_005/2" \
    "${PREFIX}${POLYA}"

write_record \
    "${OUT_DIR}/adapter_all_R1.fq" \
    "READ_006/1" \
    "${PREFIX}${POLYG}"

write_record \
    "${OUT_DIR}/adapter_all_R2.fq" \
    "READ_006/2" \
    "${PREFIX}${POLYG}"
########################################
# adapter_negative
#
# No configured adapter detection sequence.
########################################

reset_file "${OUT_DIR}/adapter_negative_R1.fq"
reset_file "${OUT_DIR}/adapter_negative_R2.fq"

write_record \
    "${OUT_DIR}/adapter_negative_R1.fq" \
    "NEGATIVE_001" \
    "${PREFIX}TGCATGCATGCATGCATGCATGCATGCATGCATGC"

write_record \
    "${OUT_DIR}/adapter_negative_R1.fq" \
    "NEGATIVE_002" \
    "${PREFIX}CACACACACACACACACACACACACACACACACACA"

write_record \
    "${OUT_DIR}/adapter_negative_R2.fq" \
    "NEGATIVE_001" \
    "${PREFIX}TGCATGCATGCATGCATGCATGCATGCATGCATGC"

write_record \
    "${OUT_DIR}/adapter_negative_R2.fq" \
    "NEGATIVE_002" \
    "${PREFIX}CACACACACACACACACACACACACACACACACACA"

########################################
# adapter_multiple
#
# Same adapter occurs twice in one read.
########################################

reset_file "${OUT_DIR}/adapter_multiple_R1.fq"
reset_file "${OUT_DIR}/adapter_multiple_R2.fq"

write_record \
    "${OUT_DIR}/adapter_multiple_R1.fq" \
    "MULTIPLE_001" \
    "${TRUSEQ_R1}AAAAAAAAAAAAAAAA${TRUSEQ_R1}"

write_record \
    "${OUT_DIR}/adapter_multiple_R2.fq" \
    "MULTIPLE_001" \
    "${TRUSEQ_R2}TTTTTTTTTTTTTTTT${TRUSEQ_R2}"

########################################
# adapter_short
#
# Reads shorter than the 12-nt detection
# sequence. No adapter should be detected.
########################################

reset_file "${OUT_DIR}/adapter_short_R1.fq"
reset_file "${OUT_DIR}/adapter_short_R2.fq"

write_record \
    "${OUT_DIR}/adapter_short_R1.fq" \
    "SHORT_001" \
    "ACGTACGTACG"

write_record \
    "${OUT_DIR}/adapter_short_R1.fq" \
    "SHORT_002" \
    "AAAAAAAAAAA"

write_record \
    "${OUT_DIR}/adapter_short_R2.fq" \
    "SHORT_001" \
    "ACGTACGTACG"

write_record \
    "${OUT_DIR}/adapter_short_R2.fq" \
    "SHORT_002" \
    "CCCCCCCCCCC"

########################################

echo "Generated:"
echo "  ${OUT_DIR}/adapter_profile_R1.fq"
echo "  ${OUT_DIR}/adapter_profile_R2.fq"
echo "  ${OUT_DIR}/partial_adapter_R1.fq"
echo "  ${OUT_DIR}/partial_adapter_R2.fq"
echo "  ${OUT_DIR}/adapter_all_R1.fq"
echo "  ${OUT_DIR}/adapter_all_R2.fq"
echo "  ${OUT_DIR}/adapter_negative_R1.fq"
echo "  ${OUT_DIR}/adapter_negative_R2.fq"
echo "  ${OUT_DIR}/adapter_multiple_R1.fq"
echo "  ${OUT_DIR}/adapter_multiple_R2.fq"
echo "  ${OUT_DIR}/adapter_short_R1.fq"
echo "  ${OUT_DIR}/adapter_short_R2.fq"