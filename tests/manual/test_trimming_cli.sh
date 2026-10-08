#!/usr/bin/env bash

set -u
set -o pipefail

# ============================================================
# NeoQC trimming CLI smoke / integration test
#
# Run from repository root:
#
#   ./tests/manual/test_trimming_cli.sh
#
# Or:
#
#   NEOQC=./build/neoqc ./tests/manual/test_trimming_cli.sh
#
# Temporary workspace:
#
#   build/manual_trimming_test/
#
# The script:
#   - creates deterministic synthetic FASTQ/TSV files;
#   - validates the fixtures before testing NeoQC;
#   - tests trimming CLI options;
#   - validates generated FASTQ;
#   - checks exact expected sequences where possible;
#   - checks PE synchronization;
#   - checks PE overlap correction;
#   - checks PE merging;
#   - checks regression without --trim.
# ============================================================

ROOT="${ROOT:-$(pwd)}"
NEOQC="${NEOQC:-${ROOT}/build/neoqc}"

WORK="${ROOT}/build/manual_trimming_test"
INPUT="${WORK}/input"
OUT="${WORK}/results"

PASS=0
FAIL=0
SKIP=0

# ============================================================
# Colors
# ============================================================

if [[ -t 1 ]]; then
    GREEN=$'\033[32m'
    RED=$'\033[31m'
    YELLOW=$'\033[33m'
    BLUE=$'\033[34m'
    BOLD=$'\033[1m'
    RESET=$'\033[0m'
else
    GREEN=""
    RED=""
    YELLOW=""
    BLUE=""
    BOLD=""
    RESET=""
fi

pass() {
    printf "%s[PASS]%s %s\n" "$GREEN" "$RESET" "$1"
    PASS=$((PASS + 1))
}

fail() {
    printf "%s[FAIL]%s %s\n" "$RED" "$RESET" "$1"
    FAIL=$((FAIL + 1))
}

skip() {
    printf "%s[SKIP]%s %s\n" "$YELLOW" "$RESET" "$1"
    SKIP=$((SKIP + 1))
}

section() {
    printf "\n%s%s=== %s ===%s\n" \
        "$BLUE" "$BOLD" "$1" "$RESET"
}

die() {
    echo "ERROR: $*" >&2
    exit 1
}

# ============================================================
# General helpers
# ============================================================

repeat_char() {
    local char="$1"
    local count="$2"

    if [[ "$count" -le 0 ]]; then
        printf ''
        return
    fi

    printf '%*s' "$count" '' | tr ' ' "$char"
}

reverse_complement() {
    local seq="$1"

    printf '%s' "$seq" |
        rev |
        tr 'ACGTNacgtn' 'TGCANtgcan'
}

write_fastq_record() {
    local file="$1"
    local id="$2"
    local seq="$3"

    local qual
    qual="$(repeat_char "I" "${#seq}")"

    {
        printf '@%s\n' "$id"
        printf '%s\n' "$seq"
        printf '+\n'
        printf '%s\n' "$qual"
    } >"$file"
}

append_fastq_record() {
    local file="$1"
    local id="$2"
    local seq="$3"

    local qual
    qual="$(repeat_char "I" "${#seq}")"

    {
        printf '@%s\n' "$id"
        printf '%s\n' "$seq"
        printf '+\n'
        printf '%s\n' "$qual"
    } >>"$file"
}

write_fastq_record_q() {
    local file="$1"
    local id="$2"
    local seq="$3"
    local qual="$4"

    if [[ ${#seq} -ne ${#qual} ]]; then
        die "Fixture '$id': sequence/quality mismatch: ${#seq} != ${#qual}"
    fi

    {
        printf '@%s\n' "$id"
        printf '%s\n' "$seq"
        printf '+\n'
        printf '%s\n' "$qual"
    } >"$file"
}

append_fastq_record_q() {
    local file="$1"
    local id="$2"
    local seq="$3"
    local qual="$4"

    if [[ ${#seq} -ne ${#qual} ]]; then
        die "Fixture '$id': sequence/quality mismatch: ${#seq} != ${#qual}"
    fi

    {
        printf '@%s\n' "$id"
        printf '%s\n' "$seq"
        printf '+\n'
        printf '%s\n' "$qual"
    } >>"$file"
}

read_text() {
    local file="$1"

    case "$file" in
        *.gz)
            gzip -cd "$file"
            ;;
        *)
            cat "$file"
            ;;
    esac
}

count_fastq_reads() {
    local file="$1"

    local lines
    lines="$(read_text "$file" | wc -l)"

    echo $((lines / 4))
}

first_header() {
    local file="$1"
    read_text "$file" | sed -n '1p'
}

first_sequence() {
    local file="$1"
    read_text "$file" | sed -n '2p'
}

first_quality() {
    local file="$1"
    read_text "$file" | sed -n '4p'
}

find_any_trimmed_fastq() {
    local dir="$1"

    find "$dir" -type f \
        \( \
        -name '*trimmed*.fastq' \
        -o \
        -name '*trimmed*.fastq.gz' \
        \) \
        | grep -vi 'merged' \
        | head -n 1
}

find_trimmed_r1() {
    local dir="$1"

    find "$dir" -type f \
        \( \
        -name '*R1*trimmed*.fastq' \
        -o \
        -name '*R1*trimmed*.fastq.gz' \
        -o \
        -name '*trimmed*R1*.fastq' \
        -o \
        -name '*trimmed*R1*.fastq.gz' \
        \) \
        | head -n 1
}

find_trimmed_r2() {
    local dir="$1"

    find "$dir" -type f \
        \( \
        -name '*R2*trimmed*.fastq' \
        -o \
        -name '*R2*trimmed*.fastq.gz' \
        -o \
        -name '*trimmed*R2*.fastq' \
        -o \
        -name '*trimmed*R2*.fastq.gz' \
        \) \
        | head -n 1
}

find_merged_fastq() {
    local dir="$1"

    find "$dir" -type f \
        \( \
        -name '*merged*.fastq' \
        -o \
        -name '*merged*.fastq.gz' \
        \) \
        | head -n 1
}

find_trimming_report() {
    local dir="$1"

    find "$dir" -type f \
        -name 'trimming_report.json' \
        | head -n 1
}

# ============================================================
# NeoQC execution helpers
# ============================================================

run_neoqc() {
    local name="$1"
    shift

    local log="${WORK}/${name}.log"

    echo
    echo ">>> $NEOQC $*"

    "$NEOQC" "$@" >"$log" 2>&1
    local rc=$?

    if [[ $rc -ne 0 ]]; then
        echo "--- log: $log ---"
        cat "$log"
    fi

    return "$rc"
}

run_expect_success() {
    local name="$1"
    shift

    if run_neoqc "$name" "$@"; then
        pass "$name"
        return 0
    fi

    fail "$name"
    return 1
}

run_expect_failure() {
    local name="$1"
    shift

    local log="${WORK}/${name}.log"

    echo
    echo ">>> EXPECT FAILURE: $NEOQC $*"

    "$NEOQC" "$@" >"$log" 2>&1
    local rc=$?

    if [[ $rc -ne 0 ]]; then
        pass "$name rejected as expected"
        return 0
    fi

    fail "$name unexpectedly succeeded"
    return 1
}

# ============================================================
# Assertions
# ============================================================

validate_fastq() {
    local file="$1"
    local name="$2"

    if [[ -z "$file" || ! -f "$file" ]]; then
        fail "$name: FASTQ does not exist"
        return 1
    fi

    local tmp="${WORK}/validate.$$.$RANDOM.fastq"

    if ! read_text "$file" >"$tmp"; then
        fail "$name: cannot read FASTQ"
        rm -f "$tmp"
        return 1
    fi

    if ! awk '
        NR % 4 == 1 {
            if (substr($0, 1, 1) != "@")
                exit 10
        }

        NR % 4 == 2 {
            seq = $0
        }

        NR % 4 == 3 {
            if (substr($0, 1, 1) != "+")
                exit 11
        }

        NR % 4 == 0 {
            if (length(seq) != length($0))
                exit 12
        }

        END {
            if (NR % 4 != 0)
                exit 13
        }
    ' "$tmp"; then
        fail "$name: invalid FASTQ structure"
        rm -f "$tmp"
        return 1
    fi

    rm -f "$tmp"

    pass "$name: valid FASTQ"
    return 0
}

expect_first_sequence() {
    local file="$1"
    local expected="$2"
    local name="$3"

    local actual
    actual="$(first_sequence "$file")"

    if [[ "$actual" == "$expected" ]]; then
        pass "$name: sequence = $expected"
    else
        fail "$name: expected '$expected', got '$actual'"
    fi
}

expect_read_count() {
    local file="$1"
    local expected="$2"
    local name="$3"

    local actual
    actual="$(count_fastq_reads "$file")"

    if [[ "$actual" -eq "$expected" ]]; then
        pass "$name: read count = $expected"
    else
        fail "$name: expected $expected reads, got $actual"
    fi
}

expect_header_contains() {
    local file="$1"
    local expected="$2"
    local name="$3"

    local header
    header="$(first_header "$file")"

    if [[ "$header" == *"$expected"* ]]; then
        pass "$name: header contains '$expected'"
    else
        fail "$name: header '$header' does not contain '$expected'"
    fi
}

expect_file_exists() {
    local file="$1"
    local name="$2"

    if [[ -n "$file" && -f "$file" ]]; then
        pass "$name: file exists"
    else
        fail "$name: file not found"
    fi
}

# ============================================================
# Environment
# ============================================================

section "Environment"

[[ -x "$NEOQC" ]] ||
    die "NeoQC executable not found: $NEOQC"

command -v gzip >/dev/null ||
    die "gzip not found"

command -v awk >/dev/null ||
    die "awk not found"

command -v find >/dev/null ||
    die "find not found"

command -v rev >/dev/null ||
    die "rev not found"

if command -v jq >/dev/null; then
    HAVE_JQ=1
    pass "jq available"
else
    HAVE_JQ=0
    skip "jq not installed; JSON syntax checks will be skipped"
fi

rm -rf "$WORK"

mkdir -p "$INPUT"
mkdir -p "$OUT"

pass "workspace created: $WORK"

# ============================================================
# Create synthetic data
# ============================================================

section "Create synthetic FASTQ"

# ------------------------------------------------------------
# clean.fastq
#
# 24 bp, high quality
# ------------------------------------------------------------

write_fastq_record \
    "${INPUT}/clean.fastq" \
    "clean_001" \
    "ACGTACGTACGTACGTACGTACGT"

# ------------------------------------------------------------
# fixed.fastq
#
# AAAA + CCCCCCCC + GGGG
#
# front 4 => CCCCCCCCGGGG
# tail  4 => AAAACCCCCCCC
# both    => CCCCCCCC
# ------------------------------------------------------------

write_fastq_record \
    "${INPUT}/fixed.fastq" \
    "fixed_001" \
    "AAAACCCCCCCCGGGG"

# ------------------------------------------------------------
# low_quality.fastq
#
# !!!! = Phred 0
# IIII = Phred 40
#
# expected with cut-front + cut-tail:
#
# ACGTACGTACGT
# ------------------------------------------------------------

write_fastq_record_q \
    "${INPUT}/low_quality.fastq" \
    "quality_001" \
    "TTTTACGTACGTACGTGGGG" \
    "!!!!IIIIIIIIIIII!!!!"

# ------------------------------------------------------------
# window.fastq
#
# 12 good bases + 8 poor-quality bases
#
# With:
#   cut-right
#   threshold 20
#   window 4
#
# expected according to current implementation:
#
# ACGTACGTACG
# ------------------------------------------------------------

write_fastq_record_q \
    "${INPUT}/window.fastq" \
    "window_001" \
    "ACGTACGTACGTTTTTTTTT" \
    "IIIIIIIIIIII!!!!!!!!"

# ------------------------------------------------------------
# adapter.fastq
#
# insert:
#   ACGTACGTACGT
#
# adapter:
#   AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
# ------------------------------------------------------------

ADAPTER="AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"

write_fastq_record \
    "${INPUT}/adapter.fastq" \
    "adapter_001" \
    "ACGTACGTACGT${ADAPTER}"

# ------------------------------------------------------------
# polyG.fastq
#
# expected:
# ACGTACGTACGT
# ------------------------------------------------------------

write_fastq_record \
    "${INPUT}/polyg.fastq" \
    "polyg_001" \
    "ACGTACGTACGTGGGGGGGGGG"

# ------------------------------------------------------------
# polyG internal
#
# Must NOT be removed by terminal polyG trimming.
# ------------------------------------------------------------

write_fastq_record \
    "${INPUT}/polyg_internal.fastq" \
    "polyg_internal_001" \
    "ACGTGGGGGGGGACGT"

# ------------------------------------------------------------
# polyX.fastq
#
# terminal A run
#
# expected:
# ACGTACGTACGT
# ------------------------------------------------------------

write_fastq_record \
    "${INPUT}/polyx.fastq" \
    "polyx_001" \
    "ACGTACGTACGTAAAAAAAAAA"

# ------------------------------------------------------------
# polyX internal
#
# Must remain unchanged.
# ------------------------------------------------------------

write_fastq_record \
    "${INPUT}/polyx_internal.fastq" \
    "polyx_internal_001" \
    "ACGTAAAAAAAACGT"

# ------------------------------------------------------------
# min_length.fastq
#
# read1 = 20 bp
# read2 = 8 bp
#
# min-length 15:
# read1 PASS
# read2 DISCARD
# ------------------------------------------------------------

write_fastq_record \
    "${INPUT}/min_length.fastq" \
    "long_001" \
    "ACGTACGTACGTACGTACGT"

append_fastq_record \
    "${INPUT}/min_length.fastq" \
    "short_001" \
    "ACGTACGT"

# ------------------------------------------------------------
# UMI
#
# UMI = AACCGGTT
# biological sequence = ACGTACGTACGT
# ------------------------------------------------------------

write_fastq_record \
    "${INPUT}/umi.fastq" \
    "umi_001" \
    "AACCGGTTACGTACGTACGT"

# ------------------------------------------------------------
# Adapter config TSV
# ------------------------------------------------------------

cat >"${INPUT}/adapters.tsv" <<EOF
IlluminaAdapter	${ADAPTER}
DummyAdapter	TTGGAATTCCGG
EOF

# ------------------------------------------------------------
# PE filtering dataset
#
# pair_001:
#   both 16 bp => survives min-length 12
#
# pair_002:
#   R1 16 bp
#   R2  8 bp
#   => whole pair excluded from paired output
# ------------------------------------------------------------

write_fastq_record \
    "${INPUT}/pe_filter_R1.fastq" \
    "pair_001/1" \
    "ACGTACGTACGTACGT"

append_fastq_record \
    "${INPUT}/pe_filter_R1.fastq" \
    "pair_002/1" \
    "TGCATGCATGCATGCA"

write_fastq_record \
    "${INPUT}/pe_filter_R2.fastq" \
    "pair_001/2" \
    "TGCATGCATGCATGCA"

append_fastq_record \
    "${INPUT}/pe_filter_R2.fastq" \
    "pair_002/2" \
    "ACGTACGT"

# ============================================================
# PE overlap fixtures
# ============================================================

# 48 bp deterministic insert.
#
# Split visually:
#
# ACGTTGCAAGTC
# GATCGTACGATG
# CTAGCTAACGTC
# AGTATCGGATCA
#
# Full overlap = 48 bp.
#
# This is intentionally > 30 bp so a typical min_overlap=30
# implementation has enough overlap.

PE_INSERT="ACGTTGCAAGTCGATCGTACGATGCTAGCTAACGTCAGTATCGGATCA"
PE_R2_PERFECT="$(reverse_complement "$PE_INSERT")"

# ------------------------------------------------------------
# Perfect-overlap pair for merging
# ------------------------------------------------------------

write_fastq_record \
    "${INPUT}/pe_merge_R1.fastq" \
    "merge_pair_001/1" \
    "$PE_INSERT"

write_fastq_record \
    "${INPUT}/pe_merge_R2.fastq" \
    "merge_pair_001/2" \
    "$PE_R2_PERFECT"

# ------------------------------------------------------------
# Correction pair
#
# R2 contains exactly one mismatch.
# That mismatch has Q=0.
#
# R1 has Q=40 everywhere.
#
# Expected after overlap correction:
#
# reverse_complement(R2) == R1
# ------------------------------------------------------------

MISMATCH_POS=10

ORIGINAL_BASE="${PE_R2_PERFECT:${MISMATCH_POS}:1}"

case "$ORIGINAL_BASE" in
    A) MUTATED_BASE="C" ;;
    C) MUTATED_BASE="G" ;;
    G) MUTATED_BASE="T" ;;
    T) MUTATED_BASE="A" ;;
    *)
        die "Unexpected nucleotide in PE fixture: $ORIGINAL_BASE"
        ;;
esac

PE_R2_MUTATED="${PE_R2_PERFECT:0:${MISMATCH_POS}}${MUTATED_BASE}${PE_R2_PERFECT:$((MISMATCH_POS + 1))}"

PE_Q1="$(repeat_char "I" "${#PE_INSERT}")"
PE_Q2="$(repeat_char "I" "${#PE_R2_MUTATED}")"

PE_Q2="${PE_Q2:0:${MISMATCH_POS}}!${PE_Q2:$((MISMATCH_POS + 1))}"

write_fastq_record_q \
    "${INPUT}/pe_correction_R1.fastq" \
    "correction_pair_001/1" \
    "$PE_INSERT" \
    "$PE_Q1"

write_fastq_record_q \
    "${INPUT}/pe_correction_R2.fastq" \
    "correction_pair_001/2" \
    "$PE_R2_MUTATED" \
    "$PE_Q2"

pass "synthetic FASTQ/TSV files created"

# ============================================================
# Validate ALL synthetic FASTQ before testing NeoQC
# ============================================================

section "Validate synthetic fixtures"

FIXTURE_ERRORS=0

while IFS= read -r fixture; do
    if validate_fastq \
        "$fixture" \
        "fixture $(basename "$fixture")"
    then
        :
    else
        FIXTURE_ERRORS=$((FIXTURE_ERRORS + 1))
    fi
done < <(
    find "$INPUT" -maxdepth 1 -type f -name '*.fastq' | sort
)

if [[ "$FIXTURE_ERRORS" -ne 0 ]]; then
    die "Synthetic fixture validation failed"
fi

# ============================================================
# CLI help
# ============================================================

section "CLI help"

if "$NEOQC" --help >"${WORK}/help.txt" 2>&1; then
    pass "--help"
else
    fail "--help"
fi

for option in \
    --trim \
    --trim-front \
    --trim-tail \
    --cut-front \
    --cut-tail \
    --cut-right \
    --quality-threshold \
    --window-size \
    --adapter-sequence \
    --adapter-config \
    --overlap-correction \
    --merge \
    --umi \
    --umi-length \
    --min-length \
    --trim-polyg \
    --trim-polyx \
    --polyg-min-length \
    --polyx-min-length
do
    if grep -q -- "$option" "${WORK}/help.txt"; then
        pass "help contains $option"
    else
        fail "help missing $option"
    fi
done

# ============================================================
# Regression without --trim
# ============================================================

section "Regression without --trim"

BASE_OUT="${OUT}/baseline"

if run_expect_success \
    "baseline_no_trim" \
    --r1 "${INPUT}/clean.fastq" \
    --sample-id baseline \
    --out "$BASE_OUT"
then
    :
fi

BASE_TRIMMED="$(
    find "$BASE_OUT" -type f \
        \( \
        -name '*trimmed*.fastq' \
        -o \
        -name '*trimmed*.fastq.gz' \
        \) \
        -print \
        -quit
)"

if [[ -n "$BASE_TRIMMED" ]]; then
    fail "without --trim: trimmed FASTQ unexpectedly exists: $BASE_TRIMMED"
else
    pass "without --trim: no trimmed FASTQ"
fi

BASE_REPORT="$(
    find "$BASE_OUT" -type f \
        -name 'trimming_report.json' \
        -print \
        -quit
)"

if [[ -n "$BASE_REPORT" ]]; then
    fail "without --trim: trimming report unexpectedly exists"
else
    pass "without --trim: no trimming report"
fi

if [[ -f "${BASE_OUT}/per_cycle_R1.tsv" ]]; then
    pass "without --trim: normal QC output exists"
else
    fail "without --trim: per_cycle_R1.tsv missing"
fi

if [[ -f "${BASE_OUT}/run_manifest.json" ]]; then
    pass "without --trim: run_manifest.json exists"
else
    fail "without --trim: run_manifest.json missing"
fi

# ============================================================
# --trim master switch
# ============================================================

section "--trim master switch"

run_expect_failure \
    "trim_front_without_trim" \
    --r1 "${INPUT}/clean.fastq" \
    --sample-id invalid \
    --out "${OUT}/invalid_trim_front" \
    --trim-front 4

run_expect_failure \
    "adapter_without_trim" \
    --r1 "${INPUT}/clean.fastq" \
    --sample-id invalid \
    --out "${OUT}/invalid_adapter" \
    --adapter-sequence "$ADAPTER"

run_expect_failure \
    "min_length_without_trim" \
    --r1 "${INPUT}/clean.fastq" \
    --sample-id invalid \
    --out "${OUT}/invalid_min_length" \
    --min-length 10

# ============================================================
# Basic --trim
# ============================================================

section "Basic trimming pipeline"

TEST_OUT="${OUT}/basic"

if run_expect_success \
    "basic_trim" \
    --r1 "${INPUT}/clean.fastq" \
    --sample-id basic \
    --out "$TEST_OUT" \
    --trim
then
    F="$(find_any_trimmed_fastq "$TEST_OUT")"

    if [[ -n "$F" ]]; then
        pass "basic trimmed FASTQ found: $F"

        validate_fastq "$F" "basic"

        expect_first_sequence \
            "$F" \
            "ACGTACGTACGTACGTACGTACGT" \
            "basic no-op trim"
    else
        fail "basic trimmed FASTQ not found"
    fi

    REPORT="$(find_trimming_report "$TEST_OUT")"

    if [[ -n "$REPORT" ]]; then
        pass "trimming report found: $REPORT"

        if [[ "$HAVE_JQ" -eq 1 ]]; then
            if jq empty "$REPORT" >/dev/null 2>&1; then
                pass "trimming report is valid JSON"
            else
                fail "trimming report is invalid JSON"
            fi
        fi
    else
        fail "trimming report not found"
    fi

    if [[ -d "${TEST_OUT}/qc/before" ]]; then
        pass "QC BEFORE directory exists"
    else
        fail "QC BEFORE directory missing"
    fi

    if [[ -d "${TEST_OUT}/qc/after" ]]; then
        pass "QC AFTER directory exists"
    else
        fail "QC AFTER directory missing"
    fi

    if [[ -f "${TEST_OUT}/qc/before/per_cycle_R1.tsv" ]]; then
        pass "QC BEFORE per-cycle output exists"
    else
        fail "QC BEFORE per_cycle_R1.tsv missing"
    fi

    if [[ -f "${TEST_OUT}/qc/after/per_cycle_R1.tsv" ]]; then
        pass "QC AFTER per-cycle output exists"
    else
        fail "QC AFTER per_cycle_R1.tsv missing"
    fi
fi

# ============================================================
# Fixed front
# ============================================================

section "Fixed front trimming"

TEST_OUT="${OUT}/fixed_front"

if run_expect_success \
    "fixed_front" \
    --r1 "${INPUT}/fixed.fastq" \
    --sample-id fixed_front \
    --out "$TEST_OUT" \
    --trim \
    --trim-front 4
then
    F="$(find_any_trimmed_fastq "$TEST_OUT")"

    if [[ -n "$F" ]]; then
        validate_fastq "$F" "fixed front"

        expect_first_sequence \
            "$F" \
            "CCCCCCCCGGGG" \
            "fixed front"
    else
        fail "fixed front output not found"
    fi
fi

# ============================================================
# Fixed tail
# ============================================================

section "Fixed tail trimming"

TEST_OUT="${OUT}/fixed_tail"

if run_expect_success \
    "fixed_tail" \
    --r1 "${INPUT}/fixed.fastq" \
    --sample-id fixed_tail \
    --out "$TEST_OUT" \
    --trim \
    --trim-tail 4
then
    F="$(find_any_trimmed_fastq "$TEST_OUT")"

    if [[ -n "$F" ]]; then
        validate_fastq "$F" "fixed tail"

        expect_first_sequence \
            "$F" \
            "AAAACCCCCCCC" \
            "fixed tail"
    else
        fail "fixed tail output not found"
    fi
fi

# ============================================================
# Fixed front + tail
# ============================================================

section "Fixed front + tail"

TEST_OUT="${OUT}/fixed_both"

if run_expect_success \
    "fixed_both" \
    --r1 "${INPUT}/fixed.fastq" \
    --sample-id fixed_both \
    --out "$TEST_OUT" \
    --trim \
    --trim-front 4 \
    --trim-tail 4
then
    F="$(find_any_trimmed_fastq "$TEST_OUT")"

    if [[ -n "$F" ]]; then
        validate_fastq "$F" "fixed both"

        expect_first_sequence \
            "$F" \
            "CCCCCCCC" \
            "fixed front + tail"
    else
        fail "fixed both output not found"
    fi
fi

# ============================================================
# Quality front/tail
# ============================================================

section "Quality front/tail"

TEST_OUT="${OUT}/quality"

if run_expect_success \
    "quality_front_tail" \
    --r1 "${INPUT}/low_quality.fastq" \
    --sample-id quality \
    --out "$TEST_OUT" \
    --trim \
    --cut-front \
    --cut-tail \
    --quality-threshold 20
then
    F="$(find_any_trimmed_fastq "$TEST_OUT")"

    if [[ -n "$F" ]]; then
        validate_fastq "$F" "quality"

        expect_first_sequence \
            "$F" \
            "ACGTACGTACGT" \
            "quality front/tail"
    else
        fail "quality output not found"
    fi
fi

# ============================================================
# Sliding-window
# ============================================================

section "Sliding-window trimming"

TEST_OUT="${OUT}/window"

if run_expect_success \
    "sliding_window" \
    --r1 "${INPUT}/window.fastq" \
    --sample-id window \
    --out "$TEST_OUT" \
    --trim \
    --cut-right \
    --quality-threshold 20 \
    --window-size 4
then
    F="$(find_any_trimmed_fastq "$TEST_OUT")"

    if [[ -n "$F" ]]; then
        validate_fastq "$F" "sliding window"

        expect_first_sequence \
            "$F" \
            "ACGTACGTACG" \
            "sliding window"
    else
        fail "sliding-window output not found"
    fi
fi

# ============================================================
# Explicit adapter
# ============================================================

section "Explicit adapter trimming"

TEST_OUT="${OUT}/adapter"

if run_expect_success \
    "explicit_adapter" \
    --r1 "${INPUT}/adapter.fastq" \
    --sample-id adapter \
    --out "$TEST_OUT" \
    --trim \
    --adapter-sequence "$ADAPTER"
then
    F="$(find_any_trimmed_fastq "$TEST_OUT")"

    if [[ -n "$F" ]]; then
        validate_fastq "$F" "explicit adapter"

        expect_first_sequence \
            "$F" \
            "ACGTACGTACGT" \
            "explicit adapter"
    else
        fail "explicit adapter output not found"
    fi
fi

# ============================================================
# Adapter config TSV
# ============================================================

section "Adapter config TSV"

TEST_OUT="${OUT}/adapter_config"

if run_expect_success \
    "adapter_config" \
    --r1 "${INPUT}/adapter.fastq" \
    --sample-id adapter_config \
    --out "$TEST_OUT" \
    --trim \
    --adapter-config "${INPUT}/adapters.tsv"
then
    F="$(find_any_trimmed_fastq "$TEST_OUT")"

    if [[ -n "$F" ]]; then
        validate_fastq "$F" "adapter config"

        expect_first_sequence \
            "$F" \
            "ACGTACGTACGT" \
            "adapter config"
    else
        fail "adapter config output not found"
    fi
fi

# ============================================================
# Minimum length
# ============================================================

section "Minimum length filtering"

TEST_OUT="${OUT}/min_length"

if run_expect_success \
    "min_length" \
    --r1 "${INPUT}/min_length.fastq" \
    --sample-id min_length \
    --out "$TEST_OUT" \
    --trim \
    --min-length 15
then
    F="$(find_any_trimmed_fastq "$TEST_OUT")"

    if [[ -n "$F" ]]; then
        validate_fastq "$F" "minimum length"

        expect_read_count \
            "$F" \
            1 \
            "minimum length"

        expect_first_sequence \
            "$F" \
            "ACGTACGTACGTACGTACGT" \
            "minimum length surviving read"
    else
        fail "minimum-length output not found"
    fi
fi

# ============================================================
# PolyG
# ============================================================

section "PolyG trimming"

TEST_OUT="${OUT}/polyg"

if run_expect_success \
    "polyg" \
    --r1 "${INPUT}/polyg.fastq" \
    --sample-id polyg \
    --out "$TEST_OUT" \
    --trim \
    --trim-polyg \
    --polyg-min-length 6
then
    F="$(find_any_trimmed_fastq "$TEST_OUT")"

    if [[ -n "$F" ]]; then
        validate_fastq "$F" "polyG"

        expect_first_sequence \
            "$F" \
            "ACGTACGTACGT" \
            "polyG"
    else
        fail "polyG output not found"
    fi
fi

# ============================================================
# PolyG internal should remain
# ============================================================

section "PolyG internal sequence"

TEST_OUT="${OUT}/polyg_internal"

if run_expect_success \
    "polyg_internal" \
    --r1 "${INPUT}/polyg_internal.fastq" \
    --sample-id polyg_internal \
    --out "$TEST_OUT" \
    --trim \
    --trim-polyg \
    --polyg-min-length 6
then
    F="$(find_any_trimmed_fastq "$TEST_OUT")"

    if [[ -n "$F" ]]; then
        expect_first_sequence \
            "$F" \
            "ACGTGGGGGGGGACGT" \
            "internal polyG untouched"
    else
        fail "internal polyG output not found"
    fi
fi

# ============================================================
# PolyX
# ============================================================

section "PolyX trimming"

TEST_OUT="${OUT}/polyx"

if run_expect_success \
    "polyx" \
    --r1 "${INPUT}/polyx.fastq" \
    --sample-id polyx \
    --out "$TEST_OUT" \
    --trim \
    --trim-polyx \
    --polyx-min-length 6
then
    F="$(find_any_trimmed_fastq "$TEST_OUT")"

    if [[ -n "$F" ]]; then
        validate_fastq "$F" "polyX"

        expect_first_sequence \
            "$F" \
            "ACGTACGTACGT" \
            "polyX"
    else
        fail "polyX output not found"
    fi
fi

# ============================================================
# PolyX internal should remain
# ============================================================

section "PolyX internal sequence"

TEST_OUT="${OUT}/polyx_internal"

if run_expect_success \
    "polyx_internal" \
    --r1 "${INPUT}/polyx_internal.fastq" \
    --sample-id polyx_internal \
    --out "$TEST_OUT" \
    --trim \
    --trim-polyx \
    --polyx-min-length 6
then
    F="$(find_any_trimmed_fastq "$TEST_OUT")"

    if [[ -n "$F" ]]; then
        expect_first_sequence \
            "$F" \
            "ACGTAAAAAAAACGT" \
            "internal polyX untouched"
    else
        fail "internal polyX output not found"
    fi
fi

# ============================================================
# UMI
# ============================================================

section "UMI preprocessing"

TEST_OUT="${OUT}/umi"

if run_expect_success \
    "umi" \
    --r1 "${INPUT}/umi.fastq" \
    --sample-id umi \
    --out "$TEST_OUT" \
    --trim \
    --umi \
    --umi-length 8
then
    F="$(find_any_trimmed_fastq "$TEST_OUT")"

    if [[ -n "$F" ]]; then
        validate_fastq "$F" "UMI"

        expect_first_sequence \
            "$F" \
            "ACGTACGTACGT" \
            "UMI extraction"

        expect_header_contains \
            "$F" \
            "UMI:AACCGGTT" \
            "UMI metadata"
    else
        fail "UMI output not found"
    fi
fi

# ============================================================
# UMI validation
# ============================================================

section "UMI validation"

run_expect_failure \
    "umi_without_length" \
    --r1 "${INPUT}/umi.fastq" \
    --sample-id invalid_umi \
    --out "${OUT}/invalid_umi" \
    --trim \
    --umi

# ============================================================
# Paired-end filtering synchronization
# ============================================================

section "Paired-end filtering synchronization"

TEST_OUT="${OUT}/pe_filter"

if run_expect_success \
    "pe_filter" \
    --r1 "${INPUT}/pe_filter_R1.fastq" \
    --r2 "${INPUT}/pe_filter_R2.fastq" \
    --sample-id pe_filter \
    --out "$TEST_OUT" \
    --trim \
    --min-length 12
then
    R1="$(find_trimmed_r1 "$TEST_OUT")"
    R2="$(find_trimmed_r2 "$TEST_OUT")"

    if [[ -n "$R1" && -n "$R2" ]]; then
        validate_fastq "$R1" "PE filtering R1"
        validate_fastq "$R2" "PE filtering R2"

        C1="$(count_fastq_reads "$R1")"
        C2="$(count_fastq_reads "$R2")"

        if [[ "$C1" -eq "$C2" ]]; then
            pass "PE R1/R2 read counts synchronized ($C1)"
        else
            fail "PE R1/R2 count mismatch: R1=$C1 R2=$C2"
        fi

        if [[ "$C1" -eq 1 ]]; then
            pass "PE pair-level filtering removed failed pair"
        else
            fail "expected 1 surviving PE pair, got $C1"
        fi

        H1="$(
            first_header "$R1" |
                sed \
                    -e 's#^@##' \
                    -e 's#/1.*##'
        )"

        H2="$(
            first_header "$R2" |
                sed \
                    -e 's#^@##' \
                    -e 's#/2.*##'
        )"

        if [[ "$H1" == "$H2" ]]; then
            pass "PE surviving IDs correspond"
        else
            fail "PE surviving IDs differ: '$H1' vs '$H2'"
        fi
    else
        fail "PE trimmed R1/R2 outputs not found"

        echo "Files produced:"
        find "$TEST_OUT" -type f -print
    fi
fi

# ============================================================
# PE overlap correction
# ============================================================

section "PE overlap correction"

TEST_OUT="${OUT}/correction"

if run_expect_success \
    "overlap_correction" \
    --r1 "${INPUT}/pe_correction_R1.fastq" \
    --r2 "${INPUT}/pe_correction_R2.fastq" \
    --sample-id correction \
    --out "$TEST_OUT" \
    --trim \
    --overlap-correction
then
    R1="$(find_trimmed_r1 "$TEST_OUT")"
    R2="$(find_trimmed_r2 "$TEST_OUT")"

    if [[ -n "$R1" && -n "$R2" ]]; then
        validate_fastq "$R1" "overlap correction R1"
        validate_fastq "$R2" "overlap correction R2"

        OUT_R1_SEQ="$(first_sequence "$R1")"
        OUT_R2_SEQ="$(first_sequence "$R2")"

        OUT_R2_RC="$(reverse_complement "$OUT_R2_SEQ")"

        if [[ "$OUT_R1_SEQ" == "$PE_INSERT" ]]; then
            pass "overlap correction: high-quality R1 preserved"
        else
            fail "overlap correction: R1 unexpectedly changed: $OUT_R1_SEQ"
        fi

        if [[ "$OUT_R2_RC" == "$OUT_R1_SEQ" ]]; then
            pass "overlap correction resolved PE mismatch"
        else
            fail "overlap correction mismatch remains"
            echo "R1       : $OUT_R1_SEQ"
            echo "RC(R2)   : $OUT_R2_RC"
        fi
    else
        fail "overlap correction R1/R2 output not found"

        echo "Files produced:"
        find "$TEST_OUT" -type f -print
    fi
fi

# ============================================================
# PE merge
# ============================================================

section "PE merging"

TEST_OUT="${OUT}/merge"

if run_expect_success \
    "merge" \
    --r1 "${INPUT}/pe_merge_R1.fastq" \
    --r2 "${INPUT}/pe_merge_R2.fastq" \
    --sample-id merge \
    --out "$TEST_OUT" \
    --trim \
    --merge
then
    MERGED="$(find_merged_fastq "$TEST_OUT")"

    if [[ -n "$MERGED" ]]; then
        pass "merged FASTQ found: $MERGED"

        validate_fastq "$MERGED" "merged FASTQ"

        expect_read_count \
            "$MERGED" \
            1 \
            "merged FASTQ"

        expect_first_sequence \
            "$MERGED" \
            "$PE_INSERT" \
            "merged full-overlap consensus"
    else
        fail "merge command succeeded but merged FASTQ was not found"

        echo "Files produced:"
        find "$TEST_OUT" -type f -print
    fi
fi

# ============================================================
# Invalid configuration
# ============================================================

section "Invalid configuration"

run_expect_failure \
    "window_zero" \
    --r1 "${INPUT}/clean.fastq" \
    --sample-id invalid_window \
    --out "${OUT}/invalid_window" \
    --trim \
    --cut-right \
    --window-size 0

run_expect_failure \
    "umi_zero" \
    --r1 "${INPUT}/clean.fastq" \
    --sample-id invalid_umi_zero \
    --out "${OUT}/invalid_umi_zero" \
    --trim \
    --umi \
    --umi-length 0

run_expect_failure \
    "missing_adapter_config" \
    --r1 "${INPUT}/clean.fastq" \
    --sample-id missing_adapter \
    --out "${OUT}/missing_adapter" \
    --trim \
    --adapter-config "${INPUT}/does_not_exist.tsv"

run_expect_failure \
    "correction_single_end" \
    --r1 "${INPUT}/clean.fastq" \
    --sample-id invalid_correction_se \
    --out "${OUT}/invalid_correction_se" \
    --trim \
    --overlap-correction

run_expect_failure \
    "merge_single_end" \
    --r1 "${INPUT}/clean.fastq" \
    --sample-id invalid_merge_se \
    --out "${OUT}/invalid_merge_se" \
    --trim \
    --merge

# ============================================================
# Inspect trimming reports
# ============================================================

section "Trimming reports"

REPORT_COUNT="$(
    find "$OUT" -type f \
        -name 'trimming_report.json' |
        wc -l
)"

echo "Found trimming reports: $REPORT_COUNT"

if [[ "$REPORT_COUNT" -gt 0 ]]; then
    pass "trimming reports generated"
else
    fail "no trimming reports generated"
fi

if [[ "$HAVE_JQ" -eq 1 ]]; then
    JSON_ERRORS=0

    while IFS= read -r report; do
        if jq empty "$report" >/dev/null 2>&1; then
            :
        else
            echo "Invalid JSON: $report"
            JSON_ERRORS=$((JSON_ERRORS + 1))
        fi
    done < <(
        find "$OUT" -type f \
            -name 'trimming_report.json' |
            sort
    )

    if [[ "$JSON_ERRORS" -eq 0 ]]; then
        pass "all trimming reports are valid JSON"
    else
        fail "$JSON_ERRORS trimming report(s) contain invalid JSON"
    fi
fi

# ============================================================
# Print compact result tree
# ============================================================

section "Generated main outputs"

find "$OUT" -type f \
    \( \
    -name '*.trimmed.fastq.gz' \
    -o \
    -name '*merged*.fastq.gz' \
    -o \
    -name 'trimming_report.json' \
    -o \
    -name 'run_manifest.json' \
    \) \
    -print |
    sort

# ============================================================
# Summary
# ============================================================

section "SUMMARY"

TOTAL=$((PASS + FAIL + SKIP))

printf "Total checks : %d\n" "$TOTAL"
printf "%sPassed       : %d%s\n" "$GREEN" "$PASS" "$RESET"
printf "%sFailed       : %d%s\n" "$RED" "$FAIL" "$RESET"
printf "%sSkipped      : %d%s\n" "$YELLOW" "$SKIP" "$RESET"

echo
echo "Workspace:"
echo "  $WORK"

echo
echo "Logs:"
echo "  ${WORK}/*.log"

echo
echo "Results:"
echo "  $OUT"

echo
echo "Useful commands after the run:"
echo
echo "  find \"$OUT\" -type f | sort"
echo
echo "  find \"$OUT\" -name 'trimming_report.json' -print"
echo
echo "  gzip -cd \"$OUT/basic/trimmed/R1.trimmed.fastq.gz\""
echo

if [[ "$FAIL" -ne 0 ]]; then
    echo
    printf "%sNeoQC trimming smoke test FAILED.%s\n" \
        "$RED" "$RESET"
    exit 1
fi

echo
printf "%sNeoQC trimming smoke test PASSED.%s\n" \
    "$GREEN" "$RESET"

exit 0
