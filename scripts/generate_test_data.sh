#!/usr/bin/env bash

set -euo pipefail

command -v gzip >/dev/null || {
    echo "gzip not found"
    exit 1
}

command -v python3 >/dev/null || {
    echo "python3 not found"
    exit 1
}

mkdir -p tests/data
rm -f tests/data/*


# 1. Корректный одиночный FASTQ: 10 reads
: > tests/data/correct_single.fq

for i in $(seq 1 10); do
    cat >> tests/data/correct_single.fq << EOF
@READ_$(printf "%03d" "$i")
AGCTTAGCCATGGCATAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGC
+
AAAAFFFFFJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ
EOF
done

gzip -f tests/data/correct_single.fq
echo "1. Correct single-end FASTQ file generated: tests/data/correct_single.fq.gz"


# 2. Корректная пара
cat > tests/data/correct_pair_R1.fq << 'EOF'
@READ_001/1
AGCTTAGCCATGGCATAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGC
+
AAAAFFFFFJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ
@READ_002/1
TGCAAGCTTAGCCATGGCATAGCTAGCTAGCTAGCTAGCTAGCTAGC
+
AAAAFFFFFJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ
EOF
gzip -f tests/data/correct_pair_R1.fq
echo "2.1 Correct pair-end FASTQ file generated: tests/data/correct_pair_R1.fq.gz"

cat > tests/data/correct_pair_R2.fq << 'EOF'
@READ_001/2
AGCTTAGCCATGGCATAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGC
+
AAAAFFFFFJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ
@READ_002/2
TGCAAGCTTAGCCATGGCATAGCTAGCTAGCTAGCTAGCTAGCTAGC
+
AAAAFFFFFJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ
EOF
gzip -f tests/data/correct_pair_R2.fq
echo "2.2 Correct pair-end FASTQ file generated: tests/data/correct_pair_R2.fq.gz"

# 3. Без @
cat > tests/data/error_no_at.fq << 'EOF'
READ_001
AGCTTAGC
+
AAAAFFFF
EOF
gzip -f tests/data/error_no_at.fq
echo "3. Error file generated (no @): tests/data/error_no_at.fq.gz"

# 4. Без +
cat > tests/data/error_no_plus.fq << 'EOF'
@READ_001
AGCTTAGC
SEQ
AAAAFFFF
EOF
gzip -f tests/data/error_no_plus.fq
echo "4. Error file generated (no +): tests/data/error_no_plus.fq.gz"

# 5. Разная длина
cat > tests/data/error_diff_len.fq << 'EOF'
@READ_001
AGCTAGCT
+
AAAAFFFFF
EOF
gzip -f tests/data/error_diff_len.fq
echo "5. Error file generated (different lengths): tests/data/error_diff_len.fq.gz"

# 6. Незаконченная
cat > tests/data/error_incomplete.fq << 'EOF'
@READ_001
AGCT
+
AAAA
@READ_002
EOF
gzip -f tests/data/error_incomplete.fq
echo "6. Error file generated (incomplete records): tests/data/error_incomplete.fq.gz"

# 7. Пустая строка
cat > tests/data/error_empty_line.fq << 'EOF'
@READ_001
AGCT

+
AAAA
EOF
gzip -f tests/data/error_empty_line.fq
echo "7. Error file generated (empty line): tests/data/error_empty_line.fq.gz"

# 8. Разное число записей
cat > tests/data/diff_counts_R1.fq << 'EOF'
@READ_001/1
AGCTTAGC
+
AAAAFFFF
@READ_002/1
TGCAAGCT
+
AAAAFFFF
@READ_003/1
GGCCATTA
+
AAAAFFFF
EOF
gzip -f tests/data/diff_counts_R1.fq
echo "8.1 Different counts file generated (R1): tests/data/diff_counts_R1.fq.gz"

cat > tests/data/diff_counts_R2.fq << 'EOF'
@READ_001/2
AGCTTAGC
+
AAAAFFFF
@READ_002/2
TGCAAGCT
+
AAAAFFFF
EOF
gzip -f tests/data/diff_counts_R2.fq
echo "8.2 Different counts file generated (R2): tests/data/diff_counts_R2.fq.gz"

# 9. Несовпадающие ID
cat > tests/data/mismatch_ids_R1.fq << 'EOF'
@READ_001/1
AGCTTAGC
+
AAAAFFFF
@READ_003/1
TGCAAGCT
+
AAAAFFFF
EOF
gzip -f tests/data/mismatch_ids_R1.fq
echo "9.1 Mismatching IDs file generated (R1): tests/data/mismatch_ids_R1.fq.gz"

cat > tests/data/mismatch_ids_R2.fq << 'EOF'
@READ_001/2
AGCTTAGC
+
AAAAFFFF
@READ_002/2
TGCAAGCT
+
AAAAFFFF
EOF
gzip -f tests/data/mismatch_ids_R2.fq
echo "9.2 Mismatching IDs file generated (R2): tests/data/mismatch_ids_R2.fq.gz"

# 10. Длинная строка >4096
python3 -c "
seq = 'A' * 5000
qual = 'I' * 5000
print('@LONG_READ')
print(seq)
print('+')
print(qual)
" > tests/data/long_line_4096.fq
gzip -f tests/data/long_line_4096.fq
echo "10. Long line file generated: tests/data/long_line_4096.fq.gz"

# 11. Windows CRLF
printf '@READ_001\r\nAGCTTAGC\r\n+\r\nAAAAFFFF\r\n' > tests/data/windows_clrf.fq
gzip -f tests/data/windows_clrf.fq
echo "11. Windows CRLF file generated: tests/data/windows_clrf.fq.gz"

# 12. Corrupted gzip (truncated)
cp tests/data/correct_single.fq.gz tests/data/corrupted_crc.fq.gz
truncate -s -8 tests/data/corrupted_crc.fq.gz

echo "12. Corrupted gzip (CRC) generated."

# 13. Пустой FASTQ
touch tests/data/empty.fq
gzip -f tests/data/empty.fq
echo "13. Empty FASTQ generated: tests/data/empty.fq.gz"

# 14. Пустая последовательность
cat > tests/data/empty_sequence.fq << 'EOF'
@READ_001

+
IIII
EOF
gzip -f tests/data/empty_sequence.fq
echo "14. Empty sequence generated: tests/data/empty_sequence.fq.gz"

# 15. Пустая качество
cat > tests/data/empty_quality.fq << 'EOF'
@READ_001
AGCT
+

EOF
gzip -f tests/data/empty_quality.fq
echo "15. Empty quality generated: tests/data/empty_quality.fq.gz"

# 16. Пустая строка между записями
cat > tests/data/blank_line_between_records.fq << 'EOF'
@READ_001
AGCT
+
IIII

@READ_002
TGCA
+
IIII
EOF
gzip -f tests/data/blank_line_between_records.fq
echo "16. Blank line between records generated."

# 17. Несовпадающие пары
cat > tests/data/pair_mismatch_R1.fq << 'EOF'
@READ1/1
ACGT
+
IIII
@READ2/1
ACGT
+
IIII
EOF
gzip -f tests/data/pair_mismatch_R1.fq

cat > tests/data/pair_mismatch_R2.fq << 'EOF'
@READ1/2
ACGT
+
IIII
@READ3/2
ACGT
+
IIII
EOF
gzip -f tests/data/pair_mismatch_R2.fq
echo "17. Mismatching pairs generated."

########################################
# 18. Invalid quality ASCII 32
########################################

python3 - << 'PY'
from pathlib import Path
import gzip

path = Path("tests/data/invalid_quality_ascii_32.fq.gz")

with gzip.open(path, "wb") as f:
    f.write(b"@READ_001\n")
    f.write(b"ACGT\n")
    f.write(b"+\n")
    f.write(bytes([32, 32, 32, 32]) + b"\n")
PY

echo "18. Invalid quality ASCII 32 generated."

########################################
# 19. Invalid quality ASCII 127
########################################

python3 - << 'PY'
from pathlib import Path
import gzip

path = Path("tests/data/invalid_quality_ascii_127.fq.gz")

with gzip.open(path, "wb") as f:
    f.write(b"@READ_001\n")
    f.write(b"ACGT\n")
    f.write(b"+\n")
    f.write(bytes([127, 127, 127, 127]) + b"\n")
PY

echo "19. Invalid quality ASCII 127 generated."

########################################
# 20. Invalid base character
########################################

cat > tests/data/invalid_base.fq << 'EOF'
@READ_001
ACGX
+
IIII
EOF

echo "20. Invalid base character generated."

########################################
# 21. Phred+33 boundary values
########################################

python3 - << 'PY'
from pathlib import Path
import gzip

path = Path("tests/data/quality_phred33_boundaries.fq.gz")

with gzip.open(path, "wb") as f:
    f.write(b"@READ_001\n")
    f.write(b"ACGT\n")
    f.write(b"+\n")
    f.write(bytes([33, 34, 125, 126]) + b"\n")
PY

echo "21. Phred+33 boundary dataset generated."

########################################
# 22. Medium dataset
########################################
echo "22. Generating medium paired-end dataset..."

python3 scripts/generate_fastq.py \
    --paired \
    --gzip \
    --reads 100000 \
    --length 150 \
    --name medium \
    --dir tests/data

########################################
# 23. Large dataset
########################################
echo "23. Generating large paired-end dataset..."

python3 scripts/generate_fastq.py \
    --paired \
    --gzip \
    --reads 1000000 \
    --length 150 \
    --name large \
    --dir tests/data

########################################
# 24. Adapter test datasets
########################################

echo
echo "Generating adapter test datasets..."

bash scripts/generate_adapter_tests.sh

echo
echo "====================================="

echo "All test datasets generated successfully."
echo "====================================="
