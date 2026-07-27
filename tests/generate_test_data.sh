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


# 1. Корректный одиночный
cat > tests/data/correct_single.fq << 'EOF'
@READ_001
AGCTTAGCCATGGCATAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGC
+
AAAAFFFFFJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ
EOF
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
echo "2.1Correct pair-end FASTQ file generated: tests/data/correct_pair_R1.fq.gz"

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

# 12. Повреждённый gzip
echo "not a gzip file" > tests/data/corrupted_gzip.fq.gz
echo "12. Corrupted gzip file generated: tests/data/corrupted_gzip.fq.gz"

echo
echo "====================================="
echo "Generated 12 FASTQ test datasets."
echo "====================================="