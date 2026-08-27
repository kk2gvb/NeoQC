if(NOT DEFINED OUTPUT_DIR)
    message(FATAL_ERROR "OUTPUT_DIR is not defined")
endif()

set(SUMMARY "${OUTPUT_DIR}/numeric_regression_R1_summary.txt")
set(PER_CYCLE "${OUTPUT_DIR}/per_cycle_R1.tsv")
set(PER_SEQUENCE_QUALITY "${OUTPUT_DIR}/per_sequence_quality_R1.tsv")
set(PER_SEQUENCE_GC "${OUTPUT_DIR}/per_sequence_gc_content_R1.tsv")
set(LENGTH_DISTRIBUTION "${OUTPUT_DIR}/sequence_length_distribution_R1.tsv")

foreach(path IN ITEMS
        "${SUMMARY}"
        "${PER_CYCLE}"
        "${PER_SEQUENCE_QUALITY}"
        "${PER_SEQUENCE_GC}"
        "${LENGTH_DISTRIBUTION}")

    if(NOT EXISTS "${path}")
        message(FATAL_ERROR "Missing regression artifact: ${path}")
    endif()

endforeach()


# ==============================================================================
# SUMMARY CONTRACT
# ==============================================================================

file(READ "${SUMMARY}" SUMMARY_CONTENT)

set(EXPECTED_SUMMARY_LINES
    "Processed reads : 4"
    "Total bases     : 16"
    "Min length      : 4"
    "Max length      : 4"
    "Avg length      : 4.00"

    "A               : 2"
    "C               : 4"
    "G               : 4"
    "T               : 2"
    "N               : 4"

    "GC content      : 50.00%"
    "%N              : 25.00%"
    "%Q20            : 50.00%"
    "%Q30            : 25.00%"
    "% with adapter  : 0.00%"
)

foreach(expected_line IN LISTS EXPECTED_SUMMARY_LINES)

    string(FIND "${SUMMARY_CONTENT}"
                "${expected_line}"
                FOUND_POSITION)

    if(FOUND_POSITION EQUAL -1)
        message(FATAL_ERROR
            "Numeric regression failed in summary.\n"
            "Expected line: ${expected_line}")
    endif()

endforeach()


# ==============================================================================
# PER-CYCLE QUALITY CONTRACT
# ==============================================================================

file(STRINGS "${PER_CYCLE}" PER_CYCLE_ROWS)

list(GET PER_CYCLE_ROWS 0 PER_CYCLE_HEADER)

if(NOT PER_CYCLE_HEADER STREQUAL
       "cycle\tmean_quality\tlower_quartile\tmedian")

    message(FATAL_ERROR
        "Unexpected per-cycle quality header: ${PER_CYCLE_HEADER}")
endif()

set(EXPECTED_PER_CYCLE_ROWS
    "1\t15.5000\t0.0000\t2.0000"
    "2\t15.5000\t0.0000\t2.0000"
    "3\t15.5000\t0.0000\t2.0000"
    "4\t15.5000\t0.0000\t2.0000"
)

foreach(expected_row IN LISTS EXPECTED_PER_CYCLE_ROWS)

    list(FIND PER_CYCLE_ROWS "${expected_row}" FOUND_INDEX)

    if(FOUND_INDEX EQUAL -1)
        message(FATAL_ERROR
            "Numeric regression failed in per-cycle quality.\n"
            "Expected row: ${expected_row}")
    endif()

endforeach()


# ==============================================================================
# PER-SEQUENCE QUALITY CONTRACT
# ==============================================================================

file(STRINGS "${PER_SEQUENCE_QUALITY}" QUALITY_ROWS)

list(GET QUALITY_ROWS 0 QUALITY_HEADER)

if(NOT QUALITY_HEADER STREQUAL
       "mean_quality\tread_count\tread_count_truncate")

    message(FATAL_ERROR
        "Unexpected per-sequence quality header: ${QUALITY_HEADER}")
endif()

set(EXPECTED_QUALITY_ROWS
    "0\t1\t1"
    "2\t1\t1"
    "20\t1\t1"
    "40\t1\t1"
)

foreach(expected_row IN LISTS EXPECTED_QUALITY_ROWS)

    list(FIND QUALITY_ROWS "${expected_row}" FOUND_INDEX)

    if(FOUND_INDEX EQUAL -1)
        message(FATAL_ERROR
            "Numeric regression failed in per-sequence quality.\n"
            "Expected row: ${expected_row}")
    endif()

endforeach()


# Проверяем, что суммарно распределение содержит ровно 4 reads.

set(TOTAL_QUALITY_READS 0)
set(TOTAL_QUALITY_READS_TRUNCATE 0)

foreach(row IN LISTS QUALITY_ROWS)

    if(row MATCHES "^[0-9]+\t[0-9]+\t[0-9]+$")

        string(REGEX REPLACE
               "^[0-9]+\t([0-9]+)\t[0-9]+$"
               "\\1"
               COUNT
               "${row}")

        string(REGEX REPLACE
               "^[0-9]+\t[0-9]+\t([0-9]+)$"
               "\\1"
               COUNT_TRUNCATE
               "${row}")

        math(EXPR TOTAL_QUALITY_READS
             "${TOTAL_QUALITY_READS} + ${COUNT}")

        math(EXPR TOTAL_QUALITY_READS_TRUNCATE
             "${TOTAL_QUALITY_READS_TRUNCATE} + ${COUNT_TRUNCATE}")

    endif()

endforeach()

if(NOT TOTAL_QUALITY_READS EQUAL 4)

    message(FATAL_ERROR
        "Per-sequence quality regression: expected 4 reads, "
        "got ${TOTAL_QUALITY_READS}")

endif()

if(NOT TOTAL_QUALITY_READS_TRUNCATE EQUAL 4)

    message(FATAL_ERROR
        "Per-sequence quality truncated regression: expected 4 reads, "
        "got ${TOTAL_QUALITY_READS_TRUNCATE}")

endif()


# ==============================================================================
# PER-SEQUENCE GC CONTRACT
# ==============================================================================

file(STRINGS "${PER_SEQUENCE_GC}" GC_ROWS)

list(GET GC_ROWS 0 GC_HEADER)

if(NOT GC_HEADER STREQUAL "gc_percent\treads")

    message(FATAL_ERROR
        "Unexpected per-sequence GC header: ${GC_HEADER}")

endif()

set(EXPECTED_GC_ROWS
    "0\t1"
    "50\t2"
    "100\t1"
)

foreach(expected_row IN LISTS EXPECTED_GC_ROWS)

    list(FIND GC_ROWS "${expected_row}" FOUND_INDEX)

    if(FOUND_INDEX EQUAL -1)
        message(FATAL_ERROR
            "Numeric regression failed in GC distribution.\n"
            "Expected row: ${expected_row}")

    endif()

endforeach()


# Проверяем общее количество reads в GC distribution.

set(TOTAL_GC_READS 0)

foreach(row IN LISTS GC_ROWS)

    if(row MATCHES "^[0-9]+\t[0-9]+$")

        string(REGEX REPLACE
               "^[0-9]+\t([0-9]+)$"
               "\\1"
               COUNT
               "${row}")

        math(EXPR TOTAL_GC_READS
             "${TOTAL_GC_READS} + ${COUNT}")

    endif()

endforeach()

if(NOT TOTAL_GC_READS EQUAL 4)

    message(FATAL_ERROR
        "GC distribution regression: expected 4 reads, "
        "got ${TOTAL_GC_READS}")

endif()


# ==============================================================================
# SEQUENCE LENGTH DISTRIBUTION CONTRACT
# ==============================================================================

file(STRINGS "${LENGTH_DISTRIBUTION}" LENGTH_ROWS)

list(GET LENGTH_ROWS 0 LENGTH_HEADER)

if(NOT LENGTH_HEADER STREQUAL "length\treads")

    message(FATAL_ERROR
        "Unexpected sequence length distribution header: ${LENGTH_HEADER}")

endif()

list(LENGTH LENGTH_ROWS LENGTH_ROW_COUNT)

if(NOT LENGTH_ROW_COUNT EQUAL 2)

    message(FATAL_ERROR
        "Sequence length regression: expected exactly one data row, "
        "got ${LENGTH_ROW_COUNT}")

endif()

list(GET LENGTH_ROWS 1 LENGTH_DATA)

if(NOT LENGTH_DATA STREQUAL "4\t4")

    message(FATAL_ERROR
        "Sequence length regression failed.\n"
        "Expected: 4<TAB>4\n"
        "Actual: ${LENGTH_DATA}")

endif()


message(STATUS
    "Numeric regression suite passed for ${OUTPUT_DIR}")