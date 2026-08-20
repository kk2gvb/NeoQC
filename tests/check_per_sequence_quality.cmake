if(NOT EXISTS "${REPORT}")
    message(FATAL_ERROR "Per-sequence quality report was not created: ${REPORT}")
endif()

file(STRINGS "${REPORT}" ROWS)

list(GET ROWS 0 HEADER)

message(STATUS "HEADER=[${HEADER}]")

if(NOT HEADER STREQUAL "mean_quality\tread_count\tread_count_truncate")
    message(FATAL_ERROR "Unexpected per-sequence quality TSV header")
endif()

set(TOTAL_READS 0)
set(TOTAL_READS_TRUNCATE 0)

foreach(ROW IN LISTS ROWS)
    if(ROW MATCHES "^[0-9]+\t[0-9]+\t[0-9]+$")
        string(REGEX REPLACE "^[0-9]+\t([0-9]+)\t[0-9]+$" "\\1" COUNT "${ROW}")
        string(REGEX REPLACE "^[0-9]+\t[0-9]+\t([0-9]+)$" "\\1" COUNT_TRUNCATE "${ROW}")

        math(EXPR TOTAL_READS "${TOTAL_READS} + ${COUNT}")
        math(EXPR TOTAL_READS_TRUNCATE "${TOTAL_READS_TRUNCATE} + ${COUNT_TRUNCATE}")
    endif()
endforeach()

message(STATUS "TOTAL_READS=[${TOTAL_READS}]")
message(STATUS "TOTAL_READS_TRUNCATE=[${TOTAL_READS_TRUNCATE}]")

if(NOT TOTAL_READS EQUAL 1)
    message(FATAL_ERROR "Expected one read in rounded distribution, got ${TOTAL_READS}")
endif()

if(NOT TOTAL_READS_TRUNCATE EQUAL 1)
    message(FATAL_ERROR "Expected one read in truncated distribution, got ${TOTAL_READS_TRUNCATE}")
endif()