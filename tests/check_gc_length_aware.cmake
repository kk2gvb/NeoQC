if(NOT DEFINED REPORT)
    message(FATAL_ERROR "REPORT variable is required")
endif()

if(NOT EXISTS "${REPORT}")
    message(FATAL_ERROR "Report does not exist: ${REPORT}")
endif()

file(STRINGS "${REPORT}" LINES)

list(LENGTH LINES LINE_COUNT)

if(NOT LINE_COUNT EQUAL 102)
    message(FATAL_ERROR
        "Unexpected number of lines in ${REPORT}: "
        "expected 102, got ${LINE_COUNT}")
endif()

set(EXPECTED_HEADER "gc_percent\traw_read_count\tfastqc_observed_count")

list(GET LINES 0 HEADER)

if(NOT HEADER STREQUAL EXPECTED_HEADER)
    message(FATAL_ERROR
        "Unexpected GC TSV header:\n"
        "expected: '${EXPECTED_HEADER}'\n"
        "actual:   '${HEADER}'")
endif()

set(EXPECTED
    "0\t2\t2"
    "1\t0\t2"
    "12\t0\t2"
    "13\t0\t1"
    "94\t0\t0.5"
    "95\t0\t1"
    "99\t0\t1"
    "100\t1\t1"
)

foreach(EXPECTED_LINE IN LISTS EXPECTED)
    string(REGEX MATCH "^([^\t]+)\t([^\t]+)\t([^\t]+)$"
           MATCH "${EXPECTED_LINE}")

    set(GC_PERCENT "${CMAKE_MATCH_1}")

    set(FOUND FALSE)

    foreach(LINE IN LISTS LINES)
        if(LINE MATCHES "^${GC_PERCENT}\t")
            set(ACTUAL_LINE "${LINE}")
            set(FOUND TRUE)
            break()
        endif()
    endforeach()

    if(NOT FOUND)
        message(FATAL_ERROR
            "GC regression: missing bin ${GC_PERCENT}%")
    endif()

    if(NOT ACTUAL_LINE STREQUAL EXPECTED_LINE)
        message(FATAL_ERROR
            "GC regression at ${GC_PERCENT}%:\n"
            "expected: '${EXPECTED_LINE}'\n"
            "actual:   '${ACTUAL_LINE}'")
    endif()
endforeach()

message(STATUS "GC length-aware contract passed")
