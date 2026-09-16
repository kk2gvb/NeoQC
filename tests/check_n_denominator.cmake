if(NOT DEFINED REPORT)
    message(FATAL_ERROR "REPORT variable is required")
endif()

if(NOT EXISTS "${REPORT}")
    message(FATAL_ERROR "Report does not exist: ${REPORT}")
endif()

file(STRINGS "${REPORT}" LINES)

set(EXPECTED
    "position\tN_percent"
    "1\t0.0000"
    "2\t0.0000"
    "3\t0.0000"
    "4\t50.0000"
    "5\t100.0000"
)

list(LENGTH EXPECTED EXPECTED_COUNT)
list(LENGTH LINES ACTUAL_COUNT)

if(NOT ACTUAL_COUNT EQUAL EXPECTED_COUNT)
    message(FATAL_ERROR
        "Unexpected number of lines in ${REPORT}: "
        "expected ${EXPECTED_COUNT}, got ${ACTUAL_COUNT}")
endif()

math(EXPR LAST_INDEX "${EXPECTED_COUNT} - 1")

foreach(INDEX RANGE ${LAST_INDEX})
    list(GET EXPECTED ${INDEX} EXPECTED_LINE)
    list(GET LINES ${INDEX} ACTUAL_LINE)

    if(NOT ACTUAL_LINE STREQUAL EXPECTED_LINE)
        message(FATAL_ERROR
            "N denominator regression at line ${INDEX}:\n"
            "expected: '${EXPECTED_LINE}'\n"
            "actual:   '${ACTUAL_LINE}'")
    endif()
endforeach()

message(STATUS "N denominator contract passed")