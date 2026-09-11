function(check_adapter_report PATH ADAPTER_NAME EXPECTED_VALUE)
    if(NOT EXISTS "${PATH}")
        message(FATAL_ERROR
            "Adapter report was not created: ${PATH}")
    endif()

    file(STRINGS "${PATH}" LINES)

    if(LINES STREQUAL "")
        message(FATAL_ERROR
            "Adapter report is empty: ${PATH}")
    endif()

    # Header
    list(GET LINES 0 HEADER)

    string(REPLACE "\t" ";" COLUMNS "${HEADER}")

    list(FIND COLUMNS "${ADAPTER_NAME}" ADAPTER_INDEX)

    if(ADAPTER_INDEX EQUAL -1)
        message(FATAL_ERROR
            "Missing ${ADAPTER_NAME} column in ${PATH}")
    endif()

    # Search adapter column for the expected value.
    set(FOUND_EXPECTED FALSE)

    list(LENGTH LINES LINE_COUNT)

    if(LINE_COUNT GREATER 1)
        math(EXPR LAST_INDEX "${LINE_COUNT} - 1")

        foreach(INDEX RANGE 1 ${LAST_INDEX})
            list(GET LINES ${INDEX} LINE)

            if(LINE STREQUAL "")
                continue()
            endif()

            string(REPLACE "\t" ";" VALUES "${LINE}")

            list(LENGTH VALUES VALUE_COUNT)

            if(VALUE_COUNT LESS_EQUAL ADAPTER_INDEX)
                message(FATAL_ERROR
                    "Malformed row ${INDEX} in ${PATH}")
            endif()

            list(GET VALUES ${ADAPTER_INDEX} VALUE)

            if(VALUE STREQUAL "${EXPECTED_VALUE}")
                set(FOUND_EXPECTED TRUE)
                break()
            endif()
        endforeach()
    endif()

    if(NOT FOUND_EXPECTED)
        message(FATAL_ERROR
            "Expected ${EXPECTED_VALUE} for ${ADAPTER_NAME} in ${PATH}")
    endif()

    message(STATUS
        "Adapter ${ADAPTER_NAME}: ${EXPECTED_VALUE} found")
endfunction()


# ==============================================================================
# R1
# ==============================================================================

check_adapter_report(
    "${R1_REPORT}"
    "TruSeq_R1"
    "16.6667"
)

check_adapter_report(
    "${R1_REPORT}"
    "SmallRNA3'"
    "16.6667"
)

check_adapter_report(
    "${R1_REPORT}"
    "SmallRNA5'"
    "16.6667"
)

check_adapter_report(
    "${R1_REPORT}"
    "Nextera"
    "16.6667"
)

check_adapter_report(
    "${R1_REPORT}"
    "PolyA"
    "16.6667"
)

check_adapter_report(
    "${R1_REPORT}"
    "PolyG"
    "16.6667"
)


# ==============================================================================
# R2
# ==============================================================================

check_adapter_report(
    "${R2_REPORT}"
    "TruSeq_R2"
    "16.6667"
)

check_adapter_report(
    "${R2_REPORT}"
    "SmallRNA3'"
    "16.6667"
)

check_adapter_report(
    "${R2_REPORT}"
    "SmallRNA5'"
    "16.6667"
)

check_adapter_report(
    "${R2_REPORT}"
    "Nextera"
    "16.6667"
)

check_adapter_report(
    "${R2_REPORT}"
    "PolyA"
    "16.6667"
)

check_adapter_report(
    "${R2_REPORT}"
    "PolyG"
    "16.6667"
)