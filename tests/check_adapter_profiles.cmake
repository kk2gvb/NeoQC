function(check_adapter_report PATH ADAPTER_NAME EXPECTED_VALUE)
    if(NOT EXISTS "${PATH}")
        message(FATAL_ERROR "Adapter report was not created: ${PATH}")
    endif()

    file(READ "${PATH}" CONTENT)

    if(NOT CONTENT MATCHES "${ADAPTER_NAME}")
        message(FATAL_ERROR "Missing ${ADAPTER_NAME} column in ${PATH}")
    endif()

    if(NOT CONTENT MATCHES "${EXPECTED_VALUE}")
        message(FATAL_ERROR
            "Expected ${EXPECTED_VALUE} adapter detection in ${PATH}")
    endif()
endfunction()

check_adapter_report("${R1_REPORT}" "TruSeq_R1" "${EXPECTED_VALUE}")
check_adapter_report("${R2_REPORT}" "TruSeq_R2" "${EXPECTED_VALUE}")
