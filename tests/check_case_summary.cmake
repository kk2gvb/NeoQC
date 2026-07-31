if(NOT EXISTS "${CASE_SUMMARY}")
    message(FATAL_ERROR "case_summary.json was not created: ${CASE_SUMMARY}")
endif()

file(READ "${CASE_SUMMARY}" SUMMARY)
foreach(expected
        "\"patient_id\": \"P001\""
        "\"status\": \"passed\""
        "\"sample_id\": \"TUMOR_DNA\""
        "\"sample_id\": \"NORMAL_DNA\""
        "\"r1_reads\": 2")
    string(FIND "${SUMMARY}" "${expected}" POSITION)
    if(POSITION EQUAL -1)
        message(FATAL_ERROR "Missing expected JSON value: ${expected}")
    endif()
endforeach()
