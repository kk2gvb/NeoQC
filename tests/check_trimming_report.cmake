if(NOT EXISTS "${REPORT}")
    message(FATAL_ERROR "Trimming report was not created: ${REPORT}")
endif()

file(READ "${REPORT}" REPORT_CONTENT)

foreach(expected
        "\"total\": 10"
        "\"passed\": 10"
        "\"discarded\": 0"
        "\"before\": 470"
        "\"after\": 470"
        "\"trimmed\": 0"
        "\"adapter_positions\": {}")
    string(FIND "${REPORT_CONTENT}" "${expected}" expected_position)
    if(expected_position EQUAL -1)
        message(FATAL_ERROR "Missing expected trimming report value: ${expected}")
    endif()
endforeach()

if(NOT EXISTS "${MANIFEST}")
    message(FATAL_ERROR "Run manifest was not created: ${MANIFEST}")
endif()

file(READ "${MANIFEST}" MANIFEST_CONTENT)
string(FIND "${MANIFEST_CONTENT}" "trimming_report.json" artifact_position)
if(artifact_position EQUAL -1)
    message(FATAL_ERROR "Run manifest does not list trimming_report.json")
endif()
