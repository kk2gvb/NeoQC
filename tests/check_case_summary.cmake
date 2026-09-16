if(NOT DEFINED CASE_SUMMARY)
    message(FATAL_ERROR "CASE_SUMMARY is not defined")
endif()

if(NOT DEFINED NEOQC)
    message(FATAL_ERROR "NEOQC is not defined")
endif()

if(NOT DEFINED SAMPLE_SHEET)
    message(FATAL_ERROR "SAMPLE_SHEET is not defined")
endif()

if(NOT DEFINED OUT_DIR)
    message(FATAL_ERROR "OUT_DIR is not defined")
endif()

if(NOT DEFINED EXPECTED_RULESET_SHA256)
    message(FATAL_ERROR "EXPECTED_RULESET_SHA256 is not defined")
endif()

if(NOT DEFINED VALID_R1)
    message(FATAL_ERROR "VALID_R1 is not defined")
endif()

if(NOT DEFINED INVALID_R1)
    message(FATAL_ERROR "INVALID_R1 is not defined")
endif()

# ------------------------------------------------------------
# Helper
# ------------------------------------------------------------

function(extract_json_string CONTENT KEY OUTPUT_VAR)
    string(
        REGEX MATCH
        "\"${KEY}\"[ \t\r\n]*:[ \t\r\n]*\"([^\"]+)\""
        MATCH
        "${CONTENT}"
    )

    if(NOT MATCH)
        message(FATAL_ERROR
            "Missing JSON field: ${KEY}"
        )
    endif()

    set(${OUTPUT_VAR}
        "${CMAKE_MATCH_1}"
        PARENT_SCOPE
    )
endfunction()

# ------------------------------------------------------------
# 1. Initial successful batch
# ------------------------------------------------------------

execute_process(
    COMMAND
        "${NEOQC}"
        --samples "${SAMPLE_SHEET}"
        --out "${OUT_DIR}"
    RESULT_VARIABLE INITIAL_RESULT
    OUTPUT_VARIABLE INITIAL_OUTPUT
    ERROR_VARIABLE INITIAL_ERROR
)

if(NOT INITIAL_RESULT EQUAL 0)
    message(FATAL_ERROR
        "Initial batch run failed.\n"
        "stdout:\n${INITIAL_OUTPUT}\n"
        "stderr:\n${INITIAL_ERROR}"
    )
endif()

if(NOT EXISTS "${CASE_SUMMARY}")
    message(FATAL_ERROR
        "case_summary.json does not exist after initial batch: "
        "${CASE_SUMMARY}"
    )
endif()

file(READ "${CASE_SUMMARY}" INITIAL_SUMMARY)

extract_json_string(
    "${INITIAL_SUMMARY}"
    "run_id"
    FIRST_RUN_ID
)

extract_json_string(
    "${INITIAL_SUMMARY}"
    "status"
    FIRST_STATUS
)

extract_json_string(
    "${INITIAL_SUMMARY}"
    "sha256"
    FIRST_SHA256
)

if(NOT FIRST_STATUS STREQUAL "passed")
    message(FATAL_ERROR
        "Initial case summary status is not passed: ${FIRST_STATUS}"
    )
endif()

if(NOT FIRST_SHA256 STREQUAL "${EXPECTED_RULESET_SHA256}")
    message(FATAL_ERROR
        "Ruleset SHA-256 mismatch.\n"
        "Expected: ${EXPECTED_RULESET_SHA256}\n"
        "Actual:   ${FIRST_SHA256}"
    )
endif()

string(LENGTH "${FIRST_SHA256}" SHA_LENGTH)

if(NOT SHA_LENGTH EQUAL 64)
    message(FATAL_ERROR
        "Ruleset SHA-256 must contain 64 hexadecimal characters, "
        "got ${SHA_LENGTH}"
    )
endif()

message(STATUS
    "Initial case summary: "
    "run_id=${FIRST_RUN_ID}, "
    "status=${FIRST_STATUS}, "
    "sha256=${FIRST_SHA256}"
)

# ------------------------------------------------------------
# 2. Repeat successful batch with changed parameters
# ------------------------------------------------------------

execute_process(
    COMMAND
        "${NEOQC}"
        --samples "${SAMPLE_SHEET}"
        --out "${OUT_DIR}"
        --skip-adapters
    RESULT_VARIABLE REPEAT_RESULT
    OUTPUT_VARIABLE REPEAT_OUTPUT
    ERROR_VARIABLE REPEAT_ERROR
)

if(NOT REPEAT_RESULT EQUAL 0)
    message(FATAL_ERROR
        "Repeated batch run failed unexpectedly.\n"
        "stdout:\n${REPEAT_OUTPUT}\n"
        "stderr:\n${REPEAT_ERROR}"
    )
endif()

if(NOT EXISTS "${CASE_SUMMARY}")
    message(FATAL_ERROR
        "case_summary.json disappeared after repeated batch"
    )
endif()

file(READ "${CASE_SUMMARY}" REPEAT_SUMMARY)

extract_json_string(
    "${REPEAT_SUMMARY}"
    "run_id"
    SECOND_RUN_ID
)

extract_json_string(
    "${REPEAT_SUMMARY}"
    "status"
    SECOND_STATUS
)

extract_json_string(
    "${REPEAT_SUMMARY}"
    "sha256"
    SECOND_SHA256
)

if(NOT SECOND_STATUS STREQUAL "passed")
    message(FATAL_ERROR
        "Repeated case summary status is not passed: ${SECOND_STATUS}"
    )
endif()

if(SECOND_RUN_ID STREQUAL FIRST_RUN_ID)
    message(FATAL_ERROR
        "case_summary.json was not regenerated: run_id did not change"
    )
endif()

if(NOT SECOND_SHA256 STREQUAL "${EXPECTED_RULESET_SHA256}")
    message(FATAL_ERROR
        "Ruleset SHA-256 changed between runs.\n"
        "Expected: ${EXPECTED_RULESET_SHA256}\n"
        "Actual:   ${SECOND_SHA256}"
    )
endif()

if(NOT REPEAT_SUMMARY MATCHES
        "\"skip_adapters\"[ \t\r\n]*:[ \t\r\n]*true")
    message(FATAL_ERROR
        "Repeated case summary does not contain skip_adapters=true"
    )
endif()

message(STATUS
    "Repeated batch regenerated summary successfully: "
    "${FIRST_RUN_ID} -> ${SECOND_RUN_ID}"
)

# ------------------------------------------------------------
# 3. Failed batch must regenerate case_summary.json
# ------------------------------------------------------------

file(READ "${SAMPLE_SHEET}" BATCH_CONTENT)

string(REPLACE
    "${VALID_R1}"
    "${INVALID_R1}"
    FAILED_BATCH_CONTENT
    "${BATCH_CONTENT}"
)

set(FAILED_SHEET
    "${OUT_DIR}/.case_summary_failure_test.csv"
)

file(WRITE
    "${FAILED_SHEET}"
    "${FAILED_BATCH_CONTENT}"
)

execute_process(
    COMMAND
        "${NEOQC}"
        --samples "${FAILED_SHEET}"
        --out "${OUT_DIR}"
    RESULT_VARIABLE FAILED_RESULT
    OUTPUT_VARIABLE FAILED_OUTPUT
    ERROR_VARIABLE FAILED_ERROR
)

file(REMOVE "${FAILED_SHEET}")

if(FAILED_RESULT EQUAL 0)
    message(FATAL_ERROR
        "Failure regression run unexpectedly succeeded"
    )
endif()

if(NOT EXISTS "${CASE_SUMMARY}")
    message(FATAL_ERROR
        "case_summary.json does not exist after failed batch"
    )
endif()

file(READ "${CASE_SUMMARY}" FAILED_SUMMARY)

extract_json_string(
    "${FAILED_SUMMARY}"
    "run_id"
    FAILED_RUN_ID
)

extract_json_string(
    "${FAILED_SUMMARY}"
    "status"
    FAILED_STATUS
)

if(NOT FAILED_STATUS STREQUAL "failed")
    message(FATAL_ERROR
        "Failed batch did not regenerate a failed case summary. "
        "Actual status: ${FAILED_STATUS}"
    )
endif()

if(FAILED_RUN_ID STREQUAL SECOND_RUN_ID)
    message(FATAL_ERROR
        "Failed batch left the previous case summary unchanged"
    )
endif()

message(STATUS
    "Failed batch regenerated case summary successfully: "
    "${SECOND_RUN_ID} -> ${FAILED_RUN_ID}"
)

message(STATUS
    "P0.2 case_summary regression test passed"
)
