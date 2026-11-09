if(NOT DEFINED NEOQC)
    message(FATAL_ERROR "NEOQC is not defined")
endif()

if(NOT DEFINED TEST_ROOT)
    message(FATAL_ERROR "TEST_ROOT is not defined")
endif()

file(MAKE_DIRECTORY "${TEST_ROOT}")
set(HEADER "patient_id,sample_id,sample_role,material,r1,r2,platform,library_type,reference\n")

set(R1 "${TEST_ROOT}/r1.fastq")
set(R2 "${TEST_ROOT}/r2.fastq")
file(WRITE "${R1}" "@r1\nA\n+\nI\n")
file(WRITE "${R2}" "@r2\nA\n+\nI\n")

set(VALID_SHEET "${TEST_ROOT}/valid.csv")
file(WRITE "${VALID_SHEET}"
    "${HEADER}"
    "Patient-01.2,Sample_1.3,tumor,dna,${R1},,illumina,wes,GRCh38\n"
    "Patient-01.2,Normal-2,normal,dna,${R2},,illumina,wes,GRCh38\n")
execute_process(
    COMMAND "${NEOQC}" --samples "${VALID_SHEET}"
    RESULT_VARIABLE VALID_RESULT
    OUTPUT_VARIABLE VALID_OUTPUT
    ERROR_VARIABLE VALID_ERROR)
if(NOT VALID_RESULT EQUAL 0)
    message(FATAL_ERROR
        "Safe path components were rejected.\nstdout:\n${VALID_OUTPUT}\nstderr:\n${VALID_ERROR}")
endif()

function(expect_invalid_component PATIENT_ID SAMPLE_ID LABEL)
    set(SHEET "${TEST_ROOT}/${LABEL}.csv")
    set(OUT_DIR "${TEST_ROOT}/${LABEL}-out")
    set(ESCAPE_DIR "${TEST_ROOT}/${LABEL}-escape")
    set(ESCAPE_IN_OUTPUT_DIR "${OUT_DIR}/${LABEL}-escape")
    file(WRITE "${SHEET}" "${HEADER}"
        "${PATIENT_ID},${SAMPLE_ID},tumor,dna,${R1},,illumina,wes,GRCh38\n")

    execute_process(
        COMMAND "${NEOQC}" --samples "${SHEET}" --out "${OUT_DIR}"
        RESULT_VARIABLE RESULT
        OUTPUT_VARIABLE OUTPUT
        ERROR_VARIABLE ERROR)
    if(RESULT EQUAL 0)
        message(FATAL_ERROR "Unsafe path component '${LABEL}' was accepted")
    endif()
    if(EXISTS "${OUT_DIR}" OR EXISTS "${ESCAPE_DIR}" OR EXISTS "${ESCAPE_IN_OUTPUT_DIR}")
        message(FATAL_ERROR
            "Unsafe path component '${LABEL}' created an output directory")
    endif()
endfunction()

expect_invalid_component("." "sample" "dot-patient")
expect_invalid_component(".." "sample" "dotdot-patient")
expect_invalid_component("../parent-patient-escape" "sample" "parent-patient")
expect_invalid_component("patient/sample" "sample" "slash-patient")
expect_invalid_component("patient\\sample" "sample" "backslash-patient")
expect_invalid_component("/absolute" "sample" "absolute-patient")
expect_invalid_component("patient" "." "dot-sample")
expect_invalid_component("patient" ".." "dotdot-sample")
expect_invalid_component("patient" "../parent-sample-escape" "parent-sample")
expect_invalid_component("patient" "sample/name" "slash-sample")
expect_invalid_component("patient" "sample\\name" "backslash-sample")

message(STATUS "Sample-sheet path component validation passed")
