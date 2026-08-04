cmake_minimum_required(VERSION 3.22)
project(neoqc LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Release-сборка по умолчанию
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# Базовые оптимизации (без OpenMP и CUDA — они вернутся, когда появится код, их использующий)
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -march=native")

find_package(ZLIB REQUIRED)

add_executable(neoqc
    src/main.cpp
    src/fastq_reader.cpp
    src/quality_analyzer.cpp
    src/plot_runner.cpp
    src/sample_sheet.cpp
)

target_include_directories(neoqc PRIVATE include)
target_link_libraries(neoqc PRIVATE ZLIB::ZLIB)

message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

enable_testing()

# Папка с тестовыми данными (пути, как вы сказали, верные)
set(TEST_DATA_DIR "${CMAKE_SOURCE_DIR}/tests/data")
set(TEST_OUT_DIR "${CMAKE_BINARY_DIR}/test_output")

# ==============================================================================
# ПОЗИТИВНЫЕ ТЕСТЫ (Программа ДОЛЖНА завершиться успешно, код 0)
# ==============================================================================

# 1. Корректный одиночный FASTQ
add_test(NAME test_correct_single
         COMMAND neoqc --r1 ${TEST_DATA_DIR}/correct_single.fq.gz --sample-id test_single --out ${TEST_OUT_DIR})

add_test(NAME test_per_sequence_quality
         COMMAND ${CMAKE_COMMAND}
                 -DREPORT=${TEST_OUT_DIR}/per_sequence_quality_R1.tsv
                 -P ${CMAKE_SOURCE_DIR}/tests/check_per_sequence_quality.cmake)
set_tests_properties(test_per_sequence_quality PROPERTIES DEPENDS test_correct_single)

# 2. Поиск адаптеров можно отключить для измерения его влияния на скорость.
add_test(NAME test_skip_adapters
         COMMAND neoqc --r1 ${TEST_DATA_DIR}/correct_single.fq.gz --sample-id test_skip_adapters --out ${TEST_OUT_DIR} --skip-adapters)

# Корректная таблица образцов (запуск QC по ней будет добавлен отдельно).
add_test(NAME test_valid_sample_sheet
         COMMAND neoqc --samples ${CMAKE_SOURCE_DIR}/examples/samples.csv)
set_tests_properties(test_valid_sample_sheet PROPERTIES WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})

foreach(sheet missing_patient duplicate_sample missing_file duplicate_file different_references wrong_material missing_normal)
    add_test(NAME test_sample_sheet_${sheet}
             COMMAND neoqc --samples ${CMAKE_SOURCE_DIR}/tests/sample_sheets/${sheet}.csv)
    set_tests_properties(test_sample_sheet_${sheet} PROPERTIES
        WILL_FAIL TRUE
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
endforeach()

add_test(NAME test_sample_sheet_missing_tumor_rna_warning
         COMMAND neoqc --samples ${CMAKE_SOURCE_DIR}/tests/sample_sheets/missing_tumor_rna.csv)
set_tests_properties(test_sample_sheet_missing_tumor_rna_warning PROPERTIES
    PASS_REGULAR_EXPRESSION "Warning: Patient P001: tumor RNA is not specified"
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})

add_test(NAME test_sample_sheet_batch_run
         COMMAND neoqc --samples ${CMAKE_SOURCE_DIR}/tests/sample_sheets/batch_run.csv
                       --out ${TEST_OUT_DIR}/batch_case)
set_tests_properties(test_sample_sheet_batch_run PROPERTIES
    PASS_REGULAR_EXPRESSION "All samples passed"
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})

add_test(NAME test_case_summary
         COMMAND ${CMAKE_COMMAND}
                 -DCASE_SUMMARY=${TEST_OUT_DIR}/batch_case/P001/case_summary.json
                 -P ${CMAKE_SOURCE_DIR}/tests/check_case_summary.cmake)
set_tests_properties(test_case_summary PROPERTIES DEPENDS test_sample_sheet_batch_run)

# 3. Корректная пара R1/R2
add_test(NAME test_correct_paired
         COMMAND neoqc --r1 ${TEST_DATA_DIR}/correct_pair_R1.fq.gz --r2 ${TEST_DATA_DIR}/correct_pair_R2.fq.gz --sample-id test_pair --out ${TEST_OUT_DIR})

add_test(NAME test_adapter_profiles
         COMMAND neoqc --r1 ${TEST_DATA_DIR}/adapter_profile_R1.fq
                       --r2 ${TEST_DATA_DIR}/adapter_profile_R2.fq
                       --sample-id adapter_profile
                       --out ${TEST_OUT_DIR}/adapter_profile)
add_test(NAME test_adapter_profiles_report
         COMMAND ${CMAKE_COMMAND}
                 -DR1_REPORT=${TEST_OUT_DIR}/adapter_profile/adapter_content_R1.tsv
                 -DR2_REPORT=${TEST_OUT_DIR}/adapter_profile/adapter_content_R2.tsv
                 -P ${CMAKE_SOURCE_DIR}/tests/check_adapter_profiles.cmake)
set_tests_properties(test_adapter_profiles_report PROPERTIES DEPENDS test_adapter_profiles)

add_test(NAME test_partial_adapter_profiles
         COMMAND neoqc --r1 ${TEST_DATA_DIR}/partial_adapter_R1.fq
                       --r2 ${TEST_DATA_DIR}/partial_adapter_R2.fq
                       --sample-id partial_adapter
                       --out ${TEST_OUT_DIR}/partial_adapter)
add_test(NAME test_partial_adapter_profiles_report
         COMMAND ${CMAKE_COMMAND}
                 -DR1_REPORT=${TEST_OUT_DIR}/partial_adapter/adapter_content_R1.tsv
                 -DR2_REPORT=${TEST_OUT_DIR}/partial_adapter/adapter_content_R2.tsv
                 -P ${CMAKE_SOURCE_DIR}/tests/check_adapter_profiles.cmake)
set_tests_properties(test_partial_adapter_profiles_report PROPERTIES DEPENDS test_partial_adapter_profiles)

# 4. Строка длиннее 4096 символов (Тут решение за вами: упадет или обработает. По умолчанию ставим WILL_FAIL)
add_test(NAME test_long_line_4096
         COMMAND neoqc --r1 ${TEST_DATA_DIR}/long_line_4096.fq.gz --sample-id test_long_line --out ${TEST_OUT_DIR})

# 5. Файл с окончаниями строк Windows (Тоже на ваше усмотрение: парсер может переварить \r, а может выдать ошибку)
add_test(NAME test_windows_clrf
         COMMAND neoqc --r1 ${TEST_DATA_DIR}/windows_clrf.fq.gz --sample-id test_clrf --out ${TEST_OUT_DIR})

# ==============================================================================
# НЕГАТИВНЫЕ ТЕСТЫ (Программа ДОЛЖНА упасть с ошибкой, WILL_FAIL)
# ==============================================================================

# 6. Неверный заголовок без @
add_test(NAME test_error_no_at
         COMMAND neoqc --r1 ${TEST_DATA_DIR}/error_no_at.fq.gz --sample-id err6 --out ${TEST_OUT_DIR})
set_tests_properties(test_error_no_at PROPERTIES WILL_FAIL TRUE)

# 7. Неверная третья строка без +
add_test(NAME test_error_no_plus
         COMMAND neoqc --r1 ${TEST_DATA_DIR}/error_no_plus.fq.gz --sample-id err7 --out ${TEST_OUT_DIR})
set_tests_properties(test_error_no_plus PROPERTIES WILL_FAIL TRUE)

# 8. Разная длина последовательности и качества
add_test(NAME test_error_diff_len
         COMMAND neoqc --r1 ${TEST_DATA_DIR}/error_diff_len.fq.gz --sample-id err8 --out ${TEST_OUT_DIR})
set_tests_properties(test_error_diff_len PROPERTIES WILL_FAIL TRUE)

# 9. Незаконченная запись
add_test(NAME test_error_incomplete
         COMMAND neoqc --r1 ${TEST_DATA_DIR}/error_incomplete.fq.gz --sample-id err9 --out ${TEST_OUT_DIR})
set_tests_properties(test_error_incomplete PROPERTIES WILL_FAIL TRUE)

# 10. Пустая строка внутри файла
add_test(NAME test_error_empty_line
         COMMAND neoqc --r1 ${TEST_DATA_DIR}/error_empty_line.fq.gz --sample-id err10 --out ${TEST_OUT_DIR})
set_tests_properties(test_error_empty_line PROPERTIES WILL_FAIL TRUE)

# 11. Разное число записей в R1 и R2
add_test(NAME test_diff_counts
         COMMAND neoqc --r1 ${TEST_DATA_DIR}/diff_counts_R1.fq.gz --r2 ${TEST_DATA_DIR}/diff_counts_R2.fq.gz --sample-id err11 --out ${TEST_OUT_DIR})
set_tests_properties(test_diff_counts PROPERTIES WILL_FAIL TRUE)

# 12. Несовпадающие идентификаторы пары
add_test(NAME test_mismatch_ids
         COMMAND neoqc --r1 ${TEST_DATA_DIR}/mismatch_ids_R1.fq.gz --r2 ${TEST_DATA_DIR}/mismatch_ids_R2.fq.gz --sample-id err12 --out ${TEST_OUT_DIR})
set_tests_properties(test_mismatch_ids PROPERTIES WILL_FAIL TRUE)

# 13. Повреждённый gzip (Библиотека ZLIB должна выкинуть ошибку при чтении)
add_test(NAME test_corrupted_gzip
         COMMAND neoqc --r1 ${TEST_DATA_DIR}/corrupted_gzip.fq.gz --sample-id err13 --out ${TEST_OUT_DIR})
set_tests_properties(test_corrupted_gzip PROPERTIES WILL_FAIL TRUE)
