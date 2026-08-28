cmake_minimum_required(VERSION 3.22)
project(neoqc LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Release-сборка по умолчанию
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# Базовые оптимизации
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -march=native")

find_package(ZLIB REQUIRED)
find_package(OpenMP REQUIRED)

add_executable(neoqc
    src/main.cpp
    src/fastq_reader.cpp
    src/quality_analyzer.cpp
    src/plot_runner.cpp
    src/sample_sheet.cpp
    src/gc_model.cpp
)

target_include_directories(neoqc PRIVATE include)
target_link_libraries(neoqc PRIVATE
    ZLIB::ZLIB
    OpenMP::OpenMP_CXX
)
target_compile_definitions(neoqc PRIVATE
    NEOQC_PLOT_SCRIPT="${CMAKE_CURRENT_SOURCE_DIR}/scripts/plot_results.py"
)

add_executable(neoqc_duplication_test
    tests/test_quality_analyzer_duplication.cpp
    src/quality_analyzer.cpp
    src/gc_model.cpp
)
target_include_directories(neoqc_duplication_test PRIVATE include)

message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

enable_testing()

add_test(NAME test_unlimited_sequence_duplication
         COMMAND neoqc_duplication_test)

find_package(Python3 COMPONENTS Interpreter)
if(Python3_Interpreter_FOUND)
    add_test(NAME test_qc_report
             COMMAND ${Python3_EXECUTABLE} ${CMAKE_SOURCE_DIR}/tests/test_qc_report.py)
    set_tests_properties(test_qc_report PROPERTIES
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})

    add_test(NAME test_qc_rules
             COMMAND ${Python3_EXECUTABLE} ${CMAKE_SOURCE_DIR}/tests/test_qc_rules.py)
    set_tests_properties(test_qc_rules PROPERTIES
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})

    execute_process(
        COMMAND ${Python3_EXECUTABLE} -c "import matplotlib"
        RESULT_VARIABLE MATPLOTLIB_IMPORT_STATUS
        OUTPUT_QUIET
        ERROR_QUIET
    )
    if(MATPLOTLIB_IMPORT_STATUS EQUAL 0)
        add_test(NAME test_plot_results
                 COMMAND ${Python3_EXECUTABLE} ${CMAKE_SOURCE_DIR}/tests/test_plot_results.py)
        set_tests_properties(test_plot_results PROPERTIES
            WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
    else()
        message(STATUS "Matplotlib not found: plot rendering test is disabled")
    endif()
endif()

# Папка с тестовыми данными
set(TEST_DATA_DIR "${CMAKE_SOURCE_DIR}/tests/data")
set(TEST_OUT_DIR "${CMAKE_BINARY_DIR}/test_output")
set(TEST_SAMPLE_SHEETS_DIR "${CMAKE_BINARY_DIR}/test_sample_sheets")

file(MAKE_DIRECTORY "${TEST_SAMPLE_SHEETS_DIR}")

configure_file(
    "${CMAKE_SOURCE_DIR}/tests/sample_sheets/batch_run.csv.in"
    "${TEST_SAMPLE_SHEETS_DIR}/batch_run.csv"
    @ONLY
)

# Генерация всех FASTQ-данных для тестов.
add_custom_target(generate_test_data ALL
    COMMAND bash
            ${CMAKE_SOURCE_DIR}/scripts/generate_test_data.sh
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    DEPENDS
        ${CMAKE_SOURCE_DIR}/scripts/generate_test_data.sh
        ${CMAKE_SOURCE_DIR}/scripts/generate_fastq.py
        ${CMAKE_SOURCE_DIR}/scripts/generate_adapter_tests.sh
    COMMENT "Generating NeoQC test data"
    VERBATIM
)

# ==============================================================================
# ПОЗИТИВНЫЕ ТЕСТЫ (Программа ДОЛЖНА завершиться успешно, код 0)
# ==============================================================================

# Корректный одиночный FASTQ
add_test(NAME test_correct_single
         COMMAND neoqc --r1 ${TEST_DATA_DIR}/correct_single.fq.gz --sample-id test_single --out ${TEST_OUT_DIR})

add_test(NAME test_per_sequence_quality
         COMMAND ${CMAKE_COMMAND}
                 -DREPORT=${TEST_OUT_DIR}/per_sequence_quality_R1.tsv
                 -P ${CMAKE_SOURCE_DIR}/tests/check_per_sequence_quality.cmake)
set_tests_properties(test_per_sequence_quality PROPERTIES DEPENDS test_correct_single)

add_test(NAME test_per_cycle_quality
         COMMAND ${CMAKE_COMMAND}
                 -DREPORT=${TEST_OUT_DIR}/per_cycle_R1.tsv
                 -P ${CMAKE_SOURCE_DIR}/tests/check_per_cycle_quality.cmake)
set_tests_properties(test_per_cycle_quality PROPERTIES DEPENDS test_correct_single)

# Проверка качества последовательностей реализована как зависимая проверка

# Пропуск поиска адаптеров для измерения влияния на скорость
add_test(NAME test_skip_adapters
         COMMAND neoqc --r1 ${TEST_DATA_DIR}/correct_single.fq.gz --sample-id test_skip_adapters --out ${TEST_OUT_DIR} --skip-adapters)

# Корректная таблица образцов (запуск QC по ней будет добавлен отдельно)
add_test(NAME test_valid_sample_sheet
         COMMAND neoqc --samples ${CMAKE_SOURCE_DIR}/examples/samples.csv)
set_tests_properties(test_valid_sample_sheet PROPERTIES WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})

# Негативные тесты таблиц образцов
foreach(sheet missing_patient duplicate_sample missing_file duplicate_file different_references wrong_material missing_normal)
    add_test(NAME test_sample_sheet_${sheet}
             COMMAND neoqc --samples ${CMAKE_SOURCE_DIR}/tests/sample_sheets/${sheet}.csv)
    set_tests_properties(test_sample_sheet_${sheet} PROPERTIES
        WILL_FAIL TRUE
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
endforeach()

# Предупреждение об отсутствии опухолевой РНК
add_test(NAME test_sample_sheet_missing_tumor_rna_warning
         COMMAND neoqc --samples ${CMAKE_SOURCE_DIR}/tests/sample_sheets/missing_tumor_rna.csv)
set_tests_properties(test_sample_sheet_missing_tumor_rna_warning PROPERTIES
    PASS_REGULAR_EXPRESSION "Warning: Patient P001: tumor RNA is not specified"
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})

# Пакетный запуск с таблицей образцов
add_test(NAME test_sample_sheet_batch_run
         COMMAND neoqc --samples ${TEST_SAMPLE_SHEETS_DIR}/batch_run.csv
                       --out ${TEST_OUT_DIR}/batch_case)
set_tests_properties(test_sample_sheet_batch_run PROPERTIES
    PASS_REGULAR_EXPRESSION "All samples passed"
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})

# Проверка итогового отчёта пакетного запуска
add_test(NAME test_case_summary
         COMMAND ${CMAKE_COMMAND}
                 -DCASE_SUMMARY=${TEST_OUT_DIR}/batch_case/P001/case_summary.json
                 -P ${CMAKE_SOURCE_DIR}/tests/check_case_summary.cmake)
set_tests_properties(test_case_summary PROPERTIES DEPENDS test_sample_sheet_batch_run)

# Корректная пара R1/R2
add_test(NAME test_correct_paired
         COMMAND neoqc
                 --r1 ${TEST_DATA_DIR}/correct_pair_R1.fq.gz
                 --r2 ${TEST_DATA_DIR}/correct_pair_R2.fq.gz
                 --sample-id test_pair
                 --out ${TEST_OUT_DIR})

# Контракт точных метрик дупликации 50-nt префиксов.
add_test(NAME test_sequence_duplication
         COMMAND neoqc
                 --r1 ${CMAKE_SOURCE_DIR}/tests/fixtures/sequence_duplication.fq
                 --sample-id duplication_contract
                 --out ${TEST_OUT_DIR}/sequence_duplication)

add_test(NAME test_numeric_regression
         COMMAND neoqc
                 --r1 ${CMAKE_SOURCE_DIR}/tests/fixtures/numeric_regression.fq
                 --sample-id numeric_regression
                 --out ${TEST_OUT_DIR}/numeric_regression)

set_tests_properties(test_numeric_regression PROPERTIES
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})

add_test(NAME test_numeric_regression_contract
         COMMAND ${CMAKE_COMMAND}
                 -DOUTPUT_DIR=${TEST_OUT_DIR}/numeric_regression
                 -P ${CMAKE_SOURCE_DIR}/tests/check_numeric_regression.cmake)

set_tests_properties(test_numeric_regression_contract PROPERTIES
    DEPENDS test_numeric_regression)

add_test(NAME test_sequence_duplication_artifacts
         COMMAND ${CMAKE_COMMAND}
                 -DOUTPUT_DIR=${TEST_OUT_DIR}/sequence_duplication
                 -P ${CMAKE_SOURCE_DIR}/tests/check_sequence_duplication.cmake)
set_tests_properties(test_sequence_duplication_artifacts PROPERTIES
    DEPENDS test_sequence_duplication)

# Профили адаптеров в полной выборке
add_test(NAME test_adapter_profiles
         COMMAND neoqc --r1 ${TEST_DATA_DIR}/adapter_profile_R1.fq
                       --r2 ${TEST_DATA_DIR}/adapter_profile_R2.fq
                       --sample-id adapter_profile
                       --out ${TEST_OUT_DIR}/adapter_profile)

# Проверка отчёта профилей адаптеров
add_test(NAME test_adapter_profiles_report
         COMMAND ${CMAKE_COMMAND}
                 -DR1_REPORT=${TEST_OUT_DIR}/adapter_profile/adapter_content_R1.tsv
                 -DR2_REPORT=${TEST_OUT_DIR}/adapter_profile/adapter_content_R2.tsv
                 -DEXPECTED_VALUE=100\\.0000
                 -P ${CMAKE_SOURCE_DIR}/tests/check_adapter_profiles.cmake)
set_tests_properties(test_adapter_profiles_report PROPERTIES DEPENDS test_adapter_profiles)

# Частичные профили адаптеров
add_test(NAME test_partial_adapter_profiles
         COMMAND neoqc --r1 ${TEST_DATA_DIR}/partial_adapter_R1.fq
                       --r2 ${TEST_DATA_DIR}/partial_adapter_R2.fq
                       --sample-id partial_adapter
                       --out ${TEST_OUT_DIR}/partial_adapter)

# Проверка отчёта частичных профилей адаптеров
add_test(NAME test_partial_adapter_profiles_report
         COMMAND ${CMAKE_COMMAND}
                 -DR1_REPORT=${TEST_OUT_DIR}/partial_adapter/adapter_content_R1.tsv
                 -DR2_REPORT=${TEST_OUT_DIR}/partial_adapter/adapter_content_R2.tsv
                 -DEXPECTED_VALUE=50\\.0000
                 -P ${CMAKE_SOURCE_DIR}/tests/check_adapter_profiles.cmake)
set_tests_properties(test_partial_adapter_profiles_report PROPERTIES DEPENDS test_partial_adapter_profiles)

# Строка длиннее 4096 символов
add_test(NAME test_long_line_4096
         COMMAND neoqc --r1 ${TEST_DATA_DIR}/long_line_4096.fq.gz --sample-id test_long_line --out ${TEST_OUT_DIR})

# Файл с окончаниями строк Windows CRLF
add_test(NAME test_windows_clrf
         COMMAND neoqc --r1 ${TEST_DATA_DIR}/windows_clrf.fq.gz --sample-id test_clrf --out ${TEST_OUT_DIR})

# Граничные допустимые символы Phred+33: ASCII 33 (!) и ASCII 126 (~)
add_test(NAME test_quality_phred33_boundaries
         COMMAND neoqc
                 --r1 ${TEST_DATA_DIR}/quality_phred33_boundaries.fq.gz
                 --sample-id quality_boundaries
                 --out ${TEST_OUT_DIR}/quality_phred33_boundaries)

# ==============================================================================
# НЕГАТИВНЫЕ ТЕСТЫ (Программа ДОЛЖНА упасть с ошибкой, WILL_FAIL)
# ==============================================================================

# Недопустимый символ качества ASCII 32 (space)
add_test(NAME test_invalid_quality_ascii_32
         COMMAND neoqc
                 --r1 ${TEST_DATA_DIR}/invalid_quality_ascii_32.fq.gz
                 --sample-id invalid_quality_32
                 --out ${TEST_OUT_DIR}/invalid_quality_ascii_32)

set_tests_properties(test_invalid_quality_ascii_32 PROPERTIES
    WILL_FAIL TRUE)


# Недопустимый символ качества ASCII 127 (DEL)
add_test(NAME test_invalid_quality_ascii_127
         COMMAND neoqc
                 --r1 ${TEST_DATA_DIR}/invalid_quality_ascii_127.fq.gz
                 --sample-id invalid_quality_127
                 --out ${TEST_OUT_DIR}/invalid_quality_ascii_127)

set_tests_properties(test_invalid_quality_ascii_127 PROPERTIES
    WILL_FAIL TRUE)

# Неверный заголовок без символа @
add_test(NAME test_error_no_at
         COMMAND neoqc --r1 ${TEST_DATA_DIR}/error_no_at.fq.gz --sample-id err6 --out ${TEST_OUT_DIR})
set_tests_properties(test_error_no_at PROPERTIES WILL_FAIL TRUE)

# Неверная третья строка без символа +
add_test(NAME test_error_no_plus
         COMMAND neoqc --r1 ${TEST_DATA_DIR}/error_no_plus.fq.gz --sample-id err7 --out ${TEST_OUT_DIR})
set_tests_properties(test_error_no_plus PROPERTIES WILL_FAIL TRUE)

# Разная длина последовательности и качества
add_test(NAME test_error_diff_len
         COMMAND neoqc --r1 ${TEST_DATA_DIR}/error_diff_len.fq.gz --sample-id err8 --out ${TEST_OUT_DIR})
set_tests_properties(test_error_diff_len PROPERTIES WILL_FAIL TRUE)

# Незаконченная запись в файле
add_test(NAME test_error_incomplete
         COMMAND neoqc --r1 ${TEST_DATA_DIR}/error_incomplete.fq.gz --sample-id err9 --out ${TEST_OUT_DIR})
set_tests_properties(test_error_incomplete PROPERTIES WILL_FAIL TRUE)

# Пустая строка внутри файла
add_test(NAME test_error_empty_line
         COMMAND neoqc --r1 ${TEST_DATA_DIR}/error_empty_line.fq.gz --sample-id err10 --out ${TEST_OUT_DIR})
set_tests_properties(test_error_empty_line PROPERTIES WILL_FAIL TRUE)

# Разное число записей в R1 и R2
add_test(NAME test_diff_counts
         COMMAND neoqc --r1 ${TEST_DATA_DIR}/diff_counts_R1.fq.gz --r2 ${TEST_DATA_DIR}/diff_counts_R2.fq.gz --sample-id err11 --out ${TEST_OUT_DIR})
set_tests_properties(test_diff_counts PROPERTIES WILL_FAIL TRUE)

# Несовпадающие идентификаторы пары
add_test(NAME test_mismatch_ids
         COMMAND neoqc --r1 ${TEST_DATA_DIR}/mismatch_ids_R1.fq.gz --r2 ${TEST_DATA_DIR}/mismatch_ids_R2.fq.gz --sample-id err12 --out ${TEST_OUT_DIR})
set_tests_properties(test_mismatch_ids PROPERTIES WILL_FAIL TRUE)

# Повреждённый gzip (ошибка CRC)
add_test(NAME test_corrupted_crc
         COMMAND neoqc
                 --r1 ${TEST_DATA_DIR}/corrupted_crc.fq.gz
                 --sample-id crc
                 --out ${TEST_OUT_DIR})
set_tests_properties(test_corrupted_crc PROPERTIES WILL_FAIL TRUE)

# Пустой FASTQ файл
add_test(NAME test_empty_fastq
         COMMAND neoqc
                 --r1 ${TEST_DATA_DIR}/empty.fq.gz
                 --sample-id empty
                 --out ${TEST_OUT_DIR})
set_tests_properties(test_empty_fastq PROPERTIES WILL_FAIL TRUE)

# Пустая последовательность в записи
add_test(NAME test_empty_sequence
         COMMAND neoqc
                 --r1 ${TEST_DATA_DIR}/empty_sequence.fq.gz
                 --sample-id empty_seq
                 --out ${TEST_OUT_DIR})
set_tests_properties(test_empty_sequence PROPERTIES WILL_FAIL TRUE)

# Пустое качество в записи
add_test(NAME test_empty_quality
         COMMAND neoqc
                 --r1 ${TEST_DATA_DIR}/empty_quality.fq.gz
                 --sample-id empty_qual
                 --out ${TEST_OUT_DIR})
set_tests_properties(test_empty_quality PROPERTIES WILL_FAIL TRUE)

# Пустая строка между записями
add_test(NAME test_blank_line_between_records
         COMMAND neoqc
                 --r1 ${TEST_DATA_DIR}/blank_line_between_records.fq.gz
                 --sample-id blank_line
                 --out ${TEST_OUT_DIR})
set_tests_properties(test_blank_line_between_records PROPERTIES WILL_FAIL TRUE)

# Несовпадающие пары R1 и R2 (разные идентификаторы)
add_test(NAME test_pair_mismatch
         COMMAND neoqc
                 --r1 ${TEST_DATA_DIR}/pair_mismatch_R1.fq.gz
                 --r2 ${TEST_DATA_DIR}/pair_mismatch_R2.fq.gz
                 --sample-id mismatch
                 --out ${TEST_OUT_DIR})
set_tests_properties(test_pair_mismatch PROPERTIES WILL_FAIL TRUE)

# Некорректные символы в последовательности
add_test(NAME test_invalid_characters
         COMMAND neoqc --r1 ${TEST_DATA_DIR}/invalid_base.fq --sample-id invalid_chars --out ${TEST_OUT_DIR})
set_tests_properties(test_invalid_characters PROPERTIES WILL_FAIL TRUE)
