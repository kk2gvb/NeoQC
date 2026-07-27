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

# 2. Корректная пара R1/R2
add_test(NAME test_correct_paired
         COMMAND neoqc --r1 ${TEST_DATA_DIR}/correct_pair_R1.fq.gz --r2 ${TEST_DATA_DIR}/correct_pair_R2.fq.gz --sample-id test_pair --out ${TEST_OUT_DIR})

# 3. Строка длиннее 4096 символов (Тут решение за вами: упадет или обработает. По умолчанию ставим WILL_FAIL)
add_test(NAME test_long_line_4096
         COMMAND neoqc --r1 ${TEST_DATA_DIR}/long_line_4096.fq.gz --sample-id err10 --out ${TEST_OUT_DIR})

# 4. Файл с окончаниями строк Windows (Тоже на ваше усмотрение: парсер может переварить \r, а может выдать ошибку)
add_test(NAME test_windows_clrf
         COMMAND neoqc --r1 ${TEST_DATA_DIR}/windows_clrf.fq.gz --sample-id err11 --out ${TEST_OUT_DIR})

# ==============================================================================
# НЕГАТИВНЫЕ ТЕСТЫ (Программа ДОЛЖНА упасть с ошибкой, WILL_FAIL)
# ==============================================================================

# 5. Неверный заголовок без @
add_test(NAME test_error_no_at
         COMMAND neoqc --r1 ${TEST_DATA_DIR}/error_no_at.fq.gz --sample-id err3 --out ${TEST_OUT_DIR})
set_tests_properties(test_error_no_at PROPERTIES WILL_FAIL TRUE)

# 6. Неверная третья строка без +
add_test(NAME test_error_no_plus
         COMMAND neoqc --r1 ${TEST_DATA_DIR}/error_no_plus.fq.gz --sample-id err4 --out ${TEST_OUT_DIR})
set_tests_properties(test_error_no_plus PROPERTIES WILL_FAIL TRUE)

# 7. Разная длина последовательности и качества
add_test(NAME test_error_diff_len
         COMMAND neoqc --r1 ${TEST_DATA_DIR}/error_diff_len.fq.gz --sample-id err5 --out ${TEST_OUT_DIR})
set_tests_properties(test_error_diff_len PROPERTIES WILL_FAIL TRUE)

# 8. Незаконченная запись
add_test(NAME test_error_incomplete
         COMMAND neoqc --r1 ${TEST_DATA_DIR}/error_incomplete.fq.gz --sample-id err6 --out ${TEST_OUT_DIR})
set_tests_properties(test_error_incomplete PROPERTIES WILL_FAIL TRUE)

# 9. Пустая строка внутри файла
add_test(NAME test_error_empty_line
         COMMAND neoqc --r1 ${TEST_DATA_DIR}/error_empty_line.fq.gz --sample-id err7 --out ${TEST_OUT_DIR})
set_tests_properties(test_error_empty_line PROPERTIES WILL_FAIL TRUE)

# 10. Разное число записей в R1 и R2
add_test(NAME test_diff_counts
         COMMAND neoqc --r1 ${TEST_DATA_DIR}/diff_counts_R1.fq.gz --r2 ${TEST_DATA_DIR}/diff_counts_R2.fq.gz --sample-id err8 --out ${TEST_OUT_DIR})
set_tests_properties(test_diff_counts PROPERTIES WILL_FAIL TRUE)

# 11. Несовпадающие идентификаторы пары
add_test(NAME test_mismatch_ids
         COMMAND neoqc --r1 ${TEST_DATA_DIR}/mismatch_ids_R1.fq.gz --r2 ${TEST_DATA_DIR}/mismatch_ids_R2.fq.gz --sample-id err9 --out ${TEST_OUT_DIR})
set_tests_properties(test_mismatch_ids PROPERTIES WILL_FAIL TRUE)

# 12. Повреждённый gzip (Библиотека ZLIB должна выкинуть ошибку при чтении)
add_test(NAME test_corrupted_gzip
         COMMAND neoqc --r1 ${TEST_DATA_DIR}/corrupted_gzip.fq.gz --sample-id err12 --out ${TEST_OUT_DIR})
set_tests_properties(test_corrupted_gzip PROPERTIES WILL_FAIL TRUE)
