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

# Тесты (подключим на шаге 5)
enable_testing()
# add_subdirectory(tests)

message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")