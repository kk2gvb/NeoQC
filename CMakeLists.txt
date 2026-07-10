cmake_minimum_required(VERSION 3.22)
project(mRNA_QC LANGUAGES CXX CUDA)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CUDA_ARCHITECTURES 86)   # RTX 3070 = sm_86

# Оптимизации
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -march=native -fopenmp")
set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -O3 -Xcompiler -fopenmp")

find_package(ZLIB REQUIRED)

add_executable(qc
    src/main.cpp
    src/fastq_reader.cpp
    src/quality_analyzer.cpp
    #src/cuda/kernels.cu
)

target_include_directories(qc PRIVATE include)
target_link_libraries(qc PRIVATE ZLIB::ZLIB)

# CUDA настройки
set_target_properties(qc PROPERTIES 
    CUDA_SEPARABLE_COMPILATION ON
    POSITION_INDEPENDENT_CODE ON
)

message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "CUDA Architectures: ${CMAKE_CUDA_ARCHITECTURES}")
