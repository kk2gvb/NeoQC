#pragma once

#include <string>

enum class ComputeBackend {
    CPU,
    CUDA
};

struct BackendSelection {
    ComputeBackend requested = ComputeBackend::CPU;
    ComputeBackend effective = ComputeBackend::CPU;
    std::string diagnostic;
};

ComputeBackend parseComputeBackend(const std::string& value);
const char* computeBackendName(ComputeBackend backend);

// Validates the requested backend before FASTQ processing starts.  CUDA runs a
// Runtime self-test before selected per-base QC metrics are submitted to GPU.
BackendSelection selectComputeBackend(ComputeBackend requested);
