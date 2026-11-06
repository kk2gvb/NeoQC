#include "compute_backend.h"

#include <stdexcept>

#ifdef NEOQC_ENABLE_CUDA
#include "cuda_backend.h"
#endif

ComputeBackend parseComputeBackend(const std::string& value)
{
    if (value == "cpu") {
        return ComputeBackend::CPU;
    }
    if (value == "cuda") {
        return ComputeBackend::CUDA;
    }
    throw std::runtime_error(
        "Invalid backend '" + value + "'. Expected 'cpu' or 'cuda'.");
}

const char* computeBackendName(ComputeBackend backend)
{
    switch (backend) {
    case ComputeBackend::CPU:
        return "CPU (OpenMP)";
    case ComputeBackend::CUDA:
        return "CUDA";
    }
    return "unknown";
}

BackendSelection selectComputeBackend(ComputeBackend requested)
{
    BackendSelection selection;
    selection.requested = requested;

    if (requested == ComputeBackend::CPU) {
        selection.diagnostic = "CPU backend selected.";
        return selection;
    }

#ifdef NEOQC_ENABLE_CUDA
    selection.diagnostic = runCudaBackendSelfTest();
    selection.effective = ComputeBackend::CUDA;
    return selection;
#else
    throw std::runtime_error(
        "CUDA backend was requested, but this NeoQC binary was built without "
        "CUDA support. Rebuild with -DNEOQC_ENABLE_CUDA=ON.");
#endif
}
