#include "cuda_backend.h"

#include <cuda_runtime.h>

#include <sstream>
#include <stdexcept>
#include <string>

namespace {

std::string cudaErrorMessage(const char* operation, cudaError_t error)
{
    std::ostringstream message;
    message << operation << " failed: " << cudaGetErrorString(error)
            << " (CUDA error " << static_cast<int>(error) << ')';
    return message.str();
}

void checkCuda(cudaError_t error, const char* operation)
{
    if (error != cudaSuccess) {
        throw std::runtime_error(cudaErrorMessage(operation, error));
    }
}

__global__ void cudaSelfTestKernel(int* value)
{
    if (blockIdx.x == 0 && threadIdx.x == 0) {
        *value = 0x4e5143; // "NQC"
    }
}

} // namespace

std::string runCudaBackendSelfTest()
{
    int deviceCount = 0;
    const cudaError_t countResult = cudaGetDeviceCount(&deviceCount);
    if (countResult != cudaSuccess) {
        throw std::runtime_error(cudaErrorMessage("cudaGetDeviceCount", countResult));
    }
    if (deviceCount == 0) {
        throw std::runtime_error("No CUDA-capable device is available.");
    }

    int device = 0;
    checkCuda(cudaGetDevice(&device), "cudaGetDevice");

    cudaDeviceProp properties{};
    checkCuda(cudaGetDeviceProperties(&properties, device), "cudaGetDeviceProperties");
    if (properties.major == 0) {
        throw std::runtime_error("Selected CUDA device reports an invalid compute capability.");
    }

    int* deviceValue = nullptr;
    try {
        checkCuda(cudaMalloc(&deviceValue, sizeof(*deviceValue)), "cudaMalloc self-test buffer");
        cudaSelfTestKernel<<<1, 1>>>(deviceValue);
        checkCuda(cudaGetLastError(), "CUDA self-test kernel launch");
        checkCuda(cudaDeviceSynchronize(), "CUDA self-test kernel execution");

        int hostValue = 0;
        checkCuda(cudaMemcpy(&hostValue, deviceValue, sizeof(hostValue),
                             cudaMemcpyDeviceToHost),
                  "cudaMemcpy self-test result");
        if (hostValue != 0x4e5143) {
            throw std::runtime_error("CUDA self-test kernel returned an unexpected result.");
        }
    }
    catch (...) {
        if (deviceValue != nullptr) {
            cudaFree(deviceValue);
        }
        throw;
    }
    checkCuda(cudaFree(deviceValue), "cudaFree self-test buffer");

    std::ostringstream message;
    message << "CUDA diagnostic passed: device " << device << " ('"
            << properties.name << "'), compute capability " << properties.major
            << '.' << properties.minor << ". Test kernel completed successfully.";
    return message.str();
}
