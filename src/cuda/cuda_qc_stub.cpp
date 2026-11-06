#include "cuda/cuda_qc.h"

#include <stdexcept>
#include <utility>

CudaQcTiming& CudaQcTiming::operator+=(const CudaQcTiming& other)
{
    packing += other.packing;
    hostToDevice += other.hostToDevice;
    kernel += other.kernel;
    deviceToHost += other.deviceToHost;
    total += other.total;
    return *this;
}

struct CudaQcProcessor::Impl {
    CudaQcTiming timing;
};

CudaQcProcessor::CudaQcProcessor()
    : impl_(std::make_unique<Impl>())
{
    throw std::runtime_error(
        "CUDA QC was requested, but this binary was built without CUDA support.");
}

CudaQcProcessor::~CudaQcProcessor() = default;
CudaQcProcessor::CudaQcProcessor(CudaQcProcessor&&) noexcept = default;
CudaQcProcessor& CudaQcProcessor::operator=(CudaQcProcessor&&) noexcept = default;

PerBaseMetrics CudaQcProcessor::process(
    const std::vector<FastqRecord>&,
    const std::vector<std::uint8_t>*)
{
    throw std::runtime_error(
        "CUDA QC was requested, but this binary was built without CUDA support.");
}

const CudaQcTiming& CudaQcProcessor::lastTiming() const noexcept
{
    return impl_->timing;
}
