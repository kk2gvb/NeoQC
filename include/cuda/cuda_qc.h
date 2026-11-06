#pragma once

#include "per_base_metrics.h"

#include <cstdint>
#include <memory>
#include <vector>

struct CudaQcTiming {
    double packing = 0.0;
    double hostToDevice = 0.0;
    double kernel = 0.0;
    double deviceToHost = 0.0;
    double total = 0.0;

    CudaQcTiming& operator+=(const CudaQcTiming& other);
};

// The public interface intentionally has no CUDA Runtime types.  It is safe
// to include from the CPU-only build; the implementation is either a CUDA
// PImpl or a CPU-only diagnostic stub selected by CMake.
class CudaQcProcessor {
public:
    CudaQcProcessor();
    ~CudaQcProcessor();

    CudaQcProcessor(const CudaQcProcessor&) = delete;
    CudaQcProcessor& operator=(const CudaQcProcessor&) = delete;
    CudaQcProcessor(CudaQcProcessor&&) noexcept;
    CudaQcProcessor& operator=(CudaQcProcessor&&) noexcept;

    // Computes exact per-base content and quality accumulators for one
    // complete batch. A null mask selects every record.
    PerBaseMetrics process(
        const std::vector<FastqRecord>& records,
        const std::vector<std::uint8_t>* selected = nullptr);

    const CudaQcTiming& lastTiming() const noexcept;

private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};
