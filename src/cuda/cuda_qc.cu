#include "cuda/cuda_qc.h"

#include <cuda_runtime.h>

#include <chrono>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <utility>

namespace {

using Clock = std::chrono::steady_clock;

std::string errorMessage(const char* operation, cudaError_t error)
{
    std::ostringstream message;
    message << operation << " failed: " << cudaGetErrorString(error)
            << " (CUDA error " << static_cast<int>(error) << ')';
    return message.str();
}

void checkCuda(cudaError_t error, const char* operation)
{
    if (error != cudaSuccess) {
        throw std::runtime_error(errorMessage(operation, error));
    }
}

std::size_t checkedAdd(std::size_t left, std::size_t right, const char* what)
{
    if (right > std::numeric_limits<std::size_t>::max() - left) {
        throw std::runtime_error(std::string("CUDA batch size overflows while calculating ") + what);
    }
    return left + right;
}

std::size_t checkedMultiply(std::size_t left, std::size_t right, const char* what)
{
    if (left != 0 && right > std::numeric_limits<std::size_t>::max() / left) {
        throw std::runtime_error(std::string("CUDA batch size overflows while calculating ") + what);
    }
    return left * right;
}

template <typename T>
class DeviceBuffer {
public:
    DeviceBuffer() = default;
    ~DeviceBuffer() { reset(); }

    DeviceBuffer(const DeviceBuffer&) = delete;
    DeviceBuffer& operator=(const DeviceBuffer&) = delete;

    T* data() noexcept { return data_; }

    void ensure(std::size_t count, const char* name)
    {
        if (count <= capacity_) {
            return;
        }
        reset();
        checkCuda(cudaMalloc(&data_, checkedMultiply(count, sizeof(T), name)), name);
        capacity_ = count;
    }

    void reset() noexcept
    {
        if (data_ != nullptr) {
            cudaFree(data_);
            data_ = nullptr;
            capacity_ = 0;
        }
    }

private:
    T* data_ = nullptr;
    std::size_t capacity_ = 0;
};

class Event {
public:
    Event() { checkCuda(cudaEventCreate(&event_), "cudaEventCreate"); }
    ~Event() { cudaEventDestroy(event_); }
    Event(const Event&) = delete;
    Event& operator=(const Event&) = delete;
    cudaEvent_t get() const noexcept { return event_; }

private:
    cudaEvent_t event_{};
};

__global__ void accumulatePerBaseKernel(
    const char* sequences,
    const char* qualities,
    const unsigned long long* offsets,
    const std::uint32_t* lengths,
    std::size_t readCount,
    unsigned long long* countA,
    unsigned long long* countC,
    unsigned long long* countG,
    unsigned long long* countT,
    unsigned long long* countN,
    unsigned long long* readsPerPosition,
    unsigned long long* qualitySum,
    unsigned long long* qualityCount,
    unsigned long long* qualityHistogram)
{
    const std::size_t position = blockIdx.x;
    const std::size_t readIndex = static_cast<std::size_t>(blockIdx.y) * blockDim.x + threadIdx.x;

    __shared__ unsigned long long sharedBases[5];
    __shared__ unsigned long long sharedReads;
    __shared__ unsigned long long sharedQualitySum;
    __shared__ unsigned long long sharedQualityCount;
    __shared__ unsigned long long sharedHistogram[94];

    for (std::size_t index = threadIdx.x; index < 5; index += blockDim.x) {
        sharedBases[index] = 0;
    }
    for (std::size_t index = threadIdx.x; index < 94; index += blockDim.x) {
        sharedHistogram[index] = 0;
    }
    if (threadIdx.x == 0) {
        sharedReads = 0;
        sharedQualitySum = 0;
        sharedQualityCount = 0;
    }
    __syncthreads();

    if (readIndex < readCount && position < lengths[readIndex]) {
        const unsigned long long offset = offsets[readIndex] + position;
        atomicAdd(&sharedReads, 1ULL);

        switch (sequences[offset]) {
        case 'A': case 'a': atomicAdd(&sharedBases[0], 1ULL); break;
        case 'C': case 'c': atomicAdd(&sharedBases[1], 1ULL); break;
        case 'G': case 'g': atomicAdd(&sharedBases[2], 1ULL); break;
        case 'T': case 't': atomicAdd(&sharedBases[3], 1ULL); break;
        case 'N': case 'n': atomicAdd(&sharedBases[4], 1ULL); break;
        default: break;
        }

        const int quality = static_cast<int>(qualities[offset]) - 33;
        if (quality >= 0) {
            atomicAdd(&sharedQualitySum, static_cast<unsigned long long>(quality));
            atomicAdd(&sharedQualityCount, 1ULL);
            if (quality < 94) {
                atomicAdd(&sharedHistogram[quality], 1ULL);
            }
        }
    }
    __syncthreads();

    if (threadIdx.x == 0) {
        atomicAdd(&countA[position], sharedBases[0]);
        atomicAdd(&countC[position], sharedBases[1]);
        atomicAdd(&countG[position], sharedBases[2]);
        atomicAdd(&countT[position], sharedBases[3]);
        atomicAdd(&countN[position], sharedBases[4]);
        atomicAdd(&readsPerPosition[position], sharedReads);
        atomicAdd(&qualitySum[position], sharedQualitySum);
        atomicAdd(&qualityCount[position], sharedQualityCount);
        for (std::size_t quality = 0; quality < 94; ++quality) {
            if (sharedHistogram[quality] != 0) {
                atomicAdd(&qualityHistogram[position * 94 + quality], sharedHistogram[quality]);
            }
        }
    }
}

} // namespace

struct CudaQcProcessor::Impl {
    DeviceBuffer<char> sequences;
    DeviceBuffer<char> qualities;
    DeviceBuffer<std::uint64_t> offsets;
    DeviceBuffer<std::uint32_t> lengths;
    DeviceBuffer<std::uint64_t> countA;
    DeviceBuffer<std::uint64_t> countC;
    DeviceBuffer<std::uint64_t> countG;
    DeviceBuffer<std::uint64_t> countT;
    DeviceBuffer<std::uint64_t> countN;
    DeviceBuffer<std::uint64_t> readsPerPosition;
    DeviceBuffer<std::uint64_t> qualitySum;
    DeviceBuffer<std::uint64_t> qualityCount;
    DeviceBuffer<std::uint64_t> qualityHistogram;
    CudaQcTiming timing;
};

CudaQcTiming& CudaQcTiming::operator+=(const CudaQcTiming& other)
{
    packing += other.packing;
    hostToDevice += other.hostToDevice;
    kernel += other.kernel;
    deviceToHost += other.deviceToHost;
    total += other.total;
    return *this;
}

CudaQcProcessor::CudaQcProcessor()
    : impl_(std::make_unique<Impl>())
{
}

CudaQcProcessor::~CudaQcProcessor() = default;
CudaQcProcessor::CudaQcProcessor(CudaQcProcessor&&) noexcept = default;
CudaQcProcessor& CudaQcProcessor::operator=(CudaQcProcessor&&) noexcept = default;

const CudaQcTiming& CudaQcProcessor::lastTiming() const noexcept
{
    return impl_->timing;
}

PerBaseMetrics CudaQcProcessor::process(
    const std::vector<FastqRecord>& records,
    const std::vector<std::uint8_t>* selected)
{
    if (selected != nullptr && selected->size() != records.size()) {
        throw std::runtime_error("CUDA QC selection mask has an invalid size.");
    }

    const auto started = Clock::now();
    impl_->timing = {};
    std::vector<char> hostSequences;
    std::vector<char> hostQualities;
    std::vector<std::uint64_t> hostOffsets;
    std::vector<std::uint32_t> hostLengths;
    std::size_t maxLength = 0;

    const auto packingStarted = Clock::now();
    for (std::size_t index = 0; index < records.size(); ++index) {
        if (selected != nullptr && !(*selected)[index]) {
            continue;
        }
        const FastqRecord& record = records[index];
        const std::size_t length = record.sequence.size();
        if (length == 0) {
            continue;
        }
        if (length != record.quality.size()) {
            throw std::runtime_error("CUDA QC received a FASTQ record with unequal sequence and quality lengths.");
        }
        if (length > std::numeric_limits<std::uint32_t>::max()) {
            throw std::runtime_error("CUDA QC read length exceeds the uint32_t packed representation.");
        }
        if (hostSequences.size() > std::numeric_limits<std::uint64_t>::max()) {
            throw std::runtime_error("CUDA QC packed sequence offset exceeds uint64_t.");
        }
        hostOffsets.push_back(static_cast<std::uint64_t>(hostSequences.size()));
        hostLengths.push_back(static_cast<std::uint32_t>(length));
        hostSequences.insert(hostSequences.end(), record.sequence.begin(), record.sequence.end());
        hostQualities.insert(hostQualities.end(), record.quality.begin(), record.quality.end());
        maxLength = std::max(maxLength, length);
    }
    impl_->timing.packing = std::chrono::duration<double>(Clock::now() - packingStarted).count();

    PerBaseMetrics result;
    if (hostOffsets.empty()) {
        impl_->timing.total = std::chrono::duration<double>(Clock::now() - started).count();
        return result;
    }

    const std::size_t positionsBytes = checkedMultiply(maxLength, sizeof(std::uint64_t), "per-position output");
    const std::size_t histogramElements = checkedMultiply(maxLength, std::size_t{94}, "quality histogram");
    const std::size_t histogramBytes = checkedMultiply(histogramElements, sizeof(std::uint64_t), "quality histogram");
    std::size_t requiredBytes = checkedAdd(hostSequences.size(), hostQualities.size(), "packed reads");
    requiredBytes = checkedAdd(requiredBytes, checkedMultiply(hostOffsets.size(), sizeof(std::uint64_t), "offsets"), "offsets");
    requiredBytes = checkedAdd(requiredBytes, checkedMultiply(hostLengths.size(), sizeof(std::uint32_t), "lengths"), "lengths");
    requiredBytes = checkedAdd(requiredBytes, checkedMultiply(positionsBytes, std::size_t{9}, "per-position outputs"), "per-position outputs");
    requiredBytes = checkedAdd(requiredBytes, histogramBytes, "quality histogram");

    std::size_t freeBytes = 0;
    std::size_t totalBytes = 0;
    checkCuda(cudaMemGetInfo(&freeBytes, &totalBytes), "cudaMemGetInfo");
    if (requiredBytes > freeBytes) {
        std::ostringstream message;
        message << "CUDA QC batch requires " << requiredBytes
                << " bytes of device memory, but only " << freeBytes << " bytes are free.";
        throw std::runtime_error(message.str());
    }

    impl_->sequences.ensure(hostSequences.size(), "cudaMalloc packed sequences");
    impl_->qualities.ensure(hostQualities.size(), "cudaMalloc packed qualities");
    impl_->offsets.ensure(hostOffsets.size(), "cudaMalloc packed offsets");
    impl_->lengths.ensure(hostLengths.size(), "cudaMalloc packed lengths");
    impl_->countA.ensure(maxLength, "cudaMalloc base A counts");
    impl_->countC.ensure(maxLength, "cudaMalloc base C counts");
    impl_->countG.ensure(maxLength, "cudaMalloc base G counts");
    impl_->countT.ensure(maxLength, "cudaMalloc base T counts");
    impl_->countN.ensure(maxLength, "cudaMalloc base N counts");
    impl_->readsPerPosition.ensure(maxLength, "cudaMalloc reads-per-position counts");
    impl_->qualitySum.ensure(maxLength, "cudaMalloc quality sums");
    impl_->qualityCount.ensure(maxLength, "cudaMalloc quality counts");
    impl_->qualityHistogram.ensure(histogramElements, "cudaMalloc quality histogram");

    Event transferStarted;
    Event transferFinished;
    checkCuda(cudaEventRecord(transferStarted.get()), "cudaEventRecord H2D start");
    checkCuda(cudaMemcpy(impl_->sequences.data(), hostSequences.data(), hostSequences.size(), cudaMemcpyHostToDevice), "cudaMemcpy sequences H2D");
    checkCuda(cudaMemcpy(impl_->qualities.data(), hostQualities.data(), hostQualities.size(), cudaMemcpyHostToDevice), "cudaMemcpy qualities H2D");
    checkCuda(cudaMemcpy(impl_->offsets.data(), hostOffsets.data(), hostOffsets.size() * sizeof(std::uint64_t), cudaMemcpyHostToDevice), "cudaMemcpy offsets H2D");
    checkCuda(cudaMemcpy(impl_->lengths.data(), hostLengths.data(), hostLengths.size() * sizeof(std::uint32_t), cudaMemcpyHostToDevice), "cudaMemcpy lengths H2D");
    checkCuda(cudaEventRecord(transferFinished.get()), "cudaEventRecord H2D finish");
    checkCuda(cudaEventSynchronize(transferFinished.get()), "cudaEventSynchronize H2D");
    float milliseconds = 0.0F;
    checkCuda(cudaEventElapsedTime(&milliseconds, transferStarted.get(), transferFinished.get()), "cudaEventElapsedTime H2D");
    impl_->timing.hostToDevice = milliseconds / 1000.0;

    checkCuda(cudaMemset(impl_->countA.data(), 0, positionsBytes), "cudaMemset base A counts");
    checkCuda(cudaMemset(impl_->countC.data(), 0, positionsBytes), "cudaMemset base C counts");
    checkCuda(cudaMemset(impl_->countG.data(), 0, positionsBytes), "cudaMemset base G counts");
    checkCuda(cudaMemset(impl_->countT.data(), 0, positionsBytes), "cudaMemset base T counts");
    checkCuda(cudaMemset(impl_->countN.data(), 0, positionsBytes), "cudaMemset base N counts");
    checkCuda(cudaMemset(impl_->readsPerPosition.data(), 0, positionsBytes), "cudaMemset reads-per-position counts");
    checkCuda(cudaMemset(impl_->qualitySum.data(), 0, positionsBytes), "cudaMemset quality sums");
    checkCuda(cudaMemset(impl_->qualityCount.data(), 0, positionsBytes), "cudaMemset quality counts");
    checkCuda(cudaMemset(impl_->qualityHistogram.data(), 0, histogramBytes), "cudaMemset quality histogram");

    Event kernelStarted;
    Event kernelFinished;
    checkCuda(cudaEventRecord(kernelStarted.get()), "cudaEventRecord kernel start");
    constexpr unsigned int threads = 256;
    const std::size_t yBlocks = (hostOffsets.size() + threads - 1) / threads;
    if (maxLength > std::numeric_limits<unsigned int>::max()
        || yBlocks > std::numeric_limits<unsigned int>::max()) {
        throw std::runtime_error("CUDA QC grid dimensions exceed CUDA limits.");
    }
    accumulatePerBaseKernel<<<dim3(static_cast<unsigned int>(maxLength), static_cast<unsigned int>(yBlocks)), threads>>>(
        impl_->sequences.data(), impl_->qualities.data(),
        reinterpret_cast<const unsigned long long*>(impl_->offsets.data()), impl_->lengths.data(),
        hostOffsets.size(),
        reinterpret_cast<unsigned long long*>(impl_->countA.data()),
        reinterpret_cast<unsigned long long*>(impl_->countC.data()),
        reinterpret_cast<unsigned long long*>(impl_->countG.data()),
        reinterpret_cast<unsigned long long*>(impl_->countT.data()),
        reinterpret_cast<unsigned long long*>(impl_->countN.data()),
        reinterpret_cast<unsigned long long*>(impl_->readsPerPosition.data()),
        reinterpret_cast<unsigned long long*>(impl_->qualitySum.data()),
        reinterpret_cast<unsigned long long*>(impl_->qualityCount.data()),
        reinterpret_cast<unsigned long long*>(impl_->qualityHistogram.data()));
    checkCuda(cudaGetLastError(), "CUDA per-base QC kernel launch");
    checkCuda(cudaEventRecord(kernelFinished.get()), "cudaEventRecord kernel finish");
    checkCuda(cudaEventSynchronize(kernelFinished.get()), "cudaEventSynchronize kernel");
    checkCuda(cudaEventElapsedTime(&milliseconds, kernelStarted.get(), kernelFinished.get()), "cudaEventElapsedTime kernel");
    impl_->timing.kernel = milliseconds / 1000.0;

    result.resize(maxLength);
    Event resultStarted;
    Event resultFinished;
    checkCuda(cudaEventRecord(resultStarted.get()), "cudaEventRecord D2H start");
    checkCuda(cudaMemcpy(result.baseCountA.data(), impl_->countA.data(), positionsBytes, cudaMemcpyDeviceToHost), "cudaMemcpy base A counts D2H");
    checkCuda(cudaMemcpy(result.baseCountC.data(), impl_->countC.data(), positionsBytes, cudaMemcpyDeviceToHost), "cudaMemcpy base C counts D2H");
    checkCuda(cudaMemcpy(result.baseCountG.data(), impl_->countG.data(), positionsBytes, cudaMemcpyDeviceToHost), "cudaMemcpy base G counts D2H");
    checkCuda(cudaMemcpy(result.baseCountT.data(), impl_->countT.data(), positionsBytes, cudaMemcpyDeviceToHost), "cudaMemcpy base T counts D2H");
    checkCuda(cudaMemcpy(result.baseCountN.data(), impl_->countN.data(), positionsBytes, cudaMemcpyDeviceToHost), "cudaMemcpy base N counts D2H");
    checkCuda(cudaMemcpy(result.readsPerPosition.data(), impl_->readsPerPosition.data(), positionsBytes, cudaMemcpyDeviceToHost), "cudaMemcpy reads-per-position D2H");
    checkCuda(cudaMemcpy(result.qualitySum.data(), impl_->qualitySum.data(), positionsBytes, cudaMemcpyDeviceToHost), "cudaMemcpy quality sums D2H");
    checkCuda(cudaMemcpy(result.qualityCount.data(), impl_->qualityCount.data(), positionsBytes, cudaMemcpyDeviceToHost), "cudaMemcpy quality counts D2H");
    checkCuda(cudaMemcpy(result.qualityHistogram.data(), impl_->qualityHistogram.data(), histogramBytes, cudaMemcpyDeviceToHost), "cudaMemcpy quality histogram D2H");
    checkCuda(cudaEventRecord(resultFinished.get()), "cudaEventRecord D2H finish");
    checkCuda(cudaEventSynchronize(resultFinished.get()), "cudaEventSynchronize D2H");
    checkCuda(cudaEventElapsedTime(&milliseconds, resultStarted.get(), resultFinished.get()), "cudaEventElapsedTime D2H");
    impl_->timing.deviceToHost = milliseconds / 1000.0;
    impl_->timing.total = std::chrono::duration<double>(Clock::now() - started).count();
    return result;
}
