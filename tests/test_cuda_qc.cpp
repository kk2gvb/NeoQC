#include "cuda/cuda_qc.h"
#include "per_base_metrics.h"

#include <cassert>
#include <exception>
#include <iostream>
#include <vector>

namespace {

FastqRecord record(std::string sequence, std::string quality, std::uint64_t number)
{
    return {"@CUDA_" + std::to_string(number), std::move(sequence), "+", std::move(quality), number};
}

void assertEqual(const PerBaseMetrics& expected, const PerBaseMetrics& actual)
{
    assert(expected.baseCountA == actual.baseCountA);
    assert(expected.baseCountC == actual.baseCountC);
    assert(expected.baseCountG == actual.baseCountG);
    assert(expected.baseCountT == actual.baseCountT);
    assert(expected.baseCountN == actual.baseCountN);
    assert(expected.readsPerPosition == actual.readsPerPosition);
    assert(expected.qualitySum == actual.qualitySum);
    assert(expected.qualityCount == actual.qualityCount);
    assert(expected.qualityHistogram == actual.qualityHistogram);
}

} // namespace

int main()
{
    const std::vector<FastqRecord> records = {
        record("ACGTN", "!5I~+", 1),
        record("ac", "~!", 2),
        record("NNNN", "IIII", 3),
    };
    const std::vector<std::uint8_t> selected = {1, 0, 1};

    try {
        CudaQcProcessor processor;
        assertEqual(collectPerBaseMetricsCpu(records), processor.process(records));
        assertEqual(collectPerBaseMetricsCpu(records, &selected), processor.process(records, &selected));
        std::cout << "CUDA per-base QC parity passed.\n";
        return 0;
    }
    catch (const std::exception& error) {
        std::cerr << "CUDA per-base QC parity could not run: " << error.what() << '\n';
        return 77;
    }
}
