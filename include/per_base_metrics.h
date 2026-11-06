#pragma once

#include "fastq_reader.h"

#include <array>
#include <cstdint>
#include <vector>

// Exact position-level accumulators shared by the CPU reference fallback and
// the CUDA implementation.  The layout matches QualityAnalyzer's existing
// per-base state; report generation therefore remains unchanged.
struct PerBaseMetrics {
    std::vector<std::uint64_t> baseCountA;
    std::vector<std::uint64_t> baseCountC;
    std::vector<std::uint64_t> baseCountG;
    std::vector<std::uint64_t> baseCountT;
    std::vector<std::uint64_t> baseCountN;
    std::vector<std::uint64_t> readsPerPosition;

    std::vector<std::uint64_t> qualitySum;
    std::vector<std::uint64_t> qualityCount;
    std::vector<std::array<std::uint64_t, 94>> qualityHistogram;

    void resize(std::size_t positions);
    void merge(const PerBaseMetrics& other);
};

// CPU reference used only when a CUDA batch cannot be submitted safely.  The
// selected mask represents a transactional batch subset, e.g. reads that
// passed trimming.  A null mask selects every record.
PerBaseMetrics collectPerBaseMetricsCpu(
    const std::vector<FastqRecord>& records,
    const std::vector<std::uint8_t>* selected = nullptr);
