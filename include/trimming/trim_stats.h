#pragma once

#include <cstdint>

// Aggregate statistics for the trimming pipeline.
struct TrimStats {
    std::uint64_t total_reads = 0;
    std::uint64_t passed_reads = 0;
    std::uint64_t discarded_reads = 0;

    std::uint64_t bases_before = 0;
    std::uint64_t bases_after = 0;
    std::uint64_t bases_trimmed = 0;

    std::uint64_t adapter_trimmed_reads = 0;
    std::uint64_t quality_trimmed_reads = 0;
    std::uint64_t too_short_reads = 0;
};
