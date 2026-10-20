#pragma once

#include <cstddef>
#include <cstdint>
#include <map>

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
    std::uint64_t polyG_trimmed_reads = 0;
    std::uint64_t polyX_trimmed_reads = 0;
    std::uint64_t too_short_reads = 0;

    // Output counters preserve read-level filtering semantics for PE while
    // making the strict pair policy explicit.
    std::uint64_t output_reads = 0;
    std::uint64_t mate_discarded_reads = 0;
    std::uint64_t total_pairs = 0;
    std::uint64_t passed_pairs = 0;
    std::uint64_t discarded_pairs = 0;
    std::uint64_t corrected_reads = 0;
    std::uint64_t corrected_bases = 0;
    std::uint64_t merged_pairs = 0;
    std::uint64_t umi_processed_reads = 0;
    std::uint64_t umi_failed_reads = 0;
    std::uint64_t umi_bases = 0;

    // Adapter positions are zero-based offsets in the read before adapter trimming.
    std::map<std::size_t, std::uint64_t> adapter_positions;

    void merge(const TrimStats& other);
};
