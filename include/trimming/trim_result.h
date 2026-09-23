#pragma once

#include <cstddef>
#include <optional>
#include <string>

// Result of trimming a single FASTQ read.
struct TrimResult {
    std::size_t original_length = 0;
    std::size_t final_length = 0;
    std::size_t trimmed_front = 0;
    std::size_t trimmed_tail = 0;

    bool quality_trimmed = false;

    bool adapter_found = false;
    std::optional<std::size_t> adapter_position;

    bool passed = true;
    std::optional<std::string> discard_reason;
};
