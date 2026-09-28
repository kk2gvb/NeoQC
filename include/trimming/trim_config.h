#pragma once

#include <cstddef>
#include <string>

// Configuration reserved for the future trimming pipeline.
// All operations are disabled by default.
struct TrimConfig {
    bool enabled = false;

    std::size_t trim_front = 0;
    std::size_t trim_tail = 0;
    bool cut_front = false;
    bool cut_tail = false;
    bool cut_right = false;

    std::size_t quality_threshold = 0;
    std::size_t window_size = 0;

    bool adapter_trimming = false;
    std::string adapter_sequence;
    std::size_t min_adapter_match = 6;

    std::size_t min_overlap = 30;
    std::size_t overlap_diff_limit = 5;
    std::size_t overlap_diff_percent_limit = 20;

    std::size_t min_length = 0;
};
