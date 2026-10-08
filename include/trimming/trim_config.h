#pragma once

#include "adapter_config.h"

#include <cstddef>
#include <string>
#include <vector>

// Trimming is disabled by default. Individual operations only take effect
// when enabled is true.
struct TrimConfig {
    bool enabled = false;

    std::size_t trim_front = 0;
    std::size_t trim_tail = 0;
    bool cut_front = false;
    bool cut_tail = false;
    bool cut_right = false;

    std::size_t quality_threshold = 20;
    std::size_t window_size = 4;

    bool adapter_trimming = false;
    std::string adapter_sequence;
    std::vector<AdapterConfigEntry> adapters;
    std::size_t min_adapter_match = 6;

    std::size_t min_overlap = 30;
    std::size_t overlap_diff_limit = 5;
    std::size_t overlap_diff_percent_limit = 20;

    std::size_t min_length = 0;

    bool trim_poly_g = false;
    bool trim_poly_x = false;
    std::size_t poly_g_min_length = 10;
    std::size_t poly_x_min_length = 10;

    bool overlap_correction = false;
    bool merge_reads = false;

    bool umi_enabled = false;
    std::size_t umi_length = 0;
};
