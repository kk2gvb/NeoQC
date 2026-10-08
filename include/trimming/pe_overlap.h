#pragma once

#include "trimming/trim_config.h"

#include <cstddef>
#include <optional>
#include <string>

struct FastqRecord;

struct PeOverlap {
    std::size_t r1_start = 0;
    std::size_t r2_reverse_start = 0;
    std::size_t length = 0;
    std::size_t mismatches = 0;

    // Read-through means that both reads cover the complete insert and carry
    // non-insert sequence at their 3' ends.
    bool isReadThrough() const noexcept {
        return r1_start == 0 && r2_reverse_start > 0;
    }
};

class PeOverlapAnalyzer {
public:
    explicit PeOverlapAnalyzer(TrimConfig config);

    std::optional<PeOverlap> find(const FastqRecord& r1,
                                  const FastqRecord& r2) const;

    static char normalizedBase(char base) noexcept;
    static char complement(char base) noexcept;
    static std::string reverseComplement(const std::string& sequence);
    static std::string reverseQuality(const std::string& quality);
    static std::size_t originalR2Index(std::size_t r2Length,
                                       std::size_t reverseIndex);

private:
    TrimConfig config_;
};
