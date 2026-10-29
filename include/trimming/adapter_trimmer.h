#pragma once

#include "trimming/trim_config.h"
#include "trimming/trim_result.h"

#include <optional>
#include <string>
#include <vector>

struct FastqRecord;
struct PeOverlap;

// Applies explicit adapter matching and paired-end overlap adapter trimming.
class AdapterTrimmer {
public:
    explicit AdapterTrimmer(TrimConfig config);

    const TrimConfig& getConfig() const noexcept;

    TrimResult trim(FastqRecord& record) const;
    void trim(FastqRecord& record, TrimResult& result) const;

    void trimPair(FastqRecord& r1,
                  FastqRecord& r2,
                  TrimResult& r1Result,
                  TrimResult& r2Result) const;
    void trimPair(FastqRecord& r1,
                  FastqRecord& r2,
                  TrimResult& r1Result,
                  TrimResult& r2Result,
                  const std::optional<PeOverlap>& overlap) const;

private:
    TrimConfig config_;
    // Normalized once at construction, so the matching hot loop only
    // normalizes bases from the incoming read.
    std::vector<std::string> adapters_;
};
