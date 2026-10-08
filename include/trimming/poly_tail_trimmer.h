#pragma once

#include "trimming/trim_config.h"
#include "trimming/trim_result.h"

struct FastqRecord;

// Removes configured terminal polyG and polyX runs after adapter trimming.
class PolyTailTrimmer {
public:
    explicit PolyTailTrimmer(TrimConfig config);

    const TrimConfig& getConfig() const noexcept;
    void trim(FastqRecord& record, TrimResult& result) const;

private:
    TrimConfig config_;
};
