#pragma once

#include "trimming/trim_config.h"
#include "trimming/trim_result.h"

struct FastqRecord;

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

private:
    TrimConfig config_;
};
