#pragma once

#include "trimming/trim_config.h"
#include "trimming/trim_result.h"

struct FastqRecord;

// Applies fixed and end-based quality trimming to one FASTQ record.
class QualityTrimmer {
public:
    explicit QualityTrimmer(TrimConfig config);

    const TrimConfig& getConfig() const noexcept;

    TrimResult trim(FastqRecord& record) const;

private:
    TrimConfig config_;
};
