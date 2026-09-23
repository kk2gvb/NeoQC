#pragma once

#include "trimming/adapter_trimmer.h"
#include "trimming/quality_trimmer.h"
#include "trimming/trim_config.h"
#include "trimming/trim_result.h"
#include "trimming/trim_stats.h"

struct FastqRecord;

// Top-level owner of trimming pipeline configuration and statistics.
class Trimmer {
public:
    explicit Trimmer(TrimConfig config);

    const TrimConfig& getConfig() const noexcept;
    const TrimStats& getStats() const noexcept;
    const QualityTrimmer& getQualityTrimmer() const noexcept;
    const AdapterTrimmer& getAdapterTrimmer() const noexcept;

    TrimResult trim(FastqRecord& record);

private:
    TrimConfig config_;
    QualityTrimmer qualityTrimmer_;
    AdapterTrimmer adapterTrimmer_;
    TrimStats stats_;
};
