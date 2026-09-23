#include "trimming/adapter_trimmer.h"
#include "trimming/quality_trimmer.h"
#include "trimming/trim_config.h"
#include "trimming/trim_result.h"
#include "trimming/trim_stats.h"
#include "trimming/trimmer.h"

#include <iostream>

int main() {
    const TrimConfig config;
    if (config.enabled) {
        std::cerr << "Trimming must be disabled by default\n";
        return 1;
    }

    const TrimStats stats;
    if (stats.total_reads != 0 || stats.passed_reads != 0
        || stats.discarded_reads != 0 || stats.bases_before != 0
        || stats.bases_after != 0 || stats.bases_trimmed != 0
        || stats.adapter_trimmed_reads != 0
        || stats.quality_trimmed_reads != 0 || stats.too_short_reads != 0) {
        std::cerr << "TrimStats counters must default to zero\n";
        return 1;
    }

    const QualityTrimmer qualityTrimmer(config);
    const AdapterTrimmer adapterTrimmer(config);
    if (qualityTrimmer.getConfig().enabled
        || adapterTrimmer.getConfig().enabled) {
        std::cerr << "Trimming components must preserve disabled configuration\n";
        return 1;
    }

    TrimConfig enabledConfig;
    enabledConfig.enabled = true;
    const Trimmer trimmer(enabledConfig);
    if (!trimmer.getConfig().enabled
        || !trimmer.getQualityTrimmer().getConfig().enabled
        || !trimmer.getAdapterTrimmer().getConfig().enabled
        || trimmer.getStats().total_reads != 0) {
        std::cerr << "Trimmer must compose its configuration and empty statistics\n";
        return 1;
    }

    const TrimResult result;
    if (result.quality_trimmed || result.adapter_found || result.adapter_position.has_value()
        || result.discard_reason.has_value()) {
        std::cerr << "TrimResult default optional state is invalid\n";
        return 1;
    }

    TrimResult adapterAtStart;
    adapterAtStart.adapter_found = true;
    adapterAtStart.adapter_position = 0;
    if (!adapterAtStart.adapter_position.has_value()
        || *adapterAtStart.adapter_position != 0) {
        std::cerr << "Adapter position zero must be distinguishable from no adapter\n";
        return 1;
    }

    return 0;
}
