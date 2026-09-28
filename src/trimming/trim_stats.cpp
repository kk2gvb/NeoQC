#include "trimming/trim_stats.h"

void TrimStats::merge(const TrimStats& other) {
    total_reads += other.total_reads;
    passed_reads += other.passed_reads;
    discarded_reads += other.discarded_reads;
    bases_before += other.bases_before;
    bases_after += other.bases_after;
    bases_trimmed += other.bases_trimmed;
    adapter_trimmed_reads += other.adapter_trimmed_reads;
    quality_trimmed_reads += other.quality_trimmed_reads;
    polyG_trimmed_reads += other.polyG_trimmed_reads;
    polyX_trimmed_reads += other.polyX_trimmed_reads;
    too_short_reads += other.too_short_reads;

    for (const auto& [position, count] : other.adapter_positions) {
        adapter_positions[position] += count;
    }
}
