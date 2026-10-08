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
    output_reads += other.output_reads;
    mate_discarded_reads += other.mate_discarded_reads;
    total_pairs += other.total_pairs;
    passed_pairs += other.passed_pairs;
    discarded_pairs += other.discarded_pairs;
    corrected_reads += other.corrected_reads;
    corrected_bases += other.corrected_bases;
    merged_pairs += other.merged_pairs;
    umi_processed_reads += other.umi_processed_reads;
    umi_failed_reads += other.umi_failed_reads;
    umi_bases += other.umi_bases;

    for (const auto& [position, count] : other.adapter_positions) {
        adapter_positions[position] += count;
    }
}
