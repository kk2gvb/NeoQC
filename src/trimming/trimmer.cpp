#include "trimming/trimmer.h"
#include "fastq_reader.h"

#include <utility>

Trimmer::Trimmer(TrimConfig config)
    : config_(std::move(config))
    , qualityTrimmer_(config_)
    , adapterTrimmer_(config_)
{
}

const TrimConfig& Trimmer::getConfig() const noexcept {
    return config_;
}

const TrimStats& Trimmer::getStats() const noexcept {
    return stats_;
}

const QualityTrimmer& Trimmer::getQualityTrimmer() const noexcept {
    return qualityTrimmer_;
}

const AdapterTrimmer& Trimmer::getAdapterTrimmer() const noexcept {
    return adapterTrimmer_;
}

TrimResult Trimmer::trim(FastqRecord& record) {
    TrimResult result = qualityTrimmer_.trim(record);

    ++stats_.total_reads;
    stats_.bases_before += result.original_length;
    stats_.bases_after += result.final_length;
    stats_.bases_trimmed += result.original_length - result.final_length;
    ++stats_.passed_reads;

    if (result.quality_trimmed) {
        ++stats_.quality_trimmed_reads;
    }

    return result;
}
