#include "trimming/trimmer.h"
#include "fastq_reader.h"

#include <stdexcept>
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
    adapterTrimmer_.trim(record, result);
    applyMinimumLengthFilter(result);

    updateStats(result);
    return result;
}

std::pair<TrimResult, TrimResult> Trimmer::trimPair(FastqRecord& r1,
                                                     FastqRecord& r2) {
    TrimResult r1Result = qualityTrimmer_.trim(r1);
    TrimResult r2Result = qualityTrimmer_.trim(r2);
    adapterTrimmer_.trimPair(r1, r2, r1Result, r2Result);
    applyMinimumLengthFilter(r1Result);
    applyMinimumLengthFilter(r2Result);

    updateStats(r1Result);
    updateStats(r2Result);
    return {std::move(r1Result), std::move(r2Result)};
}

void Trimmer::applyMinimumLengthFilter(TrimResult& result) const {
    result.passed = true;
    result.discard_reason.reset();

    if (config_.enabled && result.final_length < config_.min_length) {
        result.passed = false;
        result.discard_reason = "TOO_SHORT";
    }
}

void Trimmer::updateStats(const TrimResult& result) {
    ++stats_.total_reads;
    stats_.bases_before += result.original_length;
    stats_.bases_after += result.final_length;
    stats_.bases_trimmed += result.original_length - result.final_length;
    if (result.passed) {
        ++stats_.passed_reads;
    } else {
        ++stats_.discarded_reads;
        if (result.discard_reason == "TOO_SHORT") {
            ++stats_.too_short_reads;
        }
    }

    if (result.quality_trimmed) {
        ++stats_.quality_trimmed_reads;
    }
    if (result.adapter_found && result.adapter_position) {
        ++stats_.adapter_trimmed_reads;
        ++stats_.adapter_positions[*result.adapter_position];
    }
}
