#include "trimming/trimmer.h"
#include "fastq_reader.h"

#include <stdexcept>
#include <utility>

Trimmer::Trimmer(TrimConfig config)
    : config_(std::move(config))
    , qualityTrimmer_(config_)
    , adapterTrimmer_(config_)
    , polyTailTrimmer_(config_)
    , umiProcessor_(config_)
    , overlapAnalyzer_(config_)
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

const PolyTailTrimmer& Trimmer::getPolyTailTrimmer() const noexcept {
    return polyTailTrimmer_;
}

TrimResult Trimmer::trim(FastqRecord& record) {
    TrimResult umiResult;
    if (!umiProcessor_.process(record, umiResult)) {
        updateStats(umiResult);
        return umiResult;
    }
    TrimResult result = qualityTrimmer_.trim(record);
    result.original_length = umiResult.original_length;
    result.trimmed_front += umiResult.trimmed_front;
    result.umi_extracted = umiResult.umi_extracted;
    result.umi = std::move(umiResult.umi);
    adapterTrimmer_.trim(record, result);
    polyTailTrimmer_.trim(record, result);
    applyMinimumLengthFilter(result);

    updateStats(result);
    if (result.passed) {
        ++stats_.output_reads;
    }
    return result;
}

std::pair<TrimResult, TrimResult> Trimmer::trimPair(FastqRecord& r1,
                                                     FastqRecord& r2) {
    PairTrimResult detailed = trimPairDetailed(r1, r2);
    return {std::move(detailed.r1), std::move(detailed.r2)};
}

PairTrimResult Trimmer::trimPairDetailed(FastqRecord& r1, FastqRecord& r2) {
    PairTrimResult output;
    TrimResult r1Umi;
    const bool r1Ready = umiProcessor_.process(r1, r1Umi);

    if (r1Ready) {
        output.r1 = qualityTrimmer_.trim(r1);
        output.r1.original_length = r1Umi.original_length;
        output.r1.trimmed_front += r1Umi.trimmed_front;
        output.r1.umi_extracted = r1Umi.umi_extracted;
        output.r1.umi = std::move(r1Umi.umi);
    } else {
        output.r1 = std::move(r1Umi);
    }
    output.r2 = qualityTrimmer_.trim(r2);

    std::optional<PeOverlap> overlap;
    if (r1Ready && (config_.adapter_trimming || config_.overlap_correction
                    || config_.merge_reads)) {
        overlap = overlapAnalyzer_.find(r1, r2);
    }
    if (r1Ready && config_.overlap_correction && overlap.has_value()) {
        peProcessor_.correct(r1, r2, output.r1, output.r2, *overlap);
    }
    if (r1Ready) {
        adapterTrimmer_.trimPair(
            r1, r2, output.r1, output.r2, overlap);
        polyTailTrimmer_.trim(r1, output.r1);
    }
    polyTailTrimmer_.trim(r2, output.r2);
    applyMinimumLengthFilter(output.r1);
    applyMinimumLengthFilter(output.r2);

    if (config_.merge_reads && output.r1.passed && output.r2.passed) {
        const auto mergeOverlap = overlapAnalyzer_.find(r1, r2);
        if (mergeOverlap.has_value()) {
            output.merged = peProcessor_.merge(r1, r2, *mergeOverlap);
            ++stats_.merged_pairs;
        }
    }

    updateStats(output.r1);
    updateStats(output.r2);
    ++stats_.total_pairs;
    if (output.r1.passed && output.r2.passed) {
        ++stats_.passed_pairs;
        stats_.output_reads += 2;
    } else {
        ++stats_.discarded_pairs;
        if (output.r1.passed) {
            ++stats_.mate_discarded_reads;
        }
        if (output.r2.passed) {
            ++stats_.mate_discarded_reads;
        }
    }
    return output;
}

void Trimmer::applyMinimumLengthFilter(TrimResult& result) const {
    if (!result.passed && result.discard_reason.has_value()) {
        return;
    }
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
    if (result.polyG_trimmed) {
        ++stats_.polyG_trimmed_reads;
    }
    if (result.polyX_trimmed) {
        ++stats_.polyX_trimmed_reads;
    }
    if (result.corrected_bases > 0) {
        ++stats_.corrected_reads;
        stats_.corrected_bases += result.corrected_bases;
    }
    if (result.umi_extracted) {
        ++stats_.umi_processed_reads;
        stats_.umi_bases += result.umi.size();
    } else if (result.discard_reason == "UMI_TOO_SHORT") {
        ++stats_.umi_failed_reads;
    }
}
