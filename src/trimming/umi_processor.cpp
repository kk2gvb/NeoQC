#include "trimming/umi_processor.h"
#include "fastq_reader.h"

#include <stdexcept>
#include <utility>

UmiProcessor::UmiProcessor(TrimConfig config)
    : config_(std::move(config))
{
}

bool UmiProcessor::process(FastqRecord& record, TrimResult& result) const {
    if (record.sequence.size() != record.quality.size()) {
        throw std::invalid_argument(
            "Cannot extract UMI from different sequence and quality lengths");
    }
    result.original_length = record.sequence.size();
    result.final_length = record.sequence.size();
    if (!config_.enabled || !config_.umi_enabled) {
        return true;
    }
    if (config_.umi_length == 0) {
        throw std::invalid_argument("UMI length must be greater than zero");
    }
    if (record.sequence.size() < config_.umi_length) {
        result.passed = false;
        result.discard_reason = "UMI_TOO_SHORT";
        return false;
    }

    result.umi = record.sequence.substr(0, config_.umi_length);
    result.umi_extracted = true;
    record.header += " UMI:" + result.umi;
    record.sequence.erase(0, config_.umi_length);
    record.quality.erase(0, config_.umi_length);
    result.trimmed_front += config_.umi_length;
    result.final_length = record.sequence.size();
    return true;
}
