#include "trimming/quality_trimmer.h"
#include "fastq_reader.h"

#include <algorithm>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <utility>

namespace {

std::size_t removeFront(FastqRecord& record, std::size_t count) {
    const std::size_t removed = std::min(count, record.sequence.size());
    record.sequence.erase(0, removed);
    record.quality.erase(0, removed);
    return removed;
}

std::size_t removeTail(FastqRecord& record, std::size_t count) {
    const std::size_t removed = std::min(count, record.sequence.size());
    const std::size_t kept = record.sequence.size() - removed;
    record.sequence.resize(kept);
    record.quality.resize(kept);
    return removed;
}

bool isLowQuality(char qualityChar, std::size_t threshold) {
    const int phred =
        static_cast<int>(static_cast<unsigned char>(qualityChar)) - 33;
    return phred < 0 || static_cast<std::size_t>(phred) < threshold;
}

std::uint64_t phredScore(char qualityChar) {
    const int phred =
        static_cast<int>(static_cast<unsigned char>(qualityChar)) - 33;
    return phred < 0 ? 0 : static_cast<std::uint64_t>(phred);
}

bool windowBelowThreshold(std::uint64_t sum,
                          std::size_t windowSize,
                          std::size_t threshold) {
    const std::uint64_t window = static_cast<std::uint64_t>(windowSize);
    const std::uint64_t qualityThreshold =
        static_cast<std::uint64_t>(threshold);
    if (qualityThreshold > std::numeric_limits<std::uint64_t>::max() / window) {
        return true;
    }
    return sum < qualityThreshold * window;
}

void requireSynchronizedRead(const FastqRecord& record) {
    if (record.sequence.size() != record.quality.size()) {
        throw std::invalid_argument(
            "Cannot trim FASTQ record with different sequence and quality lengths");
    }
}

}  // namespace

QualityTrimmer::QualityTrimmer(TrimConfig config)
    : config_(std::move(config))
{
}

const TrimConfig& QualityTrimmer::getConfig() const noexcept {
    return config_;
}

TrimResult QualityTrimmer::trim(FastqRecord& record) const {
    requireSynchronizedRead(record);

    TrimResult result;
    result.original_length = record.sequence.size();

    if (!config_.enabled) {
        result.final_length = record.sequence.size();
        return result;
    }

    result.trimmed_front = removeFront(record, config_.trim_front);
    result.trimmed_tail = removeTail(record, config_.trim_tail);

    if (config_.cut_front) {
        std::size_t removed = 0;
        while (removed < record.quality.size()
               && isLowQuality(record.quality[removed],
                               config_.quality_threshold)) {
            ++removed;
        }
        if (removed > 0) {
            result.trimmed_front += removeFront(record, removed);
            result.quality_trimmed = true;
        }
    }

    if (config_.cut_tail) {
        std::size_t removed = 0;
        const std::size_t length = record.quality.size();
        while (removed < length
               && isLowQuality(record.quality[length - 1 - removed],
                               config_.quality_threshold)) {
            ++removed;
        }
        if (removed > 0) {
            result.trimmed_tail += removeTail(record, removed);
            result.quality_trimmed = true;
        }
    }

    if (config_.cut_right) {
        if (config_.window_size == 0) {
            throw std::invalid_argument(
                "Sliding-window trimming requires window_size greater than zero");
        }

        const std::size_t length = record.quality.size();
        if (config_.window_size <= length) {
            std::uint64_t rollingSum = 0;
            for (std::size_t index = 0; index < config_.window_size; ++index) {
                rollingSum += phredScore(record.quality[index]);
            }

            std::size_t trimPosition = length;
            for (std::size_t start = 0;; ++start) {
                if (windowBelowThreshold(rollingSum,
                                         config_.window_size,
                                         config_.quality_threshold)) {
                    trimPosition = start;
                    break;
                }
                if (start + config_.window_size == length) {
                    break;
                }
                rollingSum -= phredScore(record.quality[start]);
                rollingSum += phredScore(
                    record.quality[start + config_.window_size]);
            }

            if (trimPosition < length) {
                result.trimmed_tail += removeTail(record, length - trimPosition);
                result.quality_trimmed = true;
            }
        }
    }

    result.final_length = record.sequence.size();
    return result;
}
