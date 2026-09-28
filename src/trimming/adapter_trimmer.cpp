#include "trimming/adapter_trimmer.h"
#include "fastq_reader.h"

#include <algorithm>
#include <cctype>
#include <cstdint>
#include <limits>
#include <optional>
#include <stdexcept>
#include <utility>

namespace {

struct AdapterMatch {
    std::size_t position = 0;
};

struct OverlapMatch {
    std::size_t length = 0;
    std::size_t mismatches = 0;
};

void requireSynchronizedRead(const FastqRecord& record) {
    if (record.sequence.size() != record.quality.size()) {
        throw std::invalid_argument(
            "Cannot trim FASTQ record with different sequence and quality lengths");
    }
}

std::size_t removeTail(FastqRecord& record, std::size_t count) {
    const std::size_t removed = std::min(count, record.sequence.size());
    const std::size_t kept = record.sequence.size() - removed;
    record.sequence.resize(kept);
    record.quality.resize(kept);
    return removed;
}

char normalizedBase(char base) {
    return static_cast<char>(std::toupper(static_cast<unsigned char>(base)));
}

bool matchesAdapter(const FastqRecord& record,
                    std::size_t position,
                    const std::string& adapter,
                    std::size_t overlapLength) {
    const std::size_t allowedMismatches = overlapLength / 8;
    std::size_t mismatches = 0;
    for (std::size_t index = 0; index < overlapLength; ++index) {
        if (normalizedBase(record.sequence[position + index])
            != normalizedBase(adapter[index])) {
            ++mismatches;
            if (mismatches > allowedMismatches) {
                return false;
            }
        }
    }
    return true;
}

std::optional<AdapterMatch> findAdapterMatch(const FastqRecord& record,
                                             const TrimConfig& config) {
    const std::string& adapter = config.adapter_sequence;
    if (adapter.empty() || adapter.size() < config.min_adapter_match
        || record.sequence.size() < config.min_adapter_match) {
        return std::nullopt;
    }

    for (std::size_t position = 0;
         position + config.min_adapter_match <= record.sequence.size();
         ++position) {
        const std::size_t overlapLength = std::min(
            adapter.size(), record.sequence.size() - position);
        if (overlapLength >= config.min_adapter_match
            && matchesAdapter(record, position, adapter, overlapLength)) {
            return AdapterMatch{position};
        }
    }
    return std::nullopt;
}

char complement(char base) {
    switch (normalizedBase(base)) {
        case 'A': return 'T';
        case 'C': return 'G';
        case 'G': return 'C';
        case 'T': return 'A';
        default: return 'N';
    }
}

std::string reverseComplement(const std::string& sequence) {
    std::string result;
    result.reserve(sequence.size());
    for (auto iterator = sequence.rbegin(); iterator != sequence.rend(); ++iterator) {
        result.push_back(complement(*iterator));
    }
    return result;
}

std::size_t allowedOverlapMismatches(std::size_t overlapLength,
                                     const TrimConfig& config) {
    const std::uint64_t percentLimit =
        static_cast<std::uint64_t>(config.overlap_diff_percent_limit);
    const std::uint64_t overlap = static_cast<std::uint64_t>(overlapLength);
    const std::uint64_t allowedByPercent =
        percentLimit > std::numeric_limits<std::uint64_t>::max() / overlap
            ? std::numeric_limits<std::uint64_t>::max()
            : percentLimit * overlap / 100;
    return static_cast<std::size_t>(std::min<std::uint64_t>(
        config.overlap_diff_limit, allowedByPercent));
}

std::optional<OverlapMatch> findBestOverlap(const FastqRecord& r1,
                                            const FastqRecord& r2,
                                            const TrimConfig& config) {
    const std::size_t maximumOverlap =
        std::min(r1.sequence.size(), r2.sequence.size());
    if (maximumOverlap < config.min_overlap) {
        return std::nullopt;
    }

    const std::string reversedR2 = reverseComplement(r2.sequence);
    for (std::size_t overlapLength = maximumOverlap;
         overlapLength >= config.min_overlap;
         --overlapLength) {
        const std::size_t allowedMismatches =
            allowedOverlapMismatches(overlapLength, config);
        std::size_t mismatches = 0;
        const std::size_t r2Start = reversedR2.size() - overlapLength;
        for (std::size_t index = 0; index < overlapLength; ++index) {
            if (normalizedBase(r1.sequence[index])
                != reversedR2[r2Start + index]) {
                ++mismatches;
                if (mismatches > allowedMismatches) {
                    break;
                }
            }
        }
        if (mismatches <= allowedMismatches) {
            return OverlapMatch{overlapLength, mismatches};
        }
    }
    return std::nullopt;
}

void initializeResult(const FastqRecord& record, TrimResult& result) {
    if (result.original_length == 0 && result.final_length == 0
        && result.trimmed_front == 0 && result.trimmed_tail == 0) {
        result.original_length = record.sequence.size();
    }
}

void applyAdapterTrim(FastqRecord& record,
                      std::size_t position,
                      TrimResult& result) {
    const std::size_t removed = removeTail(record, record.sequence.size() - position);
    if (removed == 0) {
        return;
    }
    result.adapter_found = true;
    result.adapter_position = position;
    result.trimmed_tail += removed;
    result.final_length = record.sequence.size();
}

}  // namespace

AdapterTrimmer::AdapterTrimmer(TrimConfig config)
    : config_(std::move(config))
{
}

const TrimConfig& AdapterTrimmer::getConfig() const noexcept {
    return config_;
}

TrimResult AdapterTrimmer::trim(FastqRecord& record) const {
    TrimResult result;
    result.original_length = record.sequence.size();
    result.final_length = record.sequence.size();
    trim(record, result);
    return result;
}

void AdapterTrimmer::trim(FastqRecord& record, TrimResult& result) const {
    requireSynchronizedRead(record);
    initializeResult(record, result);
    result.final_length = record.sequence.size();

    if (!config_.enabled || !config_.adapter_trimming) {
        return;
    }
    if (config_.min_adapter_match == 0) {
        throw std::invalid_argument("Adapter trimming requires min_adapter_match greater than zero");
    }

    const std::optional<AdapterMatch> match = findAdapterMatch(record, config_);
    if (match.has_value()) {
        applyAdapterTrim(record, match->position, result);
    }
}

void AdapterTrimmer::trimPair(FastqRecord& r1,
                              FastqRecord& r2,
                              TrimResult& r1Result,
                              TrimResult& r2Result) const {
    requireSynchronizedRead(r1);
    requireSynchronizedRead(r2);
    initializeResult(r1, r1Result);
    initializeResult(r2, r2Result);
    r1Result.final_length = r1.sequence.size();
    r2Result.final_length = r2.sequence.size();

    if (!config_.enabled || !config_.adapter_trimming) {
        return;
    }
    if (config_.min_overlap == 0) {
        throw std::invalid_argument("Paired-end adapter trimming requires min_overlap greater than zero");
    }

    const std::optional<OverlapMatch> overlap = findBestOverlap(r1, r2, config_);
    if (overlap.has_value()) {
        applyAdapterTrim(r1, overlap->length, r1Result);
        applyAdapterTrim(r2, overlap->length, r2Result);
        return;
    }

    trim(r1, r1Result);
    trim(r2, r2Result);
}
