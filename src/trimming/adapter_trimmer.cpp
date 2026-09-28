#include "trimming/adapter_trimmer.h"
#include "fastq_reader.h"
#include "trimming/pe_overlap.h"

#include <algorithm>
#include <cctype>
#include <optional>
#include <stdexcept>
#include <utility>

namespace {

struct AdapterMatch {
    std::size_t position = 0;
    std::size_t mismatches = 0;
    std::size_t length = 0;
    std::size_t order = 0;
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

std::optional<std::size_t> adapterMismatches(const FastqRecord& record,
                                            std::size_t position,
                                            const std::string& adapter,
                                            std::size_t overlapLength) {
    const std::size_t allowedMismatches = overlapLength / 8;
    std::size_t mismatches = 0;
    for (std::size_t index = 0; index < overlapLength; ++index) {
        if (PeOverlapAnalyzer::normalizedBase(record.sequence[position + index])
            != PeOverlapAnalyzer::normalizedBase(adapter[index])) {
            ++mismatches;
            if (mismatches > allowedMismatches) {
                return std::nullopt;
            }
        }
    }
    return mismatches;
}

std::optional<AdapterMatch> findAdapterMatch(const FastqRecord& record,
                                             const TrimConfig& config) {
    if (record.sequence.size() < config.min_adapter_match) {
        return std::nullopt;
    }

    std::vector<std::string> adapters;
    if (!config.adapter_sequence.empty()) {
        adapters.push_back(config.adapter_sequence);
    }
    for (const TrimAdapter& adapter : config.fasta_adapters) {
        adapters.push_back(adapter.sequence);
    }

    std::optional<AdapterMatch> best;
    for (std::size_t order = 0; order < adapters.size(); ++order) {
        const std::string& adapter = adapters[order];
        if (adapter.size() < config.min_adapter_match) {
            continue;
        }
        for (std::size_t position = 0;
             position + config.min_adapter_match <= record.sequence.size();
             ++position) {
            const std::size_t length = std::min(
                adapter.size(), record.sequence.size() - position);
            const auto mismatches = adapterMismatches(
                record, position, adapter, length);
            if (!mismatches.has_value()) {
                continue;
            }
            const AdapterMatch candidate{position, *mismatches, length, order};
            if (!best.has_value()
                || candidate.position < best->position
                || (candidate.position == best->position
                    && candidate.length > best->length)
                || (candidate.position == best->position
                    && candidate.length == best->length
                    && candidate.mismatches < best->mismatches)
                || (candidate.position == best->position
                    && candidate.length == best->length
                    && candidate.mismatches == best->mismatches
                    && candidate.order < best->order)) {
                best = candidate;
            }
        }
    }
    return best;
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
    const std::optional<PeOverlap> overlap =
        config_.enabled && config_.adapter_trimming
            ? PeOverlapAnalyzer(config_).find(r1, r2)
            : std::nullopt;
    trimPair(r1, r2, r1Result, r2Result, overlap);
}

void AdapterTrimmer::trimPair(FastqRecord& r1,
                              FastqRecord& r2,
                              TrimResult& r1Result,
                              TrimResult& r2Result,
                              const std::optional<PeOverlap>& overlap) const {
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

    if (overlap.has_value() && overlap->isReadThrough()) {
        applyAdapterTrim(r1, overlap->length, r1Result);
        applyAdapterTrim(r2, overlap->length, r2Result);
        return;
    }

    trim(r1, r1Result);
    trim(r2, r2Result);
}
