#include "trimming/pe_overlap.h"
#include "fastq_reader.h"

#include <algorithm>
#include <cctype>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <utility>

namespace {

std::size_t allowedMismatches(std::size_t overlapLength,
                              const TrimConfig& config) {
    const std::uint64_t percent = config.overlap_diff_percent_limit;
    const std::uint64_t length = overlapLength;
    const std::uint64_t byPercent =
        percent > std::numeric_limits<std::uint64_t>::max() / length
            ? std::numeric_limits<std::uint64_t>::max()
            : percent * length / 100;
    return static_cast<std::size_t>(std::min<std::uint64_t>(
        config.overlap_diff_limit, byPercent));
}

bool better(const PeOverlap& candidate, const PeOverlap& current) {
    if (candidate.length != current.length) {
        return candidate.length > current.length;
    }
    if (candidate.mismatches != current.mismatches) {
        return candidate.mismatches < current.mismatches;
    }
    if (candidate.r1_start != current.r1_start) {
        return candidate.r1_start < current.r1_start;
    }
    return candidate.r2_reverse_start < current.r2_reverse_start;
}

}  // namespace

PeOverlapAnalyzer::PeOverlapAnalyzer(TrimConfig config)
    : config_(std::move(config))
{
}

char PeOverlapAnalyzer::normalizedBase(char base) noexcept {
    return static_cast<char>(std::toupper(static_cast<unsigned char>(base)));
}

char PeOverlapAnalyzer::complement(char base) noexcept {
    switch (normalizedBase(base)) {
        case 'A': return 'T';
        case 'C': return 'G';
        case 'G': return 'C';
        case 'T': return 'A';
        default: return 'N';
    }
}

std::string PeOverlapAnalyzer::reverseComplement(const std::string& sequence) {
    std::string result;
    result.reserve(sequence.size());
    for (auto iterator = sequence.rbegin(); iterator != sequence.rend(); ++iterator) {
        result.push_back(complement(*iterator));
    }
    return result;
}

std::string PeOverlapAnalyzer::reverseQuality(const std::string& quality) {
    return std::string(quality.rbegin(), quality.rend());
}

std::size_t PeOverlapAnalyzer::originalR2Index(std::size_t r2Length,
                                               std::size_t reverseIndex) {
    if (reverseIndex >= r2Length) {
        throw std::out_of_range("Reverse R2 overlap coordinate is out of range");
    }
    return r2Length - reverseIndex - 1;
}

std::optional<PeOverlap> PeOverlapAnalyzer::find(const FastqRecord& r1,
                                                 const FastqRecord& r2) const {
    if (r1.sequence.size() != r1.quality.size()
        || r2.sequence.size() != r2.quality.size()) {
        throw std::invalid_argument(
            "Cannot analyze PE overlap with different sequence and quality lengths");
    }
    if (config_.min_overlap == 0) {
        throw std::invalid_argument("PE overlap requires min_overlap greater than zero");
    }
    if (r1.sequence.size() < config_.min_overlap
        || r2.sequence.size() < config_.min_overlap) {
        return std::nullopt;
    }

    const std::string r2Reverse = reverseComplement(r2.sequence);
    const std::ptrdiff_t minimumShift =
        -static_cast<std::ptrdiff_t>(r2Reverse.size() - config_.min_overlap);
    const std::ptrdiff_t maximumShift =
        static_cast<std::ptrdiff_t>(r1.sequence.size() - config_.min_overlap);
    std::optional<PeOverlap> best;

    for (std::ptrdiff_t shift = minimumShift; shift <= maximumShift; ++shift) {
        const std::size_t r1Start = shift > 0
            ? static_cast<std::size_t>(shift) : 0;
        const std::size_t r2Start = shift < 0
            ? static_cast<std::size_t>(-shift) : 0;
        const std::size_t length = std::min(
            r1.sequence.size() - r1Start, r2Reverse.size() - r2Start);
        if (length < config_.min_overlap) {
            continue;
        }

        const std::size_t allowed = allowedMismatches(length, config_);
        std::size_t mismatches = 0;
        for (std::size_t index = 0; index < length; ++index) {
            if (normalizedBase(r1.sequence[r1Start + index])
                != normalizedBase(r2Reverse[r2Start + index])) {
                ++mismatches;
                if (mismatches > allowed) {
                    break;
                }
            }
        }
        if (mismatches > allowed) {
            continue;
        }
        const PeOverlap candidate{r1Start, r2Start, length, mismatches};
        if (!best.has_value() || better(candidate, *best)) {
            best = candidate;
        }
    }
    return best;
}
