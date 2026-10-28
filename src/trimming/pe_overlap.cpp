#include "trimming/pe_overlap.h"
#include "fastq_reader.h"

#include <algorithm>
#include <array>
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

const std::array<char, 256>& normalizedBases() {
    static const std::array<char, 256> bases = [] {
        std::array<char, 256> result{};
        for (std::size_t index = 0; index < result.size(); ++index) {
            result[index] = static_cast<char>(std::toupper(
                static_cast<unsigned char>(index)));
        }
        return result;
    }();
    return bases;
}

const std::array<char, 256>& complementedBases() {
    static const std::array<char, 256> bases = [] {
        std::array<char, 256> result{};
        const auto& normalized = normalizedBases();
        for (std::size_t index = 0; index < result.size(); ++index) {
            switch (normalized[index]) {
                case 'A': result[index] = 'T'; break;
                case 'C': result[index] = 'G'; break;
                case 'G': result[index] = 'C'; break;
                case 'T': result[index] = 'A'; break;
                default: result[index] = 'N'; break;
            }
        }
        return result;
    }();
    return bases;
}

void validateInputs(const FastqRecord& r1,
                    const FastqRecord& r2,
                    const TrimConfig& config) {
    if (r1.sequence.size() != r1.quality.size()
        || r2.sequence.size() != r2.quality.size()) {
        throw std::invalid_argument(
            "Cannot analyze PE overlap with different sequence and quality lengths");
    }
    if (config.min_overlap == 0) {
        throw std::invalid_argument("PE overlap requires min_overlap greater than zero");
    }
}

// With at most K mismatches, partitioning an alignment into K + 1 blocks
// guarantees that at least one block is exact. This filter therefore cannot
// reject an alignment accepted by the exhaustive mismatch comparison.
bool hasExactPigeonholeBlock(const FastqRecord& r1,
                             const FastqRecord& r2,
                             std::size_t r1Start,
                             std::size_t r2ReverseStart,
                             std::size_t length,
                             std::size_t allowedMismatches) {
    if (allowedMismatches == 0 || allowedMismatches >= length) {
        return true;
    }

    const std::size_t blockCount = allowedMismatches + 1;
    const std::size_t baseBlockLength = length / blockCount;
    const std::size_t largerBlockCount = length % blockCount;
    const auto& normalized = normalizedBases();
    const auto& complemented = complementedBases();
    std::size_t offset = 0;
    for (std::size_t block = 0; block < blockCount; ++block) {
        const std::size_t blockLength = baseBlockLength
            + (block < largerBlockCount ? 1 : 0);
        std::size_t index = 0;
        for (; index < blockLength; ++index) {
            const std::size_t r2Index = r2.sequence.size() - 1
                - (r2ReverseStart + offset + index);
            if (normalized[static_cast<unsigned char>(
                    r1.sequence[r1Start + offset + index])]
                != complemented[static_cast<unsigned char>(r2.sequence[r2Index])]) {
                break;
            }
        }
        if (index == blockLength) {
            return true;
        }
        offset += blockLength;
    }
    return false;
}

std::optional<PeOverlap> findExhaustive(const FastqRecord& r1,
                                        const FastqRecord& r2,
                                        const TrimConfig& config) {
    validateInputs(r1, r2, config);
    if (r1.sequence.size() < config.min_overlap
        || r2.sequence.size() < config.min_overlap) {
        return std::nullopt;
    }

    const std::string r2Reverse = PeOverlapAnalyzer::reverseComplement(r2.sequence);
    const std::ptrdiff_t minimumShift =
        -static_cast<std::ptrdiff_t>(r2Reverse.size() - config.min_overlap);
    const std::ptrdiff_t maximumShift =
        static_cast<std::ptrdiff_t>(r1.sequence.size() - config.min_overlap);
    std::optional<PeOverlap> best;

    for (std::ptrdiff_t shift = minimumShift; shift <= maximumShift; ++shift) {
        const std::size_t r1Start = shift > 0
            ? static_cast<std::size_t>(shift) : 0;
        const std::size_t r2Start = shift < 0
            ? static_cast<std::size_t>(-shift) : 0;
        const std::size_t length = std::min(
            r1.sequence.size() - r1Start, r2Reverse.size() - r2Start);
        const std::size_t allowed = allowedMismatches(length, config);
        std::size_t mismatches = 0;
        for (std::size_t index = 0; index < length; ++index) {
            if (PeOverlapAnalyzer::normalizedBase(r1.sequence[r1Start + index])
                != PeOverlapAnalyzer::normalizedBase(r2Reverse[r2Start + index])) {
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

}  // namespace

PeOverlapAnalyzer::PeOverlapAnalyzer(TrimConfig config)
    : config_(std::move(config))
{
}

char PeOverlapAnalyzer::normalizedBase(char base) noexcept {
    return normalizedBases()[static_cast<unsigned char>(base)];
}

char PeOverlapAnalyzer::complement(char base) noexcept {
    return complementedBases()[static_cast<unsigned char>(base)];
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
    validateInputs(r1, r2, config_);
    if (r1.sequence.size() < config_.min_overlap
        || r2.sequence.size() < config_.min_overlap) {
        return std::nullopt;
    }

    const std::ptrdiff_t minimumShift =
        -static_cast<std::ptrdiff_t>(r2.sequence.size() - config_.min_overlap);
    const std::ptrdiff_t maximumShift =
        static_cast<std::ptrdiff_t>(r1.sequence.size() - config_.min_overlap);
    const auto& normalized = normalizedBases();
    const auto& complemented = complementedBases();
    std::optional<PeOverlap> best;

    for (std::ptrdiff_t shift = minimumShift; shift <= maximumShift; ++shift) {
        const std::size_t r1Start = shift > 0
            ? static_cast<std::size_t>(shift) : 0;
        const std::size_t r2Start = shift < 0
            ? static_cast<std::size_t>(-shift) : 0;
        const std::size_t length = std::min(
            r1.sequence.size() - r1Start, r2.sequence.size() - r2Start);

        const std::size_t allowed = allowedMismatches(length, config_);
        if (!hasExactPigeonholeBlock(
                r1, r2, r1Start, r2Start, length, allowed)) {
            continue;
        }
        std::size_t mismatches = 0;
        for (std::size_t index = 0; index < length; ++index) {
            const std::size_t r2Index = r2.sequence.size() - 1 - (r2Start + index);
            if (normalized[static_cast<unsigned char>(r1.sequence[r1Start + index])]
                != complemented[static_cast<unsigned char>(r2.sequence[r2Index])]) {
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

std::optional<PeOverlap> PeOverlapAnalyzer::findReference(
    const FastqRecord& r1,
    const FastqRecord& r2) const {
    return findExhaustive(r1, r2, config_);
}
