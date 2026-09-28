#include "trimming/pe_processor.h"

#include <algorithm>
#include <cctype>
#include <stdexcept>
#include <string>

namespace {

bool canonical(char base) {
    switch (PeOverlapAnalyzer::normalizedBase(base)) {
        case 'A':
        case 'C':
        case 'G':
        case 'T': return true;
        default: return false;
    }
}

std::string pairIdentity(const std::string& header) {
    const std::size_t whitespace = header.find_first_of(" \t");
    std::string identity = header.substr(0, whitespace);
    if (!identity.empty() && identity.front() == '@') {
        identity.erase(identity.begin());
    }
    if (identity.size() >= 2 && identity[identity.size() - 2] == '/'
        && (identity.back() == '1' || identity.back() == '2')) {
        identity.resize(identity.size() - 2);
    }
    return identity;
}

}  // namespace

void PeProcessor::correct(FastqRecord& r1,
                          FastqRecord& r2,
                          TrimResult& r1Result,
                          TrimResult& r2Result,
                          const PeOverlap& overlap) const {
    if (r1.sequence.size() != r1.quality.size()
        || r2.sequence.size() != r2.quality.size()) {
        throw std::invalid_argument(
            "Cannot correct PE reads with different sequence and quality lengths");
    }

    for (std::size_t offset = 0; offset < overlap.length; ++offset) {
        const std::size_t r1Index = overlap.r1_start + offset;
        const std::size_t reverseIndex = overlap.r2_reverse_start + offset;
        const std::size_t r2Index = PeOverlapAnalyzer::originalR2Index(
            r2.sequence.size(), reverseIndex);
        const char r1Base = PeOverlapAnalyzer::normalizedBase(r1.sequence[r1Index]);
        const char r2Base = PeOverlapAnalyzer::complement(r2.sequence[r2Index]);
        if (r1Base == r2Base || (r1Base == 'N' && r2Base == 'N')) {
            continue;
        }

        const char r1Quality = r1.quality[r1Index];
        const char r2Quality = r2.quality[r2Index];
        if (r1Base == 'N' && canonical(r2Base) && r2Quality > r1Quality) {
            r1.sequence[r1Index] = r2Base;
            r1.quality[r1Index] = r2Quality;
            ++r1Result.corrected_bases;
        } else if (r2Base == 'N' && canonical(r1Base) && r1Quality > r2Quality) {
            r2.sequence[r2Index] = PeOverlapAnalyzer::complement(r1Base);
            r2.quality[r2Index] = r1Quality;
            ++r2Result.corrected_bases;
        } else if (canonical(r1Base) && canonical(r2Base)) {
            if (r1Quality > r2Quality) {
                r2.sequence[r2Index] = PeOverlapAnalyzer::complement(r1Base);
                r2.quality[r2Index] = r1Quality;
                ++r2Result.corrected_bases;
            } else if (r2Quality > r1Quality) {
                r1.sequence[r1Index] = r2Base;
                r1.quality[r1Index] = r2Quality;
                ++r1Result.corrected_bases;
            }
        }
    }
}

std::optional<FastqRecord> PeProcessor::merge(const FastqRecord& r1,
                                              const FastqRecord& r2,
                                              const PeOverlap& overlap) const {
    if (r1.sequence.size() != r1.quality.size()
        || r2.sequence.size() != r2.quality.size()) {
        throw std::invalid_argument(
            "Cannot merge PE reads with different sequence and quality lengths");
    }

    const std::string r2Sequence = PeOverlapAnalyzer::reverseComplement(r2.sequence);
    const std::string r2Quality = PeOverlapAnalyzer::reverseQuality(r2.quality);
    const std::ptrdiff_t shift = static_cast<std::ptrdiff_t>(overlap.r1_start)
        - static_cast<std::ptrdiff_t>(overlap.r2_reverse_start);
    const std::ptrdiff_t begin = overlap.isReadThrough()
        ? 0 : std::min<std::ptrdiff_t>(0, shift);
    const std::ptrdiff_t end = overlap.isReadThrough()
        ? static_cast<std::ptrdiff_t>(overlap.length)
        : std::max<std::ptrdiff_t>(
            static_cast<std::ptrdiff_t>(r1.sequence.size()),
            shift + static_cast<std::ptrdiff_t>(r2Sequence.size()));

    FastqRecord merged;
    merged.header = "@" + pairIdentity(r1.header) + " merged";
    merged.separator = "+";
    merged.recordNumber = r1.recordNumber;
    merged.sequence.reserve(static_cast<std::size_t>(end - begin));
    merged.quality.reserve(static_cast<std::size_t>(end - begin));

    for (std::ptrdiff_t coordinate = begin; coordinate < end; ++coordinate) {
        const bool hasR1 = coordinate >= 0
            && coordinate < static_cast<std::ptrdiff_t>(r1.sequence.size());
        const std::ptrdiff_t r2Index = coordinate - shift;
        const bool hasR2 = r2Index >= 0
            && r2Index < static_cast<std::ptrdiff_t>(r2Sequence.size());
        if (hasR1 && hasR2) {
            const std::size_t i1 = static_cast<std::size_t>(coordinate);
            const std::size_t i2 = static_cast<std::size_t>(r2Index);
            const char base1 = PeOverlapAnalyzer::normalizedBase(r1.sequence[i1]);
            const char base2 = PeOverlapAnalyzer::normalizedBase(r2Sequence[i2]);
            const char quality1 = r1.quality[i1];
            const char quality2 = r2Quality[i2];
            if (base1 == base2) {
                merged.sequence.push_back(base1);
                merged.quality.push_back(std::max(quality1, quality2));
            } else if (quality1 > quality2) {
                merged.sequence.push_back(base1);
                merged.quality.push_back(quality1);
            } else if (quality2 > quality1) {
                merged.sequence.push_back(base2);
                merged.quality.push_back(quality2);
            } else {
                merged.sequence.push_back('N');
                merged.quality.push_back(std::min(quality1, quality2));
            }
        } else if (hasR1) {
            const std::size_t index = static_cast<std::size_t>(coordinate);
            merged.sequence.push_back(r1.sequence[index]);
            merged.quality.push_back(r1.quality[index]);
        } else if (hasR2) {
            const std::size_t index = static_cast<std::size_t>(r2Index);
            merged.sequence.push_back(r2Sequence[index]);
            merged.quality.push_back(r2Quality[index]);
        }
    }
    return merged;
}
