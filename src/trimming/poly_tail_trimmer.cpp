#include "trimming/poly_tail_trimmer.h"
#include "fastq_reader.h"

#include <cctype>
#include <stdexcept>
#include <utility>

namespace {

char normalizedBase(char base) {
    return static_cast<char>(std::toupper(static_cast<unsigned char>(base)));
}

void requireSynchronizedRead(const FastqRecord& record) {
    if (record.sequence.size() != record.quality.size()) {
        throw std::invalid_argument(
            "Cannot trim FASTQ record with different sequence and quality lengths");
    }
}

std::size_t terminalRunLength(const std::string& sequence, char base) {
    std::size_t length = 0;
    for (auto iterator = sequence.rbegin(); iterator != sequence.rend(); ++iterator) {
        if (normalizedBase(*iterator) != base) {
            break;
        }
        ++length;
    }
    return length;
}

std::size_t trimTerminalRun(FastqRecord& record,
                            char base,
                            std::size_t minimumLength) {
    const std::size_t length = terminalRunLength(record.sequence, base);
    if (length < minimumLength) {
        return 0;
    }
    record.sequence.resize(record.sequence.size() - length);
    record.quality.resize(record.quality.size() - length);
    return length;
}

bool isPolyXBase(char base) {
    return base == 'A' || base == 'C' || base == 'G' || base == 'T';
}

}  // namespace

PolyTailTrimmer::PolyTailTrimmer(TrimConfig config)
    : config_(std::move(config))
{
}

const TrimConfig& PolyTailTrimmer::getConfig() const noexcept {
    return config_;
}

void PolyTailTrimmer::trim(FastqRecord& record, TrimResult& result) const {
    requireSynchronizedRead(record);
    result.final_length = record.sequence.size();

    if (!config_.enabled) {
        return;
    }

    if (config_.trim_poly_g) {
        if (config_.poly_g_min_length == 0) {
            throw std::invalid_argument(
                "PolyG trimming requires poly_g_min_length greater than zero");
        }
        const std::size_t removed = trimTerminalRun(
            record, 'G', config_.poly_g_min_length);
        if (removed > 0) {
            result.trimmed_tail += removed;
            result.polyG_trimmed = true;
        }
    }

    if (config_.trim_poly_x) {
        if (config_.poly_x_min_length == 0) {
            throw std::invalid_argument(
                "PolyX trimming requires poly_x_min_length greater than zero");
        }
        if (record.sequence.empty()) {
            result.final_length = 0;
            return;
        }
        const char base = normalizedBase(record.sequence.back());
        if (isPolyXBase(base)) {
            const std::size_t removed = trimTerminalRun(
                record, base, config_.poly_x_min_length);
            if (removed > 0) {
                result.trimmed_tail += removed;
                result.polyX_trimmed = true;
            }
        }
    }

    result.final_length = record.sequence.size();
}
