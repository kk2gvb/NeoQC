#include "quality_analyzer.h"

#include <cmath>
#include <cstdint>
#include <iostream>
#include <string>

namespace {

std::string uniqueSequence(uint64_t value) {
    constexpr char bases[] = {'A', 'C', 'G', 'T', 'N'};
    std::string sequence(DUPLICATION_PREFIX_LENGTH, 'A');
    for (std::size_t position = 0; position < sequence.size() && value > 0; ++position) {
        sequence[position] = bases[value % 5];
        value /= 5;
    }
    return sequence;
}

bool closeEnough(double left, double right) {
    return std::abs(left - right) < 1e-9;
}

}  // namespace

int main() {
    constexpr uint64_t uniqueCount = 100001;
    QualityAnalyzer analyzer;
    FastqRecord record;
    record.quality.assign(DUPLICATION_PREFIX_LENGTH, 'I');

    for (uint64_t index = 0; index < uniqueCount; ++index) {
        record.sequence = uniqueSequence(index);
        analyzer.processRecord(record);
    }
    record.sequence = uniqueSequence(0);
    analyzer.processRecord(record);

    const DuplicationStats stats = analyzer.getDuplicationStats();
    const uint64_t totalReads = uniqueCount + 1;
    const double expectedRemaining =
        100.0 * static_cast<double>(uniqueCount) / static_cast<double>(totalReads);

    if (stats.totalReads != totalReads || stats.uniqueSequences != uniqueCount) {
        std::cerr << "Unique-prefix map stopped before 100001 keys\n";
        return 1;
    }
    if (!closeEnough(stats.deduplicatedRemainingPercent, expectedRemaining)) {
        std::cerr << "Incorrect exact deduplicated percentage\n";
        return 1;
    }
    if (stats.levels.size() != 16) {
        std::cerr << "Expected 16 stable duplication bins\n";
        return 1;
    }

    const double expectedLevelOneTotal =
        100.0 * static_cast<double>(uniqueCount - 1) / static_cast<double>(totalReads);
    const double expectedLevelTwoTotal = 200.0 / static_cast<double>(totalReads);
    if (!closeEnough(stats.levels[0].totalSequencesPercent, expectedLevelOneTotal)
        || !closeEnough(stats.levels[1].totalSequencesPercent, expectedLevelTwoTotal)) {
        std::cerr << "Incorrect exact duplication-level distribution\n";
        return 1;
    }

    return 0;
}
