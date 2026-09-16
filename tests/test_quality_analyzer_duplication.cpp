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

        // Exactly 0.1%: 1 occurrence out of 1000 reads.
    {
        QualityAnalyzer thresholdAnalyzer;
        FastqRecord thresholdRecord;
        thresholdRecord.quality.assign(DUPLICATION_PREFIX_LENGTH, 'I');

        for (int i = 0; i < 999; ++i) {
            thresholdRecord.sequence = "TTTTTTTTTTTT";
            thresholdAnalyzer.processRecord(thresholdRecord);
        }

        thresholdRecord.sequence = "ACGTACGTACGT";
        thresholdAnalyzer.processRecord(thresholdRecord);

        const DuplicationStats thresholdStats =
            thresholdAnalyzer.getDuplicationStats();

        bool found = false;
        for (const auto& sequence :
             thresholdStats.overrepresentedSequences) {
            if (sequence.sequence == "ACGTACGTACGT") {
                found = true;

                if (!closeEnough(sequence.percent, 0.1)) {
                    std::cerr << "Incorrect 0.1% boundary percentage\n";
                    return 1;
                }

                if (sequence.count != 1) {
                    std::cerr << "Incorrect 0.1% boundary count\n";
                    return 1;
                }
            }
        }

        if (!found) {
            std::cerr << "Sequence at exactly 0.1% must be reported\n";
            return 1;
        }
    }

    // Just below 0.1%: 1 occurrence out of 1001 reads.
    {
        QualityAnalyzer belowThresholdAnalyzer;
        FastqRecord belowThresholdRecord;
        belowThresholdRecord.quality.assign(DUPLICATION_PREFIX_LENGTH, 'I');

        for (int i = 0; i < 1000; ++i) {
            belowThresholdRecord.sequence = "TTTTTTTTTTTT";
            belowThresholdAnalyzer.processRecord(belowThresholdRecord);
        }

        belowThresholdRecord.sequence = "ACGTACGTACGT";
        belowThresholdAnalyzer.processRecord(belowThresholdRecord);

        const DuplicationStats belowThresholdStats =
            belowThresholdAnalyzer.getDuplicationStats();

        for (const auto& sequence :
             belowThresholdStats.overrepresentedSequences) {
            if (sequence.sequence == "ACGTACGTACGT") {
                std::cerr << "Sequence below 0.1% must not be reported\n";
                return 1;
            }
        }
    }

    return 0;
}
