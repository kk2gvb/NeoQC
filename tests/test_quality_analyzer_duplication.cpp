#include "quality_analyzer.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <random>
#include <string>
#include <unordered_map>
#include <vector>

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

bool sameStats(const DuplicationStats& left, const DuplicationStats& right) {
    if (left.totalReads != right.totalReads
        || left.uniqueSequences != right.uniqueSequences
        || !closeEnough(left.deduplicatedRemainingPercent,
                        right.deduplicatedRemainingPercent)
        || left.levels.size() != right.levels.size()
        || left.overrepresentedSequences.size() != right.overrepresentedSequences.size()) {
        return false;
    }
    for (std::size_t index = 0; index < left.levels.size(); ++index) {
        if (left.levels[index].label != right.levels[index].label
            || !closeEnough(left.levels[index].totalSequencesPercent,
                            right.levels[index].totalSequencesPercent)
            || !closeEnough(left.levels[index].deduplicatedSequencesPercent,
                            right.levels[index].deduplicatedSequencesPercent)) {
            return false;
        }
    }
    for (std::size_t index = 0; index < left.overrepresentedSequences.size(); ++index) {
        const auto& a = left.overrepresentedSequences[index];
        const auto& b = right.overrepresentedSequences[index];
        if (a.sequence != b.sequence || a.count != b.count
            || !closeEnough(a.percent, b.percent)) {
            return false;
        }
    }
    return true;
}

bool testExactAggregationDifferential() {
    std::vector<std::string> sequences;
    sequences.insert(sequences.end(), 60, std::string(DUPLICATION_PREFIX_LENGTH, 'A'));
    sequences.insert(sequences.end(), 20,
                     std::string(DUPLICATION_PREFIX_LENGTH, 'C') + "suffix-one");
    sequences.insert(sequences.end(), 10,
                     std::string(DUPLICATION_PREFIX_LENGTH, 'C') + "suffix-two");
    sequences.insert(sequences.end(), 9, "ACGTNACGTN");
    sequences.insert(sequences.end(), 7, "acgtnacgtn");
    sequences.insert(sequences.end(), 5, "ACG");
    sequences.push_back(std::string(DUPLICATION_PREFIX_LENGTH - 1, 'G') + "A");
    sequences.push_back(std::string(DUPLICATION_PREFIX_LENGTH - 1, 'G') + "T");

    std::mt19937_64 random(0x4455504c49434154ULL);
    constexpr char bases[] = "ACGTNacgtn";
    std::uniform_int_distribution<std::size_t> length(1, 75);
    std::uniform_int_distribution<std::size_t> base(0, sizeof(bases) - 2);
    for (std::size_t index = 0; index < 20000; ++index) {
        std::string sequence(length(random), 'A');
        for (char& nucleotide : sequence) nucleotide = bases[base(random)];
        sequences.push_back(std::move(sequence));
    }

    std::vector<QualityAnalyzer> referenceWorkers;
    std::vector<QualityAnalyzer> optimizedWorkers;
    referenceWorkers.reserve(4);
    optimizedWorkers.reserve(4);
    for (int index = 0; index < 4; ++index) {
        referenceWorkers.emplace_back();
        optimizedWorkers.emplace_back();
    }

    FastqRecord record;
    for (std::size_t index = 0; index < sequences.size(); ++index) {
        record.sequence = sequences[index];
        record.quality.assign(record.sequence.size(), 'I');
        referenceWorkers[index % referenceWorkers.size()].processRecord(record);
        optimizedWorkers[index % optimizedWorkers.size()].processRecord(record);
    }

    std::unordered_map<DuplicationKey, uint64_t, DuplicationKeyHash> referenceCounts;
    QualityAnalyzer referenceSummary;
    for (const auto& worker : referenceWorkers) {
        referenceSummary.merge(worker);
        for (const DuplicationEntry& entry : worker.getDuplicationEntries()) {
            referenceCounts[entry.key] += entry.count;
        }
    }
    std::vector<DuplicationEntry> referenceEntries;
    referenceEntries.reserve(referenceCounts.size());
    for (const auto& [key, count] : referenceCounts) {
        referenceEntries.push_back({key, count});
    }

    QualityAnalyzer optimizedSummary;
    std::size_t expectedEntries = 0;
    for (const auto& worker : optimizedWorkers) {
        expectedEntries += worker.getDuplicationEntryCount();
    }
    optimizedSummary.reserveDuplicationEntries(expectedEntries);
    for (auto& worker : optimizedWorkers) {
        optimizedSummary.merge(worker);
        optimizedSummary.mergeDuplicationCounts(worker);
    }

    const DuplicationStats referenceStats =
        referenceSummary.getDuplicationStats(referenceEntries);
    const DuplicationStats optimizedStats = optimizedSummary.getDuplicationStats();
    if (!sameStats(referenceStats, optimizedStats)) {
        std::cerr << "Exact duplication aggregation differs from reference counts\n";
        return false;
    }
    if (optimizedStats.uniqueSequences != referenceCounts.size()) {
        std::cerr << "Exact duplication aggregation lost a unique prefix\n";
        return false;
    }
    return true;
}

bool testAdapterContentDifferential() {
    const std::vector<AdapterConfigEntry> adapterConfig = {
        {"adapter-a", "ACGTACGTACGTGGGG"},
        {"adapter-b", "TTTTCCCCAAAAGGGG"},
        {"adapter-c", "GGGGGGGGGGGGAAAA"},
        {"adapter-short", "ACGT"},
    };
    std::vector<std::string> sequences = {
        "CCCCCCCCCCCC",                         // no adapter
        "ACGTACGTACGT",                         // start
        "CCCCACGTACGTACGTGG",                   // middle
        "CCCCACGTACGTACGT",                     // end
        "ACGTACGTACGTACGTACGTACGT",             // multiple occurrences
        "TTTTCCCCAAAAGGGGACGTACGTACGT",         // multiple adapter types
        "ACGTACGTACG",                          // partial adapter
        "ACGTACGTACGA",                         // one-base difference
        "ACGTNCGTACGT",                         // N inside candidate
        "acgtacgtacgt",                         // lowercase remains case-sensitive
        "ACGT",                                 // short read
        std::string(200, 'A'),                   // long homopolymer
        "GGGGGGGGGGGGGGGGGGGGGGGG",             // repetitive sequence
    };

    std::mt19937_64 random(0x4144415054455251ULL);
    constexpr char bases[] = "ACGTNacgtn";
    std::uniform_int_distribution<std::size_t> length(1, 160);
    std::uniform_int_distribution<std::size_t> base(0, sizeof(bases) - 2);
    for (std::size_t index = 0; index < 20000; ++index) {
        std::string sequence(length(random), 'A');
        for (char& nucleotide : sequence) nucleotide = bases[base(random)];
        if (index % 5 == 0 && sequence.size() >= 12) {
            const std::size_t position = random() % (sequence.size() - 11);
            sequence.replace(position, 12, "ACGTACGTACGT");
        }
        sequences.push_back(std::move(sequence));
    }

    QualityAnalyzer reference(ReadDirection::R1, adapterConfig);
    QualityAnalyzer optimized(ReadDirection::R1, adapterConfig);
    FastqRecord record;
    for (const std::string& sequence : sequences) {
        record.sequence = sequence;
        record.quality.assign(sequence.size(), 'I');
        reference.processRecord(record);
        optimized.processRecord(record);
        reference.analyzeAdaptersReference(record);
        optimized.analyzeAdapters(record);
    }

    if (reference.adapterPosCounts != optimized.adapterPosCounts
        || !closeEnough(reference.getStats().percentWithAdapter,
                        optimized.getStats().percentWithAdapter)) {
        std::cerr << "Adapter Content optimized search differs from reference\n";
        return false;
    }
    return true;
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

    return testExactAggregationDifferential() && testAdapterContentDifferential()
        ? 0 : 1;
}
