#include "../include/quality_analyzer.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <utility>

namespace {
constexpr size_t kAdapterDetectionKmerLength = 12;
constexpr std::array<const char*, 16> kDuplicationLabels = {
    "1", "2", "3", "4", "5", "6", "7", "8", "9",
    ">10", ">50", ">100", ">500", ">1k", ">5k", ">10k+"
};

std::size_t duplicationSlot(uint64_t duplicationLevel) {
    if (duplicationLevel > 10000) return 15;
    if (duplicationLevel > 5000) return 14;
    if (duplicationLevel > 1000) return 13;
    if (duplicationLevel > 500) return 12;
    if (duplicationLevel > 100) return 11;
    if (duplicationLevel > 50) return 10;
    if (duplicationLevel > 10) return 9;
    return static_cast<std::size_t>(duplicationLevel - 1);
}

uint64_t duplicationBaseCode(char base) {
    switch (base) {
        case 'A': case 'a': return 1;
        case 'C': case 'c': return 2;
        case 'G': case 'g': return 3;
        case 'T': case 't': return 4;
        case 'N': case 'n': return 5;
        default: return 0;
    }
}

char duplicationBase(uint64_t code) {
    constexpr std::array<char, 6> bases = {'\0', 'A', 'C', 'G', 'T', 'N'};
    return code < bases.size() ? bases[code] : '\0';
}

DuplicationKey encodeDuplicationKey(const std::string& sequence) {
    DuplicationKey key;
    const std::size_t length = std::min(DUPLICATION_PREFIX_LENGTH, sequence.size());
    for (std::size_t position = 0; position < length; ++position) {
        const uint64_t code = duplicationBaseCode(sequence[position]);
        const std::size_t bit = position * 3;
        const std::size_t word = bit / 64;
        const std::size_t offset = bit % 64;
        key.words[word] |= code << offset;
        if (offset > 61) {
            key.words[word + 1] |= code >> (64 - offset);
        }
    }
    return key;
}

std::string decodeDuplicationKey(const DuplicationKey& key) {
    std::string sequence;
    sequence.reserve(DUPLICATION_PREFIX_LENGTH);
    for (std::size_t position = 0; position < DUPLICATION_PREFIX_LENGTH; ++position) {
        const std::size_t bit = position * 3;
        const std::size_t word = bit / 64;
        const std::size_t offset = bit % 64;
        uint64_t code = (key.words[word] >> offset) & 0x7;
        if (offset > 61) {
            code |= (key.words[word + 1] << (64 - offset)) & 0x7;
        }
        const char base = duplicationBase(code);
        if (base == '\0') break;
        sequence.push_back(base);
    }
    return sequence;
}

double qualityQuantile(const std::array<uint64_t, 94>& histogram,
                       uint64_t count,
                       double probability) {
    if (count == 0) return 0.0;
    const uint64_t target = std::max<uint64_t>(
        1, static_cast<uint64_t>(std::ceil(probability * count)));
    uint64_t cumulative = 0;
    for (size_t quality = 0; quality < histogram.size(); ++quality) {
        cumulative += histogram[quality];
        if (cumulative >= target) return static_cast<double>(quality);
    }
    return static_cast<double>(histogram.size() - 1);
}
}

std::size_t DuplicationKeyHash::operator()(const DuplicationKey& key) const noexcept {
    uint64_t hash = 0x9e3779b97f4a7c15ULL;
    for (uint64_t word : key.words) {
        word ^= word >> 30;
        word *= 0xbf58476d1ce4e5b9ULL;
        word ^= word >> 27;
        word *= 0x94d049bb133111ebULL;
        word ^= word >> 31;
        hash ^= word + 0x9e3779b97f4a7c15ULL + (hash << 6) + (hash >> 2);
    }
    return static_cast<std::size_t>(hash);
}

QualityAnalyzer::QualityAnalyzer(ReadDirection direction) {
    if (direction == ReadDirection::R1) {
        adapters = {
            {"TruSeq_R1", "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA", ""},
            {"SmallRNA3'", "TGGAATTCTCGGGTGCCAAGG", ""},
            {"SmallRNA5'", "GTTCAGAGTTCTACAGTCCGACGATC", ""},
            {"Nextera", "CTGTCTCTTATACACATCT", ""}
        };
    } else {
        adapters = {
            {"TruSeq_R2", "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT", ""},
            {"SmallRNA3'", "TGGAATTCTCGGGTGCCAAGG", ""},
            {"SmallRNA5'", "GTTCAGAGTTCTACAGTCCGACGATC", ""},
            {"Nextera", "CTGTCTCTTATACACATCT", ""}
        };
    }

    for (auto& adapter : adapters) {
        adapter.detectionSequence = adapter.sequence.substr(
            0, std::min(kAdapterDetectionKmerLength, adapter.sequence.size()));
    }

    adapterPosCounts.resize(adapters.size());
}

void QualityAnalyzer::processRecord(const FastqRecord& record) {
    const std::string& seq = record.sequence;
    const std::string& qual = record.quality;
    size_t len = seq.length();

    if (len == 0) return;

    totalReads++;
    totalLength += len;

    sequenceCounts[encodeDuplicationKey(seq)]++;

    if (len >= lengthDistribution.size())
    {
        lengthDistribution.resize(len + 1, 0);
    }

    lengthDistribution[len]++;

    if (len < minLength) minLength = len;
    if (len > maxLength) maxLength = len;

    // Расширение массивов качества, если нужно
    if (len > qualitySum.size()) {
        qualitySum.resize(len, 0);
        qualityCount.resize(len, 0);
        qualityHistogram.resize(len);

        baseCountA.resize(len, 0);
        baseCountC.resize(len, 0);
        baseCountG.resize(len, 0);
        baseCountT.resize(len, 0);
        baseCountN.resize(len, 0);

        readsPerPosition.resize(len, 0);
    }

    // Подсчёт оснований и качества
    uint64_t readQualSum = 0;
    uint64_t validQualityBases = 0;
    size_t gc = 0;
    for (size_t i = 0; i < len; ++i) {
        readsPerPosition[i]++;
        char c = seq[i];
        // Подсчёт оснований
        switch (c) {
            case 'A': case 'a': countA++; totalBases++; baseCountA[i]++; break;
            case 'C': case 'c': countC++; totalBases++; totalGC++; baseCountC[i]++; gc++; break;
            case 'G': case 'g': countG++; totalBases++; totalGC++; baseCountG[i]++; gc++; break;
            case 'T': case 't': countT++; totalBases++; baseCountT[i]++; break;
            case 'N': case 'n': countN++; totalBases++; baseCountN[i]++; break;
            default:
                // Неизвестный символ — считаем как N
                countN++; totalBases++;
                break;
        }

        // Качество (только если строка качества достаточно длинная)
        if (i < qual.length()) {
            int q = static_cast<int>(qual[i]) - 33; // Phred+33
            if (q >= 0) {
                qualitySum[i] += q;
                qualityCount[i]++;
                if (static_cast<size_t>(q) < qualityHistogram[i].size()) {
                    qualityHistogram[i][static_cast<size_t>(q)]++;
                }
                readQualSum += q;
                validQualityBases++;

                if (q >= 20) q20Count++;
                if (q >= 30) q30Count++;
            }
        }
    }

    int gcPercent = static_cast<int>(
        std::lround(static_cast<double>(gc) * 100.0 / len));

    gcDistribution[gcPercent]++;

    if (validQualityBases > 0) {
        const auto meanQuality = static_cast<size_t>(std::lround(
            static_cast<double>(readQualSum) / validQualityBases));
        if (meanQuality >= perSequenceQualityDistribution.size()) {
            perSequenceQualityDistribution.resize(meanQuality + 1, 0);
        }
        perSequenceQualityDistribution[meanQuality]++;
    }

}

void QualityAnalyzer::analyzeAdapters(const FastqRecord& record) {
    const std::string& seq = record.sequence;
    bool foundInRead = false;

    for (size_t aid = 0; aid < adapters.size(); ++aid) {
        const std::string& detectionSequence = adapters[aid].detectionSequence;
        const size_t position = seq.find(detectionSequence);

        if (position == std::string::npos) {
            continue;
        }

        if (adapterPosCounts[aid].size() < seq.size()) {
            adapterPosCounts[aid].resize(seq.size(), 0);
        }

        // FastQC reports a cumulative trace: once an adapter k-mer is found,
        // the read is counted from that position through its end.
        for (size_t j = position; j < seq.size(); ++j) {
            adapterPosCounts[aid][j]++;
        }

        foundInRead = true;
    }

    if (foundInRead) {
        readsWithAdapter++;
    }
}

QualityStats QualityAnalyzer::getStats() const {
    QualityStats stats;

    stats.totalReads = totalReads;
    stats.totalBases = totalBases;
    stats.minLength  = (totalReads > 0) ? minLength : 0;
    stats.maxLength  = maxLength;
    stats.avgLength  = (totalReads > 0) ? static_cast<double>(totalLength) / totalReads : 0.0;

    stats.countA = countA;
    stats.countC = countC;
    stats.countG = countG;
    stats.countT = countT;
    stats.countN = countN;

    stats.baseCountA = baseCountA;
    stats.baseCountC = baseCountC;
    stats.baseCountG = baseCountG;
    stats.baseCountT = baseCountT;
    stats.baseCountN = baseCountN;

    stats.readsPerPosition = readsPerPosition;

    stats.gcDistribution = gcDistribution;

    stats.lengthDistribution = lengthDistribution;

    stats.avgGC      = (totalBases > 0) ? static_cast<double>(totalGC) / totalBases * 100.0 : 0.0;
    stats.percentN   = (totalBases > 0) ? static_cast<double>(countN) / totalBases * 100.0 : 0.0;
    stats.percentQ20 = (totalBases > 0) ? static_cast<double>(q20Count) / totalBases * 100.0 : 0.0;
    stats.percentQ30 = (totalBases > 0) ? static_cast<double>(q30Count) / totalBases * 100.0 : 0.0;
    stats.percentWithAdapter = (totalReads > 0) ? static_cast<double>(readsWithAdapter) / totalReads * 100.0 : 0.0;

    // Качество по позициям
    stats.meanQualityPerPosition.resize(maxLength);
    stats.lowerQuartileQualityPerPosition.resize(maxLength);
    stats.medianQualityPerPosition.resize(maxLength);
    for (size_t i = 0; i < maxLength; ++i) {
        if (i < qualityCount.size() && qualityCount[i] > 0) {
            stats.meanQualityPerPosition[i] = static_cast<double>(qualitySum[i]) / qualityCount[i];
            stats.lowerQuartileQualityPerPosition[i] = qualityQuantile(
                qualityHistogram[i], qualityCount[i], 0.25);
            stats.medianQualityPerPosition[i] = qualityQuantile(
                qualityHistogram[i], qualityCount[i], 0.50);
        } else {
            stats.meanQualityPerPosition[i] = 0.0;
            stats.lowerQuartileQualityPerPosition[i] = 0.0;
            stats.medianQualityPerPosition[i] = 0.0;
        }
    }

    stats.perSequenceQualityDistribution = perSequenceQualityDistribution;

    return stats;
}

DuplicationStats QualityAnalyzer::getDuplicationStats() const {
    DuplicationStats stats;
    stats.totalReads = totalReads;
    stats.uniqueSequences = sequenceCounts.size();
    stats.levels.reserve(kDuplicationLabels.size());

    std::unordered_map<uint64_t, uint64_t> collatedCounts;
    collatedCounts.reserve(std::min<std::size_t>(sequenceCounts.size(), 16384));
    for (const auto& [sequence, count] : sequenceCounts) {
        (void)sequence;
        collatedCounts[count]++;
    }

    std::array<double, 16> rawByLevel{};
    std::array<double, 16> deduplicatedByLevel{};
    double rawTotal = 0.0;
    double deduplicatedTotal = 0.0;
    for (const auto& [level, observations] : collatedCounts) {
        const double exactCount = static_cast<double>(observations);
        const std::size_t slot = duplicationSlot(level);
        rawByLevel[slot] += exactCount * static_cast<double>(level);
        deduplicatedByLevel[slot] += exactCount;
        rawTotal += exactCount * static_cast<double>(level);
        deduplicatedTotal += exactCount;
    }

    stats.deduplicatedRemainingPercent = rawTotal > 0.0
        ? 100.0 * deduplicatedTotal / rawTotal
        : 100.0;
    for (std::size_t i = 0; i < kDuplicationLabels.size(); ++i) {
        stats.levels.push_back({
            kDuplicationLabels[i],
            rawTotal > 0.0 ? 100.0 * rawByLevel[i] / rawTotal : 0.0,
            deduplicatedTotal > 0.0
                ? 100.0 * deduplicatedByLevel[i] / deduplicatedTotal
                : 0.0,
        });
    }

    for (const auto& [sequence, count] : sequenceCounts) {
        const double percent = totalReads > 0
            ? 100.0 * static_cast<double>(count) / static_cast<double>(totalReads)
            : 0.0;
        if (percent > OVERREPRESENTED_SEQUENCE_THRESHOLD) {
            stats.overrepresentedSequences.push_back({
                decodeDuplicationKey(sequence), count, percent
            });
        }
    }
    std::sort(
        stats.overrepresentedSequences.begin(),
        stats.overrepresentedSequences.end(),
        [](const auto& left, const auto& right) {
            return left.count != right.count
                ? left.count > right.count
                : left.sequence < right.sequence;
        });

    return stats;
}

DuplicationStats QualityAnalyzer::getDuplicationStats(
    const std::vector<DuplicationEntry>& entries) const
{
    DuplicationStats stats;

    stats.totalReads = totalReads;
    stats.uniqueSequences = entries.size();

    std::unordered_map<uint64_t, uint64_t> collatedCounts;

    for (const auto& entry : entries)
    {
        collatedCounts[entry.count]++;
    }

    std::array<double, 16> rawByLevel{};
    std::array<double, 16> deduplicatedByLevel{};

    double rawTotal = 0.0;
    double deduplicatedTotal = 0.0;

    for (const auto& [level, observations] : collatedCounts)
    {
        const double exactCount =
            static_cast<double>(observations);

        const std::size_t slot = duplicationSlot(level);

        rawByLevel[slot] +=
            exactCount * static_cast<double>(level);

        deduplicatedByLevel[slot] += exactCount;

        rawTotal +=
            exactCount * static_cast<double>(level);

        deduplicatedTotal += exactCount;
    }

    stats.deduplicatedRemainingPercent =
        rawTotal > 0.0
            ? 100.0 * deduplicatedTotal / rawTotal
            : 100.0;

    stats.levels.reserve(kDuplicationLabels.size());

    for (std::size_t i = 0;
         i < kDuplicationLabels.size();
         ++i)
    {
        stats.levels.push_back({
            kDuplicationLabels[i],
            rawTotal > 0.0
                ? 100.0 * rawByLevel[i] / rawTotal
                : 0.0,
            deduplicatedTotal > 0.0
                ? 100.0 * deduplicatedByLevel[i]
                    / deduplicatedTotal
                : 0.0
        });
    }

    for (const auto& entry : entries)
    {
        const double percent =
            totalReads > 0
                ? 100.0 *
                    static_cast<double>(entry.count) /
                    static_cast<double>(totalReads)
                : 0.0;

        if (percent > OVERREPRESENTED_SEQUENCE_THRESHOLD)
        {
            stats.overrepresentedSequences.push_back({
                decodeDuplicationKey(entry.key),
                entry.count,
                percent
            });
        }
    }

    std::sort(
        stats.overrepresentedSequences.begin(),
        stats.overrepresentedSequences.end(),
        [](const auto& left, const auto& right)
        {
            return left.count != right.count
                ? left.count > right.count
                : left.sequence < right.sequence;
        });

    return stats;
}

std::vector<DuplicationEntry> QualityAnalyzer::getDuplicationEntries() const
{
    std::vector<DuplicationEntry> entries;
    entries.reserve(sequenceCounts.size());

    for (const auto& [key, count] : sequenceCounts)
    {
        entries.push_back({key, count});
    }

    return entries;
}


void QualityAnalyzer::merge(const QualityAnalyzer& other)
{
    // ---------------------------------------------------------------------
    // Простые счётчики
    // ---------------------------------------------------------------------

    totalGC += other.totalGC;
    totalBases += other.totalBases;
    totalReads += other.totalReads;

    countA += other.countA;
    countC += other.countC;
    countG += other.countG;
    countT += other.countT;
    countN += other.countN;

    totalLength += other.totalLength;

    q20Count += other.q20Count;
    q30Count += other.q30Count;

    readsWithAdapter += other.readsWithAdapter;

    minLength = std::min(minLength, other.minLength);
    maxLength = std::max(maxLength, other.maxLength);

    // ---------------------------------------------------------------------
    // Вспомогательные функции
    // ---------------------------------------------------------------------

    auto mergeVector = [](auto& lhs, const auto& rhs)
    {
        if (lhs.size() < rhs.size())
            lhs.resize(rhs.size(), 0);

        for (size_t i = 0; i < rhs.size(); ++i)
            lhs[i] += rhs[i];
    };

    // ---------------------------------------------------------------------
    // Векторы
    // ---------------------------------------------------------------------

    mergeVector(baseCountA, other.baseCountA);
    mergeVector(baseCountC, other.baseCountC);
    mergeVector(baseCountG, other.baseCountG);
    mergeVector(baseCountT, other.baseCountT);
    mergeVector(baseCountN, other.baseCountN);

    mergeVector(readsPerPosition, other.readsPerPosition);

    mergeVector(gcDistribution, other.gcDistribution);
    mergeVector(lengthDistribution, other.lengthDistribution);

    mergeVector(qualitySum, other.qualitySum);
    mergeVector(qualityCount, other.qualityCount);

    mergeVector(perSequenceQualityDistribution,
                other.perSequenceQualityDistribution);

    // ---------------------------------------------------------------------
    // Quality histogram
    // ---------------------------------------------------------------------

    if (qualityHistogram.size() < other.qualityHistogram.size())
    {
        qualityHistogram.resize(other.qualityHistogram.size());
    }

    for (size_t pos = 0; pos < other.qualityHistogram.size(); ++pos)
    {
        for (size_t q = 0; q < 94; ++q)
        {
            qualityHistogram[pos][q] += other.qualityHistogram[pos][q];
        }
    }

    // ---------------------------------------------------------------------
    // Adapter positions
    // ---------------------------------------------------------------------

    if (adapterPosCounts.size() < other.adapterPosCounts.size())
    {
        adapterPosCounts.resize(other.adapterPosCounts.size());
    }

    for (size_t adapter = 0; adapter < other.adapterPosCounts.size(); ++adapter)
    {
        mergeVector(adapterPosCounts[adapter],
                    other.adapterPosCounts[adapter]);
    }
}

    
