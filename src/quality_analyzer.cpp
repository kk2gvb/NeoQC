#include "../include/quality_analyzer.h"
#include <algorithm>
#include <cmath>

QualityAnalyzer::QualityAnalyzer() {
    adapterPosCounts.resize(adapters.size());
    qualityDistribution.assign(50, 0); // Phred scores 0..49
}

void QualityAnalyzer::processRecord(const FastqRecord& record) {
    const std::string& seq = record.sequence;
    const std::string& qual = record.quality;
    size_t len = seq.length();

    if (len == 0) return;

    // Длина
    totalReads++;
    totalLength += len;
    if (len < minLength) minLength = len;
    if (len > maxLength) maxLength = len;

    // Расширение массивов качества, если нужно
    if (len > qualitySum.size()) {
        qualitySum.resize(len, 0);
        qualityCount.resize(len, 0);
    }

    // Подсчёт оснований и качества
    uint64_t readQualSum = 0;
    uint64_t validQualityBases = 0;
    for (size_t i = 0; i < len; ++i) {
        char c = seq[i];
        // Подсчёт оснований
        switch (c) {
            case 'A': case 'a': countA++; totalBases++; break;
            case 'C': case 'c': countC++; totalBases++; totalGC++; break;
            case 'G': case 'g': countG++; totalBases++; totalGC++; break;
            case 'T': case 't': countT++; totalBases++; break;
            case 'N': case 'n': countN++; totalBases++; break;
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
                readQualSum += q;
                validQualityBases++;

                if (q < static_cast<int>(qualityDistribution.size())) {
                    qualityDistribution[q]++;
                }

                if (q >= 20) q20Count++;
                if (q >= 30) q30Count++;
            }
        }
    }

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
        const std::string& adapter = adapters[aid].sequence;
        const size_t position = seq.find(adapter);

        if (position == std::string::npos) {
            continue;
        }

        const size_t adapterLength = adapter.length();

        if (adapterPosCounts[aid].size() < position + adapterLength) {
            adapterPosCounts[aid].resize(position + adapterLength, 0);
        }

        for (size_t j = 0; j < adapterLength; ++j) {
            adapterPosCounts[aid][position + j]++;
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

    stats.avgGC      = (totalBases > 0) ? static_cast<double>(totalGC) / totalBases * 100.0 : 0.0;
    stats.percentN   = (totalBases > 0) ? static_cast<double>(countN) / totalBases * 100.0 : 0.0;
    stats.percentQ20 = (totalBases > 0) ? static_cast<double>(q20Count) / totalBases * 100.0 : 0.0;
    stats.percentQ30 = (totalBases > 0) ? static_cast<double>(q30Count) / totalBases * 100.0 : 0.0;
    stats.percentWithAdapter = (totalReads > 0) ? static_cast<double>(readsWithAdapter) / totalReads * 100.0 : 0.0;

    // Качество по позициям
    stats.meanQualityPerPosition.resize(maxLength);
    for (size_t i = 0; i < maxLength; ++i) {
        if (i < qualityCount.size() && qualityCount[i] > 0) {
            stats.meanQualityPerPosition[i] = static_cast<double>(qualitySum[i]) / qualityCount[i];
        } else {
            stats.meanQualityPerPosition[i] = 0.0;
        }
    }

    // Распределение качества
    stats.qualityDistribution = qualityDistribution;
    stats.perSequenceQualityDistribution = perSequenceQualityDistribution;

    return stats;
}
