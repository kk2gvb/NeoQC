#include "quality_analyzer.h"
#include <iostream>
#include <iomanip>
#include <algorithm>

QualityAnalyzer::QualityAnalyzer(size_t maxReadLength) {
    qualitySum.resize(maxReadLength, 0);
    qualityCount.resize(maxReadLength, 0);
}

void QualityAnalyzer::processRecord(const FastqRecord& record) {
    size_t len = record.sequence.length();
    if (len > maxLength) maxLength = len;

    // Quality
    for (size_t i = 0; i < len && i < qualitySum.size(); ++i) {
        int q = record.quality[i] - 33;  // Phred+33
        if (q >= 0) {
            qualitySum[i] += q;
            qualityCount[i]++;
        }
    }

    // GC content
    for (char c : record.sequence) {
        if (c == 'G' || c == 'C') totalGC++;
        totalBases++;
    }

    totalReads++;
}

QualityStats QualityAnalyzer::getStats() const {
    QualityStats stats;
    stats.totalReads = totalReads;
    stats.totalBases = totalBases;
    stats.avgGC = totalBases > 0 ? (double)totalGC / totalBases * 100.0 : 0.0;

    stats.meanQualityPerPosition.resize(maxLength);
    for (size_t i = 0; i < maxLength; ++i) {
        stats.meanQualityPerPosition[i] = qualityCount[i] > 0 ? 
            (double)qualitySum[i] / qualityCount[i] : 0.0;
    }
    return stats;
}

void QualityAnalyzer::printSummary() const {
    auto stats = getStats();
    std::cout << "\n=== Quality Analysis Summary ===\n";
    std::cout << "Processed reads: " << stats.totalReads << "\n";
    std::cout << "Total bases: " << stats.totalBases << "\n";
    std::cout << "Average GC: " << std::fixed << std::setprecision(2) 
              << stats.avgGC << "%\n";
    
    std::cout << "\nMean quality per position (first 50):\n";
    for (size_t i = 0; i < std::min<size_t>(50, stats.meanQualityPerPosition.size()); ++i) {
        std::cout << "Pos " << i+1 << ": " << std::fixed << std::setprecision(2) 
                  << stats.meanQualityPerPosition[i] << "\n";
    }
}