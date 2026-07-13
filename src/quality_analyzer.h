#pragma once
#include "fastq_reader.h"
#include <vector>
#include <string>

struct QualityStats {
    std::vector<double> meanQualityPerPosition;
    double avgGC = 0.0;
    size_t totalReads = 0;
    size_t totalBases = 0;
};

class QualityAnalyzer {
public:
    QualityAnalyzer(size_t maxReadLength = 150);
    void processRecord(const FastqRecord& record);
    QualityStats getStats() const;
    void printSummary() const;
    void analyzeAdapters(const FastqRecord& record);
    void printAdapterStats(const std::string& filename) const;

private:
    std::vector<long long> qualitySum;   // sum of quality scores per position
    std::vector<long long> qualityCount;
    long long totalGC = 0;
    long long totalBases = 0;
    size_t totalReads = 0;
    size_t maxLength = 0;
    std::vector<long long> adapterCounts;
    const std::string illuminaUniversal = "AGATCGGAAGAG";
};