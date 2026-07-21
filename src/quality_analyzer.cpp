#include "../include/quality_analyzer.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include "../include/utils.hpp"

QualityAnalyzer::QualityAnalyzer(){

}


void QualityAnalyzer::processRecord(const FastqRecord& record) {
    size_t len = record.sequence.length();
    if (len > qualitySum.size()) {
        qualitySum.resize(len, 0);
        qualityCount.resize(len, 0);
    }

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
    
    std::cout << "\nMean quality per position:\n";
    for (size_t i = 0; i < (size_t)(stats.meanQualityPerPosition.size()); ++i) {
        std::cout << "Pos " << i+1 << ": " << std::fixed << std::setprecision(2) 
                  << stats.meanQualityPerPosition[i] << "\n";
    }
}

void QualityAnalyzer::analyzeAdapters(const FastqRecord& record) {
    const std::string& seq = record.sequence;
    if (seq.length() < illuminaUniversal.length()) return;

    // Ищем адаптер в последовательности
    for (size_t i = 0; i <= seq.length() - illuminaUniversal.length(); ++i) {
        if (seq.substr(i, illuminaUniversal.length()) == illuminaUniversal) {
            if (i >= adapterCounts.size()) {
                adapterCounts.resize(i + 100, 0);
            }
            adapterCounts[i]++;
            break;
        }
    }
}

void QualityAnalyzer::printAdapterStats(const std::string& filename, const std::string& folder) const {
    std::cout << "\n=== Adapter Content ===\n";
    bool found = false;
    std::ofstream outputFile("../results/adapter_stats_"+Utils::trim_path(filename, folder)+".txt");

    for (size_t i = 0; i < adapterCounts.size(); ++i) {
        if (adapterCounts[i] > 0) {
            double percent = (double)adapterCounts[i] / totalReads * 100.0;
            if (outputFile.is_open()) {
                outputFile << i+1 << " " << percent << "\n";
            }
            if (percent > 0.1) {  
                std::cout << "Pos " << i+1 << ": " << percent << "% adapters\n";
                found = true;
            }
        }
    }
    outputFile.close();
    
    if (!found) {
        std::cout << "No significant adapter content detected.\n";
    }
}