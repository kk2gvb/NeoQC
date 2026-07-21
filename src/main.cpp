#include "../include/fastq_reader.h"
#include "../include/utils.hpp"
#include "../include/quality_analyzer.h"
#include <iostream>

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <folder>" << std::endl;
        return 1;
    }

    std::string folder = argv[1];

    std::string r1_path = "../data/" + folder + "/" + folder + "R1_1.fq.gz";
    std::string r2_path = "../data/" + folder + "/" + folder + "R2_2.fq.gz";

    std::cout << "Processing sample: " << folder << std::endl;
    std::cout << "R1: " << r1_path << std::endl;
    std::cout << "R2: " << r2_path << std::endl;

    QualityAnalyzer analyzer1;
    QualityAnalyzer analyzer2;

    try {
        FastqReader reader1(r1_path);
        QualityAnalyzer analyzer1;

        FastqRecord rec;
        int count = 0;

        while (reader1.readNext(rec)) {
            analyzer1.processRecord(rec);
            analyzer1.analyzeAdapters(rec);
            count++;
            if (count % 1000000 == 0) {
                std::cout << "Processed " << count << " reads...\n";
            }
        }

        analyzer1.printSummary();
        analyzer1.printAdapterStats(r1_path, folder);
        {
            std::string cmd = std::string("python3 ../scripts/plot_results.py ../results/adapter_stats_")
                + Utils::trim_path(r1_path, folder) + ".txt";
            system(cmd.c_str());
        }
        std::cout << "\nDone. Total reads: " << count << std::endl;

    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }

    try {
        FastqReader reader2(r2_path);
        QualityAnalyzer analyzer2;

        FastqRecord rec;
        int count = 0;

        while (reader2.readNext(rec)) {
            analyzer2.processRecord(rec);
            analyzer2.analyzeAdapters(rec);
            count++;
            if (count % 1000000 == 0) {
                std::cout << "Processed " << count << " reads...\n";
            }
        }

        analyzer2.printSummary();
        analyzer2.printAdapterStats(r2_path, folder);
        {
            std::string cmd = std::string("python3 ../scripts/plot_results.py ../results/adapter_stats_")
                + Utils::trim_path(r2_path, folder) + ".txt";
            system(cmd.c_str());
        }
        std::cout << "\nDone. Total reads: " << count << std::endl;

    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}