#include "fastq_reader.h"
#include "quality_analyzer.h"
#include <iostream>

int main(int argc, char* argv[]) {
    std::string filename = argc > 1 ? argv[1] : "../data/test/test10k.fastq.gz";

    try {
        FastqReader reader(filename);
        QualityAnalyzer analyzer(300);  // max expected read length

        FastqRecord rec;
        int count = 0;
        const int maxReads = 100000;  // лимит для теста

        while (reader.readNext(rec) && count < maxReads) {
            analyzer.processRecord(rec);
            count++;
            if (count % 10000 == 0) {
                std::cout << "Processed " << count << " reads...\n";
            }
        }

        analyzer.printSummary();
        std::cout << "\nDone. Total reads: " << count << std::endl;

    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    return 0;
}