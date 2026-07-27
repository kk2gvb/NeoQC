#include "../include/fastq_reader.h"
#include "../include/quality_analyzer.h"
#include "../include/plot_runner.h"

#include <iostream>
#include <string>
#include <vector>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <stdexcept>

namespace fs = std::filesystem;

// ---------------------------------------------------------------------------
// Аргументы командной строки
// ---------------------------------------------------------------------------
struct Args {
    std::string r1;
    std::string r2;          // пустая строка = single-end
    std::string sampleId;
    std::string outDir;
    bool        plot = false;
};

Args parseArgs(int argc, char* argv[]) {
    Args args;
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        auto needValue = [&](const std::string& name) -> std::string {
            if (i + 1 >= argc) {
                throw std::runtime_error("Option " + name + " requires a value");
            }
            return argv[++i];
        };

        if      (arg == "--r1")        args.r1       = needValue("--r1");
        else if (arg == "--r2")        args.r2       = needValue("--r2");
        else if (arg == "--sample-id") args.sampleId = needValue("--sample-id");
        else if (arg == "--out")       args.outDir   = needValue("--out");
        else if (arg == "--plot")      args.plot     = true;
        else if (arg == "--help" || arg == "-h") {
            throw std::runtime_error("help");
        } else {
            throw std::runtime_error("Unknown argument: " + arg);
        }
    }

    if (args.r1.empty())       throw std::runtime_error("--r1 is required");
    if (args.sampleId.empty()) throw std::runtime_error("--sample-id is required");
    if (args.outDir.empty())   throw std::runtime_error("--out is required");

    return args;
}

void printUsage(const char* progName) {
    std::cerr <<
        "NeoQC — FASTQ quality analysis\n\n"
        "Usage:\n"
        "  single-end:\n"
        "    " << progName << " --r1 <file> --sample-id <id> --out <dir> [--plot]\n\n"
        "  paired-end:\n"
        "    " << progName << " --r1 <file> --r2 <file> --sample-id <id> --out <dir> [--plot]\n\n"
        "Options:\n"
        "  --r1 <file>       Path to R1 FASTQ (plain or .gz)\n"
        "  --r2 <file>       Path to R2 FASTQ (optional, for paired-end)\n"
        "  --sample-id <id>  Sample identifier (used in output filenames)\n"
        "  --out <dir>       Output directory (created if missing)\n"
        "  --plot            Build plots via plot_results.py (optional)\n";
}

// ---------------------------------------------------------------------------
// Вывод результатов (временная реализация — позже вынесется в
// ConsoleReporter и ResultWriter)
// ---------------------------------------------------------------------------
void printConsoleSummary(const QualityStats& stats, const std::string& readName) {
    std::cout << "\n=== " << readName << " Summary ===\n";
    std::cout << "Processed reads : " << stats.totalReads << "\n";
    std::cout << "Total bases     : " << stats.totalBases << "\n";
    std::cout << "Avg length      : " << stats.avgLength << "\n";
    std::cout << "GC content      : " << std::fixed << std::setprecision(2)
              << stats.avgGC << "%\n";
    std::cout << "%N              : " << stats.percentN << "%\n";
    std::cout << "%Q20            : " << stats.percentQ20 << "%\n";
    std::cout << "%Q30            : " << stats.percentQ30 << "%\n";
    std::cout << "% with adapter  : " << stats.percentWithAdapter << "%\n";
}

void writeSummaryTxt(const QualityStats& stats,
                     const std::string& outDir,
                     const std::string& filename) {
    std::string path = outDir + "/" + filename + "_summary.txt";
    std::ofstream out(path);
    if (!out) throw std::runtime_error("Cannot write to " + path);

    out << "Processed reads: " << stats.totalReads << "\n"
        << "Total bases: "     << stats.totalBases << "\n"
        << "Avg length: "      << stats.avgLength  << "\n"
        << "GC: "              << stats.avgGC      << "%\n"
        << "%N: "              << stats.percentN   << "%\n"
        << "%Q20: "            << stats.percentQ20 << "%\n"
        << "%Q30: "            << stats.percentQ30 << "%\n"
        << "% with adapter: "  << stats.percentWithAdapter << "%\n";
}

void writeAdapterTsv(const std::vector<QualityAnalyzer::Adapter>& adapters,
                     const std::vector<std::vector<uint64_t>>& adapterPosCounts,
                     size_t totalReads,
                     size_t maxLength,              // <-- новый параметр
                     const std::string& outDir,
                     const std::string& readName) {
    std::string path = outDir + "/adapter_content_" + readName + ".tsv";
    std::ofstream out(path);
    if (!out) throw std::runtime_error("Cannot write to " + path);

    // Заголовок
    out << "pos";
    for (const auto& adapter : adapters) {
        out << "\t" << adapter.name;
    }
    out << "\n";

    // Используем maxLength, а не maxPos из счётчиков
    for (size_t i = 0; i < maxLength; ++i) {
        out << (i + 1);
        for (size_t aid = 0; aid < adapterPosCounts.size(); ++aid) {
            double percent = 0.0;
            if (i < adapterPosCounts[aid].size() && totalReads > 0) {
                percent = static_cast<double>(adapterPosCounts[aid][i])
                        / static_cast<double>(totalReads) * 100.0;
            }
            out << "\t" << std::fixed << std::setprecision(4) << percent;
        }
        out << "\n";
    }
}

// ---------------------------------------------------------------------------
// Обработка одного файла (R1 или R2)
// ---------------------------------------------------------------------------
void processOneFile(const std::string& path,
                    const std::string& readName,
                    const std::string& outDir,
                    const std::string& sampleId) {
    QualityAnalyzer analyzer;

    {
        FastqReader reader(path);
        FastqRecord rec;
        size_t count = 0;

        while (reader.readNext(rec)) {
            analyzer.processRecord(rec);
            analyzer.analyzeAdapters(rec);
            ++count;
            if (count % 1'000'000 == 0) {
                std::cout << readName << ": processed " << count << " reads...\n";
            }
        }
        std::cout << readName << ": total reads = " << count << "\n";
    }

    QualityStats stats = analyzer.getStats();

    printConsoleSummary(stats, readName);
    writeSummaryTxt(stats, outDir, sampleId + "_" + readName);
    writeAdapterTsv(analyzer.adapters,
                    analyzer.adapterPosCounts,
                    stats.totalReads,
                    stats.maxLength,             // <-- передаём maxLength
                    outDir,
                    readName);
}

// ---------------------------------------------------------------------------
// Обработка парных файлов (R1 и R2)
// ---------------------------------------------------------------------------

void processPairedFiles(const std::string& r1Path,
                        const std::string& r2Path,
                        const std::string& outDir,
                        const std::string& sampleId)
{
    FastqReader readerR1(r1Path);
    FastqReader readerR2(r2Path);

    QualityAnalyzer analyzerR1;
    QualityAnalyzer analyzerR2;

    FastqRecord rec1;
    FastqRecord rec2;

    size_t count = 0;

    while (true)
    {
        bool ok1 = readerR1.readNext(rec1);
        bool ok2 = readerR2.readNext(rec2);

        // Оба файла закончились
        if (!ok1 && !ok2)
            break;

        // R1 закончился раньше
        if (!ok1)
        {
            throw std::runtime_error(
                "FASTQ validation error:\n"
                "reason: R1 contains fewer reads than R2");
        }

        // R2 закончился раньше
        if (!ok2)
        {
            throw std::runtime_error(
                "FASTQ validation error:\n"
                "reason: R2 contains fewer reads than R1");
        }



        analyzerR1.processRecord(rec1);
        analyzerR1.analyzeAdapters(rec1);

        analyzerR2.processRecord(rec2);
        analyzerR2.analyzeAdapters(rec2);

        ++count;

        if (count % 1000000 == 0)
        {
            std::cout << "Processed "
                      << count
                      << " paired reads...\n";
        }
    }

    std::cout << "R1: total reads = "
              << analyzerR1.getStats().totalReads
              << "\n";

    std::cout << "R2: total reads = "
              << analyzerR2.getStats().totalReads
              << "\n";

    QualityStats statsR1 = analyzerR1.getStats();
    QualityStats statsR2 = analyzerR2.getStats();

    printConsoleSummary(statsR1, "R1");
    printConsoleSummary(statsR2, "R2");

    writeSummaryTxt(statsR1, outDir, sampleId + "_R1");
    writeSummaryTxt(statsR2, outDir, sampleId + "_R2");

    writeAdapterTsv(
        analyzerR1.adapters,
        analyzerR1.adapterPosCounts,
        statsR1.totalReads,
        statsR1.maxLength,
        outDir,
        "R1");

    writeAdapterTsv(
        analyzerR2.adapters,
        analyzerR2.adapterPosCounts,
        statsR2.totalReads,
        statsR2.maxLength,
        outDir,
        "R2");
}

// ---------------------------------------------------------------------------
// main
// ---------------------------------------------------------------------------
int main(int argc, char* argv[]) {
    if (argc < 2) {
        printUsage(argv[0]);
        return 1;
    }

    Args args;
    try {
        args = parseArgs(argc, argv);
    } catch (const std::exception& e) {
        std::string msg = e.what();
        if (msg == "help") {
            printUsage(argv[0]);
            return 0;
        }
        std::cerr << "Error: " << msg << "\n\n";
        printUsage(argv[0]);
        return 1;
    }

    // Проверка входных файлов
    if (!fs::exists(args.r1)) {
        std::cerr << "Error: R1 file not found: " << args.r1 << "\n";
        return 1;
    }
    if (!args.r2.empty() && !fs::exists(args.r2)) {
        std::cerr << "Error: R2 file not found: " << args.r2 << "\n";
        return 1;
    }

    // Создание выходного каталога
    {
        std::error_code ec;
        fs::create_directories(args.outDir, ec);
        if (ec) {
            std::cerr << "Error: cannot create output directory '"
                      << args.outDir << "': " << ec.message() << "\n";
            return 1;
        }
    }

    const bool isPaired = !args.r2.empty();
    std::cout << "Sample ID : " << args.sampleId << "\n"
              << "Mode      : " << (isPaired ? "paired-end" : "single-end") << "\n"
              << "R1        : " << args.r1 << "\n";
    if (isPaired) std::cout << "R2        : " << args.r2 << "\n";
    std::cout << "Output    : " << args.outDir << "\n";

    try
    {
        if (isPaired)
        {
            processPairedFiles(
                args.r1,
                args.r2,
                args.outDir,
                args.sampleId);
        }
        else
        {
            processOneFile(
                args.r1,
                "R1",
                args.outDir,
                args.sampleId);
        }
    }
    catch (const std::exception& e)
    {
        std::cerr << e.what() << '\n';
        return 1;
    }

    // Построение графиков (опционально, через PlotRunner)
    if (args.plot) {
        std::string plotDir = args.outDir + "/plots";
        PlotRunner::run(args.outDir + "/adapter_content_R1.tsv", plotDir);
        if (isPaired) {
            PlotRunner::run(args.outDir + "/adapter_content_R2.tsv", plotDir);
        }
    }

    std::cout << "\nDone.\n";
    return 0;
}