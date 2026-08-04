#include "../include/fastq_reader.h"
#include "../include/quality_analyzer.h"
#include "../include/plot_runner.h"
#include "../include/sample_sheet.h"

#include <iostream>
#include <string>
#include <vector>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <stdexcept>
#include <chrono>
#include <optional>
#include <ctime>

namespace fs = std::filesystem;
using Duration = std::chrono::nanoseconds;
using Clock = std::chrono::steady_clock;

struct PerformanceTimers {
    Duration fileOpen{};
    Duration readAndDecompress{};
    Duration validation{};
    Duration pairValidation{};
    Duration metrics{};
    Duration adapterSearch{};
    Duration reportWriting{};
    Duration plotting{};
};

struct AnalysisResult {
    PerformanceTimers timers;
    QualityStats r1Stats;
    std::optional<QualityStats> r2Stats;
};

struct BatchSampleResult {
    SampleSheetEntry entry;
    bool passed = false;
    std::optional<AnalysisResult> analysis;
    std::string error;
};

class ScopedTimer {
public:
    ScopedTimer(Duration& total, bool enabled)
        : total(total), enabled(enabled), start(enabled ? Clock::now() : Clock::time_point{}) {}
    ~ScopedTimer() { if (enabled) total += Clock::now() - start; }

private:
    Duration& total;
    bool enabled;
    Clock::time_point start;
};

void printPerformanceTimers(const PerformanceTimers& timers) {
    const auto milliseconds = [](Duration duration) {
        return std::chrono::duration<double, std::milli>(duration).count();
    };

    std::cout << "\n=== Performance timings ===\n" << std::fixed
              << std::setprecision(3)
              << "File opening                 : " << milliseconds(timers.fileOpen) << " ms\n"
              << "Reading and decompression    : " << milliseconds(timers.readAndDecompress) << " ms\n"
              << "FASTQ structure validation   : " << milliseconds(timers.validation) << " ms\n"
              << "Paired-read validation       : " << milliseconds(timers.pairValidation) << " ms\n"
              << "Metrics calculation          : " << milliseconds(timers.metrics) << " ms\n"
              << "Adapter search               : " << milliseconds(timers.adapterSearch) << " ms\n"
              << "Report writing               : " << milliseconds(timers.reportWriting) << " ms\n"
              << "Plot generation              : " << milliseconds(timers.plotting) << " ms\n"
              << "Total                        : " << milliseconds(timers.fileOpen + timers.readAndDecompress + timers.validation + timers.pairValidation + timers.metrics + timers.adapterSearch + timers.reportWriting + timers.plotting) << " ms\n";
}

// ---------------------------------------------------------------------------
// Аргументы командной строки
// ---------------------------------------------------------------------------
struct Args {
    std::string r1;
    std::string r2;          // пустая строка = single-end
    std::string sampleId;
    std::string outDir;
    std::string samples;
    bool        plot = false;
    bool        skipAdapters = false;
    bool        timings = false;
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
        else if (arg == "--samples")   args.samples  = needValue("--samples");
        else if (arg == "--plot")      args.plot     = true;
        else if (arg == "--skip-adapters") args.skipAdapters = true;
        else if (arg == "--timings")   args.timings  = true;
        else if (arg == "--help" || arg == "-h") {
            throw std::runtime_error("help");
        } else {
            throw std::runtime_error("Unknown argument: " + arg);
        }
    }

    if (!args.samples.empty()) {
        if (!args.r1.empty() || !args.r2.empty() || !args.sampleId.empty()) {
            throw std::runtime_error("--samples cannot be combined with --r1, --r2, or --sample-id");
        }
        return args;
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
        "    " << progName << " --r1 <file> --sample-id <id> --out <dir> [--plot] [--skip-adapters] [--timings]\n\n"
        "  paired-end:\n"
        "    " << progName << " --r1 <file> --r2 <file> --sample-id <id> --out <dir> [--plot] [--skip-adapters] [--timings]\n\n"
        "  validate sample sheet:\n"
        "    " << progName << " --samples <samples.csv> [--out <dir>] [--plot] [--skip-adapters] [--timings]\n\n"
        "Options:\n"
        "  --r1 <file>       Path to R1 FASTQ (plain or .gz)\n"
        "  --r2 <file>       Path to R2 FASTQ (optional, for paired-end)\n"
        "  --sample-id <id>  Sample identifier (used in output filenames)\n"
        "  --out <dir>       Output directory (created if missing); enables batch QC with --samples\n"
        "  --samples <file>  Validate a CSV table; combine with --out to run batch QC\n"
        "  --plot            Build plots via plot_results.py (optional)\n"
        "  --skip-adapters   Disable adapter search (for performance measurements)\n"
        "  --timings         Measure and print per-stage execution times\n";
}

void writeSummary(std::ostream& out,
                  const QualityStats& stats,
                  const std::string& readName,
                  bool adaptersSkipped)
{
    out << "\n=== " << readName << " Summary ===\n";

    out << "Processed reads : " << stats.totalReads << "\n";
    out << "Total bases     : " << stats.totalBases << "\n";

    out << "Min length      : " << stats.minLength << "\n";
    out << "Max length      : " << stats.maxLength << "\n";
    out << "Avg length      : "
        << std::fixed << std::setprecision(2)
        << stats.avgLength << "\n";

    out << "\nBase composition\n";

    out << "A               : " << stats.countA << "\n";
    out << "C               : " << stats.countC << "\n";
    out << "G               : " << stats.countG << "\n";
    out << "T               : " << stats.countT << "\n";
    out << "N               : " << stats.countN << "\n";

    out << "\n";

    out << "GC content      : " << stats.avgGC << "%\n";
    out << "%N              : " << stats.percentN << "%\n";
    out << "%Q20            : " << stats.percentQ20 << "%\n";
    out << "%Q30            : " << stats.percentQ30 << "%\n";
    if (adaptersSkipped) {
        out << "% with adapter  : not calculated (--skip-adapters)\n";
    } else {
        out << "% with adapter  : "
            << stats.percentWithAdapter << "%\n";
    }
}

// ---------------------------------------------------------------------------
// Вывод результатов (временная реализация — позже вынесется в
// ConsoleReporter и ResultWriter)
// ---------------------------------------------------------------------------
void printConsoleSummary(const QualityStats& stats,
                         const std::string& readName,
                         bool adaptersSkipped)
{
    writeSummary(std::cout, stats, readName, adaptersSkipped);
}


// ---------------------------------------------------------------------------
// Нормализация идентификатора считывания (удаление '@', '/1' и '/2', обрезка по пробелу)
// ---------------------------------------------------------------------------
std::string normalizeReadId(const std::string& header)
{
    std::string id = header;

    if (!id.empty() && id.front() == '@')
        id.erase(0, 1);

    auto space = id.find(' ');
    if (space != std::string::npos)
        id.erase(space);

    if (id.size() >= 2)
    {
        auto tail = id.substr(id.size() - 2);

        if (tail == "/1" || tail == "/2")
            id.erase(id.size() - 2);
    }

    return id;
}

void writeSummaryTxt(const QualityStats& stats,
                     const std::string& outDir,
                     const std::string& filename,
                     bool adaptersSkipped) {
    std::string path = outDir + "/" + filename + "_summary.txt";
    std::ofstream out(path);

    if (!out) throw std::runtime_error("Cannot write to " + path);
    writeSummary(out, stats, filename, adaptersSkipped);

}

void writePerCycleQualityTsv(const std::vector<double>& meanQuality,
                             const std::string& outDir,
                             const std::string& readName)
{
    std::string path = outDir + "/per_cycle_" + readName + ".tsv";

    std::ofstream out(path);

    if (!out)
        throw std::runtime_error("Cannot write to " + path);

    out << "cycle\tmean_quality\n";

    for (size_t i = 0; i < meanQuality.size(); ++i)
    {
        out << (i + 1)
            << "\t"
            << std::fixed
            << std::setprecision(4)
            << meanQuality[i]
            << "\n";
    }
}

void writePerBaseSequenceContentTsv(const std::vector<uint64_t>& baseCountA,
                                    const std::vector<uint64_t>& baseCountC,
                                    const std::vector<uint64_t>& baseCountG,
                                    const std::vector<uint64_t>& baseCountT,
                                    const std::vector<uint64_t>& baseCountN,
                                    const std::string& outDir,
                                    const std::string& readName)
{
    std::string path = outDir + "/per_base_sequence_content_" + readName + ".tsv";

    std::ofstream out(path);

    if (!out)
        throw std::runtime_error("Cannot write to " + path);

    out << "position\tA\tC\tG\tT\tN\n";

    for (size_t i = 0; i < baseCountA.size(); ++i)
    {
        const double total =
            baseCountA[i] +
            baseCountC[i] +
            baseCountG[i] +
            baseCountT[i] +
            baseCountN[i];

        double a = 0;
        double c = 0;
        double g = 0;
        double t = 0;
        double n = 0;

        if (total > 0)
        {
            a = baseCountA[i] * 100.0 / total;
            c = baseCountC[i] * 100.0 / total;
            g = baseCountG[i] * 100.0 / total;
            t = baseCountT[i] * 100.0 / total;
            n = baseCountN[i] * 100.0 / total;
        }

        out
            << (i + 1)
            << "\t"
            << a
            << "\t"
            << c
            << "\t"
            << g
            << "\t"
            << t
            << "\t"
            << n
            << "\n";
            }
}

void writeQualityDistributionTsv(
    const std::vector<uint64_t>& distribution,
    const std::string& outDir,
    const std::string& readName)
{
    std::string path =
        outDir + "/quality_distribution_" + readName + ".tsv";

    std::ofstream out(path);

    if (!out)
        throw std::runtime_error("Cannot write to " + path);

    out << "quality\tcount\n";

    for (size_t q = 0; q < distribution.size(); ++q)
    {
        out << q
            << "\t"
            << distribution[q]
            << "\n";
    }
}

void writePerSequenceQualityTsv(
    const std::vector<uint64_t>& distribution,
    const std::string& outDir,
    const std::string& readName)
{
    const std::string path = outDir + "/per_sequence_quality_" + readName + ".tsv";
    std::ofstream out(path);
    if (!out) throw std::runtime_error("Cannot write to " + path);

    out << "mean_quality\tread_count\n";
    for (size_t quality = 0; quality < distribution.size(); ++quality) {
        out << quality << "\t" << distribution[quality] << "\n";
    }
}

void writeAdapterTsv(const std::vector<QualityAnalyzer::Adapter>& adapters,
                     const std::vector<std::vector<uint64_t>>& adapterPosCounts,
                     size_t totalReads,
                     size_t maxLength,
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
AnalysisResult processOneFile(const std::string& path,
                              const std::string& readName,
                              const std::string& outDir,
                              const std::string& sampleId,
                              bool skipAdapters,
                              bool collectTimings) {
    AnalysisResult result;
    PerformanceTimers& timers = result.timers;
    QualityAnalyzer analyzer;

    {
        const auto openStart = collectTimings ? Clock::now() : Clock::time_point{};
        FastqReader reader(path, collectTimings);
        if (collectTimings) timers.fileOpen += Clock::now() - openStart;
        FastqRecord rec;
        size_t count = 0;

        while (reader.readNext(rec)) {
            {
                ScopedTimer timer(timers.metrics, collectTimings);
                analyzer.processRecord(rec);
            }
            if (!skipAdapters) {
                ScopedTimer timer(timers.adapterSearch, collectTimings);
                analyzer.analyzeAdapters(rec);
            }
            ++count;
            if (count % 1'000'000 == 0) {
                std::cout << readName << ": processed " << count << " reads...\n";
            }
        }
        std::cout << readName << ": total reads = " << count << "\n";
        if (collectTimings) {
            timers.readAndDecompress += reader.getTiming().readAndDecompress;
            timers.validation += reader.getTiming().validation;
        }
    }

    QualityStats stats;
    {
        ScopedTimer timer(timers.metrics, collectTimings);
        stats = analyzer.getStats();
    }

    printConsoleSummary(stats, readName, skipAdapters);
    {
        ScopedTimer timer(timers.reportWriting, collectTimings);
        writeSummaryTxt(stats, outDir, sampleId + "_" + readName, skipAdapters);
        writePerCycleQualityTsv(
            stats.meanQualityPerPosition,
            outDir,
            readName);
        writeQualityDistributionTsv(
            stats.qualityDistribution,
            outDir,
            readName);
        writePerSequenceQualityTsv(
            stats.perSequenceQualityDistribution,
            outDir,
            readName);

        writePerBaseSequenceContentTsv(
            stats.baseCountA,
            stats.baseCountC,
            stats.baseCountG,
            stats.baseCountT,
            stats.baseCountN,
            outDir,
            readName);
        if (!skipAdapters) {
            writeAdapterTsv(analyzer.adapters,
                            analyzer.adapterPosCounts,
                            stats.totalReads,
                            stats.maxLength,
                            outDir,
                            readName);
        }
    }

    result.r1Stats = stats;
    return result;
}

// ---------------------------------------------------------------------------
// Обработка парных файлов (R1 и R2)
// ---------------------------------------------------------------------------

AnalysisResult processPairedFiles(const std::string& r1Path,
                                     const std::string& r2Path,
                                     const std::string& outDir,
                                     const std::string& sampleId,
                                     bool skipAdapters,
                                     bool collectTimings)
{
    AnalysisResult result;
    PerformanceTimers& timers = result.timers;
    const auto openStart = collectTimings ? Clock::now() : Clock::time_point{};
    FastqReader readerR1(r1Path, collectTimings);
    FastqReader readerR2(r2Path, collectTimings);
    if (collectTimings) timers.fileOpen += Clock::now() - openStart;

    QualityAnalyzer analyzerR1(ReadDirection::R1);
    QualityAnalyzer analyzerR2(ReadDirection::R2);

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
        
        // Проверка идентификаторов считываний
        bool matchingReadIds;
        {
            ScopedTimer timer(timers.pairValidation, collectTimings);
            matchingReadIds = normalizeReadId(rec1.header) == normalizeReadId(rec2.header);
        }
        if (!matchingReadIds)
        {
            throw std::runtime_error(
                "FASTQ validation error:\n"
                "reason: paired read identifiers do not match\n"
                "R1: " + rec1.header + "\n"
                "R2: " + rec2.header);
        }

        {
            ScopedTimer timer(timers.metrics, collectTimings);
            analyzerR1.processRecord(rec1);
            analyzerR2.processRecord(rec2);
        }
        if (!skipAdapters) {
            ScopedTimer timer(timers.adapterSearch, collectTimings);
            analyzerR1.analyzeAdapters(rec1);
            analyzerR2.analyzeAdapters(rec2);
        }

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

    if (collectTimings) {
        timers.readAndDecompress += readerR1.getTiming().readAndDecompress;
        timers.readAndDecompress += readerR2.getTiming().readAndDecompress;
        timers.validation += readerR1.getTiming().validation;
        timers.validation += readerR2.getTiming().validation;
    }

    QualityStats statsR1;
    QualityStats statsR2;
    {
        ScopedTimer timer(timers.metrics, collectTimings);
        statsR1 = analyzerR1.getStats();
        statsR2 = analyzerR2.getStats();
    }

    printConsoleSummary(statsR1, "R1", skipAdapters);
    printConsoleSummary(statsR2, "R2", skipAdapters);

    {
        ScopedTimer timer(timers.reportWriting, collectTimings);
        writeSummaryTxt(statsR1, outDir, sampleId + "_R1", skipAdapters);
        writeSummaryTxt(statsR2, outDir, sampleId + "_R2", skipAdapters);

        writePerCycleQualityTsv(statsR1.meanQualityPerPosition, outDir, "R1");
        writePerCycleQualityTsv(statsR2.meanQualityPerPosition, outDir, "R2");

        writeQualityDistributionTsv(statsR1.qualityDistribution, outDir, "R1");
        writeQualityDistributionTsv(statsR2.qualityDistribution, outDir, "R2");

        writePerBaseSequenceContentTsv(statsR1.baseCountA, statsR1.baseCountC, statsR1.baseCountG, statsR1.baseCountT, statsR1.baseCountN, outDir, "R1");
        writePerBaseSequenceContentTsv(statsR2.baseCountA, statsR2.baseCountC, statsR2.baseCountG, statsR2.baseCountT, statsR2.baseCountN, outDir, "R2");

        writePerSequenceQualityTsv(statsR1.perSequenceQualityDistribution, outDir, "R1");
        writePerSequenceQualityTsv(statsR2.perSequenceQualityDistribution, outDir, "R2");

        if (!skipAdapters) {
            writeAdapterTsv(analyzerR1.adapters, analyzerR1.adapterPosCounts,
                            statsR1.totalReads, statsR1.maxLength, outDir, "R1");
            writeAdapterTsv(analyzerR2.adapters, analyzerR2.adapterPosCounts,
                            statsR2.totalReads, statsR2.maxLength, outDir, "R2");
        }
    }

    result.r1Stats = statsR1;
    result.r2Stats = statsR2;
    return result;
}

std::string jsonEscape(const std::string& value) {
    std::string escaped;
    escaped.reserve(value.size());
    for (const char ch : value) {
        switch (ch) {
            case '"': escaped += "\\\""; break;
            case '\\': escaped += "\\\\"; break;
            case '\n': escaped += "\\n"; break;
            case '\r': escaped += "\\r"; break;
            case '\t': escaped += "\\t"; break;
            default: escaped += ch; break;
        }
    }
    return escaped;
}

std::string currentUtcTimestamp() {
    const std::time_t now = std::time(nullptr);
    std::tm timeInfo{};
#ifdef _WIN32
    gmtime_s(&timeInfo, &now);
#else
    gmtime_r(&now, &timeInfo);
#endif
    std::ostringstream output;
    output << std::put_time(&timeInfo, "%Y-%m-%dT%H:%M:%SZ");
    return output.str();
}

void writeCaseSummary(const std::string& patientId,
                      const fs::path& caseOutputDir,
                      const std::vector<BatchSampleResult>& results,
                      const std::vector<std::string>& warnings)
{
    bool casePassed = true;
    for (const auto& result : results) {
        if (result.entry.patientId == patientId && !result.passed) casePassed = false;
    }

    const fs::path path = caseOutputDir / "case_summary.json";
    std::ofstream output(path);
    if (!output) throw std::runtime_error("Cannot write case summary: " + path.string());

    output << "{\n"
           << "  \"patient_id\": \"" << jsonEscape(patientId) << "\",\n"
           << "  \"status\": \"" << (casePassed ? "passed" : "failed") << "\",\n"
           << "  \"neoqc_version\": \"0.1\",\n"
           << "  \"run_date\": \"" << currentUtcTimestamp() << "\",\n"
           << "  \"warnings\": [";

    bool firstWarning = true;
    for (const auto& warning : warnings) {
        if (warning.find("Patient " + patientId + ":") == std::string::npos) continue;
        if (!firstWarning) output << ", ";
        output << "\"" << jsonEscape(warning) << "\"";
        firstWarning = false;
    }
    output << "],\n  \"samples\": [\n";

    bool firstSample = true;
    for (const auto& result : results) {
        if (result.entry.patientId != patientId) continue;
        if (!firstSample) output << ",\n";
        firstSample = false;

        const auto& entry = result.entry;
        output << "    {\n"
               << "      \"sample_id\": \"" << jsonEscape(entry.sampleId) << "\",\n"
               << "      \"role\": \"" << jsonEscape(entry.sampleRole) << "\",\n"
               << "      \"material\": \"" << jsonEscape(entry.material) << "\",\n"
               << "      \"r1\": \"" << jsonEscape(entry.r1) << "\",\n"
               << "      \"r2\": \"" << jsonEscape(entry.r2) << "\",\n"
               << "      \"qc_status\": \"" << (result.passed ? "passed" : "failed") << "\"";

        if (result.passed && result.analysis) {
            output << ",\n      \"r1_reads\": " << result.analysis->r1Stats.totalReads;
            if (result.analysis->r2Stats) {
                output << ",\n      \"r2_reads\": " << result.analysis->r2Stats->totalReads;
            }
        } else {
            output << ",\n      \"qc_error\": \"" << jsonEscape(result.error) << "\"";
        }
        output << "\n    }";
    }
    output << "\n  ]\n}\n";
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

    if (!args.samples.empty()) {
        try {
            const auto entries = loadAndValidateSampleSheet(args.samples);
            const auto warnings = validateCaseComposition(entries);
            std::cout << "Sample sheet is valid: " << args.samples << "\n"
                      << "Samples: " << entries.size() << "\n";
            for (const auto& warning : warnings) {
                std::cout << "Warning: " << warning << "\n";
            }

            if (args.outDir.empty()) return 0;

            std::error_code ec;
            fs::create_directories(args.outDir, ec);
            if (ec) {
                throw std::runtime_error("Cannot create output directory '" + args.outDir
                                         + "': " + ec.message());
            }

            std::vector<BatchSampleResult> results;
            bool allPassed = true;
            for (const auto& entry : entries) {
                const fs::path sampleOutDir = fs::path(args.outDir)
                                              / entry.patientId / entry.sampleId;
                fs::create_directories(sampleOutDir, ec);
                if (ec) {
                    throw std::runtime_error("Cannot create output directory '"
                                             + sampleOutDir.string() + "': " + ec.message());
                }

                std::cout << "\nCase " << entry.patientId << " — checking " << entry.sampleId
                          << "\nR1: " << entry.r1 << "\n";
                if (!entry.r2.empty()) std::cout << "R2: " << entry.r2 << "\n";

                try {
                    AnalysisResult analysis;
                    if (entry.r2.empty()) {
                        analysis = processOneFile(entry.r1, "R1", sampleOutDir.string(),
                                                  entry.sampleId, args.skipAdapters, args.timings);
                    } else {
                        analysis = processPairedFiles(entry.r1, entry.r2, sampleOutDir.string(),
                                                      entry.sampleId, args.skipAdapters, args.timings);
                    }

                    if (args.plot && !args.skipAdapters) {
                        ScopedTimer timer(analysis.timers.plotting, args.timings);
                        const std::string plotDir = (sampleOutDir / "plots").string();
                        PlotRunner::run((sampleOutDir / "adapter_content_R1.tsv").string(), plotDir);
                        if (!entry.r2.empty()) {
                            PlotRunner::run((sampleOutDir / "adapter_content_R2.tsv").string(), plotDir);
                        }
                    }
                    if (args.timings) printPerformanceTimers(analysis.timers);
                    std::cout << "Result: passed\n";
                    results.push_back({entry, true, std::move(analysis), ""});
                } catch (const std::exception& e) {
                    std::cerr << "Result: failed for " << entry.sampleId << ": " << e.what() << "\n";
                    allPassed = false;
                    results.push_back({entry, false, std::nullopt, e.what()});
                }
            }

            for (const auto& entry : entries) {
                const fs::path caseOutDir = fs::path(args.outDir) / entry.patientId;
                if (!fs::exists(caseOutDir / "case_summary.json")) {
                    writeCaseSummary(entry.patientId, caseOutDir, results, warnings);
                    std::cout << "Case summary: " << (caseOutDir / "case_summary.json") << "\n";
                }
            }

            if (!allPassed) {
                std::cerr << "One or more samples failed; the case is not successful.\n";
                return 1;
            }

            std::cout << "\nAll samples passed.\n";
            return 0;
        } catch (const std::exception& e) {
            std::cerr << e.what() << '\n';
            return 1;
        }
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
    if (args.skipAdapters) std::cout << "Adapters  : skipped\n";

    PerformanceTimers timers;
    try
    {
        if (isPaired)
        {
            timers = processPairedFiles(
                args.r1,
                args.r2,
                args.outDir,
                args.sampleId,
                args.skipAdapters,
                args.timings).timers;
        }
        else
        {
            timers = processOneFile(
                args.r1,
                "R1",
                args.outDir,
                args.sampleId,
                args.skipAdapters,
                args.timings).timers;
        }
    }
    catch (const std::exception& e)
    {
        std::cerr << e.what() << '\n';
        return 1;
    }

    // Построение графиков (опционально, через PlotRunner)
    if (args.plot && !args.skipAdapters) {
        ScopedTimer timer(timers.plotting, args.timings);
        std::string plotDir = args.outDir + "/plots";
        PlotRunner::run(args.outDir + "/adapter_content_R1.tsv", plotDir);
        if (isPaired) {
            PlotRunner::run(args.outDir + "/adapter_content_R2.tsv", plotDir);
        }
    } else if (args.plot) {
        std::cout << "Skipping adapter plots because adapter search was disabled.\n";
    }

    if (args.timings) printPerformanceTimers(timers);
    std::cout << "\nDone.\n";
    return 0;
}
