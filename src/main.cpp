#include "../include/fastq_reader.h"
#include "../include/quality_analyzer.h"
#include "../include/plot_runner.h"
#include "../include/sample_sheet.h"

#include <omp.h>
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
#include <algorithm>
#include <future>

namespace fs = std::filesystem;
using Clock = std::chrono::steady_clock;

struct AnalysisResult {
    QualityStats r1Stats;
    std::optional<QualityStats> r2Stats;
};

struct BatchSampleResult {
    SampleSheetEntry entry;
    bool passed = false;
    std::optional<AnalysisResult> analysis;
    std::string error;
};

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
        "    " << progName << " --r1 <file> --sample-id <id> --out <dir> [--plot] [--skip-adapters]\n\n"
        "  paired-end:\n"
        "    " << progName << " --r1 <file> --r2 <file> --sample-id <id> --out <dir> [--plot] [--skip-adapters]\n\n"
        "  validate sample sheet:\n"
        "    " << progName << " --samples <samples.csv> [--out <dir>] [--plot] [--skip-adapters]\n\n"
        "Options:\n"
        "  --r1 <file>       Path to R1 FASTQ (plain or .gz)\n"
        "  --r2 <file>       Path to R2 FASTQ (optional, for paired-end)\n"
        "  --sample-id <id>  Sample identifier (used in output filenames)\n"
        "  --out <dir>       Output directory (created if missing); enables batch QC with --samples\n"
        "  --samples <file>  Validate a CSV table; combine with --out to run batch QC\n"
        "  --plot            Build plots via plot_results.py (optional)\n"
        "  --skip-adapters   Disable adapter search (for performance measurements)\n";
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

BaseValidationError processRecord(
    QualityAnalyzer& analyzer,
    const FastqRecord& record,
    bool skipAdapters)
{
    BaseValidationError error =
        analyzer.processRecord(record);

    if (!error.found && !skipAdapters)
    {
        analyzer.analyzeAdapters(record);
    }

    return error;
}

void processBatchParallel(
    std::vector<QualityAnalyzer>& analyzers,
    const std::vector<FastqRecord>& batch,
    bool skipAdapters)
{
    std::vector<BaseValidationError> validationErrors(
        analyzers.size());

    #pragma omp parallel for schedule(static)
    for (int i = 0; i < static_cast<int>(batch.size()); ++i)
    {
        const int threadId = omp_get_thread_num();

        BaseValidationError error = processRecord(
            analyzers[threadId],
            batch[i],
            skipAdapters);

        if (error.found)
        {
            validationErrors[threadId] = error;
        }
    }

    for (const auto& error : validationErrors)
    {
        if (error.found)
        {
            std::ostringstream oss;

            oss << "FASTQ validation error:\n"
                << "record: " << error.recordNumber
                << "\n"
                << "reason: invalid base '"
                << error.base
                << "' at position "
                << (error.position + 1);

            throw std::runtime_error(oss.str());
        }
    }
}

bool readPairedBatch(
    FastqReader& readerR1,
    FastqReader& readerR2,
    std::vector<FastqRecord>& batchR1,
    std::vector<FastqRecord>& batchR2,
    std::size_t batchSize)
{
    batchR1.clear();
    batchR2.clear();

    if (batchR1.capacity() < batchSize)
        batchR1.reserve(batchSize);

    if (batchR2.capacity() < batchSize)
        batchR2.reserve(batchSize);

    FastqRecord rec1;
    FastqRecord rec2;

    while (batchR1.size() < batchSize)
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

        // Проверяем идентификаторы
        if (normalizeReadId(rec1.header) != normalizeReadId(rec2.header))
        {
            throw std::runtime_error(
                "FASTQ validation error:\n"
                "reason: paired read identifiers do not match\n"
                "R1: " + rec1.header + "\n"
                "R2: " + rec2.header);
        }

        batchR1.emplace_back(std::move(rec1));
        batchR2.emplace_back(std::move(rec2));

        rec1 = FastqRecord{};
        rec2 = FastqRecord{};
    }

    return !batchR1.empty();
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
                             const std::vector<double>& lowerQuartile,
                             const std::vector<double>& median,
                             const std::string& outDir,
                             const std::string& readName)
{
    std::string path = outDir + "/per_cycle_" + readName + ".tsv";

    std::ofstream out(path);

    if (!out)
        throw std::runtime_error("Cannot write to " + path);

    if (lowerQuartile.size() != meanQuality.size() ||
        median.size() != meanQuality.size()) {
        throw std::runtime_error("Per-cycle quality vectors have different sizes");
    }

    out << "cycle\tmean_quality\tlower_quartile\tmedian\n";

    for (size_t i = 0; i < meanQuality.size(); ++i)
    {
        out << (i + 1)
            << "\t"
            << std::fixed
            << std::setprecision(4)
            << meanQuality[i]
            << "\t"
            << lowerQuartile[i]
            << "\t"
            << median[i]
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

void writePerSequenceGCContentTsv(
    const std::vector<uint64_t>& gcDistribution,
    const std::string& outDir,
    const std::string& readName)
{
    std::string path =
        outDir + "/per_sequence_gc_content_" + readName + ".tsv";

    std::ofstream out(path);

    if (!out)
        throw std::runtime_error("Cannot write to " + path);

    out << "gc_percent\treads\n";

    for (size_t i = 0; i < gcDistribution.size(); ++i)
    {
        out << i
            << "\t"
            << gcDistribution[i]
            << "\n";
    }
}

void writePerBaseNContentTsv(
    const std::vector<uint64_t>& baseCountN,
    const std::vector<uint64_t>& readsPerPosition,
    const std::string& outDir,
    const std::string& readName)
{
    std::string path =
        outDir + "/per_base_n_content_" + readName + ".tsv";

    std::ofstream out(path);

    if (!out)
        throw std::runtime_error("Cannot write to " + path);

    out << "position\tN_percent\n";

    for (size_t i = 0; i < baseCountN.size(); ++i)
    {
        double percent = 0.0;

        if (readsPerPosition[i] > 0)
        {
            percent =
                static_cast<double>(baseCountN[i]) * 100.0 /
                readsPerPosition[i];
        }

        out
            << (i + 1)
            << "\t"
            << std::fixed
            << std::setprecision(4)
            << percent
            << "\n";
    }
}

void writePerSequenceQualityTsv(
    const std::vector<uint64_t>& distribution,
    const std::vector<uint64_t>& distributionTruncate,
    const std::string& outDir,
    const std::string& readName)
{
    const std::string path = outDir + "/per_sequence_quality_" + readName + ".tsv";
    std::ofstream out(path);
    if (!out) throw std::runtime_error("Cannot write to " + path);

    out << "mean_quality\tread_count\tread_count_truncate\n";
   const size_t maxSize = std::max(distribution.size(), distributionTruncate.size());

    for (size_t quality = 0; quality < maxSize; ++quality) {
        const uint64_t rounded =
            quality < distribution.size()
                ? distribution[quality]
                : 0;

        const uint64_t truncated =
            quality < distributionTruncate.size()
                ? distributionTruncate[quality]
                : 0;

        out << quality
            << "\t"
            << rounded
            << "\t"
            << truncated
            << "\n";
    }
}

void removeRetiredQualityDistributionArtifacts(
    const std::string& outDir,
    const std::string& readName)
{
    const fs::path retiredPath =
        fs::path(outDir) / ("quality_distribution_" + readName + ".tsv");
    std::error_code ec;
    fs::remove(retiredPath, ec);
    if (ec) {
        throw std::runtime_error(
            "Cannot remove retired artifact '" + retiredPath.string() + "': " + ec.message());
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

void writeSequenceLengthDistributionTsv(
    const std::vector<uint64_t>& lengthDistribution,
    const std::string& outDir,
    const std::string& readName)
{
    std::string path =
        outDir + "/sequence_length_distribution_" + readName + ".tsv";

    std::ofstream out(path);

    if (!out)
        throw std::runtime_error("Cannot write to " + path);

    out << "length\treads\n";

    for (size_t i = 0; i < lengthDistribution.size(); ++i)
    {
        if (lengthDistribution[i] == 0)
            continue;

        out << i
            << "\t"
            << lengthDistribution[i]
            << "\n";
    }
}

template <typename Writer>
void writeAtomically(const fs::path& path, Writer writer) {
    fs::path temporary = path;
    temporary += ".tmp";
    {
        std::ofstream out(temporary, std::ios::trunc);
        if (!out) throw std::runtime_error("Cannot write to " + temporary.string());
        writer(out);
        out.flush();
        if (!out) {
            std::error_code cleanupError;
            fs::remove(temporary, cleanupError);
            throw std::runtime_error("Cannot complete write to " + temporary.string());
        }
    }

    std::error_code error;
    fs::rename(temporary, path, error);
    if (error) {
        std::error_code removeError;
        fs::remove(path, removeError);
        error.clear();
        fs::rename(temporary, path, error);
    }
    if (error) {
        std::error_code cleanupError;
        fs::remove(temporary, cleanupError);
        throw std::runtime_error("Cannot publish " + path.string() + ": " + error.message());
    }
}

std::string tsvSafeFilename(const std::string& path) {
    std::string filename = fs::path(path).filename().string();
    std::replace(filename.begin(), filename.end(), '\t', '_');
    std::replace(filename.begin(), filename.end(), '\n', '_');
    std::replace(filename.begin(), filename.end(), '\r', '_');
    return filename;
}

fs::path duplicationIncompletePath(const std::string& outDir,
                                   const std::string& readName) {
    return fs::path(outDir) / ("sequence_duplication_" + readName + ".incomplete");
}

void beginDuplicationArtifacts(const std::string& outDir,
                               const std::string& readName,
                               const std::string& sourceFastq) {
    for (const char* prefix : {
             "sequence_duplication_levels_",
             "sequence_duplication_summary_",
             "overrepresented_sequences_",
         }) {
        const fs::path path = fs::path(outDir) / (prefix + readName + ".tsv");
        std::error_code error;
        fs::remove(path, error);
        if (error) {
            throw std::runtime_error("Cannot remove stale artifact " + path.string()
                                     + ": " + error.message());
        }
        fs::path temporary = path;
        temporary += ".tmp";
        error.clear();
        fs::remove(temporary, error);
        if (error) {
            throw std::runtime_error("Cannot remove stale artifact " + temporary.string()
                                     + ": " + error.message());
        }
    }

    const fs::path marker = duplicationIncompletePath(outDir, readName);
    std::error_code error;
    fs::remove(marker, error);
    if (error) {
        throw std::runtime_error("Cannot remove stale marker " + marker.string()
                                 + ": " + error.message());
    }
    fs::path temporaryMarker = marker;
    temporaryMarker += ".tmp";
    error.clear();
    fs::remove(temporaryMarker, error);
    if (error) {
        throw std::runtime_error("Cannot remove stale marker " + temporaryMarker.string()
                                 + ": " + error.message());
    }
    writeAtomically(marker, [&](std::ostream& out) {
        out << "Duplication artifacts are incomplete for "
            << tsvSafeFilename(sourceFastq) << '\n';
    });
}

void writeDuplicationArtifacts(const DuplicationStats& stats,
                               const std::string& sourceFastq,
                               const std::string& outDir,
                               const std::string& readName) {
    const fs::path root(outDir);
    writeAtomically(root / ("sequence_duplication_levels_" + readName + ".tsv"),
        [&](std::ostream& out) {
            out << "duplication_level\ttotal_sequences_percent"
                   "\tdeduplicated_sequences_percent\n";
            out << std::fixed << std::setprecision(10);
            for (const auto& row : stats.levels) {
                out << row.label << '\t'
                    << row.totalSequencesPercent << '\t'
                    << row.deduplicatedSequencesPercent << '\n';
            }
        });

    writeAtomically(root / ("overrepresented_sequences_" + readName + ".tsv"),
        [&](std::ostream& out) {
            out << "sequence\tcount\tpercentage\tpossible_source\n";
            out << std::fixed << std::setprecision(10);
            for (const auto& sequence : stats.overrepresentedSequences) {
                out << sequence.sequence << '\t'
                    << sequence.count << '\t'
                    << sequence.percent << "\tNo Hit\n";
            }
        });

    // The summary is the transaction's provenance record and is published last.
    writeAtomically(root / ("sequence_duplication_summary_" + readName + ".tsv"),
        [&](std::ostream& out) {
            out << "source_kind\talgorithm\tsource_fastq\tprefix_length"
                   "\ttotal_reads\tunique_sequences"
                   "\tdeduplicated_remaining_percent\n";
            out << "native_fastq\tneoqc-exact-prefix-v1\t"
                << tsvSafeFilename(sourceFastq) << '\t'
                << DUPLICATION_PREFIX_LENGTH << '\t'
                << stats.totalReads << '\t'
                << stats.uniqueSequences << '\t'
                << std::fixed << std::setprecision(10)
                << stats.deduplicatedRemainingPercent << '\n';
        });

    const fs::path marker = duplicationIncompletePath(outDir, readName);
    std::error_code error;
    if (!fs::remove(marker, error) || error) {
        throw std::runtime_error("Cannot complete duplication transaction "
                                 + marker.string() + ": "
                                 + (error ? error.message() : "marker is missing"));
    }
}

void writeAnalysisReports(
    const QualityStats& stats,
    const DuplicationStats& duplicationStats,
    const QualityAnalyzer& analyzer,
    const std::string& sourcePath,
    const std::string& outDir,
    const std::string& sampleId,
    const std::string& readName,
    bool skipAdapters)
{
    writeSummaryTxt(
        stats,
        outDir,
        sampleId + "_" + readName,
        skipAdapters);

    writePerCycleQualityTsv(
        stats.meanQualityPerPosition,
        stats.lowerQuartileQualityPerPosition,
        stats.medianQualityPerPosition,
        outDir,
        readName);

    writePerSequenceQualityTsv(
        stats.perSequenceQualityDistribution,
        stats.perSequenceQualityDistributionTruncate,
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

    writePerSequenceGCContentTsv(
        stats.gcDistribution,
        outDir,
        readName);

    writePerBaseNContentTsv(
        stats.baseCountN,
        stats.readsPerPosition,
        outDir,
        readName);

    writeSequenceLengthDistributionTsv(
        stats.lengthDistribution,
        outDir,
        readName);

    writeDuplicationArtifacts(
        duplicationStats,
        sourcePath,
        outDir,
        readName);

    if (!skipAdapters)
    {
        writeAdapterTsv(
            analyzer.adapters,
            analyzer.adapterPosCounts,
            stats.totalReads,
            stats.maxLength,
            outDir,
            readName);
    }
}

std::vector<DuplicationEntry> mergeSortedDuplicationEntries(
    const std::vector<DuplicationEntry>& a,
    const std::vector<DuplicationEntry>& b)
{
    std::vector<DuplicationEntry> result;
    result.reserve(a.size() + b.size());

    std::size_t i = 0;
    std::size_t j = 0;

    while (i < a.size() && j < b.size())
    {
        if (a[i].key.words < b[j].key.words)
        {
            result.push_back(a[i]);
            ++i;
        }
        else if (b[j].key.words < a[i].key.words)
        {
            result.push_back(b[j]);
            ++j;
        }
        else
        {
            result.push_back({
                a[i].key,
                a[i].count + b[j].count
            });

            ++i;
            ++j;
        }
    }

    while (i < a.size())
    {
        result.push_back(a[i]);
        ++i;
    }

    while (j < b.size())
    {
        result.push_back(b[j]);
        ++j;
    }

    return result;
}

std::vector<DuplicationEntry> mergeDuplicationEntriesTree(
    std::vector<std::vector<DuplicationEntry>> entries)
{
    if (entries.empty())
        return {};

    while (entries.size() > 1)
    {
        std::vector<std::vector<DuplicationEntry>> next;
        next.reserve((entries.size() + 1) / 2);

        for (std::size_t i = 0; i < entries.size(); i += 2)
        {
            if (i + 1 < entries.size())
            {
                next.push_back(
                    mergeSortedDuplicationEntries(
                        entries[i],
                        entries[i + 1]));
            }
            else
            {
                next.push_back(std::move(entries[i]));
            }
        }

        entries = std::move(next);
    }

    return std::move(entries[0]);
}

// ---------------------------------------------------------------------------
// Обработка одного файла (R1 или R2)
// ---------------------------------------------------------------------------
AnalysisResult processOneFile(const std::string& path,
                            const std::string& readName,
                            const std::string& outDir,
                            const std::string& sampleId,
                            bool skipAdapters) {
    AnalysisResult result;
    removeRetiredQualityDistributionArtifacts(outDir, readName);
    beginDuplicationArtifacts(outDir, readName, path);
    QualityAnalyzer analyzer;

    FastqReader reader(path);

    constexpr std::size_t BATCH_SIZE = 100000;

    const int threadCount = omp_get_max_threads();

    std::cout << "OpenMP threads: "
            << threadCount
            << "\n";

    std::vector<QualityAnalyzer> localAnalyzers;
    localAnalyzers.reserve(threadCount);

    for (int i = 0; i < threadCount; ++i)
    {
        localAnalyzers.emplace_back();
    }

    std::vector<FastqRecord> batch;

    size_t count = 0;

    while (reader.readBatch(batch, BATCH_SIZE))
    {
        processBatchParallel(
            localAnalyzers,
            batch,
            skipAdapters);

        count += batch.size();

        if (count % 1000000 == 0)
        {
            std::cout << "Processed "
                    << count
                    << " reads...\n";
        }
    }

    for (auto& localAnalyzer : localAnalyzers)
    {
        analyzer.merge(localAnalyzer);
    }

    std::vector<std::vector<DuplicationEntry>> entries;

    entries.reserve(localAnalyzers.size());

    for (const auto& localAnalyzer : localAnalyzers)
    {
        entries.push_back(localAnalyzer.getDuplicationEntries());
    }

    for (auto& localEntries : entries)
    {
        std::sort(
            localEntries.begin(),
            localEntries.end(),
            [](const DuplicationEntry& a,
            const DuplicationEntry& b)
            {
                return a.key.words < b.key.words;
            });
    }

    std::vector<DuplicationEntry> mergedEntries = mergeDuplicationEntriesTree(std::move(entries));

    QualityStats stats = analyzer.getStats();

    const DuplicationStats duplicationStats = analyzer.getDuplicationStats(mergedEntries);

    printConsoleSummary(stats, readName, skipAdapters);
    {
    writeAnalysisReports(
        stats,
        duplicationStats,
        analyzer,
        path,
        outDir,
        sampleId,
        readName,
        skipAdapters);
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
                                bool skipAdapters)
{
    AnalysisResult result;
    removeRetiredQualityDistributionArtifacts(outDir, "R1");
    removeRetiredQualityDistributionArtifacts(outDir, "R2");
    beginDuplicationArtifacts(outDir, "R1", r1Path);
    beginDuplicationArtifacts(outDir, "R2", r2Path);
    FastqReader readerR1(r1Path);
    FastqReader readerR2(r2Path);

    QualityAnalyzer analyzerR1(ReadDirection::R1);
    QualityAnalyzer analyzerR2(ReadDirection::R2);

    FastqRecord rec1;
    FastqRecord rec2;

    constexpr std::size_t BATCH_SIZE = 100000;

    const int threadCount = omp_get_max_threads();

    std::cout << "OpenMP threads: "
            << threadCount
            << "\n";

    std::vector<QualityAnalyzer> localAnalyzersR1;
    std::vector<QualityAnalyzer> localAnalyzersR2;

    localAnalyzersR1.reserve(threadCount);
    localAnalyzersR2.reserve(threadCount);

    for (int i = 0; i < threadCount; ++i)
    {
        localAnalyzersR1.emplace_back(ReadDirection::R1);
        localAnalyzersR2.emplace_back(ReadDirection::R2);
    }

    std::vector<FastqRecord> batchR1;
    std::vector<FastqRecord> batchR2;

    size_t count = 0;

    while (readPairedBatch(
        readerR1,
        readerR2,
        batchR1,
        batchR2,
        BATCH_SIZE))
    {
        processBatchParallel(
            localAnalyzersR1,
            batchR1,
            skipAdapters);

        processBatchParallel(
            localAnalyzersR2,
            batchR2,
            skipAdapters);

        count += batchR1.size();

        if (count % 1000000 == 0)
        {
            std::cout << "Processed "
                    << count
                    << " paired reads...\n";
        }
    }

    for (auto& localAnalyzer : localAnalyzersR1)
    {
        analyzerR1.merge(localAnalyzer);
    }

    for (auto& localAnalyzer : localAnalyzersR2)
    {
        analyzerR2.merge(localAnalyzer);
    }

    std::vector<std::vector<DuplicationEntry>> entriesR1(
        localAnalyzersR1.size());

    std::vector<std::vector<DuplicationEntry>> entriesR2(
        localAnalyzersR2.size());

    #pragma omp parallel for schedule(static)
    for (int i = 0;
        i < static_cast<int>(localAnalyzersR1.size());
        ++i)
    {
        entriesR1[i] =
            localAnalyzersR1[i].getDuplicationEntries();
    }

    #pragma omp parallel for schedule(static)
    for (int i = 0;
        i < static_cast<int>(localAnalyzersR2.size());
        ++i)
    {
        entriesR2[i] =
            localAnalyzersR2[i].getDuplicationEntries();
    }

    #pragma omp parallel for schedule(static)
    for (int i = 0;
        i < static_cast<int>(entriesR1.size());
        ++i)
    {
        auto& entries = entriesR1[i];

        std::sort(
            entries.begin(),
            entries.end(),
            [](const DuplicationEntry& a,
            const DuplicationEntry& b)
            {
                return a.key.words < b.key.words;
            });
    }

    #pragma omp parallel for schedule(static)
    for (int i = 0;
        i < static_cast<int>(entriesR2.size());
        ++i)
    {
        auto& entries = entriesR2[i];

        std::sort(
            entries.begin(),
            entries.end(),
            [](const DuplicationEntry& a,
            const DuplicationEntry& b)
            {
                return a.key.words < b.key.words;
            });
    }

    std::vector<DuplicationEntry> mergedR1 = mergeDuplicationEntriesTree(std::move(entriesR1));

    std::vector<DuplicationEntry> mergedR2 = mergeDuplicationEntriesTree(std::move(entriesR2));

    std::cout << "R1: total reads = "
              << analyzerR1.getTotalReads()
              << "\n";

    std::cout << "R2: total reads = "
              << analyzerR2.getTotalReads()
              << "\n";

    QualityStats statsR1 = analyzerR1.getStats();
    QualityStats statsR2 = analyzerR2.getStats();
        
    const DuplicationStats duplicationStatsR1 = analyzerR1.getDuplicationStats(mergedR1);
    const DuplicationStats duplicationStatsR2 = analyzerR2.getDuplicationStats(mergedR2);

    printConsoleSummary(statsR1, "R1", skipAdapters);
    printConsoleSummary(statsR2, "R2", skipAdapters);

    auto futureR1 = std::async(
        std::launch::async,
        [&]()
        {
            writeAnalysisReports(
                statsR1,
                duplicationStatsR1,
                analyzerR1,
                r1Path,
                outDir,
                sampleId,
                "R1",
                skipAdapters);
        });

    auto futureR2 = std::async(
        std::launch::async,
        [&]()
        {
            writeAnalysisReports(
                statsR2,
                duplicationStatsR2,
                analyzerR2,
                r2Path,
                outDir,
                sampleId,
                "R2",
                skipAdapters);
        });

    futureR1.get();
    futureR2.get();

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

    const auto start = Clock::now();

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
                        analysis = processOneFile(
                            entry.r1,
                            "R1",
                            sampleOutDir.string(),
                            entry.sampleId,
                            args.skipAdapters);
                    } else {
                        analysis = processPairedFiles(
                            entry.r1,
                            entry.r2,
                            sampleOutDir.string(),
                            entry.sampleId,
                            args.skipAdapters);
                    }

                    if (args.plot) {
                        const std::string plotDir = (sampleOutDir / "plots").string();
                        PlotOptions plotOptions;
                        plotOptions.includeAdapters = !args.skipAdapters;
                        PlotRunner::runAll(sampleOutDir.string(), plotDir, plotOptions);
                    }
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

    try
    {
        if (isPaired)
        {
            processPairedFiles(
                args.r1,
                args.r2,
                args.outDir,
                args.sampleId,
                args.skipAdapters);
        }
        else
        {
            processOneFile(
                args.r1,
                "R1",
                args.outDir,
                args.sampleId,
                args.skipAdapters);
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
        PlotOptions plotOptions;
        plotOptions.includeAdapters = !args.skipAdapters;
        PlotRunner::runAll(args.outDir, plotDir, plotOptions);
    }

    const auto end = Clock::now();

    const auto elapsed = std::chrono::duration<double>(end - start);

    std::cout << "Done.\n"
              << "Total processing time: "
              << elapsed.count()
              << " s\n";

    return 0;
}
