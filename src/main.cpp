#include "../include/fastq_reader.h"
#include "../include/quality_analyzer.h"
#include "../include/plot_runner.h"
#include "../include/sample_sheet.h"
#include "../include/trimming/trim_config.h"
#include "../include/trimming/adapter_fasta.h"
#include "../include/trimming/fastq_writer.h"
#include "../include/trimming/trim_report.h"
#include "../include/trimming/trimmer.h"

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
#include <random>
#include <sstream>
#include <charconv>
#include <cctype>
#include <system_error>
#include <memory>

#ifndef NEOQC_VERSION
#error "NEOQC_VERSION must be supplied by the build system"
#endif

namespace fs = std::filesystem;
using Clock = std::chrono::steady_clock;

struct AnalysisResult {
    QualityStats r1Stats;
    std::optional<QualityStats> r2Stats;
    std::optional<QualityStats> r1AfterStats;
    std::optional<QualityStats> r2AfterStats;
    TrimStats trimmingStats;
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
    bool        trimmingOptionSpecified = false;
    bool        adapterFastaSpecified = false;
    std::string adapterFasta;
    TrimConfig  trimConfig;
};

struct RunManifest
{
    std::string runId;
    std::string createdAt;
    std::string sampleId;

    bool paired = false;
    bool plot = false;
    bool skipAdapters = false;

    std::string r1Path;
    std::optional<std::string> r2Path;

    std::vector<std::string> artifacts;
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
        auto trimmingValue = [&](const std::string& name) -> std::string {
            args.trimmingOptionSpecified = true;
            return needValue(name);
        };
        auto parseSize = [&](const std::string& name) -> std::size_t {
            const std::string value = trimmingValue(name);
            std::size_t parsed = 0;
            const auto [end, error] = std::from_chars(
                value.data(), value.data() + value.size(), parsed);
            if (error != std::errc{} || end != value.data() + value.size()) {
                throw std::runtime_error(
                    "Option " + name + " requires a non-negative integer");
            }
            return parsed;
        };

        if      (arg == "--r1")        args.r1       = needValue("--r1");
        else if (arg == "--r2")        args.r2       = needValue("--r2");
        else if (arg == "--sample-id") args.sampleId = needValue("--sample-id");
        else if (arg == "--out")       args.outDir   = needValue("--out");
        else if (arg == "--samples")   args.samples  = needValue("--samples");
        else if (arg == "--plot")      args.plot     = true;
        else if (arg == "--skip-adapters") args.skipAdapters = true;
        else if (arg == "--trim")      args.trimConfig.enabled = true;
        else if (arg == "--trim-front") args.trimConfig.trim_front = parseSize(arg);
        else if (arg == "--trim-tail")  args.trimConfig.trim_tail = parseSize(arg);
        else if (arg == "--cut-front") {
            args.trimmingOptionSpecified = true;
            args.trimConfig.cut_front = true;
        }
        else if (arg == "--cut-tail") {
            args.trimmingOptionSpecified = true;
            args.trimConfig.cut_tail = true;
        }
        else if (arg == "--cut-right") {
            args.trimmingOptionSpecified = true;
            args.trimConfig.cut_right = true;
        }
        else if (arg == "--quality-threshold") {
            args.trimConfig.quality_threshold = parseSize(arg);
        }
        else if (arg == "--window-size") args.trimConfig.window_size = parseSize(arg);
        else if (arg == "--adapter-sequence") {
            args.trimConfig.adapter_sequence = trimmingValue(arg);
            args.trimConfig.adapter_trimming = true;
        }
        else if (arg == "--adapter-fasta") {
            args.adapterFastaSpecified = true;
            args.adapterFasta = trimmingValue(arg);
        }
        else if (arg == "--overlap-correction") {
            args.trimmingOptionSpecified = true;
            args.trimConfig.overlap_correction = true;
        }
        else if (arg == "--merge") {
            args.trimmingOptionSpecified = true;
            args.trimConfig.merge_reads = true;
        }
        else if (arg == "--umi") {
            args.trimmingOptionSpecified = true;
            args.trimConfig.umi_enabled = true;
        }
        else if (arg == "--umi-length") {
            args.trimConfig.umi_length = parseSize(arg);
        }
        else if (arg == "--min-length") args.trimConfig.min_length = parseSize(arg);
        else if (arg == "--trim-polyg") {
            args.trimmingOptionSpecified = true;
            args.trimConfig.trim_poly_g = true;
        }
        else if (arg == "--trim-polyx") {
            args.trimmingOptionSpecified = true;
            args.trimConfig.trim_poly_x = true;
        }
        else if (arg == "--polyg-min-length") {
            args.trimConfig.poly_g_min_length = parseSize(arg);
        }
        else if (arg == "--polyx-min-length") {
            args.trimConfig.poly_x_min_length = parseSize(arg);
        }
        else if (arg == "--help" || arg == "-h") {
            throw std::runtime_error("help");
        } else {
            throw std::runtime_error("Unknown argument: " + arg);
        }
    }

    if (args.trimmingOptionSpecified && !args.trimConfig.enabled) {
        throw std::runtime_error("Trimming options require --trim");
    }
    if (args.adapterFastaSpecified) {
        args.trimConfig.fasta_adapters = loadAdapterFasta(args.adapterFasta);
        args.trimConfig.adapter_trimming = true;
    }
    if (args.trimConfig.cut_right && args.trimConfig.window_size == 0) {
        throw std::runtime_error("--cut-right requires --window-size greater than zero");
    }
    if (args.trimConfig.quality_threshold > 93) {
        throw std::runtime_error("--quality-threshold must be between 0 and 93");
    }
    if (!args.trimConfig.adapter_sequence.empty()) {
        const bool valid = std::all_of(
            args.trimConfig.adapter_sequence.begin(),
            args.trimConfig.adapter_sequence.end(),
            [](char base) {
                const char normalized = static_cast<char>(
                    std::toupper(static_cast<unsigned char>(base)));
                return normalized == 'A' || normalized == 'C'
                    || normalized == 'G' || normalized == 'T'
                    || normalized == 'N';
            });
        if (!valid) {
            throw std::runtime_error(
                "--adapter-sequence may contain only A, C, G, T, or N");
        }
        if (args.trimConfig.adapter_sequence.size()
            < args.trimConfig.min_adapter_match) {
            throw std::runtime_error(
                "--adapter-sequence must contain at least "
                + std::to_string(args.trimConfig.min_adapter_match) + " bases");
        }
    }
    if (args.trimConfig.umi_enabled && args.trimConfig.umi_length == 0) {
        throw std::runtime_error("--umi requires --umi-length greater than zero");
    }
    if (!args.trimConfig.umi_enabled && args.trimConfig.umi_length != 0) {
        throw std::runtime_error("--umi-length requires --umi");
    }
    if ((args.trimConfig.adapter_trimming || args.trimConfig.overlap_correction
         || args.trimConfig.merge_reads) && args.trimConfig.min_overlap == 0) {
        throw std::runtime_error("PE overlap minimum must be greater than zero");
    }
    if (args.trimConfig.trim_poly_g && args.trimConfig.poly_g_min_length == 0) {
        throw std::runtime_error(
            "--trim-polyg requires --polyg-min-length greater than zero");
    }
    if (args.trimConfig.trim_poly_x && args.trimConfig.poly_x_min_length == 0) {
        throw std::runtime_error(
            "--trim-polyx requires --polyx-min-length greater than zero");
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
    if ((args.trimConfig.overlap_correction || args.trimConfig.merge_reads)
        && args.r2.empty()) {
        throw std::runtime_error(
            "--overlap-correction and --merge require paired-end input (--r2)");
    }

    return args;
}

void printUsage(const char* progName) {
    std::cerr <<
        "NeoQC — FASTQ quality analysis\n\n"
        "Usage:\n"
        "  single-end:\n"
        "    " << progName << " --r1 <file> --sample-id <id> --out <dir> [--plot] [--skip-adapters] [--trim]\n\n"
        "  paired-end:\n"
        "    " << progName << " --r1 <file> --r2 <file> --sample-id <id> --out <dir> [--plot] [--skip-adapters] [--trim]\n\n"
        "  validate sample sheet:\n"
        "    " << progName << " --samples <samples.csv> [--out <dir>] [--plot] [--skip-adapters] [--trim]\n\n"
        "Options:\n"
        "  --r1 <file>       Path to R1 FASTQ (plain or .gz)\n"
        "  --r2 <file>       Path to R2 FASTQ (optional, for paired-end)\n"
        "  --sample-id <id>  Sample identifier (used in output filenames)\n"
        "  --out <dir>       Output directory (created if missing); enables batch QC with --samples\n"
        "  --samples <file>  Validate a CSV table; combine with --out to run batch QC\n"
        "  --plot            Build plots via plot_results.py (optional)\n"
        "  --skip-adapters   Disable QC adapter search (for performance measurements)\n\n"
        "Trimming options (require --trim):\n"
        "  --trim                 Enable the trimming pipeline, FASTQ output, and before/after QC\n"
        "  --trim-front N         Remove N bases from 5' end (default: 0)\n"
        "  --trim-tail N          Remove N bases from 3' end (default: 0)\n"
        "  --cut-front            Remove low-quality bases from 5' end\n"
        "  --cut-tail             Remove low-quality bases from 3' end\n"
        "  --cut-right            Trim at first low-quality sliding window\n"
        "  --quality-threshold N  Phred threshold (default: 20)\n"
        "  --window-size N        Sliding-window size (default: 4)\n"
        "  --adapter-sequence S   Trim explicit adapter sequence\n"
        "  --adapter-fasta FILE   Trim adapters loaded from FASTA\n"
        "  --overlap-correction   Correct PE overlap mismatches by base quality\n"
        "  --merge                Write overlapping PE consensus reads\n"
        "  --umi                  Extract a 5' UMI from R1 (or the SE read)\n"
        "  --umi-length N         Required UMI length when --umi is enabled\n"
        "  --min-length N         Discard reads shorter than N (default: 0)\n"
        "  --trim-polyg           Trim terminal polyG runs\n"
        "  --trim-polyx           Trim terminal A/C/G/T homopolymer runs\n"
        "  --polyg-min-length N   Minimum polyG run (default: 10)\n"
        "  --polyx-min-length N   Minimum polyX run (default: 10)\n";
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
    std::vector<FastqRecord>& batch,
    bool skipAdapters,
    std::vector<Trimmer>* trimmers = nullptr)
{
    std::vector<BaseValidationError> validationErrors(
        analyzers.size());

    #pragma omp parallel for schedule(static)
    for (int i = 0; i < static_cast<int>(batch.size()); ++i)
    {
        const int threadId = omp_get_thread_num();

        if (trimmers != nullptr)
        {
            (*trimmers)[threadId].trim(batch[i]);
        }

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

void processPairedBatchParallel(
    std::vector<QualityAnalyzer>& analyzersR1,
    std::vector<QualityAnalyzer>& analyzersR2,
    std::vector<FastqRecord>& batchR1,
    std::vector<FastqRecord>& batchR2,
    bool skipAdapters,
    std::vector<Trimmer>* trimmers = nullptr)
{
    std::vector<BaseValidationError> validationErrors(analyzersR1.size());

    #pragma omp parallel for schedule(static)
    for (int i = 0; i < static_cast<int>(batchR1.size()); ++i)
    {
        const int threadId = omp_get_thread_num();

        if (trimmers != nullptr)
        {
            (*trimmers)[threadId].trimPair(batchR1[i], batchR2[i]);
        }

        BaseValidationError error = processRecord(
            analyzersR1[threadId], batchR1[i], skipAdapters);
        if (!error.found)
        {
            error = processRecord(
                analyzersR2[threadId], batchR2[i], skipAdapters);
        }
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
                << "\nreason: invalid base '" << error.base
                << "' at position " << (error.position + 1);
            throw std::runtime_error(oss.str());
        }
    }
}

void processTrimmedBatchParallel(
    std::vector<QualityAnalyzer>& beforeAnalyzers,
    std::vector<QualityAnalyzer>& afterAnalyzers,
    std::vector<Trimmer>& trimmers,
    std::vector<FastqRecord>& batch,
    std::vector<TrimResult>& trimResults,
    bool skipAdapters)
{
    trimResults.resize(batch.size());
    std::vector<BaseValidationError> validationErrors(beforeAnalyzers.size());

    #pragma omp parallel for schedule(static)
    for (int i = 0; i < static_cast<int>(batch.size()); ++i)
    {
        const int threadId = omp_get_thread_num();
        BaseValidationError error = processRecord(
            beforeAnalyzers[threadId], batch[i], skipAdapters);
        if (!error.found)
        {
            trimResults[i] = trimmers[threadId].trim(batch[i]);
            if (trimResults[i].passed)
            {
                error = processRecord(
                    afterAnalyzers[threadId], batch[i], skipAdapters);
            }
        }
        if (error.found)
        {
            validationErrors[threadId] = error;
        }
    }

    for (const auto& error : validationErrors)
    {
        if (error.found)
        {
            std::ostringstream message;
            message << "FASTQ validation error:\n"
                    << "record: " << error.recordNumber
                    << "\nreason: invalid base '" << error.base
                    << "' at position " << (error.position + 1);
            throw std::runtime_error(message.str());
        }
    }
}

void processTrimmedPairedBatchParallel(
    std::vector<QualityAnalyzer>& beforeAnalyzersR1,
    std::vector<QualityAnalyzer>& beforeAnalyzersR2,
    std::vector<QualityAnalyzer>& afterAnalyzersR1,
    std::vector<QualityAnalyzer>& afterAnalyzersR2,
    std::vector<Trimmer>& trimmers,
    std::vector<FastqRecord>& batchR1,
    std::vector<FastqRecord>& batchR2,
    std::vector<TrimResult>& trimResultsR1,
    std::vector<TrimResult>& trimResultsR2,
    std::vector<std::optional<FastqRecord>>& mergedRecords,
    bool skipAdapters)
{
    trimResultsR1.resize(batchR1.size());
    trimResultsR2.resize(batchR2.size());
    mergedRecords.resize(batchR1.size());
    std::vector<BaseValidationError> validationErrors(beforeAnalyzersR1.size());

    #pragma omp parallel for schedule(static)
    for (int i = 0; i < static_cast<int>(batchR1.size()); ++i)
    {
        const int threadId = omp_get_thread_num();
        BaseValidationError error = processRecord(
            beforeAnalyzersR1[threadId], batchR1[i], skipAdapters);
        if (!error.found)
        {
            error = processRecord(
                beforeAnalyzersR2[threadId], batchR2[i], skipAdapters);
        }
        if (!error.found)
        {
            PairTrimResult results =
                trimmers[threadId].trimPairDetailed(batchR1[i], batchR2[i]);
            trimResultsR1[i] = std::move(results.r1);
            trimResultsR2[i] = std::move(results.r2);
            mergedRecords[i] = std::move(results.merged);
            if (trimResultsR1[i].passed && trimResultsR2[i].passed)
            {
                error = processRecord(
                    afterAnalyzersR1[threadId], batchR1[i], skipAdapters);
                if (!error.found)
                {
                    error = processRecord(
                        afterAnalyzersR2[threadId], batchR2[i], skipAdapters);
                }
            }
        }
        if (error.found)
        {
            validationErrors[threadId] = error;
        }
    }

    for (const auto& error : validationErrors)
    {
        if (error.found)
        {
            std::ostringstream message;
            message << "FASTQ validation error:\n"
                    << "record: " << error.recordNumber
                    << "\nreason: invalid base '" << error.base
                    << "' at position " << (error.position + 1);
            throw std::runtime_error(message.str());
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

void writePerCycleQualityTsv(
    const std::string& outDir,
    const std::vector<PerBaseQualityGroup>& groups,
    const std::string& filename)
{
    std::string path = outDir + "/per_cycle_" + filename + ".tsv";
    std::ofstream out(path);

    if (!out) {
        throw std::runtime_error(
            "Failed to open per-cycle quality output: " + path);
    }

    out << "cycle\tmean_quality\tlower_quartile\tmedian\n";

    for (const auto& group : groups) {

        if (group.start == group.end) {
            out << group.start;
        }
        else {
            out << group.start
                << '-'
                << group.end;
        }

        out << '\t';

        if (group.evaluated) {
            out << group.mean << '\t'
                << group.lowerQuartile << '\t'
                << group.median;
        }
        else {
            /*
             * Mean может существовать, но percentile
             * недостаточно надёжны для оценки.
             */
            out << group.mean << '\t'
                << "nan\t"
                << "nan";
        }

        out << '\n';
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
        const double covered =
            baseCountA[i] +
            baseCountC[i] +
            baseCountG[i] +
            baseCountT[i] +
            baseCountN[i];
        const double canonical =
            baseCountA[i] +
            baseCountC[i] +
            baseCountG[i] +
            baseCountT[i];

        double a = 0;
        double c = 0;
        double g = 0;
        double t = 0;
        double n = 0;

        // A/C/G/T are normalized among canonical calls. N remains relative
        // to every read covering the position (also written separately in
        // per_base_n_content_R*.tsv).
        if (canonical > 0)
        {
            a = baseCountA[i] * 100.0 / canonical;
            c = baseCountC[i] * 100.0 / canonical;
            g = baseCountG[i] * 100.0 / canonical;
            t = baseCountT[i] * 100.0 / canonical;
        }
        if (covered > 0)
        {
            n = baseCountN[i] * 100.0 / covered;
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
    const std::vector<double>& gcDistributionFastQC,
    const std::string& outDir,
    const std::string& readName)
{
    std::string path =
        outDir + "/per_sequence_gc_content_" + readName + ".tsv";

    std::ofstream out(path);

    if (!out)
        throw std::runtime_error("Cannot write to " + path);

    out << "gc_percent\traw_read_count\tfastqc_observed_count\n";

    for (size_t i = 0; i < 101; ++i)
    {
        const uint64_t rawCount =
            i < gcDistribution.size() ? gcDistribution[i] : 0;

        const double fastqcCount =
            i < gcDistributionFastQC.size()
                ? gcDistributionFastQC[i]
                : 0.0;

        out << i
            << "\t"
            << rawCount
            << "\t"
            << fastqcCount
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
            out << "sequence\tcount\tpercentage\n";
            out << std::fixed << std::setprecision(10);
            for (const auto& sequence : stats.overrepresentedSequences) {
                out << sequence.sequence << '\t'
                    << sequence.count << '\t'
                    << sequence.percent << '\n';
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
        outDir,
        stats.perBaseQualityGroups,
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
        stats.gcDistributionFastQC,
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

void createDirectoryOrThrow(const fs::path& path) {
    std::error_code error;
    fs::create_directories(path, error);
    if (error) {
        throw std::runtime_error(
            "Cannot create output directory '" + path.string() + "': "
            + error.message());
    }
}

std::vector<DuplicationEntry> collectMergedDuplicationEntries(
    const std::vector<QualityAnalyzer>& analyzers)
{
    std::vector<std::vector<DuplicationEntry>> entries;
    entries.reserve(analyzers.size());
    for (const auto& analyzer : analyzers) {
        entries.push_back(analyzer.getDuplicationEntries());
    }
    for (auto& localEntries : entries) {
        std::sort(
            localEntries.begin(),
            localEntries.end(),
            [](const DuplicationEntry& left, const DuplicationEntry& right) {
                return left.key.words < right.key.words;
            });
    }
    return mergeDuplicationEntriesTree(std::move(entries));
}

void mergeAnalyzers(QualityAnalyzer& destination,
                    const std::vector<QualityAnalyzer>& sources) {
    for (const auto& source : sources) {
        destination.merge(source);
    }
}

AnalysisResult processOneFileWithTrimming(
    const std::string& path,
    const std::string& readName,
    const std::string& outDir,
    const std::string& sampleId,
    bool skipAdapters,
    const TrimConfig& trimConfig)
{
    AnalysisResult result;
    const fs::path root(outDir);
    const fs::path beforeDir = root / "qc" / "before";
    const fs::path afterDir = root / "qc" / "after";
    const fs::path trimmedDir = root / "trimmed";
    const fs::path trimmedPath = trimmedDir / (readName + ".trimmed.fastq.gz");
    createDirectoryOrThrow(beforeDir);
    createDirectoryOrThrow(afterDir);
    createDirectoryOrThrow(trimmedDir);

    beginDuplicationArtifacts(beforeDir.string(), readName, path);
    beginDuplicationArtifacts(afterDir.string(), readName, trimmedPath.string());

    FastqReader reader(path);
    FastqWriter writer(trimmedPath);
    constexpr std::size_t BATCH_SIZE = 100000;
    const int threadCount = omp_get_max_threads();
    std::cout << "OpenMP threads: " << threadCount << "\n";

    std::vector<QualityAnalyzer> beforeAnalyzers;
    std::vector<QualityAnalyzer> afterAnalyzers;
    std::vector<Trimmer> trimmers;
    beforeAnalyzers.reserve(threadCount);
    afterAnalyzers.reserve(threadCount);
    trimmers.reserve(threadCount);
    for (int i = 0; i < threadCount; ++i) {
        beforeAnalyzers.emplace_back();
        afterAnalyzers.emplace_back();
        trimmers.emplace_back(trimConfig);
    }

    std::vector<FastqRecord> batch;
    std::vector<TrimResult> trimResults;
    std::size_t count = 0;
    while (reader.readBatch(batch, BATCH_SIZE)) {
        processTrimmedBatchParallel(
            beforeAnalyzers,
            afterAnalyzers,
            trimmers,
            batch,
            trimResults,
            skipAdapters);
        for (std::size_t i = 0; i < batch.size(); ++i) {
            if (trimResults[i].passed) {
                writer.write(batch[i]);
            }
        }
        count += batch.size();
        if (count % 1000000 == 0) {
            std::cout << "Processed " << count << " reads...\n";
        }
    }
    writer.close();

    QualityAnalyzer beforeAnalyzer;
    QualityAnalyzer afterAnalyzer;
    mergeAnalyzers(beforeAnalyzer, beforeAnalyzers);
    mergeAnalyzers(afterAnalyzer, afterAnalyzers);
    for (const auto& trimmer : trimmers) {
        result.trimmingStats.merge(trimmer.getStats());
    }

    const auto beforeEntries = collectMergedDuplicationEntries(beforeAnalyzers);
    const auto afterEntries = collectMergedDuplicationEntries(afterAnalyzers);
    const QualityStats beforeStats = beforeAnalyzer.getStats();
    const QualityStats afterStats = afterAnalyzer.getStats();
    const DuplicationStats beforeDuplication =
        beforeAnalyzer.getDuplicationStats(beforeEntries);
    const DuplicationStats afterDuplication =
        afterAnalyzer.getDuplicationStats(afterEntries);

    printConsoleSummary(beforeStats, readName + " BEFORE", skipAdapters);
    printConsoleSummary(afterStats, readName + " AFTER", skipAdapters);
    writeAnalysisReports(
        beforeStats,
        beforeDuplication,
        beforeAnalyzer,
        path,
        beforeDir.string(),
        sampleId,
        readName,
        skipAdapters);
    writeAnalysisReports(
        afterStats,
        afterDuplication,
        afterAnalyzer,
        trimmedPath.string(),
        afterDir.string(),
        sampleId,
        readName,
        skipAdapters);

    result.r1Stats = beforeStats;
    result.r1AfterStats = afterStats;
    return result;
}

// ---------------------------------------------------------------------------
// Обработка одного файла (R1 или R2)
// ---------------------------------------------------------------------------
AnalysisResult processOneFile(const std::string& path,
                            const std::string& readName,
                            const std::string& outDir,
                            const std::string& sampleId,
                            bool skipAdapters,
                            const TrimConfig* trimConfig = nullptr) {
    if (trimConfig != nullptr) {
        return processOneFileWithTrimming(
            path, readName, outDir, sampleId, skipAdapters, *trimConfig);
    }
    AnalysisResult result;
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

    std::vector<Trimmer> localTrimmers;
    if (trimConfig != nullptr)
    {
        localTrimmers.reserve(threadCount);
        for (int i = 0; i < threadCount; ++i)
        {
            localTrimmers.emplace_back(*trimConfig);
        }
    }

    std::vector<FastqRecord> batch;

    size_t count = 0;

    while (reader.readBatch(batch, BATCH_SIZE))
    {
        processBatchParallel(
            localAnalyzers,
            batch,
            skipAdapters,
            trimConfig != nullptr ? &localTrimmers : nullptr);

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
    for (const auto& localTrimmer : localTrimmers)
    {
        result.trimmingStats.merge(localTrimmer.getStats());
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


AnalysisResult processPairedFilesWithTrimming(
    const std::string& r1Path,
    const std::string& r2Path,
    const std::string& outDir,
    const std::string& sampleId,
    bool skipAdapters,
    const TrimConfig& trimConfig)
{
    AnalysisResult result;
    const fs::path root(outDir);
    const fs::path beforeDir = root / "qc" / "before";
    const fs::path afterDir = root / "qc" / "after";
    const fs::path trimmedDir = root / "trimmed";
    const fs::path trimmedR1Path = trimmedDir / "R1.trimmed.fastq.gz";
    const fs::path trimmedR2Path = trimmedDir / "R2.trimmed.fastq.gz";
    const fs::path mergedPath = trimmedDir / "merged.fastq.gz";
    createDirectoryOrThrow(beforeDir);
    createDirectoryOrThrow(afterDir);
    createDirectoryOrThrow(trimmedDir);

    beginDuplicationArtifacts(beforeDir.string(), "R1", r1Path);
    beginDuplicationArtifacts(beforeDir.string(), "R2", r2Path);
    beginDuplicationArtifacts(afterDir.string(), "R1", trimmedR1Path.string());
    beginDuplicationArtifacts(afterDir.string(), "R2", trimmedR2Path.string());

    FastqReader readerR1(r1Path);
    FastqReader readerR2(r2Path);
    FastqWriter writerR1(trimmedR1Path);
    FastqWriter writerR2(trimmedR2Path);
    std::unique_ptr<FastqWriter> mergedWriter;
    if (trimConfig.merge_reads) {
        mergedWriter = std::make_unique<FastqWriter>(mergedPath);
    }
    constexpr std::size_t BATCH_SIZE = 100000;
    const int threadCount = omp_get_max_threads();
    std::cout << "OpenMP threads: " << threadCount << "\n";

    std::vector<QualityAnalyzer> beforeAnalyzersR1;
    std::vector<QualityAnalyzer> beforeAnalyzersR2;
    std::vector<QualityAnalyzer> afterAnalyzersR1;
    std::vector<QualityAnalyzer> afterAnalyzersR2;
    std::vector<Trimmer> trimmers;
    beforeAnalyzersR1.reserve(threadCount);
    beforeAnalyzersR2.reserve(threadCount);
    afterAnalyzersR1.reserve(threadCount);
    afterAnalyzersR2.reserve(threadCount);
    trimmers.reserve(threadCount);
    for (int i = 0; i < threadCount; ++i) {
        beforeAnalyzersR1.emplace_back(ReadDirection::R1);
        beforeAnalyzersR2.emplace_back(ReadDirection::R2);
        afterAnalyzersR1.emplace_back(ReadDirection::R1);
        afterAnalyzersR2.emplace_back(ReadDirection::R2);
        trimmers.emplace_back(trimConfig);
    }

    std::vector<FastqRecord> batchR1;
    std::vector<FastqRecord> batchR2;
    std::vector<TrimResult> trimResultsR1;
    std::vector<TrimResult> trimResultsR2;
    std::vector<std::optional<FastqRecord>> mergedRecords;
    std::size_t count = 0;
    while (readPairedBatch(
        readerR1, readerR2, batchR1, batchR2, BATCH_SIZE))
    {
        processTrimmedPairedBatchParallel(
            beforeAnalyzersR1,
            beforeAnalyzersR2,
            afterAnalyzersR1,
            afterAnalyzersR2,
            trimmers,
            batchR1,
            batchR2,
            trimResultsR1,
            trimResultsR2,
            mergedRecords,
            skipAdapters);
        for (std::size_t i = 0; i < batchR1.size(); ++i) {
            if (trimResultsR1[i].passed && trimResultsR2[i].passed) {
                writerR1.write(batchR1[i]);
                writerR2.write(batchR2[i]);
                if (mergedWriter != nullptr && mergedRecords[i].has_value()) {
                    mergedWriter->write(*mergedRecords[i]);
                }
            }
        }
        count += batchR1.size();
        if (count % 1000000 == 0) {
            std::cout << "Processed " << count << " paired reads...\n";
        }
    }
    writerR1.close();
    writerR2.close();
    if (mergedWriter != nullptr) {
        mergedWriter->close();
    }

    QualityAnalyzer beforeAnalyzerR1(ReadDirection::R1);
    QualityAnalyzer beforeAnalyzerR2(ReadDirection::R2);
    QualityAnalyzer afterAnalyzerR1(ReadDirection::R1);
    QualityAnalyzer afterAnalyzerR2(ReadDirection::R2);
    mergeAnalyzers(beforeAnalyzerR1, beforeAnalyzersR1);
    mergeAnalyzers(beforeAnalyzerR2, beforeAnalyzersR2);
    mergeAnalyzers(afterAnalyzerR1, afterAnalyzersR1);
    mergeAnalyzers(afterAnalyzerR2, afterAnalyzersR2);
    for (const auto& trimmer : trimmers) {
        result.trimmingStats.merge(trimmer.getStats());
    }

    const auto beforeEntriesR1 = collectMergedDuplicationEntries(beforeAnalyzersR1);
    const auto beforeEntriesR2 = collectMergedDuplicationEntries(beforeAnalyzersR2);
    const auto afterEntriesR1 = collectMergedDuplicationEntries(afterAnalyzersR1);
    const auto afterEntriesR2 = collectMergedDuplicationEntries(afterAnalyzersR2);
    const QualityStats beforeStatsR1 = beforeAnalyzerR1.getStats();
    const QualityStats beforeStatsR2 = beforeAnalyzerR2.getStats();
    const QualityStats afterStatsR1 = afterAnalyzerR1.getStats();
    const QualityStats afterStatsR2 = afterAnalyzerR2.getStats();

    printConsoleSummary(beforeStatsR1, "R1 BEFORE", skipAdapters);
    printConsoleSummary(beforeStatsR2, "R2 BEFORE", skipAdapters);
    printConsoleSummary(afterStatsR1, "R1 AFTER", skipAdapters);
    printConsoleSummary(afterStatsR2, "R2 AFTER", skipAdapters);

    writeAnalysisReports(
        beforeStatsR1,
        beforeAnalyzerR1.getDuplicationStats(beforeEntriesR1),
        beforeAnalyzerR1,
        r1Path,
        beforeDir.string(),
        sampleId,
        "R1",
        skipAdapters);
    writeAnalysisReports(
        beforeStatsR2,
        beforeAnalyzerR2.getDuplicationStats(beforeEntriesR2),
        beforeAnalyzerR2,
        r2Path,
        beforeDir.string(),
        sampleId,
        "R2",
        skipAdapters);
    writeAnalysisReports(
        afterStatsR1,
        afterAnalyzerR1.getDuplicationStats(afterEntriesR1),
        afterAnalyzerR1,
        trimmedR1Path.string(),
        afterDir.string(),
        sampleId,
        "R1",
        skipAdapters);
    writeAnalysisReports(
        afterStatsR2,
        afterAnalyzerR2.getDuplicationStats(afterEntriesR2),
        afterAnalyzerR2,
        trimmedR2Path.string(),
        afterDir.string(),
        sampleId,
        "R2",
        skipAdapters);

    result.r1Stats = beforeStatsR1;
    result.r2Stats = beforeStatsR2;
    result.r1AfterStats = afterStatsR1;
    result.r2AfterStats = afterStatsR2;
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
                                const TrimConfig* trimConfig = nullptr)
{
    if (trimConfig != nullptr) {
        return processPairedFilesWithTrimming(
            r1Path, r2Path, outDir, sampleId, skipAdapters, *trimConfig);
    }
    AnalysisResult result;
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

    std::vector<Trimmer> localTrimmers;
    if (trimConfig != nullptr)
    {
        localTrimmers.reserve(threadCount);
        for (int i = 0; i < threadCount; ++i)
        {
            localTrimmers.emplace_back(*trimConfig);
        }
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
        processPairedBatchParallel(
            localAnalyzersR1,
            localAnalyzersR2,
            batchR1,
            batchR2,
            skipAdapters,
            trimConfig != nullptr ? &localTrimmers : nullptr);

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
    for (const auto& localTrimmer : localTrimmers)
    {
        result.trimmingStats.merge(localTrimmer.getStats());
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

std::vector<std::string> collectArtifacts(
    const fs::path& directory)
{
    std::vector<std::string> artifacts;

    if (!fs::exists(directory))
    {
        return artifacts;
    }

    for (const auto& entry :
         fs::recursive_directory_iterator(directory))
    {
        if (!entry.is_regular_file())
        {
            continue;
        }
        if (entry.path().filename() == "run_manifest.json")
        {
            continue;
        }
        artifacts.push_back(
            fs::relative(
                entry.path(),
                directory
            ).generic_string());
    }

    std::sort(
        artifacts.begin(),
        artifacts.end());

    return artifacts;
}

void writeRunManifest(
    const fs::path& outputDir,
    const RunManifest& manifest)
{
    const fs::path manifestPath =
        outputDir / "run_manifest.json";

    const std::vector<std::string> artifacts =
        collectArtifacts(outputDir);

    std::ofstream output(manifestPath);

    if (!output)
    {
        throw std::runtime_error(
            "Cannot write run manifest: " +
            manifestPath.string());
    }

    output << "{\n";
    output << "  \"schema_version\": 1,\n";

    output << "  \"run_id\": \""
           << jsonEscape(manifest.runId)
           << "\",\n";

    output << "  \"created_at\": \""
           << jsonEscape(manifest.createdAt)
           << "\",\n";

    output << "  \"neoqc_version\": \""
           << NEOQC_VERSION
           << "\",\n";

    output << "  \"run\": {\n";

    output << "    \"sample_id\": \""
           << jsonEscape(manifest.sampleId)
           << "\",\n";

    output << "    \"mode\": \""
           << (manifest.paired
               ? "paired-end"
               : "single-end")
           << "\",\n";

    output << "    \"plot\": "
           << (manifest.plot ? "true" : "false")
           << ",\n";

    output << "    \"skip_adapters\": "
           << (manifest.skipAdapters
               ? "true"
               : "false")
           << "\n";

    output << "  },\n";

    output << "  \"reads\": [";

    output << "\"R1\"";

    if (manifest.paired)
    {
        output << ", \"R2\"";
    }

    output << "],\n";

    output << "  \"inputs\": {\n";

    output << "    \"R1\": \""
           << jsonEscape(manifest.r1Path)
           << "\"";

    if (manifest.r2Path)
    {
        output << ",\n";

        output << "    \"R2\": \""
               << jsonEscape(*manifest.r2Path)
               << "\"";
    }

    output << "\n";
    output << "  },\n";

    output << "  \"artifacts\": [\n";

    for (std::size_t i = 0;
         i < artifacts.size();
         ++i)
    {
        output << "    \""
               << jsonEscape(artifacts[i])
               << "\"";

        if (i + 1 < artifacts.size())
        {
            output << ",";
        }

        output << "\n";
    }

    output << "  ]\n";
    output << "}\n";
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

//================================================================================
// Run ID generation and staging directory creation
//================================================================================

std::string generateRunId()
{
    std::random_device randomDevice;
    std::mt19937_64 generator(randomDevice());
    const uint64_t randomValue = generator();

    std::ostringstream output;
    output << currentUtcTimestamp()
           << "-"
           << std::hex
           << std::setw(16)
           << std::setfill('0')
           << randomValue;

    return output.str();
}

fs::path createStagingDirectory(
    const fs::path& outputDir,
    const std::string& runId)
{
    const fs::path parent = outputDir.parent_path();

    if (!parent.empty())
    {
        std::error_code ec;
        fs::create_directories(parent, ec);

        if (ec)
        {
            throw std::runtime_error(
                "Cannot create output parent directory '" +
                parent.string() + "': " +
                ec.message());
        }
    }

    const std::string outputName = outputDir.filename().string();

    if (outputName.empty())
    {
        throw std::runtime_error(
            "Output directory must have a valid directory name");
    }

    const fs::path stagingDir =
        parent / ("." + outputName + ".neoqc-tmp-" + runId);

    std::error_code ec;
    fs::create_directory(stagingDir, ec);

    if (ec)
    {
        throw std::runtime_error(
            "Cannot create staging directory '" +
            stagingDir.string() + "': " +
            ec.message());
    }

    return stagingDir;
}

bool pathIsInside(const fs::path& path, const fs::path& directory) {
    std::error_code pathError;
    const fs::path normalizedPath = fs::weakly_canonical(path, pathError);
    if (pathError) {
        throw std::runtime_error(
            "Cannot resolve input path '" + path.string() + "': "
            + pathError.message());
    }
    std::error_code directoryError;
    const fs::path normalizedDirectory =
        fs::weakly_canonical(directory, directoryError);
    if (directoryError) {
        throw std::runtime_error(
            "Cannot resolve output path '" + directory.string() + "': "
            + directoryError.message());
    }
    const fs::path relative = normalizedPath.lexically_relative(normalizedDirectory);
    return !relative.empty()
        && *relative.begin() != "..";
}

void ensureInputsOutsideOutput(const std::string& r1,
                               const std::string& r2,
                               const fs::path& outputDir) {
    if (pathIsInside(r1, outputDir)
        || (!r2.empty() && pathIsInside(r2, outputDir))) {
        throw std::runtime_error(
            "Trimming output directory must not contain an input FASTQ: "
            + outputDir.string());
    }
}

void removeStagingDirectory(const fs::path& stagingDir)
{
    std::error_code ec;
    fs::remove_all(stagingDir, ec);

    if (ec)
    {
        throw std::runtime_error(
            "Cannot remove staging directory '" +
            stagingDir.string() + "': " +
            ec.message());
    }
}

void publishRun(
    const fs::path& stagingDir,
    const fs::path& outputDir)
{
    const fs::path parent = outputDir.parent_path();

    if (!fs::exists(stagingDir))
    {
        throw std::runtime_error(
            "Staging directory does not exist: " +
            stagingDir.string());
    }

    const std::string runId = generateRunId();

    const fs::path backupDir =
        parent /
        ("." + outputDir.filename().string() +
         ".neoqc-backup-" + runId);

    std::error_code ec;

    bool hadPreviousOutput = fs::exists(outputDir);

    if (hadPreviousOutput)
    {
        fs::rename(outputDir, backupDir, ec);

        if (ec)
        {
            throw std::runtime_error(
                "Cannot move previous output directory '" +
                outputDir.string() +
                "' to backup '" +
                backupDir.string() +
                "': " +
                ec.message());
        }
    }

    fs::rename(stagingDir, outputDir, ec);

    if (ec)
    {
        if (hadPreviousOutput)
        {
            std::error_code restoreError;
            fs::rename(backupDir, outputDir, restoreError);

            if (restoreError)
            {
                throw std::runtime_error(
                    "Cannot publish new NeoQC result and cannot restore "
                    "previous output directory. "
                    "Original publish error: " +
                    ec.message() +
                    "; restore error: " +
                    restoreError.message());
            }
        }

        throw std::runtime_error(
            "Cannot publish NeoQC result '" +
            stagingDir.string() +
            "' to '" +
            outputDir.string() +
            "': " +
            ec.message());
    }

    if (hadPreviousOutput)
    {
        std::error_code cleanupError;
        fs::remove_all(backupDir, cleanupError);

        if (cleanupError)
        {
            throw std::runtime_error(
                "NeoQC result was published successfully, but previous "
                "output could not be removed: " +
                cleanupError.message());
        }
    }
}

AnalysisResult runSampleTransaction(
    const std::string& r1,
    const std::string& r2,
    const std::string& sampleId,
    const fs::path& outputDir,
    bool plot,
    bool skipAdapters,
    const TrimConfig* trimConfig = nullptr)
{
    const bool paired = !r2.empty();
    if (trimConfig != nullptr
        && (trimConfig->overlap_correction || trimConfig->merge_reads)
        && !paired) {
        throw std::runtime_error(
            "PE overlap correction and merging require paired-end input");
    }
    if (trimConfig != nullptr) {
        ensureInputsOutsideOutput(r1, r2, outputDir);
    }
    const std::string runId = generateRunId();
    fs::path stagingDir = createStagingDirectory(outputDir, runId);
    try
    {
        AnalysisResult analysis = paired
            ? processPairedFiles(r1, r2, stagingDir.string(), sampleId, skipAdapters, trimConfig)
            : processOneFile(r1, "R1", stagingDir.string(), sampleId, skipAdapters, trimConfig);

        if (trimConfig != nullptr && trimConfig->enabled)
        {
            writeTrimmingReport(
                stagingDir / "trimming_report.json",
                analysis.trimmingStats);
        }

        const RunManifest manifest{
            .runId = runId,
            .createdAt = currentUtcTimestamp(),
            .sampleId = sampleId,
            .paired = paired,
            .plot = plot,
            .skipAdapters = skipAdapters,
            .r1Path = r1,
            .r2Path = paired ? std::optional<std::string>(r2) : std::nullopt,
            .artifacts = {},
        };
        // Plot generation uses the evaluation engine, which is manifest-only.
        writeRunManifest(stagingDir, manifest);
        if (plot)
        {
            PlotOptions options;
            options.includeAdapters = !skipAdapters;
            if (trimConfig != nullptr) {
                PlotRunner::runAll(
                    (stagingDir / "qc" / "before").string(),
                    (stagingDir / "plots" / "before").string(),
                    options);
                PlotRunner::runAll(
                    (stagingDir / "qc" / "after").string(),
                    (stagingDir / "plots" / "after").string(),
                    options);
            } else {
                PlotRunner::runAll(
                    stagingDir.string(),
                    (stagingDir / "plots").string(),
                    options);
            }
        }
        writeRunManifest(stagingDir, manifest);
        publishRun(stagingDir, outputDir);
        return analysis;
    }
    catch (...)
    {
        if (fs::exists(stagingDir))
        {
            std::error_code ignored;
            fs::remove_all(stagingDir, ignored);
        }
        throw;
    }
}

void writeCaseSummary(
    const std::string& patientId,
    const fs::path& caseOutputDir,
    const std::vector<BatchSampleResult>& results,
    const std::vector<std::string>& warnings,
    const std::string& batchRunId,
    bool plot,
    bool skipAdapters)
{
    bool casePassed = true;

    for (const auto& result : results)
    {
        if (result.entry.patientId == patientId &&
            !result.passed)
        {
            casePassed = false;
        }
    }

    const fs::path path =
        caseOutputDir / "case_summary.json";

    writeAtomically(path, [&](std::ostream& output)
    {
        output << "{\n"
               << "  \"patient_id\": \""
               << jsonEscape(patientId)
               << "\",\n"
               << "  \"status\": \""
               << (casePassed ? "passed" : "failed")
               << "\",\n"
               << "  \"neoqc_version\": \""
               << NEOQC_VERSION
               << "\",\n"
               << "  \"run_id\": \""
               << jsonEscape(batchRunId)
               << "\",\n"
               << "  \"parameters\": {\"plot\": "
               << (plot ? "true" : "false")
               << ", \"skip_adapters\": "
               << (skipAdapters ? "true" : "false")
               << "},\n"
               << "  \"ruleset\": {"
                  "\"id\": \"fastqc-compatible-v1\", "
                  "\"version\": \"1.0.0\", "
                  "\"sha256\": \"" NEOQC_RULESET_SHA256 "\"},\n"
               << "  \"run_date\": \""
               << currentUtcTimestamp()
               << "\",\n"
               << "  \"warnings\": [";

        bool firstWarning = true;

        for (const auto& warning : warnings)
        {
            if (warning.find(
                    "Patient " + patientId + ":")
                == std::string::npos)
            {
                continue;
            }

            if (!firstWarning)
                output << ", ";

            output << "\""
                   << jsonEscape(warning)
                   << "\"";

            firstWarning = false;
        }

        output << "],\n"
               << "  \"samples\": [\n";

        bool firstSample = true;

        for (const auto& result : results)
        {
            if (result.entry.patientId != patientId)
                continue;

            if (!firstSample)
                output << ",\n";

            firstSample = false;

            const auto& entry = result.entry;

            output << "    {\n"
                   << "      \"sample_id\": \""
                   << jsonEscape(entry.sampleId)
                   << "\",\n"
                   << "      \"role\": \""
                   << jsonEscape(entry.sampleRole)
                   << "\",\n"
                   << "      \"material\": \""
                   << jsonEscape(entry.material)
                   << "\",\n"
                   << "      \"r1\": \""
                   << jsonEscape(entry.r1)
                   << "\",\n"
                   << "      \"r2\": \""
                   << jsonEscape(entry.r2)
                   << "\",\n"
                   << "      \"qc_status\": \""
                   << (result.passed
                           ? "passed"
                           : "failed")
                   << "\"";

            if (result.passed &&
                result.analysis)
            {
                output
                    << ",\n      \"r1_reads\": "
                    << result.analysis->r1Stats.totalReads;

                if (result.analysis->r2Stats)
                {
                    output
                        << ",\n      \"r2_reads\": "
                        << result.analysis->r2Stats->totalReads;
                }
            }
            else
            {
                output
                    << ",\n      \"qc_error\": \""
                    << jsonEscape(result.error)
                    << "\"";
            }

            output << "\n    }";
        }

        output << "\n  ]\n}\n";
    });
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

    // Safe integration point for the future trimming pipeline.  The object
    // only owns configuration and component structure until algorithms exist.
    const Trimmer trimmer(args.trimConfig);
    const TrimConfig* activeTrimConfig = args.trimConfig.enabled
        ? &args.trimConfig
        : nullptr;

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

            const std::string batchRunId = generateRunId();
            std::vector<BatchSampleResult> results;
            bool allPassed = true;
            for (const auto& entry : entries) {
                const fs::path sampleOutDir = fs::path(args.outDir)
                                              / entry.patientId / entry.sampleId;
                std::cout << "\nCase " << entry.patientId << " — checking " << entry.sampleId
                          << "\nR1: " << entry.r1 << "\n";
                if (!entry.r2.empty()) std::cout << "R2: " << entry.r2 << "\n";

                try {
                    AnalysisResult analysis = runSampleTransaction(
                        entry.r1, entry.r2, entry.sampleId, sampleOutDir,
                        args.plot, args.skipAdapters, activeTrimConfig);
                    std::cout << "Result: passed\n";
                    results.push_back({entry, true, std::move(analysis), ""});
                } catch (const std::exception& e) {
                    std::cerr << "Result: failed for " << entry.sampleId << ": " << e.what() << "\n";
                    allPassed = false;
                    results.push_back({entry, false, std::nullopt, e.what()});
                }
            }
            
            {
                std::vector<std::string> patientIds;
                for (const auto& entry : entries)
                {
                    if (std::find(patientIds.begin(), patientIds.end(), entry.patientId) == patientIds.end())
                        patientIds.push_back(entry.patientId);
                }
                for (const auto& patientId : patientIds)
                {
                    const fs::path caseOutDir = fs::path(args.outDir) / patientId;
                    writeCaseSummary(patientId, caseOutDir, results, warnings,
                                     batchRunId, args.plot, args.skipAdapters);
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


    const bool isPaired = !args.r2.empty();
    if (activeTrimConfig != nullptr) {
        try {
            ensureInputsOutsideOutput(args.r1, args.r2, fs::path(args.outDir));
        } catch (const std::exception& error) {
            std::cerr << "Error: " << error.what() << '\n';
            return 1;
        }
    }
    std::cout << "Sample ID : " << args.sampleId << "\n"
              << "Mode      : " << (isPaired ? "paired-end" : "single-end") << "\n"
              << "R1        : " << args.r1 << "\n";
    if (isPaired) std::cout << "R2        : " << args.r2 << "\n";
    std::cout << "Output    : " << args.outDir << "\n";
    if (args.skipAdapters) std::cout << "Adapters  : skipped\n";
    if (trimmer.getConfig().enabled) {
        std::cout << "Trimming  : enabled\n";
    }

    const std::string runId = generateRunId();

    fs::path stagingDir;
    RunManifest manifest;

    try
    {
        stagingDir = createStagingDirectory(
            fs::path(args.outDir),
            runId);

        std::cout << "Run ID    : " << runId << "\n";
        std::cout << "Staging   : " << stagingDir << "\n";

        AnalysisResult analysis;
        if (isPaired)
        {
            analysis = processPairedFiles(
                args.r1,
                args.r2,
                stagingDir.string(),
                args.sampleId,
                args.skipAdapters,
                activeTrimConfig);
        }
        else
        {
            analysis = processOneFile(
                args.r1,
                "R1",
                stagingDir.string(),
                args.sampleId,
                args.skipAdapters,
                activeTrimConfig);
        }

        if (activeTrimConfig != nullptr)
        {
            writeTrimmingReport(
                stagingDir / "trimming_report.json",
                analysis.trimmingStats);
        }

        manifest = RunManifest{
            .runId = runId,
            .createdAt = currentUtcTimestamp(),
            .sampleId = args.sampleId,
            .paired = isPaired,
            .plot = args.plot,
            .skipAdapters = args.skipAdapters,
            .r1Path = args.r1,
            .r2Path = args.r2.empty()
                ? std::nullopt
                : std::optional<std::string>(args.r2),
            .artifacts = {},
        };

        writeRunManifest(
            stagingDir,
            manifest);
    }
    catch (const std::exception& e)
    {
        std::cerr << e.what() << '\n';

        if (!stagingDir.empty())
        {
            try
            {
                removeStagingDirectory(stagingDir);
            }
            catch (const std::exception& cleanupError)
            {
                std::cerr
                    << "Warning: failed to remove staging directory: "
                    << cleanupError.what()
                    << '\n';
            }
        }

        return 1;
    }

    // Построение графиков (опционально, через PlotRunner)
    if (args.plot)
    {
        try
        {
            PlotOptions plotOptions;
            plotOptions.includeAdapters = !args.skipAdapters;
            if (activeTrimConfig != nullptr) {
                PlotRunner::runAll(
                    (stagingDir / "qc" / "before").string(),
                    (stagingDir / "plots" / "before").string(),
                    plotOptions);
                PlotRunner::runAll(
                    (stagingDir / "qc" / "after").string(),
                    (stagingDir / "plots" / "after").string(),
                    plotOptions);
            } else {
                PlotRunner::runAll(
                    stagingDir.string(),
                    (stagingDir / "plots").string(),
                    plotOptions);
            }
        }
        catch (const std::exception& e)
        {
            std::cerr << e.what() << '\n';

            try
            {
                removeStagingDirectory(stagingDir);
            }
            catch (const std::exception& cleanupError)
            {
                std::cerr
                    << "Warning: failed to remove staging directory: "
                    << cleanupError.what()
                    << '\n';
            }

            return 1;
        }
    }

    writeRunManifest(
        stagingDir,
        manifest);

    try
    {
        publishRun(
            stagingDir,
            fs::path(args.outDir));
    }
    catch (const std::exception& e)
    {
        std::cerr << e.what() << '\n';

        if (fs::exists(stagingDir))
        {
            try
            {
                removeStagingDirectory(stagingDir);
            }
            catch (const std::exception& cleanupError)
            {
                std::cerr
                    << "Warning: failed to remove staging directory: "
                    << cleanupError.what()
                    << '\n';
            }
        }

        return 1;
    }


    const auto end = Clock::now();

    const auto elapsed = std::chrono::duration<double>(end - start);

    std::cout << "Done.\n"
              << "Total processing time: "
              << elapsed.count()
              << " s\n";

    return 0;
}
