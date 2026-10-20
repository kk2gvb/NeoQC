#include "fastq_reader.h"
#include "trimming/trim_report.h"
#include "trimming/trimmer.h"

#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>

namespace {

FastqRecord record(const std::string& sequence,
                   const std::string& quality = "") {
    FastqRecord value;
    value.header = "@READ";
    value.separator = "+";
    value.sequence = sequence;
    value.quality = quality.empty() ? std::string(sequence.size(), '?') : quality;
    return value;
}

bool check(bool condition, const char* message) {
    if (!condition) std::cerr << message << '\n';
    return condition;
}

std::string readFile(const std::filesystem::path& path) {
    std::ifstream input(path);
    return {std::istreambuf_iterator<char>(input), std::istreambuf_iterator<char>()};
}

}  // namespace

int main() {
    TrimConfig config;
    config.enabled = true;
    config.trim_front = 2;
    config.min_length = 3;
    Trimmer trimmer(config);

    FastqRecord passed = record("ABCDE");
    const TrimResult passedResult = trimmer.trim(passed);
    if (!check(passed.sequence == "CDE" && passed.quality == "???",
               "filtering changed a passing read incorrectly")
        || !check(passedResult.passed && !passedResult.discard_reason,
                  "length equal to minimum must pass")) return 1;

    FastqRecord discarded = record("ABCD");
    const TrimResult discardedResult = trimmer.trim(discarded);
    if (!check(discarded.sequence == "CD" && discarded.quality == "??",
               "filtering must not change a discarded read")
        || !check(!discardedResult.passed
                  && discardedResult.discard_reason == "TOO_SHORT",
                  "short read discard result")) return 1;

    FastqRecord shorterThanMinimum = record("AB");
    const TrimResult tooShortResult = trimmer.trim(shorterThanMinimum);
    if (!check(!tooShortResult.passed
               && tooShortResult.discard_reason == "TOO_SHORT"
               && shorterThanMinimum.sequence.empty()
               && shorterThanMinimum.quality.empty(),
               "shorter-than-minimum read must be discarded after trimming")) return 1;

    TrimConfig filterOnlyConfig;
    filterOnlyConfig.enabled = true;
    filterOnlyConfig.min_length = 5;
    Trimmer filterOnlyTrimmer(filterOnlyConfig);
    FastqRecord untrimmedDiscard = record("ABCD");
    const std::string untrimmedSequence = untrimmedDiscard.sequence;
    const std::string untrimmedQuality = untrimmedDiscard.quality;
    const TrimResult untrimmedDiscardResult = filterOnlyTrimmer.trim(untrimmedDiscard);
    if (!check(!untrimmedDiscardResult.passed
               && untrimmedDiscardResult.discard_reason == "TOO_SHORT"
               && untrimmedDiscard.sequence == untrimmedSequence
               && untrimmedDiscard.quality == untrimmedQuality
               && untrimmedDiscard.sequence.size() == untrimmedDiscard.quality.size(),
               "filtering must not mutate sequence or quality")) return 1;

    TrimConfig zeroConfig;
    zeroConfig.enabled = true;
    Trimmer zeroTrimmer(zeroConfig);
    FastqRecord empty = record("");
    if (!check(zeroTrimmer.trim(empty).passed, "min_length zero must retain empty reads")) {
        return 1;
    }

    TrimConfig emptyDiscardConfig;
    emptyDiscardConfig.enabled = true;
    emptyDiscardConfig.min_length = 1;
    Trimmer emptyDiscardTrimmer(emptyDiscardConfig);
    FastqRecord emptyDiscard = record("");
    if (!check(!emptyDiscardTrimmer.trim(emptyDiscard).passed,
               "empty read must be discarded when min_length is positive")) return 1;

    const TrimStats& stats = trimmer.getStats();
    if (!check(stats.total_reads == 3 && stats.passed_reads == 1
               && stats.discarded_reads == 2 && stats.too_short_reads == 2
               && stats.bases_before == 11 && stats.bases_after == 5
               && stats.bases_trimmed == 6
               && stats.total_reads == stats.passed_reads + stats.discarded_reads
               && stats.bases_before == stats.bases_after + stats.bases_trimmed,
               "filtering statistics")) return 1;

    TrimConfig qualityConfig;
    qualityConfig.enabled = true;
    qualityConfig.cut_tail = true;
    qualityConfig.quality_threshold = 20;
    Trimmer qualityTrimmer(qualityConfig);
    FastqRecord qualityRead = record("ABCD", std::string(2, '?') + "!!");
    qualityTrimmer.trim(qualityRead);
    if (!check(qualityRead.sequence == "AB"
               && qualityTrimmer.getStats().quality_trimmed_reads == 1,
               "quality trimming counter")) return 1;

    TrimConfig adapterConfig;
    adapterConfig.enabled = true;
    adapterConfig.adapter_trimming = true;
    adapterConfig.adapter_sequence = "TTTT";
    adapterConfig.min_adapter_match = 4;
    adapterConfig.min_overlap = 10;
    Trimmer adapterTrimmer(adapterConfig);
    FastqRecord adapterR1 = record("AACCTTTT");
    FastqRecord adapterR2 = record("GGTTTT");
    adapterTrimmer.trimPair(adapterR1, adapterR2);
    const TrimStats& adapterStats = adapterTrimmer.getStats();
    if (!check(adapterStats.total_reads == 2
               && adapterStats.adapter_trimmed_reads == 2
               && adapterStats.adapter_positions.at(2) == 1
               && adapterStats.adapter_positions.at(4) == 1,
               "paired adapter positions must count mates separately")) return 1;

    TrimStats positionsA;
    positionsA.adapter_positions[0] = 1;
    positionsA.adapter_positions[5] = 2;
    TrimStats positionsB;
    positionsB.adapter_positions[0] = 3;
    positionsB.adapter_positions[9] = 1;
    positionsA.merge(positionsB);
    if (!check(positionsA.adapter_positions[0] == 4
               && positionsA.adapter_positions[5] == 2
               && positionsA.adapter_positions[9] == 1,
               "adapter position merge")) return 1;

    const std::filesystem::path reportPath =
        std::filesystem::temp_directory_path() / "neoqc_trim_filter_report_test.json";
    writeTrimmingReport(reportPath, adapterStats);
    const std::string report = readFile(reportPath);
    std::filesystem::remove(reportPath);
    const std::size_t position2 = report.find("\"2\": 1");
    const std::size_t position4 = report.find("\"4\": 1");
    if (!check(report.find("\"total\": 2") != std::string::npos
               && report.find("\"adapter_trimmed_reads\": 2") != std::string::npos
               && position2 != std::string::npos && position4 != std::string::npos
               && position2 < position4,
               "trimming report contents and deterministic adapter order")) return 1;

    const std::filesystem::path zeroReportPath =
        std::filesystem::temp_directory_path() / "neoqc_trim_zero_report_test.json";
    writeTrimmingReport(zeroReportPath, TrimStats{});
    const std::string zeroReport = readFile(zeroReportPath);
    std::filesystem::remove(zeroReportPath);
    if (!check(zeroReport.find("\"total\": 0") != std::string::npos
               && zeroReport.find("\"adapter_positions\": {}") != std::string::npos,
               "zero-read report")) return 1;

    TrimStats invalidStats;
    invalidStats.total_reads = 1;
    bool rejectedInvalidStats = false;
    try {
        writeTrimmingReport(reportPath, invalidStats);
    } catch (const std::logic_error&) {
        rejectedInvalidStats = true;
    }
    return check(rejectedInvalidStats, "report must reject invalid statistics") ? 0 : 1;
}
