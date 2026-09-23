#include "fastq_reader.h"
#include "trimming/quality_trimmer.h"
#include "trimming/trimmer.h"

#include <initializer_list>
#include <iostream>
#include <string>

namespace {

FastqRecord makeRecord(const std::string& sequence,
                       std::initializer_list<int> phredScores) {
    FastqRecord record;
    record.header = "@READ_001/1";
    record.separator = "+";
    record.recordNumber = 42;
    record.sequence = sequence;
    for (const int score : phredScores) {
        record.quality.push_back(static_cast<char>(score + 33));
    }
    return record;
}

bool require(bool condition, const std::string& message) {
    if (!condition) {
        std::cerr << message << '\n';
    }
    return condition;
}

bool verifyRecord(const FastqRecord& record,
                  const std::string& expectedSequence,
                  const std::string& expectedQuality,
                  const std::string& message) {
    return require(record.sequence == expectedSequence, message + ": sequence")
        && require(record.quality == expectedQuality, message + ": quality")
        && require(record.sequence.size() == record.quality.size(),
                   message + ": sequence/quality lengths diverged")
        && require(record.header == "@READ_001/1" && record.separator == "+"
                   && record.recordNumber == 42,
                   message + ": FASTQ metadata changed");
}

bool verifyResult(const TrimResult& result,
                  std::size_t originalLength,
                  std::size_t finalLength,
                  std::size_t trimmedFront,
                  std::size_t trimmedTail,
                  bool qualityTrimmed,
                  const std::string& message) {
    return require(result.original_length == originalLength,
                   message + ": original length")
        && require(result.final_length == finalLength, message + ": final length")
        && require(result.trimmed_front == trimmedFront, message + ": trimmed front")
        && require(result.trimmed_tail == trimmedTail, message + ": trimmed tail")
        && require(result.quality_trimmed == qualityTrimmed,
                   message + ": quality-trim marker");
}

bool testFixedFront() {
    TrimConfig config;
    config.enabled = true;

    FastqRecord record = makeRecord("ACGT", {30, 31, 32, 33});
    QualityTrimmer trimmer(config);
    TrimResult result = trimmer.trim(record);
    if (!verifyRecord(record, "ACGT", "?@AB", "fixed front zero")
        || !verifyResult(result, 4, 4, 0, 0, false, "fixed front zero")) return false;

    config.trim_front = 1;
    record = makeRecord("ACGT", {30, 31, 32, 33});
    result = QualityTrimmer(config).trim(record);
    if (!verifyRecord(record, "CGT", "@AB", "fixed front one")
        || !verifyResult(result, 4, 3, 1, 0, false, "fixed front one")) return false;

    config.trim_front = 3;
    record = makeRecord("ACGT", {30, 31, 32, 33});
    result = QualityTrimmer(config).trim(record);
    if (!verifyRecord(record, "T", "B", "fixed front N")
        || !verifyResult(result, 4, 1, 3, 0, false, "fixed front N")) return false;

    config.trim_front = 4;
    record = makeRecord("ACGT", {30, 31, 32, 33});
    result = QualityTrimmer(config).trim(record);
    if (!verifyRecord(record, "", "", "fixed front equal length")
        || !verifyResult(result, 4, 0, 4, 0, false,
                         "fixed front equal length")) return false;

    config.trim_front = 10;
    record = makeRecord("ACGT", {30, 31, 32, 33});
    result = QualityTrimmer(config).trim(record);
    return verifyRecord(record, "", "", "fixed front exceeds length")
        && verifyResult(result, 4, 0, 4, 0, false,
                        "fixed front exceeds length");
}

bool testFixedTail() {
    TrimConfig config;
    config.enabled = true;

    FastqRecord record = makeRecord("ACGT", {30, 31, 32, 33});
    TrimResult result = QualityTrimmer(config).trim(record);
    if (!verifyRecord(record, "ACGT", "?@AB", "fixed tail zero")
        || !verifyResult(result, 4, 4, 0, 0, false, "fixed tail zero")) return false;

    config.trim_tail = 1;
    record = makeRecord("ACGT", {30, 31, 32, 33});
    result = QualityTrimmer(config).trim(record);
    if (!verifyRecord(record, "ACG", "?@A", "fixed tail one")
        || !verifyResult(result, 4, 3, 0, 1, false, "fixed tail one")) return false;

    config.trim_tail = 3;
    record = makeRecord("ACGT", {30, 31, 32, 33});
    result = QualityTrimmer(config).trim(record);
    if (!verifyRecord(record, "A", "?", "fixed tail N")
        || !verifyResult(result, 4, 1, 0, 3, false, "fixed tail N")) return false;

    config.trim_tail = 4;
    record = makeRecord("ACGT", {30, 31, 32, 33});
    result = QualityTrimmer(config).trim(record);
    if (!verifyRecord(record, "", "", "fixed tail equal length")
        || !verifyResult(result, 4, 0, 0, 4, false,
                         "fixed tail equal length")) return false;

    config.trim_tail = 10;
    record = makeRecord("ACGT", {30, 31, 32, 33});
    result = QualityTrimmer(config).trim(record);
    return verifyRecord(record, "", "", "fixed tail exceeds length")
        && verifyResult(result, 4, 0, 0, 4, false,
                        "fixed tail exceeds length");
}

bool testFixedFrontAndTail() {
    TrimConfig config;
    config.enabled = true;
    config.trim_front = 3;
    config.trim_tail = 2;

    FastqRecord record = makeRecord("ABCDEFGHIJ", {30, 30, 30, 30, 30,
                                                      30, 30, 30, 30, 30});
    TrimResult result = QualityTrimmer(config).trim(record);
    if (!verifyRecord(record, "DEFGH", "?????", "fixed front and tail")
        || !verifyResult(result, 10, 5, 3, 2, false,
                         "fixed front and tail")) return false;

    config.trim_front = 8;
    config.trim_tail = 5;
    record = makeRecord("ABCDEFGHIJ", {30, 30, 30, 30, 30,
                                        30, 30, 30, 30, 30});
    result = QualityTrimmer(config).trim(record);
    return verifyRecord(record, "", "", "fixed front and tail exceed length")
        && verifyResult(result, 10, 0, 8, 2, false,
                        "fixed front and tail exceed length");
}

bool testQualityFront() {
    TrimConfig config;
    config.enabled = true;
    config.cut_front = true;
    config.quality_threshold = 20;

    FastqRecord record = makeRecord("ABC", {30, 30, 30});
    TrimResult result = QualityTrimmer(config).trim(record);
    if (!verifyRecord(record, "ABC", "???", "quality front all good")
        || !verifyResult(result, 3, 3, 0, 0, false,
                         "quality front all good")) return false;

    record = makeRecord("ABC", {5, 30, 30});
    result = QualityTrimmer(config).trim(record);
    if (!verifyRecord(record, "BC", "??", "quality front one low")
        || !verifyResult(result, 3, 2, 1, 0, true,
                         "quality front one low")) return false;

    record = makeRecord("ABCDE", {5, 8, 12, 30, 30});
    result = QualityTrimmer(config).trim(record);
    if (!verifyRecord(record, "DE", "??", "quality front contiguous lows")
        || !verifyResult(result, 5, 2, 3, 0, true,
                         "quality front contiguous lows")) return false;

    record = makeRecord("ABC", {20, 5, 5});
    result = QualityTrimmer(config).trim(record);
    if (!verifyRecord(record, "ABC", "5&&", "quality front threshold boundary")
        || !verifyResult(result, 3, 3, 0, 0, false,
                         "quality front threshold boundary")) return false;

    record = makeRecord("ABCDE", {30, 5, 5, 30, 30});
    result = QualityTrimmer(config).trim(record);
    if (!verifyRecord(record, "ABCDE", "?&&??", "quality front middle lows")
        || !verifyResult(result, 5, 5, 0, 0, false,
                         "quality front middle lows")) return false;

    record = makeRecord("ABCDE", {30, 30, 30, 5, 5});
    result = QualityTrimmer(config).trim(record);
    if (!verifyRecord(record, "ABCDE", "???&&", "quality front tail lows")
        || !verifyResult(result, 5, 5, 0, 0, false,
                         "quality front tail lows")) return false;

    record = makeRecord("ABCDE", {5, 5, 30, 5, 5});
    result = QualityTrimmer(config).trim(record);
    if (!verifyRecord(record, "CDE", "?&&", "quality front stops at good")
        || !verifyResult(result, 5, 3, 2, 0, true,
                         "quality front stops at good")) return false;

    record = makeRecord("ABC", {5, 5, 5});
    result = QualityTrimmer(config).trim(record);
    if (!verifyRecord(record, "", "", "quality front all low")
        || !verifyResult(result, 3, 0, 3, 0, true,
                         "quality front all low")) return false;

    record = makeRecord("A", {5});
    result = QualityTrimmer(config).trim(record);
    if (!verifyRecord(record, "", "", "quality front one base")
        || !verifyResult(result, 1, 0, 1, 0, true,
                         "quality front one base")) return false;

    record = makeRecord("A", {30});
    result = QualityTrimmer(config).trim(record);
    if (!verifyRecord(record, "A", "?", "quality front one good base")
        || !verifyResult(result, 1, 1, 0, 0, false,
                         "quality front one good base")) return false;

    config.trim_front = 1;
    record = makeRecord("ABCD", {30, 5, 30, 30});
    result = QualityTrimmer(config).trim(record);
    if (!verifyRecord(record, "CD", "??", "fixed and quality front")
        || !verifyResult(result, 4, 2, 2, 0, true,
                         "fixed and quality front")) return false;

    config.trim_front = 0;
    record = makeRecord("", {});
    result = QualityTrimmer(config).trim(record);
    return verifyRecord(record, "", "", "quality front empty")
        && verifyResult(result, 0, 0, 0, 0, false, "quality front empty");
}

bool testQualityTail() {
    TrimConfig config;
    config.enabled = true;
    config.cut_tail = true;
    config.quality_threshold = 20;

    FastqRecord record = makeRecord("ABC", {30, 30, 30});
    TrimResult result = QualityTrimmer(config).trim(record);
    if (!verifyRecord(record, "ABC", "???", "quality tail all good")
        || !verifyResult(result, 3, 3, 0, 0, false,
                         "quality tail all good")) return false;

    record = makeRecord("ABC", {30, 30, 5});
    result = QualityTrimmer(config).trim(record);
    if (!verifyRecord(record, "AB", "??", "quality tail one low")
        || !verifyResult(result, 3, 2, 0, 1, true,
                         "quality tail one low")) return false;

    record = makeRecord("ABCDE", {30, 30, 12, 8, 5});
    result = QualityTrimmer(config).trim(record);
    if (!verifyRecord(record, "AB", "??", "quality tail contiguous lows")
        || !verifyResult(result, 5, 2, 0, 3, true,
                         "quality tail contiguous lows")) return false;

    record = makeRecord("ABC", {5, 5, 20});
    result = QualityTrimmer(config).trim(record);
    if (!verifyRecord(record, "ABC", "&&5", "quality tail threshold boundary")
        || !verifyResult(result, 3, 3, 0, 0, false,
                         "quality tail threshold boundary")) return false;

    record = makeRecord("ABCDE", {5, 5, 30, 30, 30});
    result = QualityTrimmer(config).trim(record);
    if (!verifyRecord(record, "ABCDE", "&&???", "quality tail front lows")
        || !verifyResult(result, 5, 5, 0, 0, false,
                         "quality tail front lows")) return false;

    record = makeRecord("ABCDE", {30, 30, 5, 5, 30});
    result = QualityTrimmer(config).trim(record);
    if (!verifyRecord(record, "ABCDE", "??&&?", "quality tail middle lows")
        || !verifyResult(result, 5, 5, 0, 0, false,
                         "quality tail middle lows")) return false;

    record = makeRecord("ABCDE", {5, 5, 30, 5, 5});
    result = QualityTrimmer(config).trim(record);
    if (!verifyRecord(record, "ABC", "&&?", "quality tail stops at good")
        || !verifyResult(result, 5, 3, 0, 2, true,
                         "quality tail stops at good")) return false;

    record = makeRecord("ABC", {5, 5, 5});
    result = QualityTrimmer(config).trim(record);
    if (!verifyRecord(record, "", "", "quality tail all low")
        || !verifyResult(result, 3, 0, 0, 3, true,
                         "quality tail all low")) return false;

    record = makeRecord("A", {5});
    result = QualityTrimmer(config).trim(record);
    if (!verifyRecord(record, "", "", "quality tail one base")
        || !verifyResult(result, 1, 0, 0, 1, true,
                         "quality tail one base")) return false;

    record = makeRecord("A", {30});
    result = QualityTrimmer(config).trim(record);
    if (!verifyRecord(record, "A", "?", "quality tail one good base")
        || !verifyResult(result, 1, 1, 0, 0, false,
                         "quality tail one good base")) return false;

    config.trim_tail = 1;
    record = makeRecord("ABCD", {30, 30, 5, 30});
    result = QualityTrimmer(config).trim(record);
    if (!verifyRecord(record, "AB", "??", "fixed and quality tail")
        || !verifyResult(result, 4, 2, 0, 2, true,
                         "fixed and quality tail")) return false;

    config.trim_tail = 0;
    record = makeRecord("", {});
    result = QualityTrimmer(config).trim(record);
    return verifyRecord(record, "", "", "quality tail empty")
        && verifyResult(result, 0, 0, 0, 0, false, "quality tail empty");
}

bool testTrimmerStatsAndPairedReads() {
    TrimConfig config;
    config.enabled = true;
    config.cut_front = true;
    config.cut_tail = true;
    config.quality_threshold = 20;

    Trimmer trimmer(config);
    FastqRecord r1 = makeRecord("ABCD", {5, 30, 30, 5});
    FastqRecord r2 = makeRecord("WXYZ", {5, 30, 30, 5});
    const TrimResult r1Result = trimmer.trim(r1);
    const TrimResult r2Result = trimmer.trim(r2);

    const TrimStats& stats = trimmer.getStats();
    return verifyRecord(r1, "BC", "??", "paired R1")
        && verifyRecord(r2, "XY", "??", "paired R2")
        && verifyResult(r1Result, 4, 2, 1, 1, true, "paired R1")
        && verifyResult(r2Result, 4, 2, 1, 1, true, "paired R2")
        && require(stats.total_reads == 2 && stats.passed_reads == 2
                   && stats.discarded_reads == 0 && stats.bases_before == 8
                   && stats.bases_after == 4 && stats.bases_trimmed == 4
                   && stats.quality_trimmed_reads == 2
                   && stats.too_short_reads == 0,
                   "Trimmer statistics for paired records");
}

bool testDisabledConfiguration() {
    TrimConfig config;
    config.trim_front = 2;
    config.trim_tail = 2;
    config.cut_front = true;
    config.cut_tail = true;
    config.quality_threshold = 20;

    Trimmer trimmer(config);
    FastqRecord record = makeRecord("ABCD", {5, 5, 5, 5});
    const TrimResult result = trimmer.trim(record);
    const TrimStats& stats = trimmer.getStats();

    return verifyRecord(record, "ABCD", "&&&&", "disabled configuration")
        && verifyResult(result, 4, 4, 0, 0, false, "disabled configuration")
        && require(stats.total_reads == 1 && stats.bases_before == 4
                   && stats.bases_after == 4 && stats.bases_trimmed == 0
                   && stats.quality_trimmed_reads == 0,
                   "disabled configuration statistics");
}

bool testFixedTrimmingDoesNotCountAsQualityTrimming() {
    TrimConfig config;
    config.enabled = true;
    config.trim_front = 1;
    config.trim_tail = 1;
    config.cut_front = true;
    config.cut_tail = true;
    config.quality_threshold = 20;

    Trimmer trimmer(config);
    FastqRecord record = makeRecord("ABCD", {30, 30, 30, 30});
    const TrimResult result = trimmer.trim(record);
    const TrimStats& stats = trimmer.getStats();

    return verifyRecord(record, "BC", "??", "fixed trimming statistics")
        && verifyResult(result, 4, 2, 1, 1, false,
                        "fixed trimming statistics")
        && require(stats.bases_trimmed == 2 && stats.quality_trimmed_reads == 0,
                   "fixed trimming must not count as quality trimming");
}

}  // namespace

int main() {
    return testFixedFront()
        && testFixedTail()
        && testFixedFrontAndTail()
        && testQualityFront()
        && testQualityTail()
        && testTrimmerStatsAndPairedReads()
        && testDisabledConfiguration()
        && testFixedTrimmingDoesNotCountAsQualityTrimming()
        ? 0 : 1;
}
