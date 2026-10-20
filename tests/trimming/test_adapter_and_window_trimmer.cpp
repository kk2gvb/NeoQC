#include "fastq_reader.h"
#include "trimming/adapter_trimmer.h"
#include "trimming/quality_trimmer.h"
#include "trimming/trimmer.h"

#include <initializer_list>
#include <iostream>
#include <stdexcept>
#include <string>

namespace {

FastqRecord makeRecord(const std::string& sequence,
                       std::initializer_list<int> qualities = {}) {
    FastqRecord record;
    record.header = "@PAIR_001/1";
    record.separator = "+";
    record.recordNumber = 17;
    record.sequence = sequence;
    if (qualities.size() == 0) {
        record.quality.assign(sequence.size(), '?');
    } else {
        for (const int quality : qualities) {
            record.quality.push_back(static_cast<char>(quality + 33));
        }
    }
    return record;
}

bool check(bool condition, const std::string& message) {
    if (!condition) std::cerr << message << '\n';
    return condition;
}

bool checkRecord(const FastqRecord& record,
                 const std::string& sequence,
                 const std::string& message) {
    return check(record.sequence == sequence, message + ": sequence")
        && check(record.sequence.size() == record.quality.size(),
                 message + ": sequence/quality synchronization")
        && check(record.header == "@PAIR_001/1" && record.separator == "+"
                 && record.recordNumber == 17,
                 message + ": FASTQ metadata");
}

bool checkResult(const TrimResult& result,
                 std::size_t original,
                 std::size_t final,
                 std::size_t front,
                 std::size_t tail,
                 bool quality,
                 bool adapter,
                 std::optional<std::size_t> adapterPosition,
                 const std::string& message) {
    return check(result.original_length == original, message + ": original length")
        && check(result.final_length == final, message + ": final length")
        && check(result.trimmed_front == front, message + ": front trimming")
        && check(result.trimmed_tail == tail, message + ": tail trimming")
        && check(result.quality_trimmed == quality, message + ": quality marker")
        && check(result.adapter_found == adapter, message + ": adapter marker")
        && check(result.adapter_position == adapterPosition,
                 message + ": adapter position");
}

std::string reverseComplement(const std::string& sequence) {
    std::string result;
    result.reserve(sequence.size());
    for (auto iterator = sequence.rbegin(); iterator != sequence.rend(); ++iterator) {
        switch (*iterator) {
            case 'A': result.push_back('T'); break;
            case 'C': result.push_back('G'); break;
            case 'G': result.push_back('C'); break;
            case 'T': result.push_back('A'); break;
            default: result.push_back('N'); break;
        }
    }
    return result;
}

bool testSlidingWindow() {
    TrimConfig config;
    config.enabled = true;
    config.cut_right = true;
    config.window_size = 4;
    config.quality_threshold = 20;

    FastqRecord record = makeRecord("ABCDEFGHI", {30, 30, 30, 30, 30, 30, 30, 30, 30});
    TrimResult result = QualityTrimmer(config).trim(record);
    if (!checkRecord(record, "ABCDEFGHI", "sliding all good")
        || !checkResult(result, 9, 9, 0, 0, false, false, std::nullopt,
                        "sliding all good")) return false;

    record = makeRecord("ABCDEFGHI", {10, 10, 10, 10, 30, 30, 30, 30, 30});
    result = QualityTrimmer(config).trim(record);
    if (!checkRecord(record, "", "sliding first window low")
        || !checkResult(result, 9, 0, 0, 9, true, false, std::nullopt,
                        "sliding first window low")) return false;

    record = makeRecord("ABCDEFGHI", {30, 30, 30, 30, 30, 10, 10, 10, 10});
    result = QualityTrimmer(config).trim(record);
    if (!checkRecord(record, "ABCD", "sliding middle window low")
        || !checkResult(result, 9, 4, 0, 5, true, false, std::nullopt,
                        "sliding middle window low")) return false;

    record = makeRecord("ABCDEFGHI", {30, 30, 30, 30, 30, 30, 10, 10, 10});
    result = QualityTrimmer(config).trim(record);
    if (!checkRecord(record, "ABCDE", "sliding final window low")
        || !checkResult(result, 9, 5, 0, 4, true, false, std::nullopt,
                        "sliding final window low")) return false;

    record = makeRecord("ABCD", {20, 20, 20, 20});
    result = QualityTrimmer(config).trim(record);
    if (!checkRecord(record, "ABCD", "sliding threshold boundary")
        || !checkResult(result, 4, 4, 0, 0, false, false, std::nullopt,
                        "sliding threshold boundary")) return false;

    config.window_size = 1;
    record = makeRecord("ABC", {30, 10, 30});
    result = QualityTrimmer(config).trim(record);
    if (!checkRecord(record, "A", "sliding one-base window")
        || !checkResult(result, 3, 1, 0, 2, true, false, std::nullopt,
                        "sliding one-base window")) return false;

    config.window_size = 3;
    record = makeRecord("ABC", {10, 10, 10});
    result = QualityTrimmer(config).trim(record);
    if (!checkRecord(record, "", "sliding whole-read window")
        || !checkResult(result, 3, 0, 0, 3, true, false, std::nullopt,
                        "sliding whole-read window")) return false;

    config.window_size = 4;
    record = makeRecord("ABC", {10, 10, 10});
    result = QualityTrimmer(config).trim(record);
    if (!checkRecord(record, "ABC", "sliding window exceeds read")
        || !checkResult(result, 3, 3, 0, 0, false, false, std::nullopt,
                        "sliding window exceeds read")) return false;

    record = makeRecord("", {});
    result = QualityTrimmer(config).trim(record);
    if (!checkRecord(record, "", "sliding empty read")
        || !checkResult(result, 0, 0, 0, 0, false, false, std::nullopt,
                        "sliding empty read")) return false;

    config.window_size = 1;
    record = makeRecord("A", {10});
    result = QualityTrimmer(config).trim(record);
    if (!checkRecord(record, "", "sliding one-base read")
        || !checkResult(result, 1, 0, 0, 1, true, false, std::nullopt,
                        "sliding one-base read")) return false;

    config.window_size = 0;
    bool rejectedInvalidWindow = false;
    try {
        record = makeRecord("A", {30});
        QualityTrimmer(config).trim(record);
    } catch (const std::invalid_argument&) {
        rejectedInvalidWindow = true;
    }
    return check(rejectedInvalidWindow, "sliding zero window must be rejected");
}

bool testSlidingStats() {
    TrimConfig config;
    config.enabled = true;
    config.cut_front = true;
    config.cut_right = true;
    config.quality_threshold = 20;
    config.window_size = 1;
    Trimmer trimmer(config);
    FastqRecord record = makeRecord("ABCD", {5, 30, 10, 30});
    const TrimResult result = trimmer.trim(record);
    const TrimStats& stats = trimmer.getStats();
    return checkRecord(record, "B", "combined quality operations")
        && checkResult(result, 4, 1, 1, 2, true, false, std::nullopt,
                       "combined quality operations")
        && check(stats.quality_trimmed_reads == 1 && stats.bases_trimmed == 3,
                 "sliding quality statistics must not double-count reads");
}

bool testAdapterMatching() {
    TrimConfig config;
    config.enabled = true;
    config.adapter_trimming = true;
    config.adapter_sequence = "AGATCGGA";
    AdapterTrimmer trimmer(config);

    FastqRecord record = makeRecord("TTTTAGATCGGA");
    TrimResult result = trimmer.trim(record);
    if (!checkRecord(record, "TTTT", "full adapter")
        || !checkResult(result, 12, 4, 0, 8, false, true, 4,
                        "full adapter")) return false;

    record = makeRecord("TTTTAGATCG");
    result = trimmer.trim(record);
    if (!checkRecord(record, "TTTT", "partial adapter")
        || !checkResult(result, 10, 4, 0, 6, false, true, 4,
                        "partial adapter")) return false;

    record = makeRecord("TTTTCCCCGGGG");
    result = trimmer.trim(record);
    if (!checkRecord(record, "TTTTCCCCGGGG", "no adapter")
        || !checkResult(result, 12, 12, 0, 0, false, false, std::nullopt,
                        "no adapter")) return false;

    record = makeRecord("AGATCGGA");
    result = trimmer.trim(record);
    if (!checkRecord(record, "", "adapter at start")
        || !checkResult(result, 8, 0, 0, 8, false, true, 0,
                        "adapter at start")) return false;

    record = makeRecord("TTTTAGATCG");
    result = trimmer.trim(record);
    if (!checkRecord(record, "TTTT", "exact minimum adapter match")
        || !checkResult(result, 10, 4, 0, 6, false, true, 4,
                        "exact minimum adapter match")) return false;

    record = makeRecord("TTTTAGATC");
    result = trimmer.trim(record);
    if (!checkRecord(record, "TTTTAGATC", "short adapter prefix")
        || !checkResult(result, 9, 9, 0, 0, false, false, std::nullopt,
                        "short adapter prefix")) return false;

    record = makeRecord("TTTTAGATCGTA");
    result = trimmer.trim(record);
    if (!checkRecord(record, "TTTT", "one mismatch at length eight")
        || !checkResult(result, 12, 4, 0, 8, false, true, 4,
                        "one mismatch at length eight")) return false;

    config.adapter_sequence = "ACGTACGTACGTACGT";
    trimmer = AdapterTrimmer(config);
    record = makeRecord("TTTTACGTTCGTTCGTACGT");
    result = trimmer.trim(record);
    if (!checkRecord(record, "TTTT", "allowed mismatch limit")
        || !check(result.adapter_found, "allowed mismatch limit marker")) return false;

    config.adapter_sequence = "AAAAAA";
    trimmer = AdapterTrimmer(config);
    record = makeRecord("CCCCCC");
    result = trimmer.trim(record);
    if (!checkRecord(record, "CCCCCC", "excess mismatch limit")
        || !check(!result.adapter_found, "excess mismatch limit marker")) return false;

    config.adapter_sequence.clear();
    trimmer = AdapterTrimmer(config);
    record = makeRecord("AGATCGGA");
    result = trimmer.trim(record);
    if (!checkRecord(record, "AGATCGGA", "empty adapter")
        || !check(!result.adapter_found, "empty adapter marker")) return false;

    config.adapter_sequence = "AGATCGGA";
    trimmer = AdapterTrimmer(config);
    record = makeRecord("", {});
    result = trimmer.trim(record);
    return checkRecord(record, "", "empty adapter read")
        && check(!result.adapter_found, "empty adapter read marker");
}

bool testAdapterStats() {
    TrimConfig config;
    config.enabled = true;
    config.adapter_trimming = true;
    config.adapter_sequence = "AGATCGGA";
    Trimmer trimmer(config);
    FastqRecord record = makeRecord("TTTTAGATCGGA");
    trimmer.trim(record);
    const TrimStats& stats = trimmer.getStats();
    return check(stats.adapter_trimmed_reads == 1 && stats.total_reads == 1
                 && stats.bases_trimmed == 8,
                 "adapter statistics");
}

bool testPairedOverlap() {
    const std::string insert = "ACGTTGCA";
    const std::string reverseInsert = reverseComplement(insert);
    TrimConfig config;
    config.enabled = true;
    config.adapter_trimming = true;
    config.min_overlap = 8;
    config.overlap_diff_limit = 1;
    config.overlap_diff_percent_limit = 20;

    Trimmer trimmer(config);
    FastqRecord r1 = makeRecord(insert + "AAAA");
    FastqRecord r2 = makeRecord(reverseInsert + "CCCC");
    const auto results = trimmer.trimPair(r1, r2);
    if (!checkRecord(r1, insert, "paired short insert R1")
        || !checkRecord(r2, reverseInsert, "paired short insert R2")
        || !checkResult(results.first, 12, 8, 0, 4, false, true, 8,
                        "paired short insert R1")
        || !checkResult(results.second, 12, 8, 0, 4, false, true, 8,
                        "paired short insert R2")) return false;

    Trimmer normalTrimmer(config);
    r1 = makeRecord("ACGTTGCATGCA");
    r2 = makeRecord(reverseComplement(r1.sequence));
    const auto normalResults = normalTrimmer.trimPair(r1, r2);
    if (!checkRecord(r1, "ACGTTGCATGCA", "normal paired reads R1")
        || !checkRecord(r2, reverseComplement("ACGTTGCATGCA"), "normal paired reads R2")
        || !check(!normalResults.first.adapter_found && !normalResults.second.adapter_found,
                  "normal paired reads must not trim")) return false;

    Trimmer noOverlapTrimmer(config);
    r1 = makeRecord("AAAAAAAAAAAA");
    r2 = makeRecord("AAAAAAAAAAAA");
    const auto noOverlapResults = noOverlapTrimmer.trimPair(r1, r2);
    if (!check(!noOverlapResults.first.adapter_found && !noOverlapResults.second.adapter_found,
               "non-overlapping reads must remain unchanged")) return false;

    r1 = makeRecord(insert + "AAAA");
    r2 = makeRecord(reverseInsert + "CCCA");
    const auto mismatchResults = trimmer.trimPair(r1, r2);
    if (!check(mismatchResults.first.adapter_found && mismatchResults.second.adapter_found,
               "permitted overlap mismatch")) return false;

    r1 = makeRecord(insert + "AAAA");
    r2 = makeRecord("AA" + reverseInsert.substr(2) + "CCCC");
    const auto excessiveMismatchResults = trimmer.trimPair(r1, r2);
    if (!check(!excessiveMismatchResults.first.adapter_found
               && !excessiveMismatchResults.second.adapter_found,
               "excess overlap mismatches")) return false;

    TrimConfig percentageConfig = config;
    percentageConfig.overlap_diff_limit = 5;
    percentageConfig.overlap_diff_percent_limit = 10;
    Trimmer percentageTrimmer(percentageConfig);
    r1 = makeRecord(insert + "AAAA");
    r2 = makeRecord("A" + reverseInsert.substr(1) + "CCCC");
    const auto percentageResults = percentageTrimmer.trimPair(r1, r2);
    if (!check(!percentageResults.first.adapter_found
               && !percentageResults.second.adapter_found,
               "overlap mismatch percentage")) return false;

    TrimConfig shortOverlapConfig = config;
    shortOverlapConfig.min_overlap = 9;
    Trimmer shortOverlapTrimmer(shortOverlapConfig);
    r1 = makeRecord(insert + "AAAA");
    r2 = makeRecord(reverseInsert + "CCCC");
    const auto shortOverlapResults = shortOverlapTrimmer.trimPair(r1, r2);
    if (!check(!shortOverlapResults.first.adapter_found
               && !shortOverlapResults.second.adapter_found,
               "overlap shorter than configured minimum")) return false;

    TrimConfig fallbackConfig = config;
    fallbackConfig.min_overlap = 30;
    fallbackConfig.adapter_sequence = "AGATCGGA";
    Trimmer fallbackTrimmer(fallbackConfig);
    r1 = makeRecord("TTTTAGATCGGA");
    r2 = makeRecord("GGGGCCCCCCCC");
    const auto r1FallbackResults = fallbackTrimmer.trimPair(r1, r2);
    if (!checkRecord(r1, "TTTT", "paired fallback only R1")
        || !checkRecord(r2, "GGGGCCCCCCCC", "paired fallback only R1 R2")
        || !check(r1FallbackResults.first.adapter_found
                  && !r1FallbackResults.second.adapter_found,
                  "paired fallback only R1 marker")) return false;

    r1 = makeRecord("TTTTCCCCCCCC");
    r2 = makeRecord("GGGGAGATCGGA");
    const auto r2FallbackResults = fallbackTrimmer.trimPair(r1, r2);
    if (!checkRecord(r1, "TTTTCCCCCCCC", "paired fallback only R2 R1")
        || !checkRecord(r2, "GGGG", "paired fallback only R2")
        || !check(!r2FallbackResults.first.adapter_found
                  && r2FallbackResults.second.adapter_found,
                  "paired fallback only R2 marker")) return false;

    r1 = makeRecord("TTTTAGATCGGA");
    r2 = makeRecord("GGGGAGATCGGA");
    const auto fallbackResults = fallbackTrimmer.trimPair(r1, r2);
    if (!checkRecord(r1, "TTTT", "paired fallback R1")
        || !checkRecord(r2, "GGGG", "paired fallback R2")
        || !check(fallbackResults.first.adapter_found && fallbackResults.second.adapter_found,
                  "paired fallback markers")
        || !check(fallbackTrimmer.getStats().adapter_trimmed_reads == 4,
                  "paired adapter statistics")
        || !check(r1.header == "@PAIR_001/1" && r2.header == "@PAIR_001/1",
                  "paired read identifiers")) return false;

    Trimmer emptyPairTrimmer(config);
    r1 = makeRecord("", {});
    r2 = makeRecord("", {});
    const auto emptyResults = emptyPairTrimmer.trimPair(r1, r2);
    if (!checkRecord(r1, "", "empty paired R1")
        || !checkRecord(r2, "", "empty paired R2")
        || !check(!emptyResults.first.adapter_found && !emptyResults.second.adapter_found,
                  "empty paired adapter markers")) return false;

    TrimConfig tieConfig = config;
    tieConfig.min_overlap = 4;
    tieConfig.overlap_diff_limit = 0;
    tieConfig.overlap_diff_percent_limit = 0;
    Trimmer tieTrimmer(tieConfig);
    r1 = makeRecord("AAAAAAAACCCC");
    r2 = makeRecord("TTTTTTTTGGGG");
    const auto tieResults = tieTrimmer.trimPair(r1, r2);
    return checkRecord(r1, "AAAAAAAA", "overlap tie-breaking R1")
        && checkRecord(r2, "TTTTTTTT", "overlap tie-breaking R2")
        && check(tieResults.first.adapter_position == 8
                 && tieResults.second.adapter_position == 8,
                 "overlap tie-breaking must choose the largest overlap");
}

}  // namespace

int main() {
    return testSlidingWindow()
        && testSlidingStats()
        && testAdapterMatching()
        && testAdapterStats()
        && testPairedOverlap()
        ? 0 : 1;
}
