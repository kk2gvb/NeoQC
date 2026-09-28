#include "fastq_reader.h"
#include "trimming/adapter_fasta.h"
#include "trimming/adapter_trimmer.h"
#include "trimming/pe_overlap.h"
#include "trimming/pe_processor.h"
#include "trimming/quality_trimmer.h"
#include "trimming/trimmer.h"
#include "trimming/umi_processor.h"

#include <filesystem>
#include <fstream>
#include <iostream>
#include <optional>
#include <stdexcept>
#include <string>

namespace {

FastqRecord record(std::string header,
                   std::string sequence,
                   std::string quality = {}) {
    if (quality.empty()) quality.assign(sequence.size(), 'I');
    return {std::move(header), std::move(sequence), "+", std::move(quality), 7};
}

bool check(bool condition, const std::string& message) {
    if (!condition) std::cerr << message << '\n';
    return condition;
}

TrimConfig peConfig() {
    TrimConfig config;
    config.enabled = true;
    config.min_overlap = 4;
    config.overlap_diff_limit = 2;
    config.overlap_diff_percent_limit = 25;
    return config;
}

bool testOverlapCoordinatesAndCorrection() {
    TrimConfig config = peConfig();
    FastqRecord r1 = record("@PAIR/1 metadata", "AAAACCCC", "IIII!III");
    const std::string r2Reverse = "TCCCGGGG";
    FastqRecord r2 = record("@PAIR/2 metadata",
                            PeOverlapAnalyzer::reverseComplement(r2Reverse));
    const auto overlap = PeOverlapAnalyzer(config).find(r1, r2);
    if (!check(overlap.has_value(), "partial overlap must be found")
        || !check(overlap->r1_start == 4 && overlap->r2_reverse_start == 0
                  && overlap->length == 4 && overlap->mismatches == 1,
                 "partial overlap coordinates")) return false;

    // reverse index 0 maps to the last original R2 base.
    r2.quality.back() = 'I';
    TrimResult result1;
    TrimResult result2;
    PeProcessor().correct(r1, r2, result1, result2, *overlap);
    if (!check(r1.sequence == "AAAATCCC", "higher-quality R2 correction")
        || !check(r1.quality[4] == 'I', "corrected quality comes from R2")
        || !check(result1.corrected_bases == 1 && result2.corrected_bases == 0,
                 "correction counters")
        || !check(r1.header == "@PAIR/1 metadata" && r2.header == "@PAIR/2 metadata",
                 "correction preserves IDs")
        || !check(r1.sequence.size() == r1.quality.size()
                  && r2.sequence.size() == r2.quality.size(),
                 "correction preserves FASTQ lengths")) return false;

    r1 = record("@PAIR/1", "AAAACCCC");
    r2 = record("@PAIR/2", PeOverlapAnalyzer::reverseComplement(r2Reverse));
    r2.quality.back() = '!';
    result1 = {};
    result2 = {};
    PeProcessor().correct(r1, r2, result1, result2, *overlap);
    if (!check(PeOverlapAnalyzer::reverseComplement(r2.sequence) == "CCCCGGGG",
               "higher-quality R1 correction at reverse R2 coordinate")
        || !check(result2.corrected_bases == 1, "R2 correction count")) return false;

    r1 = record("@PAIR/1", "AAAACCCC");
    r2 = record("@PAIR/2", PeOverlapAnalyzer::reverseComplement(r2Reverse));
    result1 = {};
    result2 = {};
    PeProcessor().correct(r1, r2, result1, result2, *overlap);
    if (!check(r1.sequence == "AAAACCCC"
               && PeOverlapAnalyzer::reverseComplement(r2.sequence) == r2Reverse,
               "equal qualities must not correct")) return false;

    r1 = record("@PAIR/1", "AAAANCCC", "IIII!III");
    r2 = record("@PAIR/2", PeOverlapAnalyzer::reverseComplement(r2Reverse));
    result1 = {};
    result2 = {};
    PeProcessor().correct(r1, r2, result1, result2, *overlap);
    if (!check(r1.sequence == "AAAATCCC", "N to canonical correction")) return false;

    TrimConfig twoMismatchConfig = peConfig();
    twoMismatchConfig.min_overlap = 8;
    r1 = record("@PAIR/1", "ACGTACGT", "!III!III");
    r2 = record("@PAIR/2", PeOverlapAnalyzer::reverseComplement("TCGTTCGT"));
    const auto twoMismatchOverlap = PeOverlapAnalyzer(twoMismatchConfig).find(r1, r2);
    result1 = {};
    result2 = {};
    PeProcessor().correct(r1, r2, result1, result2, *twoMismatchOverlap);
    if (!check(r1.sequence == "TCGTTCGT" && result1.corrected_bases == 2,
               "multiple mismatch correction")) return false;

    r1 = record("@PAIR/1", "AAAANCCC");
    r2 = record("@PAIR/2", PeOverlapAnalyzer::reverseComplement("NCCCGGGG"));
    result1 = {};
    result2 = {};
    const auto nOverlap = PeOverlapAnalyzer(config).find(r1, r2);
    PeProcessor().correct(r1, r2, result1, result2, *nOverlap);
    if (!check(result1.corrected_bases == 0 && result2.corrected_bases == 0,
               "N versus N is not corrected")) return false;

    config.overlap_correction = true;
    Trimmer trimmer(config);
    r1 = record("@PAIR/1", "AAAACCCC", "IIII!III");
    r2 = record("@PAIR/2", PeOverlapAnalyzer::reverseComplement(r2Reverse));
    trimmer.trimPair(r1, r2);
    if (!check(trimmer.getStats().corrected_reads == 1
               && trimmer.getStats().corrected_bases == 1,
               "aggregate correction statistics")) return false;

    FastqRecord noR1 = record("@NONE/1", "AAAAAAAA");
    FastqRecord noR2 = record("@NONE/2", "AAAAAAAA");
    Trimmer noOverlapTrimmer(config);
    noOverlapTrimmer.trimPair(noR1, noR2);
    return check(noOverlapTrimmer.getStats().corrected_bases == 0,
                 "invalid overlap forbids correction");
}

bool testMerge() {
    TrimConfig config = peConfig();
    FastqRecord r1 = record("@PAIR/1 source", "AAAACCCC");
    FastqRecord r2 = record("@PAIR/2 source",
                            PeOverlapAnalyzer::reverseComplement("CCCCGGGG"));
    const auto overlap = PeOverlapAnalyzer(config).find(r1, r2);
    const auto merged = PeProcessor().merge(r1, r2, *overlap);
    if (!check(merged->sequence == "AAAACCCCGGGG", "partial merge sequence")
        || !check(merged->quality == std::string(12, 'I'), "partial merge quality")
        || !check(merged->header == "@PAIR merged", "deterministic merged ID")
        || !check(merged->sequence.size() == merged->quality.size(),
                  "merged FASTQ lengths")
        || !check(r1.sequence == "AAAACCCC", "merge does not mutate R1")) return false;

    r1 = record("@FULL/1", "ACGT");
    r2 = record("@FULL/2", PeOverlapAnalyzer::reverseComplement("ACGT"));
    const auto fullOverlap = PeOverlapAnalyzer(config).find(r1, r2);
    if (!check(PeProcessor().merge(r1, r2, *fullOverlap)->sequence == "ACGT",
               "full overlap merge")) return false;

    r1 = record("@MISMATCH/1", "AAAACCCC");
    r2 = record("@MISMATCH/2", PeOverlapAnalyzer::reverseComplement("TCCCGGGG"));
    const auto mismatchOverlap = PeOverlapAnalyzer(config).find(r1, r2);
    if (!check(PeProcessor().merge(r1, r2, *mismatchOverlap)->sequence
                   == "AAAANCCCGGGG",
               "equal-quality disagreement becomes N")) return false;

    r1.quality[4] = '!';
    const auto qualityConsensus = PeProcessor().merge(r1, r2, *mismatchOverlap);
    if (!check(qualityConsensus->sequence == "AAAATCCCGGGG"
               && qualityConsensus->quality[4] == 'I',
               "merge selects the higher-quality base and quality")) return false;

    r1 = record("@SHORT/1", "ACGTACAAAA");
    r2 = record("@SHORT/2", "GTACGTCCCC");
    const auto readThroughOverlap = PeOverlapAnalyzer(config).find(r1, r2);
    if (!check(readThroughOverlap.has_value() && readThroughOverlap->isReadThrough()
               && PeProcessor().merge(r1, r2, *readThroughOverlap)->sequence == "ACGTAC",
               "read-through merge excludes adapter tails")) return false;

    FastqRecord noR1 = record("@NONE/1", "AAAAAAAA");
    FastqRecord noR2 = record("@NONE/2", "AAAAAAAA");
    if (!check(!PeOverlapAnalyzer(config).find(noR1, noR2).has_value(),
               "invalid overlap must not merge")) return false;

    config.merge_reads = true;
    Trimmer trimmer(config);
    r1 = record("@PAIR/1", "AAAACCCC");
    r2 = record("@PAIR/2", PeOverlapAnalyzer::reverseComplement("CCCCGGGG"));
    PairTrimResult output = trimmer.trimPairDetailed(r1, r2);
    return check(output.merged.has_value()
                 && trimmer.getStats().merged_pairs == 1,
                 "Trimmer merge orchestration and statistics");
}

bool testUmi() {
    TrimConfig config;
    config.enabled = true;
    config.umi_enabled = true;
    config.umi_length = 4;
    UmiProcessor processor(config);
    FastqRecord read = record("@READ metadata", "ACGTGATT", "12345678");
    TrimResult result;
    if (!check(processor.process(read, result), "valid UMI")
        || !check(read.sequence == "GATT" && read.quality == "5678",
                 "UMI removes synchronized prefix")
        || !check(read.header == "@READ metadata UMI:ACGT", "UMI header format")
        || !check(result.umi == "ACGT" && result.trimmed_front == 4,
                 "UMI result metadata")) return false;

    config.umi_length = 8;
    read = record("@READ", "ACGTGATT");
    result = {};
    if (!check(UmiProcessor(config).process(read, result) && read.sequence.empty()
               && read.quality.empty(), "UMI equal to read length")) return false;

    config.umi_length = 9;
    read = record("@READ", "ACGTGATT");
    result = {};
    if (!check(!UmiProcessor(config).process(read, result)
               && result.discard_reason == "UMI_TOO_SHORT"
               && read.sequence == "ACGTGATT", "oversized UMI fails without mutation")) {
        return false;
    }

    bool zeroRejected = false;
    try {
        config.umi_length = 0;
        UmiProcessor(config).process(read, result);
    } catch (const std::invalid_argument&) {
        zeroRejected = true;
    }
    if (!check(zeroRejected, "zero UMI length rejected")) return false;

    config.umi_length = 1;
    config.trim_front = 1;
    Trimmer trimmer(config);
    read = record("@ORDER", "AACG", "IIII");
    result = trimmer.trim(read);
    if (!check(read.sequence == "CG" && result.umi == "A"
               && result.trimmed_front == 2,
               "UMI runs before fixed trimming")) return false;

    config.trim_front = 0;
    Trimmer pairTrimmer(config);
    FastqRecord r1 = record("@UMI_PAIR/1", "AACG");
    FastqRecord r2 = record("@UMI_PAIR/2", "TGCA");
    pairTrimmer.trimPair(r1, r2);
    return check(r1.sequence == "ACG" && r1.header == "@UMI_PAIR/1 UMI:A"
                 && r2.sequence == "TGCA" && r2.header == "@UMI_PAIR/2",
                 "PE UMI is extracted from R1 only");
}

bool testAdapterFasta() {
    const std::filesystem::path path =
        std::filesystem::temp_directory_path() / "neoqc-adapters-test.fa";
    {
        std::ofstream output(path);
        output << ">first\nAGAT\nCGGA\n>second\nACGTACGT\n";
    }
    const auto adapters = loadAdapterFasta(path);
    std::filesystem::remove(path);
    if (!check(adapters.size() == 2 && adapters[0].sequence == "AGATCGGA",
               "multiline adapter FASTA")) return false;

    TrimConfig config;
    config.enabled = true;
    config.adapter_trimming = true;
    config.fasta_adapters = adapters;
    FastqRecord read = record("@A", "TTTTACGTACGT");
    const TrimResult result = AdapterTrimmer(config).trim(read);
    if (!check(read.sequence == "TTTT" && result.adapter_position == 4,
               "multiple FASTA adapters")) return false;

    config.adapter_sequence = "CCCCCCCC";
    read = record("@A", "GGGGAGATCGGA");
    if (!check(AdapterTrimmer(config).trim(read).adapter_found
               && read.sequence == "GGGG",
               "explicit and FASTA adapters form one collection")) return false;

    bool malformedRejected = false;
    {
        std::ofstream output(path);
        output << "ACGT\n";
    }
    try {
        (void)loadAdapterFasta(path);
    } catch (const std::runtime_error&) {
        malformedRejected = true;
    }
    std::filesystem::remove(path);
    bool missingRejected = false;
    try {
        (void)loadAdapterFasta(path);
    } catch (const std::runtime_error&) {
        missingRejected = true;
    }
    bool emptyRejected = false;
    {
        std::ofstream output(path);
    }
    try {
        (void)loadAdapterFasta(path);
    } catch (const std::runtime_error&) {
        emptyRejected = true;
    }
    std::filesystem::remove(path);
    bool alphabetRejected = false;
    {
        std::ofstream output(path);
        output << ">invalid\nACGTZ\n";
    }
    try {
        (void)loadAdapterFasta(path);
    } catch (const std::runtime_error&) {
        alphabetRejected = true;
    }
    std::filesystem::remove(path);
    return check(malformedRejected && missingRejected && emptyRejected
                 && alphabetRejected,
                 "malformed, empty, invalid, and missing adapter FASTA rejected");
}

bool testQualityBoundaries() {
    TrimConfig config;
    config.enabled = true;
    config.cut_front = true;
    config.quality_threshold = 0;
    FastqRecord read = record("@Q0", "A", "!");
    TrimResult result = QualityTrimmer(config).trim(read);
    if (!check(read.sequence == "A" && result.final_length == 1,
               "quality threshold zero includes Q0")) return false;

    config.quality_threshold = 93;
    read = record("@Q93", "A", "~");
    result = QualityTrimmer(config).trim(read);
    if (!check(read.sequence == "A", "maximum Phred threshold boundary")) return false;
    read = record("@Q92", "A", "}");
    result = QualityTrimmer(config).trim(read);
    return check(read.sequence.empty() && result.quality_trimmed,
                 "base below maximum Phred threshold is trimmed");
}

}  // namespace

int main() {
    return testOverlapCoordinatesAndCorrection()
        && testMerge()
        && testUmi()
        && testAdapterFasta()
        && testQualityBoundaries()
        ? 0 : 1;
}
