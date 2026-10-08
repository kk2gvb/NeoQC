#include "fastq_reader.h"
#include "trimming/trimmer.h"

#include <iostream>
#include <string>

namespace {

FastqRecord makeRecord(const std::string& sequence) {
    FastqRecord record;
    record.header = "@READ";
    record.sequence = sequence;
    record.separator = "+";
    record.quality.assign(sequence.size(), 'I');
    return record;
}

bool check(bool condition, const std::string& message) {
    if (!condition) {
        std::cerr << message << '\n';
    }
    return condition;
}

bool checkPolyG() {
    TrimConfig config;
    config.enabled = true;
    config.trim_poly_g = true;
    config.poly_g_min_length = 4;

    struct Case {
        std::string input;
        std::string expected;
        bool trimmed;
    };
    const Case cases[] = {
        {"ACGT", "ACGT", false},
        {"ACGTGGG", "ACGTGGG", false},
        {"ACGTGGGG", "ACGT", true},
        {"ACGTGGGGG", "ACGT", true},
        {"ACGGGGTA", "ACGGGGTA", false},
        {"GGGGGG", "", true},
        {"", "", false},
        {"G", "G", false},
    };

    Trimmer trimmer(config);
    for (const auto& test : cases) {
        FastqRecord record = makeRecord(test.input);
        const TrimResult result = trimmer.trim(record);
        if (!check(record.sequence == test.expected,
                   "unexpected polyG sequence for " + test.input)
            || !check(record.quality.size() == record.sequence.size(),
                      "polyG sequence/quality desynchronization")
            || !check(result.polyG_trimmed == test.trimmed,
                      "unexpected polyG TrimResult flag")
            || !check(!result.polyX_trimmed,
                      "polyG-only run set polyX flag")) {
            return false;
        }
    }
    return check(trimmer.getStats().polyG_trimmed_reads == 3,
                 "polyG statistics count");
}

bool checkPolyX() {
    for (const char base : std::string("ACGT")) {
        TrimConfig config;
        config.enabled = true;
        config.trim_poly_x = true;
        config.poly_x_min_length = 4;
        Trimmer trimmer(config);
        const std::string prefix = base == 'T' ? "ACGA" : "ACGT";

        FastqRecord exact = makeRecord(prefix + std::string(4, base));
        const TrimResult exactResult = trimmer.trim(exact);
        if (!check(exact.sequence == prefix && exactResult.polyX_trimmed,
                   std::string("polyX threshold case failed for ") + base)) {
            return false;
        }

        FastqRecord longRun = makeRecord(prefix + std::string(5, base));
        if (!check(trimmer.trim(longRun).polyX_trimmed
                   && longRun.sequence == prefix,
                   std::string("polyX long run failed for ") + base)) {
            return false;
        }

        FastqRecord shortRun = makeRecord(prefix + std::string(3, base));
        if (!check(!trimmer.trim(shortRun).polyX_trimmed,
                   std::string("polyX short run trimmed for ") + base)) {
            return false;
        }
    }

    TrimConfig config;
    config.enabled = true;
    config.trim_poly_x = true;
    config.poly_x_min_length = 4;
    Trimmer trimmer(config);
    FastqRecord middle = makeRecord("ACAAAAGT");
    FastqRecord nTail = makeRecord("ACGTNNNN");
    FastqRecord empty = makeRecord("");
    if (!check(!trimmer.trim(middle).polyX_trimmed
               && middle.sequence == "ACAAAAGT",
               "middle homopolymer was trimmed")
        || !check(!trimmer.trim(nTail).polyX_trimmed
                  && nTail.sequence == "ACGTNNNN",
                  "N tail was treated as polyX")
        || !check(!trimmer.trim(empty).polyX_trimmed,
                  "empty read was marked polyX-trimmed")) {
        return false;
    }
    return true;
}

bool checkCombinedAndFiltering() {
    TrimConfig config;
    config.enabled = true;
    config.trim_poly_g = true;
    config.trim_poly_x = true;
    config.poly_g_min_length = 4;
    config.poly_x_min_length = 4;
    config.min_length = 3;
    Trimmer trimmer(config);

    FastqRecord gTail = makeRecord("ACGTGGGG");
    const TrimResult gResult = trimmer.trim(gTail);
    if (!check(gResult.polyG_trimmed && !gResult.polyX_trimmed,
               "polyG/polyX double counting")
        || !check(gTail.sequence == "ACGT", "combined polyG output")) {
        return false;
    }

    FastqRecord tooShort = makeRecord("ACGGGG");
    const TrimResult shortResult = trimmer.trim(tooShort);
    return check(shortResult.polyG_trimmed
                 && !shortResult.passed
                 && shortResult.discard_reason == "TOO_SHORT"
                 && tooShort.sequence == "AC",
                 "minimum length was not applied after polyG");
}

}  // namespace

int main() {
    return checkPolyG() && checkPolyX() && checkCombinedAndFiltering() ? 0 : 1;
}
