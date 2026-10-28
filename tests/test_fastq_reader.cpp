#include "fastq_reader.h"

#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

void writePlain(const std::filesystem::path& path, const std::string& content) {
    std::ofstream output(path, std::ios::binary);
    output.write(content.data(), static_cast<std::streamsize>(content.size()));
}

void writeGzip(const std::filesystem::path& path, const std::string& content) {
    gzFile output = gzopen(path.string().c_str(), "wb");
    if (output == nullptr) {
        throw std::runtime_error("cannot create gzip FASTQ fixture");
    }
    const int written = gzwrite(
        output, content.data(), static_cast<unsigned int>(content.size()));
    const int closeStatus = gzclose(output);
    if (written != static_cast<int>(content.size()) || closeStatus != Z_OK) {
        throw std::runtime_error("cannot create gzip FASTQ fixture");
    }
}

bool check(bool condition, const std::string& message) {
    if (!condition) std::cerr << message << '\n';
    return condition;
}

bool checkValid(const std::filesystem::path& path,
                const std::vector<std::string>& expectedSequences) {
    FastqReader reader(path.string());
    for (const std::string& expected : expectedSequences) {
        FastqRecord record;
        if (!reader.readNext(record) || record.sequence != expected
            || record.quality.size() != expected.size()) {
            return false;
        }
    }
    FastqRecord end;
    return !reader.readNext(end);
}

bool checkInvalid(const std::filesystem::path& path, const std::string& reason) {
    try {
        FastqReader reader(path.string());
        FastqRecord record;
        reader.readNext(record);
    } catch (const std::runtime_error& error) {
        return std::string(error.what()).find(reason) != std::string::npos;
    }
    return false;
}

}  // namespace

int main() {
    const std::filesystem::path root = std::filesystem::temp_directory_path();
    const std::filesystem::path plain = root / "neoqc_reader_valid.fastq";
    const std::filesystem::path gzip = root / "neoqc_reader_valid.fastq.gz";
    const std::filesystem::path crlf = root / "neoqc_reader_crlf.fastq";
    const std::filesystem::path longLine = root / "neoqc_reader_long.fastq";
    const std::filesystem::path invalid = root / "neoqc_reader_invalid.fastq";

    const std::string valid = "@READ_1\nACGT\n+ metadata\nIIII\n"
        "@READ_2\nTGCA\n+\n####";
    writePlain(plain, valid);
    writeGzip(gzip, valid);
    writePlain(crlf, "@CRLF\r\nACGT\r\n+\r\nIIII\r\n");
    writePlain(longLine, "@LONG\n" + std::string(20000, 'A') + "\n+\n"
        + std::string(20000, 'I') + "\n");

    bool success = check(checkValid(plain, {"ACGT", "TGCA"}), "plain FASTQ")
        && check(checkValid(gzip, {"ACGT", "TGCA"}), "gzip FASTQ")
        && check(checkValid(crlf, {"ACGT"}), "CRLF FASTQ")
        && check(checkValid(longLine, {std::string(20000, 'A')}), "long FASTQ line");

    const std::vector<std::pair<std::string, std::string>> invalidCases = {
        {"", "FASTQ file is empty"},
        {"@READ\nACGT\n+\n", "missing quality"},
        {"@READ\nACGT\n+\nIII\n", "sequence and quality lengths differ"},
        {"@READ\nACGT\n+\nII I\n", "invalid Phred+33 quality character"},
        {"READ\nACGT\n+\nIIII\n", "header does not start with '@'"},
        {"@READ\nACGT\n-\nIIII\n", "separator line does not start with '+'"},
    };
    for (const auto& [content, reason] : invalidCases) {
        writePlain(invalid, content);
        success = check(checkInvalid(invalid, reason), "FASTQ validation: " + reason)
            && success;
    }

    std::filesystem::remove(plain);
    std::filesystem::remove(gzip);
    std::filesystem::remove(crlf);
    std::filesystem::remove(longLine);
    std::filesystem::remove(invalid);
    return success ? 0 : 1;
}
