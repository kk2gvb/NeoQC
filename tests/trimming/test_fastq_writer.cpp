#include "fastq_reader.h"
#include "trimming/fastq_writer.h"

#include <filesystem>
#include <iostream>
#include <string>

namespace {

FastqRecord makeRecord(const std::string& header,
                       const std::string& sequence,
                       const std::string& quality) {
    FastqRecord record;
    record.header = header;
    record.sequence = sequence;
    record.separator = "+ original metadata";
    record.quality = quality;
    return record;
}

}  // namespace

int main() {
    const std::filesystem::path path =
        std::filesystem::temp_directory_path() / "neoqc_fastq_writer_test.fastq.gz";
    const FastqRecord expected1 = makeRecord("@READ_001/1 metadata", "ACGT", "IIII");
    const FastqRecord expected2 = makeRecord("@READ_002/1", "GGA", "HHH");

    {
        FastqWriter writer(path);
        writer.write(expected1);
        writer.write(expected2);
        writer.close();
    }

    FastqReader reader(path.string());
    FastqRecord actual1;
    FastqRecord actual2;
    FastqRecord end;
    const bool valid = reader.readNext(actual1)
        && reader.readNext(actual2)
        && !reader.readNext(end)
        && actual1.header == expected1.header
        && actual1.sequence == expected1.sequence
        && actual1.separator == expected1.separator
        && actual1.quality == expected1.quality
        && actual2.header == expected2.header
        && actual2.sequence == expected2.sequence
        && actual2.quality == expected2.quality;
    if (!valid) {
        std::cerr << "gzip FASTQ round-trip failed\n";
        return 1;
    }

    const std::filesystem::path bufferedPath =
        std::filesystem::temp_directory_path() / "neoqc_fastq_writer_buffered_test.fastq.gz";
    {
        FastqWriter writer(bufferedPath);
        for (std::size_t index = 0; index < 5003; ++index) {
            const std::string id = std::to_string(index);
            writer.write(makeRecord("@ORDER_" + id, "ACGTACGT", "IIIIIIII"));
        }
        writer.close();
    }
    FastqReader bufferedReader(bufferedPath.string());
    for (std::size_t index = 0; index < 5003; ++index) {
        FastqRecord actual;
        const std::string id = std::to_string(index);
        if (!bufferedReader.readNext(actual)
            || actual.header != "@ORDER_" + id
            || actual.sequence != "ACGTACGT"
            || actual.separator != "+ original metadata"
            || actual.quality != "IIIIIIII") {
            std::cerr << "buffered gzip FASTQ order or formatting failed\n";
            return 1;
        }
    }
    FastqRecord bufferedEnd;
    if (bufferedReader.readNext(bufferedEnd)) {
        std::cerr << "buffered gzip FASTQ had unexpected extra records\n";
        return 1;
    }

    const std::filesystem::path emptyPath =
        std::filesystem::temp_directory_path() / "neoqc_fastq_writer_empty_test.fastq.gz";
    {
        FastqWriter writer(emptyPath);
        writer.close();
    }
    gzFile empty = gzopen(emptyPath.string().c_str(), "rb");
    char byte = 0;
    const bool emptyOutput = empty != nullptr && gzread(empty, &byte, 1) == 0;
    if (empty != nullptr) {
        gzclose(empty);
    }

    const std::filesystem::path fullPath("/dev/full");
    bool propagatedWriteError = false;
    if (std::filesystem::exists(fullPath)) {
        try {
            FastqWriter writer(fullPath);
            writer.write(expected1);
            writer.close();
        } catch (const std::runtime_error&) {
            propagatedWriteError = true;
        }
    } else {
        propagatedWriteError = true;
    }

    std::filesystem::remove(path);
    std::filesystem::remove(bufferedPath);
    std::filesystem::remove(emptyPath);
    if (!emptyOutput || !propagatedWriteError) {
        std::cerr << "buffered gzip FASTQ empty output or error propagation failed\n";
        return 1;
    }
    return 0;
}
