#include "fastq_reader.h"
#include "trimming/fastq_writer.h"

#include <filesystem>
#include <iostream>

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
    std::filesystem::remove(path);
    if (!valid) {
        std::cerr << "gzip FASTQ round-trip failed\n";
        return 1;
    }
    return 0;
}
